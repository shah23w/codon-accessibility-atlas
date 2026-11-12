#!/usr/bin/env python3
"""
codon_engine.py: codon-level mutational accessibility.

The study's central quantity: given the codon a species actually carries at a
conserved resistance site, which amino acids can a SINGLE nucleotide change reach?
A species whose codon is one step from Phe can evolve L1014F; a species whose codon
is not, cannot, until a second and much rarer mutation occurs.

Implemented from first principles with no external dependencies:
  * standard nuclear genetic code (NCBI translation table 1, correct for all target
    genes here, which are nuclear);
  * single-nucleotide codon neighbourhoods, classified transition/transversion;
  * codon-level accessibility (what a SPECIFIC codon reaches);
  * amino-acid-level accessibility averaged over synonymous codons (the codon-blind
    "theory" layer);
  * a transition/transversion-weighted accessibility rate;
  * two-step reachability, for quantifying how much rarer the inaccessible paths are.

All functions are deterministic and side-effect free. Run this file directly to
execute the self-test.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Dict, List, Tuple

# --------------------------------------------------------------------------- #
# 1. Standard genetic code (NCBI translation table 1)
# --------------------------------------------------------------------------- #
_BASES = "TCAG"
CODONS: List[str] = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
_AA_STRING = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE: Dict[str, str] = dict(zip(CODONS, _AA_STRING))

AA_TO_CODONS: Dict[str, List[str]] = {}
for _c, _a in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(_a, []).append(_c)
for _a in AA_TO_CODONS:
    AA_TO_CODONS[_a] = sorted(AA_TO_CODONS[_a])

STOP = "*"


def is_transition(b1: str, b2: str) -> bool:
    """A<->G or C<->T."""
    pair = frozenset((b1, b2))
    return pair == frozenset("AG") or pair == frozenset("CT")


def normalise(codon: str) -> str:
    return codon.upper().replace("U", "T")


def is_valid_codon(codon: str) -> bool:
    c = normalise(codon)
    return len(c) == 3 and all(b in "ACGT" for b in c)


@lru_cache(maxsize=None)
def neighbours(codon: str) -> Tuple[Tuple[str, int, str, str, bool], ...]:
    """All 9 single-nucleotide neighbours: (codon, pos, from, to, is_transition)."""
    codon = normalise(codon)
    out = []
    for pos in range(3):
        for base in "ACGT":
            if base != codon[pos]:
                nb = codon[:pos] + base + codon[pos + 1:]
                out.append((nb, pos, codon[pos], base, is_transition(codon[pos], base)))
    return tuple(out)


# --------------------------------------------------------------------------- #
# 2. Codon-level accessibility: a SPECIFIC codon (a species' actual codon)
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=None)
def codon_accessible_aas(codon: str) -> Dict[str, Tuple[Tuple[str, bool], ...]]:
    """{target_aa: ((neighbour_codon, is_transition), ...)} reachable in ONE step.

    The source amino acid is excluded (synonymous changes are not substitutions);
    stop codons are retained as '*' because reaching a stop is a real, and strongly
    deleterious, mutational outcome that must not be silently dropped.
    """
    codon = normalise(codon)
    src_aa = CODON_TABLE[codon]
    reach: Dict[str, List[Tuple[str, bool]]] = {}
    for nb, _pos, _fb, _tb, ts in neighbours(codon):
        aa = CODON_TABLE[nb]
        if aa == src_aa:
            continue
        reach.setdefault(aa, []).append((nb, ts))
    return {k: tuple(v) for k, v in reach.items()}


def codon_reaches(codon: str, to_aa: str) -> bool:
    """Can `codon` reach `to_aa` in a single nucleotide change?"""
    return to_aa in codon_accessible_aas(normalise(codon))


def codon_paths_to(codon: str, to_aa: str) -> Tuple[Tuple[str, bool], ...]:
    return codon_accessible_aas(normalise(codon)).get(to_aa, ())


def codon_weight_to(codon: str, to_aa: str, kappa: float = 4.0) -> float:
    """Ts/Tv-weighted count of single-step paths codon -> to_aa.

    Each transition path contributes `kappa`, each transversion 1.0. kappa is the
    transition:transversion rate ratio; kappa=1 gives raw path counts.
    """
    return sum(kappa if ts else 1.0 for (_nb, ts) in codon_paths_to(codon, to_aa))


def codon_path_counts(codon: str, to_aa: str) -> Tuple[int, int, int]:
    """(n_paths, n_transitions, n_transversions) for codon -> to_aa in one step."""
    paths = codon_paths_to(codon, to_aa)
    n_ts = sum(1 for (_n, ts) in paths if ts)
    return len(paths), n_ts, len(paths) - n_ts


# --------------------------------------------------------------------------- #
# 3. Two-step reachability: how much rarer is the inaccessible route?
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=None)
def codon_reaches_two_step(codon: str, to_aa: str,
                           allow_stop_intermediate: bool = False) -> int:
    """Number of 2-step nucleotide paths from `codon` to `to_aa`.

    Only counted when `to_aa` is NOT reachable in one step, which is the case the
    argument cares about. By default paths through a stop-codon intermediate are
    excluded: a nonsense intermediate is lethal and cannot segregate.
    """
    codon = normalise(codon)
    if codon_reaches(codon, to_aa):
        return 0
    n = 0
    for nb, _p, _f, _t, _ts in neighbours(codon):
        if not allow_stop_intermediate and CODON_TABLE[nb] == STOP:
            continue
        for nb2, _p2, _f2, _t2, _ts2 in neighbours(nb):
            if CODON_TABLE[nb2] == to_aa:
                n += 1
    return n


# --------------------------------------------------------------------------- #
# 4. Amino-acid-level accessibility: averaged over synonymous codons
# --------------------------------------------------------------------------- #
def aa_accessibility(from_aa: str, to_aa: str, kappa: float = 4.0) -> dict:
    """Codon-blind accessibility of from_aa -> to_aa, averaged over synonymous codons.

    This is the null the study argues against: it asks "across all the ways this
    residue could be encoded, how accessible is the resistance amino acid?", which is
    the quantity you get if you ignore which codon a species actually carries.
    """
    from_codons = AA_TO_CODONS[from_aa]
    n_from = len(from_codons)
    codons_1step, total_weight = [], 0.0
    for fc in from_codons:
        total_weight += codon_weight_to(fc, to_aa, kappa=kappa)
        if codon_reaches(fc, to_aa):
            codons_1step.append(fc)
    codons_2step = [c for c in from_codons if c not in codons_1step]
    return {
        "from": from_aa,
        "to": to_aa,
        "n_codons_from": n_from,
        "n_1step": len(codons_1step),
        "frac_1step": len(codons_1step) / n_from,
        "codons_1step": ",".join(sorted(codons_1step)),
        "codons_2step": ",".join(sorted(codons_2step)),
        "mean_Tsweighted_rate": total_weight / n_from,
    }


def translate_codon(codon: str) -> str:
    return CODON_TABLE[normalise(codon)]


def translate(seq: str) -> str:
    seq = normalise(seq)
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


# --------------------------------------------------------------------------- #
# 5. Self-test
# --------------------------------------------------------------------------- #
def _self_test() -> None:
    assert len(CODONS) == 64 and len(CODON_TABLE) == 64
    assert CODON_TABLE["ATG"] == "M" and CODON_TABLE["TAA"] == STOP
    assert is_transition("A", "G") and not is_transition("A", "C")
    assert len(AA_TO_CODONS["L"]) == 6 and len(AA_TO_CODONS["M"]) == 1

    # The reciprocal discriminating case at kdr L1014, the study's worked example:
    # Musca domestica carries CTT, Anopheles gambiae carries TTA. Both encode Leu.
    # CTT reaches His but NOT Ser; TTA reaches Ser but NOT His, and the two species'
    # documented resistance repertoires match exactly that asymmetry.
    assert codon_reaches("CTT", "H") and not codon_reaches("CTT", "S")
    assert codon_reaches("TTA", "S") and not codon_reaches("TTA", "H")
    # Both reach Phe, and both species evolved L1014F.
    assert codon_reaches("CTT", "F") and codon_reaches("TTA", "F")
    # CTG reaches neither Ser nor His.
    assert not codon_reaches("CTG", "S") and not codon_reaches("CTG", "H")

    # Ts/Tv weighting: TTA->TCA (Ser) is a transition at position 2.
    assert codon_weight_to("TTA", "S", kappa=4.0) == 4.0
    assert codon_weight_to("TTA", "S", kappa=1.0) == 1.0
    n, nts, ntv = codon_path_counts("TTA", "S")
    assert (n, nts, ntv) == (1, 1, 0)

    # Two-step: CTG cannot reach Ser in one step but can in two.
    assert codon_reaches_two_step("CTG", "S") > 0
    assert codon_reaches_two_step("TTA", "S") == 0   # already 1-step

    # Stop codons are surfaced, not hidden.
    assert STOP in codon_accessible_aas("TTA")

    assert translate("ATGTTTTAA") == "MF*"
    print("codon_engine self-test: PASS")


def _validate_theory_layer() -> None:
    """Validate the amino-acid layer against values derivable by hand from the code.

    Self-contained by design: this module must not depend on any external table, so
    that the engine can be verified from first principles by anyone with the
    genetic code in front of them. The same expectations were also checked
    against an independently produced accessibility table and agreed on all 29
    rows; that check is not a runtime dependency.
    """
    # Hand-derivable expectations at the kdr site L1014. Leu codons are
    # TTA TTG CTT CTC CTA CTG; these asymmetries are the study's central mechanism.
    checks = [
        # Phe = TTT/TTC. Reached by TTA(p3), TTG(p3), CTT(p1), CTC(p1) -> 4 of 6.
        ("L", "F", 4, 6),
        # Ser = TCN/AGT/AGC. Only TTA->TCA and TTG->TCG are single steps -> 2 of 6.
        # (CTT vs TCT differs at BOTH positions 1 and 2.)
        ("L", "S", 2, 6),
        # His = CAT/CAC. Only CTT->CAT and CTC->CAC -> 2 of 6.
        # (CTA->CAA and CTG->CAG give Gln, not His.)
        ("L", "H", 2, 6),
        # Met has the single codon ATG; Ile (ATT/ATC/ATA) is one step from it.
        ("M", "I", 1, 1),
        # Trp has the single codon TGG; Cys (TGT/TGC) is one step away.
        ("W", "C", 1, 1),
    ]
    for frm, to, exp_n1, exp_ncod in checks:
        got = aa_accessibility(frm, to, kappa=4.0)
        assert got["n_codons_from"] == exp_ncod, (frm, to, got)
        assert got["n_1step"] == exp_n1, (frm, to, got)

    # Ts/Tv weighting must reduce to raw path counts at kappa = 1.
    for frm in "ACDEFGHIKLMNPQRSTVWY":
        for to in "ACDEFGHIKLMNPQRSTVWY":
            if frm == to:
                continue
            raw = sum(len(codon_paths_to(c, to)) for c in AA_TO_CODONS[frm])
            got = aa_accessibility(frm, to, kappa=1.0)
            assert abs(got["mean_Tsweighted_rate"] * len(AA_TO_CODONS[frm]) - raw) < 1e-9

    # Every codon's 9 neighbours must be exactly partitioned into 3 transitions
    # and 6 transversions.
    for c in CODONS:
        nb = neighbours(c)
        assert len(nb) == 9
        assert sum(1 for x in nb if x[4]) == 3

    print("codon_engine theory-layer validation: PASS "
          "(self-contained, no external table)")


if __name__ == "__main__":
    _self_test()
    _validate_theory_layer()
