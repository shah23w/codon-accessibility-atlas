#!/usr/bin/env python3
"""
enum01  Where in the genetic code can codon accessibility ever matter?

This module uses no data. It is arithmetic on the standard genetic code, and it
produces the paper's generality result: a statement, derivable before any
species is looked at, about which amino acid substitutions CAN show a codon
accessibility effect in any organism, and which cannot.

THE DEFINITION. A site whose wild-type residue is X gives an informative test
for a target residue Y when the synonymous codons of X disagree about whether
they reach Y in one nucleotide change: some do, some do not. If every codon of X
reaches Y, every species can reach it and accessibility predicts nothing. If no
codon of X reaches Y, no species can and accessibility again predicts nothing.
Only disagreement among synonymous codons makes a species' own codon
informative.

    informative(X) = { Y : 0 < |{c in codons(X) : Y in reach(c)}| < |codons(X)| }

WHY THE RULE IS NOT DEGENERACY. "Every two-fold degenerate residue contributes
exactly zero" is false. Lysine is two-fold degenerate (AAA, AAG) yet
informative, because a second-position change gives ATA (Ile) from AAA but ATG
(Met) from AAG -- the ATN block is itself split between two amino acids.
Glutamate and glutamine are two-fold degenerate over the same purines and are
NOT informative, because their neighbour blocks GTN/GCN/GGN and CTN/CCN/CGN are
unsplit. The rule is checked here rather than asserted:

    informativeness requires EITHER the synonymous family of X to be split
    across codon blocks, OR a one-step neighbour block to be split.

The module then scores every site in this project's own reference table, and the
cross-kingdom sites from the literature, against the enumeration.

writes
    results/enum01_informative_targets.csv     per amino acid
    results/enum01_site_predictions.csv        per site, ours and literature
    logs/enum01_informativeness.log
"""
import csv, os, sys, time
from itertools import product

T0 = time.time()
DIR_RES = "results"
DIR_LOG = "logs"
REF = os.path.join("data", "refs", "reference_sites.csv")

os.makedirs(DIR_LOG, exist_ok=True)
_lines = []
def log(s=""):
    print(("[%6.1fs] %s" % (time.time() - T0, s)) if s else "")
    _lines.append(s)

FAILURES = []
def chk(label, condition, detail=""):
    if condition:
        log("  [PASS] %s%s" % (label, ("   " + detail) if detail else ""))
    else:
        log("  [FAIL] %s%s" % (label, ("   " + detail) if detail else ""))
        FAILURES.append(label)
    return condition

# --------------------------------------------------------------------------- #
#  the standard genetic code, written out rather than imported, so that this
#  module has no dependency and a reader can check every line of it
# --------------------------------------------------------------------------- #
BASES = "TCAG"
AAS = (
    "FFLLSSSSYY**CC*W"
    "LLLLCCCCHHQQRRRR".replace("CCCC", "PPPP")
    + "IIIMTTTTNNKKSSRR"
    + "VVVVAAAADDEEGGGG"
)
CODON_TABLE = {}
i = 0
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            CODON_TABLE[b1 + b2 + b3] = AAS[i]
            i += 1

CODONS_OF = {}
for c, a in CODON_TABLE.items():
    CODONS_OF.setdefault(a, []).append(c)

AA20 = sorted(a for a in CODONS_OF if a != "*")


def neighbours(codon):
    """every codon one nucleotide substitution away"""
    out = []
    for pos in range(3):
        for b in BASES:
            if b != codon[pos]:
                out.append(codon[:pos] + b + codon[pos + 1:])
    return out


def reach(codon):
    """amino acids reachable from this codon in one change, stops excluded"""
    return {CODON_TABLE[n] for n in neighbours(codon) if CODON_TABLE[n] != "*"}


def informative_targets(x):
    """targets Y for which the synonymous codons of X disagree about reachability"""
    cs = CODONS_OF[x]
    out = []
    for y in AA20:
        if y == x:
            continue
        hits = sum(1 for c in cs if y in reach(c))
        if 0 < hits < len(cs):
            out.append(y)
    return sorted(out)


log("=" * 74)
log("enum01  which substitutions can ever show a codon accessibility effect")
log("=" * 74)

# sanity on the table itself before anything is built on it
chk("the genetic code table has 64 codons", len(CODON_TABLE) == 64,
    "%d" % len(CODON_TABLE))
chk("it encodes 20 amino acids and a stop", len(CODONS_OF) == 21,
    "%d symbols" % len(CODONS_OF))
chk("ATG is the only methionine codon", CODONS_OF["M"] == ["ATG"],
    "%s" % CODONS_OF["M"])
chk("TGG is the only tryptophan codon", CODONS_OF["W"] == ["TGG"],
    "%s" % CODONS_OF["W"])
chk("leucine has six codons", len(CODONS_OF["L"]) == 6,
    "%d" % len(CODONS_OF["L"]))
chk("serine has six codons in two separated blocks", len(CODONS_OF["S"]) == 6,
    "%s" % " ".join(sorted(CODONS_OF["S"])))
chk("there are three stop codons", len(CODONS_OF["*"]) == 3,
    "%s" % " ".join(sorted(CODONS_OF["*"])))

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("1  THE ENUMERATION")
log("-" * 74)
log("   For each wild-type residue, the targets whose reachability DIFFERS")
log("   between its synonymous codons. Those are the only substitutions where")
log("   a species' own codon can predict anything.")
log()

TABLE = {x: informative_targets(x) for x in AA20}
ordered = sorted(AA20, key=lambda x: (-len(TABLE[x]), x))

log("   %-4s %-7s %-8s  %s" % ("res", "codons", "n_infor", "informative targets"))
for x in ordered:
    log("   %-4s %-7d %-8d  %s"
        % (x, len(CODONS_OF[x]), len(TABLE[x]),
           " ".join(TABLE[x]) if TABLE[x] else "-- none, dead by arithmetic --"))

live = [x for x in AA20 if TABLE[x]]
dead = [x for x in AA20 if not TABLE[x]]
log()
log("   informative residues (%d) : %s" % (len(live), " ".join(sorted(live))))
log("   dead residues       (%d) : %s" % (len(dead), " ".join(sorted(dead))))

EXPECTED = {"R": 11, "S": 11, "L": 9, "G": 5, "I": 5, "T": 5, "V": 5,
            "A": 2, "K": 2, "P": 2}
log()
log("   Against an independently derived table:")
ok = True
for x, n in sorted(EXPECTED.items()):
    got = len(TABLE[x])
    if got != n:
        ok = False
        log("     %s expected %d, computed %d   <-- MISMATCH" % (x, n, got))
chk("every residue's informative-target count matches the published enumeration",
    ok and len(live) == len(EXPECTED),
    "%d informative residues, expected %d" % (len(live), len(EXPECTED)))
chk("the other ten residues are dead by arithmetic",
    set(dead) == set(AA20) - set(EXPECTED),
    " ".join(sorted(dead)))

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("2  THE CORRECTION: 'TWO-FOLD DEGENERATE MEANS ZERO' IS FALSE")
log("-" * 74)

twofold = [x for x in AA20 if len(CODONS_OF[x]) == 2]
log("   two-fold degenerate residues : %s" % " ".join(sorted(twofold)))
log()
for x in sorted(twofold):
    log("   %s  %-9s -> %s"
        % (x, "/".join(sorted(CODONS_OF[x])),
           (" ".join(TABLE[x]) if TABLE[x] else "none")))

chk("lysine is two-fold degenerate yet informative",
    "K" in twofold and len(TABLE["K"]) > 0,
    "K reaches %s" % " ".join(TABLE["K"]))
chk("glutamate is two-fold degenerate and NOT informative",
    "E" in twofold and len(TABLE["E"]) == 0)
chk("glutamine is two-fold degenerate and NOT informative",
    "Q" in twofold and len(TABLE["Q"]) == 0)
chk("so the slogan 'every two-fold degenerate residue contributes zero' is refuted",
    len(TABLE["K"]) > 0 and len(TABLE["E"]) == 0,
    "K is a counterexample; E and Q are not")

log()
log("   Why K differs from E and Q, codon by codon:")
for x in ("K", "E", "Q"):
    for c in sorted(CODONS_OF[x]):
        mid = {c[:1] + b + c[2:]: CODON_TABLE[c[:1] + b + c[2:]]
               for b in BASES if b != c[1]}
        log("     %s %s  second-position changes -> %s"
            % (x, c, "  ".join("%s=%s" % (k, v) for k, v in sorted(mid.items()))))
log()
log("   K: AAA gives ATA (Ile) but AAG gives ATG (Met). The ATN block is split")
log("      between two amino acids, so the two lysine codons disagree.")
log("   E: GAA gives GTA and GAG gives GTG, both Val. The GTN block is unsplit,")
log("      so the two glutamate codons agree and nothing is informative.")
log()
log("   CORRECT RULE: informativeness requires EITHER the synonymous family to")
log("   be split across codon blocks, OR a one-step neighbour block to be split.")

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("3  THIS PROJECT'S OWN SITES, SCORED AGAINST THE ENUMERATION")
log("-" * 74)

SROWS = []
if not os.path.exists(REF):
    log("   FATAL: %s not found" % REF)
    FAILURES.append("reference_sites.csv missing")
else:
    seen = set()
    rows = list(csv.DictReader(open(REF)))
    log()
    log("   %-10s %-9s %-4s %-4s %-9s %s"
        % ("gene", "site", "wt", "res", "capable", "why"))
    for r in rows:
        key = (r["gene"], r["site_label"], r["resistance_allele"])
        if key in seen:
            continue
        seen.add(key)
        wt, res = r["wild_type_aa"], r["resistance_allele"]
        if wt not in TABLE:
            continue
        capable = res in TABLE[wt]
        if TABLE[wt]:
            why = ("target is informative for %s" % wt if capable
                   else "%s is informative only for %s" % (wt, " ".join(TABLE[wt])))
        else:
            why = "%s is dead by arithmetic" % wt
        log("   %-10s %-9s %-4s %-4s %-9s %s"
            % (r["gene"], r["site_label"], wt, res,
               "YES" if capable else "no", why))
        SROWS.append(dict(source="this study", system="arthropod",
                          gene=r["gene"], site=r["site_label"],
                          wild_type=wt, resistance_allele=res,
                          wt_informative_targets=" ".join(TABLE[wt]),
                          predicted_informative=int(capable)))

    n_cap = sum(d["predicted_informative"] for d in SROWS)
    log()
    log("   %d of %d curated site x allele pairs are predicted informative"
        % (n_cap, len(SROWS)))

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("4  THE CROSS-KINGDOM PREDICTION, MADE WITHOUT DATA")
log("-" * 74)
log("   These sites are typed from the literature, not computed from our atlas,")
log("   and are labelled as such. They are the sites that killed the three")
log("   cross-kingdom arms of the original project. The enumeration says in")
log("   advance which of them could ever have worked.")

# (system, gene, site, wild-type, documented resistance alleles, citation)
LIT = [
    ("fungal", "CYP51/erg11", "Y132", "Y", ["F", "H"],
     "Sterol 14a-demethylase; Y132F/H in Candida and plant pathogens"),
    ("fungal", "SDH B", "H277", "H", ["Y", "R"],
     "Succinate dehydrogenase B H277Y/R"),
    ("fungal", "beta-tubulin", "E198", "E", ["A", "K", "G"],
     "Benzimidazole resistance E198A/K/G"),
    ("fungal", "beta-tubulin", "F200", "F", ["Y"],
     "Benzimidazole resistance F200Y"),
    ("fungal", "cytb", "F129", "F", ["L"],
     "QoI resistance F129L"),
    ("fungal", "cytb", "G143", "G", ["A"],
     "QoI resistance G143A"),
    ("weed", "ALS/AHAS", "W574", "W", ["L"],
     "ALS inhibitor resistance W574L"),
    ("weed", "ALS/AHAS", "P197", "P", ["H", "Q", "S", "T", "L", "A"],
     "ALS inhibitor resistance; P197 is the most substituted ALS residue"),
    ("weed", "EPSPS", "P106", "P", ["S", "T", "A", "L"],
     "Glyphosate resistance P106S/T/A/L"),
]

log()
log("   %-8s %-14s %-7s %-4s %-18s %s"
    % ("system", "gene", "site", "wt", "alleles capable", "verdict"))
for sysname, gene, site, wt, alleles, note in LIT:
    cap = [a for a in alleles if a in TABLE.get(wt, [])]
    if not TABLE.get(wt):
        verdict = "DEAD: %s is dead by arithmetic" % wt
    elif cap:
        verdict = "CAN work via %s" % " ".join(cap)
    else:
        verdict = "dead here: %s informative only for %s" % (wt, " ".join(TABLE[wt]))
    log("   %-8s %-14s %-7s %-4s %-18s %s"
        % (sysname, gene, site, wt, " ".join(cap) if cap else "none", verdict))
    for a in alleles:
        SROWS.append(dict(source="literature", system=sysname, gene=gene,
                          site=site, wild_type=wt, resistance_allele=a,
                          wt_informative_targets=" ".join(TABLE.get(wt, [])),
                          predicted_informative=int(a in TABLE.get(wt, [])),
                          note=note))

log()
log("   Where the surveyed alleles stand.")
log("     Of the fungal (E198, F200, F129, Y132) and weed (W574, P197)")
log("     alleles, the fungal four and W574 cannot show the effect. P197 can.")
log()
p_inf = TABLE["P"]
log("     Proline is informative for %s." % " and ".join(p_inf))
for c in sorted(CODONS_OF["P"]):
    log("       %s -> %s" % (c, " ".join(sorted(reach(c) & set(p_inf))) or "neither"))
log("     P197H and P197Q are both documented ALS-resistance alleles, so P197 is")
log("     a site where the test COULD be run in weeds: W574 is dead and P197")
log("     is live.")

chk("proline is informative, so P197 admits the test",
    len(p_inf) == 2 and "H" in p_inf and "Q" in p_inf,
    "P informative for %s" % " ".join(p_inf))
chk("tryptophan is uninformative, so W574 cannot admit the test",
    len(TABLE["W"]) == 0)
chk("the four fungal residues E, F, Y and the G143 glycine behave as claimed",
    len(TABLE["E"]) == 0 and len(TABLE["F"]) == 0 and len(TABLE["Y"]) == 0,
    "E, F, Y all dead; G is live (%s) so G143A needs checking separately"
    % " ".join(TABLE["G"]))

g143 = "A" in TABLE["G"]
log()
log("     Note G143A: glycine IS informative (for %s), so whether G143A can"
    % " ".join(TABLE["G"]))
log("     show the effect depends on whether alanine is among those targets: %s."
    % ("YES, it can" if g143 else "no, it cannot"))

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("5  DOES THE A PRIORI PREDICTION SURVIVE CONTACT WITH THE DATA?")
log("-" * 74)
log("   Sections 1 to 4 use no data at all. This section is the test of them.")
log()
log("   The enumeration says a (wild-type, target) pair can only be informative")
log("   if the synonymous codons of the wild-type residue disagree about")
log("   reaching the target. If that is right, then the events our core test")
log("   found to be CODON-DISCRIMINATING -- where a species' own codon reaches")
log("   some but not all of the alleles curated at its site -- must sit at sites")
log("   carrying at least one predicted-informative pair. An event at a site")
log("   where every pair is dead by arithmetic cannot discriminate, and finding")
log("   one would refute the enumeration.")

EV = os.path.join(DIR_RES, "core_test_events_primary.csv")
if not os.path.exists(EV):
    log("   skipped: %s not found" % EV)
else:
    log()
    log("   THE RIGHT QUANTITY TO PREDICT. An event is supply-informative when")
    log("   some but not all panel species could reach that allele at that site,")
    log("   i.e. 0 < site_supply_fraction < 1. That is between-species variation,")
    log("   which is what the enumeration is about. It is NOT the same as the")
    log("   core test's codon-discriminating set, which asks whether one species'")
    log("   own codon separates the alleles in its site repertoire;")
    log("   canonical_quantities.csv warns that the two are different questions")
    log("   and are not in a subset relation.")
    log()
    log("   The hard direction, which sampling cannot excuse:")
    log("     if the synonymous codons of the wild-type residue AGREE about")
    log("     reaching the allele, every wild-type species behaves the same way,")
    log("     so the supply fraction MUST be exactly 0 or exactly 1.")
    log("   The soft direction: a pair the enumeration calls informative may")
    log("   still show a fraction of 0 or 1 if the panel happens to carry only")
    log("   one of the relevant codons. That is sampling, not a failed rule.")

    rows = list(csv.DictReader(open(EV)))
    tab = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
    violations = []
    for r in rows:
        try:
            f = float(r["site_supply_fraction"])
        except (KeyError, ValueError):
            continue
        wt, al = r["wild_type_aa"], r["resistance_allele"]
        pred = int(al in TABLE.get(wt, []))
        obs = int(0.0 < f < 1.0)
        tab[(pred, obs)] += 1
        if obs and not pred:
            violations.append((r["species"], r["gene"], r["site_label"],
                               wt, al, r.get("codon", "?"), f))

    log()
    log("   %-34s %12s %12s" % ("", "fraction 0/1", "strictly between"))
    log("   %-34s %12d %12d" % ("pair predicted NOT informative",
                                tab[(0, 0)], tab[(0, 1)]))
    log("   %-34s %12d %12d" % ("pair predicted informative",
                                tab[(1, 0)], tab[(1, 1)]))
    log()
    log("   events scored : %d" % sum(tab.values()))

    chk("no uninformative pair shows between-species variation (the hard direction)",
        len(violations) == 0,
        "%d violation(s)" % len(violations))
    if violations:
        log()
        log("   VIOLATIONS -- these refute the enumeration and must be explained:")
        for v in violations[:15]:
            log("     %-26s %-6s %-7s %s->%s codon %s fraction %.4f"
                % (v[0][:26], v[1], v[2], v[3], v[4], v[5], v[6]))

    n_inf = tab[(0, 1)] + tab[(1, 1)]
    # The canonical count is read from the register, not typed here: a check
    # written against a literal stops being a check the moment the pipeline
    # moves and the literal does not.
    _canon = {}
    with open(os.path.join(DIR_RES, "canonical_quantities.csv"),
              encoding="utf-8-sig") as _fh:
        for _r in csv.DictReader(_fh):
            _canon[_r["quantity"]] = _r["value"]
    n_canon = int(float(_canon["n_events_supply_informative"]))
    log()
    log("   supply-informative events found here : %d" % n_inf)
    log("   canonical n_events_supply_informative : %d" % n_canon)
    chk("this reproduces the canonical count of supply-informative events",
        n_inf == n_canon, "%d against %d" % (n_inf, n_canon))

    # Figure 5c reads this table rather than recovering the counts from the
    # prose above.
    with open(os.path.join(DIR_RES, "enum01_predicted_vs_observed.csv"), "w",
              newline="") as _fh:
        _w = csv.writer(_fh, lineterminator="\n")
        _w.writerow(["pair_predicted_informative", "observed_reachable_fraction",
                     "n_events"])
        for _pred in (0, 1):
            for _obs in (0, 1):
                _w.writerow([_pred,
                             "strictly between 0 and 1" if _obs else "exactly 0 or 1",
                             tab[(_pred, _obs)]])
        _w.writerow(["", "total", sum(tab.values())])
    log("  wrote results/enum01_predicted_vs_observed.csv")

log()
log("   THE STATEMENT FOR THE PAPER. It is not the residue that has to be")
log("   informative, it is the (wild-type, target) PAIR. Glycine is informative")
log("   for C, D, E, S and W, so G143A in fungal cytb is dead even though")
log("   glycine is live; proline is informative for H and Q, so EPSPS P106S is")
log("   dead while ALS P197H is live. Stating the rule at residue level gets")
log("   these wrong.")

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("WRITING")
log("-" * 74)
with open(os.path.join(DIR_RES, "enum01_informative_targets.csv"), "w",
          newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["wild_type_aa", "n_codons", "n_informative_targets",
                "informative_targets"])
    for x in ordered:
        w.writerow([x, len(CODONS_OF[x]), len(TABLE[x]), " ".join(TABLE[x])])
log("  wrote results/enum01_informative_targets.csv")

flds = ["source", "system", "gene", "site", "wild_type", "resistance_allele",
        "wt_informative_targets", "predicted_informative", "note"]
with open(os.path.join(DIR_RES, "enum01_site_predictions.csv"), "w",
          newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=flds, lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    for r in SROWS:
        w.writerow(r)
log("  wrote results/enum01_site_predictions.csv")

log()
log("=" * 74)
if FAILURES:
    log("RESULT: %d CHECK(S) FAILED" % len(FAILURES))
    for f in FAILURES:
        log("   - %s" % f)
else:
    log("RESULT: every check passed")
log("=" * 74)

with open(os.path.join(DIR_LOG, "enum01_informativeness.log"), "w") as fh:
    fh.write("\n".join(_lines) + "\n")

sys.exit(1 if FAILURES else 0)
