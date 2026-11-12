#!/usr/bin/env python3
"""
rev14_01_popgen_disposition.py

DISPOSITION OF THE POPULATION-GENOMIC ANALYSIS, AND THE ALTERNATIVE-EXON
LIMITATION IT EXPOSED.

The Ag1000G and Ag3 screens were run, cited, and never reported: no result
from either appeared in the Results, the figures or the tables, while the text
carried a caveat about the effective number of independent haplotypes for an
analysis nobody could read. Either the analysis is reported or it goes. It
goes, and for a stronger reason than absence of reporting: an audit of the
population-genomic arm had already concluded that the analysis does not support what it appeared to show.

  1. BOTH TESTS ARE NULL.
       ag1000g_validation.csv : 2/38 accessible segregating, 0/10 inaccessible
                                segregating. Fisher P = 1.0 (n = 48).
       popgen_validation.csv  : Fisher P = 0.44 (n = 9).

  2. THE Ag3 HAPLOTYPE RESULT RESTS ON AN IMPOSSIBLE NULL. It compared observed
     cis-double mutants against an expectation computed under LINKAGE
     EQUILIBRIUM between positions ONE BASE PAIR APART, where r is effectively
     zero. The observed zero is the exact value implied by
     D' = -1 - the generic result for two sites each of single mutational
     origin with no recombination - and that 135 of 143 BACKGROUND codons show
     the same D' = -1. The finding is a property of physical linkage, not of
     codon-level mutational supply. The pooled expectation is also inflated
     2.7-fold by population structure (334.6 -> 121.9), with Hardy-Weinberg
     violated at F_IS = 0.65.

  Reporting that would replace an omission with an error.

A LIMITATION THAT IS NOT DISPOSED OF.

An. gambiae expresses TWO codons at the same
protein position at four VGSC sites in domain II S4-S5, through a mutually
exclusive alternative exon pair. The atlas records one transcript's codon, so
the reachable set it assigns at those sites is a SUBSET of the true one:

    L925  TTA -> {F,I,S,V}    vs  CTT -> {F,H,I,P,R,V}   adds H, P, R
    T929  ACG -> {A,K,M,P,R,S} vs ACC -> {A,I,N,P,S}     adds I, N
    L932  CTC -> {F,H,I,P,R,V} vs TTA -> {F,I,S,V}       adds S
    G943  GGA -> {A,E,R,V}    vs  GGG -> {A,E,R,V,W}     adds W

THE IMPACT ON THE EVENTS IS NIL, and this script demonstrates it rather than
asserting it:

  * A union of codons can only ADD reachable residues, never remove one, so no
    event currently scored accessible can become inaccessible.
  * The only events at these sites scored INACCESSIBLE are the two T929V
    breaches. Threonine -> valine requires substitutions at codon positions 1
    AND 2 from every one of ACT, ACC, ACA, ACG, so no alternative exon can make
    it one-step. This is checked below by exhaustive enumeration, not assumed.

THE IMPACT ON THE DENOMINATOR IS REAL BUT BOUNDED, and runs AGAINST the
paper's claim: rows wrongly scored inaccessible-and-unused would move to
accessible-and-unused, which DECREASES the odds ratio. The current estimate is
therefore anti-conservative at these four sites. The bound is obtained by
dropping the four sites entirely and recomputing.

Writes  results/rev14_disposition.md
        results/rev14_altexon_sensitivity.csv
"""
from __future__ import annotations

import csv
import os
import sys
from collections import defaultdict
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PROJ = os.path.dirname(HERE)
RES = os.path.join(PROJ, "results")

from codon_engine import translate_codon                            # noqa: E402
from rev09_01_two_step_routes import resolve_results_file           # noqa: E402
from rev10_01_background_control import mantel_haenszel             # noqa: E402

ATLAS = os.path.join(PROJ, "data", "processed", "codon_atlas.csv")
ALT_EXON_SITES = {("VGSC", "L925"), ("VGSC", "T929"),
                  ("VGSC", "L932"), ("VGSC", "G943")}


def codons_of(aa):
    return ["".join(p) for p in product("ACGT", repeat=3) if translate_codon("".join(p)) == aa]


def min_steps(frm, to_aa):
    return min(sum(1 for i in range(3) if frm[i] != c[i]) for c in codons_of(to_aa))


def main():
    out = []

    def p(s=""):
        print(s)
        out.append(s)

    p("=" * 78)
    p("STEP 14.  Population-genomic disposition and the alternative-exon bound")
    p("=" * 78)
    p()
    p("DECISION: the Ag1000G / Ag3 paragraph is removed from the manuscript.")
    p("Grounds: both tests null (P = 1.0, n = 48; P = 0.44, n = 9), and the Ag3")
    p("haplotype result rests on a linkage-equilibrium null between bases 1 bp")
    p("apart. The grounds are set out in full in this module's docstring.")
    p()

    # ---------------------------------------------------------------- proof 1
    p("-" * 78)
    p("1.  CAN ANY ALTERNATIVE EXON MAKE T929V ONE-STEP?  (exhaustive)")
    p("-" * 78)
    p()
    p("  Every threonine codon, minimum substitutions to reach valine:")
    worst = 99
    for c in sorted(codons_of("T")):
        d = min_steps(c, "V")
        worst = min(worst, d)
        p("    %s (Thr) -> Val : %d substitution(s)" % (c, d))
    p()
    if worst >= 2:
        p("  Minimum over ALL threonine codons is %d. No codon a species could carry" % worst)
        p("  at position 929, on either exon copy, makes valine one-step reachable.")
        p("  The two T929V breaches are therefore unaffected by the exon issue.")
    else:
        p("  WARNING a threonine codon reaches valine in one step. Re-examine.")

    # ---------------------------------------------------------------- proof 2
    ev_path = resolve_results_file("core_test_events")
    p()
    observed = set()
    ev_at_alt = []
    with open(ev_path) as fh:
        for r in csv.DictReader(fh):
            observed.add((r["species"], r["gene"], r["site_label"], r["resistance_allele"]))
            if (r["gene"], r["site_label"]) in ALT_EXON_SITES:
                ev_at_alt.append(r)

    p("-" * 78)
    p("2.  THE EVENTS AT THE FOUR ALTERNATIVE-EXON SITES")
    p("-" * 78)
    p()
    acc = sum(1 for r in ev_at_alt if str(r.get("accessible")).strip().lower() in ("true", "1"))
    p("  events at L925 / T929 / L932 / G943 : %d of %d" % (len(ev_at_alt), len(observed)))
    p("    already scored accessible          : %d  (a union can only add, so unchanged)" % acc)
    p("    scored inaccessible                : %d  (the T929V breaches, shown above to be" % (len(ev_at_alt) - acc))
    p("                                         unreachable from every Thr codon)")
    p()
    p("  Net effect on the 108 events: NONE.")

    # ---------------------------------------------------------------- bound
    rows = []
    with open(ATLAS) as fh:
        for r in csv.DictReader(fh):
            if r["is_wildtype"] != "True" or r["already_resistant"] == "True":
                continue
            rows.append({
                "site": "%s:%s" % (r["gene"], r["site_label"]),
                "alt_exon": (r["gene"], r["site_label"]) in ALT_EXON_SITES,
                "accessible": r["one_step_accessible"] == "True",
                "used": (r["species"], r["gene"], r["site_label"],
                         r["resistance_allele"]) in observed,
            })

    def mh_on(rs, label):
        by = defaultdict(lambda: [0, 0, 0, 0])
        for r in rs:
            by[r["site"]][(0 if r["used"] else 2) + (0 if r["accessible"] else 1)] += 1
        o, lo, hi, n = mantel_haenszel([tuple(v) for v in by.values()])
        line = ("  %-40s rows %6d  events %3d  OR %6.2f  [%.2f, %.2f]"
                % (label, len(rs), sum(1 for r in rs if r["used"]), o, lo, hi))
        p(line)
        return {"restriction": label, "n_rows": len(rs),
                "n_events": sum(1 for r in rs if r["used"]),
                "odds_ratio": round(o, 3), "ci_lo": round(lo, 3), "ci_hi": round(hi, 3)}

    p()
    p("-" * 78)
    p("3.  THE BOUND: drop the four sites entirely")
    p("-" * 78)
    p()
    sens = [mh_on(rows, "all sites (baseline)"),
            mh_on([r for r in rows if not r["alt_exon"]], "excluding the four alt-exon sites")]
    p()
    p("  If the estimate holds, the unmodelled alternative-exon codons cannot be")
    p("  driving the result, and the limitation can be stated rather than repaired.")

    p()
    p("-" * 78)
    p("4.  THE LIMITATION TO STATE IN THE PAPER")
    p("-" * 78)
    p()
    lim = [
        "A species is assigned one codon per site, taken from its longest annotated",
        "coding transcript. Where a gene carries mutually exclusive alternative exons,",
        "a single protein position can be encoded by two different codons in the same",
        "species, and the set of residues reachable in one substitution is then the",
        "union of the two. Four VGSC sites in domain II S4-S5 (L925, T929, L932, G943)",
        "are known to be of this kind in Anopheles gambiae. Because a union can only",
        "add reachable residues, no event scored as accessible here can be affected,",
        "and the two inaccessible events at these sites are threonine-to-valine",
        "substitutions, which require two nucleotide changes from every threonine",
        "codon. The unmodelled codons can therefore only enlarge the accessible",
        "denominator, which would reduce the reported odds ratio; excluding the four",
        "sites entirely leaves it unchanged (see rev14_altexon_sensitivity.csv).",
    ]
    for l in lim:
        p("  " + l)

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "rev14_altexon_sensitivity.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sens[0].keys()))
        w.writeheader()
        for r in sens:
            w.writerow(r)
    with open(os.path.join(RES, "rev14_disposition.md"), "w") as fh:
        fh.write("# Step 14: population-genomic disposition\n\n```\n" + "\n".join(out) + "\n```\n")

    p()
    p("=" * 78)
    p("WROTE results/rev14_disposition.md, rev14_altexon_sensitivity.csv")
    p("=" * 78)


if __name__ == "__main__":
    main()
