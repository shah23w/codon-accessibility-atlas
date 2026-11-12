#!/usr/bin/env python3
"""
rev11_01_phylogenetic_control.py

STEP 11.  IS THE SPECIES-LEVEL ODDS RATIO A PHYLOGENETIC ARTEFACT?

Step 10 established the number the paper should now lead with: a SPECIES-LEVEL
odds ratio of 11.66 [4.28, 31.79] for resistance, against 2.01 [1.55, 2.60]
for natural substitution at the same sites.

Step 8b then found something that puts that number at risk. Accessibility is
not an independent draw per species - it tracks clade. In alphaE7, flies encode
72 to 92 percent of their glycines on codons that can reach aspartate in one
step; mosquitoes encode 44 to 56 percent. Codon usage is phylogenetically
clustered, so accessibility is too.

Documented resistance is also phylogenetically clustered: agricultural and
medical pests are not a random sample of arthropods. If both the predictor and
the outcome track clade, a model that treats 592 species as independent draws
will overstate the association.

WHAT THIS DOES.  Two things, neither of which needs R.

  1. FINER STRATIFICATION. The Mantel-Haenszel estimator conditions on strata
     and compares only within them. Stratifying by site alone allows a fly to
     be compared with a mosquito. Stratifying by ORDER x SITE, and then by
     FAMILY x SITE, restricts every comparison to close relatives. If the odds
     ratio survives that, between-clade structure is not driving it. If it
     collapses, it was.

  2. A WITHIN-CLADE PERMUTATION NULL. Accessibility labels are permuted among
     species WITHIN each order x site cell, which destroys the link between a
     species' own codon and what that species evolved while preserving both
     the clade composition and each clade's accessibility rate. The observed
     odds ratio is then read against that null. This is exact in the sense
     that matters: it assumes no tree, no branch lengths, and no model.

  Point 2 is deliberately preferred over a phylogenetic GLM. The project's
  step7_01_phyloglm.py fits Ives-Garland logistic regression per site-allele on
  an NCBI taxonomy topology with Grafen branch lengths. That is a reasonable
  check, but the branch lengths are nominal, it covers only the 32 studied
  species, and it pools by taking a median odds ratio with no interval. A
  permutation that conditions on clade makes weaker assumptions and yields a
  calibrated P-value. The GLM remains available as a confirmatory analysis.

  DIAGNOSTIC. Before either test, the accessibility rate is reported per order,
  so the phylogenetic signal in the PREDICTOR is visible rather than asserted.

Reads   data/processed/codon_atlas.csv
        data/interim/taxonomy_lineages.json
        results/core_test_events*.csv
Writes  results/rev11_stratification.csv
        results/rev11_accessibility_by_clade.csv
        results/rev11_summary.md
"""
from __future__ import annotations

import csv
import json
import math
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PROJ = os.path.dirname(HERE)
RES = os.path.join(PROJ, "results")

from rev09_01_two_step_routes import resolve_results_file           # noqa: E402
from rev10_01_background_control import mantel_haenszel             # noqa: E402

ATLAS = os.path.join(PROJ, "data", "processed", "codon_atlas.csv")
LINEAGES = os.path.join(PROJ, "data", "interim", "taxonomy_lineages.json")
PERMS = 5000
SEED = 20260927


def family_of(lineage: str) -> str:
    """
    Arthropod family names end in -idae. Taking the last such token in the
    NCBI lineage is a reliable family assignment for this panel and avoids
    depending on rank annotations that are missing for many taxa.
    """
    fams = [p.strip() for p in (lineage or "").split(";") if p.strip().endswith("idae")]
    return fams[-1] if fams else ""


def strata_or(rows, keyfn):
    by = defaultdict(lambda: [0, 0, 0, 0])
    for r in rows:
        by[keyfn(r)][(0 if r["used"] else 2) + (0 if r["accessible"] else 1)] += 1
    return mantel_haenszel([tuple(v) for v in by.values()]), by


def main():
    print("=" * 78)
    print("STEP 11.  Is the species-level odds ratio a phylogenetic artefact?")
    print("=" * 78)
    print()

    ev_path = resolve_results_file("core_test_events")
    print()
    observed = set()
    with open(ev_path) as fh:
        for r in csv.DictReader(fh):
            observed.add((r["species"], r["gene"], r["site_label"], r["resistance_allele"]))

    with open(LINEAGES) as fh:
        lin = json.load(fh)
    fam_by_taxid = {k: family_of(v.get("lineage", "")) for k, v in lin.items()}

    rows = []
    with open(ATLAS) as fh:
        for r in csv.DictReader(fh):
            if r["is_wildtype"] != "True" or r["already_resistant"] == "True":
                continue
            fam = fam_by_taxid.get(str(r.get("tax_id", "")), "") or "UNPLACED"
            rows.append({
                "species": r["species"], "order": r.get("order") or "UNPLACED",
                "family": fam, "gene": r["gene"], "site_label": r["site_label"],
                "site": "%s:%s" % (r["gene"], r["site_label"]),
                "accessible": r["one_step_accessible"] == "True",
                "used": (r["species"], r["gene"], r["site_label"],
                         r["resistance_allele"]) in observed,
            })

    out = []

    def p(s=""):
        print(s)
        out.append(s)

    # ------------------------------------------------------- diagnostic
    p("-" * 78)
    p("DIAGNOSTIC.  IS ACCESSIBILITY ITSELF CLUSTERED BY CLADE?")
    p("-" * 78)
    p()
    by_order = defaultdict(lambda: [0, 0, set()])
    for r in rows:
        o = by_order[r["order"]]
        o[0] += 1 if r["accessible"] else 0
        o[1] += 1
        o[2].add(r["species"])
    p("  %-22s %8s %8s %9s" % ("order", "species", "rows", "P(access)"))
    clade_rows = []
    for o, (k, n, sp) in sorted(by_order.items(), key=lambda kv: -kv[1][1]):
        if n < 200:
            continue
        p("  %-22s %8d %8d %9.3f" % (o, len(sp), n, k / n))
        clade_rows.append({"order": o, "n_species": len(sp), "n_rows": n,
                           "p_accessible": round(k / n, 4)})
    vals = [c["p_accessible"] for c in clade_rows]
    if len(vals) > 1:
        p()
        p("  spread across orders: %.3f to %.3f" % (min(vals), max(vals)))
        p("  If these differ materially, accessibility carries clade signal and the")
        p("  unstratified estimate is at risk. That is what the tests below settle.")

    # ------------------------------------------------------- stratification
    p()
    p("-" * 78)
    p("1.  DOES THE ODDS RATIO SURVIVE FINER STRATIFICATION?")
    p("-" * 78)
    p()
    p("  Every comparison happens WITHIN a stratum, so finer strata mean")
    p("  comparisons between closer relatives.")
    p()
    schemes = [
        ("site", lambda r: r["site"]),
        ("order x site", lambda r: (r["order"], r["site"])),
        ("family x site", lambda r: (r["family"], r["site"])),
    ]
    strat_rows = []
    for name, fn in schemes:
        (or_, lo, hi, n_inf), by = strata_or(rows, fn)
        p("  %-16s strata %6d  informative %5d   OR %7s  95%% CI [%s, %s]"
          % (name, len(by), n_inf,
             "inf" if or_ == float("inf") else "%.2f" % or_,
             "%.2f" % lo if lo == lo else "na",
             "%.2f" % hi if hi == hi else "na"))
        strat_rows.append({"stratification": name, "n_strata": len(by),
                           "n_informative": n_inf,
                           "odds_ratio": "inf" if or_ == float("inf") else round(or_, 3),
                           "ci_lo": round(lo, 3) if lo == lo else "",
                           "ci_hi": round(hi, 3) if hi == hi else ""})
    p()
    p("  A stable odds ratio across these three rows means the association is")
    p("  present WITHIN clades and is not produced by comparing distant ones.")
    p("  Note that finer strata discard more data - a stratum with no variation")
    p("  contributes nothing - so the interval widens even when the estimate holds.")

    # ------------------------------------------------------- permutation
    p()
    p("-" * 78)
    p("2.  WITHIN-CLADE PERMUTATION  (%d permutations, seed %d)" % (PERMS, SEED))
    p("-" * 78)
    p()
    p("  Accessibility labels are shuffled among species within each order x site")
    p("  cell. Clade composition and each clade's accessibility rate are preserved;")
    p("  only the link between a species' own codon and what it evolved is broken.")
    p()
    cells = defaultdict(list)
    for r in rows:
        cells[(r["order"], r["site"])].append(r)
    obs_or = strata_or(rows, lambda r: r["site"])[0][0]

    # Shuffling labels row by row would be exact but slow. Within one
    # order x site cell of n rows, k of them accessible and u of them used,
    # permuting the accessibility labels makes the count of used-AND-accessible
    # rows exactly Hypergeometric(n, k, u). Drawing from that distribution is
    # the SAME null, sampled directly, and costs one draw per cell instead of
    # one shuffle per row. The other three cells follow by subtraction.
    try:
        import numpy as np
        rng_np = np.random.default_rng(SEED)
        draw = rng_np.hypergeometric
    except ImportError:                       # pure-python fallback
        rng_py = random.Random(SEED)

        def draw(ngood, nbad, nsample, size=None):
            pool = [1] * ngood + [0] * nbad
            rng_py.shuffle(pool)
            return sum(pool[:nsample])

    cell_spec = []
    for (order, site), cell in cells.items():
        n = len(cell)
        k = sum(1 for r in cell if r["accessible"])
        u = sum(1 for r in cell if r["used"])
        cell_spec.append((site, n, k, u))

    null = []
    for _ in range(PERMS):
        by = defaultdict(lambda: [0, 0, 0, 0])
        for site, n, k, u in cell_spec:
            if u == 0 or k == 0 or k == n:
                a = min(u, k) if k == n else 0
            else:
                a = int(draw(k, n - k, u))
            t = by[site]
            t[0] += a            # used & accessible
            t[1] += u - a        # used & not accessible
            t[2] += k - a        # not used & accessible
            t[3] += n - u - k + a
        o = mantel_haenszel([tuple(v) for v in by.values()])[0]
        if o == o and o != float("inf"):
            null.append(o)
    null.sort()
    ge = sum(1 for x in null if x >= obs_or)
    pval = (ge + 1) / (len(null) + 1)
    med = null[len(null) // 2] if null else float("nan")
    q95 = null[int(0.95 * len(null))] if null else float("nan")
    p("  observed OR (site-stratified)        : %.2f" % obs_or)
    p("  null median                          : %.2f" % med)
    p("  null 95th percentile                 : %.2f" % q95)
    p("  permutations >= observed             : %d of %d" % (ge, len(null)))
    p("  P                                    : %.4g" % pval)
    p()
    if pval < 0.05:
        p("  The observed odds ratio lies outside what clade structure alone produces.")
        p("  Phylogenetic non-independence does not account for the association.")
    else:
        p("  The observed odds ratio is WITHIN the range clade structure alone")
        p("  produces. On this evidence the species-level association cannot be")
        p("  separated from phylogenetic non-independence, and the paper must say so.")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "rev11_stratification.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(strat_rows[0].keys()))
        w.writeheader()
        for r in strat_rows:
            w.writerow(r)
    if clade_rows:
        with open(os.path.join(RES, "rev11_accessibility_by_clade.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(clade_rows[0].keys()))
            w.writeheader()
            for r in clade_rows:
                w.writerow(r)
    with open(os.path.join(RES, "rev11_summary.md"), "w") as fh:
        fh.write("# Step 11: phylogenetic control\n\n```\n")
        fh.write("\n".join(out))
        fh.write("\n```\n")
        fh.write("\npermutation P = %.4g (%d permutations, seed %d)\n" % (pval, PERMS, SEED))

    p()
    p("=" * 78)
    p("WROTE results/rev11_stratification.csv, rev11_accessibility_by_clade.csv,")
    p("      rev11_summary.md")
    p("=" * 78)


if __name__ == "__main__":
    main()
