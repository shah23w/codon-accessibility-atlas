#!/usr/bin/env python3
"""
rev12_01_tolerability_control.py

STEP 12.  DOES ACCESSIBILITY SURVIVE ADJUSTMENT FOR TOLERABILITY,
          AT THE LEVEL THE PAPER NOW CLAIMS?

THE CONFOUND.  The genetic code is error-minimising: codons one nucleotide
apart tend to encode biochemically similar residues. So a substitution that is
mutationally accessible is also, on average, biochemically conservative - and
conservative substitutions are more likely to be tolerated by the protein.
Accessibility and tolerability are therefore correlated BY CONSTRUCTION, and an
association between accessibility and which allele is used could be an
association with tolerability wearing a disguise. step6_03_multivariate.py
states this plainly in its own docstring.

WHAT THE MANUSCRIPT DID.  It adjusted for ESM-2 variant-effect scores. That
analysis is sound as far as it goes, but it has three problems:

  1. UNDOCUMENTED. The specification of the scoring exists, in
     step6_01_esm_variant_effects.py and in the columns of
     esm_variant_effects.csv, but was not written out anywhere a reader could
     find it: how the model is applied, which checkpoint, which parameters,
     and what its limitations are. This script prints it.

  2. WRONG LEVEL. ESM was scored once per gene x site x candidate residue
     against a single REFERENCE protein. It has no species dimension. The
     headline estimate is a SPECIES-LEVEL odds ratio (Step 10), so a
     site-level covariate cannot formally adjust it.

  3. MODEL-DEPENDENT. A protein language model's score is not a measurement.
     Resting the tolerability control on one model, from one family, is a
     single point of failure.

WHAT THIS DOES.  Adjusts the species-level odds ratio for tolerability using
three covariates that are independent of one another, and reports all of them.

  A. GRANTHAM DISTANCE. Physicochemical distance between the wild-type and
     candidate residue, from composition, polarity and molecular volume
     (Grantham 1974 Science 185:862-864). Chosen deliberately over BLOSUM: a
     BLOSUM score is estimated FROM observed substitutions, which already
     embeds mutational accessibility, so conditioning on it would be partly
     circular. Grantham is purely physicochemical and carries no information
     about what evolution actually did. It is computed here from Grantham's
     own formula and property table and self-validated against five published
     pairwise distances before use.

  B. ESM-2 SCORE. The existing masked-marginal scores, joined in. Site-level,
     and labelled as such.

  C. PANEL RESIDUE FREQUENCY. How often the candidate residue is actually seen
     at that site across the 594-species panel. This is an empirical
     tolerability measure taken from the data itself: a residue that no
     arthropod carries at a position is unlikely to be tolerated there. Unlike
     ESM it needs no model, and unlike ESM it is available for every row.

  ADJUSTMENT IS BY STRATIFICATION, not regression. Mantel-Haenszel conditions
  on strata exactly and makes no functional-form assumption about how
  tolerability enters. Each covariate is cut into tertiles and crossed with
  site, so accessible and inaccessible alleles are compared only within the
  same site AND the same tolerability band.

  READING IT. If the odds ratio holds across all three adjustments, the
  association is not tolerability in disguise. If it collapses under any of
  them, that must be reported - and the covariate that kills it tells you which
  mechanism was really operating.

Reads   data/processed/codon_atlas.csv
        results/esm_variant_effects*.csv
        results/core_test_events*.csv
Writes  results/rev12_adjusted.csv
        results/rev12_rows.csv
        results/rev12_methods_esm.md
        results/rev12_summary.md
"""
from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PROJ = os.path.dirname(HERE)
RES = os.path.join(PROJ, "results")

from rev09_01_two_step_routes import resolve_results_file           # noqa: E402
from rev10_01_background_control import mantel_haenszel             # noqa: E402

ATLAS = os.path.join(PROJ, "data", "processed", "codon_atlas.csv")

# --------------------------------------------------------------------------- #
# Grantham (1974) Science 185:862-864.
#   D(i,j) = rho * [ alpha(c_i-c_j)^2 + beta(p_i-p_j)^2 + gamma(v_i-v_j)^2 ]^0.5
# Properties: c = composition, p = polarity, v = molecular volume.
# Computing from the formula rather than transcribing 190 table entries removes
# the transcription-error risk; the self-test below checks it against published
# values, including the matrix maximum D(Cys,Trp)=215.
# --------------------------------------------------------------------------- #
GRANTHAM_PROPS = {
    "A": (0.00, 8.1, 31.0),  "R": (0.65, 10.5, 124.0), "N": (1.33, 11.6, 56.0),
    "D": (1.38, 13.0, 54.0), "C": (2.75, 5.5, 55.0),   "Q": (0.89, 10.5, 85.0),
    "E": (0.92, 12.3, 83.0), "G": (0.74, 9.0, 3.0),    "H": (0.58, 10.4, 96.0),
    "I": (0.00, 5.2, 111.0), "L": (0.00, 4.9, 111.0),  "K": (0.33, 11.3, 119.0),
    "M": (0.00, 5.7, 105.0), "F": (0.00, 5.2, 132.0),  "P": (0.39, 8.0, 32.5),
    "S": (1.42, 9.2, 32.0),  "T": (0.71, 8.6, 61.0),   "W": (0.13, 5.4, 170.0),
    "Y": (0.20, 6.2, 136.0), "V": (0.00, 5.9, 84.0),
}
ALPHA, BETA, GAMMA, RHO = 1.833, 0.1018, 0.000399, 50.723


def grantham(a, b):
    if a not in GRANTHAM_PROPS or b not in GRANTHAM_PROPS:
        return None
    (c1, p1, v1), (c2, p2, v2) = GRANTHAM_PROPS[a], GRANTHAM_PROPS[b]
    return RHO * math.sqrt(ALPHA * (c1 - c2) ** 2
                           + BETA * (p1 - p2) ** 2
                           + GAMMA * (v1 - v2) ** 2)


def validate_grantham():
    """Published values from Grantham's Table 2. Refuse to run if any is off."""
    known = [("S", "R", 110), ("L", "I", 5), ("C", "W", 215),
             ("I", "V", 29), ("G", "W", 184), ("D", "E", 45)]
    print("  Grantham self-test (computed vs published):")
    bad = 0
    for a, b, ref in known:
        got = grantham(a, b)
        ok = abs(got - ref) <= 1.0
        bad += 0 if ok else 1
        print("    D(%s,%s) = %6.1f   published %3d   %s"
              % (a, b, got, ref, "OK" if ok else "MISMATCH"))
    if bad:
        sys.exit("  Grantham implementation is wrong. Stop.")
    print("  Grantham validated.")


def tertile_cut(values):
    """Return two cut points splitting sorted values into thirds."""
    v = sorted(values)
    if len(v) < 3:
        return (float("inf"), float("inf"))
    return (v[len(v) // 3], v[2 * len(v) // 3])


def band(x, cuts):
    if x is None:
        return "NA"
    return "lo" if x < cuts[0] else ("mid" if x < cuts[1] else "hi")


def run(rows, keyfn, label, out):
    by = defaultdict(lambda: [0, 0, 0, 0])
    for r in rows:
        by[keyfn(r)][(0 if r["used"] else 2) + (0 if r["accessible"] else 1)] += 1
    or_, lo, hi, n_inf = mantel_haenszel([tuple(v) for v in by.values()])
    line = ("  %-34s strata %6d  informative %4d   OR %7s  [%s, %s]"
            % (label, len(by), n_inf,
               "inf" if or_ == float("inf") else "%.2f" % or_,
               "%.2f" % lo if lo == lo else "na",
               "%.2f" % hi if hi == hi else "na"))
    print(line)
    out.append(line)
    return {"adjustment": label, "n_strata": len(by), "n_informative": n_inf,
            "odds_ratio": "inf" if or_ == float("inf") else round(or_, 3),
            "ci_lo": round(lo, 3) if lo == lo else "",
            "ci_hi": round(hi, 3) if hi == hi else ""}


def main():
    print("=" * 78)
    print("STEP 12.  Accessibility adjusted for tolerability, at species level")
    print("=" * 78)
    print()
    validate_grantham()
    print()

    ev_path = resolve_results_file("core_test_events")
    esm_path = resolve_results_file("esm_variant_effects")
    print()

    observed = set()
    with open(ev_path) as fh:
        for r in csv.DictReader(fh):
            observed.add((r["species"], r["gene"], r["site_label"], r["resistance_allele"]))

    esm = {}
    with open(esm_path) as fh:
        for r in csv.DictReader(fh):
            try:
                esm[(r["gene"], r["site_label"], r["candidate_aa"])] = float(r["esm_score"])
            except (ValueError, KeyError):
                pass

    # panel residue frequency: how often each residue is actually seen at a site
    site_counts = defaultdict(lambda: defaultdict(int))
    site_total = defaultdict(int)
    rowbuf = []
    with open(ATLAS) as fh:
        for r in csv.DictReader(fh):
            k = (r["gene"], r["site_label"])
            site_counts[k][r["species_aa"]] += 1
            site_total[k] += 1
            rowbuf.append(r)

    rows = []
    for r in rowbuf:
        if r["is_wildtype"] != "True" or r["already_resistant"] == "True":
            continue
        gene, site, allele = r["gene"], r["site_label"], r["resistance_allele"]
        k = (gene, site)
        g = grantham(r["wild_type_aa"], allele)
        freq = site_counts[k].get(allele, 0) / site_total[k] if site_total[k] else None
        rows.append({
            "species": r["species"], "gene": gene, "site_label": site,
            "site": "%s:%s" % (gene, site),
            "wild_type_aa": r["wild_type_aa"], "candidate_aa": allele,
            "accessible": r["one_step_accessible"] == "True",
            "used": (r["species"], gene, site, allele) in observed,
            "grantham": g, "esm_score": esm.get((gene, site, allele)),
            "panel_freq": freq,
        })

    n = len(rows)
    n_g = sum(1 for r in rows if r["grantham"] is not None)
    n_e = sum(1 for r in rows if r["esm_score"] is not None)
    n_f = sum(1 for r in rows if r["panel_freq"] is not None)
    print("-" * 78)
    print("COVARIATE COVERAGE")
    print("-" * 78)
    print("  rows (species x site x candidate allele) : %d" % n)
    print("  with Grantham distance                   : %d (%.1f%%)" % (n_g, 100 * n_g / n))
    print("  with an ESM-2 score                      : %d (%.1f%%)" % (n_e, 100 * n_e / n))
    print("  with panel residue frequency             : %d (%.1f%%)" % (n_f, 100 * n_f / n))
    if n_e < n:
        print()
        print("  ESM does not cover every row, because it was scored per gene x site x")
        print("  candidate residue on a reference protein rather than per species.")
        print("  The ESM-adjusted row below is therefore computed on a SUBSET; the")
        print("  Grantham and panel-frequency rows use all rows.")

    cuts_g = tertile_cut([r["grantham"] for r in rows if r["grantham"] is not None])
    cuts_e = tertile_cut([r["esm_score"] for r in rows if r["esm_score"] is not None])
    cuts_f = tertile_cut([r["panel_freq"] for r in rows if r["panel_freq"] is not None])

    out = []
    print()
    print("-" * 78)
    print("THE SPECIES-LEVEL ODDS RATIO, ADJUSTED")
    print("-" * 78)
    print()
    res = [
        run(rows, lambda r: r["site"], "site only (Step 10 baseline)", out),
        run(rows, lambda r: (r["site"], band(r["grantham"], cuts_g)),
            "site x Grantham tertile", out),
        run([r for r in rows if r["esm_score"] is not None],
            lambda r: (r["site"], band(r["esm_score"], cuts_e)),
            "site x ESM-2 tertile (subset)", out),
        run(rows, lambda r: (r["site"], band(r["panel_freq"], cuts_f)),
            "site x panel-frequency tertile", out),
        run(rows, lambda r: (r["site"], band(r["grantham"], cuts_g),
                             band(r["panel_freq"], cuts_f)),
            "site x Grantham x panel-frequency", out),
    ]

    base = res[0]["odds_ratio"]
    print()
    print("  Baseline is %s. An adjusted estimate that stays near it means the" % base)
    print("  covariate does not explain the association. Intervals widen as strata")
    print("  multiply, because strata without variation drop out; read the estimate.")

    # ---------------------------------------------------------------- methods
    meth = []
    meth.append("## ESM-2 variant-effect scoring, as actually performed")
    meth.append("")
    meth.append("Model: **esm2_t33_650M_UR50D** (ESM-2, 33 layers, 650M parameters,")
    meth.append("trained on UniRef50D), loaded with `esm.pretrained.load_model_and_alphabet`.")
    meth.append("")
    meth.append("Scoring: **masked marginals**. For a site i in protein P the residue at i")
    meth.append("is replaced by the mask token and the model is run once; the variant effect is")
    meth.append("")
    meth.append("    score(mut) = log p(mut | P_masked) - log p(wt | P_masked)")
    meth.append("")
    meth.append("so a negative score means the model considers the substitution less likely")
    meth.append("than the wild-type residue at that position. Both terms come from the same")
    meth.append("forward pass, which is why no sampling or ensembling is involved.")
    meth.append("")
    meth.append("Inputs: the curated REFERENCE protein for each gene, truncated to the window")
    meth.append("recorded per row in `esm_variant_effects.csv` (`window_start`, `window_end`).")
    meth.append("")
    meth.append("Coverage: %d scored gene x site x candidate-residue combinations." % len(esm))
    meth.append("")
    meth.append("### Limitations, stated")
    meth.append("")
    meth.append("1. The score is **site-level**: one reference protein per gene, with no")
    meth.append("   species dimension. It therefore cannot formally adjust a species-level")
    meth.append("   estimate, which is why Grantham distance and panel residue frequency are")
    meth.append("   reported alongside it here.")
    meth.append("2. A language-model likelihood is a statement about sequence plausibility")
    meth.append("   under the training distribution, not a measurement of fitness or of")
    meth.append("   protein function. It is used as a proxy and described as one.")
    meth.append("3. The model has no representation of codons, which is the property that")
    meth.append("   makes it useful here: it cannot see the variable the analysis is testing.")
    meth.append("4. Resistance-associated residues are themselves present in the training")
    meth.append("   data, so the model is not blind to the outcome being studied.")

    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "rev12_adjusted.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0].keys()))
        w.writeheader()
        for r in res:
            w.writerow(r)
    with open(os.path.join(RES, "rev12_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(RES, "rev12_methods_esm.md"), "w") as fh:
        fh.write("\n".join(meth) + "\n")
    with open(os.path.join(RES, "rev12_summary.md"), "w") as fh:
        fh.write("# Step 12: tolerability adjustment\n\n```\n" + "\n".join(out) + "\n```\n")

    print()
    print("=" * 78)
    print("WROTE")
    print("  results/rev12_adjusted.csv        the adjusted odds ratios")
    print("  results/rev12_rows.csv            every row with all three covariates")
    print("  results/rev12_methods_esm.md      the ESM specification, for the appendix")
    print("  results/rev12_summary.md")
    print("=" * 78)


if __name__ == "__main__":
    main()
