#!/usr/bin/env python3
"""
step6_03_multivariate.py: does accessibility predict which allele is observed AFTER
adjusting for predicted tolerability?

The question this exists to answer
----------------------------------
Steps 3-5 show accessibility predicts the observed allele. Step 6.2 shows ESM-predicted
tolerability ALSO separates observed from unobserved alleles, and that the two axes are
correlated (Pearson r = +0.304 over 893 substitutions). Two univariate results plus a
correlation is not an answer: a reviewer will ask whether accessibility carries any
information once tolerability is accounted for. That is the only form of the question
that matters, and it is answered here.

Why the two axes are correlated at all is worth stating plainly: the genetic code is
error-minimising, so single-nucleotide-accessible substitutions are biochemically
conservative more often than chance, and conservative substitutions score better under a
protein language model. The confound is a property of the code, not of resistance.

Design
------
Unit of analysis: one candidate substitution at one verified site (site x candidate
amino acid). Outcome: is that substitution documented as a resistance allele in ANY
species. Predictors: mutational accessibility, and ESM-2 score as the tolerability axis.

Four analyses, deliberately overlapping, because each fails differently:

  1. NECESSARY-NOT-SUFFICIENT 2x2. P(observed | accessible) vs P(observed |
     inaccessible), with an odds ratio and Fisher exact test. Quantifies the claim
     that has so far only been asserted.
  2. ADJUSTED LOGISTIC REGRESSION. observed ~ accessible + esm_score, fitted by GEE
     with sites as clusters so standard errors are robust to site-level dependence
     (alleles at the same site are not independent observations).
  3. LIKELIHOOD-RATIO TEST. Does adding accessibility improve fit over an
     ESM-only model? This is the direct test of independent information.
  4. ESM-STRATIFIED MANTEL-HAENSZEL. A nonparametric check that makes no
     functional-form assumption: within strata of similar predicted tolerability,
     is accessibility still associated with being observed? If the association is a
     fitness artefact it disappears inside strata.

Effect attenuation (unadjusted OR vs adjusted OR) is reported: how much of the
accessibility signal is explained away by tolerability is itself the finding.

Outputs
-------
results/multivariate_dataset.csv
results/multivariate_results.csv
logs/step6_03_multivariate.log
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2
import statsmodels.api as sm
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.cov_struct import Independence
from statsmodels.stats.contingency_tables import StratifiedTable

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import DIR_PROC, DIR_REFS, DIR_RESULTS, Tee, ensure_dirs  # noqa: E402
from codon_engine import codon_reaches  # noqa: E402

AAS = "ACDEFGHIKLMNPQRSTVWY"


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default="codon_atlas_extended.csv")
    ap.add_argument("--record", default="resistance_record_merged.csv")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    TAG = args.tag

    ensure_dirs()
    log = Tee("step6_03_multivariate" + TAG)
    log("=== Step 6.3  accessibility adjusted for predicted tolerability ===\n")
    log(f"atlas  : {args.atlas}")
    log(f"record : {args.record}\n")

    # ---- assemble the substitution-level table ---------------------------
    esm = {(r["gene"], r["site_label"], r["candidate_aa"]): float(r["esm_score"])
           for r in read_csv(os.path.join(DIR_RESULTS, "esm_variant_effects.csv"))}

    documented = defaultdict(set)
    for r in read_csv(os.path.join(DIR_REFS, "reference_sites.csv")):
        documented[(r["gene"], r["site_label"])].add(r["resistance_allele"])
    for r in read_csv(os.path.join(DIR_PROC, args.record)):
        key = (r["gene"], r["site_label"])
        if key in documented:
            documented[key].add(r["resistance_allele"])

    codons_at, wt_of = defaultdict(list), {}
    for r in read_csv(os.path.join(DIR_PROC, args.atlas)):
        if r["is_wildtype"] != "True":
            continue
        key = (r["gene"], r["site_label"])
        codons_at[key].append(r["codon"])
        wt_of[key] = r["wild_type_aa"]

    rows = []
    for key, codons in sorted(codons_at.items()):
        gene, site = key
        wt = wt_of[key]
        doc = documented.get(key, set())
        for aa in AAS:
            if aa == wt:
                continue
            score = esm.get((gene, site, aa))
            if score is None:
                continue
            n_reach = sum(1 for c in codons if codon_reaches(c, aa))
            rows.append({
                "gene": gene, "site_label": site, "site": f"{gene}:{site}",
                "wild_type_aa": wt, "candidate_aa": aa,
                "observed": int(aa in doc),
                "accessible": int(n_reach > 0),
                "frac_reaching": n_reach / len(codons),
                "n_codons_at_site": len(codons),
                "esm_score": score,
            })
    df = pd.DataFrame(rows)
    log(f"substitutions: {len(df)}  across {df['site'].nunique()} sites")
    log(f"observed     : {int(df.observed.sum())}")
    log(f"accessible   : {int(df.accessible.sum())}")
    log(f"Pearson r(frac_reaching, esm_score) = "
        f"{df.frac_reaching.corr(df.esm_score):+.4f}")

    p0 = os.path.join(DIR_RESULTS, "multivariate_dataset" + TAG + ".csv")
    df.to_csv(p0, index=False)
    log(f"wrote {p0}")

    out = []

    # ---- 1. necessary but not sufficient ---------------------------------
    log("\n" + "=" * 70)
    log("1  NECESSARY BUT NOT SUFFICIENT")
    log("=" * 70)
    a = int(((df.accessible == 1) & (df.observed == 1)).sum())
    b = int(((df.accessible == 1) & (df.observed == 0)).sum())
    c = int(((df.accessible == 0) & (df.observed == 1)).sum())
    d = int(((df.accessible == 0) & (df.observed == 0)).sum())
    orr, pf = fisher_exact([[a, b], [c, d]])
    log(f"                     observed   not observed")
    log(f"  accessible         {a:8d}   {b:12d}")
    log(f"  not accessible     {c:8d}   {d:12d}")
    log(f"\n  P(observed | accessible)     = {a/(a+b):.4f}  ({a}/{a+b})")
    log(f"  P(observed | not accessible) = {c/(c+d):.4f}  ({c}/{c+d})")
    log(f"  odds ratio = {orr:.2f}   Fisher exact p = {pf:.3e}")
    log(f"\n  NECESSARY?     {c} of {a+c} observed alleles are inaccessible "
        f"({c/(a+c):.1%}) -> accessibility is not strictly required")
    log(f"  SUFFICIENT?    {b} of {a+b} accessible alleles are never observed "
        f"({b/(a+b):.1%}) -> plainly not sufficient")
    out.append({"analysis": "2x2_necessary_sufficient", "estimate": orr,
                "p_value": pf, "n": len(df),
                "detail": f"P(obs|acc)={a/(a+b):.4f} P(obs|inacc)={c/(c+d):.4f}"})
    # the raw cells, so a figure never has to recover them from log prose
    for nm, val in (("2x2_accessible_observed", a), ("2x2_accessible_total", a + b),
                    ("2x2_inaccessible_observed", c),
                    ("2x2_inaccessible_total", c + d)):
        out.append({"analysis": nm, "estimate": val, "p_value": "", "n": len(df),
                    "detail": "cell count of the necessary/sufficient 2x2"})

    # ---- 2. adjusted logistic regression ---------------------------------
    log("\n" + "=" * 70)
    log("2  LOGISTIC REGRESSION, GEE with site clusters")
    log("=" * 70)
    df["site_id"] = pd.factorize(df["site"])[0]
    df = df.sort_values("site_id").reset_index(drop=True)

    specs = []

    def fit(formula_cols, label, data=None):
        d = df if data is None else data
        X = sm.add_constant(d[formula_cols].astype(float), has_constant="add")
        m = GEE(d["observed"].astype(float), X, groups=d["site_id"],
                family=Binomial(), cov_struct=Independence()).fit()
        log(f"\n  {label}")
        for name in X.columns:
            coef, se, p = m.params[name], m.bse[name], m.pvalues[name]
            lo, hi = coef - 1.96 * se, coef + 1.96 * se
            if name == "const":
                log(f"    {name:16} beta={coef:+.4f}  p={p:.3g}")
            else:
                log(f"    {name:16} beta={coef:+.4f}  OR={np.exp(coef):7.3f}  "
                    f"95% CI [{np.exp(lo):.3f}, {np.exp(hi):.3f}]  p={p:.3g}")
                specs.append({"specification": label, "term": name,
                              "odds_ratio": round(float(np.exp(coef)), 4),
                              "ci_lo": round(float(np.exp(lo)), 4),
                              "ci_hi": round(float(np.exp(hi)), 4),
                              "p_value": float(p), "n": len(d)})
        return m

    m_acc = fit(["accessible"], "accessibility alone")
    m_esm = fit(["esm_score"], "tolerability alone")
    m_both = fit(["accessible", "esm_score"], "BOTH (the question)")

    or_un = float(np.exp(m_acc.params["accessible"]))
    or_adj = float(np.exp(m_both.params["accessible"]))
    p_adj = float(m_both.pvalues["accessible"])
    log(f"\n  accessibility OR unadjusted = {or_un:.3f}")
    log(f"  accessibility OR adjusted   = {or_adj:.3f}   p = {p_adj:.3g}")
    att = 100 * (1 - (np.log(or_adj) / np.log(or_un))) if or_un > 1 else float("nan")
    log(f"  attenuation of log-OR after adjustment: {att:.1f}%")
    out.append({"analysis": "gee_adjusted_accessible", "estimate": or_adj,
                "p_value": p_adj, "n": len(df),
                "detail": f"unadjusted OR={or_un:.3f}; attenuation={att:.1f}%"})
    out.append({"analysis": "gee_adjusted_esm",
                "estimate": float(np.exp(m_both.params["esm_score"])),
                "p_value": float(m_both.pvalues["esm_score"]), "n": len(df),
                "detail": "per 1 log-unit of ESM score, adjusted for accessibility"})

    # ESM-2 scores one reference protein per gene, and acetylcholinesterase's
    # is Torpedo californica P04058: the frame the ace-1 literature numbers
    # against, and the only reference in this paper that is not an arthropod.
    # Refitting without that gene bounds what a vertebrate orthologue's scores
    # can be doing to the adjusted estimate. It is an exclusion of rows, not a
    # rescoring, so it costs no model inference.
    _no_ache = df[df["gene"] != "AChE"]
    if len(_no_ache) == len(df) or _no_ache.empty:
        sys.exit("FATAL: excluding gene == 'AChE' removes %d of %d rows. The "
                 "vertebrate-reference sensitivity fit quoted in the Methods "
                 "cannot be built." % (len(df) - len(_no_ache), len(df)))
    m_no_ache = fit(["accessible", "esm_score"],
                    "BOTH, dropping AChE (vertebrate ESM-2 reference)",
                    data=_no_ache)
    _excl = "AChE is scored on the Torpedo californica reference; %d of %d "
    _excl = _excl % (len(df) - len(_no_ache), len(df))
    out.append({"analysis": "gee_adjusted_accessible_no_ache",
                "estimate": float(np.exp(m_no_ache.params["accessible"])),
                "p_value": float(m_no_ache.pvalues["accessible"]),
                "n": len(_no_ache),
                "detail": _excl + "candidate rows excluded"})
    out.append({"analysis": "gee_adjusted_esm_no_ache",
                "estimate": float(np.exp(m_no_ache.params["esm_score"])),
                "p_value": float(m_no_ache.pvalues["esm_score"]),
                "n": len(_no_ache),
                "detail": _excl + "candidate rows excluded; per 1 log-unit "
                                  "of ESM score, adjusted for accessibility"})

    # ---- 3. likelihood-ratio test ----------------------------------------
    log("\n" + "=" * 70)
    log("3  LIKELIHOOD-RATIO TEST  (ordinary logit; does accessibility add fit?)")
    log("=" * 70)
    y = df["observed"].astype(float)
    X1 = sm.add_constant(df[["esm_score"]].astype(float), has_constant="add")
    X2 = sm.add_constant(df[["esm_score", "accessible"]].astype(float),
                         has_constant="add")
    r1 = sm.Logit(y, X1).fit(disp=0)
    r2 = sm.Logit(y, X2).fit(disp=0)
    lr = 2 * (r2.llf - r1.llf)
    p_lr = float(chi2.sf(lr, 1))
    log(f"  ESM only          log-likelihood = {r1.llf:.3f}")
    log(f"  ESM + accessible  log-likelihood = {r2.llf:.3f}")
    log(f"  LR chi2(1) = {lr:.3f}   p = {p_lr:.3e}")
    log(f"  pseudo-R2: {r1.prsquared:.4f} -> {r2.prsquared:.4f}")
    out.append({"analysis": "likelihood_ratio_accessible", "estimate": lr,
                "p_value": p_lr, "n": len(df),
                "detail": f"pseudoR2 {r1.prsquared:.4f}->{r2.prsquared:.4f}"})

    # ---- 4. ESM-stratified Mantel-Haenszel -------------------------------
    log("\n" + "=" * 70)
    log("4  MANTEL-HAENSZEL, stratified by predicted tolerability")
    log("=" * 70)
    df["esm_stratum"] = pd.qcut(df["esm_score"], 5, labels=False, duplicates="drop")
    tables = []
    log(f"  {'stratum':8} {'ESM range':>22} {'acc/obs':>9} {'acc/unobs':>10} "
        f"{'inacc/obs':>10} {'inacc/unobs':>12}")
    for s, g in df.groupby("esm_stratum"):
        aa_ = int(((g.accessible == 1) & (g.observed == 1)).sum())
        bb = int(((g.accessible == 1) & (g.observed == 0)).sum())
        cc = int(((g.accessible == 0) & (g.observed == 1)).sum())
        dd = int(((g.accessible == 0) & (g.observed == 0)).sum())
        rng = f"[{g.esm_score.min():.2f},{g.esm_score.max():.2f}]"
        log(f"  {int(s):8d} {rng:>22} {aa_:>9d} {bb:>10d} {cc:>10d} {dd:>12d}")
        if (aa_ + cc) > 0 and (bb + dd) > 0:
            tables.append(np.array([[aa_, bb], [cc, dd]]))
    if len(tables) >= 2:
        st = StratifiedTable(tables)
        log(f"\n  Mantel-Haenszel pooled OR = {st.oddsratio_pooled:.3f}")
        ci = st.oddsratio_pooled_confint()
        log(f"  95% CI [{ci[0]:.3f}, {ci[1]:.3f}]")
        tp = st.test_null_odds()
        # statsmodels can underflow this p-value to exactly 0; recompute from the
        # statistic so the reported value is a real number, not a display artefact.
        p_mh = float(chi2.sf(tp.statistic, 1))
        shown = f"{p_mh:.3g}" if p_mh > 0 else "<1e-300"
        log(f"  test of OR=1: chi2 = {tp.statistic:.3f}  p = {shown}")
        hom = st.test_equal_odds()
        log(f"  homogeneity across strata: chi2 = {hom.statistic:.3f} "
            f"p = {hom.pvalue:.3g}")
        out.append({"analysis": "mantel_haenszel_esm_stratified",
                    "estimate": float(st.oddsratio_pooled),
                    "p_value": p_mh, "n": len(df),
                    "detail": f"95% CI [{ci[0]:.3f},{ci[1]:.3f}]; "
                              f"{len(tables)} strata"})
    else:
        log("  too few usable strata")

    p1 = os.path.join(DIR_RESULTS, "multivariate_results" + TAG + ".csv")
    with open(p1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    sp = os.path.join(DIR_RESULTS, "multivariate_specifications" + TAG + ".csv")
    with open(sp, "w", newline="") as fh:
        w2 = csv.DictWriter(fh, fieldnames=list(specs[0].keys()))
        w2.writeheader(); w2.writerows(specs)
    log(f"\nwrote {p1}")
    log(f"wrote {sp}  ({len(specs)} model terms)")
    log("\nStep 6.3 complete.")


if __name__ == "__main__":
    main()
