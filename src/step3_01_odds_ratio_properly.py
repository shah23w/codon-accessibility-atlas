#!/usr/bin/env python3
"""
step3_01_odds_ratio_properly.py
===============================
The accessibility odds ratio, computed five ways.

WHY FIVE ESTIMATORS
-------------------
The manuscript reports an accessibility odds ratio of 19.33 [8.17, 45.72] after
adjustment for ESM-2, from a pooled logistic model over 893 candidate
substitutions. The project's own statistical audit (finding M5) established
that:

  * the estimate is not an artefact - there is no separation, and Firth
    correction moves it only -7.5% - but
  * the whole thing rests on SIX observations: 42 of 47 sites have zero
    inaccessible-and-observed candidates, so moving that cell by one changes
    the odds ratio by about 20%;
  * a site-cluster bootstrap gives 95% CI [10.07, 97.98], so the published
    upper limit of ~46 is roughly two-fold too tight;
  * sites entered the analysis BECAUSE they carry at least one documented
    allele, and Spearman(pool size, n documented alleles) = +0.69, so the
    marginal (unstratified) estimate is confounded between sites and the
    site-conditional estimate is the correct target.

WHAT THIS SCRIPT DOES INSTEAD
-----------------------------
  A. Reproduces the crude 2x2 as a check.
  B. EXACT STRATIFIED TEST. Conditioning on each site's margins, the
     accessible-and-observed cell is hypergeometric. Summing over sites gives an
     exact distribution by convolution - no asymptotics, no sparse-data bias, and
     immune to the between-site confound because it conditions on site margins.
     It is a stronger claim than the GEE, and the figure to lead with.
  C. Mantel-Haenszel odds ratio stratified by site, with the
     Robins-Breslow-Greenland interval.
  D. Conditional (fixed-effects) logistic regression, strata = site, if
     statsmodels is available. This is the site-conditional estimate to quote.
  E. Site-cluster bootstrap, to get an honest interval.
  F. Leave-one-site-out, to show which sites carry the estimate.
  G. Sparse-cell diagnostics, including the +-1 sensitivity on the six
     observations the estimate depends on.
  H. The between-site confound, quantified.

INPUT
  results/multivariate_dataset_primary.csv
      893 rows: gene, site_label, site, wild_type_aa, candidate_aa,
      observed, accessible, frac_reaching, n_codons_at_site, esm_score

OUTPUTS
  results/step3_odds_ratio_estimators.csv
  results/step3_site_strata.csv
  results/step3_leave_one_site_out.csv
  results/step3_summary.md

Standard library only, except that statsmodels - if present - adds the
conditional logistic regression. Everything else runs anywhere.
"""

import csv
import math
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
DATA = os.path.join(RES, "multivariate_dataset_primary.csv")

OUT_EST = os.path.join(RES, "step3_odds_ratio_estimators.csv")
OUT_STRATA = os.path.join(RES, "step3_site_strata.csv")
OUT_LOSO = os.path.join(RES, "step3_leave_one_site_out.csv")
OUT_MD = os.path.join(RES, "step3_summary.md")

BOOTSTRAP_REPS = 5000
SEED = 20260926

if not os.path.exists(DATA):
    sys.exit("ERROR: missing input:\n  %s" % DATA)


def tf(v):
    return str(v).strip().lower() in ("true", "1", "yes", "t")


# ======================================================================
# load and build per-site 2x2 tables
# ======================================================================
with open(DATA, newline="") as fh:
    rows = list(csv.DictReader(fh))

for r in rows:
    r["_obs"] = tf(r["observed"])
    r["_acc"] = tf(r["accessible"])
    r["_site"] = r.get("site") or "%s:%s" % (r["gene"], r["site_label"])

print("=" * 76)
print("STEP 3  The accessibility odds ratio, done properly")
print("=" * 76)
print()
print("Candidate substitutions: %d   sites: %d"
      % (len(rows), len({r["_site"] for r in rows})))
print()

# ----------------------------------------------------------------------
# A. crude 2x2
# ----------------------------------------------------------------------
A = sum(1 for r in rows if r["_acc"] and r["_obs"])
B = sum(1 for r in rows if r["_acc"] and not r["_obs"])
C = sum(1 for r in rows if not r["_acc"] and r["_obs"])
D = sum(1 for r in rows if not r["_acc"] and not r["_obs"])
crude_or = (A / B) / (C / D) if B and C and D else float("nan")

print("-" * 76)
print("A.  CRUDE 2x2  (reproduction check)")
print("-" * 76)
print("                     observed    not observed     total")
print("  accessible       %8d    %12d  %8d" % (A, B, A + B))
print("  not accessible   %8d    %12d  %8d" % (C, D, C + D))
print("  total            %8d    %12d  %8d" % (A + C, B + D, len(rows)))
print()
print("  crude odds ratio = %.4f" % crude_or)
print("  P(observed | accessible)     = %.4f" % (A / (A + B)))
print("  P(observed | not accessible) = %.4f" % (C / (C + D)))
print()

# ----------------------------------------------------------------------
# build strata
# ----------------------------------------------------------------------
strata = {}
for r in rows:
    s = strata.setdefault(r["_site"], {"a": 0, "b": 0, "c": 0, "d": 0,
                                       "gene": r["gene"],
                                       "site_label": r["site_label"],
                                       "n_codons_at_site": r.get("n_codons_at_site", "")})
    if r["_acc"] and r["_obs"]:
        s["a"] += 1
    elif r["_acc"]:
        s["b"] += 1
    elif r["_obs"]:
        s["c"] += 1
    else:
        s["d"] += 1

for s in strata.values():
    s["N"] = s["a"] + s["b"] + s["c"] + s["d"]
    s["n_obs"] = s["a"] + s["c"]          # documented alleles at this site
    s["K_acc"] = s["a"] + s["b"]          # accessible candidates at this site
    # a stratum is informative for the conditional analysis only if both
    # margins allow the cell to vary
    lo = max(0, s["n_obs"] - (s["N"] - s["K_acc"]))
    hi = min(s["n_obs"], s["K_acc"])
    s["cell_lo"], s["cell_hi"] = lo, hi
    s["informative"] = hi > lo

n_inf = sum(1 for s in strata.values() if s["informative"])
zero_c = sum(1 for s in strata.values() if s["c"] == 0)

print("-" * 76)
print("G.  SPARSE-CELL DIAGNOSTIC  <-- read this before believing any interval")
print("-" * 76)
print("  sites total                                        : %d" % len(strata))
print("  sites informative for the site-conditional analysis : %d" % n_inf)
print("  sites with ZERO inaccessible-and-observed candidates: %d" % zero_c)
print("  total inaccessible-and-observed observations        : %d   <-- the whole"
      % C)
print("                                                            estimate rests")
print("                                                            on these")
print()


# ----------------------------------------------------------------------
# B. EXACT stratified test by hypergeometric convolution
# ----------------------------------------------------------------------
def hypergeom_pmf(N, K, n):
    """Distribution of (K-type count) in a sample of n drawn from N without
    replacement, returned as {x: probability}."""
    lo = max(0, n - (N - K))
    hi = min(n, K)
    denom = math.comb(N, n)
    return {x: math.comb(K, x) * math.comb(N - K, n - x) / denom
            for x in range(lo, hi + 1)}


exp_total = 0.0
var_total = 0.0
dist = {0: 1.0}
for s in strata.values():
    N, K, n = s["N"], s["K_acc"], s["n_obs"]
    if n == 0 or K == 0 or K == N or n == N:
        # margins force the cell; it contributes deterministically
        forced = s["a"]
        dist = {x + forced: p for x, p in dist.items()}
        exp_total += forced
        continue
    pmf = hypergeom_pmf(N, K, n)
    exp_total += n * K / N
    var_total += n * (K / N) * ((N - K) / N) * ((N - n) / (N - 1)) if N > 1 else 0.0
    nxt = defaultdict(float)
    for x, px in dist.items():
        for y, py in pmf.items():
            nxt[x + y] += px * py
    dist = dict(nxt)

p_exact = sum(p for x, p in dist.items() if x >= A)
sd_total = math.sqrt(var_total) if var_total > 0 else float("nan")
z = (A - exp_total) / sd_total if sd_total and not math.isnan(sd_total) else float("nan")

print("-" * 76)
print("B.  EXACT STRATIFIED TEST  (multiple hypergeometric, conditioned on site")
print("    margins)  <-- LEAD WITH THIS")
print("-" * 76)
print("  observed accessible-and-documented : %d" % A)
print("  expected under the exact null      : %.2f" % exp_total)
print("  standard deviation                 : %.2f" % sd_total)
print("  z                                  : %.2f" % z)
print("  one-sided exact P                  : %.4g" % p_exact)
print()
print("  Assumption-free: no asymptotics, no sparse-data bias, and immune to the")
print("  between-site confound because it conditions on each site's margins.")
print()


# ----------------------------------------------------------------------
# C. Mantel-Haenszel
# ----------------------------------------------------------------------
def mantel_haenszel(strat_list):
    R = S = 0.0
    for s in strat_list:
        N = s["N"]
        if N == 0:
            continue
        R += s["a"] * s["d"] / N
        S += s["b"] * s["c"] / N
    if S == 0 or R == 0:
        return float("nan"), (float("nan"), float("nan"))
    or_mh = R / S
    # Robins-Breslow-Greenland variance of ln(OR_MH)
    t1 = t2 = t3 = 0.0
    for s in strat_list:
        N = s["N"]
        if N == 0:
            continue
        P = (s["a"] + s["d"]) / N
        Q = (s["b"] + s["c"]) / N
        Ri = s["a"] * s["d"] / N
        Si = s["b"] * s["c"] / N
        t1 += P * Ri
        t2 += P * Si + Q * Ri
        t3 += Q * Si
    var_ln = t1 / (2 * R * R) + t2 / (2 * R * S) + t3 / (2 * S * S)
    se = math.sqrt(var_ln)
    return or_mh, (math.exp(math.log(or_mh) - 1.959964 * se),
                   math.exp(math.log(or_mh) + 1.959964 * se))


slist = list(strata.values())
or_mh, (mh_lo, mh_hi) = mantel_haenszel(slist)

print("-" * 76)
print("C.  MANTEL-HAENSZEL, stratified by site")
print("-" * 76)
print("  OR_MH = %.2f   95%% CI (Robins-Breslow-Greenland) [%.2f, %.2f]"
      % (or_mh, mh_lo, mh_hi))
print()

# ----------------------------------------------------------------------
# D. conditional logistic regression (optional)
# ----------------------------------------------------------------------
cond_or = cond_lo = cond_hi = float("nan")
cond_note = "statsmodels not available - skipped"
try:
    import numpy as np
    import pandas as pd
    from statsmodels.discrete.conditional_models import ConditionalLogit

    df = pd.DataFrame({
        "observed": [1 if r["_obs"] else 0 for r in rows],
        "accessible": [1.0 if r["_acc"] else 0.0 for r in rows],
        "esm": [float(r["esm_score"]) for r in rows],
        "site": [r["_site"] for r in rows],
    })
    m = ConditionalLogit(df["observed"], df[["accessible", "esm"]],
                         groups=df["site"]).fit(disp=False)
    b = m.params["accessible"]
    se = m.bse["accessible"]
    cond_or = math.exp(b)
    cond_lo = math.exp(b - 1.959964 * se)
    cond_hi = math.exp(b + 1.959964 * se)
    cond_note = ("conditional (fixed-effects) logit, strata = site, "
                 "adjusted for ESM-2; ESM-2 OR = %.4f" % math.exp(m.params["esm"]))
    print("-" * 76)
    print("D.  CONDITIONAL LOGISTIC REGRESSION, strata = site  <-- QUOTE THIS")
    print("-" * 76)
    print("  accessibility OR = %.2f   95%% CI [%.2f, %.2f]" % (cond_or, cond_lo, cond_hi))
    print("  ESM-2 OR         = %.4f  (per log-unit)" % math.exp(m.params["esm"]))
    print()
    print("  Sites entered the analysis because they carry a documented allele, so")
    print("  the marginal intercept is not meaningful and this is the correct target.")
    print()
except Exception as exc:  # noqa: BLE001
    cond_note = "not computed: %s" % exc
    print("-" * 76)
    print("D.  CONDITIONAL LOGISTIC REGRESSION  -- %s" % cond_note)
    print("-" * 76)
    print()

# ----------------------------------------------------------------------
# E. site-cluster bootstrap on OR_MH
# ----------------------------------------------------------------------
random.seed(SEED)
keys = list(strata.keys())
boots, undefined = [], 0
for _ in range(BOOTSTRAP_REPS):
    sample = [strata[random.choice(keys)] for _ in keys]
    o, _ci = mantel_haenszel(sample)
    if math.isnan(o):
        undefined += 1
    else:
        boots.append(o)
boots.sort()


def pctl(v, q):
    if not v:
        return float("nan")
    i = min(len(v) - 1, max(0, int(round(q * (len(v) - 1)))))
    return v[i]


b_lo, b_med, b_hi = pctl(boots, 0.025), pctl(boots, 0.5), pctl(boots, 0.975)

print("-" * 76)
print("E.  SITE-CLUSTER BOOTSTRAP on OR_MH  (%d resamples, seed %d)"
      % (BOOTSTRAP_REPS, SEED))
print("-" * 76)
print("  median %.2f   95%% interval [%.2f, %.2f]" % (b_med, b_lo, b_hi))
print("  resamples with an undefined odds ratio: %d of %d (%.1f%%)"
      % (undefined, BOOTSTRAP_REPS, 100.0 * undefined / BOOTSTRAP_REPS))
if undefined:
    print("  Those are resamples in which no site contributed a discordant pair.")
    print("  That they occur at all is itself evidence of how thin the data are.")
print()
print("  Compare with the published interval [8.17, 45.72]. If this upper bound is")
print("  much larger, the published interval was too tight and must be replaced.")
print()

# ----------------------------------------------------------------------
# F. leave-one-site-out
# ----------------------------------------------------------------------
loso = []
for k in keys:
    sub = [strata[j] for j in keys if j != k]
    o, (l, h) = mantel_haenszel(sub)
    loso.append({"site_dropped": k, "gene": strata[k]["gene"],
                 "or_mh": round(o, 4) if not math.isnan(o) else "",
                 "ci_lo": round(l, 4) if not math.isnan(l) else "",
                 "ci_hi": round(h, 4) if not math.isnan(h) else "",
                 "a": strata[k]["a"], "b": strata[k]["b"],
                 "c": strata[k]["c"], "d": strata[k]["d"]})
vals = [x["or_mh"] for x in loso if x["or_mh"] != ""]
print("-" * 76)
print("F.  LEAVE-ONE-SITE-OUT on OR_MH")
print("-" * 76)
print("  range across all %d drops: %.2f to %.2f" % (len(vals), min(vals), max(vals)))
movers = sorted([x for x in loso if x["or_mh"] != ""],
                key=lambda x: abs(x["or_mh"] - or_mh), reverse=True)[:5]
print("  sites that move the estimate most:")
for m in movers:
    print("    drop %-22s -> OR %6.2f   (its cells a=%d b=%d c=%d d=%d)"
          % (m["site_dropped"], m["or_mh"], m["a"], m["b"], m["c"], m["d"]))
print()

# ----------------------------------------------------------------------
# G2. the +-1 sensitivity on the six observations
# ----------------------------------------------------------------------
print("-" * 76)
print("G2. SENSITIVITY: move one observation")
print("-" * 76)
sens = []
for k in keys:
    s = strata[k]
    if s["b"] > 0:
        pert = dict(s)
        pert["c"] = s["c"] + 1
        pert["b"] = s["b"] - 1
        sub = [strata[j] if j != k else pert for j in keys]
        o, _ = mantel_haenszel(sub)
        if not math.isnan(o):
            sens.append((k, o))
if sens:
    os_ = sorted(o for _, o in sens)
    print("  Reclassifying ONE accessible-but-unobserved candidate as observed-and-")
    print("  inaccessible, at each site in turn, moves OR_MH to between %.2f and %.2f"
          % (os_[0], os_[-1]))
    print("  against a point estimate of %.2f - a swing of up to %.0f%%."
          % (or_mh, 100.0 * max(abs(os_[0] - or_mh), abs(os_[-1] - or_mh)) / or_mh))
    print("  This is why four significant figures cannot be justified.")
print()

# ----------------------------------------------------------------------
# H. the between-site confound
# ----------------------------------------------------------------------
def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    rk = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for t in range(i, j + 1):
            rk[order[t]] = avg
        i = j + 1
    return rk


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else float("nan")


pool, ndoc = [], []
for s in strata.values():
    try:
        pool.append(float(s["n_codons_at_site"]))
    except (TypeError, ValueError):
        continue
    ndoc.append(float(s["n_obs"]))
rho = pearson(rank(pool), rank(ndoc)) if len(pool) == len(ndoc) and pool else float("nan")

random.seed(SEED + 1)
if not math.isnan(rho):
    perm = 0
    for _ in range(20000):
        shuffled = ndoc[:]
        random.shuffle(shuffled)
        if abs(pearson(rank(pool), rank(shuffled))) >= abs(rho):
            perm += 1
    rho_p = (perm + 1) / 20001.0
else:
    rho_p = float("nan")

print("-" * 76)
print("H.  THE BETWEEN-SITE CONFOUND")
print("-" * 76)
print("  Spearman(codon pool size at site, number of documented alleles at site)")
print("    rho = %.3f   permutation P = %.4g  (n = %d sites, 20,000 shuffles)"
      % (rho, rho_p, len(pool)))
print()
print("  Sites entered the candidate table because they carry at least one")
print("  documented allele. A positive correlation here means the marginal")
print("  (unstratified) odds ratio mixes a within-site effect with a between-site")
print("  one, so the site-conditional estimate is the correct target.")
print()

# ======================================================================
# outputs
# ======================================================================
est_rows = [
    {"estimator": "crude_2x2_unadjusted", "odds_ratio": round(crude_or, 4),
     "ci_lo": "", "ci_hi": "", "p_value": "",
     "note": "a=%d b=%d c=%d d=%d; not adjusted, not stratified" % (A, B, C, D)},
    {"estimator": "exact_stratified_hypergeometric", "odds_ratio": "",
     "ci_lo": "", "ci_hi": "", "p_value": "%.6g" % p_exact,
     "note": "observed %d vs expected %.2f, sd %.2f, z %.2f; assumption-free, "
             "conditions on site margins" % (A, exp_total, sd_total, z)},
    {"estimator": "mantel_haenszel_by_site", "odds_ratio": round(or_mh, 4),
     "ci_lo": round(mh_lo, 4), "ci_hi": round(mh_hi, 4), "p_value": "",
     "note": "Robins-Breslow-Greenland interval"},
    {"estimator": "conditional_logit_by_site_adj_esm",
     "odds_ratio": "" if math.isnan(cond_or) else round(cond_or, 4),
     "ci_lo": "" if math.isnan(cond_lo) else round(cond_lo, 4),
     "ci_hi": "" if math.isnan(cond_hi) else round(cond_hi, 4),
     "p_value": "", "note": cond_note},
    {"estimator": "site_cluster_bootstrap_of_mh", "odds_ratio": round(b_med, 4),
     "ci_lo": round(b_lo, 4), "ci_hi": round(b_hi, 4), "p_value": "",
     "note": "%d resamples, seed %d, %d undefined" % (BOOTSTRAP_REPS, SEED, undefined)},
]
# No row here carries a typed odds ratio. A snapshot row, held as literals for
# comparison, is read sooner or later as the paper's own estimate, which freezes
# it at the value it had when it was typed. The adjusted odds ratio is read from
# multivariate_specifications, which step6_03 writes, and the row is gone so
# that nothing can read it again.
with open(OUT_EST, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["estimator", "odds_ratio", "ci_lo", "ci_hi",
                                       "p_value", "note"])
    w.writeheader()
    w.writerows(est_rows)

with open(OUT_STRATA, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["site", "gene", "site_label", "a_acc_obs", "b_acc_notobs",
                "c_inacc_obs", "d_inacc_notobs", "N", "n_documented",
                "n_accessible", "cell_min", "cell_max", "informative"])
    for k in sorted(strata):
        s = strata[k]
        w.writerow([k, s["gene"], s["site_label"], s["a"], s["b"], s["c"], s["d"],
                    s["N"], s["n_obs"], s["K_acc"], s["cell_lo"], s["cell_hi"],
                    s["informative"]])

with open(OUT_LOSO, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["site_dropped", "gene", "or_mh", "ci_lo",
                                       "ci_hi", "a", "b", "c", "d"])
    w.writeheader()
    w.writerows(loso)

with open(OUT_MD, "w") as fh:
    fh.write("# Step 3: the accessibility odds ratio\n\n")
    fh.write("Generated by `src/step3_01_odds_ratio_properly.py`.\n\n")
    fh.write("## Why the published figure had to be replaced\n\n")
    fh.write("The estimate rests on **%d observations**: %d of %d sites have zero "
             "inaccessible-and-observed candidates. Reclassifying a single "
             "candidate moves the Mantel-Haenszel odds ratio by up to %.0f%%. "
             "Four significant figures cannot be justified.\n\n"
             % (C, zero_c, len(strata),
                100.0 * max(abs(os_[0] - or_mh), abs(os_[-1] - or_mh)) / or_mh
                if sens else 0.0))
    fh.write("## Estimators\n\n")
    fh.write("| Estimator | OR | 95% CI | Note |\n|---|---|---|---|\n")
    fh.write("| Crude 2x2 | %.2f | - | a=%d b=%d c=%d d=%d |\n" % (crude_or, A, B, C, D))
    fh.write("| Mantel-Haenszel by site | %.2f | %.2f-%.2f | RBG interval |\n"
             % (or_mh, mh_lo, mh_hi))
    if not math.isnan(cond_or):
        fh.write("| **Conditional logit by site, adj. ESM-2** | **%.2f** | "
                 "**%.2f-%.2f** | the estimate to quote |\n"
                 % (cond_or, cond_lo, cond_hi))
    fh.write("| Site-cluster bootstrap of MH | %.2f | %.2f-%.2f | %d resamples |\n"
             % (b_med, b_lo, b_hi, BOOTSTRAP_REPS))
    fh.write("| Published GEE (for comparison) | 19.33 | 8.17-45.72 | currently in the MS |\n\n")
    fh.write("## The assumption-free test\n\n")
    fh.write("Conditioning on each site's margins, the accessible-and-documented "
             "count is a sum of independent hypergeometric variables, so its "
             "distribution is exact by convolution. **Observed %d against %.2f "
             "expected (sd %.2f, z = %.2f), one-sided exact *P* = %.3g.** This is "
             "immune to sparse-data bias and to the between-site confound, and it "
             "is a stronger claim than the regression.\n\n"
             % (A, exp_total, sd_total, z, p_exact))
    fh.write("## The between-site confound\n\n")
    fh.write("Spearman(codon pool size, documented alleles per site) = **%.3f** "
             "(permutation *P* = %.3g, n = %d sites). Sites entered the table "
             "because they carry a documented allele, so the marginal odds ratio "
             "mixes within- and between-site effects.\n\n" % (rho, rho_p, len(pool)))
    fh.write("## Leave-one-site-out\n\n")
    fh.write("OR_MH ranges from %.2f to %.2f across single-site drops.\n"
             % (min(vals), max(vals)))

print("=" * 76)
print("WROTE")
for p in (OUT_EST, OUT_STRATA, OUT_LOSO, OUT_MD):
    print("  %s" % os.path.relpath(p, ROOT))
print("=" * 76)
