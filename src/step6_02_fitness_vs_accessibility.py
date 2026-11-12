#!/usr/bin/env python3
"""
step6_02_fitness_vs_accessibility.py: Aim 3: is it accessibility, or is it fitness?

Runs in the MAIN environment (stdlib + scipy only); consumes the scores written by
step6_01, which runs in the isolated ESM venv.

The alternative hypothesis
--------------------------
Everything in Steps 3-5 is compatible with a simpler story: accessible alleles are just
the fitter ones, and mutational accessibility is incidental. This step is built to
distinguish those, and it can come out either way.

Every candidate substitution at every verified site is placed in exactly one class:

  OBSERVED              documented as a resistance allele in at least one species
  ACCESSIBLE_UNOBSERVED reachable in one step from >=1 species' codon, never documented
  INACCESSIBLE          documented nowhere AND not reachable in one step from any codon
                        segregating at that site

and scored on an independent fitness axis (ESM-2), which never sees a codon.

The discriminating comparison
-----------------------------
  * If FITNESS drives which allele appears, OBSERVED alleles should score markedly higher
    than ACCESSIBLE_UNOBSERVED ones - the accessible ones exist as mutational options but
    are rejected because they are worse.
  * If ACCESSIBILITY drives it, OBSERVED and ACCESSIBLE_UNOBSERVED should be
    indistinguishable - both are viable, and what separates them is only whether a species
    happened to carry a codon that could reach them.

A null result in that comparison is therefore EVIDENCE FOR the mutation-supply model, not
an absence of evidence. To keep that interpretation honest the test is two-sided,
equivalence is assessed explicitly (TOST) rather than inferred from p > 0.05, and the
positive control below must succeed.

POSITIVE CONTROL: the same ESM scores must separate documented resistance alleles from
random amino acids at the same sites. If they cannot, the fitness axis is uninformative
and the main comparison means nothing. This guards against concluding "no difference"
merely because the predictor is noise.

Outputs
-------
results/aim3_fitness_classes.csv       every candidate with its class and score
results/aim3_discrimination.csv        the statistical comparisons
logs/step6_02_fitness_vs_accessibility.log
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import DIR_PROC, DIR_REFS, DIR_RESULTS, Tee, ensure_dirs  # noqa: E402
from codon_engine import codon_reaches  # noqa: E402

SEED = 20260729
N_BOOT_EQ = 10000        # bootstrap resamples for the equivalence CI
N_PERM = 100_000


def mannwhitney(a, b):
    """Two-sided Mann-Whitney U with a normal approximation; returns (U, p, rank-biserial)."""
    if not a or not b:
        return None
    comb = [(v, 0) for v in a] + [(v, 1) for v in b]
    comb.sort(key=lambda t: t[0])
    ranks, i = [0.0] * len(comb), 0
    while i < len(comb):
        j = i
        while j + 1 < len(comb) and comb[j + 1][0] == comb[i][0]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[k] = r
        i = j + 1
    r1 = sum(rk for rk, (_v, g) in zip(ranks, comb) if g == 0)
    n1, n2 = len(a), len(b)
    u1 = r1 - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (u1 - mu) / sd if sd else 0.0
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return u1, p, (2 * u1 / (n1 * n2)) - 1


def median(v):
    s = sorted(v)
    n = len(s)
    return (s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2) if n else float("nan")


def main():
    # Optional: treat alleles documented in a resistance record as OBSERVED too, not
    # only those in reference_sites.csv. This matters because a genuinely observed
    # allele misclassified as ACCESSIBLE_UNOBSERVED biases the headline comparison
    # TOWARD the accessibility model, which is the direction that would flatter the
    # hypothesis. Default is reference-only, reproducing the original run.
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", default="", help="extra source of documented alleles")
    ap.add_argument("--atlas", default="codon_atlas.csv")
    ap.add_argument("--tag", default="", help="suffix for output files and log")
    args = ap.parse_args()
    TAG = args.tag
    ensure_dirs()
    log = Tee("step6_02_fitness_vs_accessibility" + TAG)
    log("=== Step 6.2  Aim 3: accessibility vs fitness ===\n")
    rng = random.Random(SEED)

    esm_path = os.path.join(DIR_RESULTS, "esm_variant_effects.csv")
    if not os.path.exists(esm_path):
        log("esm_variant_effects.csv not found - run step6_01 in the ESM venv first")
        return
    esm = {}
    with open(esm_path) as fh:
        for r in csv.DictReader(fh):
            esm[(r["gene"], r["site_label"], r["candidate_aa"])] = float(r["esm_score"])
    log(f"ESM scores loaded: {len(esm)}")

    # documented resistance alleles per site
    documented = defaultdict(set)
    with open(os.path.join(DIR_REFS, "reference_sites.csv")) as fh:
        for r in csv.DictReader(fh):
            documented[(r["gene"], r["site_label"])].add(r["resistance_allele"])
    n_ref = sum(len(v) for v in documented.values())
    log(f"documented alleles from reference_sites.csv: {n_ref}")
    if args.record:
        added = []
        with open(os.path.join(DIR_PROC, args.record)) as fh:
            for r in csv.DictReader(fh):
                key = (r["gene"], r["site_label"])
                if key in documented and r["resistance_allele"] not in documented[key]:
                    documented[key].add(r["resistance_allele"])
                    added.append(f"{r['gene']} {r['site_label']}->{r['resistance_allele']}")
        for a in sorted(set(added)):
            log(f"  + documented by {args.record}: {a}")
        log(f"documented alleles after merge: {sum(len(v) for v in documented.values())}")

    # codons actually segregating at each site across the panel
    codons_at = defaultdict(set)
    wt_of = {}
    with open(os.path.join(DIR_PROC, args.atlas)) as fh:
        for r in csv.DictReader(fh):
            if r["is_wildtype"] != "True":
                continue
            key = (r["gene"], r["site_label"])
            codons_at[key].add(r["codon"])
            wt_of[key] = r["wild_type_aa"]
    log(f"sites with codon data: {len(codons_at)}")

    rows = []
    for key, codons in sorted(codons_at.items()):
        gene, site = key
        wt = wt_of[key]
        doc = documented.get(key, set())
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            if aa == wt:
                continue
            score = esm.get((gene, site, aa))
            if score is None:
                continue
            n_reaching = sum(1 for c in codons if codon_reaches(c, aa))
            # A documented allele is not necessarily accessible: six of them are
            # obligate two-step. Splitting the class here makes the discriminating test
            # accessibility-matched by construction rather than by assertion, and keeps
            # the six visible as a class of their own.
            if aa in doc and n_reaching > 0:
                cls = "OBSERVED_ACCESSIBLE"
            elif aa in doc:
                cls = "OBSERVED_INACCESSIBLE"
            elif n_reaching > 0:
                cls = "ACCESSIBLE_UNOBSERVED"
            else:
                cls = "INACCESSIBLE"
            rows.append({
                "gene": gene, "site_label": site, "wild_type_aa": wt,
                "candidate_aa": aa, "class": cls,
                "n_codons_reaching": n_reaching,
                "n_codons_at_site": len(codons),
                "frac_codons_reaching": round(n_reaching / len(codons), 4),
                "esm_score": score,
            })

    out = os.path.join(DIR_RESULTS, "aim3_fitness_classes" + TAG + ".csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    log(f"wrote {out}  ({len(rows)} candidate substitutions)")

    by_cls = defaultdict(list)
    for r in rows:
        by_cls[r["class"]].append(r["esm_score"])
    log("\n--- ESM score by class (higher = better tolerated) ---")
    log(f"  {'class':24s} {'n':>5s} {'median':>9s} {'mean':>9s}")
    for c in ("OBSERVED_ACCESSIBLE", "OBSERVED_INACCESSIBLE",
              "ACCESSIBLE_UNOBSERVED", "INACCESSIBLE"):
        v = by_cls.get(c, [])
        if v:
            log(f"  {c:24s} {len(v):5d} {median(v):9.3f} {sum(v)/len(v):9.3f}")

    results = []

    # ---- POSITIVE CONTROL --------------------------------------------------
    log("\n" + "=" * 70)
    log("POSITIVE CONTROL  can ESM separate documented alleles from random ones?")
    log("=" * 70)
    obs_acc = by_cls.get("OBSERVED_ACCESSIBLE", [])
    obs_inacc = by_cls.get("OBSERVED_INACCESSIBLE", [])
    obs = obs_acc + obs_inacc          # all documented alleles, for the positive control
    others = by_cls.get("ACCESSIBLE_UNOBSERVED", []) + by_cls.get("INACCESSIBLE", [])
    if obs and others:
        u, p, rb = mannwhitney(obs, others)
        log(f"  OBSERVED (n={len(obs)}) median {median(obs):.3f}  vs  "
            f"all others (n={len(others)}) median {median(others):.3f}")
        log(f"  Mann-Whitney p = {p:.3g}   rank-biserial r = {rb:+.3f}")
        ok = p < 0.05 and rb > 0
        log(f"  -> fitness axis is {'INFORMATIVE' if ok else 'NOT clearly informative'}")
        results.append({"comparison": "positive_control_observed_vs_all_others",
                        "n_a": len(obs), "n_b": len(others),
                        "median_a": round(median(obs), 4), "median_b": round(median(others), 4),
                        "p_value": p, "rank_biserial": round(rb, 4),
                        "verdict": "informative" if ok else "uninformative"})

    # ---- THE DISCRIMINATING TEST ------------------------------------------
    log("\n" + "=" * 70)
    log("DISCRIMINATING TEST  OBSERVED_ACCESSIBLE vs ACCESSIBLE_UNOBSERVED")
    log("=" * 70)
    log("  Restricted to documented alleles that are themselves accessible, so the two")
    log("  classes are matched on accessibility and differ only in whether the allele")
    log(f"  was used. The {len(obs_inacc)} documented but inaccessible alleles are")
    log("  excluded here and reported separately below.")
    au = by_cls.get("ACCESSIBLE_UNOBSERVED", [])
    obs = obs_acc
    if obs and au:
        u, p, rb = mannwhitney(obs, au)
        log(f"  OBSERVED_ACCESSIBLE   n={len(obs):4d}  median {median(obs):7.3f}")
        log(f"  ACCESSIBLE_UNOBSERVED n={len(au):4d}  median {median(au):7.3f}")
        log(f"  difference in medians : {median(obs)-median(au):+.3f}")
        log(f"  Mann-Whitney p = {p:.4g}   rank-biserial r = {rb:+.3f}")
        log("")
        if p >= 0.05:
            log("  -> NO detectable fitness difference.")
            log("     Consistent with mutation supply: the two classes are equally viable,")
            log("     and what separates them is whether a codon could reach them.")
        else:
            log("  -> A fitness difference IS detectable.")
            log("     Fitness cannot be excluded as a driver; report both effects.")
        results.append({"comparison": "observed_vs_accessible_unobserved",
                        "n_a": len(obs), "n_b": len(au),
                        "median_a": round(median(obs), 4), "median_b": round(median(au), 4),
                        "p_value": p, "rank_biserial": round(rb, 4),
                        "verdict": "no_fitness_difference" if p >= 0.05 else "fitness_difference"})

        # EQUIVALENCE, done properly.
        #
        # The previous block here was labelled "TOST-style" but was arithmetically THE SAME
        # TEST as the Mann-Whitney reported a few lines above. Rank-biserial is returned as
        # (2*u1/(n1*n2)) - 1, a monotone affine function of U, so permuting the labels and
        # counting |r_perm| >= |r_obs| IS the exact permutation Mann-Whitney. It was
        # reproduced to the digit: the permutation gave 4.29996e-04 and the stored value was
        # 0.00042999550004499957.
        #
        # That tests H0: NO DIFFERENCE. Equivalence testing requires the opposite null -
        # H0: the difference is AT LEAST as large as the bound - so the old block could
        # never have demonstrated equivalence, whatever it returned. Its "verdict" was also
        # a bare point estimate compared with 0.2, with no interval, logically disconnected
        # from the p-value printed beside it.
        #
        # A real equivalence claim needs a CONFIDENCE INTERVAL on the effect size: the
        # effect is equivalent-to-zero only if the whole interval lies inside the bound.
        bound = 0.2                     # |r| < 0.2 treated as no meaningful difference
        boots = []
        for _ in range(N_BOOT_EQ):
            ra = [obs[rng.randrange(len(obs))] for _ in range(len(obs))]
            rb_ = [au[rng.randrange(len(au))] for _ in range(len(au))]
            _u3, _p3, rbb = mannwhitney(ra, rb_)
            boots.append(rbb)
        boots.sort()
        lo90 = boots[int(0.05 * len(boots))]
        hi90 = boots[int(0.95 * len(boots))]
        # TOST at alpha=0.05 <=> the 90% CI lies entirely within (-bound, +bound).
        # A binary equivalent/not-equivalent verdict is wrong: an interval that STRADDLES
        # the bound establishes neither equivalence nor a difference exceeding it. The
        # verdict is three-way.
        equivalent = (lo90 > -bound) and (hi90 < bound)
        exceeds = (lo90 > bound) or (hi90 < -bound)
        verdict = ("equivalent" if equivalent
                   else "difference_exceeds_bound" if exceeds
                   else "inconclusive")
        log(f"\n  rank-biserial r = {rb:+.3f}")
        log(f"  90% bootstrap CI [{lo90:+.3f}, {hi90:+.3f}]  ({N_BOOT_EQ} resamples)")
        log(f"  equivalence bound |r| < {bound}")
        log(f"  -> {verdict.upper()}")
        if verdict == "equivalent":
            log("     the interval lies entirely inside the bound")
        elif verdict == "difference_exceeds_bound":
            log("     the interval lies entirely outside the bound, so the two classes")
            log("     differ by more than the pre-specified negligible effect")
        else:
            log("     the interval SPANS the bound: these data neither establish")
            log("     equivalence nor establish a difference exceeding it")
        results.append({"comparison": "observed_vs_accessible_unobserved_equivalence",
                        "n_a": len(obs), "n_b": len(au),
                        "median_a": "", "median_b": "",
                        "p_value": "", "rank_biserial": round(rb, 4),
                        "ci90_low": round(lo90, 4), "ci90_high": round(hi90, 4),
                        "verdict": verdict})

    # ---- accessibility vs fitness as competing predictors -------------------
    log("\n" + "=" * 70)
    log("ARE ACCESSIBILITY AND FITNESS CORRELATED AT ALL?")
    log("=" * 70)
    xs = [r["frac_codons_reaching"] for r in rows]
    ys = [r["esm_score"] for r in rows]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = math.sqrt(sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys))
    r_pear = num / den if den else float("nan")
    log(f"  Pearson r(codon accessibility, ESM score) = {r_pear:+.4f}  (n={n})")
    # The interpretation follows the number rather than printing regardless.
    if abs(r_pear) < 0.10:
        log("  -> near zero: the two axes are effectively independent")
    elif abs(r_pear) < 0.35:
        log(f"  -> r = {r_pear:+.3f} is SMALL BUT NOT ZERO. The axes are partly")
        log("     confounded, which is expected: the genetic code is error-minimising, so")
        log("     one-step-accessible substitutions are conservative more often than")
        log("     chance and score better under a protein language model. Accessibility")
        log("     must therefore be tested ADJUSTED for tolerability (step6_03), not")
        log("     asserted to be independent of it.")
    else:
        log(f"  -> r = {r_pear:+.3f}: the axes are substantially confounded; univariate")
        log("     comparisons on either axis alone are not interpretable")
    results.append({"comparison": "correlation_accessibility_vs_esm",
                    "n_a": n, "n_b": "", "median_a": "", "median_b": "",
                    "p_value": "", "rank_biserial": round(r_pear, 4),
                    "verdict": "independent" if abs(r_pear) < 0.2 else "correlated"})

    rp = os.path.join(DIR_RESULTS, "aim3_discrimination" + TAG + ".csv")
    with open(rp, "w", newline="") as fh:
        # union of keys across rows: the equivalence row carries CI columns the others
        # do not, and taking fieldnames from results[0] alone raises on write.
        fields = list(dict.fromkeys(k for r in results for k in r))
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in fields})
    log(f"\nwrote {rp}")
    log.close()


if __name__ == "__main__":
    main()
