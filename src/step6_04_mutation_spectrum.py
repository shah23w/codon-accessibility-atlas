#!/usr/bin/env python3
"""
step6_04_mutation_spectrum.py: is the signal mutational SUPPLY, or amino-acid
similarity wearing a mutational disguise?

The problem with everything so far
----------------------------------
Accessibility and tolerability are correlated (r = +0.31) because the genetic code is
error-minimising. Step 6.3 showed accessibility survives adjustment for tolerability
(OR ~20, attenuation 6.3%), which is strong, but a determined sceptic can still argue
that "one nucleotide away" is a proxy for "biochemically similar" and that the whole
result is a restatement of protein tolerance.

This step closes that argument with a variable the sceptic's explanation cannot touch.

The test
--------
Transitions (A<->G, C<->T) occur several-fold more often than transversions in real
genomes. That is a property of DNA replication and repair chemistry -- of the MUTATION
PROCESS -- and it is invisible to protein-level reasoning. A protein language model
scores amino acids; it has no representation of which nucleotide change produced them.

So, holding accessibility CONSTANT (restricting to substitutions already reachable in
one step), ask:

    among alleles a species could reach, are the ones actually observed
    preferentially those reachable by a TRANSITION?

  * If the pattern is mutational supply, yes: transition-reachable alleles arise more
    often per generation and should be over-represented among observed alleles.
  * If the pattern is protein tolerance, there is no reason for a transition-reachable
    allele to be favoured over a transversion-reachable one of equal tolerability.

The prediction is directional, mechanistic, and orthogonal to the confound.

Analyses
--------
  1. Substitution level (site x candidate amino acid, accessible only): transition
     fraction of observed vs never-observed alleles.
  2. GEE logistic regression adjusted for ESM score, sites as clusters.
  3. Event level (one row per documented species x site x allele): the same contrast
     using each species' OWN codon, compared against the alleles that species could
     have reached but did not.
  4. Within-site permutation, in the style used throughout this project: shuffle which
     accessible allele is "observed" within each site and recompute.

Outputs
-------
results/mutation_spectrum_substitutions.csv
results/mutation_spectrum_results.csv
logs/step6_04_mutation_spectrum.log
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.genmod.cov_struct import Independence
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.generalized_estimating_equations import GEE

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import DIR_PROC, DIR_REFS, DIR_RESULTS, Tee, ensure_dirs  # noqa: E402
from codon_engine import codon_path_counts, codon_reaches  # noqa: E402

AAS = "ACDEFGHIKLMNPQRSTVWY"
N_PERM = 20000
SEED = 20260730


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default="codon_atlas_extended.csv")
    ap.add_argument("--record", default="resistance_record_merged.csv")
    ap.add_argument("--events", default="core_test_events_primary.csv")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    TAG = args.tag
    rng = random.Random(SEED)

    ensure_dirs()
    log = Tee("step6_04_mutation_spectrum" + TAG)
    log("=== Step 6.4  mutation spectrum: transitions vs transversions ===\n")

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

    # ---------- 1. substitution level, accessible only --------------------
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
            ts = tv = n_reach = 0
            for c in codons:
                n, t1, t2 = codon_path_counts(c, aa)
                if n:
                    n_reach += 1
                    ts += t1
                    tv += t2
            if n_reach == 0:
                continue                      # accessibility held constant
            rows.append({
                "gene": gene, "site_label": site, "site": f"{gene}:{site}",
                "wild_type_aa": wt, "candidate_aa": aa,
                "observed": int(aa in doc),
                "n_transition_paths": ts, "n_transversion_paths": tv,
                "ts_fraction": ts / (ts + tv) if (ts + tv) else 0.0,
                "transition_only": int(tv == 0),
                "frac_codons_reaching": n_reach / len(codons),
                "esm_score": score,
            })
    df = pd.DataFrame(rows)
    log(f"accessible substitutions: {len(df)}  observed: {int(df.observed.sum())}")
    log(f"corr(ts_fraction, esm_score) = {df.ts_fraction.corr(df.esm_score):+.4f}"
        "   <- near zero means the two axes are nearly independent")

    obs = df[df.observed == 1]
    unobs = df[df.observed == 0]
    log(f"\n  mean transition fraction")
    log(f"    observed alleles      : {obs.ts_fraction.mean():.4f}  (n={len(obs)})")
    log(f"    never-observed alleles: {unobs.ts_fraction.mean():.4f}  (n={len(unobs)})")
    u, p_mw = mannwhitneyu(obs.ts_fraction, unobs.ts_fraction, alternative="greater")
    log(f"    Mann-Whitney (observed > unobserved) p = {p_mw:.4g}")

    a = int(obs.transition_only.sum()); b = len(obs) - a
    c = int(unobs.transition_only.sum()); d = len(unobs) - c
    orr, p_f = fisher_exact([[a, b], [c, d]])
    log(f"\n  reachable ONLY by transition")
    log(f"    observed      : {a}/{len(obs)} ({a/len(obs):.1%})")
    log(f"    never observed: {c}/{len(unobs)} ({c/len(unobs):.1%})")
    log(f"    OR = {orr:.2f}  Fisher exact p = {p_f:.4g}")

    results = [
        {"analysis": "ts_fraction_mannwhitney", "estimate": float(obs.ts_fraction.mean()
                                                                 - unobs.ts_fraction.mean()),
         "p_value": float(p_mw), "n": len(df), "detail": "observed minus unobserved"},
        {"analysis": "transition_only_fisher", "estimate": float(orr),
         "p_value": float(p_f), "n": len(df), "detail": f"{a}/{len(obs)} vs {c}/{len(unobs)}"},
        {"analysis": "ts_fraction_mean_observed", "estimate": float(obs.ts_fraction.mean()),
         "p_value": "", "n": len(obs), "detail": "mean transition fraction, observed"},
        {"analysis": "ts_fraction_mean_unobserved",
         "estimate": float(unobs.ts_fraction.mean()), "p_value": "", "n": len(unobs),
         "detail": "mean transition fraction, never observed"},
        {"analysis": "corr_ts_fraction_esm", "estimate": float(df.ts_fraction.corr(df.esm_score)), "p_value": "",
         "n": len(df), "detail": "Pearson r between transition fraction and ESM score"},
    ]

    # ---------- 2. adjusted for tolerability ------------------------------
    log("\n" + "=" * 70)
    log("ADJUSTED FOR PREDICTED TOLERABILITY  (GEE, sites as clusters)")
    log("=" * 70)
    df["site_id"] = pd.factorize(df["site"])[0]
    df = df.sort_values("site_id").reset_index(drop=True)

    # Every specification is recorded as a row, not only printed, so no figure
    # has to recover these numbers from the spacing of a print statement.
    specs = []

    def fit(cols, label):
        X = sm.add_constant(df[cols].astype(float), has_constant="add")
        m = GEE(df["observed"].astype(float), X, groups=df["site_id"],
                family=Binomial(), cov_struct=Independence()).fit()
        log(f"\n  {label}")
        for cname in cols:
            se = m.bse[cname]
            lo, hi = m.params[cname] - 1.96 * se, m.params[cname] + 1.96 * se
            log(f"    {cname:22} OR={np.exp(m.params[cname]):7.3f}  "
                f"95% CI [{np.exp(lo):.3f}, {np.exp(hi):.3f}]  p={m.pvalues[cname]:.4g}")
            specs.append({"specification": label, "term": cname,
                          "odds_ratio": round(float(np.exp(m.params[cname])), 4),
                          "ci_lo": round(float(np.exp(lo)), 4),
                          "ci_hi": round(float(np.exp(hi)), 4),
                          "p_value": float(m.pvalues[cname]), "n": len(df)})
        return m

    fit(["ts_fraction"], "transition fraction alone")
    fit(["esm_score"], "tolerability alone")
    m_both = fit(["ts_fraction", "esm_score"], "BOTH")
    m_all = fit(["ts_fraction", "esm_score", "frac_codons_reaching"],
                "plus codon supply")
    results.append({"analysis": "gee_ts_adjusted_for_esm",
                    "estimate": float(np.exp(m_both.params["ts_fraction"])),
                    "p_value": float(m_both.pvalues["ts_fraction"]), "n": len(df),
                    "detail": "transition fraction, adjusted for ESM score"})
    results.append({"analysis": "gee_ts_adjusted_full",
                    "estimate": float(np.exp(m_all.params["ts_fraction"])),
                    "p_value": float(m_all.pvalues["ts_fraction"]), "n": len(df),
                    "detail": "adjusted for ESM and codon supply"})

    # ---------- 3. event level -------------------------------------------
    log("\n" + "=" * 70)
    log("EVENT LEVEL  each species' own codon, observed vs its foregone options")
    log("=" * 70)
    ev = [r for r in read_csv(os.path.join(DIR_RESULTS, args.events))
          if r["accessible"] == "True"]
    ev_obs, ev_alt = [], []
    for r in ev:
        codon, wt = r["codon"], r["wild_type_aa"]
        key = (r["gene"], r["site_label"])
        n, ts, tv = codon_path_counts(codon, r["resistance_allele"])
        ev_obs.append(ts / (ts + tv) if (ts + tv) else 0.0)
        for aa in AAS:
            if aa == wt or aa in documented.get(key, set()):
                continue
            n2, t1, t2 = codon_path_counts(codon, aa)
            if n2:
                ev_alt.append(t1 / (t1 + t2) if (t1 + t2) else 0.0)
    log(f"  documented events (accessible): {len(ev_obs)}")
    log(f"  foregone accessible options   : {len(ev_alt)}")
    log(f"    mean ts fraction, chosen  : {np.mean(ev_obs):.4f}")
    log(f"    mean ts fraction, foregone: {np.mean(ev_alt):.4f}")
    u2, p2 = mannwhitneyu(ev_obs, ev_alt, alternative="greater")
    log(f"    Mann-Whitney p = {p2:.4g}")
    results.append({"analysis": "event_level_ts", "estimate":
                    float(np.mean(ev_obs) - np.mean(ev_alt)), "p_value": float(p2),
                    "n": len(ev_obs), "detail": "chosen vs foregone, per-species codon"})

    # ---------- 4. within-site permutation --------------------------------
    log("\n" + "=" * 70)
    log("WITHIN-SITE PERMUTATION  (shuffle which accessible allele is observed)")
    log("=" * 70)
    by_site = defaultdict(list)
    for _, r in df.iterrows():
        by_site[r["site"]].append((int(r["observed"]), float(r["ts_fraction"])))
    usable = {s: v for s, v in by_site.items()
              if 0 < sum(o for o, _ in v) < len(v)}
    log(f"  sites with both observed and unobserved accessible alleles: {len(usable)}")
    obs_stat = float(np.mean([t for s, v in usable.items() for o, t in v if o]))
    ge = 0
    for _ in range(N_PERM):
        vals = []
        for s, v in usable.items():
            k = sum(o for o, _ in v)
            ts_vals = [t for _, t in v]
            vals.extend(rng.sample(ts_vals, k))
        if np.mean(vals) >= obs_stat:
            ge += 1
    p_perm = (ge + 1) / (N_PERM + 1)
    log(f"  observed mean ts fraction of observed alleles : {obs_stat:.4f}")
    log(f"  permutation p ({N_PERM} shuffles)             : {p_perm:.5f}")
    results.append({"analysis": "within_site_permutation_ts", "estimate": obs_stat,
                    "p_value": p_perm, "n": len(usable),
                    "detail": f"{N_PERM} within-site shuffles"})

    p0 = os.path.join(DIR_RESULTS, "mutation_spectrum_substitutions" + TAG + ".csv")
    df.to_csv(p0, index=False)
    p1 = os.path.join(DIR_RESULTS, "mutation_spectrum_results" + TAG + ".csv")
    with open(p1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader(); w.writerows(results)
    p2 = os.path.join(DIR_RESULTS, "mutation_spectrum_specifications" + TAG + ".csv")
    with open(p2, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(specs[0].keys()))
        w.writeheader(); w.writerows(specs)
    log(f"\nwrote {p0}")
    log(f"wrote {p1}")
    log(f"wrote {p2}  ({len(specs)} model terms across "
        f"{len({r['specification'] for r in specs})} specifications)")
    log("\nStep 6.4 complete.")


if __name__ == "__main__":
    main()
