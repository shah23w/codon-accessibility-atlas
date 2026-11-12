#!/usr/bin/env python3
"""
step3_04_controls.py: phylogenetic, effort, and mutational-step controls.

The core test (step3_03) is vulnerable to three specific criticisms. Each is addressed
here with an analysis designed to FAIL if the criticism is correct.

CONTROL A - phylogenetic non-independence
-----------------------------------------
Codon identity is phylogenetically clustered, and so is resistance reporting. A naive
association could therefore be an artefact of shared ancestry. The permutation null is
re-run with codons resampled ONLY WITHIN the same taxonomic order, so any signal driven
by between-clade differences in codon usage is removed by construction. If the effect
is phylogenetic, the order-stratified test collapses.

CONTROL B - obligate two-step alleles
-------------------------------------
An orthogonal prediction that uses no per-species information. Some resistance amino
acids are unreachable in one nucleotide step from ANY synonymous codon of the wild-type
residue (e.g. Cys from Leu, Val from Thr): they are obligate two-step substitutions.
Under a mutation-supply model these should be reported in strictly fewer species than
one-step-accessible alleles at the same site. Under a fitness-only model there is no
reason for them to be rarer.

CONTROL C - within-site supply vs prevalence
--------------------------------------------
For each site, correlate how many atlas species can reach an allele (empirical codon
supply, n up to 541) against how many species are reported to carry it. Computed WITHIN
site and combined across sites, because pooling across sites is invalid: sites differ
enormously in study effort and effect size, and a site whose wild-type residue has a
single codon (M918, Met = ATG) has zero accessibility variance by construction.

CONTROL D - detection effort
----------------------------
Species differ by orders of magnitude in how intensively they are studied. Effort is
proxied by the number of distinct resistance records for that species, and the core
test is repeated restricted to well-studied species, where "not observed" is
informative rather than merely unsampled.

Outputs
-------
results/control_phylogenetic.csv
results/control_two_step.csv
results/control_supply_vs_prevalence.csv
results/control_effort.csv
logs/step3_04_controls.log
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import (  # noqa: E402
    DIR_PROC, DIR_RESULTS, Tee, cluster_bootstrap_ci, ensure_dirs,
)
from codon_engine import AA_TO_CODONS, codon_reaches, codon_reaches_two_step  # noqa: E402

N_PERM = 100_000
SEED = 20260729


def spearman(x, y):
    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[s[j + 1]] == v[s[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[s[k]] = avg
            i = j + 1
        return r
    if len(x) < 3:
        return float("nan")
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def main():
    # Same parameterisation as step3_03 so the controls can be run against the
    # original VGSC-only record and the merged record independently. Defaults
    # reproduce the original run exactly.
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", default="resistance_record.csv")
    ap.add_argument("--atlas", default="codon_atlas.csv")
    ap.add_argument("--reconciliation", default="species_reconciliation.csv")
    ap.add_argument("--events", default="core_test_events.csv")
    ap.add_argument("--tag", default="", help="suffix for output files and log")
    args = ap.parse_args()
    TAG = args.tag

    ensure_dirs()
    log = Tee("step3_04_controls" + TAG)
    log("=== Step 3.4  phylogenetic, two-step, supply and effort controls ===\n")
    log(f"record         : {args.record}")
    log(f"events         : {args.events}")
    log(f"reconciliation : {args.reconciliation}")
    log(f"output tag     : '{TAG}'\n")
    rng = random.Random(SEED)

    # ---- load ---------------------------------------------------------------
    order_of, codon, site_wt = {}, {}, {}
    pool_by_site = defaultdict(list)
    pool_by_site_order = defaultdict(list)
    with open(os.path.join(DIR_PROC, args.atlas)) as fh:
        for r in csv.DictReader(fh):
            order_of[r["species"]] = r["order"]
            if r["is_wildtype"] != "True":
                continue
            key = (r["gene"], r["site_label"])
            codon[(r["species"], *key)] = r["codon"]
            site_wt[key] = r["wild_type_aa"]
    for (sp, g, s), c in codon.items():
        pool_by_site[(g, s)].append(c)
        pool_by_site_order[(g, s, order_of.get(sp, ""))].append(c)

    with open(os.path.join(DIR_RESULTS, args.events)) as fh:
        events = list(csv.DictReader(fh))
    for e in events:
        e["accessible"] = e["accessible"] == "True"
        e["order"] = order_of.get(e["species"], "")
    log(f"events: {len(events)};  orders represented: "
        f"{len({e['order'] for e in events})}")

    with open(os.path.join(DIR_PROC, args.record)) as fh:
        rec = list(csv.DictReader(fh))

    # ================= CONTROL A: order-stratified permutation ==============
    log("\n" + "=" * 70)
    log("CONTROL A  order-stratified permutation (phylogenetic control)")
    log("=" * 70)
    usable = [e for e in events
              if len(set(pool_by_site_order.get((e["gene"], e["site_label"], e["order"]), []))) > 0]
    log(f"  events with an order-matched codon pool: {len(usable)}/{len(events)}")
    obs = sum(1 for e in usable if e["accessible"])
    ge, tot = 0, 0
    for _ in range(N_PERM):
        s = 0
        for e in usable:
            p = pool_by_site_order[(e["gene"], e["site_label"], e["order"])]
            s += codon_reaches(rng.choice(p), e["resistance_allele"])
        tot += s
        if s >= obs:
            ge += 1
    exp = tot / N_PERM
    pA = (ge + 1) / (N_PERM + 1)
    log(f"  observed accessible : {obs}/{len(usable)} ({obs/len(usable):.1%})")
    log(f"  null (within-order) : {exp:.2f} ({exp/len(usable):.1%})")
    log(f"  enrichment          : {obs/exp:.2f}x")
    log(f"  p (order-stratified): {pA:.6f}")
    log("  -> if the effect were driven by clade structure this would collapse to ~1x")
    # Interval over whole species, for the same reason as everywhere else: a
    # species with several events is one draw, not several.
    _pA_ci = cluster_bootstrap_ci(
        [e["species"] for e in usable],
        [e["accessible"] for e in usable],
        [sum(1 for c in pool_by_site_order[(e["gene"], e["site_label"],
                                            e["order"])]
             if codon_reaches(c, e["resistance_allele"]))
         / len(pool_by_site_order[(e["gene"], e["site_label"], e["order"])])
         for e in usable])
    log(f"  95% CI (species clusters): {_pA_ci[0]} to {_pA_ci[1]}")
    with open(os.path.join(DIR_RESULTS, "control_phylogenetic" + TAG + ".csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["stratification", "n_events", "n_accessible", "expected_null",
                    "enrichment", "p_permutation", "ci_lo", "ci_hi"])
        w.writerow(["taxonomic_order", len(usable), obs, round(exp, 3),
                    round(obs / exp, 3) if exp else "", pA,
                    _pA_ci[0], _pA_ci[1]])

    # also: does the signal hold within Diptera alone, and outside Diptera?
    for label, sel in (("Diptera only", lambda e: e["order"] == "Diptera"),
                       ("non-Diptera", lambda e: e["order"] != "Diptera")):
        sub = [e for e in events if sel(e)]
        if len(sub) < 4:
            continue
        o = sum(1 for e in sub if e["accessible"])
        ge2, t2 = 0, 0
        for _ in range(20000):
            s = sum(1 for e in sub
                    if codon_reaches(rng.choice(pool_by_site[(e["gene"], e["site_label"])]),
                                     e["resistance_allele"]))
            t2 += s
            if s >= o:
                ge2 += 1
        ex = t2 / 20000
        log(f"    {label:14s} n={len(sub):3d} obs={o:3d} null={ex:5.2f} "
            f"enrich={o/ex if ex else float('nan'):4.2f}x p={(ge2+1)/20001:.5f}")

    # ================= CONTROL B: obligate two-step alleles =================
    log("\n" + "=" * 70)
    log("CONTROL B  obligate two-step alleles are reported in fewer species")
    log("=" * 70)
    sites_in_atlas = set(pool_by_site)   # live from the atlas, never a stale flag
    rows, one_counts, two_counts = [], [], []
    for key in sorted(sites_in_atlas):
        wt = site_wt.get(key)
        if not wt:
            continue
        wt_codons = AA_TO_CODONS[wt]
        by_allele = defaultdict(set)
        for r in rec:
            if (r["gene"], r["site_label"]) == key:
                by_allele[r["resistance_allele"]].add(r["species"])
        for allele, sp in sorted(by_allele.items()):
            n_reach = sum(1 for c in wt_codons if codon_reaches(c, allele))
            obligate2 = n_reach == 0
            rows.append({"gene": key[0], "site_label": key[1], "wild_type_aa": wt,
                         "resistance_allele": allele, "n_species_reported": len(sp),
                         "n_wt_codons": len(wt_codons), "n_codons_reaching": n_reach,
                         "obligate_two_step": obligate2,
                         "min_two_step_paths": min(codon_reaches_two_step(c, allele)
                                                   for c in wt_codons) if obligate2 else 0})
            (two_counts if obligate2 else one_counts).append(len(sp))
            log(f"  {key[1]:7s} {wt}->{allele}  reported in {len(sp):3d} species;  "
                f"{n_reach}/{len(wt_codons)} wild-type codons reach it"
                f"{'   [OBLIGATE 2-STEP]' if obligate2 else ''}")
    if one_counts and two_counts:
        m1 = sum(one_counts) / len(one_counts)
        m2 = sum(two_counts) / len(two_counts)
        # exact Mann-Whitney U via permutation
        allv = one_counts + two_counts
        obs_stat = m1 - m2
        ge = 0
        for _ in range(N_PERM):
            rng.shuffle(allv)
            a = allv[:len(one_counts)]
            b = allv[len(one_counts):]
            if (sum(a) / len(a) - sum(b) / len(b)) >= obs_stat:
                ge += 1
        pB = (ge + 1) / (N_PERM + 1)
        log(f"\n  one-step-accessible alleles : n={len(one_counts):2d}, "
            f"mean {m1:.2f} species (median {sorted(one_counts)[len(one_counts)//2]})")
        log(f"  obligate two-step alleles   : n={len(two_counts):2d}, "
            f"mean {m2:.2f} species (median {sorted(two_counts)[len(two_counts)//2]})")
        log(f"  permutation p (one > two)   : {pB:.5f}")
    with open(os.path.join(DIR_RESULTS, "control_two_step" + TAG + ".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # ================= CONTROL C: supply vs prevalence, within site =========
    log("\n" + "=" * 70)
    log("CONTROL C  within-site codon supply vs species prevalence")
    log("=" * 70)
    crows = []
    rhos = []
    for key in sorted(sites_in_atlas):
        pool = pool_by_site.get(key, [])
        if not pool:
            continue
        by_allele = defaultdict(set)
        for r in rec:
            if (r["gene"], r["site_label"]) == key:
                by_allele[r["resistance_allele"]].add(r["species"])
        if len(by_allele) < 3:
            continue
        xs, ys, labs = [], [], []
        for allele, sp in sorted(by_allele.items()):
            supply = sum(1 for c in pool if codon_reaches(c, allele)) / len(pool)
            xs.append(supply); ys.append(len(sp)); labs.append(allele)
            crows.append({"gene": key[0], "site_label": key[1], "resistance_allele": allele,
                          "atlas_species_pool": len(pool),
                          "frac_pool_reaching": round(supply, 4),
                          "n_species_reported": len(sp)})
        rho = spearman(xs, ys)
        rhos.append((key, rho, len(xs)))
        log(f"  {key[1]:7s} n_alleles={len(xs)}  rho={rho:+.3f}   "
            + "  ".join(f"{a}:{s:.2f}->{n}" for a, s, n in zip(labs, xs, ys)))
    if rhos:
        valid = [r for _k, r, _n in rhos if not math.isnan(r)]
        log(f"\n  sites with >=3 alleles: {len(valid)};  mean within-site rho = "
            f"{sum(valid)/len(valid):+.3f}")
        log("  (pooling ACROSS sites is invalid - see docstring)")
    with open(os.path.join(DIR_RESULTS, "control_supply_vs_prevalence" + TAG + ".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(crows[0].keys()))
        w.writeheader(); w.writerows(crows)

    # ================= CONTROL D: detection effort ==========================
    log("\n" + "=" * 70)
    log("CONTROL D  detection effort")
    log("=" * 70)
    effort = Counter(r["species"] for r in rec)
    with open(os.path.join(DIR_PROC, args.reconciliation)) as fh:
        resolve = {r["record_species"]: r["resolved_species"]
                   for r in csv.DictReader(fh) if r["matched"] == "True"}
    eff_res = Counter()
    for sp, n in effort.items():
        if sp in resolve:
            eff_res[resolve[sp]] += n
    for e in events:
        e["effort"] = eff_res.get(e["species"], 0)
    med_eff = sorted(e["effort"] for e in events)[len(events) // 2]
    log(f"  median records per species: {med_eff}")
    erows = []
    for label, sel in (("well-studied (>= median)", lambda e: e["effort"] >= med_eff),
                       ("less-studied (< median)", lambda e: e["effort"] < med_eff)):
        sub = [e for e in events if sel(e)]
        if len(sub) < 4:
            continue
        o = sum(1 for e in sub if e["accessible"])
        ge2, t2 = 0, 0
        for _ in range(20000):
            s = sum(1 for e in sub
                    if codon_reaches(rng.choice(pool_by_site[(e["gene"], e["site_label"])]),
                                     e["resistance_allele"]))
            t2 += s
            if s >= o:
                ge2 += 1
        ex = t2 / 20000
        p = (ge2 + 1) / 20001
        log(f"  {label:26s} n={len(sub):3d} obs={o:3d} ({o/len(sub):5.1%})  "
            f"null={ex:5.2f}  enrich={o/ex if ex else float('nan'):4.2f}x  p={p:.5f}")
        _ci = cluster_bootstrap_ci(
            [e["species"] for e in sub],
            [e["accessible"] for e in sub],
            [sum(1 for c in pool_by_site[(e["gene"], e["site_label"])]
                 if codon_reaches(c, e["resistance_allele"]))
             / len(pool_by_site[(e["gene"], e["site_label"])]) for e in sub])
        erows.append({"stratum": label, "n_events": len(sub), "n_accessible": o,
                      "ci_lo": _ci[0], "ci_hi": _ci[1],
                      "expected_null": round(ex, 3),
                      "enrichment": round(o / ex, 3) if ex else "", "p_permutation": p})
    log("  -> the effect should NOT be confined to well-studied species")
    with open(os.path.join(DIR_RESULTS, "control_effort" + TAG + ".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(erows[0].keys()))
        w.writeheader(); w.writerows(erows)

    log.close()


if __name__ == "__main__":
    main()
