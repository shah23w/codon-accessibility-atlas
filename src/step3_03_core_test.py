#!/usr/bin/env python3
"""
step3_03_core_test.py: does a species' OWN codon predict which allele it evolved?

This is the study's confirmatory test. Three designs are run, deliberately, because
each is vulnerable to a different criticism and they fail independently.

TEST 1 - Codon-permutation test (primary)
-----------------------------------------
For every documented event (species i evolved allele a at site S), score 1 if a is
reachable in one nucleotide step from that species' own codon c_i.
Statistic = number of accessible events.

Null: the observed pairing of codons to species is broken by resampling each event's
codon from the EMPIRICAL codon distribution at that site across all wild-type species in
the atlas (n up to 541), holding the observed allele fixed.

Why this null and not the proposal's original one: drawing the allele from the site's
"known repertoire" is circular, because that repertoire is itself compiled from the
species that evolved those alleles, and is therefore already accessibility-enriched.
Permuting CODONS instead leaves the allele distribution untouched and asks only whether
the codon a species happens to carry is non-randomly matched to the allele it evolved.
The atlas supplies the codon distribution, so the null is empirical, not assumed.

TEST 2 - Within-species conditional test
----------------------------------------
Restricted to species carrying a codon that can reach SOME but not ALL alleles
documented at a site. These species are individually informative: they discriminate the
hypothesis from any site-level or species-level effect, because the prediction differs
between alleles WITHIN the same species and the same site. Detection effort, insecticide
use, and phylogeny are all held constant by construction.

TEST 3 - Reciprocal-accessibility test
--------------------------------------
The strongest qualitative case. For allele pairs (a1, a2) at one site, find species where
a1 is accessible and a2 is not, and species where the reverse holds. If mutational
supply drives the pattern, each species should show only its accessible allele.
The Musca CTT (reaches His, not Ser) versus Anopheles TTA (reaches Ser, not His) case is
the worked example; this generalises it across the whole atlas.

Outputs
-------
results/core_test_permutation.csv     per-site and pooled permutation results
results/core_test_events.csv          every event with its accessibility call
results/core_test_reciprocal.csv      reciprocal-accessibility contrasts
results/core_test_conditional.csv     within-species discriminating events
logs/step3_03_core_test.log
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import DIR_PROC, DIR_RESULTS, Tee, ensure_dirs  # noqa: E402
from codon_engine import codon_reaches, codon_reaches_two_step  # noqa: E402

N_PERM = 100_000
SEED = 20260729          # fixed for reproducibility; permutation only, never data generation


def main():
    # Optional inputs so the test can be run on the original VGSC-only record and on
    # the merged record without either overwriting the other. Defaults reproduce the
    # original run exactly.
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", default="resistance_record.csv")
    ap.add_argument("--atlas", default="codon_atlas.csv")
    ap.add_argument("--reconciliation", default="species_reconciliation.csv")
    ap.add_argument("--tag", default="", help="suffix for output files and log")
    args = ap.parse_args()
    TAG = args.tag

    ensure_dirs()
    log = Tee("step3_03_core_test" + TAG)
    log("=== Step 3.3  core test: does a species' own codon predict its allele? ===\n")
    log(f"record         : {args.record}")
    log(f"reconciliation : {args.reconciliation}")
    log(f"output tag     : '{TAG}'\n")

    # ---- load ---------------------------------------------------------------
    with open(os.path.join(DIR_PROC, args.reconciliation)) as fh:
        resolve = {r["record_species"]: r["resolved_species"]
                   for r in csv.DictReader(fh) if r["matched"] == "True"}
    log(f"reconciled record species: {len(resolve)}")

    # codon per species per site, restricted to species carrying the wild type
    codon = {}
    site_codon_pool = defaultdict(list)
    site_wt = {}
    with open(os.path.join(DIR_PROC, args.atlas)) as fh:
        for r in csv.DictReader(fh):
            if r["is_wildtype"] != "True":
                continue
            key = (r["gene"], r["site_label"])
            codon[(r["species"], *key)] = r["codon"]
            site_wt[key] = r["wild_type_aa"]
    for (sp, g, s), c in codon.items():
        site_codon_pool[(g, s)].append(c)
    log(f"atlas wild-type codon calls: {len(codon)}")

    with open(os.path.join(DIR_PROC, args.record)) as fh:
        rec = list(csv.DictReader(fh))

    # Site membership is determined LIVE from the atlas, not from a precomputed flag in
    # the record. The flag written by step3_01 reflects whatever reference set existed
    # when the record was parsed, and goes stale as soon as the site set is expanded
    # (step1_05) - which silently discarded 49 events in an earlier run. Deriving it here
    # removes that circular dependency entirely.
    atlas_sites = {(g, s) for (_sp, g, s) in codon}

    # ---- assemble events ----------------------------------------------------
    events, skipped = [], Counter()
    for r in rec:
        if (r["gene"], r["site_label"]) not in atlas_sites:
            skipped["site not in atlas"] += 1
            continue
        sp = resolve.get(r["species"])
        if not sp:
            skipped["species not in atlas"] += 1
            continue
        key = (r["gene"], r["site_label"])
        c = codon.get((sp, *key))
        if not c:
            skipped["species not wild-type / no codon at site"] += 1
            continue
        allele = r["resistance_allele"]
        events.append({
            "species": sp, "record_species": r["species"], "gene": r["gene"],
            "site_label": r["site_label"], "wild_type_aa": r["wild_type_aa"],
            "codon": c, "resistance_allele": allele,
            "accessible": codon_reaches(c, allele),
            "n_2step": codon_reaches_two_step(c, allele),
            "in_combination": r["is_part_of_combination"],
            "primary_ref": r["primary_ref"],
        })
    log(f"\nusable events: {len(events)}")
    for k, v in skipped.most_common():
        log(f"  skipped - {k}: {v}")

    if not events:
        raise RuntimeError("no usable events - aborting")

    n_acc = sum(1 for e in events if e["accessible"])
    log(f"\nobserved: {n_acc}/{len(events)} events accessible in one step "
        f"({n_acc/len(events):.1%})")

    # ---- TEST 1: codon permutation -----------------------------------------
    log("\n" + "=" * 70)
    log("TEST 1  codon-permutation test")
    log("=" * 70)
    rng = random.Random(SEED)
    by_site = defaultdict(list)
    for e in events:
        by_site[(e["gene"], e["site_label"])].append(e)

    rows = []
    for key in sorted(by_site, key=lambda k: -len(by_site[k])):
        evs = by_site[key]
        pool = site_codon_pool.get(key, [])
        if not pool:
            continue
        obs = sum(1 for e in evs if e["accessible"])
        ge = 0
        null_tot = 0
        for _ in range(N_PERM):
            stat = sum(1 for e in evs
                       if codon_reaches(rng.choice(pool), e["resistance_allele"]))
            null_tot += stat
            if stat >= obs:
                ge += 1
        exp = null_tot / N_PERM
        p = (ge + 1) / (N_PERM + 1)
        rows.append({"gene": key[0], "site_label": key[1], "n_events": len(evs),
                     "n_accessible": obs, "frac_accessible": round(obs / len(evs), 4),
                     "expected_null": round(exp, 3),
                     "enrichment": round(obs / exp, 3) if exp else "",
                     "p_permutation": p, "codon_pool_n": len(pool),
                     "n_distinct_codons": len(set(pool))})
        log(f"  {key[0]:6s} {key[1]:7s} n={len(evs):3d}  accessible {obs:3d} "
            f"({obs/len(evs):5.1%})  null {exp:6.2f}  "
            f"enrich {obs/exp if exp else float('nan'):5.2f}x  p={p:.5f}")

    # pooled across sites
    obs_all = n_acc
    ge = 0
    null_tot = 0
    for _ in range(N_PERM):
        stat = 0
        for e in events:
            pool = site_codon_pool[(e["gene"], e["site_label"])]
            stat += codon_reaches(rng.choice(pool), e["resistance_allele"])
        null_tot += stat
        if stat >= obs_all:
            ge += 1
    exp_all = null_tot / N_PERM
    p_all = (ge + 1) / (N_PERM + 1)
    rows.append({"gene": "ALL", "site_label": "POOLED", "n_events": len(events),
                 "n_accessible": obs_all, "frac_accessible": round(obs_all / len(events), 4),
                 "expected_null": round(exp_all, 3),
                 "enrichment": round(obs_all / exp_all, 3) if exp_all else "",
                 "p_permutation": p_all, "codon_pool_n": "", "n_distinct_codons": ""})
    log(f"\n  POOLED  n={len(events)}  accessible {obs_all} ({obs_all/len(events):.1%})  "
        f"null {exp_all:.2f} ({exp_all/len(events):.1%})  "
        f"enrichment {obs_all/exp_all:.2f}x  p={p_all:.6f}")

    # ---- TEST 1b: informative subset ---------------------------------------
    # A site x allele combination carries information ONLY if the site's codon pool
    # contains both codons that reach the allele and codons that do not. Where every
    # wild-type codon reaches the allele (e.g. any substitution from Met, which has the
    # single codon ATG), accessibility has zero variance and the event is uninformative
    # by construction: it can neither support nor refute the hypothesis, but it does
    # inflate the denominator and dilute the apparent effect size.
    #
    # This stratification is defined by the DESIGN, not by the outcome, and both strata
    # are reported.
    log("\n" + "=" * 70)
    log("TEST 1b  informative vs uninformative site x allele combinations")
    log("=" * 70)
    informative, uninformative = [], []
    for e in events:
        pool = site_codon_pool.get((e["gene"], e["site_label"]), [])
        frac = (sum(1 for c in pool if codon_reaches(c, e["resistance_allele"])) / len(pool)
                if pool else 0.0)
        e["site_supply_fraction"] = round(frac, 4)
        (informative if 0.0 < frac < 1.0 else uninformative).append(e)
    log(f"  informative   (0 < supply < 1): {len(informative)} events")
    log(f"  uninformative (supply = 0 or 1): {len(uninformative)} events")

    for label, sub in (("informative", informative), ("uninformative", uninformative)):
        if len(sub) < 3:
            continue
        obs = sum(1 for e in sub if e["accessible"])
        ge, tot = 0, 0
        for _ in range(N_PERM):
            s = sum(1 for e in sub
                    if codon_reaches(rng.choice(site_codon_pool[(e["gene"], e["site_label"])]),
                                     e["resistance_allele"]))
            tot += s
            if s >= obs:
                ge += 1
        exp = tot / N_PERM
        log(f"  {label:14s} n={len(sub):3d}  obs {obs:3d} ({obs/len(sub):5.1%})  "
            f"null {exp:6.2f} ({exp/len(sub):5.1%})  "
            f"enrich {obs/exp if exp else float('nan'):4.2f}x  p={(ge+1)/(N_PERM+1):.6f}")
        rows.append({"gene": "ALL", "site_label": f"POOLED_{label}", "n_events": len(sub),
                     "n_accessible": obs, "frac_accessible": round(obs / len(sub), 4),
                     "expected_null": round(exp, 3),
                     "enrichment": round(obs / exp, 3) if exp else "",
                     "p_permutation": (ge + 1) / (N_PERM + 1),
                     "codon_pool_n": "", "n_distinct_codons": ""})

    with open(os.path.join(DIR_RESULTS, "core_test_permutation" + TAG + ".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    with open(os.path.join(DIR_RESULTS, "core_test_events" + TAG + ".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(events[0].keys()))
        w.writeheader(); w.writerows(sorted(events, key=lambda e: (e["site_label"], e["species"])))

    # ---- TEST 2: within-species conditional --------------------------------
    log("\n" + "=" * 70)
    log("TEST 2  within-species conditional (codon discriminates BETWEEN alleles)")
    log("=" * 70)
    repertoire = defaultdict(set)
    for r in rec:
        if (r["gene"], r["site_label"]) in atlas_sites:
            repertoire[(r["gene"], r["site_label"])].add(r["resistance_allele"])

    cond, n_disc, n_hit = [], 0, 0
    for e in events:
        key = (e["gene"], e["site_label"])
        rep = repertoire[key]
        reach = {a for a in rep if codon_reaches(e["codon"], a)}
        if not reach or reach == rep:
            continue                      # codon does not discriminate at this site
        n_disc += 1
        hit = e["resistance_allele"] in reach
        n_hit += hit
        cond.append({**{k: e[k] for k in ("species", "gene", "site_label", "codon",
                                          "resistance_allele")},
                     "site_repertoire": ",".join(sorted(rep)),
                     "accessible_subset": ",".join(sorted(reach)),
                     "n_repertoire": len(rep), "n_accessible": len(reach),
                     "prob_by_chance": round(len(reach) / len(rep), 4),
                     "observed_allele_accessible": hit})
    if n_disc:
        exp = sum(c["prob_by_chance"] for c in cond)
        log(f"  discriminating events: {n_disc}")
        log(f"  observed allele accessible: {n_hit}/{n_disc} ({n_hit/n_disc:.1%})")
        log(f"  expected by chance        : {exp:.2f}/{n_disc} ({exp/n_disc:.1%})")
        # exact binomial-style permutation using per-event probabilities
        ge = 0
        for _ in range(N_PERM):
            s = sum(1 for c in cond if rng.random() < c["prob_by_chance"])
            if s >= n_hit:
                ge += 1
        p2 = (ge + 1) / (N_PERM + 1)
        log(f"  p (Poisson-binomial permutation) = {p2:.6f}")
        with open(os.path.join(DIR_RESULTS, "core_test_conditional" + TAG + ".csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(cond[0].keys()))
            w.writeheader(); w.writerows(cond)
    else:
        log("  no discriminating events")

    # ---- TEST 3: reciprocal accessibility -----------------------------------
    log("\n" + "=" * 70)
    log("TEST 3  reciprocal accessibility")
    log("=" * 70)
    recip = []
    obs_by_sp = defaultdict(set)
    for e in events:
        obs_by_sp[(e["species"], e["gene"], e["site_label"])].add(e["resistance_allele"])

    for key, rep in sorted(repertoire.items()):
        alleles = sorted(rep)
        for i in range(len(alleles)):
            for j in range(i + 1, len(alleles)):
                a1, a2 = alleles[i], alleles[j]
                grp1, grp2 = [], []
                for (sp, g, s), got in obs_by_sp.items():
                    if (g, s) != key:
                        continue
                    c = codon.get((sp, g, s))
                    if not c:
                        continue
                    r1, r2 = codon_reaches(c, a1), codon_reaches(c, a2)
                    if r1 and not r2:
                        grp1.append((sp, c, sorted(got)))
                    elif r2 and not r1:
                        grp2.append((sp, c, sorted(got)))
                if not grp1 or not grp2:
                    continue
                # among species where only a1 is accessible, how many show a1 vs a2?
                g1_a1 = sum(1 for _s, _c, got in grp1 if a1 in got)
                g1_a2 = sum(1 for _s, _c, got in grp1 if a2 in got)
                g2_a1 = sum(1 for _s, _c, got in grp2 if a1 in got)
                g2_a2 = sum(1 for _s, _c, got in grp2 if a2 in got)
                recip.append({"gene": key[0], "site_label": key[1],
                              "allele_1": a1, "allele_2": a2,
                              "n_species_only_a1_accessible": len(grp1),
                              "n_species_only_a2_accessible": len(grp2),
                              "a1_observed_in_a1group": g1_a1, "a2_observed_in_a1group": g1_a2,
                              "a1_observed_in_a2group": g2_a1, "a2_observed_in_a2group": g2_a2,
                              "concordant": g1_a1 + g2_a2, "discordant": g1_a2 + g2_a1,
                              "species_a1group": ";".join(f"{s}({c})" for s, c, _ in grp1),
                              "species_a2group": ";".join(f"{s}({c})" for s, c, _ in grp2)})
                log(f"  {key[1]:7s} {a1} vs {a2}:  "
                    f"only-{a1}-accessible n={len(grp1)} -> {a1}:{g1_a1} {a2}:{g1_a2}   |   "
                    f"only-{a2}-accessible n={len(grp2)} -> {a1}:{g2_a1} {a2}:{g2_a2}")
    if recip:
        tot_c = sum(r["concordant"] for r in recip)
        tot_d = sum(r["discordant"] for r in recip)
        log(f"\n  TOTAL concordant (allele appears where it IS accessible)   : {tot_c}")
        log(f"  TOTAL discordant (allele appears where it is NOT accessible): {tot_d}")
        with open(os.path.join(DIR_RESULTS, "core_test_reciprocal" + TAG + ".csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(recip[0].keys()))
            w.writeheader(); w.writerows(recip)
    else:
        log("  no reciprocal contrasts available")

    log.close()


if __name__ == "__main__":
    main()
