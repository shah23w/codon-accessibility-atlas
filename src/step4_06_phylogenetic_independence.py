#!/usr/bin/env python3
"""
step4_06_phylogenetic_independence.py: how much of the result survives if the events
are not treated as independent?

The objection
-------------
The core test treats 143 species x site x allele events as independent observations.
They are not. L1014F is documented in many species; if some of those share the allele by
descent rather than by repeated origin, the effective sample size is smaller than 143 and
every p-value is optimistic. CONTROL A (order-stratified permutation) removes
between-clade differences in CODON composition, but it does not address non-independence
of the EVENTS themselves.

The test
--------
Collapse events to one per (gene, site, allele) within each taxonomic rank, so that a
resistance allele shared across a genus, family or order counts once rather than many
times. This is deliberately conservative: kdr alleles are known to have arisen
repeatedly and independently, so collapsing genuinely independent origins throws away
real signal. If the effect survives collapsing to ORDER, the most brutal version, where
all Diptera carrying L1014F count as a single event, it cannot be an artefact of
pseudoreplication.

Under each collapsing scheme the same permutation null is used as step3_03: resample
each event's codon from the empirical distribution of codons at that site across the
atlas, holding the observed allele fixed.

Outputs
-------
results/phylogenetic_independence.csv
logs/step4_06_phylogenetic_independence.log
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import (  # noqa: E402
    DIR_PROC, DIR_RESULTS, Tee, cluster_bootstrap_ci, ensure_dirs,
)
from codon_engine import codon_reaches  # noqa: E402
from step9_01_corrected_statistics import poisson_binomial_sf  # noqa: E402

N_PERM = 20000
SEED = 20260730


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def exact_accessibility_probs(events, pool):
    """Per-event accessibility probability under the SAME null the permutation
    samples from: one codon drawn uniformly from that site's pool.

    P(accessible) is therefore the fraction of the pool that reaches the allele,
    the accessible count is a sum of independent Bernoullis, and its tail is
    exact by Poisson-binomial convolution. The permutation is a Monte Carlo
    estimate of that same quantity and reports the floor 1/(n_perm + 1) when
    there are no exceedances, which is what every row of this file did at
    20,000 permutations. sum(ps) must agree with the permutation's expected
    value to within Monte Carlo error.
    """
    ps = []
    for a, _c, site in events:
        codons = pool.get(site)
        if not codons:
            ps.append(0.0)
            continue
        ps.append(sum(1 for c2 in codons if codon_reaches(c2, a)) / len(codons))
    return ps


def permutation_test(events, pool, rng, n_perm=N_PERM):
    """events: list of (allele, codon). pool: {site_key: [codons]}"""
    obs = sum(1 for a, c, _ in events if codon_reaches(c, a))
    ge = 0
    exp_total = 0.0
    for _ in range(n_perm):
        k = 0
        for a, _c, site in events:
            c2 = rng.choice(pool[site])
            if codon_reaches(c2, a):
                k += 1
        exp_total += k
        if k >= obs:
            ge += 1
    exp = exp_total / n_perm
    ps = exact_accessibility_probs(events, pool)
    return (obs, exp, (ge + 1) / (n_perm + 1),
            sum(ps), poisson_binomial_sf(ps, obs), 1.0 / (n_perm + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default="codon_atlas_extended.csv")
    ap.add_argument("--events", default="core_test_events_primary.csv")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    ensure_dirs()
    log = Tee("step4_06_phylogenetic_independence" + args.tag)
    log("=== Step 4.6  are the events phylogenetically independent? ===\n")
    rng = random.Random(SEED)

    events = read_csv(os.path.join(DIR_RESULTS, args.events))
    log(f"events: {len(events)}")

    # taxonomy for each species
    panel = {r["species"]: r for r in
             read_csv(os.path.join(DIR_PROC, "species_panel.csv"))}
    atlas_tax = {}
    pool = defaultdict(list)
    for r in read_csv(os.path.join(DIR_PROC, args.atlas)):
        if r["is_wildtype"] == "True":
            pool[(r["gene"], r["site_label"])].append(r["codon"])
            atlas_tax[r["species"]] = r["order"]

    def rank_of(sp, rank):
        p = panel.get(sp)
        if p and p.get(rank):
            return p[rank]
        if rank == "order":
            return atlas_tax.get(sp, "")
        if rank == "genus":
            return sp.split()[0]
        return ""

    missing = Counter()
    for r in events:
        for rank in ("genus", "family", "order"):
            if not rank_of(r["species"], rank):
                missing[rank] += 1
    log(f"species lacking rank information: {dict(missing)}")

    schemes = [("none (every event)", None), ("genus", "genus"),
               ("family", "family"), ("order", "order")]
    out = []
    for label, rank in schemes:
        # The cluster each collapsed event belongs to is carried alongside,
        # because the interval below must resample the unit the stratum is
        # defined on: the species when nothing is collapsed, otherwise the
        # genus, family or order.
        seen, collapsed, labels = set(), [], []
        for r in events:
            if rank is None:
                key = (r["species"], r["gene"], r["site_label"],
                       r["resistance_allele"])
            else:
                grp = rank_of(r["species"], rank)
                if not grp:
                    grp = f"__unranked__{r['species']}"
                key = (grp, r["gene"], r["site_label"], r["resistance_allele"])
            if key in seen:
                continue
            seen.add(key)
            collapsed.append((r["resistance_allele"], r["codon"],
                              (r["gene"], r["site_label"])))
            labels.append(key[0])
        pairs = [(e, lab) for e, lab in zip(collapsed, labels)
                 if pool.get(e[2])]
        usable = [e for e, _ in pairs]
        obs, exp, p, exp_exact, p_exact, p_floor = permutation_test(
            usable, pool, rng)
        enr = obs / exp if exp else float("nan")
        log(f"\n  collapse by {label}")
        log(f"    independent events : {len(usable)}  (was {len(events)})")
        log(f"    accessible         : {obs} ({obs/len(usable):.1%})")
        log(f"    permutation null   : {exp:.2f} ({exp/len(usable):.1%})")
        log(f"    enrichment         : {enr:.2f}x")
        log(f"    p                  : {p:.5f}")
        ci_lo, ci_hi = cluster_bootstrap_ci(
            [lab for _, lab in pairs],
            [codon_reaches(e[1], e[0]) for e in usable],
            exact_accessibility_probs(usable, pool))
        log(f"    95% CI             : {ci_lo} to {ci_hi}")
        out.append({"collapse_rank": label, "n_events": len(usable),
                    "n_accessible": obs, "ci_lo": ci_lo, "ci_hi": ci_hi,
                    "null_expected": round(exp, 3),
                    "enrichment": round(enr, 4), "p_value": p,
                    "expected_exact": round(exp_exact, 4),
                    "p_exact": p_exact,
                    "n_perm": N_PERM,
                    "p_permutation_floor": p_floor})

    log("\n  reading: if the effect were an artefact of counting the same")
    log("  evolutionary origin many times, enrichment would fall toward 1.0 as")
    log("  events are collapsed to coarser ranks.")

    p1 = os.path.join(DIR_RESULTS, "phylogenetic_independence" + args.tag + ".csv")
    with open(p1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    log(f"\nwrote {p1}")
    log("\nStep 4.6 complete.")


if __name__ == "__main__":
    main()
