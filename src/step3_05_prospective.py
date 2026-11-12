#!/usr/bin/env python3
"""
step3_05_prospective.py: out-of-sample tests and prospective predictions.

The within-sample tests (step3_03) show that a species' codon is non-randomly matched to
the allele it evolved. They cannot show that the relationship is PREDICTIVE. These
analyses address that, and they are the ones that decide whether the claim is a
description or a forecast.

TEST 4 - temporal replication
-----------------------------
Publication years are parsed from the primary reference of each record. Events are split
at a cutoff year and the accessibility test is re-run independently in each stratum.
Events first reported AFTER the cutoff are out-of-sample with respect to any pattern
visible before it. If the effect is an artefact of the early, heavily-studied literature
(kdr in mosquitoes and houseflies), it should weaken or vanish in the later stratum.

TEST 5 - leave-one-site-out
---------------------------
The accessibility rule is not fitted, so conventional cross-validation does not apply.
The equivalent check is whether the effect is carried by a single site: the test is
re-run with each site removed in turn. A result driven entirely by kdr L1014 would
collapse when L1014 is dropped.

PROSPECTIVE PREDICTIONS
-----------------------
For every species x site in the atlas where the species carries the wild-type residue,
enumerate resistance alleles that are (i) documented at that site in some other species,
and (ii) reachable in ONE nucleotide step from this species' own codon, but (iii) not yet
reported in this species. These are the alleles expected to appear next, and they are
ranked by mutational supply (Ts/Tv-weighted rate). The complementary list - alleles that
are documented elsewhere but NOT reachable in one step - is the negative prediction, and
is the more informative half: those should stay absent.

Outputs
-------
results/temporal_replication.csv
results/leave_one_site_out.csv
results/prospective_predictions.csv
results/prospective_negative_predictions.csv
logs/step3_05_prospective.log
"""
from __future__ import annotations

import csv
import os
import random
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import (  # noqa: E402
    DIR_PROC, DIR_RESULTS, Tee, cluster_bootstrap_ci, ensure_dirs,
)
from step9_01_corrected_statistics import poisson_binomial_sf  # noqa: E402
from codon_engine import (  # noqa: E402
    codon_reaches, codon_reaches_two_step, codon_weight_to,
)

N_PERM = 50_000
SEED = 20260729
YEAR_RE = re.compile(r"(19[5-9]\d|20[0-2]\d)")



def exact_accessibility_probs(events, pools):
    """Per-event accessibility probability under the SAME null the permutation
    samples from: one codon drawn uniformly from that site's wild-type pool.

    P(accessible) is therefore the fraction of the pool that reaches the allele,
    the accessible count is a sum of independent Bernoullis, and its tail is
    exact by Poisson-binomial convolution. The permutation above is a Monte
    Carlo estimate of that same quantity, so it cannot be more accurate and it
    reports the floor 1/(n_perm + 1) whenever there are no exceedances.
    sum(ps) must agree with the permutation's expected value to within Monte
    Carlo error; that agreement is the check that the two constructions are the
    same null.
    """
    ps = []
    for e in events:
        pool = pools.get((e["gene"], e["site_label"]))
        if not pool:
            ps.append(0.0)
            continue
        ps.append(sum(1 for c in pool if codon_reaches(c, e["resistance_allele"]))
                  / len(pool))
    return ps

def perm_test(events, pools, rng, n_perm=N_PERM):
    """Permutation test: resample codons from the site's empirical pool."""
    if not events:
        return None
    obs = sum(1 for e in events if e["accessible"])
    ge, tot = 0, 0
    for _ in range(n_perm):
        s = 0
        for e in events:
            p = pools.get((e["gene"], e["site_label"]))
            if p:
                s += codon_reaches(rng.choice(p), e["resistance_allele"])
        tot += s
        if s >= obs:
            ge += 1
    exp = tot / n_perm
    ps = exact_accessibility_probs(events, pools)
    # Species-cluster interval, so the temporal and leave-one-site-out strata
    # in Supplementary Fig. 2 carry the same uncertainty the headline row does.
    ci_lo, ci_hi = cluster_bootstrap_ci([e["species"] for e in events],
                                        [e["accessible"] for e in events], ps)
    return {"n": len(events), "observed": obs, "expected": exp,
            "ci_lo": ci_lo, "ci_hi": ci_hi,
            "frac_obs": obs / len(events), "frac_exp": exp / len(events),
            "enrichment": obs / exp if exp else float("nan"),
            "p": (ge + 1) / (n_perm + 1),
            "expected_exact": sum(ps),
            "p_exact": poisson_binomial_sf(ps, obs),
            "n_perm": n_perm,
            "p_permutation_floor": 1.0 / (n_perm + 1)}


def main():
    import argparse as _ap
    _p = _ap.ArgumentParser()
    _p.add_argument("--atlas", default="codon_atlas.csv")
    _p.add_argument("--events", default="core_test_events.csv")
    _p.add_argument("--record", default="resistance_record.csv")
    _p.add_argument("--reconciliation", default="species_reconciliation.csv")
    _p.add_argument("--tag", default="")
    _a = _p.parse_args()
    ensure_dirs()
    log = Tee("step3_05_prospective" + _a.tag)
    log("=== Step 3.5  out-of-sample tests and prospective predictions ===\n")
    rng = random.Random(SEED)

    # ---- load atlas ---------------------------------------------------------
    codon, pools, wt_of = {}, defaultdict(list), {}
    order_of, tier_of = {}, {}
    with open(os.path.join(DIR_PROC, _a.atlas)) as fh:
        for r in csv.DictReader(fh):
            order_of[r["species"]] = r["order"]
            tier_of[r["species"]] = r["tier"]
            if r["is_wildtype"] != "True":
                continue
            key = (r["gene"], r["site_label"])
            codon[(r["species"], *key)] = r["codon"]
            wt_of[key] = r["wild_type_aa"]
    for (sp, g, s), c in codon.items():
        pools[(g, s)].append(c)
    log(f"atlas wild-type codon calls: {len(codon)} over {len(pools)} sites")

    with open(os.path.join(DIR_RESULTS, _a.events)) as fh:
        events = list(csv.DictReader(fh))
    for e in events:
        e["accessible"] = e["accessible"] == "True"
        yrs = YEAR_RE.findall(e.get("primary_ref", ""))
        e["year"] = int(yrs[-1]) if yrs else None
    dated = [e for e in events if e["year"]]
    log(f"events: {len(events)};  with a parseable year: {len(dated)}")

    # ---- TEST 4: temporal replication --------------------------------------
    log("\n" + "=" * 70)
    log("TEST 4  temporal replication")
    log("=" * 70)
    trows = []
    if dated:
        yrs = sorted(e["year"] for e in dated)
        log(f"  year range {yrs[0]}-{yrs[-1]}, median {yrs[len(yrs)//2]}")
        for cutoff in (2005, 2008, 2010):
            early = [e for e in dated if e["year"] <= cutoff]
            late = [e for e in dated if e["year"] > cutoff]
            if len(early) < 5 or len(late) < 5:
                log(f"  cutoff {cutoff}: too few events either side "
                    f"({len(early)}/{len(late)}) - skipped")
                continue
            for label, sub in (("<=%d" % cutoff, early), (">%d" % cutoff, late)):
                r = perm_test(sub, pools, rng)
                log(f"  cutoff {cutoff}  {label:>7s}  n={r['n']:3d}  "
                    f"obs {r['observed']:3d} ({r['frac_obs']:5.1%})  "
                    f"null {r['expected']:6.2f} ({r['frac_exp']:5.1%})  "
                    f"enrich {r['enrichment']:4.2f}x  p={r['p']:.5f}")
                trows.append({"cutoff": cutoff, "stratum": label, **r})
        with open(os.path.join(DIR_RESULTS, "temporal_replication" + _a.tag + ".csv"), "w", newline="") as fh:
            if trows:
                w = csv.DictWriter(fh, fieldnames=list(trows[0].keys()))
                w.writeheader(); w.writerows(trows)
        log("  -> the LATER stratum is out-of-sample; the effect should persist there")

    # ---- TEST 5: leave-one-site-out ----------------------------------------
    log("\n" + "=" * 70)
    log("TEST 5  leave-one-site-out")
    log("=" * 70)
    sites = sorted({(e["gene"], e["site_label"]) for e in events},
                   key=lambda k: -sum(1 for e in events if (e["gene"], e["site_label"]) == k))
    full = perm_test(events, pools, rng)
    log(f"  ALL sites          n={full['n']:3d}  enrich {full['enrichment']:4.2f}x  p={full['p']:.6f}")
    lrows = [{"held_out": "none", **full}]
    for key in sites:
        sub = [e for e in events if (e["gene"], e["site_label"]) != key]
        if len(sub) < 5:
            continue
        r = perm_test(sub, pools, rng)
        log(f"  drop {key[1]:8s}      n={r['n']:3d}  enrich {r['enrichment']:4.2f}x  p={r['p']:.6f}")
        lrows.append({"held_out": f"{key[0]}:{key[1]}", **r})
    with open(os.path.join(DIR_RESULTS, "leave_one_site_out" + _a.tag + ".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(lrows[0].keys()))
        w.writeheader(); w.writerows(lrows)
    log("  -> no single site should account for the whole effect")

    # ---- PROSPECTIVE PREDICTIONS -------------------------------------------
    log("\n" + "=" * 70)
    log("PROSPECTIVE PREDICTIONS")
    log("=" * 70)
    with open(os.path.join(DIR_PROC, _a.record)) as fh:
        rec = list(csv.DictReader(fh))
    with open(os.path.join(DIR_PROC, _a.reconciliation)) as fh:
        resolve = {r["record_species"]: r["resolved_species"]
                   for r in csv.DictReader(fh) if r["matched"] == "True"}

    repertoire = defaultdict(set)
    reported = defaultdict(set)          # (species, gene, site) -> alleles already reported
    atlas_sites = set(pools)            # live from the atlas, never a stale flag
    for r in rec:
        if (r["gene"], r["site_label"]) not in atlas_sites:
            continue
        key = (r["gene"], r["site_label"])
        repertoire[key].add(r["resistance_allele"])
        sp = resolve.get(r["species"])
        if sp:
            reported[(sp, *key)].add(r["resistance_allele"])

    pos_rows, neg_rows = [], []
    for (sp, g, s), c in sorted(codon.items()):
        key = (g, s)
        rep = repertoire.get(key)
        if not rep:
            continue
        already = reported.get((sp, g, s), set())
        for allele in sorted(rep):
            if allele in already:
                continue
            reach = codon_reaches(c, allele)
            row = {"species": sp, "order": order_of.get(sp, ""), "tier": tier_of.get(sp, ""),
                   "gene": g, "site_label": s, "wild_type_aa": wt_of.get(key, ""),
                   "codon": c, "candidate_allele": allele,
                   "one_step_accessible": reach,
                   "ts_tv_weighted_rate": round(codon_weight_to(c, allele), 3),
                   "n_2step_paths": codon_reaches_two_step(c, allele) if not reach else 0,
                   "n_species_with_allele": len({resolve.get(r["species"])
                                                 for r in rec
                                                 if (r["gene"], r["site_label"]) == key
                                                 and r["resistance_allele"] == allele
                                                 and resolve.get(r["species"])}),
                   "already_reported_in_species": ";".join(sorted(already))}
            (pos_rows if reach else neg_rows).append(row)

    pos_rows.sort(key=lambda r: (-r["ts_tv_weighted_rate"], -r["n_species_with_allele"], r["species"]))
    neg_rows.sort(key=lambda r: (r["species"], r["gene"], r["site_label"]))

    p1 = os.path.join(DIR_RESULTS, "prospective_predictions" + _a.tag + ".csv")
    with open(p1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pos_rows[0].keys()))
        w.writeheader(); w.writerows(pos_rows)
    p2 = os.path.join(DIR_RESULTS, "prospective_negative_predictions" + _a.tag + ".csv")
    with open(p2, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(neg_rows[0].keys()))
        w.writeheader(); w.writerows(neg_rows)

    log(f"  POSITIVE (accessible, unreported): {len(pos_rows):6d}  -> expected to appear")
    log(f"  NEGATIVE (inaccessible, unreported): {len(neg_rows):6d}  -> expected to stay absent")
    log(f"  species covered: {len({r['species'] for r in pos_rows})}")
    log(f"\n  wrote {p1}")
    log(f"  wrote {p2}")

    log("\n  --- highest-supply predictions at kdr L1014 ---")
    n = 0
    for r in pos_rows:
        if r["site_label"] != "L1014":
            continue
        log(f"    {r['species'][:34]:34s} ({r['order'][:12]:12s}) codon {r['codon']} "
            f"-> L1014{r['candidate_allele']}  rate {r['ts_tv_weighted_rate']:.1f}  "
            f"(allele known in {r['n_species_with_allele']} spp.)")
        n += 1
        if n >= 15:
            break

    by_order = Counter(r["order"] for r in pos_rows)
    log("\n  predictions per order (top 12):")
    for o, k in by_order.most_common(12):
        log(f"    {o[:22]:22s} {k}")
    log.close()


if __name__ == "__main__":
    main()
