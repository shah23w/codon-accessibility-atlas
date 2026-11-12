#!/usr/bin/env python3
"""
step4_07_matched_null_contrast.py: score the rapid-adaptation and deep-time classes
under the SAME nulls, so the timescale contrast can be read honestly.

Why this step exists
--------------------
The two classes were never compared under a common null. The rapid class (step3_03) uses
a CODON-randomisation null: hold the observed residue fixed, draw the species' codon from
the pool of wild-type codons at that site. The deep-time class (step4_01) uses a
RESIDUE-randomisation null: hold the codon fixed, draw the derived residue from the viable
set at that site. Those are opposite randomisations of different objects, and a ratio
computed under one is not comparable with a ratio computed under the other.

Worse, step4_01 computes the residue null two ways - weighted by panel frequency (1.010x)
and uniform (1.166x) - and only the first was reported. Its own docstring anticipated the
problem: frequency weighting partly encodes accessibility, because residues that are easy
to reach become common in the panel. The uniform value sits inside the rapid class's
bootstrap CI, so the reported contrast may be an artefact of the null choice.

What this does
--------------
Both classes are reduced to a common representation - (species, gene, site, own codon(s),
derived residue) - and scored under all three nulls:

  codon_pool        draw the codon from the wild-type codons at that site (the rapid
                    class's own null, now applied to both)
  residue_uniform   draw the derived residue uniformly from the viable set at that site
  residue_weighted  the same, weighted by how often each residue occurs in the panel

Only same-null (matched) comparisons between classes are interpretable. Every cell gets an
exact Poisson-binomial tail and a species-cluster bootstrap CI, and everything is repeated
on the subset of species and sites the two classes share.

Output: results/matched_null_contrast.csv
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
from codon_engine import codon_reaches, translate_codon  # noqa: E402
from step9_01_corrected_statistics import poisson_binomial_sf  # noqa: E402

N_BOOT = 2000
SEED = 20260731


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


class Obs:
    """One observation, in a form both classes share."""

    __slots__ = ("species", "gene", "site", "codons", "derived", "ancestral_aas",
                 "accessible")

    def __init__(self, species, gene, site, codons, derived, ancestral_aas):
        self.species = species
        self.gene = gene
        self.site = site
        self.codons = tuple(codons)
        self.derived = derived
        self.ancestral_aas = set(ancestral_aas)
        self.accessible = any(codon_reaches(c, derived) for c in self.codons)

    @property
    def key(self):
        return (self.gene, self.site)


# ------------------------------------------------------------------ the nulls
def null_codon_pool(o, codon_pool, viable):
    """Randomise the codon; hold the observed residue fixed."""
    pool = codon_pool.get(o.key)
    if not pool:
        return None
    return sum(1 for c in pool if codon_reaches(c, o.derived)) / len(pool)


def _residue_pool(o, viable):
    """Residues this observation could plausibly have moved TO at this site."""
    v = viable.get(o.key)
    if not v:
        return None
    pool = {a: n for a, n in v.items() if a not in o.ancestral_aas}
    return pool or None


def null_residue_uniform(o, codon_pool, viable):
    pool = _residue_pool(o, viable)
    if not pool:
        return None
    hits = sum(1 for a in pool if any(codon_reaches(c, a) for c in o.codons))
    return hits / len(pool)


def null_residue_weighted(o, codon_pool, viable):
    pool = _residue_pool(o, viable)
    if not pool:
        return None
    tot = sum(pool.values())
    hits = sum(n for a, n in pool.items()
               if any(codon_reaches(c, a) for c in o.codons))
    return hits / tot if tot else None


NULLS = [("codon_pool", null_codon_pool),
         ("residue_uniform", null_residue_uniform),
         ("residue_weighted", null_residue_weighted)]


# ------------------------------------------------------------------ estimation
def estimate(obs, probs, rng):
    """Observed/expected with an exact tail and a species-cluster bootstrap CI."""
    n = len(obs)
    o = sum(1 for x in obs if x.accessible)
    e = sum(probs)
    ratio = o / e if e else float("nan")
    p = poisson_binomial_sf(probs, o) if 0 < e < n else float("nan")

    by_sp = defaultdict(list)
    for x, pr in zip(obs, probs):
        by_sp[x.species].append((x.accessible, pr))
    species = list(by_sp)
    boots = []
    for _ in range(N_BOOT):
        oo = ee = 0.0
        for _ in range(len(species)):
            for acc, pr in by_sp[species[rng.randrange(len(species))]]:
                oo += acc
                ee += pr
        if ee:
            boots.append(oo / ee)
    boots.sort()
    lo = boots[int(0.025 * len(boots))] if boots else float("nan")
    hi = boots[int(0.975 * len(boots)) - 1] if boots else float("nan")
    return dict(n=n, observed=o, expected=round(e, 3), ratio=round(ratio, 4),
                p_exact=p, ci_lo=round(lo, 4), ci_hi=round(hi, 4),
                n_species=len(species))


def validate(deep, panel_only, log):
    """Reproduce step4_01's pooled deep-time nulls exactly, or fail loudly.

    step4_01 reports 489 usable substitutions, 358 accessible, expected 354.601
    (frequency-weighted) and 307.037 (uniform). If this script cannot reproduce those,
    its matched-null numbers cannot be trusted either.
    """
    log("-" * 78)
    log("VALIDATION  reproduce step4_01's deep-time nulls with the panel-only universe")
    log("-" * 78)
    n = obs = 0
    exp_w = exp_u = 0.0
    for o in deep:
        pool = {a: c for a, c in panel_only.get(o.key, {}).items()
                if a not in o.ancestral_aas}
        if not pool:
            continue
        n += 1
        obs += o.accessible
        hit = [a for a in pool if any(codon_reaches(c, a) for c in o.codons)]
        exp_u += len(hit) / len(pool)
        exp_w += sum(pool[a] for a in hit) / sum(pool.values())
    log(f"  usable substitutions : {n:>8}   (step4_01: 489)")
    log(f"  accessible           : {obs:>8}   (step4_01: 358)")
    log(f"  weighted expected    : {exp_w:>8.3f}   (step4_01: 354.601)  "
        f"ratio {obs / exp_w:.3f} (1.010)")
    log(f"  uniform  expected    : {exp_u:>8.3f}   (step4_01: 307.037)  "
        f"ratio {obs / exp_u:.3f} (1.166)")
    ok = (n == 489 and obs == 358
          and abs(exp_w - 354.601) < 0.05 and abs(exp_u - 307.037) < 0.05)
    log(f"\n  REPRODUCED: {ok}")
    if not ok:
        raise SystemExit("validation FAILED - the matched-null numbers cannot be trusted")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default="codon_atlas.csv")
    ap.add_argument("--events", default="core_test_events_primary.csv")
    ap.add_argument("--substitutions", default="ancestral_substitutions.csv")
    ap.add_argument("--validate", action="store_true",
                    help="reproduce step4_01's two deep-time nulls and assert they match")
    args = ap.parse_args()

    ensure_dirs()
    log = Tee("step4_07_matched_null_contrast")
    log("=== Step 4.7  matched-null contrast: rapid adaptation vs deep time ===\n")
    rng = random.Random(SEED)

    # ALL species, deduped to one codon per species per site - NOT wild-type only.
    # step4_01 does the same, and for a reason: a deep-time substitution is by
    # definition a species that does NOT carry the reference residue, so filtering to
    # wild-type erases exactly the variation the deep-time class is built from. (An
    # earlier version of this script filtered to wild-type and produced a viable set of
    # one residue at sites where step4_01 recorded three.)
    seen = set()
    codon_pool = defaultdict(list)
    n_wt_calls = 0
    for r in read_csv(os.path.join(DIR_PROC, args.atlas)):
        k = (r["gene"], r["site_label"], r["species"])
        if k in seen:
            continue
        seen.add(k)
        if r["is_wildtype"] == "True":
            n_wt_calls += 1
        codon_pool[(r["gene"], r["site_label"])].append(r["codon"])
    n_all_calls = sum(len(v) for v in codon_pool.values())
    # The codon-randomised ratio reported here for the rapid class is NOT the
    # pooled enrichment in Tables S1 and S2, and the difference is entirely the
    # pool it draws from: this script draws from every species' codon at the
    # site, step9_01 draws only from the wild-type codons. The wider pool
    # contains the derived residues a deep-time substitution is defined by, so
    # restricting it to wild-type would erase the deep-time class. Both numbers
    # are written out, so that the Supplementary Information states the reason
    # from a file rather than asserting it.
    _pp = os.path.join(DIR_RESULTS, "matched_null_codon_pool.csv")
    with open(_pp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "value", "definition"])
        w.writerow(["pool_all_species", n_all_calls,
                    "codon calls in the pool this script randomises over: one per "
                    "(species, gene, site), wild-type or not"])
        w.writerow(["pool_wild_type_only", n_wt_calls,
                    "the subset of that pool carrying the wild-type residue, which "
                    "is the pool step9_01 uses for the pooled enrichment in "
                    "Tables S1 and S2"])
        w.writerow(["pool_difference", n_all_calls - n_wt_calls,
                    "calls present here and absent there: species that already "
                    "carry a non-reference residue at the site"])
    log(f"  codon pool: {n_all_calls:,} calls (all species), of which "
        f"{n_wt_calls:,} are wild type")
    log(f"  wrote {_pp}")

    # ---- the viable residue universe, shared by both classes ---------------
    # Panel residues at the site PLUS every documented resistance allele there. A
    # panel-only universe is not usable for the rapid class, because a resistance allele
    # need not occur in any wild-type panel species; using it would silently drop events.
    viable = defaultdict(Counter)
    for key, cs in codon_pool.items():
        for c in cs:
            viable[key][translate_codon(c)] += 1
    panel_only = {k: Counter(v) for k, v in viable.items()}   # for --validate
    events = read_csv(os.path.join(DIR_RESULTS, args.events))
    for e in events:
        k = (e["gene"], e["site_label"])
        if e["resistance_allele"] not in viable[k]:
            viable[k][e["resistance_allele"]] = 0
    # a residue with zero panel count still needs positive weight in the weighted null,
    # or it is unreachable by construction; give every member of the universe >= 1
    for k in viable:
        for a in viable[k]:
            viable[k][a] = max(viable[k][a], 1)

    # ---- build the two classes ---------------------------------------------
    rapid = [Obs(e["species"], e["gene"], e["site_label"], [e["codon"]],
                 e["resistance_allele"], [e["wild_type_aa"]]) for e in events]
    deep = []
    for s in read_csv(os.path.join(DIR_RESULTS, args.substitutions)):
        deep.append(Obs(s["species"], s["gene"], s["site_label"],
                        s["ancestral_codons"].split(","), s["tip_aa"],
                        s["ancestral_aas"].split(",")))
    log(f"  rapid class     : {len(rapid)} events, "
        f"{len({o.species for o in rapid})} species, {len({o.key for o in rapid})} sites")
    log(f"  deep-time class : {len(deep)} substitutions, "
        f"{len({o.species for o in deep})} species, {len({o.key for o in deep})} sites")

    shared_sp = {o.species for o in rapid} & {o.species for o in deep}
    shared_site = {o.key for o in rapid} & {o.key for o in deep}
    log(f"  shared species  : {len(shared_sp)}")
    log(f"  shared sites    : {len(shared_site)}\n")

    if args.validate:
        validate(deep, panel_only, log)
        return

    subsets = [
        ("all", lambda o: True),
        ("shared_sites", lambda o: o.key in shared_site),
        ("shared_species_and_sites",
         lambda o: o.species in shared_sp and o.key in shared_site),
    ]

    rows = []
    for sub_name, keep in subsets:
        log("-" * 78)
        log(f"SUBSET: {sub_name}")
        log("-" * 78)
        log(f"  {'null':18} {'class':10} {'n':>5} {'obs':>5} {'exp':>8} "
            f"{'ratio':>7} {'95% CI':>17} {'P':>10}")
        for null_name, fn in NULLS:
            cells = {}
            for cls_name, cls in (("rapid", rapid), ("deep_time", deep)):
                keptobs, probs = [], []
                for o in cls:
                    if not keep(o):
                        continue
                    p = fn(o, codon_pool, viable)
                    if p is None:
                        continue
                    keptobs.append(o)
                    probs.append(p)
                if not keptobs:
                    continue
                est = estimate(keptobs, probs, rng)
                cells[cls_name] = est
                rows.append(dict(subset=sub_name, null=null_name, cls=cls_name, **est))
                log(f"  {null_name:18} {cls_name:10} {est['n']:5d} {est['observed']:5d} "
                    f"{est['expected']:8.1f} {est['ratio']:7.3f} "
                    f"[{est['ci_lo']:6.3f},{est['ci_hi']:6.3f}] {est['p_exact']:10.2e}")
            if len(cells) == 2:
                r, d = cells["rapid"], cells["deep_time"]
                sep = "SEPARATED" if r["ci_lo"] > d["ci_hi"] else "OVERLAPPING"
                log(f"  {'':18} -> rapid {r['ratio']:.3f} vs deep {d['ratio']:.3f}; "
                    f"bootstrap CIs {sep}")
                rows.append(dict(subset=sub_name, null=null_name, cls="CONTRAST",
                                 n="", observed="", expected="",
                                 ratio=round(r["ratio"] / d["ratio"], 4)
                                 if d["ratio"] else "",
                                 p_exact="", ci_lo="", ci_hi="",
                                 n_species=sep))
        log("")

    out = os.path.join(DIR_RESULTS, "matched_null_contrast.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["subset", "null", "cls", "n", "observed",
                                           "expected", "ratio", "p_exact", "ci_lo",
                                           "ci_hi", "n_species"])
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {out}")
    log("\nOnly SAME-NULL comparisons between classes are interpretable. A contrast that "
        "holds\nunder every null is robust; one that holds under only the "
        "frequency-weighted null is\nnot a timescale result, because that null partly "
        "encodes accessibility.")
    log("\nStep 4.7 complete.")


if __name__ == "__main__":
    main()
