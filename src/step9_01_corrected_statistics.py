#!/usr/bin/env python3
"""
step9_01_corrected_statistics.py: the manuscript's headline statistics, computed the way
the project's six internal audits require.

Every number the paper quotes is produced HERE, from the corrected atlas and the corrected
record, so there is one authoritative source and no opportunity to re-type a value by hand
(REVIEW 06 found two numbers in RESULTS_SUMMARY.md that exist in no file at all).

What this fixes, and why each matters
-------------------------------------
1. EXACT p-values instead of permutation floors. Every "p = 1e-5" in the old output was
   literally 1/(N_perm+1) with zero exceedances, so statements like "the p-value is
   unchanged across collapse levels" carried no information - it could not go lower. The
   null statistic is a sum of INDEPENDENT Bernoullis, so the exact tail is computable by
   Poisson-binomial dynamic programming. The real values are far smaller, and they now
   move when the data move.

2. The reciprocal test, computed IN-PIPELINE and correctly. The headline
   `[[13,0],[0,6]], Fisher p = 3.7e-05` was never computed by any script in src/ - it was
   done ad hoc at a shell - and summing three 2x2 tables into one Fisher test is invalid.
   The correct test of a common odds ratio is the exact conditional test.

3. DOUBLE COUNTING in that test. The two L1014 contrasts share their entire a1 group, so
   four species were counted twice. 19 observations -> 15.

4. PHYLOGENETIC COLLAPSE of the reciprocal test, which had none, although step4_06 applies
   exactly that logic elsewhere.

5. MULTIPLE TESTING. ~44 per-site tests with no correction. RyR G4946 (p = 0.027) does not
   survive BH, so the claim that a second gene carries an individually significant site is
   withdrawn; the pooled six-gene test is the argument that survives.

6. AN ASSUMPTION-FREE effect estimate. The GEE odds ratio rests on 6 discordant
   observations and moves 15 -> 27 depending on estimator. Conditioning on each site's own
   margins and convolving the resulting hypergeometrics gives an exact test that needs no
   model at all.

7. UNCERTAINTY on the effect size, which was reported nowhere.

8. BOTH NULLS. The atlas-pool null asks "is this codon unusual among arthropods"; a
   within-record label permutation asks "did each species take the allele ITS codon
   reaches". They answer different questions and differ materially. Both are reported.

Outputs
-------
results/corrected_statistics.csv     one row per statistic, with method and provenance
logs/step9_01_corrected_statistics.log
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import DIR_PROC, DIR_RESULTS, Tee, ensure_dirs  # noqa: E402
from codon_engine import codon_reaches  # noqa: E402

SEED = 20260731
N_BOOT = 5000


def sp_name(grp_key):
    """'species_a1group' -> 'first arm'. Only for the log."""
    return "first arm" if grp_key.endswith("a1group") else "second arm"


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


# --------------------------------------------------------------------------- #
# exact distributions
# --------------------------------------------------------------------------- #
def poisson_binomial_pmf(ps):
    """Exact pmf of a sum of independent Bernoullis, by DP. O(n^2)."""
    dist = [1.0]
    for p in ps:
        nxt = [0.0] * (len(dist) + 1)
        for k, v in enumerate(dist):
            if v == 0.0:
                continue
            nxt[k] += v * (1.0 - p)
            nxt[k + 1] += v * p
        dist = nxt
    return dist


def poisson_binomial_sf(ps, obs):
    """P(X >= obs) exactly."""
    d = poisson_binomial_pmf(ps)
    return sum(d[obs:]) if obs < len(d) else 0.0


def hypergeom_pmf_vector(N, K, n):
    """pmf of Hypergeometric(N, K, n) over k = 0..min(K, n), by ratio recurrence."""
    lo, hi = max(0, n - (N - K)), min(K, n)
    if hi < lo:
        return lo, [1.0]
    # unnormalised, then normalise; recurrence avoids large factorials
    vals = [1.0]
    for k in range(lo + 1, hi + 1):
        prev = vals[-1]
        num = (K - k + 1) * (n - k + 1)
        den = k * (N - K - n + k)
        vals.append(prev * num / den if den else 0.0)
    tot = sum(vals)
    return lo, [v / tot for v in vals] if tot else [1.0]


def convolve_offset(dists):
    """Convolve [(offset, pmf), ...] -> (offset, pmf)."""
    off_tot, cur = 0, [1.0]
    for off, pmf in dists:
        off_tot += off
        nxt = [0.0] * (len(cur) + len(pmf) - 1)
        for i, a in enumerate(cur):
            if a == 0.0:
                continue
            for j, b in enumerate(pmf):
                nxt[i + j] += a * b
        cur = nxt
    return off_tot, cur


def exact_common_or_test(tables):
    """Exact conditional test of a common odds ratio = 1 across 2x2 tables.

    Under H0 each table's a-cell is central hypergeometric given its margins. The exact
    null of sum(a) is the convolution. Returns (observed_sum, expected, p_upper).
    Correct replacement for summing tables into one Fisher test, which is invalid.
    """
    dists, obs, exp = [], 0, 0.0
    for a, b, c, d in tables:
        N, K, n = a + b + c + d, a + c, a + b
        lo, pmf = hypergeom_pmf_vector(N, K, n)
        dists.append((lo, pmf))
        obs += a
        exp += sum((lo + i) * p for i, p in enumerate(pmf))
    off, pmf = convolve_offset(dists)
    p = sum(p for i, p in enumerate(pmf) if off + i >= obs)
    return obs, exp, min(1.0, p)


def benjamini_hochberg(pvals):
    """Return BH q-values in the input order."""
    idx = sorted(range(len(pvals)), key=lambda i: pvals[i])
    m = len(pvals)
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(idx), start=1):
        k = m - rank + 1
        val = min(prev, pvals[i] * m / k)
        q[i] = val
        prev = val
    return q


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default="codon_atlas.csv")
    ap.add_argument("--record", default="resistance_record_corrected.csv")
    ap.add_argument("--reconciliation", default="species_reconciliation_final.csv")
    ap.add_argument("--events", default="core_test_events_primary.csv")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    rng = random.Random(SEED)

    ensure_dirs()
    log = Tee("step9_01_corrected_statistics" + args.tag)
    log("=== Step 9.1  corrected headline statistics ===\n")
    log(f"atlas  : {args.atlas}")
    log(f"record : {args.record}")
    log(f"events : {args.events}\n")

    out = []

    def emit(name, value, p, method, detail="", q=""):
        out.append({"statistic": name, "value": value, "p_value": p,
                    "q_value": q, "method": method, "detail": detail})

    # ---- codon pools, deduplicated to ONE codon per species -----------------
    codon, wt_of = {}, {}
    for r in read_csv(os.path.join(DIR_PROC, args.atlas)):
        if r["is_wildtype"] != "True":
            continue
        codon[(r["species"], r["gene"], r["site_label"])] = r["codon"]
        wt_of[(r["gene"], r["site_label"])] = r["wild_type_aa"]
    pool = defaultdict(list)
    for (sp, g, s), c in codon.items():
        pool[(g, s)].append(c)
    log(f"codon pools: {len(pool)} sites, {len(codon)} species-site calls")

    events = read_csv(os.path.join(DIR_RESULTS, args.events))
    log(f"events: {len(events)}")

    # ---- 1. EXACT pooled test ----------------------------------------------
    log("\n" + "=" * 70)
    log("1  POOLED ACCESSIBILITY: exact, replacing the permutation floor")
    log("=" * 70)
    ps, obs = [], 0
    for e in events:
        key = (e["gene"], e["site_label"])
        cs = pool.get(key)
        if not cs:
            continue
        a = e["resistance_allele"]
        p_i = sum(1 for c in cs if codon_reaches(c, a)) / len(cs)
        ps.append(p_i)
        obs += 1 if e["accessible"] == "True" else 0
    exp = sum(ps)
    p_exact = poisson_binomial_sf(ps, obs)
    log(f"  events            : {len(ps)}")
    log(f"  accessible        : {obs} ({obs/len(ps):.1%})")
    log(f"  expected (null)   : {exp:.6f} ({exp/len(ps):.1%})")
    log(f"  enrichment        : {obs/exp:.6f}x    [= {obs} / {exp:.6f}]")
    log(f"  EXACT p (Poisson-binomial) = {p_exact:.4g}")
    log("    (the old permutation reported 1e-5, which was simply 1/(N+1) with 0 exceedances)")
    # The expectation is emitted as a quantity in its own right, at full
    # precision. Recorded only as a rounded string, a real disagreement between
    # this script and step2_01 is indistinguishable from display rounding and
    # the ratio cannot be reproduced from the numbers the file reports.
    emit("expected_null_pooled", round(exp, 6), "",
         "Poisson-binomial exact expectation, sum(p_i) over all events",
         f"{len(ps)} events; sum of per-event accessibility probabilities")
    emit("pooled_enrichment", round(obs / exp, 6), p_exact,
         "Poisson-binomial exact tail",
         f"{obs}/{len(ps)} accessible, expected {exp:.6f} (exact; not rounded)")

    # ---- 2. informative / uninformative ------------------------------------
    log("\n" + "=" * 70)
    log("2  STRATA: and why the uninformative stratum proves nothing")
    log("=" * 70)
    inf_p, inf_o, uninf_n, uninf_o = [], 0, 0, 0
    for e, p_i in zip([e for e in events if (e["gene"], e["site_label"]) in pool], ps):
        acc = e["accessible"] == "True"
        if 0.0 < p_i < 1.0:
            inf_p.append(p_i)
            inf_o += acc
        else:
            uninf_n += 1
            uninf_o += acc
    if inf_p:
        ie = sum(inf_p)
        ip = poisson_binomial_sf(inf_p, inf_o)
        log(f"  informative   n={len(inf_p):4d}  obs {inf_o} ({inf_o/len(inf_p):.1%})  "
            f"exp {ie:.2f}  enrich {inf_o/ie:.3f}x  EXACT p = {ip:.4g}")
        emit("informative_enrichment", round(inf_o / ie, 6), ip,
             "Poisson-binomial exact tail",
             f"{inf_o}/{len(inf_p)}, expected {ie:.6f}")
    # Concentration check. VGSC L1014 supplies about a third of the informative events,
    # so the obvious objection is that one site carries the confirmatory test. Recompute
    # it with that site, and then the whole of VGSC, removed.
    inf_rows = [(e, p_i) for e, p_i in
                zip([e for e in events if (e["gene"], e["site_label"]) in pool], ps)
                if 0.0 < p_i < 1.0]
    for label, keep in (
            ("informative_excl_VGSC_L1014",
             lambda e: not (e["gene"] == "VGSC" and e["site_label"] == "L1014")),
            ("informative_excl_VGSC", lambda e: e["gene"] != "VGSC")):
        sub = [(e, p_i) for e, p_i in inf_rows if keep(e)]
        if not sub:
            continue
        o = sum(1 for e, _ in sub if e["accessible"] == "True")
        ex = sum(p_i for _, p_i in sub)
        pv = poisson_binomial_sf([p_i for _, p_i in sub], o)
        log(f"  {label:30} n={len(sub):4d}  obs {o}  exp {ex:.2f}  "
            f"enrich {o/ex:.3f}x  EXACT p = {pv:.4g}")
        emit(label, round(o / ex, 4), pv, "Poisson-binomial exact tail",
             f"{o}/{len(sub)}")

    # ---- 2b. the per-site conditional formulation --------------------------
    # The test above conditions on the allele that evolved: the null draws a
    # codon and asks whether it reaches THAT allele. An allele no panel codon
    # reaches in one step therefore has p_i = 0, contributes no likelihood,
    # and cannot lower the exact tail however often it occurs. Conditioning
    # instead on the repertoire of alleles reported at the site makes a
    # species informative whenever its codon reaches some but not all of that
    # repertoire, whichever allele it acquired, so an obligate two-step allele
    # enters as a failure. step3_03 builds that table and gives it a
    # Monte-Carlo P; this gives it the same exact tail as every other test
    # reported here, which is what the main text needs in order to quote it.
    pool_events = [e for e in events if (e["gene"], e["site_label"]) in pool]
    # step3_03 names its two outputs in step, so the conditional table is
    # found from whichever events file this run was pointed at rather than
    # from --tag, which does not apply to step3_03's filenames.
    if "core_test_events" not in args.events:
        sys.exit("FATAL: --events is %r, from which the matching "
                 "core_test_conditional table cannot be named."
                 % args.events)
    cond_path = os.path.join(
        DIR_RESULTS, args.events.replace("core_test_events",
                                         "core_test_conditional"))
    if not os.path.exists(cond_path):
        sys.exit("FATAL: %s is missing, so the per-site conditional estimand "
                 "cannot be computed. Run src/step3_03_core_test.py first."
                 % cond_path)
    cond = read_csv(cond_path)
    cps, cobs = [], 0
    for c in cond:
        n_rep, n_acc = int(c["n_repertoire"]), int(c["n_accessible"])
        p_c = n_acc / n_rep
        # The stored probability IS the null of this test, so it is recomputed
        # from the two counts rather than trusted.
        if abs(p_c - float(c["prob_by_chance"])) > 5e-5:
            sys.exit("FATAL: %s, %s %s %s: prob_by_chance is %s but %d/%d = "
                     "%.6f. The conditional null is not what the table says."
                     % (os.path.basename(cond_path), c["species"], c["gene"],
                        c["site_label"], c["prob_by_chance"], n_acc, n_rep,
                        p_c))
        cps.append(p_c)
        cobs += c["observed_allele_accessible"] == "True"
    cexp = sum(cps)
    cp = poisson_binomial_sf(cps, cobs)
    # Which obligate two-step events this formulation can penalise. They are
    # the events the primary test scored inaccessible, and they appear here
    # only where the species' codon reaches part of its site's repertoire;
    # where it reaches none of it the event is uninformative here too.
    two_step = {(e["species"], e["gene"], e["site_label"],
                 e["resistance_allele"]) for e in pool_events
                if e["accessible"] != "True"}
    weighted = {(c["species"], c["gene"], c["site_label"],
                 c["resistance_allele"]) for c in cond
                if c["observed_allele_accessible"] != "True"}
    if not weighted <= two_step:
        sys.exit("FATAL: the per-site conditional table scores %d event(s) "
                 "inaccessible that the primary test scores accessible: %s. "
                 "The two tables disagree about one-step accessibility."
                 % (len(weighted - two_step), sorted(weighted - two_step)))
    log(f"  per-site cond n={len(cond):4d}  obs {cobs}  exp {cexp:.2f}  "
        f"enrich {cobs/cexp:.3f}x  EXACT p = {cp:.4g}")
    log(f"    obligate two-step events carrying failure weight: "
        f"{len(weighted)} of {len(two_step)}")
    pc_path = os.path.join(DIR_RESULTS,
                           f"step9_01_per_site_conditional{args.tag}.csv")
    with open(pc_path, "w", newline="", encoding="utf-8") as fh:
        _w = csv.writer(fh)
        _w.writerow(["n", "observed", "expected", "ratio", "mean_p",
                     "exact_p", "two_step_weighted", "two_step_total"])
        _w.writerow([len(cond), cobs, round(cexp, 4), round(cobs / cexp, 4),
                     round(cexp / len(cond), 4), "%.6g" % cp,
                     len(weighted), len(two_step)])
    log(f"  wrote {os.path.basename(pc_path)}")
    emit("per_site_conditional_enrichment", round(cobs / cexp, 6), cp,
         "Poisson-binomial exact tail; the allele is drawn from the "
         "repertoire reported at the site, so an obligate two-step allele "
         "counts as a failure",
         f"{cobs}/{len(cond)} accessible, expected {cexp:.6f}; "
         f"{len(weighted)} of {len(two_step)} obligate two-step events carry "
         f"failure weight")

    log(f"  uninformative n={uninf_n:4d}  obs {uninf_o}")
    log("    Every event here has p_i in {0,1}, so observed EQUALS expected by arithmetic")
    log("    and p = 1.0 by construction. It cannot corroborate anything, and the pooled")
    log("    and informative tests are ONE test, not two - the uninformative events")
    log("    contribute deterministically to both numerator and null.")

    # ---- 3. BOTH NULLS ------------------------------------------------------
    log("\n" + "=" * 70)
    log("3  TWO NULLS, TWO QUESTIONS")
    log("=" * 70)
    by_site_events = defaultdict(list)
    for e in events:
        by_site_events[(e["gene"], e["site_label"])].append(e)
    tot_obs = sum(1 for e in events if e["accessible"] == "True")
    exp_rec, n_rec = 0.0, 0
    # Keep the per-event probabilities, not just their sum: the exact tail for
    # this null is a Poisson-binomial over exactly these p_i, the same as for
    # the atlas-pool null. Without them the row went out with an empty P
    # column, and a downstream figure script supplied a typed literal instead.
    rec_ps, rec_by_species = [], defaultdict(list)
    for key, evs in by_site_events.items():
        codons_here = [e["codon"] for e in evs]
        if not codons_here:
            continue
        for e in evs:
            a = e["resistance_allele"]
            p_i = sum(1 for c in codons_here if codon_reaches(c, a)) / len(codons_here)
            exp_rec += p_i
            rec_ps.append(p_i)
            rec_by_species[e["species"]].append((e["accessible"] == "True", p_i))
            n_rec += 1
    p_within_record = poisson_binomial_sf(rec_ps, tot_obs)
    # Same species-cluster bootstrap as the atlas-pool null, so that the three
    # null constructions can be drawn with comparable intervals rather than as
    # bare points (Supplementary Fig. 4c).
    # Its own generator, seeded independently. Sharing one stream means that
    # adding a bootstrap anywhere shifts every bootstrap after it: the pooled
    # interval moved in its third decimal the first time this block was added.
    _rng_wr = random.Random(SEED + 1)
    _rsp = list(rec_by_species)
    _rr = []
    for _ in range(N_BOOT):
        _o = _t = 0.0
        for _ in _rsp:
            for _acc, _pi in rec_by_species[_rng_wr.choice(_rsp)]:
                _o += _acc
                _t += _pi
        if _t > 0:
            _rr.append(_o / _t)
    _rr.sort()
    wr_lo, wr_hi = _rr[int(.025 * len(_rr))], _rr[int(.975 * len(_rr))]
    log(f"  atlas-pool null      : expected {exp:.2f}/{len(ps)}  -> {obs/exp:.3f}x")
    log(f"  within-record null   : expected {exp_rec:.2f}/{n_rec}  -> {tot_obs/exp_rec:.3f}x"
        f"   95% CI [{wr_lo:.3f}, {wr_hi:.3f}]")
    log("    The first asks 'is this codon unusual among arthropods generally'; the second")
    log("    asks 'did each species take the allele ITS OWN codon reaches', permuting only")
    log("    among species that actually appear in the record. Record species already")
    log("    carry more accessibility-friendly codons than arthropods at large, so part of")
    log("    the atlas-pool effect is a property of which species get studied.")
    emit("within_record_null_enrichment", round(tot_obs / exp_rec, 6),
         p_within_record,
         "Poisson-binomial exact tail, permuting within record species",
         f"{tot_obs}/{n_rec} accessible, expected {exp_rec:.6f}")
    emit("within_record_null_enrichment_CI", f"[{wr_lo:.3f}, {wr_hi:.3f}]", "",
         "species-cluster bootstrap", f"{N_BOOT} resamples over {len(_rsp)} species")

    # ---- 4. bootstrap CI on the effect size ---------------------------------
    log("\n" + "=" * 70)
    log("4  UNCERTAINTY ON THE EFFECT SIZE (species-cluster bootstrap)")
    log("=" * 70)
    by_species = defaultdict(list)
    for e, p_i in zip([e for e in events if (e["gene"], e["site_label"]) in pool], ps):
        by_species[e["species"]].append((e["accessible"] == "True", p_i))
    sp_list = list(by_species)
    ratios = []
    for _ in range(N_BOOT):
        o = t = 0.0
        for _ in sp_list:
            for acc, p_i in by_species[rng.choice(sp_list)]:
                o += acc
                t += p_i
        if t > 0:
            ratios.append(o / t)
    ratios.sort()
    lo_, hi_ = ratios[int(.025 * len(ratios))], ratios[int(.975 * len(ratios))]
    log(f"  pooled enrichment {obs/exp:.3f}x   95% CI [{lo_:.3f}, {hi_:.3f}]  "
        f"({N_BOOT} resamples over {len(sp_list)} species)")
    emit("pooled_enrichment_CI", f"[{lo_:.3f}, {hi_:.3f}]", "",
         "species-cluster bootstrap", f"{N_BOOT} resamples")

    # ---- 5. per-site tests with BH ------------------------------------------
    log("\n" + "=" * 70)
    log("5  PER-SITE TESTS WITH MULTIPLE-TESTING CORRECTION")
    log("=" * 70)
    site_rows = []
    for key, evs in sorted(by_site_events.items()):
        cs = pool.get(key)
        if not cs:
            continue
        pl = [sum(1 for c in cs if codon_reaches(c, e["resistance_allele"])) / len(cs)
              for e in evs]
        o = sum(1 for e in evs if e["accessible"] == "True")
        if not pl:
            continue
        pv = poisson_binomial_sf(pl, o)
        degenerate = all(x in (0.0, 1.0) for x in pl)
        site_rows.append({"gene": key[0], "site": key[1], "n": len(evs), "obs": o,
                          "exp": sum(pl), "p": pv, "degenerate": degenerate})
    nondeg = [r for r in site_rows if not r["degenerate"]]
    qs = benjamini_hochberg([r["p"] for r in nondeg])
    for r, q in zip(nondeg, qs):
        r["q"] = q
    log(f"  sites tested: {len(site_rows)};  non-degenerate (p<1 attainable): {len(nondeg)}")
    log(f"  {'gene':9} {'site':7} {'n':>3} {'obs':>4} {'exp':>7} {'p':>10} {'BH q':>10}")
    for r in sorted(nondeg, key=lambda x: x["p"]):
        star = "  *" if r["q"] < 0.05 else ""
        log(f"  {r['gene']:9} {r['site']:7} {r['n']:>3} {r['obs']:>4} {r['exp']:>7.2f} "
            f"{r['p']:>10.3g} {r['q']:>10.3g}{star}")
        emit(f"site_{r['gene']}_{r['site']}", f"{r['obs']}/{r['n']}", r["p"],
             "Poisson-binomial exact, P unadjusted; q is Benjamini-Hochberg "
             f"over the {len(nondeg)} non-degenerate sites",
             f"expected {r['exp']:.3f} of {r['n']}", q=round(r["q"], 6))
    surv = [r for r in nondeg if r["q"] < 0.05]
    log(f"\n  sites surviving BH q<0.05: {len(surv)} "
        f"({', '.join(r['gene']+' '+r['site'] for r in surv) or 'none'})")
    log("  NOTE: the claim that a second gene carries an individually significant site is")
    log("  WITHDRAWN if RyR G4946 is not in that list. The pooled multi-gene test is the")
    log("  argument that survives correction.")

    # ---- 6. reciprocal, done properly ---------------------------------------
    log("\n" + "=" * 70)
    log("6  RECIPROCAL DESIGN: exact conditional test, de-duplicated, collapsed")
    log("=" * 70)
    # run_all passes no --tag here, so this resolved to the SUPERSEDED _t1 file.
    # Prefer the canonical _primary and fall back only to older tags.
    recf = ""
    for _t in ([args.tag] if args.tag else []) + ["_primary", "_t1", "_final", ""]:
        _c = os.path.join(DIR_RESULTS, f"core_test_reciprocal{_t}.csv")
        if os.path.exists(_c):
            recf = _c
            break
    log(f"  reciprocal input: {os.path.basename(recf) if recf else 'NONE FOUND'}"),
    if os.path.exists(recf):
        rec = read_csv(recf)

        def _mem(s):
            return [z.split("(")[0] for z in (s or "").split(";") if z]

        # Which alleles each species actually carries at each site. The earlier
        # de-duplication took the FIRST n_obs members of each arm instead,
        # which gave the right count by luck and the wrong species.
        alle_of = defaultdict(set)
        for e in events:
            alle_of[(e["gene"], e["site_label"],
                     e["species"])].add(e["resistance_allele"])

        tables, seen_obs = [], {}
        for x in rec:
            a = int(x["a1_observed_in_a1group"]); b = int(x["a2_observed_in_a1group"])
            c = int(x["a1_observed_in_a2group"]); d = int(x["a2_observed_in_a2group"])
            tables.append((a, b, c, d))
            log(f"  {x['gene']} {x['site_label']} {x['allele_1']} vs {x['allele_2']}: "
                f"[[{a},{b}],[{c},{d}]]")
            # One observation per (site, allele, species) that actually carries
            # the allele. The naive pooled total double-counts because the two
            # L1014 contrasts share their first arm.
            for grp_key, al, n_obs in (
                    ("species_a1group", x["allele_1"], a),
                    ("species_a2group", x["allele_2"], d)):
                hits = [sp_ for sp_ in _mem(x.get(grp_key))
                        if al in alle_of[(x["gene"], x["site_label"], sp_)]]
                if len(hits) != n_obs:
                    sys.exit("FATAL: %s %s %s arm reports %d observation(s) of "
                             "%s but %d member species carry it (%s). The "
                             "reciprocal table and the event record disagree."
                             % (x["gene"], x["site_label"], grp_key, n_obs, al,
                                len(hits), ", ".join(hits)))
                log(f"      {al} observed in {sp_name(grp_key)}: "
                    f"{', '.join(hits) if hits else 'none'}")
                for sp_ in hits:
                    seen_obs[(x["gene"], x["site_label"], al, sp_)] = 1
        naive = sum(t[0] for t in tables) + sum(t[3] for t in tables)
        o, e_, p = exact_common_or_test(tables)
        seen_species = {k[3] for k in seen_obs}
        log(f"\n  naive pooled 2x2 concordant count: {naive}")
        log(f"  DISTINCT concordant (site, allele, species) observations: {len(seen_obs)}")
        log(f"  DISTINCT species among them                             : "
            f"{len(seen_species)}")
        log(f"    {naive - len(seen_obs)} observation(s) are duplicates: the two L1014")
        log("    contrasts share their ENTIRE first arm, so each species observed")
        log("    carrying that allele enters both tables.")
        log(f"    {len(seen_obs) - len(seen_species)} species contributes two different")
        log("    alleles, which is why the species count is lower again.")
        log(f"\n  exact conditional test of common OR = 1:")
        log(f"    observed sum(a) = {o}, expected {e_:.3f}, one-sided p = {p:.4g}")
        log("    (summing the tables into a single Fisher test is invalid and inflates this)")

        # The de-duplication above corrects the COUNT that is reported; it does
        # not make the strata independent. Entering both L1014 tables still
        # reuses the shared first arm inside the conditional likelihood.
        # Collapsing each site to one stratum - its shared first allele against
        # the union of the others - makes the strata independent, at the cost of
        # the within-site resolution between S and W. Both are reported and the
        # independent one is the test that assumes nothing it does not have.
        by_site = defaultdict(list)
        for x in rec:
            by_site[(x["gene"], x["site_label"])].append(x)
        ind_tables, ind_labels = [], []
        for (g_, s_), xs in sorted(by_site.items()):
            first = {x["allele_1"] for x in xs}
            if len(first) != 1:
                sys.exit("FATAL: the contrasts at %s %s do not share a single "
                         "first allele (%s), so they cannot be collapsed into "
                         "one stratum." % (g_, s_, ", ".join(sorted(first))))
            al1 = first.pop()
            others = {x["allele_2"] for x in xs}
            arm_a, arm_b = set(), set()
            for x in xs:
                arm_a |= set(_mem(x["species_a1group"]))
                arm_b |= set(_mem(x["species_a2group"]))
            if arm_a & arm_b:
                sys.exit("FATAL: %s %s has %d species in both arms (%s); the "
                         "collapsed stratum would double-count them."
                         % (g_, s_, len(arm_a & arm_b),
                            ", ".join(sorted(arm_a & arm_b))))
            t = (sum(1 for sp_ in arm_a if al1 in alle_of[(g_, s_, sp_)]),
                 sum(1 for sp_ in arm_a if alle_of[(g_, s_, sp_)] & others),
                 sum(1 for sp_ in arm_b if al1 in alle_of[(g_, s_, sp_)]),
                 sum(1 for sp_ in arm_b if alle_of[(g_, s_, sp_)] & others))
            ind_tables.append(t)
            ind_labels.append("%s %s %s vs %s (%d and %d species)"
                              % (g_, s_, al1, "/".join(sorted(others)),
                                 len(arm_a), len(arm_b)))
            log(f"  collapsed {ind_labels[-1]}: [[{t[0]},{t[1]}],[{t[2]},{t[3]}]]")
        o2, e2, p2 = exact_common_or_test(ind_tables)
        log(f"\n  exact conditional test over {len(ind_tables)} INDEPENDENT strata:")
        log(f"    observed sum(a) = {o2}, expected {e2:.3f}, one-sided p = {p2:.4g}")
        # The cell counts themselves, so the appendix can print them and the
        # three totals (naive, distinct observations, distinct species) can be
        # reconciled on the page instead of asserted.
        rp = os.path.join(DIR_RESULTS,
                          f"step9_01_reciprocal_accounting{args.tag}.csv")
        with open(rp, "w", newline="", encoding="utf-8") as fh:
            rw = csv.writer(fh)
            rw.writerow(["kind", "gene", "site", "allele_1", "allele_2",
                         "a", "b", "c", "d", "n_arm_1", "n_arm_2", "note"])
            for x, t in zip(rec, tables):
                rw.writerow(["contrast", x["gene"], x["site_label"],
                             x["allele_1"], x["allele_2"], t[0], t[1], t[2],
                             t[3], len(_mem(x["species_a1group"])),
                             len(_mem(x["species_a2group"])),
                             "as published; the two L1014 rows share arm 1"])
            for lab, t in zip(ind_labels, ind_tables):
                rw.writerow(["independent_stratum", "", lab, "", "", t[0],
                             t[1], t[2], t[3], "", "",
                             "one stratum per site, arms disjoint"])
            rw.writerow(["total", "", "", "", "", naive, len(seen_obs),
                         len(seen_species), "", "", "",
                         "a = naive pooled concordant, b = distinct "
                         "(site, allele, species), c = distinct species"])
            rw.writerow(["overlapping_total", "", "", "", "", o,
                         round(e_, 3), "%.6g" % p, "", "", "",
                         "a = observed, b = expected, c = one-sided P, "
                         "three contrasts sharing the L1014 first arm"])
            rw.writerow(["independent_total", "", "", "", "", o2,
                         round(e2, 3), "%.6g" % p2, "", "", "",
                         "a = observed, b = expected, c = one-sided P, "
                         "one non-overlapping stratum per site"])
        log(f"  wrote {os.path.basename(rp)}")
        emit("reciprocal_independent_strata", f"{o2}/{o2}", p2,
             "exact conditional test of a common odds ratio over "
             f"{len(ind_tables)} non-overlapping strata, one per site",
             f"{o2} concordant of {o2}; expected {e2:.3f} under a common odds "
             f"ratio of 1; each site enters once, its shared first allele "
             f"against the union of the others, so no arm is reused; "
             f"strata: {'; '.join(ind_labels)}")
        _sites = sorted({f"{x['gene']} {x['site_label']}" for x in rec})
        _pairs = ", ".join(f"{x['gene']} {x['site_label']} {x['allele_1']} vs "
                           f"{x['allele_2']}" for x in rec)
        emit("reciprocal_concordant_observations", f"{o}/{o}", p,
             "exact conditional test of a common odds ratio across "
             f"{len(tables)} contrasts at {len(_sites)} sites",
             f"{o} concordant of {o} informative; expected {e_:.3f} under a "
             f"common odds ratio of 1; {len(seen_obs)} distinct "
             f"(site, allele, species) observations, not {naive}; contrasts: "
             f"{_pairs}")
    else:
        log("  reciprocal file not found - skipped")

    p1 = os.path.join(DIR_RESULTS, f"corrected_statistics{args.tag}.csv")
    with open(p1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["statistic", "value", "p_value",
                                           "q_value", "method", "detail"])
        w.writeheader(); w.writerows(out)
    log(f"\nwrote {p1}")
    log("\nStep 9.1 complete.")


if __name__ == "__main__":
    main()
