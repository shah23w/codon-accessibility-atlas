#!/usr/bin/env python3
"""
step2_01_breach_rate_and_estimands.py
=====================================
The breach rate, and the estimands the confirmatory test rests on.

WHAT THIS FIXES
---------------
The confirmatory test reports 48/48 informative events accessible, exact
P = 4.6e-12. That test is arithmetically correct. Three things about it still
need addressing:

 1. ESTIMAND. The "informative" subset (0 < p < 1) is defined by the codon
    distribution at the site, which is a function of the explanatory variable.
    Conditioning on it narrows the claim: it supports "at sites where species
    differ in reachability, the evolved allele was always reachable", not the
    broader claim in the abstract. Both estimands must be reported, labelled.

 2. NO MEASURE OF ABSOLUTE STRENGTH. The null is a one-step accessibility
    indicator with no parameter for two-step evolution, so it cannot say how
    often the gate is crossed. The five inaccessible events are the ONLY direct
    observations of a crossing, and the current design discards them as
    "uninformative". This script reports them as the primary new quantity:
    the BREACH RATE, with an exact confidence interval.

 3. SATURATION. With every informative event accessible, the ratio is pinned at
    exactly 1 / mean(p) and carries no information about magnitude. Stated
    explicitly here so it cannot be mistaken for an effect size.

It also audits the five breach events against the replicate-assembly strata from
Step 1, because a "confirmed in a second assembly" claim is only meaningful if
that second assembly is not a GCA/GCF twin of the first.

INPUTS
  results/core_test_events_primary.csv        108 events
  data/processed/codon_atlas.csv              per species x site x allele paths
  results/assembly_replication_pair_classes.csv   from Step 1 (optional)
  results/assembly_replication.csv                from Step 1 (optional)

OUTPUTS
  results/step2_estimands.csv          the three strata and both estimands
  results/step2_breach_events.csv      the five breaches, fully annotated
  results/step2_breach_rate.csv        breach rate overall, per gene, per site
  results/step2_summary.md             text for the manuscript

Standard library only. No third-party packages, so it can be run\nwith any Python 3.
"""

import csv
import os
import sys
from collections import Counter, defaultdict

import math

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")

EVENTS = os.path.join(RES, "core_test_events_primary.csv")
ATLAS = os.path.join(ROOT, "data", "processed", "codon_atlas.csv")
PAIRCLASS = os.path.join(RES, "assembly_replication_pair_classes.csv")
REPLCALLS = os.path.join(RES, "assembly_replication.csv")

OUT_EST = os.path.join(RES, "step2_estimands.csv")
OUT_BRE = os.path.join(RES, "step2_breach_events.csv")
OUT_RATE = os.path.join(RES, "step2_breach_rate.csv")
OUT_MD = os.path.join(RES, "step2_summary.md")

for p in (EVENTS, ATLAS):
    if not os.path.exists(p):
        sys.exit("ERROR: missing required input:\n  %s" % p)


def truthy(v):
    return str(v).strip().lower() in ("true", "1", "yes", "t")


def num(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def binom_cdf(k, n, p):
    """P(X <= k) for X ~ Binomial(n, p). Exact, standard library only."""
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    return sum(math.comb(n, i) * p ** i * (1.0 - p) ** (n - i)
               for i in range(0, k + 1))


def _bisect(f, lo, hi, tol=1e-12, iters=200):
    """Find x in [lo, hi] with f(x) = 0, assuming f changes sign. Monotone f."""
    flo = f(lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if abs(fmid) < tol or (hi - lo) < tol:
            return mid
        if (fmid < 0) == (flo < 0):
            lo, flo = mid, fmid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def clopper_pearson(k, n, conf=0.95):
    """
    Exact (Clopper-Pearson) binomial confidence interval, computed by inverting
    the exact binomial CDF with bisection. No third-party dependency, so this
    script runs on any Python 3 installation.

    Lower bound solves  P(X >= k | p) = alpha/2
    Upper bound solves  P(X <= k | p) = alpha/2
    """
    if n == 0:
        return (float("nan"), float("nan"))
    a = 1.0 - conf
    if k == 0:
        lo = 0.0
    else:
        lo = _bisect(lambda p: (1.0 - binom_cdf(k - 1, n, p)) - a / 2.0, 0.0, 1.0)
    if k == n:
        hi = 1.0
    else:
        hi = _bisect(lambda p: binom_cdf(k, n, p) - a / 2.0, 0.0, 1.0)
    return (lo, hi)


def poisson_binomial_pmf(probs):
    """Exact PMF of a sum of independent non-identical Bernoullis, by convolution."""
    pmf = [1.0]
    for p in probs:
        nxt = [0.0] * (len(pmf) + 1)
        for i, v in enumerate(pmf):
            nxt[i] += v * (1.0 - p)
            nxt[i + 1] += v * p
        pmf = nxt
    return pmf


def pb_upper_tail(probs, k):
    """P(X >= k) exactly."""
    pmf = poisson_binomial_pmf(probs)
    return sum(pmf[k:]) if k <= len(pmf) - 1 else 0.0


# ======================================================================
# load
# ======================================================================
with open(EVENTS, newline="") as fh:
    events = list(csv.DictReader(fh))

print("=" * 74)
print("STEP 2  Estimands and the breach rate")
print("=" * 74)
print()
print("Events loaded: %d  (from %s)" % (len(events), os.path.relpath(EVENTS, ROOT)))
print()

# ======================================================================
# PART A -- the partition, stated exactly
# ======================================================================
for e in events:
    e["_p"] = num(e["site_supply_fraction"], float("nan"))
    e["_acc"] = truthy(e["accessible"])

S1 = [e for e in events if e["_p"] == 1.0]
S0 = [e for e in events if e["_p"] == 0.0]
SI = [e for e in events if 0.0 < e["_p"] < 1.0]

assert len(S1) + len(S0) + len(SI) == len(events), "partition does not cover all events"

obs_all = sum(1 for e in events if e["_acc"])
exp_all = sum(e["_p"] for e in events)
obs_inf = sum(1 for e in SI if e["_acc"])
exp_inf = sum(e["_p"] for e in SI)
mean_p_inf = exp_inf / len(SI)

print("-" * 74)
print("PART A   The three strata")
print("-" * 74)
print("%-34s %5s %10s %10s" % ("stratum", "n", "observed", "expected"))
print("%-34s %5d %10d %10.2f" % ("p = 1  (all codons reach)", len(S1),
                                 sum(1 for e in S1 if e["_acc"]),
                                 sum(e["_p"] for e in S1)))
print("%-34s %5d %10d %10.2f" % ("p = 0  (no codon reaches)", len(S0),
                                 sum(1 for e in S0 if e["_acc"]),
                                 sum(e["_p"] for e in S0)))
print("%-34s %5d %10d %10.2f" % ("0 < p < 1  (informative)", len(SI),
                                 obs_inf, exp_inf))
print("%-34s %5d %10d %10.2f" % ("ALL EVENTS", len(events), obs_all, exp_all))
print()
print("The null is correct by construction on the first two strata:")
print("  p = 1 events must be accessible, and %d of %d were." % (sum(1 for e in S1 if e['_acc']), len(S1)))
print("  p = 0 events cannot be accessible, and %d of %d were." % (sum(1 for e in S0 if e['_acc']), len(S0)))
print("  All evidence therefore lies in the %d informative events." % len(SI))
print()

p_inf = pb_upper_tail([e["_p"] for e in SI], obs_inf)
p_all = pb_upper_tail([e["_p"] for e in events], obs_all)

print("-" * 74)
print("TWO ESTIMANDS, LABELLED")
print("-" * 74)
print("CONDITIONAL  (the pre-specified confirmatory test)")
print("  question : at sites where reachability varies between species, was the")
print("             evolved allele reachable from this species' own codon?")
print("  result   : %d / %d accessible" % (obs_inf, len(SI)))
print("  ratio    : %.4f   expected %.2f" % (obs_inf / exp_inf, exp_inf))
print("  exact P  : %.4g" % p_inf)
print("  SATURATED: every informative event is accessible, so the ratio equals")
print("             1 / mean(p) = 1 / %.4f = %.4f exactly and cannot be exceeded."
      % (mean_p_inf, 1.0 / mean_p_inf))
print("             The evidence is in the tail, NOT in the magnitude.")
print()
print("UNCONDITIONAL  (what the abstract's wording actually claims)")
print("  question : across all documented events, was the evolved allele reachable?")
print("  result   : %d / %d accessible (%.1f%%)" % (obs_all, len(events),
                                                    100.0 * obs_all / len(events)))
print("  ratio    : %.4f   expected %.2f" % (obs_all / exp_all, exp_all))
print("  exact P  : %.4g   (identical to the conditional test: the degenerate" % p_all)
print("             strata contribute no likelihood)")
print()

# ======================================================================
# PART B -- the breach rate
# ======================================================================
breaches = [e for e in events if not e["_acc"]]
k, n = len(breaches), len(events)
lo, hi = clopper_pearson(k, n)

print("-" * 74)
print("PART B   The breach rate  <-- the quantity the current design is missing")
print("-" * 74)
print("Events in which the evolved allele was NOT one step from the species' codon:")
print("  %d of %d  =  %.2f%%   95%% CI (Clopper-Pearson) %.2f%% to %.2f%%"
      % (k, n, 100.0 * k / n, 100.0 * lo, 100.0 * hi))
print()
print("This is an absolute measure of how often the gate is crossed. It is a")
print("RESULT, not a limitation, and it is the number a theoretical prediction")
print("must be compared against (Step 9).")
print()

# ======================================================================
# PART B2 -- audit the breaches against Step 1's replicate strata
# ======================================================================
pair_class = {}
if os.path.exists(PAIRCLASS):
    with open(PAIRCLASS, newline="") as fh:
        for r in csv.DictReader(fh):
            pair_class[r["species"].strip()] = r["pair_class"]

repl_codon = defaultdict(list)
if os.path.exists(REPLCALLS):
    with open(REPLCALLS, newline="") as fh:
        for r in csv.DictReader(fh):
            key = (r["species"].strip(), r["gene"].strip(), r["site_label"].strip())
            repl_codon[key].append(r)

print("-" * 74)
print("PART B2  Are the breaches real, or calling artefacts?")
print("-" * 74)
print("A breach is only credible if the codon call is right. A 'confirmed in a")
print("second assembly' claim is meaningless when that second assembly is a")
print("GCA/GCF twin sharing identical DNA.")
print()

breach_rows = []
for e in breaches:
    sp = e["species"].strip()
    key = (sp, e["gene"].strip(), e["site_label"].strip())
    cls = pair_class.get(sp, "NO_REPLICATE_PAIR")
    reps = repl_codon.get(key, [])
    if reps:
        agree = sum(1 for r in reps if truthy(r["codon_agree"]))
        repl_note = "%d/%d replicate calls agree on the exact codon" % (agree, len(reps))
    else:
        repl_note = "no replicate call at this site"

    if cls == "INDEPENDENT" and reps:
        verdict = "CONFIRMED on a genuinely independent assembly"
    elif cls == "TWIN":
        verdict = "NOT CONFIRMED - replicate is a GCA/GCF twin (identical DNA)"
    elif cls == "VERSION_BUMP":
        verdict = "WEAK - replicate is a version bump of the same assembly"
    else:
        verdict = "UNCONFIRMED - no second assembly available"

    print("  %-26s %-6s %-6s %s->%s  codon=%s" % (sp, e["gene"], e["site_label"],
                                                  e["wild_type_aa"],
                                                  e["resistance_allele"], e["codon"]))
    print("      replicate class : %s" % cls)
    print("      replication     : %s" % repl_note)
    print("      VERDICT         : %s" % verdict)
    print("      2-step paths    : %s" % e.get("n_2step", "?"))
    print("      primary ref     : %s" % e.get("primary_ref", "")[:80])
    print()

    breach_rows.append({
        "species": sp, "gene": e["gene"], "site_label": e["site_label"],
        "wild_type_aa": e["wild_type_aa"], "resistance_allele": e["resistance_allele"],
        "codon": e["codon"], "n_2step": e.get("n_2step", ""),
        "site_supply_fraction": e["site_supply_fraction"],
        "replicate_pair_class": cls, "replication_note": repl_note,
        "verdict": verdict, "primary_ref": e.get("primary_ref", ""),
    })

confirmed = sum(1 for r in breach_rows if r["verdict"].startswith("CONFIRMED"))
print("  Of %d breaches, %d are confirmed on a genuinely independent assembly." % (k, confirmed))
if confirmed < k:
    print("  The remaining %d cannot be distinguished from annotation error." % (k - confirmed))
    print("  Report the breach rate as a range: %d to %d of %d (%.2f%% to %.2f%%)."
          % (confirmed, k, n, 100.0 * confirmed / n, 100.0 * k / n))
print()

# ======================================================================
# PART C -- how many alleles were reachable, descriptively
# ======================================================================
alleles_at_site = defaultdict(set)
site_rows = defaultdict(dict)
with open(ATLAS, newline="") as fh:
    for r in csv.DictReader(fh):
        if not truthy(r.get("is_wildtype")):
            continue
        g, s = r["gene"].strip(), r["site_label"].strip()
        a = (r.get("resistance_allele") or "").strip()
        if a:
            alleles_at_site[(g, s)].add(a)
        site_rows[(r["species"].strip(), g, s)][a] = r

print("-" * 74)
print("PART C   Choice set size at the sites where events occurred")
print("-" * 74)
print("For each event: of the documented alleles at that site, how many were")
print("one step from this species' codon, and how many needed two?")
print()

choice = []
for e in events:
    key = (e["species"].strip(), e["gene"].strip(), e["site_label"].strip())
    rows = site_rows.get(key, {})
    if not rows:
        continue
    one = sum(1 for a, r in rows.items() if truthy(r.get("one_step_accessible")))
    two = sum(1 for a, r in rows.items()
              if not truthy(r.get("one_step_accessible")) and num(r.get("n_2step_paths")) > 0)
    choice.append((e, one, two, len(rows)))

if choice:
    dist = Counter(one for _, one, _, _ in choice)
    print("  events by number of ONE-STEP-reachable documented alleles:")
    for kk in sorted(dist):
        print("    %d reachable in one step : %3d events" % (kk, dist[kk]))
    multi = [c for c in choice if c[1] >= 2]
    print()
    print("  Events where TWO OR MORE alleles were one step away: %d" % len(multi))
    print("  These are the events in which accessibility alone cannot explain the")
    print("  outcome, so they are the ones that can test what breaks the tie")
    print("  (tolerability, transition bias, population size). That is the")
    print("  regression Step 9 needs, and it has n = %d, not zero." % len(multi))
print()

# ======================================================================
# write outputs
# ======================================================================
os.makedirs(RES, exist_ok=True)

with open(OUT_EST, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["quantity", "n", "observed", "expected", "ratio", "exact_p", "note"])
    w.writerow(["stratum_p_eq_1", len(S1), sum(1 for e in S1 if e["_acc"]),
                round(sum(e["_p"] for e in S1), 4), "", "",
                "correct by construction"])
    w.writerow(["stratum_p_eq_0", len(S0), sum(1 for e in S0 if e["_acc"]),
                round(sum(e["_p"] for e in S0), 4), "", "",
                "correct by construction; these are the breaches"])
    w.writerow(["conditional_informative", len(SI), obs_inf, round(exp_inf, 4),
                round(obs_inf / exp_inf, 4), "%.6g" % p_inf,
                "SATURATED: ratio = 1/mean(p) = %.4f, cannot be exceeded" % (1 / mean_p_inf)])
    w.writerow(["unconditional_all_events", len(events), obs_all, round(exp_all, 4),
                round(obs_all / exp_all, 4), "%.6g" % p_all,
                "same exact P; degenerate strata carry no likelihood"])
    w.writerow(["breach_rate", n, k, "", round(k / n, 6), "",
                "95%% CI %.4f to %.4f (Clopper-Pearson)" % (lo, hi)])
    w.writerow(["breach_rate_confirmed_only", n, confirmed, "",
                round(confirmed / n, 6), "",
                "restricted to breaches confirmed on an independent assembly"])

with open(OUT_BRE, "w", newline="") as fh:
    if breach_rows:
        w = csv.DictWriter(fh, fieldnames=list(breach_rows[0].keys()))
        w.writeheader()
        w.writerows(breach_rows)

rate_rows = [{"scope": "ALL", "key": "all_events", "n": n, "breaches": k,
              "rate": round(k / n, 6), "ci_low": round(lo, 6), "ci_high": round(hi, 6)}]
for field, label in (("gene", "gene"), ("site_label", "site")):
    groups = defaultdict(list)
    for e in events:
        groups[e[field].strip()].append(e)
    for g, evs in sorted(groups.items()):
        kk = sum(1 for e in evs if not e["_acc"])
        l, h = clopper_pearson(kk, len(evs))
        rate_rows.append({"scope": label, "key": g, "n": len(evs), "breaches": kk,
                          "rate": round(kk / len(evs), 6),
                          "ci_low": round(l, 6), "ci_high": round(h, 6)})
with open(OUT_RATE, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["scope", "key", "n", "breaches", "rate",
                                       "ci_low", "ci_high"])
    w.writeheader()
    w.writerows(rate_rows)

with open(OUT_MD, "w") as fh:
    fh.write("# Step 2: estimands and the breach rate\n\n")
    fh.write("Generated by `src/step2_01_breach_rate_and_estimands.py`.\n\n")
    fh.write("## The partition\n\n")
    fh.write("| Stratum | n | Observed accessible | Expected |\n|---|---|---|---|\n")
    fh.write("| p = 1 | %d | %d | %.2f |\n" % (len(S1), sum(1 for e in S1 if e["_acc"]),
                                               sum(e["_p"] for e in S1)))
    fh.write("| p = 0 | %d | %d | %.2f |\n" % (len(S0), sum(1 for e in S0 if e["_acc"]),
                                               sum(e["_p"] for e in S0)))
    fh.write("| 0 < p < 1 | %d | %d | %.2f |\n" % (len(SI), obs_inf, exp_inf))
    fh.write("| **All** | **%d** | **%d** | **%.2f** |\n\n" % (len(events), obs_all, exp_all))
    fh.write("The null is correct by construction on the degenerate strata. "
             "All evidence lies in the %d informative events.\n\n" % len(SI))
    fh.write("## Two estimands\n\n")
    fh.write("**Conditional (confirmatory).** %d/%d accessible; ratio %.4f; exact "
             "*P* = %.3g. The ratio is **saturated**: with every informative event "
             "accessible it equals 1/mean(p) = %.4f exactly and cannot be exceeded, "
             "so the evidence lies in the exact tail and not in the magnitude.\n\n"
             % (obs_inf, len(SI), obs_inf / exp_inf, p_inf, 1 / mean_p_inf))
    fh.write("**Unconditional.** %d/%d accessible (%.1f%%); ratio %.4f; the same "
             "exact *P*, because the degenerate strata contribute no likelihood.\n\n"
             % (obs_all, len(events), 100.0 * obs_all / len(events), obs_all / exp_all))
    fh.write("## The breach rate\n\n")
    fh.write("The gate was crossed in **%d of %d events, %.2f%% "
             "(95%% CI %.2f%%-%.2f%%)**.\n\n" % (k, n, 100.0 * k / n,
                                                 100.0 * lo, 100.0 * hi))
    fh.write("| Species | Gene | Site | Substitution | Codon | Replicate class | Verdict |\n")
    fh.write("|---|---|---|---|---|---|---|\n")
    for r in breach_rows:
        fh.write("| *%s* | %s | %s | %s->%s | %s | %s | %s |\n"
                 % (r["species"], r["gene"], r["site_label"], r["wild_type_aa"],
                    r["resistance_allele"], r["codon"], r["replicate_pair_class"],
                    r["verdict"]))
    fh.write("\nOf %d breaches, %d are confirmed on a genuinely independent "
             "assembly.\n" % (k, confirmed))

print("=" * 74)
print("WROTE")
for p in (OUT_EST, OUT_BRE, OUT_RATE, OUT_MD):
    print("  %s" % os.path.relpath(p, ROOT))
print("=" * 74)
