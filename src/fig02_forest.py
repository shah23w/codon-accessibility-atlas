#!/usr/bin/env python3
"""
fig02_forest.py - FIGURE 2.  The effect, site by site and choice by choice.

WHY THIS FIGURE IS NOT A CONVENTIONAL FOREST PLOT
-------------------------------------------------
A forest plot shows point estimates and the heterogeneity among them. Tabulating
the 47 site strata shows that 15 are informative for the Mantel-Haenszel
estimate, and that at nearly all of those 15 the cell b - events landing on a
codon that could NOT reach the resistance allele in one step - is EMPTY.

That is the paper's finding stated in its rawest form: at almost every site,
every documented resistance event used an allele the species' own codon could
already reach. But it means the exact conditional odds ratio is unbounded above
at nearly every site, so a forest of point estimates would be a row of arrows
leaving the panel - a figure whose visual content is an artefact of a zero cell
rather than a comparison the reader can make.

The exact count is COMPUTED and printed at run time rather than stated here,
because a number written into a docstring goes stale silently: this one said
"13 of those 15" while the script had been printing 12 of 15.

So panel a plots the comparison that does NOT degenerate, on a probability axis:

    grey mark   the accessible fraction of the panel at that site, (a+c)/n.
                This is a census of every wild-type species in the atlas at that
                site, not a sample, so it carries no interval.
    blue mark   the fraction of OBSERVED events at that site that landed on an
                accessible codon, a/(a+b), with an exact Clopper-Pearson interval.
                Small event counts give wide intervals, and they are drawn wide.

The distance between the two marks is the effect at that site, and the sites are
ordered within each gene by the grey mark, so the reader can see that the blue
marks sit at or near 1.0 whether the background was 0.9 or 0.4. The exact
conditional odds ratio and its interval are still printed for every site, in the
right-hand column and in full in results/fig02_per_site_or.csv - the forest
information is kept, it is simply not asked to carry the picture.

Panel b answers the second fair objection - that the pooled number is an artefact
of one analytic choice - by drawing every control from Steps 10 to 14 on one log
axis, against a band showing the primary estimate.

CHECKS
------
The pooled Mantel-Haenszel estimate is recomputed here from the same two inputs
Step 11 used, and is then compared against results/rev11_stratification.csv. A
mismatch is fatal, not cosmetic: a figure that quietly disagrees with the text is
worse than no figure. The per-site intervals are exact conditional (Cornfield)
intervals; no continuity correction is applied anywhere, because a 0.5 added to
an empty cell is a number nobody measured.

Reads   data/processed/codon_atlas.csv
        results/core_test_events*.csv          (via the loud resolver)
        results/rev10_comparison.csv
        results/rev11_stratification.csv
        results/rev12_adjusted.csv
        results/rev13_sensitivity.csv
        results/rev14_altexon_sensitivity.csv
Writes  results/fig02_per_site_or.csv
        results/figures/figure2_forest.pdf
        results/figures/figure2_forest.png
"""
from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402
from matplotlib.lines import Line2D                                  # noqa: E402
from scipy.stats import beta as _beta                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PROJ = os.path.dirname(HERE)
RES = os.path.join(PROJ, "results")
FIGDIR = os.path.join(RES, "figures")
ATLAS = os.path.join(PROJ, "data", "processed", "codon_atlas.csv")

import fig00_style as S                                              # noqa: E402
from rev09_01_two_step_routes import resolve_results_file            # noqa: E402
from rev10_01_background_control import mantel_haenszel              # noqa: E402

ALPHA = 0.05
# Vertical offset of the Clopper-Pearson interval from the row's own baseline,
# in row units. Small enough to read as the same row, large enough that the
# interval and the background-to-observed connector do not share a line.
CP_OFFSET = 0.20
XLIM_OR = (0.04, 400.0)


# ---------------------------------------------------------------------------
# exact conditional inference for a 2x2 table
# ---------------------------------------------------------------------------
def _support(a, b, c, d):
    m1, n1, N = a + b, a + c, a + b + c + d
    return max(0, m1 - (N - n1)), min(m1, n1), m1, n1, N


def _logw(k, m1, n1, N):
    """log of the unnormalised conditional weight at psi = 1."""
    return (math.lgamma(n1 + 1) - math.lgamma(k + 1) - math.lgamma(n1 - k + 1)
            + math.lgamma(N - n1 + 1) - math.lgamma(m1 - k + 1)
            - math.lgamma(N - n1 - m1 + k + 1))


def _probs(t, ks, base):
    lw = [b + k * t for k, b in zip(ks, base)]
    mx = max(lw)
    w = [math.exp(x - mx) for x in lw]
    s = sum(w)
    return [x / s for x in w]


def _solve(fn, lo=-45.0, hi=45.0, tol=1e-10, _widen=6):
    """
    Bisection on log(OR); fn must be monotone with a sign change in [lo, hi].

    The bracket widens rather than giving up. Returning None on a failed
    bracket was indistinguishable, downstream, from a genuinely open interval:
    exact_or turns None into 0.0 or inf, so a numerical failure and a true
    one-sided bound printed identically. Widening first means None now means
    what it says, and the caller checks whether the cell count is interior
    before accepting it.
    """
    flo, fhi = fn(lo), fn(hi)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    tries = 0
    while flo * fhi > 0 and tries < _widen:
        lo, hi = lo * 2.0, hi * 2.0
        flo, fhi = fn(lo), fn(hi)
        tries += 1
    if flo * fhi > 0:
        return None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        fm = fn(mid)
        if fm == 0.0 or (hi - lo) < tol:
            return mid
        if flo * fm < 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def exact_or(a, b, c, d, alpha=ALPHA):
    """
    Conditional MLE of the odds ratio and its exact (Cornfield) interval.

    Open bounds come back as 0.0 or math.inf. Nothing is continuity-corrected.
    """
    k_lo, k_hi, m1, n1, N = _support(a, b, c, d)
    if k_hi < k_lo or m1 == 0 or n1 == 0 or n1 == N or m1 == N:
        return float("nan"), float("nan"), float("nan")
    ks = list(range(k_lo, k_hi + 1))
    base = [_logw(k, m1, n1, N) for k in ks]

    if a >= k_hi:
        est = float("inf")
    elif a <= k_lo:
        est = 0.0
    else:
        t = _solve(lambda t: sum(k * p for k, p in zip(ks, _probs(t, ks, base))) - a)
        est = math.exp(t) if t is not None else float("nan")

    # An interior cell count has a finite bound on both sides, so a failure to
    # bracket there is a numerical problem, not an open interval. Silently
    # substituting 0 or infinity made the two look identical in the figure.
    if a <= k_lo:
        lo_b = 0.0
    else:
        t = _solve(lambda t: sum(p for k, p in zip(ks, _probs(t, ks, base))
                                 if k >= a) - alpha / 2)
        if t is None:
            raise SystemExit(
                "FATAL: the lower exact bound did not bracket for an interior "
                "table (a=%d, b=%d, c=%d, d=%d). This is a numerical failure, "
                "not a one-sided interval; it must not be reported as 0."
                % (a, b, c, d))
        lo_b = math.exp(t)
    if a >= k_hi:
        hi_b = float("inf")
    else:
        t = _solve(lambda t: sum(p for k, p in zip(ks, _probs(t, ks, base))
                                 if k <= a) - alpha / 2)
        if t is None:
            raise SystemExit(
                "FATAL: the upper exact bound did not bracket for an interior "
                "table (a=%d, b=%d, c=%d, d=%d). This is a numerical failure, "
                "not a one-sided interval; it must not be reported as "
                "infinite." % (a, b, c, d))
        hi_b = math.exp(t)
    return est, lo_b, hi_b


def clopper_pearson(k, n, alpha=ALPHA):
    """Exact binomial interval. Degenerate ends return a hard 0 or 1, not a nudge."""
    if n == 0:
        return float("nan"), float("nan")
    lo = 0.0 if k == 0 else _beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else _beta.ppf(1 - alpha / 2, k + 1, n - k)
    return float(lo), float(hi)


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load_rows():
    """Rebuild exactly the rows Step 11 stratified - same filters, same keys."""
    ev_path = resolve_results_file("core_test_events")
    observed = set()
    with open(ev_path) as fh:
        for r in csv.DictReader(fh):
            observed.add((r["species"], r["gene"], r["site_label"], r["resistance_allele"]))
    rows = []
    with open(ATLAS) as fh:
        for r in csv.DictReader(fh):
            if r["is_wildtype"] != "True" or r["already_resistant"] == "True":
                continue
            rows.append({
                "species": r["species"], "gene": r["gene"],
                "site_label": r["site_label"],
                "site": "%s:%s" % (r["gene"], r["site_label"]),
                "accessible": r["one_step_accessible"] == "True",
                "used": (r["species"], r["gene"], r["site_label"],
                         r["resistance_allele"]) in observed,
            })
    print("  events file : %s" % os.path.basename(ev_path))
    print("  atlas rows  : %d  (wild-type, not already resistant)" % len(rows))
    return rows


def tabulate(rows):
    cells = defaultdict(lambda: [0, 0, 0, 0])
    species = defaultdict(set)
    meta = {}
    for r in rows:
        idx = (0 if r["used"] else 2) + (0 if r["accessible"] else 1)
        cells[r["site"]][idx] += 1
        species[r["site"]].add(r["species"])
        meta[r["site"]] = (r["gene"], r["site_label"])
    return cells, species, meta


def informative(a, b, c, d):
    """A stratum informs the MH estimate when r = ad/n or s = bc/n is non-zero."""
    return (a + b + c + d) > 0 and not (a * d == 0 and b * c == 0)


def read_table(name):
    with open(os.path.join(RES, name)) as fh:
        return list(csv.DictReader(fh))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("FIGURE 2.  Per-site accessibility and the specification ladder")
    print("=" * 78)
    print()

    rows = load_rows()
    cells, species, meta = tabulate(rows)

    pooled = mantel_haenszel([tuple(v) for v in cells.values()])
    print("\n  recomputed pooled MH : OR %.3f [%.3f, %.3f]  (%d informative of %d strata)"
          % (pooled[0], pooled[1], pooled[2], pooled[3], len(cells)))

    ref = {r["stratification"]: r for r in read_table("rev11_stratification.csv")}["site"]
    want = (fnum(ref["odds_ratio"]), fnum(ref["ci_lo"]), fnum(ref["ci_hi"]),
            int(ref["n_informative"]))
    print("  Step 11 recorded     : OR %.3f [%.3f, %.3f]  (%d informative)" % want)
    for got, exp, lab in zip(pooled, want, ("OR", "CI lower", "CI upper", "n informative")):
        tol = 0.01 if lab != "n informative" else 0
        if abs(got - exp) > tol:
            sys.exit("FATAL: figure disagrees with results/rev11_stratification.csv on %s "
                     "(%.4f vs %.4f). The figure is not plotting the published analysis."
                     % (lab, got, exp))
    print("  CHECK PASSED - the figure plots the same analysis the text reports.\n")

    # ---- per-site table --------------------------------------------------
    per_site = []
    for site, (a, b, c, d) in cells.items():
        if not informative(a, b, c, d):
            continue
        gene, label = meta[site]
        n = a + b + c + d
        est, lo, hi = exact_or(a, b, c, d)
        p_obs = a / (a + b)
        cp_lo, cp_hi = clopper_pearson(a, a + b)
        per_site.append({
            "site": site, "gene": gene, "site_label": label,
            "a": a, "b": b, "c": c, "d": d,
            "n_rows": n, "n_species": len(species[site]), "n_events": a + b,
            "background_accessible": (a + c) / n,
            "observed_accessible": p_obs,
            "observed_lo": cp_lo, "observed_hi": cp_hi,
            "or_exact": est, "or_lo": lo, "or_hi": hi,
            "mh_weight": (a * d / n) if n else 0.0,
        })
    # Two different senses of "informs the estimate" meet in this panel and
    # the legend keeps them apart. A site is PLOTTED when its 2x2
    # table has a within-site contrast, which is what the Mantel-Haenszel
    # estimator needs. A site admits a per-site EXACT test only when at least
    # one of its events has a null probability strictly between 0 and 1. The
    # second is the stricter condition, so some plotted sites carry no test.
    _tested = {r["statistic"][5:] for r in read_table("corrected_statistics.csv")
               if r["statistic"].startswith("site_")}
    _plotted_keys = {"%s_%s" % (r["gene"], r["site_label"]) for r in per_site}
    _plot_no_test = sorted(_plotted_keys - _tested)
    _ev_no_test = sum(r["n_events"] for r in per_site
                      if "%s_%s" % (r["gene"], r["site_label"]) in _plot_no_test)
    _unplotted = [s for s in cells if s not in {r["site"] for r in per_site}]
    _unplot_with_ev = sum(1 for s in _unplotted if (cells[s][0] + cells[s][1]) > 0)
    _unplot_no_ev = len(_unplotted) - _unplot_with_ev
    if len(_plotted_keys) - len(_plotted_keys & _tested) != len(_plot_no_test):
        sys.exit("FATAL: the plotted and tested site sets do not reconcile.")
    print("  plotted sites with no per-site exact test : %d (%s), %d events"
          % (len(_plot_no_test), ", ".join(_plot_no_test), _ev_no_test))
    print("  unplotted sites: %d carry events, %d carry none"
          % (_unplot_with_ev, _unplot_no_ev))

    n_uninf = len(cells) - len(per_site)
    print("  informative sites : %d      no within-site contrast (not plotted) : %d"
          % (len(per_site), n_uninf))
    n_zero_b = sum(1 for r in per_site if r["b"] == 0)
    print("  sites where NO event used an inaccessible codon (b = 0) : %d of %d"
          % (n_zero_b, len(per_site)))

    with open(os.path.join(RES, "fig02_per_site_or.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_site[0].keys()))
        w.writeheader()
        for r in sorted(per_site, key=lambda r: (r["gene"], r["site"])):
            w.writerow(r)
    print("  wrote results/fig02_per_site_or.csv")

    # Grouped by GENE, not by gene label. GluCl_Px and GluCl_Tu are two sites
    # on one GluCl gene (ortho01: the same protein accession in 236 of 569
    # species), so using the raw label as a group heading printed "GluCl_Tu" as
    # if it were a gene name and would have shown seven genes where the paper,
    # canonical_quantities.csv and Figures 1, 4 and S1 all say six.
    by_gene = defaultdict(list)
    for r in per_site:
        g = "GluCl" if r["gene"].startswith("GluCl") else r["gene"]
        by_gene[g].append(r)
    gene_order = sorted(by_gene, key=lambda g: (-len(by_gene[g]), g))

    # ---- panel b: the specification ladder -------------------------------
    strat = {r["stratification"]: r for r in read_table("rev11_stratification.csv")}
    adj = {r["adjustment"]: r for r in read_table("rev12_adjusted.csv")}
    par = {r["restriction"]: r for r in read_table("rev13_sensitivity.csv")}
    alt = {r["restriction"]: r for r in read_table("rev14_altexon_sensitivity.csv")}
    arm = {r["arm"]: r for r in read_table("rev10_comparison.csv")}
    armN = arm["ARM N  natural substitutions"]

    # The ascertainment control belongs on this axis more than anything else
    # on it. Optional, because the figure must still build before APRD has
    # been run.
    aprd = None
    try:
        rows = read_table("aprd02_restricted_odds_ratio.csv")
        aprd = {r["set"]: r for r in rows}.get("APRD phenotypically resistant")
    except SystemExit:
        pass
    except FileNotFoundError:
        pass

    def spec(label, row, colour=S.C_OBS, primary=False):
        return {"label": label, "est": fnum(row.get("odds_ratio", row.get("mh_or"))),
                "lo": fnum(row["ci_lo"]), "hi": fnum(row["ci_hi"]),
                "colour": colour, "primary": primary}

    ladder = [
        ("Stratification", [
            spec("by site", strat["site"], primary=True),
            spec("by order × site", strat["order x site"]),
            spec("by family × site", strat["family x site"]),
        ]),
        ("Substitution tolerability", [
            spec("× Grantham tertile", adj["site x Grantham tertile"]),
            spec("× ESM-2 tertile", adj["site x ESM-2 tertile (subset)"]),
            spec("× panel-frequency tertile", adj["site x panel-frequency tertile"]),
            spec("× Grantham × panel frequency",
                 adj["site x Grantham x panel-frequency"]),
        ]),
        ("Paralogy and site mapping", [
            spec("flagged species × gene dropped",
                 par["dropping every FLAGGED species x gene"]),
            spec("site context identity ≥ 0.50", par["site context identity >= 0.50"]),
            spec("site context identity ≥ 0.70", par["site context identity >= 0.70"]),
            spec("single mapped candidate",
                 par["only ONE candidate mapped at that site"]),
        ]),
        ("Alternative exons", [
            spec("four alt-exon sites excluded",
                 alt["excluding the four alt-exon sites"]),
        ]),
        # Not a robustness check on the same quantity: it changes the unit of
        # analysis, and the figure says so by giving it its own heading. The
        # estimand everywhere else on this axis is one row per species, site
        # and candidate allele. Collapsing to one row per species and site
        # merges the candidate alleles, which is precisely the contrast the
        # hypothesis is about, and at a five-allele site almost every species
        # reaches at least one allele, so the collapsed table is close to
        # degenerate by construction. It is drawn because the reader is
        # entitled to see how much of the estimate rests on the expansion.
        ("Unit of analysis", ([
            # Its own colour, not the contrast arm's. Both were orange, so
            # the two rows that are not specifications of the resistance
            # estimate were indistinguishable from each other, and the legend
            # could only name one of them. Aqua is directly labelled, as
            # fig00_style requires of it.
            spec("collapsed to one row per species × site",
                 par["collapsed to one row per species x site"],
                 colour=S.C_THIRD),
        ] if "collapsed to one row per species x site" in par else [])),
        ("Ascertainment", ([
            {"label": "phenotypically resistant species only",
             "est": fnum(aprd["odds_ratio"]), "lo": fnum(aprd["ci_lo"]),
             "hi": fnum(aprd["ci_hi"]), "colour": S.C_OBS, "primary": False},
        ] if aprd else [])),
        ("Contrast arm", [
            {"label": "natural substitutions, same estimator",
             "est": fnum(armN["mh_or"]),
             "lo": fnum(armN["ci_lo"]), "hi": fnum(armN["ci_hi"]),
             "colour": S.C_ALT, "primary": False},
        ]),
    ]

    # ---- row plans -------------------------------------------------------
    rows_a = []
    for g in gene_order:
        rows_a.append(("head", g))
        for r in sorted(by_gene[g], key=lambda r: r["background_accessible"]):
            rows_a.append(("site", r))

    # ---- the pooled result, as a summary row ------------------------------
    # The 15 rows above are the sites informative for the Mantel-Haenszel
    # estimate. They do NOT carry the paper's headline comparison, which is
    # unconditional and runs over all 108 curated events: 103 of them used an
    # allele the species' own codon could already reach, against 84.69 expected
    # under the atlas-pool null. That number had no home in the five-figure
    # lineup. It belongs here, on this axis, because this axis is already a
    # percentage of alleles on an accessible codon -- and a forest plot's own
    # idiom for a pooled estimate is a summary marker at the foot of the panel.
    est = {r["quantity"]: r for r in csv.DictReader(open(
        os.path.join(RES, "step2_estimands.csv")))}
    pool = est["unconditional_all_events"]
    n_pool = int(float(pool["n"]))
    obs_pool = int(float(pool["observed"]))
    exp_pool = float(pool["expected"])
    p_pool = float(pool["exact_p"])

    # Clopper-Pearson on the observed count, derived from the counts read above
    # rather than written down. No continuity correction, as everywhere else.
    cp_lo = float(_beta.ppf(ALPHA / 2, obs_pool, n_pool - obs_pool + 1)) \
        if obs_pool else 0.0
    cp_hi = float(_beta.ppf(1 - ALPHA / 2, obs_pool + 1, n_pool - obs_pool)) \
        if obs_pool < n_pool else 1.0

    # Only the two keys this row depends on are parsed as numbers: several
    # canonical values are deliberately prose (profile intervals, provenance
    # notes), so a blanket float() over the file fails.
    canon = {r["quantity"]: r["value"] for r in csv.DictReader(open(
        os.path.join(RES, "canonical_quantities.csv")))}
    for got, key, what in ((n_pool, "n_events_total", "event count"),
                           (exp_pool, "expected_null_pooled", "null expectation")):
        if abs(got - float(canon[key])) > 1e-3:
            sys.exit("FATAL: the summary row's %s is %s, but "
                     "canonical_quantities.csv declares %s (%s)"
                     % (what, got, canon[key], key))

    rows_a.append(("gap", None))
    rows_a.append(("pooled", {
        "label": "All %d curated events" % n_pool,
        "observed": obs_pool / n_pool,
        "expected": exp_pool / n_pool,
        "lo": cp_lo, "hi": cp_hi,
        "events": "%d / %d" % (obs_pool, n_pool),
        "ratio": float(pool["ratio"]),
        "p": p_pool,
    }))
    rows_b = []
    for title, items in ladder:
        if not items:                       # a control that has not been run yet
            continue
        rows_b.append(("head", title))
        for it in items:
            rows_b.append(("spec", it))

    nA, nB = len(rows_a), len(rows_b)

    # ---- layout, in explicit figure coordinates --------------------------
    # Gridspec fractions are hard to reason about when text columns hang off the
    # side of an axes, so both panels are placed by hand and every column has a
    # measured budget. Widths are fractions of W2 = 183 mm.
    S.apply()
    # 0.115, not 0.120: the summary row and its spacer added two rows, which at
    # the old height pushed the figure to 173.7 mm against NEE's 170 mm
    # ceiling. Trimming 4% off the row pitch keeps BOTH panels in the main
    # figure rather than exiling the specification ladder to Extended Data,
    # which is the one thing in this figure that answers "is the result an
    # artefact of one analytic choice".
    # 0.110 and a 0.50 in gap, not 0.115 and 0.60: adding the unit-of-analysis
    # arm to panel b added a heading and a row, which took the figure to
    # 174 mm against NEE's 170 mm ceiling. Both panels stay in the main figure,
    # because panel b is the one thing here that answers "is the result an
    # artefact of one analytic choice".
    ROW_H = 0.110                      # inches per row
    GAP = 0.50                         # inches between the panels
    TOP_PAD, BOT_PAD = 0.60, 0.46      # TOP_PAD carries panel a's legend band
    hA, hB = nA * ROW_H, nB * ROW_H
    height = TOP_PAD + hA + GAP + hB + BOT_PAD
    if height > S.H_MAX:
        sys.exit("FATAL: figure would be %.2f in (%.0f mm) tall, above the "
                 "%.2f in (%.0f mm) page ceiling, with %d rows. ROW_H would "
                 "have to fall to %.4f in, or a row has to go."
                 % (height, height * 25.4, S.H_MAX, S.H_MAX * 25.4, nA + nB,
                    (S.H_MAX - TOP_PAD - GAP - BOT_PAD) / (nA + nB)))
    print("  canvas : %.2f x %.2f in  (%.0f x %.0f mm)"
          % (S.W2, height, S.W2 * 25.4, height * 25.4))

    fig = plt.figure(figsize=(S.W2, height))
    yA = 1.0 - (TOP_PAD + hA) / height
    yB = BOT_PAD / height
    axA = fig.add_axes([0.088, yA, 0.400, hA / height])
    axB = fig.add_axes([0.262, yB, 0.455, hB / height])

    # ================= panel a =================
    axA.set_xlim(-0.02, 1.02)
    axA.set_ylim(nA - 0.5, -0.5)
    axA.set_yticks([])
    axA.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    axA.set_xticklabels(["0", "25", "50", "75", "100"])
    axA.set_xlabel("Resistance alleles on a one-step-accessible codon (%)")
    axA.spines["left"].set_visible(False)
    axA.tick_params(axis="y", length=0)
    axA.set_axisbelow(True)
    # No gridlines. These four ran at 25/50/75/100% straight through every
    # interval in the panel, and the specification says to avoid them. The
    # x tick marks carry the scale.
    axA.tick_params(axis="x", length=2.6, width=0.6, color=S.MUTED)

    tr = axA.get_yaxis_transform()
    # Columns start well clear of the axis. At 1.075 the widest value, "17 / 19",
    # reached back to 0.99 and sat on top of the 100% gridline, and the header
    # "Events on accessible" overran "Species" entirely.
    CX = {"ev": 1.14, "sp": 1.28, "or": 1.35}
    for i, (kind, payload) in enumerate(rows_a):
        if kind == "head":
            axA.text(-0.175, i, payload, transform=tr, ha="left", va="center",
                     fontsize=7, fontweight="bold", color=S.INK, clip_on=False)
            continue
        if kind == "gap":
            # A rule, not blank space: it says the row beneath is a different
            # kind of quantity from the site rows above, not another site.
            axA.plot([-0.175, CX["or"]], [i, i], transform=tr, color=S.RULE,
                     linewidth=0.6, clip_on=False, zorder=1)
            continue
        if kind == "pooled":
            p = payload
            axA.plot([p["expected"], p["observed"]], [i, i], color=S.MUTED,
                     linewidth=1.3, zorder=1, alpha=0.65, solid_capstyle="butt")
            axA.plot([p["lo"], p["hi"]], [i, i], color=S.C_POOLED,
                     linewidth=0.9, alpha=0.45, zorder=2, solid_capstyle="butt")
            axA.scatter([p["expected"]], [i], s=19, facecolor=S.SURFACE,
                        edgecolor=S.MUTED, linewidths=1.0, zorder=3, marker="o")
            # A diamond, not a circle: the standard mark for a pooled estimate,
            # and it distinguishes this row from the per-site rows at a glance.
            axA.scatter([p["observed"]], [i], s=34, color=S.C_POOLED, zorder=4,
                        linewidths=0, marker="D")
            axA.text(-0.150, i, p["label"], transform=tr, ha="left",
                     va="center", fontsize=6.5, fontweight="bold",
                     color=S.INK, clip_on=False)
            # The third column carries an EFFECT SIZE in every row: an exact
            # odds ratio for the sites, the enrichment ratio here. Not the P,
            # which would change the column's meaning mid-table and is in any
            # case unprintable here: mathtext renders a superscript at
            # ~70% of the base size, so at any size inside NEE's 5-7 pt band
            # the exponent of "4.6x10^-12" fell below the 5 pt floor (4.6 pt),
            # and Arial has no superscript glyphs beyond U+00B9/B2/B3, so a
            # literal "10" with superscript minus twelve would have printed as
            # empty boxes. The exact P belongs in the caption, where it is set
            # as text rather than drawn.
            S.text_column(axA, i,
                          [p["events"], "", "%.2f x expected" % p["ratio"]],
                          [CX["ev"], CX["sp"], CX["or"]],
                          ha=["right", "right", "left"], color=S.INK,
                          weight="bold")
            continue
        r = payload
        bg, ob = r["background_accessible"], r["observed_accessible"]
        # The gap between the two marks IS the effect at that site, so the
        # connector is drawn with more weight than the interval. At equal
        # weight a 1-of-1 site's exact interval, which spans almost the whole
        # axis, swallows the grey mark and reads as though the estimate itself
        # were spread across the panel.
        axA.plot([bg, ob], [i, i], color=S.MUTED, linewidth=1.1, zorder=1,
                 alpha=0.55, solid_capstyle="butt")
        # The interval sits on its own baseline a fifth of a row below the
        # connector. At twelve of the fifteen sites every observed event is
        # accessible, so the interval runs leftward from 1.0 along exactly the
        # row the grey connector occupies; drawn on one line the two could not
        # be told apart, and the upper cap disappeared under the filled mark.
        yi = i + CP_OFFSET
        axA.plot([r["observed_lo"], r["observed_hi"]], [yi, yi],
                 color=S.C_OBS, linewidth=0.7, alpha=0.38, zorder=2,
                 solid_capstyle="butt")
        # End caps, so a wide interval reads as an interval rather than a rule.
        for xe in (r["observed_lo"], r["observed_hi"]):
            axA.plot([xe, xe], [yi - 0.16, yi + 0.16], color=S.C_OBS,
                     linewidth=0.7, alpha=0.38, zorder=2)
        axA.scatter([bg], [i], s=17, facecolor=S.SURFACE, edgecolor=S.MUTED,
                    linewidths=0.9, zorder=3, marker="o")
        axA.scatter([ob], [i], s=19, color=S.C_OBS, zorder=4, linewidths=0,
                    marker="o")
        axA.text(-0.150, i, r["site_label"], transform=tr, ha="left", va="center",
                 fontsize=6.5, color=S.INK2, clip_on=False)
        # Counts only, with no third column for the odds ratio and its
        # interval. None of the NEE research Articles surveyed prints a
        # confidence interval inside a figure: the point and the whisker against
        # the log axis carry the estimate, and the numbers belong in the text
        # and in Dataset S6. The open-bounded rows read particularly badly as
        # text, since seven of them printed as "inf [0.03, inf]".
        S.text_column(axA, i,
                      ["%d / %d" % (r["a"], r["n_events"]),
                       "%d" % r["n_species"]],
                      [CX["ev"], CX["sp"]],
                      ha=["right", "right"], color=S.INK2)

    for x, lab, ha in ((CX["ev"], "Events\nacc. / all", "right"),
                       (CX["sp"], "Species", "right")):
        axA.text(x, -1.55, lab, transform=tr, ha=ha, va="center", fontsize=6,
                 color=S.MUTED, clip_on=False, linespacing=1.25)

    handles = [
        Line2D([], [], marker="o", linestyle="none", markerfacecolor=S.SURFACE,
               markeredgecolor=S.MUTED, markeredgewidth=0.9, markersize=3.7,
               label="Accessible fraction of the whole panel at that site"),
        Line2D([], [], marker="o", linestyle="none", color=S.C_OBS, markersize=3.9,
               label="Observed resistance events, with exact 95% interval"),
    ]
    axA.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.175, 1.055),
               frameon=False, handletextpad=0.45, borderpad=0.0,
               labelspacing=0.40, borderaxespad=0.0, fontsize=6.5,
               labelcolor=S.INK2)

    # The summary row is a different quantity from the rows above it, and its
    # third column is a P rather than an odds ratio. Say so on the figure: a
    # reader should not have to infer from the caption why one column changes
    # meaning in one row.
    # Placed in AXES fractions below the x-axis title, not in data
    # coordinates just under the last row, where it ran straight through the
    # tick labels and the axis title. One line, so it stays inside the 15 mm
    # gap between the two panels.
    # Two lines, not one: at 5.8 pt a single line of this ran past the right
    # edge of the page. Panel a's axes is only 0.400 of the figure wide, so the
    # note has to wrap. y = -0.15 of the axes height puts it clear of the tick
    # labels and the axis title while staying inside the 15 mm gap above panel
    # b -- a third line would not fit.
    # The summary row is explained in the legend rather than on the canvas,
    # which is where interpretation belongs.

    # ================= panel b =================
    S.forest_axis(axB, XLIM_OR, ticks=[0.1, 1, 10, 100])
    axB.set_ylim(nB - 0.5, -0.5)
    axB.set_yticks([])
    # The unit is one (species, site, candidate allele) row, and the axis says
    # so rather than saying "species-level". Accessibility is a property of the
    # codon-allele pair and not of the species: a species' codon can reach one
    # documented allele at a site and not another, which is the contrast the
    # hypothesis is about, so the atlas is long over alleles and so is this
    # odds ratio. "Species-level" would be wrong in a way that matters: at a
    # five-allele site each wild-type species contributes five rows, the
    # Species column counts it once, and a reader could not then recover the
    # plotted percentages from the two columns.
    # Figure 3's other odds-ratio family is coarser still, counting one row per
    # site and candidate residue with no species dimension at all.
    axB.set_xlabel("Odds ratio, resistance allele reachable in one step\n"
                   "(species x candidate allele at each site, log scale)")
    axB.axvspan(want[1], want[2], color=S.BLUE, alpha=0.07, zorder=0, linewidth=0)
    axB.axvline(want[0], color=S.BLUE, alpha=0.45, linewidth=0.7, zorder=0,
                linestyle=(0, (3, 2)))

    trB = axB.get_yaxis_transform()
    for i, (kind, payload) in enumerate(rows_b):
        if kind == "head":
            axB.text(-0.565, i, payload, transform=trB, ha="left", va="center",
                     fontsize=7, fontweight="bold", color=S.INK, clip_on=False)
            continue
        it = payload
        S.forest_row(axB, i, it["est"], it["lo"], it["hi"], XLIM_OR,
                     size=17, color=it["colour"])
        axB.text(-0.540, i, S.cap(it["label"]), transform=trB, ha="left", va="center",
                 fontsize=6.5, color=S.INK if it["primary"] else S.INK2,
                 fontweight="bold" if it["primary"] else "normal", clip_on=False)
    # No numeric column here either, for the same reason as panel a. The
    # specification names are the information; the estimates are in Dataset S6.

    S.fig_label(fig, "a", 0.006, 1.0 - 0.012)
    # Below panel a's footnote, not above it. At -0.030 the letter b sat level
    # with panel a's axis title and ABOVE the note explaining panel a's summary
    # row, so the note read as though it belonged to panel b.
    # Clear of the first group heading. At -0.075 the letter sat on the "S"
    # of "Stratification"; at -0.030 it sat level with panel a's axis title.
    # This puts it one row above the first group heading.
    S.fig_label(fig, "b", 0.006, yA - 0.048)

    S.save(fig, "figure2_forest", FIGDIR,
           expect_mm=(S.W2 * 25.4, height * 25.4))
    plt.close(fig)

    print("\n  CAPTION NUMBERS TO QUOTE")
    print("  ------------------------")
    print("  informative site strata plotted : %d across %d genes (%s)"
          % (len(per_site), len(gene_order), ", ".join(gene_order)))
    print("  strata with no within-site contrast, not plotted : %d" % n_uninf)
    print("  sites at which every event used an accessible codon : %d of %d"
          % (n_zero_b, len(per_site)))
    print("  pooled : OR %.2f [%.2f, %.2f]" % (pooled[0], pooled[1], pooled[2]))
    print("  contrast arm : OR %.2f [%.2f, %.2f]"
          % (fnum(armN["mh_or"]), fnum(armN["ci_lo"]), fnum(armN["ci_hi"])))
    _sd = []
    for r in sorted(per_site, key=lambda r: (r["gene"], r["site_label"])):
        _sd.append(("a", "%s %s background accessible" % (r["gene"], r["site_label"]),
                    r["background_accessible"], "", "",
                    "%d of %d rows" % (r["a"] + r["c"], r["n_rows"])))
        _sd.append(("a", "%s %s observed accessible" % (r["gene"], r["site_label"]),
                    r["observed_accessible"], r["observed_lo"], r["observed_hi"],
                    "%d of %d events; %d species"
                    % (r["a"], r["n_events"], r["n_species"])))
        _sd.append(("a", "%s %s exact OR" % (r["gene"], r["site_label"]),
                    r["or_exact"], r["or_lo"], r["or_hi"], "conditional exact"))
    _sd.append(("a", "pooled, all curated events", float(pool["ratio"]),
                cp_lo, cp_hi,
                "%d of %d accessible; %.4f expected; exact P %.6g"
                % (obs_pool, n_pool, exp_pool, p_pool)))
    for _k, _it in rows_b:
        if _k == "spec":
            _sd.append(("b", _it["label"], _it["est"], _it["lo"], _it["hi"],
                        "primary" if _it.get("primary") else ""))
    S.write_source_data("Fig2", _sd)

    S.write_legend_data("figure2_forest", {
        "n_sites_plotted": len(per_site),
        "n_sites_uninformative": n_uninf,
        "n_sites_all_accessible": n_zero_b,
        "pooled_mh_or": pooled[0],
        "pooled_mh_lo": pooled[1],
        "pooled_mh_hi": pooled[2],
        "armN_or": fnum(armN["mh_or"]),
        "armN_lo": fnum(armN["ci_lo"]),
        "armN_hi": fnum(armN["ci_hi"]),
        "n_events_pooled": n_pool,
        "n_events_accessible": obs_pool,
        "expected_accessible": exp_pool,
        "exact_p": p_pool,
        # The event arithmetic the plotted rows and the summary row imply, so
        # that the legend states it instead of leaving the reader to sum a
        # column and find 62 against a summary of 108.
        "n_events_at_plotted_sites": sum(r["n_events"] for r in per_site),
        "n_events_accessible_at_plotted_sites": sum(r["a"] for r in per_site),
        "n_events_at_unplotted_sites":
            n_pool - sum(r["n_events"] for r in per_site),
        "n_sites_with_exact_test": len(_plotted_keys & _tested),
        "n_sites_plotted_without_exact_test": len(_plot_no_test),
        "sites_plotted_without_exact_test": ", ".join(
            s.replace("_", " ", 1) for s in _plot_no_test),
        "n_events_at_plotted_sites_without_exact_test": _ev_no_test,
        "n_unplotted_sites_with_events": _unplot_with_ev,
        "n_unplotted_sites_without_events": _unplot_no_ev,
        "n_specifications": sum(1 for k, _ in rows_b if k == "spec"),
        # Counted on the row's own colour, one key per kind, so that
        # recolouring a row cannot silently move a count the legend quotes.
        "n_specifications_resistance":
            sum(1 for k, it in rows_b
                if k == "spec" and it["colour"] == S.C_OBS),
        "n_specifications_unit":
            sum(1 for k, it in rows_b
                if k == "spec" and it["colour"] == S.C_THIRD),
        "n_specifications_contrast":
            sum(1 for k, it in rows_b
                if k == "spec" and it["colour"] == S.C_ALT),
        "collapsed_or": next((it["est"] for k, it in rows_b if k == "spec"
                              and it["label"].startswith("collapsed")), None),
        "collapsed_lo": next((it["lo"] for k, it in rows_b if k == "spec"
                              and it["label"].startswith("collapsed")), None),
        "collapsed_hi": next((it["hi"] for k, it in rows_b if k == "spec"
                              and it["label"].startswith("collapsed")), None),
    })
    # A reader who sums the Events column must land on the summary row, so the
    # arithmetic is asserted rather than hoped for.
    _at = sum(r["n_events"] for r in per_site)
    if _at > n_pool:
        sys.exit("FATAL: the plotted sites carry %d events but the summary row "
                 "pools %d. One of the two is wrong." % (_at, n_pool))

    print("\n  PER-SITE ROWS")
    print("  %-14s %7s %7s %6s  %s" % ("site", "bg %", "obs %", "ev", "exact OR [95% CI]"))
    for r in sorted(per_site, key=lambda r: (r["gene"], r["site"])):
        print("  %-14s %6.1f%% %6.1f%% %3d/%-3d  %s"
              % (r["site"], 100 * r["background_accessible"],
                 100 * r["observed_accessible"], r["a"], r["n_events"],
                 S.fmt_or(r["or_exact"], r["or_lo"], r["or_hi"])))
    print("\nDone.")


if __name__ == "__main__":
    main()
