#!/usr/bin/env python3
"""
step8_04_figure2.py: Figure 2: necessary but not sufficient, and separable from
protein-level tolerability.

Panel A  both halves of the claim in one view. Accessibility raises the chance an
         allele is ever observed ~19-fold, and yet four in five accessible alleles
         are never observed.
Panel B  the independent fitness axis. ESM-2 sees only protein sequence - never a
         codon - and it does separate the classes, so the axis is informative.
Panel C  the question that matters: does accessibility survive adjustment for that
         axis? Odds ratios on a log scale, unadjusted and adjusted.

Every value is parsed from a result file or log at run time.
"""
from __future__ import annotations

import csv
import os
import re
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import DIR_FIGS, DIR_RESULTS, Tee, ensure_dirs  # noqa: E402
import fig00_style as st  # noqa: E402

TAGS = ("_primary", "_t1", "_t5", "_final", "")  # one order everywhere


def read_csv(p):
    try:
        with open(p, newline="", encoding="utf-8-sig") as fh:
            return list(csv.DictReader(fh))
    except OSError:
        return []


def first_csv(stem):
    for t in TAGS:
        p = os.path.join(DIR_RESULTS, f"{stem}{t}.csv")
        if os.path.exists(p):
            return read_csv(p)
    return []


def _tag(fig, ax, letter):
    """
    Panel letter in FIGURE coordinates, taken from the axes box.

    The previous call placed it at an axes-fraction dx of -0.44. Panel a's
    axes start at figure x = 0.055 and are about 0.19 wide, so that put the
    letter at figure x = -0.03 -- off the page, where it was silently clipped
    and the figure went out with no panel a at all. Anchoring to the axes box
    and clamping to the page means a layout change can move the letter but
    cannot delete it.
    """
    pos = ax.get_position()
    st.fig_label(fig, letter,
                 max(pos.x0 - 0.030, 0.005),
                 min(pos.y1 + 0.055, 0.980))


def main():
    ensure_dirs()
    st.apply()
    log = Tee("step8_04_figure2")
    log("=== Figure 3 ===")

    # Every value below comes from a results table, not from the run log.
    res_rows = {r["analysis"]: r for r in first_csv("multivariate_results")}
    spec_rows = first_csv("multivariate_specifications")
    if not res_rows or not spec_rows:
        raise SystemExit("multivariate results/specifications not found; run step6_03")

    # ---- Panel A: the 2x2 --------------------------------------------------
    a = int(float(res_rows["2x2_accessible_observed"]["estimate"]))
    na = int(float(res_rows["2x2_accessible_total"]["estimate"]))
    c = int(float(res_rows["2x2_inaccessible_observed"]["estimate"]))
    nc = int(float(res_rows["2x2_inaccessible_total"]["estimate"]))
    b, d = na - a, nc - c
    or2 = float(res_rows["2x2_necessary_sufficient"]["estimate"])
    log(f"  2x2: accessible {a}/{na}, inaccessible {c}/{nc}, OR {or2:.2f}")

    # ---- Panel B: ESM by class --------------------------------------------
    cls = first_csv("aim3_fitness_classes")
    groups = {}
    for r in cls:
        groups.setdefault(r["class"], []).append(float(r["esm_score"]))
    order = ["OBSERVED_ACCESSIBLE", "OBSERVED_INACCESSIBLE",
             "ACCESSIBLE_UNOBSERVED", "INACCESSIBLE"]
    pretty = {"OBSERVED_ACCESSIBLE": "Observed,\naccessible",
              "OBSERVED_INACCESSIBLE": "Observed,\ntwo-step",
              "ACCESSIBLE_UNOBSERVED": "Accessible,\nunobserved",
              "INACCESSIBLE": "Inaccessible,\nunobserved"}
    # n goes in the tick label; as free text inside the axes it landed on the whiskers
    log(f"  ESM classes: { {k: len(v) for k, v in groups.items()} }")

    # ---- Panel C: odds ratios ---------------------------------------------
    # The declared headline is not in multivariate_specifications; it is the
    # conditional logit with site fixed effects adjusted for ESM-2, and it
    # lives in step3_odds_ratio_estimators.csv. Read it from there, which is
    # the file canonical_quantities.csv names as its owner.
    _est = {r["estimator"]: r for r in first_csv("step3_odds_ratio_estimators")
            if r.get("estimator") and r["estimator"] != "estimator"}
    _HEAD = "conditional_logit_by_site_adj_esm"
    if _HEAD not in _est:
        raise SystemExit(
            "step3_odds_ratio_estimators.csv has no %r row. Figure 3 cannot "
            "draw the declared headline odds ratio." % _HEAD)

    def headline():
        r = _est[_HEAD]
        return (float(r["odds_ratio"]), float(r["ci_lo"]), float(r["ci_hi"]),
                float(r["p_value"]) if r.get("p_value") else float("nan"))

    def spec(label, term):
        r = next((x for x in spec_rows
                  if x["specification"] == label and x["term"] == term), None)
        return (float(r["odds_ratio"]), float(r["ci_lo"]), float(r["ci_hi"]),
                float(r["p_value"])) if r else None

    rows = [
        # Three lines, not two. On two lines the first read "Accessibility,
        # site fixed effects", 33 characters, and panel c's tick labels reached
        # 30.7 mm left of its own axes -- 11.9 mm inside panel b's box, so the
        # label sat on panel b and looked as though it belonged to it.
        ("Accessibility, site\nfixed effects\nadjusted for ESM", headline(),
         st.C_OBS),
        ("Accessibility\nunadjusted", spec("accessibility alone", "accessible"),
         st.C_OBS),
        ("Accessibility\nadjusted for ESM", spec("BOTH (the question)", "accessible"),
         st.C_OBS),
        ("ESM tolerability\nunadjusted", spec("tolerability alone", "esm_score"),
         st.C_ALT),
        ("ESM tolerability\nadjusted", spec("BOTH (the question)", "esm_score"),
         st.C_ALT),
    ]
    rows = [r for r in rows if r[1]]
    log(f"  ORs from multivariate_specifications: {len(rows)}")

    # -------------------------------------------------------------- figure
    # Three different forms: proportion columns, distributions, and an effect-size
    # forest. Panel titles are omitted; the legend and Results carry the interpretation.
    fig = plt.figure(figsize=(st.W2, 2.95))
    # left = 0.055 left only 10 mm for panel a's y tick labels, which now sit
    # outside the axes instead of under it. 0.105 gives 19 mm.
    gs = fig.add_gridspec(1, 3, width_ratios=[0.70, 1.22, 1.02], wspace=0.46,
                          left=0.105, right=0.985, top=0.80, bottom=0.22)

    # ---- A: the two proportions, with exact intervals ----------------------
    # A dot and interval on an axis that stops at 28%, not two bars against a
    # 100% ceiling: the second proportion is 1.1% and would draw as a sliver
    # one pixel high, so the contrast the panel exists to show would be the
    # one thing it could not display. This also shows how much better
    # determined the small proportion is.
    axA = fig.add_subplot(gs[0])
    from scipy.stats import beta as _beta
    obs_n, tot_n = [a, c], [a + b, c + d]
    pct = [100 * o / t for o, t in zip(obs_n, tot_n)]
    cis = []
    for o, t in zip(obs_n, tot_n):
        lo = 100 * float(_beta.ppf(0.025, o, t - o + 1)) if o else 0.0
        hi = 100 * float(_beta.ppf(0.975, o + 1, t - o)) if o < t else 100.0
        cis.append((lo, hi))
    yy = [1, 0]
    for i, y in enumerate(yy):
        col = st.C_OBS if i == 0 else st.MUTED
        axA.plot([cis[i][0], cis[i][1]], [y, y], color=col, lw=1.1,
                 solid_capstyle="round", zorder=2)
        axA.plot(pct[i], y, "o", ms=5.0, color=col, mec=st.SURFACE, mew=0.8,
                 zorder=3)
    axA.set_yticks(yy)
    # Short noun phrases only. The counts are in the legend and in the
    # source-data file; printed here they made the labels wider than the
    # margin and the page cut them.
    axA.set_yticklabels(["One\nmutation\naway", "Not one\nmutation\naway"],
                        fontsize=6.2)
    axA.set_ylim(-0.6, 1.6)
    axA.set_xlim(0, 28)
    axA.set_xticks([0, 10, 20])
    axA.set_xlabel("Ever observed as resistance (%)", fontsize=6.8)
    st.recessive_grid(axA, "x")
    axA.spines["left"].set_visible(False)
    axA.tick_params(axis="y", length=0)
    # No OR bracket here: panel c carries that odds ratio with its interval,
    # and a quantity that appears twice in one figure will eventually disagree
    # with itself.
    _tag(fig, axA, "a")

    # ---- B: score distributions -------------------------------------------
    axB = fig.add_subplot(gs[1])
    # Shifted left of where the gridspec puts it. The three panels are very
    # different shapes, so an even wspace leaves panel b crowding panel c's row
    # labels while a wide gap opens between panels a and b.
    _pb = axB.get_position()
    axB.set_position([_pb.x0 - 0.0355, _pb.y0, _pb.width, _pb.height])
    present = [k for k in order if k in groups]
    data = [groups[k] for k in present]
    labels = [f"{pretty[k]}\n(n = {len(groups[k])})" for k in present]
    cols = [st.C_OBS if k.startswith("OBSERVED") else st.MUTED for k in present]
    # a violin over six points is not a distribution; that class is drawn as its points
    SMALL = 12
    # not d: d is the inaccessible-and-never-observed cell of the 2x2 set
    # above, and must not be shadowed by a loop variable here.
    big = [i for i, vals in enumerate(data) if len(vals) >= SMALL]
    parts = axB.violinplot([data[i] for i in big], positions=big, widths=0.78,
                           showextrema=False, showmedians=False)
    for pc, i in zip(parts["bodies"], big):
        pc.set_facecolor(cols[i])
        pc.set_alpha(0.45 if cols[i] == st.C_OBS else 0.35)
        pc.set_edgecolor("none")
    bp = axB.boxplot([data[i] for i in big], positions=big, widths=0.13, vert=True,
                     showfliers=False, patch_artist=True,
                     medianprops=dict(color=st.SURFACE, lw=1.1),
                     whiskerprops=dict(color=st.INK2, lw=0.7),
                     capprops=dict(color=st.INK2, lw=0.7), boxprops=dict(lw=0))
    for patch, i in zip(bp["boxes"], big):
        patch.set_facecolor(cols[i])
    rng = __import__("random").Random(7)
    for i, vals in enumerate(data):
        if len(vals) >= SMALL:
            continue
        for v in vals:
            axB.plot(i + rng.uniform(-0.13, 0.13), v, "o", ms=3.4, color=cols[i],
                     mec=st.SURFACE, mew=0.5, zorder=3)
    axB.set_xticks(range(len(labels)))
    axB.set_xticklabels(labels, fontsize=5.6)
    axB.set_ylabel("ESM-2 variant effect score", fontsize=6.8)
    axB.set_xlim(-0.65, len(data) - 0.35)
    st.recessive_grid(axB, "y")
    _tag(fig, axB, "b")

    # ---- C: odds ratios ----------------------------------------------------
    axC = fig.add_subplot(gs[2])
    axC.set_xscale("log")
    st.null_line(axC, 1.0)
    for i, (name, (o, lo, hi, p), col) in enumerate(rows):
        axC.plot([lo, hi], [i, i], color=col, lw=1.1, solid_capstyle="butt", zorder=2)
        for x in (lo, hi):
            axC.plot([x, x], [i - 0.13, i + 0.13], color=col, lw=1.0, zorder=2)
        axC.plot(o, i, "o", ms=4.6, color=col, mec=st.SURFACE, mew=0.7, zorder=3)
        axC.text(hi * 1.18, i, f"{o:.2f}", va="center", fontsize=6.0, color=st.INK)
    axC.set_yticks(range(len(rows)))
    axC.set_yticklabels([r[0] for r in rows], fontsize=6.2)
    axC.invert_yaxis()
    axC.set_xlabel("Odds ratio (log scale, 95% CI)", fontsize=6.8)
    axC.set_xlim(0.75, 240)
    axC.set_xticks([1, 3, 10, 30, 100])
    axC.set_xticklabels(["1", "3", "10", "30", "100"])
    st.recessive_grid(axC, "x")
    axC.legend(handles=[
        Line2D([], [], marker="o", ls="-", ms=4.6, lw=1.1, color=st.C_OBS,
               label="Codon accessibility"),
        Line2D([], [], marker="o", ls="-", ms=4.6, lw=1.1, color=st.C_ALT,
               label="ESM tolerability")],
        # ncol=1, not 2: see patch07. A two-column key anchored to the left
        # edge grows rightwards off the page whenever the panel narrows.
        loc="lower left", bbox_to_anchor=(0.0, 1.00), ncol=1, handletextpad=0.3,
        borderpad=0.2, labelspacing=0.30, fontsize=6.2)
    _tag(fig, axC, "c")

    _sd = [
        ("a", "accessible, observed as resistance", a, "", "", ""),
        ("a", "accessible, never observed", na - a, "", "", ""),
        ("a", "inaccessible, observed as resistance", c, "", "", ""),
        ("a", "inaccessible, never observed", nc - c, "", "", ""),
        ("a", "percent observed | accessible", 100.0 * a / na, "", "",
         "%d of %d" % (a, na)),
        ("a", "percent observed | inaccessible", 100.0 * c / nc, "", "",
         "%d of %d" % (c, nc)),
    ]
    for _k in present:
        for _v in groups[_k]:
            _sd.append(("b", pretty.get(_k, _k).replace("\n", " "), _v, "", "",
                        "ESM-2 variant effect score"))
    for _lab, _v, _c in rows:
        _p = "" if _v[3] != _v[3] else "P = %.6g" % _v[3]
        _sd.append(("c", _lab.replace("\n", " "), _v[0], _v[1], _v[2], _p))
    st.write_source_data("Fig3", _sd)

    # The unit of this figure is one (site, candidate residue) pair, with no
    # species dimension at all, and the legend has to be able to say so: the
    # universe is every verified site crossed with every non-wild-type residue
    # carrying an ESM-2 score. It read "site-species combinations" until
    # 2026-11-04, which named the wrong unit and made the six inaccessible
    # observed pairs look like they should have been the five species-level
    # two-step events.
    _mvd = first_csv("multivariate_dataset")
    _sites = sorted({r["site"] for r in _mvd})
    _per_site = sorted({sum(1 for r in _mvd if r["site"] == st_) for st_ in _sites})
    _breach = sorted("%s %s>%s" % (r["site"], r["wild_type_aa"], r["candidate_aa"])
                     for r in _mvd
                     if r["accessible"] == "0" and r["observed"] == "1")
    if len(_mvd) != a + (na - a) + c + (nc - c):
        sys.exit("FATAL: the candidate table has %d rows but the two-by-two "
                 "totals %d. One of the two is stale."
                 % (len(_mvd), a + (na - a) + c + (nc - c)))
    st.write_legend_data("figure2_necessary_not_sufficient", {
        "n_candidate_rows": len(_mvd),
        "n_sites": len(_sites),
        "n_candidates_per_site": (_per_site[0] if len(_per_site) == 1
                                  else [_per_site[0], _per_site[-1]]),
        "inaccessible_observed_pairs": _breach,
        "obs_accessible": a,
        "unobs_accessible": na - a,
        "obs_inaccessible": c,
        "unobs_inaccessible": nc - c,
        # na and nc are the totals as the result file states them. a + b
        # and c + d reconstruct the same thing and depend on b and d still
        # meaning what they meant at line 84, which is what broke here.
        "pct_accessible": 100.0 * a / na,
        "pct_inaccessible": 100.0 * c / nc,
        "or_2x2": or2,
        # The figure's own labels, in the figure's own order: a legend that
        # calls a class something the axis does not is a legend the reader
        # has to reconcile.
        "esm_classes": [[pretty.get(k, k).replace("\n", " "), len(groups[k])]
                        for k in order if k in groups],
        "specifications": [
            {"label": lab.replace("\n", " "), "or": v[0], "lo": v[1],
             "hi": v[2], "p": (None if v[3] != v[3] else v[3])}
            for lab, v, _col in rows],
        "headline_label": "Accessibility, site fixed effects adjusted for ESM",
    })

    for ext in ("pdf", "png"):
        p = os.path.join(DIR_FIGS, f"figure2_necessary_not_sufficient.{ext}")
        fig.savefig(p)
        if ext == "pdf":
            st.assert_no_type3(p)
            st.audit(fig, p)
        log(f"  wrote {p}")
    plt.close(fig)
    log("Figure 3 complete.")


if __name__ == "__main__":
    main()
