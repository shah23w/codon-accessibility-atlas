#!/usr/bin/env python3
"""
step8_08_figureS2.py: Supplementary Figure S2: every robustness check as one forest.

The pooled effect is modest (1.22x), so the question a reviewer will ask is not "is it
big" but "does it survive". This figure answers that in one view: the same estimate
recomputed under each threat to validity - pseudoreplication, codon-call error,
detection effort, time, and dependence on the single canonical site.

Every value is read from a result file at run time; nothing is typed in.
"""
from __future__ import annotations

import csv
import math
import os
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import DIR_FIGS, DIR_RESULTS, Tee, ensure_dirs  # noqa: E402
import fig00_style as st  # noqa: E402


# One order, the same in every figure script. Three scripts searched three
# different orders, so with more than one tagged copy of a table in results/
# two figures could be drawn from two different versions of it.
TAGS = ("_primary", "_t1", "_t5", "_final", "")


def read_csv(stem):
    for tag in TAGS:
        p = os.path.join(DIR_RESULTS, f"{stem}{tag}.csv")
        if os.path.exists(p):
            with open(p, newline="", encoding="utf-8-sig") as fh:
                return list(csv.DictReader(fh))
    raise SystemExit(f"missing input: {stem}*.csv")


def main():
    ensure_dirs()
    st.apply()
    log = Tee("step8_08_figureS2")
    log("=== Figure S2 ===")

    loo = {r["held_out"]: r for r in read_csv("leave_one_site_out")}
    ind = {r["collapse_rank"]: r for r in read_csv("phylogenetic_independence")}
    rob = {r["stratum"]: r for r in read_csv("robustness_replicated_subset")}
    eff = {r["stratum"]: r for r in read_csv("control_effort")}
    phy = read_csv("control_phylogenetic")[0]
    tmp = {r["stratum"]: r for r in read_csv("temporal_replication")
           if r["cutoff"] == "2008"}

    def row(label, d, n_key, e_key, p_key):
        # NOTE on the P value. It is carried in the tuple and is NOT drawn:
        # this figure plots the enrichment and n only. That is deliberate and
        # it is also why the permutation columns here are harmless. Several of
        # them are permutation FLOORS rather than estimates -- a run of n draws
        # cannot report below 1/(n+1), so "p = 9.99e-06" means only "nothing
        # smaller appeared in 100,000 draws", where the exact Poisson-binomial
        # value for the headline row is 4.59e-12. The result files now carry
        # p_exact alongside, and if this figure is ever changed to print a P,
        # it must read d["p_exact"], never p_key, or it would reintroduce in a
        # figure the defect contradiction 4 fixed in the analysis.
        # The interval is required, not optional: without it this figure plots
        # point estimates alone, and a stratum that lost its interval would
        # silently go back to doing that.
        if not d.get("ci_lo") or not d.get("ci_hi"):
            sys.exit("FATAL: stratum %r carries no species-cluster interval. "
                     "Rerun the step that writes it (step3_04, step3_05, "
                     "step4_04 or step4_06) before building this figure."
                     % label)
        return ("row", label, int(float(d[n_key])), float(d[e_key]),
                float(d[p_key]), float(d["ci_lo"]), float(d["ci_hi"]))

    items = [
        ("head", "Headline"),
        row("all events, atlas-pool null", loo["none"], "n", "enrichment", "p"),
        ("head", "Phylogenetic non-independence"),
        row("collapsed to one event per genus", ind["genus"],
            "n_events", "enrichment", "p_value"),
        row("collapsed to one event per family", ind["family"],
            "n_events", "enrichment", "p_value"),
        row("collapsed to one event per order", ind["order"],
            "n_events", "enrichment", "p_value"),
        row("permuted within taxonomic order", phy,
            "n_events", "enrichment", "p_permutation"),
        ("head", "Codon-call error"),
        row("call replicated in a 2nd assembly", rob["has independent replicate"],
            "n_events", "enrichment", "p_permutation"),
        row("amino acid confirmed", rob["amino acid replicated"],
            "n_events", "enrichment", "p_permutation"),
        row("exact codon confirmed (strictest)", rob["CODON replicated (strictest)"],
            "n_events", "enrichment", "p_permutation"),
        # The row above counts GCA/GCF twins and version bumps as
        # confirmation; this one does not. step2_01 classifies the pairs and
        # step4_04 applies the classification.
        # Kept to the length of the row above it: the panel is 120 mm and a
        # longer label runs off the left edge. The legend says what
        # "independent" excludes.
        row("exact codon, independent only",
            rob["CODON replicated, independent assemblies only"],
            "n_events", "enrichment", "p_permutation"),
        ("head", "Detection effort"),
        row("well-studied species", eff["well-studied (>= median)"],
            "n_events", "enrichment", "p_permutation"),
        row("less-studied species", eff["less-studied (< median)"],
            "n_events", "enrichment", "p_permutation"),
        ("head", "Time, and the canonical site"),
        row("first reported 2008 or earlier", tmp["<=2008"],
            "n", "enrichment", "p"),
        row("first reported after 2008", tmp[">2008"], "n", "enrichment", "p"),
        row("VGSC L1014 removed entirely", loo["VGSC:L1014"],
            "n", "enrichment", "p"),
    ]

    rows = [i for i in items if i[0] == "row"]
    pooled = rows[0][3]
    worst_p = max(r[4] for r in rows)
    # round the bound UP: plain %.2g turns 0.00104 into "0.001", which claims a
    # tighter bound than the data supports
    _mag = 10 ** (math.floor(math.log10(worst_p)) - 1)
    p_bound = math.ceil(worst_p / _mag) * _mag
    log(f"  {len(rows)} strata, pooled {pooled:.3f}x, worst P = {worst_p:.2g}")
    assert all(r[3] > 1.0 for r in rows), "a stratum fell below 1.0 - check the caption"

    # -------------------------------------------------------------- figure
    # 3.70, not 4.05: the two-line title is gone from the canvas.
    fig = plt.figure(figsize=(st.W15, 3.70))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.415, right=0.965, top=0.975, bottom=0.125)

    ax.axvline(1.0, color=st.INK2, lw=0.7, zorder=1)
    ax.axvline(pooled, color=st.MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)

    ylabels, y = [], 0
    for it in items:
        if it[0] == "head":
            ylabels.append(it[1])
        else:
            _, lab, n, e, p, lo, hi = it
            # An interval with end-caps, not a lollipop stem to the reference
            # line: the stem carried no information and read as uncertainty.
            ax.plot([lo, hi], [y, y], color=st.C_OBS, lw=1.0,
                    solid_capstyle="butt", zorder=2)
            for x in (lo, hi):
                ax.plot([x, x], [y - 0.17, y + 0.17], color=st.C_OBS, lw=1.0,
                        solid_capstyle="butt", zorder=2)
            ax.plot(e, y, "o", ms=4.4, color=st.C_OBS, mec=st.SURFACE, mew=0.7,
                    zorder=3)
            # the pooled reference line runs behind several of these labels; a
            # surface backing keeps both readable
            ax.text(hi + 0.006, y,
                    f"{e:.3f}× [{lo:.2f}, {hi:.2f}]",
                    va="center", ha="left",
                    fontsize=6.0, color=st.INK, zorder=4,
                    bbox=dict(facecolor=st.SURFACE, edgecolor="none",
                              boxstyle="square,pad=0.12"))
            ylabels.append("   %s  n = %d" % (st.cap(lab), n))
        y += 1

    ax.set_yticks(range(len(items)))
    ax.set_yticklabels(ylabels, fontsize=6.4)
    # group headings are flush left and bold; the strata under them are indented
    for tick, it in zip(ax.get_yticklabels(), items):
        if it[0] == "head":
            tick.set_fontweight("bold")
            tick.set_color(st.INK)
        else:
            tick.set_color(st.INK2)
    ax.invert_yaxis()
    ax.set_ylim(len(items) - 0.4, -0.6)
    # Room on the right for the widest "1.308x [1.20, 1.43]". The page check
    # in fig99 is the arbiter: it reported the previous label 21.9 pt off the
    # edge, which is why n moved to the row label and this margin grew.
    ax.set_xlim(min(0.985, min(r[5] for r in rows) - 0.02),
                max(r[6] for r in rows) * 1.14)
    ax.set_xlabel("Observed / expected accessibility  (atlas-pool null)")
    st.recessive_grid(ax, "x")
    # The title belongs to the FIGURE, not to the axes. Set on the axes with
    # loc="left" it started at 41.5% of the page -- where the axes start,
    # because the stratum labels take the left two fifths -- and the second
    # line ran off the right edge. Matplotlib wrote the glyphs anyway and the
    # page simply cut them, so the figure went out reading "every stratum
    # P <= 0.001" with "(permutation)" gone and no sign anything was missing.
    # No title on the canvas. "The effect is modest, and it survives every
    # stratification" is the conclusion the panel is evidence FOR, and the
    # pooled value and the P bound are both stated in the legend.

    st.write_source_data("FigS2", [
        ("", r[1], r[3], "", "", "n = %d; permutation P = %.6g" % (r[2], r[4]))
        for r in rows])

    st.write_legend_data("figureS2_robustness", {
        "n_strata": len(rows),
        "pooled_enrichment": pooled,
        "worst_p": worst_p,
        "p_bound": p_bound,
        "strata": [{"label": r[1], "n": r[2], "enrichment": r[3], "p": r[4],
                    "ci_lo": r[5], "ci_hi": r[6]} for r in rows],
        "n_strata_interval_excludes_1": sum(1 for r in rows if r[5] > 1.0),
    })

    for ext in ("pdf", "png"):
        p = os.path.join(DIR_FIGS, f"figureS2_robustness.{ext}")
        fig.savefig(p)
        if ext == "pdf":
            st.assert_no_type3(p)
            st.audit(fig, p)
        log(f"  wrote {p}")
    plt.close(fig)
    log("Figure S2 complete.")


if __name__ == "__main__":
    main()
