#!/usr/bin/env python3
"""
step8_05_figure3.py: Figure 3: the mutation spectrum, and how the effect behaves
under progressively stricter nulls.

Panel A  the mechanistic test. Among substitutions that are EQUALLY reachable in one
         step, the ones actually observed are preferentially reachable by a
         transition. A protein language model has no representation of nucleotides,
         so this cannot be tolerability in disguise.
Panel B  that effect across model specifications, so the reader can see it is not an
         artefact of one adjustment set.
Panel C  the honest summary of the whole paper: the effect shrinks as the null gets
         stricter, and survives all three.

Every value is read from a result file or log at run time.
"""
from __future__ import annotations

import csv
import os
import re
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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


def main():
    ensure_dirs()
    st.apply()
    log = Tee("step8_05_figure3")
    log("=== Supplementary Figure 4 ===")

    subs = first_csv("mutation_spectrum_substitutions")
    res_rows = first_csv("mutation_spectrum_results")

    # ---- Panel A -----------------------------------------------------------
    obs = [float(r["ts_fraction"]) for r in subs if r["observed"] == "1"]
    unobs = [float(r["ts_fraction"]) for r in subs if r["observed"] == "0"]
    ts_only_o = sum(1 for r in subs if r["observed"] == "1" and r["transition_only"] == "1")
    ts_only_u = sum(1 for r in subs if r["observed"] == "0" and r["transition_only"] == "1")
    _mw = next((r for r in res_rows
                if r["analysis"] == "ts_fraction_mannwhitney"), None)
    p_mw = float(_mw["p_value"]) if _mw else None
    log(f"  panel A: observed n={len(obs)}, never-observed n={len(unobs)}")

    # ---- Panel B: ts OR across specifications ------------------------------
    # Read from the results table, not from the run log: recovering these three
    # odds ratios from log prose would let a change to a print statement drop a
    # bar silently.
    spec_rows = first_csv("mutation_spectrum_specifications")
    if not spec_rows:
        raise SystemExit("mutation_spectrum_specifications*.csv not found; run step6_04")
    pretty = {"transition fraction alone": "Transition fraction\nalone",
              "BOTH": "+ ESM tolerability",
              "plus codon supply": "+ ESM + codon supply"}
    specs = []
    for key in ("transition fraction alone", "BOTH", "plus codon supply"):
        r = next((x for x in spec_rows
                  if x["specification"] == key and x["term"] == "ts_fraction"), None)
        if r:
            specs.append((pretty[key], float(r["odds_ratio"]), float(r["ci_lo"]),
                          float(r["ci_hi"]), float(r["p_value"])))
    log(f"  panel B: {len(specs)} specifications (from "
        f"mutation_spectrum_specifications*.csv)")

    # ---- Panel C: the three nulls ------------------------------------------
    cs = {r["statistic"]: r for r in
          read_csv(os.path.join(DIR_RESULTS, "corrected_statistics.csv"))}
    def _p_of(table, key):
        """
        The P value for a row, or a hard stop.

        A typed literal standing in for a P the results file does not carry
        would flow into the source-data file as though the pipeline had
        produced it. A number no step computed does not belong in a figure, so
        a missing P stops the build.
        """
        row = table.get(key)
        if row is None:
            raise SystemExit("FATAL: corrected_statistics.csv has no row %r" % key)
        v = (row.get("p_value") or "").strip()
        if not v:
            raise SystemExit(
                "FATAL: corrected_statistics.csv row %r has no P value. Rerun "
                "src/step9_01_corrected_statistics.py, which computes it as an "
                "exact Poisson-binomial tail." % key)
        return float(v)

    sel = read_csv(os.path.join(DIR_RESULTS, "codon_selection_summary.csv"))
    # Matched exactly, not by prefix. A prefix match took the first row whose
    # name began "(a)", so adding a per-gene breakdown silently moved this
    # panel onto one gene while the Results went on quoting the pooled arm.
    WP_NAME = "(a) species' OWN documented allele"
    wp = next((r for r in sel if r["name"] == WP_NAME), None)
    if wp is None:
        sys.exit("FATAL: codon_selection_summary.csv has no row named %r, "
                 "which this panel plots. Rerun "
                 "src/step9_04_codon_selection.py." % WP_NAME)
    def _ci(txt):
        """'[1.159, 1.288]' -> (1.159, 1.288). Empty or malformed -> (None, None)."""
        try:
            a, b = txt.strip().lstrip("[").rstrip("]").split(",")
            return float(a), float(b)
        except (AttributeError, ValueError):
            return None, None

    _ap_lo, _ap_hi = _ci(cs["pooled_enrichment_CI"]["value"])
    _wr_lo, _wr_hi = _ci(cs["within_record_null_enrichment_CI"]["value"]) \
        if "within_record_null_enrichment_CI" in cs else (None, None)
    _wp_lo = float(wp["ci_lo"]) if wp and wp.get("ci_lo") else None
    _wp_hi = float(wp["ci_hi"]) if wp and wp.get("ci_hi") else None
    nulls = [
        ("Atlas pool\nall wild-type species",
         float(cs["pooled_enrichment"]["value"]),
         float(cs["pooled_enrichment"]["p_value"]), _ap_lo, _ap_hi),
        ("Within record\nstudied species only",
         float(cs["within_record_null_enrichment"]["value"]),
         _p_of(cs, "within_record_null_enrichment"), _wr_lo, _wr_hi),
        ("Within protein\nsame gene + residue",
         float(wp["ratio"]) if wp else None,
         float(wp["p_excess"]) if wp else None, _wp_lo, _wp_hi),
    ]
    nulls = [n for n in nulls if n[1]]
    # fig00_style's rule is that uncertainty is drawn, never implied. All three
    # nulls carry a species-cluster bootstrap interval, and a missing interval
    # is a failure rather than a blank.
    _missing = [n[0].split("\n")[0] for n in nulls if n[3] is None or n[4] is None]
    if _missing:
        sys.exit("FATAL: no confidence interval for %s. Rerun "
                 "src/step9_01_corrected_statistics.py and "
                 "src/step9_04_codon_selection.py so that the bootstrap "
                 "intervals are written out." % ", ".join(_missing))
    log(f"  panel c: {len(nulls)} nulls, each with a bootstrap interval")

    # -------------------------------------------------------------- figure
    # Grouped columns, a forest, and a slope. The null ladder is a slope because the
    # substance of the panel is that the estimate falls as the null tightens.
    def _tag(ax, letter):
        """
        Panel letter in FIGURE coordinates, taken from the axes box.

        An axes-fraction offset puts the letter a proportion of the PANEL
        outside the panel, so a narrow panel throws it further. Panel a is
        0.28 of the figure wide, so the old dx of -0.32 landed at figure
        x = -0.034 and the page cut it off.
        """
        pos = ax.get_position()
        st.fig_label(fig, letter,
                     max(pos.x0 - 0.030, 0.005),
                     min(pos.y1 + 0.055, 0.980))

    fig = plt.figure(figsize=(st.W2, 2.78))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.00, 1.02, 0.94], wspace=0.52,
                          left=0.055, right=0.985, top=0.85, bottom=0.233)

    # ---- A: transition fraction, grouped columns --------------------------
    axA = fig.add_subplot(gs[0])
    means = [100 * sum(unobs) / len(unobs), 100 * sum(obs) / len(obs)]
    only = [100 * ts_only_u / len(unobs), 100 * ts_only_o / len(obs)]
    xx = [0, 1]
    w = 0.30
    axA.bar([x - w / 2 - 0.015 for x in xx], means, width=w,
            color=st.C_OBS, linewidth=0)
    axA.bar([x + w / 2 + 0.015 for x in xx], only, width=w,
            color=st.C_OBS, alpha=0.45, linewidth=0)
    axA.set_xticks(xx)
    axA.set_xticklabels([f"Never\nobserved\n(n = {len(unobs)})",
                         f"Observed\n(n = {len(obs)})"], fontsize=6.2)
    axA.set_xlim(-0.62, 1.62)
    axA.set_ylim(0, 58)
    axA.set_ylabel("Proportion (%); denominators differ, see legend",
                   fontsize=6.4)
    st.recessive_grid(axA, "y")
    for i in xx:
        axA.text(i - w / 2 - 0.015, means[i] + 1.2, f"{means[i]:.0f}", ha="center",
                 va="bottom", fontsize=6.0, color=st.INK)
        axA.text(i + w / 2 + 0.015, only[i] + 1.2, f"{only[i]:.0f}", ha="center",
                 va="bottom", fontsize=6.0, color=st.INK2)
    # Shade encodes the measure and position the group, which the axis names,
    # so one swatch per measure keys every bar.
    axA.legend(
        handles=[Patch(color=st.C_OBS, lw=0),
                 Patch(color=st.C_OBS, alpha=0.45, lw=0)],
        labels=["Paths that are transitions (mean)",
                "Reachable only by a transition"],
        handlelength=1.1, handleheight=0.9,
        loc="lower left", bbox_to_anchor=(-0.02, 1.00), ncol=1, handletextpad=0.4,
        borderpad=0.2, labelspacing=0.3, fontsize=6.0)
    if p_mw:
        axA.text(1.58, 52, st.fmt_p(p_mw), ha="right", va="center", fontsize=6.2,
                 color=st.INK)
    _tag(axA, "a")

    # ---- B: transition OR across specifications ---------------------------
    axB = fig.add_subplot(gs[1])
    st.null_line(axB, 1.0)
    for i, (lab, o, lo, hi, p) in enumerate(specs):
        axB.plot([lo, hi], [i, i], color=st.C_OBS, lw=1.1, solid_capstyle="butt",
                 zorder=2)
        for x in (lo, hi):
            axB.plot([x, x], [i - 0.12, i + 0.12], color=st.C_OBS, lw=1.0, zorder=2)
        axB.plot(o, i, "o", ms=4.6, color=st.C_OBS, mec=st.SURFACE, mew=0.7, zorder=3)
        axB.text(hi + 0.14, i, f"{o:.2f}", va="center", fontsize=6.0, color=st.INK)
    axB.set_yticks(range(len(specs)))
    axB.set_yticklabels([s[0] for s in specs], fontsize=6.2)
    axB.invert_yaxis()
    axB.set_xlabel("Odds ratio for transition-reachability\n(95% CI)", fontsize=6.8)
    axB.set_xlim(0.8, max(s[3] for s in specs) * 1.34)
    st.recessive_grid(axB, "x")
    _tag(axB, "b")

    # ---- C: the three nulls, unconnected ----------------------------------
    # Deliberately NOT a slope or a connected line. The three nulls are alternative
    # constructions, not an ordered series, and the estimates are not monotone
    # (1.216, 1.070, 1.243). A connecting line would assert a progression the data
    # do not show.
    axC = fig.add_subplot(gs[2])
    xs = list(range(len(nulls)))
    st.null_line(axC, 1.0, axis="y")
    for x, (lab, ratio, p, lo, hi) in zip(xs, nulls):
        axC.plot([x, x], [lo, hi], color=st.C_OBS, lw=1.1,
                 solid_capstyle="butt", zorder=2)
        for y in (lo, hi):
            axC.plot([x - 0.085, x + 0.085], [y, y], color=st.C_OBS, lw=1.0,
                     zorder=2)
        axC.plot(x, ratio, "o", ms=5.4, color=st.C_OBS, mec=st.SURFACE, mew=1.0,
                 zorder=3)
        axC.text(x, hi + 0.012, f"{ratio:.3f}", ha="center", va="bottom",
                 fontsize=6.2, color=st.INK)
    axC.set_xticks(xs)
    axC.set_xticklabels(["Atlas\npool", "Within\nrecord", "Within\nprotein"],
                        fontsize=6.2)
    axC.set_xlim(-0.58, len(nulls) - 0.42)
    _ylo = min(min(n[3] for n in nulls), 1.0) - 0.020
    axC.set_ylim(_ylo, max(n[4] for n in nulls) + 0.065)
    axC.set_ylabel("Observed / expected", fontsize=6.8)
    st.recessive_grid(axC, "y")
    _tag(axC, "c")

    st.write_source_data("FigS4", (
        [("a", "never observed, mean transition fraction", means[0], "", "",
          "n = %d" % len(unobs)),
         ("a", "observed, mean transition fraction", means[1], "", "",
          "n = %d" % len(obs)),
         ("a", "never observed, reachable only by a transition", only[0], "", "",
          "per cent of substitutions"),
         ("a", "observed, reachable only by a transition", only[1], "", "",
          "per cent of substitutions")]
        + [("b", lab, o, lo, hi, "P = %.6g" % p) for lab, o, lo, hi, p in specs]
        + [("c", lab, ratio, lo, hi, "P = %.6g" % p)
           for lab, ratio, p, lo, hi in nulls]))

    st.write_legend_data("figureS4_mutation_spectrum", {
        "n_observed": len(obs),
        "n_never_observed": len(unobs),
        "transition_fraction_observed": means[1],
        "transition_fraction_never_observed": means[0],
        "transition_only_observed": only[1],
        "transition_only_never_observed": only[0],
        "p_mann_whitney": p_mw,
        "specifications": [{"label": lab, "or": o, "lo": lo, "hi": hi, "p": p}
                           for lab, o, lo, hi, p in specs],
        "nulls": [{"label": lab, "ratio": ratio, "p": p, "lo": lo, "hi": hi}
                  for lab, ratio, p, lo, hi in nulls],
    })

    for ext in ("pdf", "png"):
        p = os.path.join(DIR_FIGS, f"figureS4_mutation_spectrum.{ext}")
        fig.savefig(p)
        if ext == "pdf":
            st.assert_no_type3(p)
            st.audit(fig, p)
        log(f"  wrote {p}")
    plt.close(fig)
    log("Supplementary Figure S4 complete.")


if __name__ == "__main__":
    main()
