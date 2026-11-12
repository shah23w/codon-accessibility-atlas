#!/usr/bin/env python3
"""
step8_10_figureS3_matched_nulls.py: Supplementary Figure S3: the timescale
contrast under matched nulls.

A claim that codon accessibility constrains rapid adaptation but leaves no
signature on deep-time substitution at the same residues would be unsound if
the two classes were scored under different nulls. The rapid class
randomised the CODON holding the residue fixed; the deep-time class randomised the RESIDUE
holding the codon fixed, weighted by how common each residue is in the panel - and that
weighting partly encodes accessibility, because residues that are easy to reach become
common.

This figure reports the corrected analysis. Both classes are scored under all three nulls.
The contrast survives only under the frequency-weighted residue null, i.e. only under the
null with the known circularity; it is absent under the other two, and on the cleanest
subset - the species and sites the two classes actually share - it is absent under every
null. The claim was therefore dropped.

Drawn entirely from results/matched_null_contrast.csv.
"""
from __future__ import annotations

import csv
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import DIR_FIGS, DIR_RESULTS, Tee, ensure_dirs  # noqa: E402
import fig00_style as st  # noqa: E402

NULL_LABEL = {
    "codon_pool": "Codon randomised\n(rapid class's own null)",
    "residue_uniform": "Residue randomised,\nuniform",
    "residue_weighted": "Residue randomised,\nweighted by panel frequency",
}
NULL_ORDER = ["codon_pool", "residue_uniform", "residue_weighted"]
PANELS = [("all", "All observations"),
          ("shared_species_and_sites", "Shared species and sites")]


def main():
    ensure_dirs()
    st.apply()
    log = Tee("step8_10_figureS3_matched_nulls")
    log("=== Figure S3: matched-null timescale contrast ===")

    rows = []
    with open(os.path.join(DIR_RESULTS, "matched_null_contrast.csv"),
              newline="", encoding="utf-8-sig") as fh:
        rows = [r for r in csv.DictReader(fh) if r["cls"] in ("rapid", "deep_time")]
    idx = {(r["subset"], r["null"], r["cls"]): r for r in rows}

    fig = plt.figure(figsize=(st.W2, 2.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.30,
                          left=0.235, right=0.90, top=0.74, bottom=0.20)

    for pi, (subset, ptitle) in enumerate(PANELS):
        ax = fig.add_subplot(gs[pi])
        ax.axvline(1.0, color=st.INK2, lw=0.7, zorder=1)
        yt, ylab = [], []
        xmax = 1.0
        for i, nl in enumerate(NULL_ORDER):
            base = i * 1.0
            yt.append(base)
            ylab.append(NULL_LABEL[nl])
            for j, (cls, col, dy) in enumerate((("rapid", st.C_OBS, -0.19),
                                                ("deep_time", st.MUTED, 0.19))):
                r = idx.get((subset, nl, cls))
                if not r:
                    continue
                lo, hi = float(r["ci_lo"]), float(r["ci_hi"])
                ratio = float(r["ratio"])
                xmax = max(xmax, hi)
                ax.plot([lo, hi], [base + dy] * 2, color=col, lw=1.2,
                        solid_capstyle="butt", zorder=2)
                # End-caps, so an interval cannot be mistaken for the stem
                # style used elsewhere, and so its ends are unambiguous where
                # two intervals nearly touch.
                for _x in (lo, hi):
                    ax.plot([_x, _x], [base + dy - 0.085, base + dy + 0.085],
                            color=col, lw=1.2, solid_capstyle="butt", zorder=2)
                ax.plot(ratio, base + dy, "o", ms=4.6, color=col, mec=st.SURFACE,
                        mew=0.7, zorder=3)
                # The bounds are printed as well as drawn: the legend and the
                # main text both turn on whether these intervals separate, and
                # a reader should not have to measure the figure to check.
                ax.text(hi + 0.016, base + dy,
                        f"{ratio:.3f} [{lo:.3f}, {hi:.3f}]", va="center",
                        fontsize=5.8, color=st.INK, zorder=4,
                        bbox=dict(facecolor=st.SURFACE, edgecolor="none",
                                  boxstyle="square,pad=0.10"))
            # Whether the two intervals separate under each null is a reading
            # of the panel, and the legend states which null separates them.
            # It was drawn here as the italic words "separated"/"overlapping".
        ax.set_yticks(yt)
        ax.set_yticklabels(ylab if pi == 0 else [""] * len(ylab), fontsize=6.0)
        ax.set_title(ptitle, pad=5, loc="left", fontsize=6.8, color=st.INK)
        ax.invert_yaxis()
        ax.set_ylim(len(NULL_ORDER) - 0.35, -0.62)
        # Wider than before: the labels now carry the interval bounds too.
        ax.set_xlim(0.93, xmax * 1.30)
        ax.set_xlabel("Observed / expected accessibility")
        st.recessive_grid(ax, "x")
        st.panel_label(ax, "AB"[pi], dx=-0.62 if pi == 0 else -0.10, dy=1.16)

    fig.legend(handles=[
        Line2D([], [], marker="o", ls="-", ms=4.6, lw=1.2, color=st.C_OBS,
               label="Rapid adaptation (resistance events)"),
        Line2D([], [], marker="o", ls="-", ms=4.6, lw=1.2, color=st.MUTED,
               label="Deep time (inferred substitutions)")],
        loc="upper left", bbox_to_anchor=(0.235, 1.0), ncol=2, frameon=False,
        handletextpad=0.4, columnspacing=1.6, fontsize=6.4)

    _sep = {}
    for _subset, _ in PANELS:
        for _nl in NULL_ORDER:
            _r, _d = idx.get((_subset, _nl, "rapid")), idx.get((_subset, _nl, "deep_time"))
            if _r and _d:
                # "ratio", not "enrichment": this file names the column
                # differently from the other robustness outputs, and the
                # figure itself reads float(r["ratio"]).
                _sep["%s|%s" % (_subset, _nl)] = {
                    "rapid": float(_r["ratio"]),
                    "rapid_lo": float(_r["ci_lo"]),
                    "rapid_hi": float(_r["ci_hi"]),
                    "rapid_n": int(float(_r["n"])),
                    "rapid_p_exact": float(_r["p_exact"]),
                    "deep_time": float(_d["ratio"]),
                    "deep_lo": float(_d["ci_lo"]),
                    "deep_hi": float(_d["ci_hi"]),
                    "deep_n": int(float(_d["n"])),
                    "deep_p_exact": float(_d["p_exact"]),
                    # Two-sided. The earlier test asked only whether rapid
                    # adaptation sat above deep time, so the codon-randomised
                    # contrast - where deep time is the higher of the two and
                    # the intervals do not overlap (1.388 > 1.366) - was
                    # counted as overlapping. The question is whether the two
                    # classes differ, not whether one exceeds the other.
                    "separated": bool(float(_r["ci_lo"]) > float(_d["ci_hi"])
                                      or float(_d["ci_lo"]) > float(_r["ci_hi"])),
                    "direction": ("rapid higher"
                                  if float(_r["ci_lo"]) > float(_d["ci_hi"])
                                  else "deep time higher"
                                  if float(_d["ci_lo"]) > float(_r["ci_hi"])
                                  else "overlapping"),
                }
    _sd = []
    for _pi, (_subset, _) in enumerate(PANELS):
        _pan = "ab"[_pi]
        for _nl in NULL_ORDER:
            for _cls in ("rapid", "deep_time"):
                _r = idx.get((_subset, _nl, _cls))
                if not _r:
                    continue
                _sd.append((_pan, "%s | %s" % (_nl, _cls), float(_r["ratio"]),
                            float(_r["ci_lo"]), float(_r["ci_hi"]),
                            "n = %s; exact P = %s"
                            % (_r["n"], _r["p_exact"])))
    st.write_source_data("FigS3", _sd)

    st.write_legend_data("figureS3_matched_null_contrast", {
        "nulls": list(NULL_ORDER),
        "subsets": [s for s, _ in PANELS],
        "contrasts": _sep,
        "n_separating": sum(1 for v in _sep.values() if v["separated"]),
    })

    for ext in ("pdf", "png"):
        p = os.path.join(DIR_FIGS, f"figureS3_matched_null_contrast.{ext}")
        fig.savefig(p)
        if ext == "pdf":
            st.assert_no_type3(p)
            st.audit(fig, p)
        log(f"  wrote {p}")
    plt.close(fig)
    log("Figure S3 complete.")


if __name__ == "__main__":
    main()
