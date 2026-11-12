#!/usr/bin/env python3
"""
Figure 5.  Where in the genetic code codon accessibility can ever matter.

The paper's generality panel. Everything drawn here is derived from the standard
genetic code alone, with no data, and is then shown to predict the data exactly.

  A  the enumeration: for each wild-type residue, how many target residues its
     synonymous codons DISAGREE about reaching. Ten residues can ever give an
     informative test; ten are dead by arithmetic.
  B  the correction. Lysine is two-fold degenerate yet informative, because the
     ATN block it reaches is split between Ile and Met, while glutamate's GTN
     block is not, which the residue-level form of the rule gets wrong.
  C  the prediction meeting the data: every one of the 108 curated events falls
     where the enumeration says it must, 60 uninformative and 48 informative,
     with no exceptions in either direction.
  D  the a priori cross-kingdom prediction, and the one place it overturns a
     claim made about our own system.

Sources, all CURRENT:
    results/enum01_informative_targets.csv
    results/enum01_site_predictions.csv
    results/enum01_predicted_vs_observed.csv   (panel c's contingency table)

No number is typed in. Panel c's counts are read from the enumeration's own
CSV export and asserted against results/canonical_quantities.csv before
anything is drawn, rather than recovered from its prose log: a log
line is written for a reader, and reformatting one silently changes a figure.
"""
import csv, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig00_style as S

DIR_RES = "results"
DIR_FIG = os.path.join("results", "figures")
os.makedirs(DIR_FIG, exist_ok=True)

TARGETS = os.path.join(DIR_RES, "enum01_informative_targets.csv")
SITES = os.path.join(DIR_RES, "enum01_site_predictions.csv")
PVO = os.path.join(DIR_RES, "enum01_predicted_vs_observed.csv")
CANON = os.path.join(DIR_RES, "canonical_quantities.csv")

for p in (TARGETS, SITES, PVO):
    if not os.path.exists(p):
        sys.exit("FATAL: missing %s -- run src/enum01_informativeness.py first" % p)

# --------------------------------------------------------------------------- #
#  inputs, parsed never typed
# --------------------------------------------------------------------------- #
rows = list(csv.DictReader(open(TARGETS)))
enum = [(r["wild_type_aa"], int(r["n_codons"]), int(r["n_informative_targets"]),
         r["informative_targets"]) for r in rows]
enum.sort(key=lambda t: (-t[2], t[0]))

# panel c contingency, read from the enumeration's own export
_pvo = {}
for _r in csv.DictReader(open(PVO)):
    if _r["pair_predicted_informative"] == "":
        continue
    _pvo[(int(_r["pair_predicted_informative"]),
          _r["observed_reachable_fraction"])] = int(_r["n_events"])
_BETW = "strictly between 0 and 1"
_ZERO = "exactly 0 or 1"
try:
    NOT_01, NOT_BETWEEN = _pvo[(0, _ZERO)], _pvo[(0, _BETW)]
    INF_01, INF_BETWEEN = _pvo[(1, _ZERO)], _pvo[(1, _BETW)]
except KeyError as _exc:
    sys.exit("FATAL: %s is missing cell %s -- rerun "
             "src/enum01_informativeness.py" % (os.path.basename(PVO), _exc))

canon = {r["quantity"]: r["value"] for r in csv.DictReader(open(CANON))}
n_inf_canon = int(float(canon["n_events_supply_informative"]))
n_tot_canon = int(float(canon["n_events_total"]))

# Rule 1: the panel asserts its own claim and aborts if it stops being true.
assert INF_BETWEEN == n_inf_canon, (
    "panel C would misreport: enumeration gives %d informative events, "
    "canonical_quantities.csv declares %d" % (INF_BETWEEN, n_inf_canon))
assert NOT_01 + NOT_BETWEEN + INF_01 + INF_BETWEEN == n_tot_canon, (
    "panel C contingency does not sum to the canonical event count")
assert NOT_BETWEEN == 0 and INF_01 == 0, (
    "panel C claims perfect separation but the off-diagonal is not empty")

lit = [r for r in csv.DictReader(open(SITES)) if r["source"] == "literature"]

# --------------------------------------------------------------------------- #
S.apply()
# NEE double column is 183 mm exactly; S.W2 carries it so the width is
# declared in one place and checked against the written media box below.
FIG_W = S.W2
# Margins are generous on the left and bottom for two reasons, both of which
# were invisible until the figure was written at its true size: panel c's row
# labels are two long words, and panel d closes with a three-line note. The
# figure is not saved with bbox_inches="tight", which would grow the canvas
# past the 183 mm column so that nothing appeared to overflow while the panel
# letters fell off the left edge at the real width.
# 0.60 and hspace 0.30: four titles removed, so the gutter that carried them
# goes too. The second row grows relative to the first because panel b filled
# only its top half while the 2x2 and the allele list were cramped.
fig = plt.figure(figsize=(FIG_W, FIG_W * 0.60))
gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1.0], height_ratios=[1.0, 1.15],
                      hspace=0.30, wspace=0.30,
                      left=0.125, right=0.985, top=0.925, bottom=0.135)

# ---- A  the enumeration --------------------------------------------------- #
axA = fig.add_subplot(gs[0, 0])
labels = [e[0] for e in enum]
counts = [e[2] for e in enum]
cols = [S.C_OBS if c else S.MUTED for c in counts]
axA.bar(range(len(labels)), counts, color=cols, width=0.74, linewidth=0)
# Ten residues score zero, and a bar of height zero draws nothing, so the half
# of the result that says "ten cannot" was invisible. Give each of them a
# visible rule on the baseline: the reader should be able to see the absence,
# not infer it from a gap between tick labels.
for _i, _c in enumerate(counts):
    if not _c:
        axA.plot([_i - 0.37, _i + 0.37], [0, 0], color=S.MUTED, lw=1.6,
                 solid_capstyle="butt", zorder=3)
axA.set_xticks(range(len(labels)))
axA.set_xticklabels(labels)
axA.set_ylabel("Informative targets")
axA.set_xlabel("Wild-type residue")
S.recessive_grid(axA, "y")
n_live = sum(1 for c in counts if c)
# No title and no "dead by arithmetic" callout. Ten bars at zero beside ten
# that are not IS the finding, and the legend names it. Nature's figures do
# not argue on the canvas.
S.panel_label(axA, "a")

# ---- B  the lysine correction --------------------------------------------- #
axB = fig.add_subplot(gs[0, 1])
axB.axis("off")
blocks = [
    ("K", "AAA", "ATA", "Ile", True),
    ("K", "AAG", "ATG", "Met", True),
    ("E", "GAA", "GTA", "Val", False),
    ("E", "GAG", "GTG", "Val", False),
]
# Two groups of two, each with its verdict underneath rather than beside it:
# this panel is the narrow column and a right-hand label column collides with
# the residue names.
axB.add_patch(Rectangle((-0.02, 0.595), 1.03, 0.345, transform=axB.transAxes,
                        facecolor=S.C_OBS, alpha=0.07, linewidth=0))
YS = [0.865, 0.725, 0.470, 0.340]
for i, (aa, src, dst, dstaa, split) in enumerate(blocks):
    y = YS[i]
    c = S.C_OBS if split else S.MUTED
    axB.text(0.02, y, aa, fontsize=7, color=S.INK, va="center", weight="bold")
    axB.text(0.13, y, src, fontsize=7, color=S.INK, va="center",
             family="monospace")
    axB.annotate("", xy=(0.47, y), xytext=(0.34, y),
                 arrowprops=dict(arrowstyle="->", color=c, lw=1.1))
    axB.text(0.52, y, dst, fontsize=7, color=c, va="center",
             family="monospace")
    axB.text(0.80, y, dstaa, fontsize=6.8, color=c, va="center")
# These two are KEYS, not commentary: they say what blue and grey mean,
# which is exactly what the specification asks a key to do. The explanation
# of why the ATN block splits and the GTN block does not is in the legend.
axB.text(0.02, 0.625, "Informative", fontsize=6.6, color=S.C_OBS, va="center")
axB.text(0.02, 0.252, "Not informative", fontsize=6.6, color=S.MUTED,
         va="center")
S.panel_label(axB, "b")

# ---- C  prediction against data ------------------------------------------- #
axC = fig.add_subplot(gs[1, 0])
cells = [[NOT_01, NOT_BETWEEN], [INF_01, INF_BETWEEN]]
for i in range(2):
    for j in range(2):
        on_diag = (i == j)
        v = cells[i][j]
        # Row 0 is "predicted not informative", column 1 is "strictly
        # between". That combination cannot occur: if every synonymous codon
        # of a residue agrees about reaching the target, the reachable
        # fraction in any sample is 0 or 1. Hatched, so the panel does not
        # present an arithmetic certainty as an empirical result.
        forced = (i == 0 and j == 1)
        axC.add_patch(Rectangle((j, 1 - i), 1, 1,
                                facecolor=(S.C_OBS if on_diag and v else S.GRID),
                                alpha=0.85 if on_diag and v else 0.5,
                                edgecolor="white", linewidth=2))
        if forced:
            axC.add_patch(Rectangle((j, 1 - i), 1, 1, facecolor="none",
                                    hatch="/////", edgecolor=S.MUTED,
                                    linewidth=0, alpha=0.9))
            axC.text(j + 0.5, 1 - i + 0.26, "forbidden\nby the code",
                     # 5.2, not 4.9: the NEE specification floor is 5 pt and
                     # fig99 rejects anything under it. The cell is about
                     # 120 pt wide, so the label has room at this size.
                     ha="center", va="center", fontsize=5.2, color=S.MUTED,
                     linespacing=1.15)
        axC.text(j + 0.5, 1.5 - i, str(v), ha="center", va="center",
                 fontsize=7,
                 color="white" if (on_diag and v) else S.MUTED,
                 weight="bold" if v else "normal")
axC.set_xlim(0, 2); axC.set_ylim(0, 2)
axC.set_xticks([0.5, 1.5])
axC.set_xticklabels(["fraction\n0 or 1", "strictly\nbetween"])
axC.set_yticks([1.5, 0.5])
axC.set_yticklabels(["predicted\nnot informative", "predicted\ninformative"])
axC.set_xlabel("observed between-species variation in reachability")
for s in axC.spines.values():
    s.set_visible(False)
axC.tick_params(length=0)
# The title asserted the result and the caption underneath restated it; the
# caption was also the run that clipped 0.5 pt off the bottom of the page.
# Two empty cells on the off-diagonal are legible without being narrated.
S.panel_label(axC, "c")

# ---- D  the cross-kingdom prediction -------------------------------------- #
axD = fig.add_subplot(gs[1, 1])
axD.axis("off")
# Pin the panel to the unit square and stop autoscaling. Everything here is
# placed by hand in what is effectively axes fractions, and an artist that
# triggers autoscale makes data coordinates drift away from them - see the
# note at the top of patch06.
axD.set_xlim(0, 1)
axD.set_ylim(0, 1)
axD.set_autoscale_on(False)
seen, order, wt_targets, targets = {}, [], {}, {}
for r in lit:
    k = (r["system"], r["gene"], r["site"], r["wild_type"])
    if k not in seen:
        seen[k] = []
        order.append(k)
    seen[k].append(int(r["predicted_informative"]))
    # The surveyed target residues, so the row names the (wild-type, target)
    # pair the enumeration is actually evaluating rather than the position
    # alone. Informativeness is a property of the pair, which is the whole
    # claim of this figure.
    targets.setdefault(k, []).append(r["resistance_allele"])
    # Empty means the wild-type residue reaches no informative target at all,
    # which is a different kind of dead from a residue that has informative
    # targets none of which is the substitution under surveillance.
    wt_targets[k] = (r.get("wt_informative_targets") or "").split()
y = 0.90
for k in order:
    sysname, gene, site, wt = k
    live = any(seen[k])
    c = S.C_OBS if live else S.MUTED
    axD.text(0.00, y, S.cap(sysname), fontsize=6.6, color=S.MUTED, va="center")
    axD.text(0.17, y, "%s %s%s" % (gene, site, "/".join(targets[k])),
             fontsize=6.8, color=S.INK, va="center")
    # A mark, not a sentence. "can show it" / "cannot" set in words made this
    # panel a table rendered as a figure; a filled against an open circle
    # carries the same two-level distinction and reads at a glance.
    axD.plot(0.88, y, "o", ms=4.6,
             mfc=(S.C_OBS if live else "none"), mec=c, mew=1.0, clip_on=False)
    y -= 0.108
axD.plot(0.055, y + 0.028, "o", ms=4.6,
         mfc=S.C_OBS, mec=S.C_OBS, mew=1.0, clip_on=False)
axD.text(0.10, y + 0.028, "Informative", fontsize=6.2, color=S.INK,
         va="center")
axD.plot(0.465, y + 0.028, "o", ms=4.6,
         mfc="none", mec=S.MUTED, mew=1.0, clip_on=False)
axD.text(0.51, y + 0.028, "dead by arithmetic", fontsize=6.2, color=S.INK,
         va="center")
S.panel_label(axD, "d")

# --------------------------------------------------------------------------- #
S.write_legend_data("figure5_enumeration", {
    # Panel b's two worked examples, so the legend names the codons the panel
    # draws rather than a pair typed into the legend builder. The published
    # legend said "the two glutamate codons both reach valine, GTA and GTG",
    # which reads as though GTA and GTG were the glutamate codons.
    "worked_examples": [{"wild_type": aa, "from_codon": src,
                         "to_codon": dst, "to_residue": dstaa,
                         "informative": bool(split)}
                        for aa, src, dst, dstaa, split in blocks],
    "n_residues_informative": sum(1 for e in enum if e[2]),
    "n_residues_dead": sum(1 for e in enum if not e[2]),
    "targets_by_residue": {e[0]: int(e[2]) for e in enum},
    "contingency": {"predicted_not_informative__fraction_0_or_1": NOT_01,
                    "predicted_not_informative__strictly_between": NOT_BETWEEN,
                    "predicted_informative__fraction_0_or_1": INF_01,
                    "predicted_informative__strictly_between": INF_BETWEEN},
    "n_events_predicted": n_tot_canon,
    # A site is not an allele. These two counts were published as
    # "n_literature_alleles" and "n_literature_informative" and counted SITES,
    # in a paper whose own claim is that informativeness attaches to the
    # (wild-type residue, target residue) PAIR. Both grains are now named for
    # what they are, and the old keys are gone so that nothing can read a
    # site count as an allele count again.
    "n_literature_sites": len(order),
    "n_literature_sites_informative": sum(1 for k in order if any(seen[k])),
    "n_literature_pairs": sum(len(seen[k]) for k in order),
    "n_literature_pairs_informative": sum(sum(seen[k]) for k in order),
    "n_literature_sites_residue_dead": sum(
        1 for k in order if not any(seen[k]) and not wt_targets[k]),
    "n_literature_sites_pair_dead": sum(
        1 for k in order if not any(seen[k]) and wt_targets[k]),
})

for ext in ("pdf", "png"):
    p = os.path.join(DIR_FIG, "figure5_enumeration.%s" % ext)
    # no dpi override: fig00_style sets the house value
    fig.savefig(p)
    print("wrote %s" % p)
S.assert_no_type3(os.path.join(DIR_FIG, "figure5_enumeration.pdf"))
S.audit(fig, os.path.join(DIR_FIG, "figure5_enumeration.pdf"))
print("PASS no Type 3 fonts")

src = os.path.join(DIR_FIG, "source_data", "Fig5_source_data.csv")
os.makedirs(os.path.dirname(src), exist_ok=True)
with open(src, "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(list(S.SOURCE_HEADER))
    for aa, nc, ni, tg in enum:
        w.writerow(["a", aa, ni, "", "", "%d codons; targets: %s" % (nc, tg)])
    for aa, s_, d_, da, sp in blocks:
        w.writerow(["b", "%s %s" % (aa, s_), d_, "", "", "%s; %s" %
                    (da, "codons disagree" if sp else "codons agree")])
    w.writerow(["c", "not_informative_fraction_0_or_1", NOT_01, "", "", ""])
    w.writerow(["c", "not_informative_strictly_between", NOT_BETWEEN, "", "", ""])
    w.writerow(["c", "informative_fraction_0_or_1", INF_01, "", "", ""])
    w.writerow(["c", "informative_strictly_between", INF_BETWEEN, "", "",
                "equals canonical n_events_supply_informative"])
    for k in order:
        w.writerow(["d", " ".join(k), int(any(seen[k])), "", "",
                    "1 = can show the effect"])
print("wrote %s" % src)
print("SUMMARY figure 5 built from enum01 outputs; %d live residues, "
      "%d of %d events informative" % (n_live, INF_BETWEEN, n_tot_canon))
