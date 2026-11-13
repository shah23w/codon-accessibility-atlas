#!/usr/bin/env python3
"""
Figure 4.  Mutational supply orders the alleles the gate admits.

The paper's mechanism figure, and the one that has to be honest about its own
limits. Its central estimate comes from a single locus, and the figure says
so.

  A  why the estimate is confined to one gene. The sodium channel carries 34 of
     the 47 informative strata; the five other genes supply 13 between them.
  B  beta at the sodium channel against the three published values, and against
     the two theoretical predictions that our interval cannot separate.
  C  the benchmark. Our estimator applied to the published data of Cano et al.
     returns their published coefficients, so the comparison in B is between
     like and like.
  D  what more data would buy: the multiplier at which beta = 1 becomes
     separable from beta = 0.6.

Sources, all CURRENT:
    results/step9_02_gene_information.csv
    results/step9_02_beta_context.csv
    results/step9_02_power.csv
    results/step9_03_benchmark_cano.csv

No number is typed in. The three literature values in panel B are read from
the benchmark file, which carries them as published, not retyped here.
"""
import csv, os, re, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig00_style as S

DIR_RES = "results"
DIR_FIG = os.path.join("results", "figures")
os.makedirs(DIR_FIG, exist_ok=True)

GENEINFO = os.path.join(DIR_RES, "step9_02_gene_information.csv")
CTX = os.path.join(DIR_RES, "step9_02_beta_context.csv")
POWER = os.path.join(DIR_RES, "step9_02_power.csv")
BENCH = os.path.join(DIR_RES, "step9_03_benchmark_cano.csv")
for p in (GENEINFO, CTX, POWER, BENCH):
    if not os.path.exists(p):
        sys.exit("FATAL: missing %s -- run step9_02 and step9_03 first" % p)

genes = list(csv.DictReader(open(GENEINFO)))
# ortho01 established that GluCl_Px and GluCl_Tu are two sites on ONE gene.
# Drawing them as two bars would contradict the paper's own gene count.
_merged, _gl = [], None
for g in genes:
    if g["gene"].startswith("GluCl"):
        if _gl is None:
            _gl = dict(g); _gl["gene"] = "GluCl"
        else:
            for k in ("strata", "candidates", "observed"):
                _gl[k] = str(int(_gl[k]) + int(g[k]))
    else:
        _merged.append(g)
if _gl:
    _merged.append(_gl)
genes = _merged
ctx = {r["quantity"]: r for r in csv.DictReader(open(CTX))}
power = list(csv.DictReader(open(POWER)))
bench = list(csv.DictReader(open(BENCH)))


def interval(note):
    m = re.search(r"\[\s*([-0-9.]+),\s*([-0-9.]+)\s*\]", note)
    if not m:
        sys.exit("FATAL: could not parse an interval from %r" % note)
    return float(m.group(1)), float(m.group(2))


# The pooled beta and its interval. The point estimate is in the context
# table; the profile interval is only in the register, so it is read from
# there rather than restated here.
_CANON = os.path.join(DIR_RES, "canonical_quantities.csv")
if not os.path.exists(_CANON):
    raise SystemExit("FATAL: %s is missing; the pooled beta interval cannot "
                     "be read." % _CANON)
_canon = {r["quantity"]: r["value"]
          for r in csv.DictReader(open(_CANON, encoding="utf-8-sig"))}
for _k in ("beta_mutational_supply", "beta_CI95"):
    if _k not in _canon:
        raise SystemExit("FATAL: canonical_quantities.csv has no %r." % _k)
B_POOLED = float(_canon["beta_mutational_supply"])
B_POOLED_LO, B_POOLED_HI = interval(_canon["beta_CI95"])
B_VGSC = float(ctx["beta_vgsc_only"]["value"])
B_VGSC_LO, B_VGSC_HI = interval(ctx["beta_vgsc_only"]["note"])
B_REST = float(ctx["beta_without_vgsc"]["value"])
B_REST_LO, B_REST_HI = interval(ctx["beta_without_vgsc"]["note"])
P_BOOT = float(ctx["beta_gene_interaction_P_bootstrap"]["value"])
SHARE = float(ctx["vgsc_share_of_informative_strata"]["value"])
MULT = ctx["multiplier_for_80pc_power"]["value"]

# the figure asserts its own claim and aborts if it stops being true
n_strata = {g["gene"]: int(g["strata"]) for g in genes}
tot = sum(n_strata.values())
assert abs(n_strata.get("VGSC", 0) / tot - SHARE) < 1e-3, (
    "panel A would misreport the VGSC share: %d/%d against declared %.4f"
    % (n_strata.get("VGSC", 0), tot, SHARE))
assert B_VGSC_LO > B_REST_LO and B_VGSC > B_REST, (
    "panel B assumes the VGSC estimate is the higher one")

# --------------------------------------------------------------------------- #
S.apply()
# NEE double column is 183 mm exactly; see src/fig00_style.py.
FIG_W = S.W2
# 0.50, not 0.56: the four panel titles are gone and the gutter they needed
# with them.
fig = plt.figure(figsize=(FIG_W, FIG_W * 0.50))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.25], height_ratios=[1.0, 1.0],
                      hspace=0.34, wspace=0.34,
                      left=0.085, right=0.975, top=0.915, bottom=0.105)

# ---- A  where the information is ------------------------------------------ #
axA = fig.add_subplot(gs[0, 0])
order = sorted(genes, key=lambda g: -int(g["strata"]))
names = [g["gene"] for g in order]
vals = [int(g["strata"]) for g in order]
cols = [S.C_OBS if n == "VGSC" else S.MUTED for n in names]
axA.barh(range(len(names))[::-1], vals, color=cols, height=0.68, linewidth=0)
axA.set_yticks(range(len(names))[::-1])
axA.set_yticklabels(names)
axA.set_xlabel("Informative strata")
S.recessive_grid(axA, "x")
axA.text(vals[0] - 1.2, len(names) - 1, "%d of %d" % (vals[0], tot),
         ha="right", va="center", color="white", fontsize=7.0, weight="bold")
S.panel_label(axA, "a")

# ---- B  beta in context ---------------------------------------------------- #
axB = fig.add_subplot(gs[0, 1])
rows = []
for b in bench:
    rows.append((b["species"], float(b["beta_published"]),
                 float(b["ci_lo_published"]), float(b["ci_hi_published"]),
                 S.MUTED, "clonal"))
# The pooled estimate across all informative strata. The legend's claim that
# "our interval contains both and separates neither" is about THIS interval,
# and it was the one row the panel did not draw.
# One line per row: the colour separates this study (blue, orange) from the
# published estimates (grey), and the legend says so.
rows.append(("All genes pooled", B_POOLED, B_POOLED_LO, B_POOLED_HI,
             S.C_OBS, "sexual"))
rows.append(("Sodium channel", B_VGSC, B_VGSC_LO, B_VGSC_HI,
             S.C_OBS, "sexual"))
rows.append(("Five other genes", B_REST, B_REST_LO, B_REST_HI,
             S.C_ALT, "sexual"))

axB.axvspan(0.1, 0.6, color=S.MUTED, alpha=0.10, linewidth=0)
axB.axvline(1.0, color=S.INK, lw=0.8, ls=(0, (4, 3)))
# the band labels need their own strip of space, or they land on the lowest
# interval; the axis is extended downward rather than the labels squeezed in
# Deeper than the annotation needs, so the two band labels clear the lowest
# interval rather than sitting a hair beneath it.
axB.set_ylim(-1.45, len(rows) - 0.45)
axB.text(0.35, -1.25, "Asexual\nsupply curve", ha="center", va="bottom",
         fontsize=6.6, color=S.MUTED)
axB.text(1.04, -1.25, "Sexual\norigin-fixation", ha="left", va="bottom",
         fontsize=6.6, color=S.INK)

for i, (lab, b, lo, hi, c, _) in enumerate(rows):
    y = len(rows) - 1 - i
    axB.plot([lo, hi], [y, y], color=c, lw=1.6, solid_capstyle="round")
    axB.plot([b], [y], "o", color=c, ms=5.0, mec="white", mew=0.8)
axB.set_yticks(range(len(rows))[::-1])
# The three benchmark rows are binomials and are set in italic; the rows this
# study contributes are not species names and stay roman.
_btl = axB.set_yticklabels([r[0] for r in rows], fontsize=6.8)
for _lab in _btl:
    if "." in _lab.get_text().split("\n")[0] and "(" not in _lab.get_text():
        _lab.set_fontstyle("italic")
axB.set_xlabel("Mutational influence $\\beta$")
S.null_line(axB, 0.0)
axB.set_xlim(-1.1, 2.0)
S.recessive_grid(axB, "x")
S.panel_label(axB, "b")

# ---- C  the benchmark ------------------------------------------------------ #
axC = fig.add_subplot(gs[1, 0])
for b in bench:
    x, y = float(b["beta_published"]), float(b["beta_ours"])
    axC.plot([x], [y], "o", color=S.C_OBS, ms=6, mec="white", mew=0.8)
    # One fixed offset put every label on the same side of the identity line.
    # E. coli is the only benchmark above the line, so its label was written
    # across the line and onto a neighbouring marker. Offset each label away
    # from the line, on the side its own point already falls, and set the
    # binomials in italic.
    above = y > x
    axC.annotate(b["species"], xy=(x, y),
                 xytext=(-6, 7) if above else (5, -11),
                 textcoords="offset points", fontsize=6.8, color=S.INK,
                 fontstyle="italic",
                 ha="right" if above else "left",
                 va="bottom" if above else "top")
lim = (0.6, 1.25)
axC.plot(lim, lim, color=S.MUTED, lw=0.8, ls=(0, (4, 3)))
axC.set_xlim(*lim); axC.set_ylim(*lim)
axC.set_xlabel("$\\beta$ published by Cano et al.")
axC.set_ylabel("$\\beta$ from our estimator")
axC.set_aspect("equal", adjustable="box")
S.recessive_grid(axC, "both")
# worst is kept: the source-data file and the legend both report it.
# "largest difference 0.089" is a derived statistic, not a label on the data,
# so the legend states it, where a reader can see what it came from.
worst = max(float(b["abs_difference"]) for b in bench)
S.panel_label(axC, "c")

# ---- D  power -------------------------------------------------------------- #
axD = fig.add_subplot(gs[1, 1])
mx = [int(r["multiplier"]) for r in power]
my = [float(r["power"]) for r in power]
ms_ = [int(r["strata"]) for r in power]
axD.plot(mx, my, "-o", color=S.C_OBS, lw=1.5, ms=4.5, mec="white", mew=0.8)
axD.axhline(0.80, color=S.INK, lw=0.8, ls=(0, (4, 3)))
axD.text(mx[-1], 0.815, "80%", ha="right", va="bottom", fontsize=6.8,
         color=S.INK)
try:
    k = int(MULT)
    j = mx.index(k)
    axD.plot([k], [my[j]], "o", color=S.C_ALT, ms=8, mec="white", mew=1.0,
             zorder=5)
    # the marked point is identified in the legend, not on the canvas
    axD.annotate("%dx" % k, xy=(k, my[j]), xytext=(0, 9),
                 textcoords="offset points", ha="center",
                 fontsize=6.5, color=S.C_ALT)
except (ValueError, TypeError):
    pass
axD.set_xlabel("Multiple of the present informative strata")
axD.set_ylabel("Power to separate\n$\\beta=1$ from $\\beta=0.6$")
axD.set_ylim(0, 1.05)
axD.set_xticks(mx)
S.recessive_grid(axD, "y")
S.panel_label(axD, "d")

# --------------------------------------------------------------------------- #
S.write_legend_data("figure4_mutational_supply", {
    "strata_by_gene": {g["gene"]: int(g["strata"]) for g in order},
    "strata_total": tot,
    "largest_benchmark_difference": worst,
    "beta_pooled": B_POOLED, "beta_pooled_lo": B_POOLED_LO,
    "beta_pooled_hi": B_POOLED_HI,
    "beta_vgsc": B_VGSC, "beta_vgsc_lo": B_VGSC_LO, "beta_vgsc_hi": B_VGSC_HI,
    "beta_rest": B_REST, "beta_rest_lo": B_REST_LO, "beta_rest_hi": B_REST_HI,
    "benchmarks": [{"species": b["species"],
                    "beta_published": float(b["beta_published"]),
                    "beta_ours": float(b["beta_ours"])} for b in bench],
    "power_multiples": [int(v) for v in mx],
    "power_values": [float(v) for v in my],
    "power_strata": [int(v) for v in ms_],
    "power_highlight_multiple": int(MULT),
})

for ext in ("pdf", "png"):
    p = os.path.join(DIR_FIG, "figure4_mutational_supply.%s" % ext)
    # no dpi override: fig00_style sets the house value, and 400 here
    # silently produced a different raster resolution from every other figure
    fig.savefig(p)
    print("wrote %s" % p)
S.assert_no_type3(os.path.join(DIR_FIG, "figure4_mutational_supply.pdf"))
S.audit(fig, os.path.join(DIR_FIG, "figure4_mutational_supply.pdf"))
print("PASS no Type 3 fonts")

src = os.path.join(DIR_FIG, "source_data", "Fig4_source_data.csv")
os.makedirs(os.path.dirname(src), exist_ok=True)
with open(src, "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["panel", "key", "value", "low", "high", "note"])
    for g in order:
        w.writerow(["a", g["gene"], g["strata"], "", "",
                    "%s candidates, %s observed" % (g["candidates"],
                                                    g["observed"])])
    for lab, b, lo, hi, colour, repro in rows:
        # the panel tells this study from the published estimates by colour;
        # the table has to say it in words
        if colour != S.MUTED:
            lab += " (this study)"
        w.writerow(["b", lab.replace("\n", " "), b, lo, hi, repro])
    for b in bench:
        w.writerow(["c", b["species"], b["beta_ours"], "", "",
                    "published %s, difference %s"
                    % (b["beta_published"], b["abs_difference"])])
    for r in power:
        w.writerow(["d", "multiplier %s" % r["multiplier"], r["power"], "", "",
                    "%s strata, %s simulations" % (r["strata"], r["nsim"])])
print("wrote %s" % src)
print("SUMMARY figure 4: beta_VGSC %.2f [%.2f, %.2f] on %d of %d strata; "
      "interaction bootstrap P %.3f; %sx for 80%% power"
      % (B_VGSC, B_VGSC_LO, B_VGSC_HI, n_strata["VGSC"], tot, P_BOOT, MULT))
