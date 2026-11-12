#!/usr/bin/env python3
"""
fig01_dataset_and_context.py
============================
Figure 1. The dataset, its phylogenetic context, and the mechanism.

  a  fan phylogeny of 463 arthropods on a dated tree, branches coloured by order,
     with two annotation rings: the species' codon class at VGSC L1014, and
     whether that species appears in the resistance record
  b  panel composition: species per order
  c  the mechanism: which leucine codons can reach which resistance alleles at
     VGSC L1014 in one nucleotide substitution
  d  curated events per gene, split by whether the observed allele was accessible

NOTHING HERE IS DRAWN BY HAND. Panel c is computed from the standard genetic code,
not from a diagram someone sketched: the script enumerates every single-nucleotide
neighbour of every leucine codon and reports which resistance alleles they encode.
If the genetic code says a codon cannot reach an allele, the panel says so.

DESIGN

The panel presents evidence; the caption interprets it. No sentences inside the
axes, no arrows pointing at conclusions, no callouts. Labels are noun phrases.

Colour follows the four-job rule. Order identity in panel a is CATEGORICAL, eight
slots in fixed order, and each is directly labelled by an outer arc -- identity is
never carried by colour alone. Codon class in the inner ring is a three-level
ORDINAL distinction (reaches the allele / does not / already carries it) and uses
three well-separated hues rather than a ramp, because "already Phe" is not a
midpoint between the other two. Panel d uses one hue for accessible and a neutral
grey for breaches: grey recedes, which is correct, because breaches are the
exception being counted, not a competing category.

WHY VGSC L1014

It is the canonical knockdown-resistance site, it has the widest coverage in the
panel (432 of 463 tree species carry a codon), and it has the most documented
species (15). It is also the clearest case of the paper's claim: two of the six
leucine codons cannot reach phenylalanine in one step at all.

OUTPUTS
    results/figures/fig01_dataset_and_context.pdf
    results/figures/fig01_dataset_and_context.png   (600 dpi)

"""
import collections, csv, math, os, re, sys

TREE  = "data/refs/timetree_arthropods_collapsed.nwk"
ATLAS = "data/processed/codon_atlas.csv"
CORE  = "results/core_test_events_primary.csv"
SEL   = "results/codon_selection_per_site.csv"
OUTDIR = os.path.join("results", "figures")
FOCAL_GENE, FOCAL_SITE = "VGSC", "L1014"

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Wedge, Rectangle
    from matplotlib.lines import Line2D
except ImportError:
    sys.exit("needs numpy + matplotlib:  micromamba install -y numpy matplotlib")

for p in (TREE, ATLAS, CORE):
    if not os.path.exists(p):
        sys.exit("FATAL: missing %s" % p)

# --------------------------------------------------------------------------- #
#  palette
# --------------------------------------------------------------------------- #
# The palette comes from fig00_style, so this figure sits on the same white
# ground and uses the same greys as every other. Two values that differ by one
# hex
# digit (#8a8983 against #8a8985). The constants below are now bound to the one
# style module, so the names used throughout this script still work.
import fig00_style as S

S.apply(verbose=False)

SURFACE, INK, INK2, INK3 = S.SURFACE, S.INK, S.INK2, S.MUTED
# Eight hues are needed here and nowhere else: the fan phylogeny colours seven
# named orders plus a neutral "Other". The first three are the validated slots.
ORDER_HUES = [S.BLUE, S.ORANGE, S.AQUA, "#eda100",
              "#e87ba4", "#008300", "#4a3aa7", "#b9b8b2"]   # 8th = Other, neutral
REACH, NOREACH, ALREADY = S.BLUE, "#a9a8a1", S.ORANGE
# Species the atlas has no VGSC L1014 codon for. Drawn in their own shade
# rather than left blank, so a missing genome cannot be read as a negative
# result.
NOCALL = "#eceae4"
ACCESSIBLE, BREACH = S.BLUE, "#b9b8b2"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})

# --------------------------------------------------------------------------- #
#  Newick (same parser as phylo02/phylo03/perm01 -- change all together)
# --------------------------------------------------------------------------- #
class Newick:
    def __init__(self, text):
        self.s = text.replace("\n", "").replace("\r", "").replace("\t", "")
        i = self.s.rfind(";")
        if i >= 0: self.s = self.s[:i]
        self.i = 0
    def parse(self): return self._node()
    def _node(self):
        if self.i < len(self.s) and self.s[self.i] == "(":
            self.i += 1; kids = []
            while True:
                kids.append(self._node())
                if self.s[self.i] == ",": self.i += 1; continue
                if self.s[self.i] == ")": self.i += 1; break
                raise ValueError("unexpected %r" % self.s[self.i])
            return dict(kids=kids, label=self._label(), length=self._length())
        return dict(kids=[], label=self._label(), length=self._length())
    def _label(self):
        if self.i < len(self.s) and self.s[self.i] in "'\"":
            q = self.s[self.i]; self.i += 1; st = self.i
            while self.i < len(self.s) and self.s[self.i] != q: self.i += 1
            lab = self.s[st:self.i]; self.i += 1; return lab.strip()
        st = self.i
        while self.i < len(self.s) and self.s[self.i] not in "(),:;": self.i += 1
        return self.s[st:self.i].strip()
    def _length(self):
        if self.i < len(self.s) and self.s[self.i] == ":":
            self.i += 1; st = self.i
            while self.i < len(self.s) and (self.s[self.i].isdigit()
                                           or self.s[self.i] in ".eE+-"): self.i += 1
            try: return float(self.s[st:self.i])
            except ValueError: return None
        return None

def norm(s): return " ".join(s.replace("_", " ").split()).strip()

root = Newick(open(TREE).read()).parse()

# ladderise so clades read cleanly, then collect tips left to right
def count_tips(n):
    n["_nt"] = 1 if not n["kids"] else sum(count_tips(k) for k in n["kids"])
    return n["_nt"]
count_tips(root)
def ladderise(n):
    if n["kids"]:
        n["kids"].sort(key=lambda k: k["_nt"])
        for k in n["kids"]: ladderise(k)
ladderise(root)

tips = []
def collect(n, depth=0.0):
    d = depth + (n["length"] or 0.0)
    n["_r"] = d
    if n["kids"]:
        for k in n["kids"]: collect(k, d)
    else:
        tips.append(n)
collect(root)
NT = len(tips)
ROOT_AGE = max(t["_r"] for t in tips)

# angles: leave a wedge open so the start and end of the circle are distinguishable
GAP_DEG = 8.0
span = 360.0 - GAP_DEG
for i, t in enumerate(tips):
    t["_th"] = math.radians(90.0 - GAP_DEG / 2.0 - span * (i + 0.5) / NT)
def set_angles(n):
    if not n["kids"]: return n["_th"]
    a = [set_angles(k) for k in n["kids"]]
    n["_th"] = 0.5 * (min(a) + max(a)); return n["_th"]
set_angles(root)

# --------------------------------------------------------------------------- #
#  data
# --------------------------------------------------------------------------- #
istrue = lambda v: str(v).strip().lower() in ("true", "1", "yes")
atlas = list(csv.DictReader(open(ATLAS)))
tipnames = {norm(t["label"]) for t in tips}

order_of, focal_codon = {}, {}
for r in atlas:
    sp = norm(r["species"])
    if sp not in tipnames: continue
    order_of.setdefault(sp, (r.get("order") or "Other").strip())
    if r["gene"] == FOCAL_GENE and r["site_label"] == FOCAL_SITE and r.get("codon"):
        focal_codon[sp] = r["codon"].strip().upper()

documented = set()
if os.path.exists(SEL):
    for r in csv.DictReader(open(SEL)):
        if istrue(r.get("documented_here")):
            sp = norm(r["species"])
            if sp in tipnames: documented.add(sp)

ocount = collections.Counter(order_of.values())
TOP = [o for o, _ in ocount.most_common(7)]
ORDER_COLOUR = {o: ORDER_HUES[i] for i, o in enumerate(TOP)}
def ocol(sp): return ORDER_COLOUR.get(order_of.get(sp, ""), ORDER_HUES[7])

# --------------------------------------------------------------------------- #
#  genetic code -- panel c is computed, not drawn
# --------------------------------------------------------------------------- #
BASES = "TCAG"
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODE = {}
i = 0
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            CODE[b1 + b2 + b3] = AAS[i]; i += 1

def one_step_aas(codon):
    out = set()
    for pos in range(3):
        for b in "ACGT":
            if b == codon[pos]: continue
            m = codon[:pos] + b + codon[pos + 1:]
            out.add(CODE[m])
    return out

LEU = sorted([c for c in CODE if CODE[c] == "L"])
focal_alleles = sorted({r["resistance_allele"] for r in atlas
                        if r["gene"] == FOCAL_GENE and r["site_label"] == FOCAL_SITE})
REACH_MAP = {c: one_step_aas(c) for c in LEU}

TARGET = "F"          # the canonical kdr allele, for the ring classification
def codon_class(sp):
    c = focal_codon.get(sp)
    if not c or c not in CODE: return None
    if CODE[c] == TARGET: return "already"
    if CODE[c] != "L": return None
    return "reach" if TARGET in REACH_MAP.get(c, set()) else "noreach"

CLASS_COLOUR = {"reach": REACH, "noreach": NOREACH, "already": ALREADY}

# --------------------------------------------------------------------------- #
#  figure
# --------------------------------------------------------------------------- #
# NEE double column is 183 mm exactly; 158 mm tall is inside the 170 mm
# ceiling that leaves room for the legend.
fig = plt.figure(figsize=(S.W2, 158 / 25.4))

def panel_tag(ax, s, x=-0.02, y=1.04):
    ax.text(x, y, str(s).lower(), transform=ax.transAxes,
            fontsize=S.PANEL_LETTER_PT, fontweight="bold",
            color=INK, va="bottom", ha="left")

# ---- a: fan tree --------------------------------------------------------- #
# Three things were wrong with the first version and are fixed here.
#
#  1. The open wedge was widened to 26 degrees to make room for age labels, which
#     broke the ring open so it no longer read as a circle. It is back to 8, and
#     the time scale has moved out of the tree entirely (see the scale bar below).
#     Labelling the guide circles could never have worked: near the root every
#     branch converges on the centre, so ANY label there is crossed by branches.
#  2. Only the tips carried colour, so the tree read as grey scaffolding with a
#     coloured fringe. Every branch is now coloured by its clade -- a node whose
#     descendants are all one order takes that order's hue -- so the groups are
#     visible as shapes, which is what makes a fan tree worth drawing.
#  3. Labels on the left half of the circle were printed upside down.
axa = fig.add_axes([-0.015, 0.015, 0.70, 0.965])
axa.set_aspect("equal"); axa.axis("off"); axa.patch.set_visible(False)
RT = 1.0
def rad(r): return RT * (r / ROOT_AGE)

# which order does each clade belong to? (None where a clade spans several)
def clade_orders(n):
    if not n["kids"]:
        n["_o"] = {order_of.get(norm(n["label"]), "")}
    else:
        n["_o"] = set()
        for k in n["kids"]:
            clade_orders(k); n["_o"] |= k["_o"]
    return n["_o"]
clade_orders(root)
BRANCH_GREY = "#d6d5cf"
def bcol(n):
    if len(n["_o"]) == 1:
        o = next(iter(n["_o"]))
        if o in ORDER_COLOUR: return ORDER_COLOUR[o]
    return BRANCH_GREY

# faint depth rings, no labels on them
for age in range(100, int(ROOT_AGE) + 1, 100):
    rr = rad(ROOT_AGE - age)
    aa = np.linspace(math.radians(90 - GAP_DEG / 2 - span),
                     math.radians(90 - GAP_DEG / 2), 240)
    axa.plot(rr * np.cos(aa), rr * np.sin(aa), color="#eeece7", lw=0.35, zorder=0)

def draw(n):
    if not n["kids"]: return
    th = [k["_th"] for k in n["kids"]]
    r0 = rad(n["_r"])
    a = np.linspace(min(th), max(th), 80)
    axa.plot(r0 * np.cos(a), r0 * np.sin(a), color=bcol(n), lw=0.5,
             solid_capstyle="round", zorder=2)
    for k in n["kids"]:
        r1 = rad(k["_r"])
        axa.plot([r0 * math.cos(k["_th"]), r1 * math.cos(k["_th"])],
                 [r0 * math.sin(k["_th"]), r1 * math.sin(k["_th"])],
                 color=bcol(k), lw=(0.75 if not k["kids"] else 0.5),
                 solid_capstyle="round", zorder=2)
        draw(k)
draw(root)

dth = math.radians(span / NT)
R1I, R1O = RT + 0.030, RT + 0.078          # ring 1: codon class
R2I, R2O = RT + 0.090, RT + 0.122          # ring 2: in the record
for t in tips:
    sp = norm(t["label"]); th = t["_th"]
    cl = codon_class(sp)
    axa.add_patch(Wedge((0, 0), R1O, math.degrees(th - dth / 2),
                        math.degrees(th + dth / 2), width=R1O - R1I,
                        facecolor=CLASS_COLOUR.get(cl, NOCALL),
                        edgecolor="none", zorder=3))
    if sp in documented:
        axa.add_patch(Wedge((0, 0), R2O, math.degrees(th - dth / 2),
                            math.degrees(th + dth / 2), width=R2O - R2I,
                            facecolor=INK, edgecolor="none", zorder=4))

# outer arcs, one label per order at its longest run, never upside down
RA = RT + 0.146
runs, cur = [], None
for i, t in enumerate(tips):
    o = order_of.get(norm(t["label"]), "")
    o = o if o in TOP else "Other"
    if cur and cur[0] == o: cur[2] = i
    else:
        if cur: runs.append(cur)
        cur = [o, i, i]
if cur: runs.append(cur)
best = {}
for o, i0, i1 in runs:
    if o not in best or (i1 - i0) > (best[o][2] - best[o][1]):
        best[o] = (o, i0, i1)
for o, i0, i1 in runs:
    if (i1 - i0 + 1) < 4: continue
    t0, t1 = tips[i0]["_th"] + dth / 2, tips[i1]["_th"] - dth / 2
    a = np.linspace(t1, t0, 60)
    col = ORDER_COLOUR.get(o, "#b9b8b2")
    axa.plot(RA * np.cos(a), RA * np.sin(a), color=col, lw=2.2,
             solid_capstyle="butt", zorder=4)
    if best.get(o) != (o, i0, i1): continue
    tm = 0.5 * (t0 + t1); deg = math.degrees(tm) % 360
    flip = 90 < deg <= 270                       # left half: flip so text is upright
    axa.text((RA + 0.028) * math.cos(tm), (RA + 0.028) * math.sin(tm),
             ("Other orders" if o == "Other" else o),
             rotation=(deg + 180 if flip else deg), rotation_mode="anchor",
             ha=("right" if flip else "left"), va="center",
             fontsize=6.2, color=INK2)

LIM = RA + 0.30
axa.set_xlim(-LIM, LIM); axa.set_ylim(-LIM, LIM)

# time scale bar, outside the tree where nothing can cross it
sb = rad(ROOT_AGE) - rad(ROOT_AGE - 100.0)
# bottom RIGHT: the legend occupies the bottom left of this panel
xb, yb = LIM - 0.06 - sb, -LIM + 0.22
axa.plot([xb, xb + sb], [yb, yb], color=INK2, lw=1.0, solid_capstyle="butt")
for xx in (xb, xb + sb):
    axa.plot([xx, xx], [yb - 0.018, yb + 0.018], color=INK2, lw=1.0)
axa.text(xb + sb / 2.0, yb + 0.030, "100 Myr", ha="center", va="bottom",
         fontsize=6.0, color=INK2)
axa.text(0, -LIM + 0.06, "%d species  ·  root %.0f Ma" % (NT, ROOT_AGE),
         ha="center", va="bottom", fontsize=6.6, color=INK2)
# axa starts at figure x = -0.015, so an axes-fraction of 0.015 put the
# letter at figure x = -0.0045, just off the page.
panel_tag(axa, "a", x=0.040, y=0.955)

nreach = sum(1 for t in tips if codon_class(norm(t["label"])) == "reach")
nno = sum(1 for t in tips if codon_class(norm(t["label"])) == "noreach")
nal = sum(1 for t in tips if codon_class(norm(t["label"])) == "already")
nnc = sum(1 for t in tips if codon_class(norm(t["label"])) is None)
# the four classes partition the tips: if they ever stop doing so, the legend
# is describing a different figure from the one on the page
assert nreach + nno + nal + nnc == NT, (nreach, nno, nal, nnc, NT)
axa.legend(handles=[
    Line2D([], [], marker="s", ls="", ms=5, mfc=REACH, mec="none",
           label="Reaches L1014F in one step (%d)" % nreach),
    Line2D([], [], marker="s", ls="", ms=5, mfc=NOREACH, mec="none",
           label="Cannot reach it (%d)" % nno),
    Line2D([], [], marker="s", ls="", ms=5, mfc=ALREADY, mec="none",
           label="Already F (%d)" % nal),
    Line2D([], [], marker="s", ls="", ms=5, mfc=NOCALL, mec="none",
           label="No VGSC codon call (%d)" % nnc),
    Line2D([], [], marker="s", ls="", ms=5, mfc=INK, mec="none",
           label="In the resistance record (%d)" % len(documented))],
    # Three columns in the strip ABOVE the fan, not a single stacked column in
    # the upper-left corner. Stacked, the fifth entry ran into the tree: the
    # corner is only empty near the very top, and by the fifth line the circle
    # has widened back under the text. Measured at the fifth entry, the nearest
    # corner sat at r^2 = 0.944 of the outer radius, i.e. inside it. The strip
    # above the fan is free of both the circle and the panel letter.
    loc="lower left", bbox_to_anchor=(0.105, 0.900), frameon=False, ncol=3,
    fontsize=5.8, handletextpad=0.4, columnspacing=1.0, labelspacing=0.30,
    labelcolor=INK2)

# ---- b: species per order ------------------------------------------------ #
axb = fig.add_axes([0.775, 0.700, 0.205, 0.255])
labs = TOP + ["Other"]
vals = [ocount[o] for o in TOP] + [sum(v for k, v in ocount.items() if k not in TOP)]
ypos = np.arange(len(labs))[::-1]
axb.barh(ypos, vals, height=0.62,
         color=[ORDER_COLOUR.get(o, ORDER_HUES[7]) for o in labs])
for yy, v in zip(ypos, vals):
    axb.text(v + max(vals) * 0.03, yy, str(v), va="center", ha="left",
             fontsize=6.0, color=INK2)
axb.set_yticks(ypos); axb.set_yticklabels(labs, fontsize=6.2, color=INK)
axb.set_xlim(0, max(vals) * 1.22)
axb.set_xticks([]); axb.tick_params(length=0)
for s in ("top", "right", "bottom", "left"): axb.spines[s].set_visible(False)
axb.patch.set_visible(False)
axb.set_xlabel("Species", fontsize=6.4, color=INK2, labelpad=1)
panel_tag(axb, "b", x=-0.50, y=1.04)

# ---- c: the mechanism, computed from the genetic code -------------------- #
# Narrower than axb and axd on purpose: the per-codon species counts are
# drawn to the RIGHT of the matrix, and at a width of 0.205 they started at
# figure x = 0.986 and ran off the page.
axc = fig.add_axes([0.775, 0.345, 0.170, 0.270])
axc.set_xlim(0, len(focal_alleles)); axc.set_ylim(0, len(LEU))
axc.invert_yaxis(); axc.axis("off"); axc.patch.set_visible(False)
counts = collections.Counter(focal_codon.values())
for i, c in enumerate(LEU):
    for j, al in enumerate(focal_alleles):
        ok = al in REACH_MAP[c]
        axc.add_patch(Rectangle((j + 0.06, i + 0.10), 0.88, 0.80,
                                facecolor=(REACH if ok else NOREACH),
                                edgecolor="none"))
        if ok:
            # Drawn, not typed. Arial and Helvetica have no U+2713 CHECK MARK,
            # so "\u2713" rendered as an empty box here. y is inverted on this
            # axis, so the vertex is the LARGEST y of the three points.
            axc.plot([j + 0.30, j + 0.44, j + 0.72],
                     [i + 0.48, i + 0.66, i + 0.30],
                     color="#ffffff", lw=1.0, solid_capstyle="round",
                     solid_joinstyle="miter", zorder=5)
    axc.text(-0.12, i + 0.5, c, ha="right", va="center", fontsize=6.6,
             color=INK, family="monospace")
    axc.text(len(focal_alleles) + 0.12, i + 0.5, "%d" % counts.get(c, 0),
             ha="left", va="center", fontsize=6.0, color=INK2)
for j, al in enumerate(focal_alleles):
    axc.text(j + 0.5, -0.12, al, ha="center", va="bottom", fontsize=6.8, color=INK)
axc.text(-0.12, -0.12, "Codon", ha="right", va="bottom", fontsize=6.2,
         color=INK2, style="italic")
axc.text(len(focal_alleles) + 0.12, -0.12, "spp.", ha="left", va="bottom",
         fontsize=6.2, color=INK2, style="italic")
axc.text(len(focal_alleles) / 2.0, len(LEU) + 0.45,
         "%s %s → resistance allele" % (FOCAL_GENE, FOCAL_SITE),
         ha="center", va="top", fontsize=6.4, color=INK2)
panel_tag(axc, "c", x=-0.50, y=1.06)

# ---- d: curated events per gene ------------------------------------------ #
# Horizontal bars, not vertical. The counts run from 71 to 1, and in a vertical
# bar chart the three single-event genes were invisible slivers with orphaned
# labels floating above them. Horizontal gives every label a place to sit and
# matches panel b's idiom, so the two read as one system.
axd = fig.add_axes([0.775, 0.045, 0.205, 0.205])
core = list(csv.DictReader(open(CORE)))
# GluCl_Px and GluCl_Tu are two sites on ONE GluCl gene (ortho01: the same
# protein accession in 236 of 569 species). This panel is per GENE, so the two
# labels are merged -- otherwise Figure 1 would show seven genes while the
# paper, canonical_quantities.csv and Figures 4 and S1 all say six.
per = collections.defaultdict(lambda: [0, 0])          # gene -> [accessible, breach]
for e in core:
    g = "GluCl" if e["gene"].startswith("GluCl") else e["gene"]
    per[g][0 if istrue(e.get("accessible")) else 1] += 1
genes = sorted(per, key=lambda g: (per[g][0] + per[g][1]))
ypos = np.arange(len(genes))
acc = np.array([per[g][0] for g in genes], float)
brc = np.array([per[g][1] for g in genes], float)
tot = acc + brc
axd.barh(ypos, acc, height=0.62, color=ACCESSIBLE)
axd.barh(ypos, brc, height=0.62, left=acc, color=BREACH)
for yy, a_, b_, t_ in zip(ypos, acc, brc, tot):
    lab = "%d" % t_ + ("  (%d)" % b_ if b_ else "")
    axd.text(t_ + tot.max() * 0.035, yy, lab, va="center", ha="left",
             fontsize=6.0, color=INK2)
axd.set_yticks(ypos); axd.set_yticklabels(genes, fontsize=6.2, color=INK)
axd.set_xlim(0, tot.max() * 1.26)
axd.set_xticks([]); axd.tick_params(length=0)
for sp_ in ("top", "right", "bottom", "left"): axd.spines[sp_].set_visible(False)
axd.patch.set_visible(False)
axd.set_xlabel("Curated events", fontsize=6.4, color=INK2, labelpad=1)
axd.legend(handles=[
    Line2D([], [], marker="s", ls="", ms=5, mfc=ACCESSIBLE, mec="none",
           label="Accessible"),
    Line2D([], [], marker="s", ls="", ms=5, mfc=BREACH, mec="none",
           label="Breach (n)")],
    loc="lower right", bbox_to_anchor=(1.03, 0.98), frameon=False, ncol=2,
    fontsize=6.0, handletextpad=0.4, columnspacing=0.9, labelcolor=INK2)
panel_tag(axd, "d", x=-0.50, y=1.10)

# --------------------------------------------------------------------------- #
#  SOURCE DATA
#
#  eLife and Nature Communications both expect a source data file per figure: the
#  numbers behind every panel, so a reader can rebuild it without the code. It is
#  written here, from the same variables the panels are drawn from, so it cannot
#  drift from what is plotted. One sheet per panel if openpyxl is available,
#  otherwise one CSV per panel.
# --------------------------------------------------------------------------- #
SD = os.path.join(OUTDIR, "source_data")
os.makedirs(SD, exist_ok=True)

panel_a = [["tip_label", "order", "codon_VGSC_L1014", "codon_class",
            "in_resistance_record", "root_to_tip_Myr"]]
for t in tips:
    spn = norm(t["label"])
    panel_a.append([t["label"], order_of.get(spn, ""), focal_codon.get(spn, ""),
                    codon_class(spn) or "", int(spn in documented), round(t["_r"], 4)])
panel_b = [["order", "n_species"]] + [[o, ocount[o]] for o in TOP] + \
          [["Other", sum(v for k, v in ocount.items() if k not in TOP)]]
panel_c = [["leucine_codon", "n_species"] + ["reaches_" + a for a in focal_alleles]]
for c in LEU:
    panel_c.append([c, counts.get(c, 0)] +
                   [int(a in REACH_MAP[c]) for a in focal_alleles])
panel_d = [["gene", "n_events", "n_accessible", "n_breach"]]
for g in sorted(per, key=lambda g: -(per[g][0] + per[g][1])):
    panel_d.append([g, per[g][0] + per[g][1], per[g][0], per[g][1]])

SHEETS = [("a_phylogeny_and_rings", panel_a), ("b_species_per_order", panel_b),
          ("c_codon_reachability", panel_c), ("d_events_per_gene", panel_d)]
wrote = []
try:
    from openpyxl import Workbook
    wb = Workbook(); wb.remove(wb.active)
    for name, data in SHEETS:
        ws = wb.create_sheet(title=name[:31])
        for row in data: ws.append(row)
    fp = os.path.join(SD, "Fig1_source_data.xlsx")
    wb.save(fp); wrote.append(fp)
except ImportError:
    for name, data in SHEETS:
        fp = os.path.join(SD, "Fig1_source_data_%s.csv" % name)
        with open(fp, "w", newline="") as fh:
            csv.writer(fh).writerows(data)
        wrote.append(fp)

S.write_legend_data("fig01_dataset_and_context", {
    "n_tips": NT,
    "root_age_ma": ROOT_AGE,
    "focal_gene": FOCAL_GENE,
    "focal_site": FOCAL_SITE,
    "target_residue": TARGET,
    "ring1_reach": nreach,
    "ring1_noreach": nno,
    # Panel c tallies codons over the species in the TREE that carry a
    # wild-type codon here, which is fewer than the whole panel does, because
    # TimeTree resolves only part of the atlas. Both counts are emitted so the
    # legend can name the difference instead of leaving two numbers for the
    # same-looking quantity in two figures.
    "n_wildtype_focal_tree": nreach + nno,
    "n_wildtype_focal_panel": len({r["species"] for r in atlas
                                   if r["gene"] == FOCAL_GENE
                                   and r["site_label"] == FOCAL_SITE
                                   and r["is_wildtype"] == "True"}),
    "ring1_already": nal,
    "ring1_nocall": nnc,
    # Named for what it is: species present anywhere in the record, across
    # all sites. It is NOT the count at the focal site.
    "ring2_in_record_any_site": len(documented),
    "orders_shown": {o: int(ocount[o]) for o in TOP},
    "orders_other": int(sum(v for k, v in ocount.items() if k not in TOP)),
    "n_orders_total": len(ocount),
})

os.makedirs(OUTDIR, exist_ok=True)
for ext, kw in (("pdf", {}), ("png", dict(dpi=600))):
    p = os.path.join(OUTDIR, "fig01_dataset_and_context." + ext)
    # bbox_inches must stay None: "tight" discards the declared figsize and
    # wrote this figure at 189.2 mm against a 183 mm column.
    fig.savefig(p, facecolor=SURFACE, bbox_inches=None, pad_inches=0.0, **kw)
    if ext == "pdf":
        S.assert_no_type3(p)
        S.audit(fig, p)
    print("wrote %s" % os.path.abspath(p))
for f in wrote: print("wrote %s" % os.path.abspath(f))
print()
print("For the caption (these are NOT printed on the figure):")
print("  %d curated events, %d species, %d references"
      % (len(core), len({e["species"] for e in core}),
         len({(e.get("primary_ref") or "").strip() for e in core})))

print()
print("Panel c, computed from the genetic code (verify against the figure):")
for c in LEU:
    r = sorted(REACH_MAP[c] & set(focal_alleles))
    print("   %s (%3d spp.)  reaches: %s" % (c, counts.get(c, 0),
                                             ", ".join(r) if r else "none"))
print()
print("Ring 1 classes: reach %d, cannot reach %d, already %s %d, no call %d"
      % (nreach, nno, TARGET, nal, nnc))
print("Ring 2: %d species in the resistance record" % len(documented))
print("Tree: %d tips, root %.2f Myr" % (NT, ROOT_AGE))
