#!/usr/bin/env python3
"""
fig00_style.py: the shared figure style and the specification checks.

CANVAS
------
Every figure is laid out for the Nature-family canvas:

    Nature family   89 mm (1 col)   183 mm (2 col)   max height 170 mm
    Nat Commun      88 mm           180 mm
    eLife           89 mm           183 mm

A figure supplied at one journal's width and then rescaled by production is a
figure whose type size is no longer the size you chose. So the widths below are
set to the NARROWEST of the three, and every figure is built at final size:

    W1  = 89 mm    fits all three as a single-column figure
    W15 = 120 mm   eLife 1.5-column; comfortably inside Nature's 120 mm
    W2  = 183 mm   fits all three as a double-column figure

Nothing is scaled after export. Submitting the 183 mm file to Nature leaves a
9 mm gutter and costs nothing; submitting a 183 mm file to eLife would shrink
7 pt type to 6.6 pt without anyone deciding to.

TYPE
----
All three journals require sans-serif and a minimum size AFTER reduction of
5 pt (Nature) to 6 pt (eLife). Body text here is 7 pt, tick labels 6.5 pt, and
nothing drops below 6 pt. Fonts are embedded as TrueType (pdf.fonttype 42) so
text stays selectable and editable in production.

Arial is the house face. The cluster may not have it, and matplotlib falls back
to DejaVu Sans SILENTLY - a figure that looks right to the script and wrong to
the journal. resolve_font() therefore reports which face was actually used and
says so out loud when it had to fall back. It does not guess from the request.

PALETTE
-------
The palette was checked rather than trusted: worst normal-vision OKLab dE 24.0
(floor 15), worst CVD dE 9.2
across protan / deutan / tritan (target 8). Aqua sits at 2.82:1 against white,
below the 3:1 contrast gate, so anything drawn in aqua carries a direct label -
identity is never left to colour alone.

CONVENTIONS, so every panel reads as one system
-----------------------------------------------
  * grid and spines recessive; no top/right spines; no chartjunk
  * one axis per panel - never a second y-scale
  * text uses an ink token, never a series colour; a coloured mark beside a
    label carries the identity instead
  * uncertainty is always drawn, never implied
  * effect sizes on a log axis, with the null line at 1 drawn but not labelled
    as if it were data
  * selective direct labels, never a number on every mark
"""
from __future__ import annotations

import os

import matplotlib as mpl
from matplotlib.patches import Polygon

MM = 1 / 25.4

# ---- canvas widths, in inches ----------------------------------------------
# Nature Ecology & Evolution follows the Nature research-figure specification:
#   single column   89 mm
#   1.5 column     120 mm
#   double column  183 mm
#   maximum height 170 mm, leaving room for the legend beneath
# A figure drawn narrower than the specification is scaled up at typesetting,
# which shrinks every font, so these are the specification widths exactly.
W1 = 89 * MM             # 3.50"  single column
W15 = 120 * MM           # 4.72"  1.5 column
W2 = 183 * MM            # 7.20"  double column
H_MAX = 170 * MM         # 6.69"  page-height ceiling for one figure

# ---- validated categorical slots (light surface) ---------------------------
BLUE = "#2a78d6"         # slot 1
ORANGE = "#eb6834"       # slot 2
AQUA = "#1baf7a"         # slot 3 - 2.82:1, direct-label it

# ---- ink and surface -------------------------------------------------------
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e4e1"
RULE = "#c9c9c5"
SURFACE = "#ffffff"

# semantic roles, so a colour means the same thing in every panel
C_OBS = BLUE             # what was observed
C_EXP = MUTED            # what the null expects
C_ALT = ORANGE           # the contrasting group
C_THIRD = AQUA
C_POOLED = INK           # a pooled / summary estimate is ink, not a series colour

FONT_STACK = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans",
              "FreeSans", "DejaVu Sans"]

# Fonts kept with the project, so a figure does not depend on what happens to
# be installed on whichever machine draws it. assets/fonts/ holds Liberation
# Sans, which is metrically compatible with Arial and OFL-licensed.
import os as _os
_FONTDIR = _os.path.join(_os.path.dirname(_os.path.dirname(
    _os.path.abspath(__file__))), "assets", "fonts")


def register_local_fonts(verbose=True):
    """Make any font shipped in assets/fonts/ visible to matplotlib."""
    if not _os.path.isdir(_FONTDIR):
        return 0
    from matplotlib import font_manager
    n = 0
    for fn in sorted(_os.listdir(_FONTDIR)):
        if fn.lower().endswith((".ttf", ".otf")):
            try:
                font_manager.fontManager.addfont(_os.path.join(_FONTDIR, fn))
                n += 1
            except Exception as e:                           # noqa: BLE001
                if verbose:
                    print("[fig00_style] could not load %s: %s" % (fn, e))
    if n and verbose:
        print("[fig00_style] loaded %d font file(s) from assets/fonts/" % n)
    return n


def resolve_font(verbose=True):
    """
    Say which sans-serif face matplotlib will ACTUALLY use, and warn on fallback.

    matplotlib resolves font.sans-serif silently: ask for Arial on a machine
    without it and you get DejaVu Sans with no message. DejaVu is wider than
    Arial at the same point size, so a panel laid out against it will not match
    the journal's metrics. This checks the installed set rather than the request.
    """
    register_local_fonts(verbose=verbose)
    from matplotlib import font_manager
    installed = {f.name for f in font_manager.fontManager.ttflist}
    for name in FONT_STACK:
        if name in installed:
            if verbose:
                if name in ("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans"):
                    note = "OK" if name in ("Arial", "Helvetica") else \
                           "OK - metrically compatible with Arial"
                else:
                    note = ("FALLBACK - not Arial-metric; type will not match "
                            "the journal spec. Put an Arial-metric font in "
                            "assets/fonts/")
                print("[fig00_style] font in use: %s  (%s)" % (name, note))
            return name
    if verbose:
        print("[fig00_style] WARNING: none of %s is installed; matplotlib default in use"
              % ", ".join(FONT_STACK))
    return None


def apply(verbose=True):
    resolve_font(verbose=verbose)
    mpl.rcParams.update({
        "figure.dpi": 300,
        "savefig.dpi": 600,          # raster companion well above the 300 dpi floor
        # NOT "tight". A tight bbox re-crops the canvas to whatever the artists
        # happen to occupy, so the written file is no longer the width that was
        # chosen - which is the one thing a figure supplied at final size must
        # be. Panels are placed explicitly instead; see save().
        "savefig.bbox": None,
        "savefig.pad_inches": 0.0,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "font.size": 7,
        "axes.labelsize": 7,
        # Nature's specification is that EVERY font in the figure sits between
        # 5 pt and 7 pt. 7.5 was over that ceiling.
        "axes.titlesize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.6,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.minor.size": 1.4,
        "ytick.minor.size": 1.4,
        "text.color": INK,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "legend.frameon": False,
        "lines.linewidth": 1.1,
        "lines.solid_capstyle": "round",
        "pdf.fonttype": 42,          # embed TrueType; production rejects Type 3
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


# --------------------------------------------------------------------------
# export
# --------------------------------------------------------------------------
def save(fig, stem, outdir, formats=("pdf", "png"), expect_mm=None):
    """
    Write the figure, then CHECK the file rather than assume it is right.

    bbox_inches is forced to None so the PDF media box is exactly figsize. If
    expect_mm is given as (width, height) in millimetres, the written PDF's own
    MediaBox is measured and compared - because "I set figsize" is a request and
    the media box is the evidence, and a figure that reaches production 4 mm
    wide of the column is a figure production will rescale.
    """
    os.makedirs(outdir, exist_ok=True)
    written = []
    for ext in formats:
        path = os.path.join(outdir, "%s.%s" % (stem, ext))
        fig.savefig(path, bbox_inches=None, pad_inches=0.0)
        written.append(path)
        if ext == "pdf":
            assert_no_type3(path)
            # from the live figure, before it is discarded: see "verification"
            audit(fig, path)
            w_mm, h_mm = pdf_size_mm(path)
            print("  wrote %s   media box %.1f x %.1f mm" % (path, w_mm, h_mm))
            if expect_mm:
                for got, exp, lab in ((w_mm, expect_mm[0], "width"),
                                      (h_mm, expect_mm[1], "height")):
                    if abs(got - exp) > 0.5:
                        raise SystemExit(
                            "%s: %s is %.2f mm, expected %.2f mm. The file is not "
                            "at final size." % (path, lab, got, exp))
        else:
            print("  wrote %s" % path)
    return written


def pdf_size_mm(path):
    """Width and height of the first MediaBox, in millimetres (1 pt = 1/72 in)."""
    import re
    with open(path, "rb") as fh:
        blob = fh.read()
    m = re.search(rb"/MediaBox\s*\[\s*([\d.+-]+)\s+([\d.+-]+)\s+"
                  rb"([\d.+-]+)\s+([\d.+-]+)\s*\]", blob)
    if not m:
        raise SystemExit("%s: no /MediaBox found; cannot verify the page size" % path)
    x0, y0, x1, y1 = (float(v) for v in m.groups())
    return (x1 - x0) * 25.4 / 72.0, (y1 - y0) * 25.4 / 72.0


def assert_no_type3(path):
    """
    Production at all three journals rejects Type 3 fonts.

    rcParams sets pdf.fonttype 42, which yields Type0/CIDFontType2 with an
    embedded TrueType program. This asserts the file that was actually written,
    because the rcParam is a request and the file is the evidence.
    """
    import re
    with open(path, "rb") as fh:
        blob = fh.read()
    if re.search(rb"/Subtype\s*/Type3", blob):
        raise SystemExit("%s: contains Type 3 fonts - production will reject it" % path)
    if not re.search(rb"/FontFile2", blob):
        raise SystemExit("%s: no embedded TrueType font program found" % path)
    return True


# --------------------------------------------------------------------------
# Nature Ecology & Evolution compliance
# --------------------------------------------------------------------------
# The specification, from Nature's research-figure guide, which NEE follows:
#   width   89 mm (single), 120 mm (1.5), or 183 mm (double column)
#   height  170 mm maximum, leaving room for the legend
#   fonts   every font between 5 pt and 7 pt, Arial or Helvetica, embedded as
#           TrueType (type 2 or 42), never outlined
#   format  vector PDF; .png and .tiff are not accepted for main figures
#
# Panel letters are the one documented departure. The guide fixes body text at
# 5-7 pt and is silent on panel letters, and Nature's own published figures set
# them in bold noticeably larger than the axis text. They are allowed here at a
# single declared size, which must be passed explicitly so it can never drift
# unnoticed.
NEE_WIDTHS_MM = (89.0, 120.0, 183.0)
NEE_MAX_HEIGHT_MM = 170.0
NEE_MIN_PT = 5.0
NEE_MAX_PT = 7.0
PANEL_LETTER_PT = 8.0


# --------------------------------------------------------------------------
# verification
#
# These checks run against the LIVE FIGURE at save time, not against the
# written PDF. That is a correction, not a convenience.
#
# The earlier implementation read the PDF with PyMuPDF, and it could not work.
# PyMuPDF clips text extraction at the page boundary: a run of text that falls
# off the page comes back TRUNCATED, with a bounding box that stops exactly at
# the edge. A title overflowing by 191 pt and one sitting 0.1 pt proud of the
# margin are indistinguishable in the extracted spans, and the words that were
# lost leave no trace at all. Measured on a deliberate 265 px overflow, the
# reported box was 0.1 pt over the page and the string was silently cut from
# "INSIDE_A WAY_OUTSIDE_THE_PAGE_B..." to "INSIDE_A WAY_OUTS". Figure S2
# passed every check while printing "every stratum P <= 0.001" and dropping
# "(permutation)" entirely.
#
# The font-size and glyph checks read the same spans, so they shared the
# blindness: text that overflows is text they never saw, and the overflowing
# text is exactly the text most likely to be wrong.
#
# Matplotlib knows all three facts exactly, before anything is written. Every
# Text carries its own point size; get_window_extent gives its true rendered
# box in figure pixels however far outside the canvas it lands; and FT2Font
# answers directly whether the resolved face can draw a character.
#
# audit() records the answers in a sidecar beside the PDF so fig99 --verify
# can still report on files built earlier. A missing or stale sidecar is a
# FAILURE, never a pass: a check that cannot run must not be able to look like
# a check that ran and found nothing.
# --------------------------------------------------------------------------
CHECK_SUFFIX = ".check.json"


def _renderer(fig):
    """
    A renderer with the figure fully laid out.

    Axis labels, offset text and tick labels are POSITIONED during the draw,
    not when they are created, so measuring against a renderer that has not
    drawn reports stale coordinates. Measured that way an x label sitting
    comfortably inside a 60 mm page read as 0.5 pt off the bottom. A check
    that cries wolf is worse than no check, because it gets switched off.
    Drawing through Agg costs a few milliseconds and makes every extent exact.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return canvas.get_renderer()


def _live_texts(fig):
    """
    Every Text the reader will actually see.

    Tick labels belonging to ticks outside the current view limits are still
    live artists, still report get_visible() as True, and still carry real
    positions far off the canvas - 64 px off, in the case that prompted this.
    Matplotlib does not draw them. Counting them would make every figure with
    a set_xlim look like it was losing text, and a check that cries wolf gets
    switched off.
    """
    from matplotlib.text import Text as _Text
    skip = set()
    for ax in fig.get_axes():
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            eps = (hi - lo) * 1e-9
            for ticks, locs in (
                    (axis.get_major_ticks(), axis.get_majorticklocs()),
                    (axis.get_minor_ticks(), axis.get_minorticklocs())):
                for tick, loc in zip(ticks, locs):
                    if lo - eps <= loc <= hi + eps:
                        continue
                    skip.add(id(tick.label1))
                    skip.add(id(tick.label2))
    out = []
    for t in fig.findobj(_Text):
        if id(t) in skip or not t.get_visible():
            continue
        if not str(t.get_text()).strip():
            continue
        out.append(t)
    return out


def figure_text_sizes(fig):
    """Point size -> the text drawn at it, from the figure itself."""
    import collections as _c
    out = _c.defaultdict(list)
    for t in _live_texts(fig):
        out[round(float(t.get_fontsize()), 1)].append(str(t.get_text()))
    return dict(out)


def figure_overflow(fig, tol_pt=0.3):
    """
    Text that falls outside the page, as (text, side, points_over).

    Saving at an exact figsize is what makes this meaningful: with
    bbox_inches="tight" matplotlib grows the canvas until nothing overflows,
    hiding the problem by changing the page size instead of the layout. Once
    the width is pinned to the column, anything that does not fit is simply
    cut, and a reader sees "nformative" where the figure says "predicted
    informative" - or, worse, sees a sentence that ends early and reads as
    though it were complete.
    """
    r = _renderer(fig)
    W, H = fig.bbox.x1, fig.bbox.y1
    tol_px = tol_pt * fig.dpi / 72.0
    out = []
    for t in _live_texts(fig):
        bb = t.get_window_extent(r)
        for over, side in ((-bb.x0, "left"), (bb.x1 - W, "right"),
                           (-bb.y0, "bottom"), (bb.y1 - H, "top")):
            if over > tol_px:
                out.append((str(t.get_text()), side, over * 72.0 / fig.dpi))
                break
    return out


def _resolved_face():
    """The FT2Font matplotlib will actually draw with, not the request."""
    import matplotlib.font_manager as _fm
    from matplotlib.ft2font import FT2Font
    return FT2Font(_fm.findfont(
        _fm.FontProperties(family=mpl.rcParams["font.sans-serif"])))


def figure_missing_glyphs(fig):
    """
    Characters the figure asks for that the resolved font cannot draw.

    A missing glyph does not fail, warn, or look wrong in every viewer: it
    prints as an empty box. Figure 1 carried U+2713 CHECK MARK, which Arial
    and Helvetica do not contain, so every "reachable" cell showed a blank
    square. DejaVu Sans has it, which is why it looked right on a machine that
    fell through to DejaVu.
    """
    font = _resolved_face()
    chars = set()
    for t in _live_texts(fig):
        chars.update(str(t.get_text()))
    return sorted(c for c in chars
                  if ord(c) > 127 and font.get_char_index(ord(c)) == 0)


def _sidecar_path(pdf_path):
    return (pdf_path[:-4] if pdf_path.lower().endswith(".pdf")
            else pdf_path) + CHECK_SUFFIX


def audit(fig, pdf_path, verbose=False):
    """
    Check a figure the moment it is drawn, and record the result beside it.

    Call this immediately after writing the PDF, with the same figure object.
    Overflow and missing glyphs raise at once, because they are defects the
    reader sees and no later stage can recover the lost text. The sidecar is
    written FIRST, so a failing build still leaves the evidence on disk.
    """
    import json as _json
    over = figure_overflow(fig)
    glyphs = figure_missing_glyphs(fig)
    sizes = figure_text_sizes(fig)
    w_mm, h_mm = [v * 25.4 for v in fig.get_size_inches()]

    rec = {
        "pdf": os.path.basename(pdf_path),
        "page_mm": [round(w_mm, 2), round(h_mm, 2)],
        # the face actually resolved, not the stack that was requested
        "font": _resolved_face().family_name,
        "n_texts": sum(len(v) for v in sizes.values()),
        "sizes": {"%.1f" % k: sorted(set(v))[:6] for k, v in sizes.items()},
        "overflow": [[t, s, round(pt, 2)] for t, s, pt in over],
        "missing_glyphs": ["U+%04X %s" % (ord(c), c) for c in glyphs],
    }
    with open(_sidecar_path(pdf_path), "w", encoding="utf-8") as fh:
        _json.dump(rec, fh, indent=1, sort_keys=True, ensure_ascii=False)

    if over:
        raise SystemExit(
            "%s: %d text run(s) fall outside the page:\n    %s"
            % (pdf_path, len(over),
               "\n    ".join("%.1f pt off the %-6s %r" % (pt, s, t[:60])
                             for t, s, pt in over[:8])))
    if glyphs:
        raise SystemExit(
            "%s: the figure font has no glyph for %s - these print as empty "
            "boxes. Draw the mark instead of typing the character, or choose "
            "one Arial contains."
            % (pdf_path, ", ".join("%r (U+%04X)" % (c, ord(c)) for c in glyphs)))
    if verbose:
        print("  audit %s: %d texts, %.1f-%.1f pt, no overflow, no tofu"
              % (os.path.basename(pdf_path), rec["n_texts"],
                 min(sizes) if sizes else 0, max(sizes) if sizes else 0))
    return rec


def read_audit(path):
    """
    The sidecar for a written PDF, or a hard failure explaining why there
    is none. Never returns a silent empty result.
    """
    import json as _json
    side = _sidecar_path(path)
    if not os.path.exists(side):
        raise SystemExit(
            "%s has no %s beside it. The clipping, font-size and glyph checks "
            "are made from the live figure at save time and cannot be "
            "recovered from the PDF alone - PyMuPDF truncates any text that "
            "runs off the page, which is precisely the text in question. "
            "Rebuild the figure." % (path, CHECK_SUFFIX))
    if os.path.getmtime(side) < os.path.getmtime(path) - 1.0:
        raise SystemExit(
            "%s is older than the PDF beside it, so it describes an earlier "
            "build. Rebuild the figure." % side)
    with open(side, encoding="utf-8") as fh:
        rec = _json.load(fh)
    w_pdf, h_pdf = pdf_size_mm(path)
    w_rec, h_rec = rec.get("page_mm", [0, 0])
    if abs(w_pdf - w_rec) > 0.6 or abs(h_pdf - h_rec) > 0.6:
        raise SystemExit(
            "%s records a %.1f x %.1f mm page but the PDF media box is "
            "%.1f x %.1f mm. The sidecar is not from this file."
            % (side, w_rec, h_rec, w_pdf, h_pdf))
    return rec


def text_sizes(path):
    """Point size -> sample text, from the sidecar written at save time."""
    return {float(k): v for k, v in read_audit(path)["sizes"].items()}


def assert_no_overflow(path, tol_mm=0.3):
    bad = read_audit(path)["overflow"]
    if bad:
        raise SystemExit(
            "%s: %d text run(s) fall outside the page:\n    %s"
            % (path, len(bad),
               "\n    ".join("%.1f pt off the %-6s %r" % (pt, s, t[:60])
                             for t, s, pt in bad[:8])))
    return True


def assert_glyphs_present(path):
    bad = read_audit(path)["missing_glyphs"]
    if bad:
        raise SystemExit(
            "%s: the figure font has no glyph for %s - these print as empty "
            "boxes. Draw the mark instead of typing the character, or choose "
            "one Arial contains." % (path, ", ".join(bad)))
    return True


def assert_nee_compliant(path, panel_letter_pt=PANEL_LETTER_PT, verbose=True):
    """
    Assert a written figure meets the NEE specification, and say exactly what
    fails. Checks the FILE, not the intent: size, height, font embedding and
    every rendered text size.
    """
    problems = []

    w, h = pdf_size_mm(path)
    if not any(abs(w - target) <= 0.6 for target in NEE_WIDTHS_MM):
        problems.append(
            "width %.1f mm is not one of the permitted widths %s"
            % (w, ", ".join("%g" % t for t in NEE_WIDTHS_MM)))
    if h > NEE_MAX_HEIGHT_MM + 0.6:
        problems.append("height %.1f mm exceeds the %g mm maximum"
                        % (h, NEE_MAX_HEIGHT_MM))

    assert_no_type3(path)

    allowed = {round(float(panel_letter_pt), 1)} if panel_letter_pt else set()
    sizes = text_sizes(path)
    for pt in sorted(sizes):
        if pt in allowed:
            continue
        if pt < NEE_MIN_PT or pt > NEE_MAX_PT:
            sample = ", ".join(repr(s) for s in sorted(set(sizes[pt]))[:4])
            problems.append(
                "text at %.1f pt is outside %g-%g pt (%d run(s), e.g. %s)"
                % (pt, NEE_MIN_PT, NEE_MAX_PT, len(sizes[pt]), sample))

    if problems:
        raise SystemExit(
            "%s is not NEE compliant:\n    %s"
            % (path, "\n    ".join(problems)))
    if verbose:
        print("  PASS NEE compliant: %.1f x %.1f mm, text %.1f-%.1f pt"
              % (w, h,
                 min(s for s in sizes if s not in allowed) if sizes else 0,
                 max(s for s in sizes if s not in allowed) if sizes else 0))
    return True


# --------------------------------------------------------------------------
# panel furniture
# --------------------------------------------------------------------------
LEGEND_DIR = os.path.join("results", "figures", "legend_data")
SOURCE_DIR = os.path.join("results", "figures", "source_data")
SOURCE_HEADER = ("panel", "key", "value", "low", "high", "note")


def write_source_data(figure_label, rows):
    """
    The numbers behind one figure, named by FIGURE NUMBER.

    figure_label is "Fig2", "FigS1" and so on, never the script's output stem:
    figure2_necessary_not_sufficient.pdf is Figure 3 and figure2_forest.pdf is
    Figure 2, so anyone pairing source data to figures by stem pairs the wrong
    two.

    Rows are (panel, key, value, low, high, note); short rows are padded, so a
    panel with no interval can write four fields and a panel with one can
    write six.
    """
    import csv as _csv
    os.makedirs(SOURCE_DIR, exist_ok=True)
    path = os.path.join(SOURCE_DIR, "%s_source_data.csv" % figure_label)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh, lineterminator="\n")
        w.writerow(SOURCE_HEADER)
        for r in rows:
            r = list(r) + [""] * (len(SOURCE_HEADER) - len(r))
            w.writerow(r[:len(SOURCE_HEADER)])
    return path


def write_legend_data(stem, mapping):
    """
    Record the numbers this figure's legend will quote.

    A legend is a document, and this project does not type numbers into
    documents: the step that computes a number writes it down and the builder
    reads it. Several of these exist nowhere else - a figure script computes
    them, prints them to its log and discards them - so without this the
    legend would have to restate them from memory, which is how a caption
    comes to disagree with the panel above it.

    Values are written as given. Rounding is the legend builder's business,
    so that the same quantity rounds the same way everywhere it appears.
    """
    import json as _json
    os.makedirs(LEGEND_DIR, exist_ok=True)
    path = os.path.join(LEGEND_DIR, "%s.json" % stem)
    with open(path, "w", encoding="utf-8") as fh:
        _json.dump(mapping, fh, indent=1, sort_keys=True, ensure_ascii=False)
    return path


def panel_label(ax, letter, dx=-0.16, dy=1.04):
    """
    A bold lowercase panel letter at the one declared size.

    The house style is lowercase, and the case is FORCED here rather than left
    to each call site. The size is PANEL_LETTER_PT, the single value
    assert_nee_compliant() exempts from the 5-7 pt rule.
    """
    ax.text(dx, dy, str(letter).lower(), transform=ax.transAxes,
            fontsize=PANEL_LETTER_PT, fontweight="bold", va="top", ha="left",
            color=INK)


def fig_label(fig, letter, x, y):
    """Panel letter placed in FIGURE coordinates - for panels of unequal width,
    where axes-relative offsets would not line the letters up."""
    fig.text(x, y, str(letter).lower(), fontsize=PANEL_LETTER_PT,
             fontweight="bold", va="top", ha="left", color=INK)


def recessive_grid(ax, axis="y"):
    """
    Ticks instead of a grid.

    Nature's figure specification says plainly: "Avoid background gridlines."
    Every panel in this project drew one. A grid is a reading aid for looking
    up values, and these panels are read for pattern and effect size, not for
    lookup -- the exact numbers are in the source-data file. Longer tick marks
    on the relevant axis give the same orientation without laying a mesh over
    the data.

    The call sites are kept so the intent stays visible at each panel, and so
    the axis that WOULD have carried the grid still gets its ticks.
    """
    ax.set_axisbelow(True)
    ax.grid(False)
    for which in (("x", "y") if axis == "both" else (axis,)):
        ax.tick_params(axis=which, length=2.6, width=0.6, color=MUTED)


def strip_axes(ax):
    """A bare canvas for a text column sitting beside a plot."""
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.patch.set_alpha(0)


def fmt_p(p):
    p = float(p)
    if p < 1e-100:
        return "P < 1e-100"
    if p < 0.001:
        m, e = ("%.1e" % p).split("e")
        return "P = %s×10$^{%d}$" % (m, int(e))
    if p < 0.01:
        return "P = %.3f" % p
    return "P = %.2f" % p


def fmt_or(est, lo, hi):
    """An estimate the way a forest plot prints it, open bounds included."""
    import math

    def one(v):
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return "–"
        if math.isinf(v):
            return "∞"
        return "%.2f" % v if v < 100 else "%.0f" % v
    return "%s [%s, %s]" % (one(est), one(lo), one(hi))


# --------------------------------------------------------------------------
# forest-plot primitives
# --------------------------------------------------------------------------
def forest_axis(ax, xlim, null=1.0, ticks=None):
    """A log effect-size axis with the null drawn as furniture, not as data."""
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.axvline(null, color=RULE, linewidth=0.7, zorder=0)
    if ticks is None:
        ticks = [t for t in (0.1, 0.25, 0.5, 1, 2, 5, 10, 25, 50, 100, 250)
                 if xlim[0] <= t <= xlim[1]]
    ax.set_xticks(ticks)
    ax.set_xticklabels([("%g" % t) for t in ticks])
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def forest_row(ax, y, est, lo, hi, xlim, size=14, color=C_OBS, alpha=1.0):
    """
    One study row: interval, marker, and an ARROW where a bound is open.

    An unbounded interval is drawn as an arrow to the panel edge rather than
    clipped to it, because a clipped whisker reads as a finite bound.
    """
    import math
    lo_open = (lo is None) or (isinstance(lo, float) and (math.isnan(lo) or lo <= xlim[0]))
    hi_open = (hi is None) or (isinstance(hi, float) and (math.isinf(hi) or math.isnan(hi)
                                                          or hi >= xlim[1]))
    x0 = xlim[0] * 1.02 if lo_open else lo
    x1 = xlim[1] * 0.98 if hi_open else hi
    ax.plot([x0, x1], [y, y], color=color, linewidth=1.0, alpha=alpha,
            solid_capstyle="butt", zorder=2)
    for openb, xe, d in ((lo_open, x0, -1), (hi_open, x1, +1)):
        if openb:
            ax.annotate("", xy=(xe * (1.18 ** d), y), xytext=(xe, y),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.9,
                                        mutation_scale=4, shrinkA=0, shrinkB=0),
                        zorder=2)
        else:
            ax.plot([xe, xe], [y - 0.16, y + 0.16], color=color,
                    linewidth=0.9, alpha=alpha, zorder=2)
    if est is not None and math.isfinite(est):
        ax.scatter([est], [y], s=size, color=color, alpha=alpha,
                   zorder=3, linewidths=0, marker="s")
    else:
        ax.scatter([xlim[1] * 0.98], [y], s=size, color=color, alpha=alpha,
                   zorder=3, linewidths=0, marker=">")


def forest_diamond(ax, y, est, lo, hi, height=0.32, color=C_POOLED):
    """The pooled estimate, drawn as a diamond so it cannot be read as a study."""
    ax.add_patch(Polygon([[lo, y], [est, y + height], [hi, y], [est, y - height]],
                         closed=True, facecolor=color, edgecolor="none", zorder=4))


def cap(s):
    """Capitalise a label's first letter and leave the rest alone.

    str.capitalize() lowercases the remainder, which would turn "ESM-2
    tertile" into "Esm-2 tertile". Across the ten NEE research Articles
    surveyed, 89% of figure labels (2,249 of 2,535) begin with a capital; the
    lowercase remainder are symbols and units such as "log10 GFP" and
    "r = 0.9", not prose. Labels built from result files arrive lowercase, so
    they are capitalised here at the point of drawing rather than edited in
    the data.
    """
    if not s:
        return s
    return s[:1].upper() + s[1:]


def text_column(ax, y, parts, xs, ha=None, color=INK, size=6.5, weight="normal"):
    """
    Aligned numeric columns beside a forest - the table half of a forest plot.

    x positions are AXES fractions (0-1) so the columns stay put whatever the
    data range; y is in DATA units so a row lines up with its interval.
    """
    ha = ha or ["left"] * len(parts)
    trans = ax.get_yaxis_transform()          # x: axes fraction, y: data
    for txt, x, a in zip(parts, xs, ha):
        ax.text(x, y, txt, ha=a, va="center", fontsize=size, color=color,
                fontweight=weight, transform=trans, clip_on=False)
