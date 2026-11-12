#!/usr/bin/env python3
"""
fig99_build_and_verify.py - rebuild every figure, then check the written files.

    python3 src/fig99_build_and_verify.py            # build, then verify
    python3 src/fig99_build_and_verify.py --verify   # verify what is on disk

WHERE TO RUN IT. Either machine works, but all nine figures must be built on
the SAME one. NEE requires Arial or Helvetica. macOS resolves Arial from
/System/Library/Fonts/Supplemental; the cluster has no Arial but the project
ships Liberation Sans in assets/fonts/, which is Arial-metric and which
fig00_style.register_local_fonts() loads, so the cluster renders at the correct
metrics under a different font NAME. Mixing the two across a figure set would
embed two different fonts in one paper.

What must NOT happen is a silent fall-through to DejaVu Sans: it is wider than
Arial at the same point size, so line breaks and label overlaps differ from
what the journal will print, and it would mask a missing-glyph defect because
DejaVu contains characters Arial does not. fig00_style.resolve_font() reports
the font actually in use; read that line, do not assume.

Figures are light work either way: every script here runs in one to three
seconds. The step that cannot be automated is a person LOOKING at the rendered
PNGs -- the empty-box ticks and the clipped "Hymenoptera" below passed every
mechanical check until someone looked.

THREE CHECKS. Page size and font embedding are read from the written PDF.
Clipping, font sizes and glyph coverage are measured from the LIVE FIGURE at
save time and recorded in a .check.json beside each PDF, which this script
reads back.

  spec      183 mm (or 89 / 120) wide, at most 170 mm tall, every font between
            5 pt and 7 pt apart from the one declared panel-letter size,
            TrueType embedded and no Type 3.
  clipping  no text outside the page. Only meaningful because the figures are
            saved at an exact figsize: bbox_inches="tight" grows the canvas
            until nothing overflows, hiding the problem by changing the page
            size instead of the layout.
  glyphs    no character the resolved font cannot draw.

WHY NOT READ THE PDF FOR ALL THREE. The first version did, with PyMuPDF, and
the clipping check could not work. PyMuPDF clips text extraction at the page
boundary: text that runs off the page comes back truncated, with a bounding
box that stops exactly at the edge. On a deliberate 265 px overflow it
reported 0.1 pt and returned a string cut in half. Font sizes and glyphs read
the same spans, so all three were blind to overflowing text -- which is the
text most likely to be wrong. Two real defects survived a clean 9/9 run:
Figure S2 printed "every stratum P <= 0.001" with "(permutation)" cut off,
and Figure 3's panel letter a was laid at figure x = -0.03 and vanished.
Both are caught now, before anything is written.

A missing or stale .check.json is a FAILURE here, not a pass. A check that
cannot run must not look like a check that ran and found nothing.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.chdir(ROOT)

import fig00_style as S                                        # noqa: E402

# script -> the PDF stem it writes. Main figures first, in NEE order.
FIGURES = [
    ("Fig 1", "fig01_dataset_and_context", "fig01_dataset_and_context"),
    ("Fig 2", "fig02_forest", "figure2_forest"),
    ("Fig 3", "step8_04_figure2", "figure2_necessary_not_sufficient"),
    ("Fig 4", "step8_12_figure4_supply", "figure4_mutational_supply"),
    ("Fig 5", "step8_11_figure5_enumeration", "figure5_enumeration"),
    ("S1", "step8_06_figureS1", "figureS1_atlas_composition"),
    ("S2", "step8_08_figureS2", "figureS2_robustness"),
    ("S3", "step8_10_figureS3_matched_nulls", "figureS3_matched_null_contrast"),
    ("S4", "step8_05_figure3", "figureS4_mutation_spectrum"),
]
DIR_FIG = os.path.join("results", "figures")
# figure slot -> the source-data file that must exist beside it. Named by
# FIGURE NUMBER: the script stems disagree with the numbering (Fig 3's stem is
# "figure2_necessary_not_sufficient"), and pairing by stem pairs the wrong two.
SOURCE_DATA = {
    "Fig 1": "Fig1_source_data.xlsx",
    "Fig 2": "Fig2_source_data.csv",
    "Fig 3": "Fig3_source_data.csv",
    "Fig 4": "Fig4_source_data.csv",
    "Fig 5": "Fig5_source_data.csv",
    "S1": "FigS1_source_data.csv",
    "S2": "FigS2_source_data.csv",
    "S3": "FigS3_source_data.csv",
    "S4": "FigS4_source_data.csv",
}


def build():
    print("=" * 78)
    print("BUILD")
    print("=" * 78)
    bad = []
    for slot, script, _ in FIGURES:
        sys.stdout.write("  %-6s %-34s " % (slot, script))
        sys.stdout.flush()
        r = subprocess.run([sys.executable, os.path.join("src", script + ".py")],
                           capture_output=True, text=True)
        if r.returncode == 0:
            print("ok")
        else:
            print("FAILED")
            tail = (r.stdout + r.stderr).strip().splitlines()
            for line in tail[-6:]:
                print("            %s" % line[:100])
            bad.append(slot)
    return bad


def verify():
    print()
    print("=" * 78)
    print("VERIFY  (the written PDFs, not the scripts)")
    print("=" * 78)
    S.apply(verbose=False)
    # macOS writes AppleDouble sidecars on this external drive; they are not
    # PDFs, they break any reader that opens them, and they must never reach
    # the deposit.
    for f in os.listdir(DIR_FIG):
        if f.startswith("._"):
            os.remove(os.path.join(DIR_FIG, f))

    bad = []
    for slot, _, stem in FIGURES:
        path = os.path.join(DIR_FIG, stem + ".pdf")
        if not os.path.exists(path):
            print("  %-6s %-40s MISSING" % (slot, stem))
            bad.append(slot)
            continue
        problems = []
        for check in (S.assert_nee_compliant, S.assert_no_overflow,
                      S.assert_glyphs_present):
            try:
                if check is S.assert_nee_compliant:
                    check(path, verbose=False)
                else:
                    check(path)
            except SystemExit as exc:
                problems += [ln.strip() for ln in str(exc).split("\n")[1:]] \
                            or [str(exc).split(":", 1)[-1].strip()]
        src = os.path.join(DIR_FIG, "source_data", SOURCE_DATA[slot])
        if not os.path.exists(src):
            problems.append("no source-data file: %s. Every legend ends "
                            "\"Source data are provided with this paper\"."
                            % os.path.basename(src))
        elif os.path.getsize(src) < 80:
            problems.append("source-data file %s is empty"
                            % os.path.basename(src))
        w, h = S.pdf_size_mm(path)
        if problems:
            bad.append(slot)
            print("  %-6s %-40s %5.1f x %5.1f mm   FAIL" % (slot, stem, w, h))
            for p in problems:
                print("            %s" % p[:104])
        else:
            sizes = [s for s in S.text_sizes(path) if s != S.PANEL_LETTER_PT]
            print("  %-6s %-40s %5.1f x %5.1f mm   pass  text %.1f-%.1f pt"
                  % (slot, stem, w, h, min(sizes), max(sizes)))
    return bad


def typed_numbers():
    """
    No figure may type a value it reports. Run before the build, because a
    typed number produces a figure that looks right and a source-data file
    that is wrong, and neither the page checks nor a reader can tell.
    """
    print("=" * 78)
    print("PROVENANCE  (no reported value is typed into a figure script)")
    print("=" * 78)
    r = subprocess.run([sys.executable,
                        os.path.join("src", "fig97_no_typed_numbers.py")],
                       capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    if r.returncode != 0:
        sys.stdout.write(r.stderr)
    return [] if r.returncode == 0 else ["typed values"]


def main():
    verify_only = "--verify" in sys.argv
    failed = [] if verify_only else (typed_numbers() + build())
    failed += verify()
    print()
    print("=" * 78)
    if failed:
        print("RESULT: %d figure(s) not ready: %s"
              % (len(set(failed)), ", ".join(sorted(set(failed)))))
    else:
        print("RESULT: %d figures build and meet the NEE specification"
              % len(FIGURES))
    print("=" * 78)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
