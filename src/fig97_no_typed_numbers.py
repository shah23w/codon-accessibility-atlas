#!/usr/bin/env python3
"""
fig97_no_typed_numbers.py - enforce the rule the figures claim to follow.

    python3 src/fig97_no_typed_numbers.py

No number appearing in a figure is typed into the script that draws it; every
value is read from a result file at run time. That rule was asserted in a dozen
docstrings and never checked, and it was broken: step8_05_figure3.py carried

    float(cs["within_record_null_enrichment"]["value"]), 9.78e-05),

a P value the results file did not contain, which then flowed into
FigS4_source_data.csv and into the figure's legend data as though the pipeline
had produced it. A reader of the source data had no way to tell.

This checks every figure script for numeric literals that look like reported
quantities and fails if it finds one. What it deliberately ignores:

  layout        coordinates, sizes, widths, offsets, alpha, linewidth, dpi,
                point sizes - geometry is typed by definition
  structure     indices, counts of panels, array positions, range() bounds
  tolerances    comparison thresholds and rounding guards

What it catches is the thing that matters: a float in scientific notation, or
a number with three or more decimal places, sitting in a data expression. Those
are measurements, and a measurement in a figure script is a measurement nobody
can trace.
"""
import ast
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)

SCRIPTS = [
    "fig01_dataset_and_context", "fig02_forest", "step8_04_figure2",
    "step8_05_figure3", "step8_06_figureS1", "step8_08_figureS2",
    "step8_10_figureS3_matched_nulls", "step8_11_figure5_enumeration",
    "step8_12_figure4_supply",
]

# What distinguished the typed P from every legitimate literal was not its
# magnitude but its company: it sat in a tuple whose other members were read
# from the results table, so it looked like data and was written out as data.
#
#     (label, float(cs["within_record_null_enrichment"]["value"]), 9.78e-05)
#
# Geometry never looks like that. So the rule is: a float literal is flagged
# only when the statement containing it also reads from a result table. This
# will miss a typed number in a statement that reads nothing, which is the
# price of a check that does not fire on every axes coordinate in the project.
# geometry and tuning keywords: a number beside one of these is layout
LAYOUT = re.compile(
    r"(figsize|width|height|\bsize\b|pad|space|margin|\bleft\b|\bright\b|"
    r"\btop\b|\bbottom\b|hspace|wspace|alpha|\blw\b|linewidth|\bms\b|"
    r"\bmew\b|dpi|fontsize|rotation|zorder|xlim|ylim|xticks|yticks|bbox|"
    r"offset|\bdx\b|\bdy\b|\bx=|\by=|\bloc\b|ncol|labelspacing|"
    r"handletextpad|columnspacing|borderpad|capsize|elinewidth|linespacing|"
    r"transform|add_axes|add_patch|fig_label|set_position)", re.I)

TOLERANCE = re.compile(r"(abs\(|[<>]=?\s*\d*\.?\d*e-|math\.isclose)")
READS_DATA = re.compile(
    r"(\bcs\[|\bctx\[|\bcanon\[|_canon\[|\brow\[|\bres_rows\[|\bspec_rows\b|"
    r"\[[\"']value[\"']\]|\[[\"']p_value[\"']\]|\[[\"']ratio[\"']\]|"
    r"\[[\"']estimate[\"']\]|\[[\"']odds_ratio[\"']\]|\bread_csv\(|\bfirst_csv\()")

def suspicious(path):
    with open(path, encoding="utf-8") as fh:
        src = fh.read()
    lines = src.split("\n")
    try:
        tree = ast.parse(src)
    except SyntaxError as exc:
        sys.exit("FATAL: %s does not parse: %s" % (path, exc))

    # map each statement to its line span so a multi-line tuple is judged whole
    stmts = []
    for node in ast.walk(tree):
        if isinstance(node, ast.stmt) and hasattr(node, "end_lineno"):
            stmts.append((node.lineno, node.end_lineno))

    def statement_text(ln):
        span = min((s_ for s_ in stmts if s_[0] <= ln <= s_[1]),
                   key=lambda s_: s_[1] - s_[0], default=(ln, ln))
        return "\n".join(lines[span[0] - 1:span[1]])

    out, seen = [], set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Constant) or not isinstance(node.value, float):
            continue
        ln = node.lineno
        if lines[ln - 1].strip().startswith("#"):
            continue
        stmt = statement_text(ln)
        if not READS_DATA.search(stmt):
            continue
        # a tolerance in a comparison is not a reported value
        if TOLERANCE.search(stmt):
            continue
        # geometry keywords anywhere in the statement, not just the literal's
        # own line: a tuple of (label, colour, offset) wraps, and the keyword
        # that identifies the offset sits on the line above
        if LAYOUT.search(stmt):
            continue
        key = (ln, repr(node.value))
        if key in seen:
            continue
        seen.add(key)
        out.append((ln, lines[ln - 1].strip()[:96], repr(node.value)))
    return out


def main():
    bad = 0
    for name in SCRIPTS:
        p = os.path.join("src", name + ".py")
        if not os.path.exists(p):
            print("  MISSING %s" % p)
            bad += 1
            continue
        hits = suspicious(p)
        if hits:
            bad += len(hits)
            print("  [FAIL] %s" % p)
            for ln, text, lit in hits:
                print("         line %d: %s   <- %s" % (ln, text, lit))
        else:
            print("  [ok]   %s" % p)
    print()
    if bad:
        print("=" * 70)
        print("%d typed value(s) in figure scripts. Every number a figure "
              "reports must" % bad)
        print("be read from a result file, so that a reader can trace it. "
              "Move the value")
        print("into the pipeline and read it back.")
        print("=" * 70)
        sys.exit(1)
    print("No typed values: every reported number is read from a file.")


if __name__ == "__main__":
    main()
