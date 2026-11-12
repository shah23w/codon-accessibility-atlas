#!/usr/bin/env python3
"""
phylo03_collapse_subspecies.py
------------------------------
Collapse TimeTree's subspecies tips back to the species we asked for.

WHY

phylo02 found 69 tips in the TimeTree output that were never submitted. They were
not substitutions -- every one is a SUBSPECIES (a trinomial) of a species that WAS
submitted, and the parent binomial is correspondingly absent from the tree:

    submitted   Apis mellifera
    returned    Apis_mellifera_ligustica, Apis_mellifera_carnica, ... (13 tips)

So the tree is right; it is resolved one rank below what we need, and the join by
species name fails. Collapsing recovers 13 species that would otherwise be dropped,
including Apis mellifera and Culex pipiens.

HOW

For each parent species, keep ONE tip and prune the rest, then rename the survivor
to the binomial. The representative is the nominotypical subspecies when present
(Apis_mellifera_mellifera, Culex_pipiens_pipiens), because that is the type form;
otherwise the first alphabetically, which is arbitrary but recorded.

This is defensible at the scale that matters here. Subspecies of one species
diverged within the last few million years at most, against a root age of ~559 Ma,
so which subspecies represents the species changes its position in the tree by a
rounding error. Every choice is written to a CSV so the decision is
inspectable rather than implicit.

Pruning cannot distort the surviving tips' ages: removing a sibling does not change
a tip's own path to the root. Internal nodes left with a single child are collapsed
and their branch lengths summed, which is the standard operation and again leaves
root-to-tip distances untouched. phylo02 re-run on the output must still report the
tree ultrametric with the same root age -- that is the check that this script did
not quietly damage anything.

USAGE
    python src/phylo03_collapse_subspecies.py [tree_in] [species_list] [tree_out]
  defaults:
    data/refs/timetree_arthropods.nwk
    results/phylo01_species_list.txt
    data/refs/timetree_arthropods_collapsed.nwk

OUTPUTS
    the collapsed Newick
    results/phylo03_collapse_audit.csv
    logs/phylo03_collapse_subspecies.log
"""
import csv, os, sys

TREE_IN  = sys.argv[1] if len(sys.argv) > 1 else "data/refs/timetree_arthropods.nwk"
LIST     = sys.argv[2] if len(sys.argv) > 2 else "results/phylo01_species_list.txt"
TREE_OUT = sys.argv[3] if len(sys.argv) > 3 else "data/refs/timetree_arthropods_collapsed.nwk"
DIR_RES, DIR_LOG = "results", "logs"

_out = []
def say(s=""):
    print(s); _out.append(s)

# --------------------------------------------------------------------------- #
#  Newick parse / write.  Deliberately duplicated from phylo02_verify_tree.py:
#  phylo02 executes at import time, so it cannot be imported as a module, and a
#  40-line parser is cheaper than a refactor that would force phylo02 to be
#  re-verified. If either is changed, change both.
# --------------------------------------------------------------------------- #
class Newick:
    def __init__(self, text):
        self.s = text.replace("\n", "").replace("\r", "").replace("\t", "")
        i = self.s.rfind(";")
        if i >= 0: self.s = self.s[:i]
        self.i = 0
    def parse(self):
        return self._node()
    def _node(self):
        if self.i < len(self.s) and self.s[self.i] == "(":
            self.i += 1; kids = []
            while True:
                kids.append(self._node())
                if self.s[self.i] == ",": self.i += 1; continue
                if self.s[self.i] == ")": self.i += 1; break
                raise ValueError("unexpected %r at %d" % (self.s[self.i], self.i))
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

def write_newick(n, is_root=True):
    if n["kids"]:
        inner = ",".join(write_newick(k, False) for k in n["kids"])
        s = "(" + inner + ")" + (n["label"] or "")
    else:
        s = n["label"].replace(" ", "_")
    if not is_root and n["length"] is not None:
        s += ":%.10g" % n["length"]
    return s + (";" if is_root else "")

def tips_of(n, acc=None, depth=0.0):
    acc = [] if acc is None else acc
    d = depth + (n["length"] or 0.0)
    if n["kids"]:
        for k in n["kids"]: tips_of(k, acc, d)
    else:
        acc.append((n, d))
    return acc

def norm(s): return " ".join(s.replace("_", " ").split()).strip().lower()

def prune(n, keep):
    """return a copy containing only kept tips; collapse unary nodes, summing
    lengths. Returns None if nothing under this node is kept."""
    if not n["kids"]:
        return dict(n) if id(n) in keep else None
    kids = [prune(k, keep) for k in n["kids"]]
    kids = [k for k in kids if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        c = kids[0]
        c["length"] = (c["length"] or 0.0) + (n["length"] or 0.0)
        return c
    return dict(kids=kids, label=n["label"], length=n["length"])

# --------------------------------------------------------------------------- #
say("=" * 78)
say("phylo03  collapse subspecies tips to the submitted species")
say("=" * 78)
for p in (TREE_IN, LIST):
    if not os.path.exists(p):
        say("FATAL: missing %s" % p); sys.exit(2)

root = Newick(open(TREE_IN).read()).parse()
tips = tips_of(root)
say("  input tree : %s" % TREE_IN)
say("  tips in    : %d" % len(tips))

submitted = [ln.strip() for ln in open(LIST) if ln.strip()]
sub_norm = {norm(x): x for x in submitted}
say("  submitted  : %d species" % len(submitted))

direct, groups = {}, {}
for node, _d in tips:
    k = norm(node["label"])
    if k in sub_norm:
        direct[k] = node
        continue
    w = k.split()
    if len(w) >= 3 and " ".join(w[:2]) in sub_norm:
        groups.setdefault(" ".join(w[:2]), []).append(node)

say()
say("  tips matching a submitted species directly : %d" % len(direct))
say("  parent species represented only by subspecies : %d" % len(groups))

# choose one representative per group
audit, keep = [], set(id(n) for n in direct.values())
for parent, members in sorted(groups.items()):
    members = sorted(members, key=lambda n: norm(n["label"]))
    nomino = [m for m in members
              if norm(m["label"]).split()[-1] == parent.split()[1]]
    rep = nomino[0] if nomino else members[0]
    keep.add(id(rep))
    for m in members:
        audit.append(dict(parent_species=sub_norm[parent],
                          subspecies_tip=m["label"],
                          action=("KEPT as representative" if m is rep else "pruned"),
                          rule=("nominotypical subspecies" if nomino
                                else "first alphabetically (no nominotypical form present)")))
    say("    %-32s %2d tips -> %s" % (sub_norm[parent], len(members), rep["label"]))

dropped_tips = [n["label"] for n, _ in tips if id(n) not in keep]
say()
say("  tips kept    : %d" % len(keep))
say("  tips pruned  : %d" % len(dropped_tips))

new_root = prune(root, keep)
if new_root is None:
    say("FATAL: pruning removed everything"); sys.exit(2)
new_root["length"] = None

# rename representatives to the binomial
renamed = 0
for node, _d in tips_of(new_root):
    k = norm(node["label"])
    if k in sub_norm:
        node["label"] = sub_norm[k]; continue
    w = k.split()
    if len(w) >= 3 and " ".join(w[:2]) in sub_norm:
        node["label"] = sub_norm[" ".join(w[:2])]; renamed += 1
say("  tips renamed to binomial : %d" % renamed)

# ---- the safety check: ages must not have moved ------------------------- #
d_in  = [d for _, d in tips]
d_out = [d for _, d in tips_of(new_root)]
say()
say("-" * 78)
say("SAFETY CHECK: did pruning move any ages?")
say("-" * 78)
say("  root-to-tip before : min %.6f  max %.6f" % (min(d_in), max(d_in)))
say("  root-to-tip after  : min %.6f  max %.6f" % (min(d_out), max(d_out)))
# Tolerance is 1e-3 Ma -- one thousand years against a root of ~559 million.
# Not 1e-6: branch lengths are re-serialised as text, and that rounding accumulates
# down a deep tree, so a tighter bound would fail on arithmetic noise rather than on
# any real change. 1e-3 is far below anything that could affect an analysis and far
# above the noise floor.
TOL = 1e-3
shift = abs(max(d_out) - max(d_in))
bad = shift > TOL or (max(d_out) - min(d_out)) > TOL
if bad:
    say("  FAILED: root age moved by %.8f Ma (tolerance %.0e), or the tree is no" % (shift, TOL))
    say("  longer ultrametric: spread after pruning is %.8f." % (max(d_out) - min(d_out)))
    say("  The output has NOT been written. Send me this log.")
    sys.exit(1)
say("  root age unchanged to within %.0e Ma (%.6f), still ultrametric." % (TOL, max(d_out)))
say("  Pruning was clean: no surviving tip's age moved.")

os.makedirs(os.path.dirname(TREE_OUT) or ".", exist_ok=True)
os.makedirs(DIR_RES, exist_ok=True); os.makedirs(DIR_LOG, exist_ok=True)
open(TREE_OUT, "w").write(write_newick(new_root) + "\n")
with open(os.path.join(DIR_RES, "phylo03_collapse_audit.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, ["parent_species", "subspecies_tip", "action", "rule"])
    w.writeheader(); w.writerows(audit)
open(os.path.join(DIR_LOG, "phylo03_collapse_subspecies.log"), "w").write("\n".join(_out) + "\n")

say()
say("  wrote %s" % TREE_OUT)
say("  wrote results/phylo03_collapse_audit.csv   (every kept/pruned decision)")
say()
say("  NEXT: re-run phylo02 on the collapsed tree. It must report 0 unexpected")
say("  tips, the same root age, and EXIT=0.")
