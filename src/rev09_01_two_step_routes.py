#!/usr/bin/env python3
"""
rev09_01_two_step_routes.py

STEP 9.  ARE THE BREACHES REALLY BREACHES?

Five curated events use a resistance allele that is NOT one nucleotide
substitution from the species' own codon. They are the paper's exceptions and
they set the breach rate at 4.63%.

This asks whether they are exceptions at all, or second steps.

  THE IDEA.  A two-step substitution is not one improbable event; it is two
  ordinary ones. If the intermediate residue is ITSELF a documented resistance
  allele at that same site, then the population had a reason to be sitting on
  the intermediate codon, and the "two-step" move is two sequential one-step
  moves, each accessible at the moment it was taken. Accessibility would then
  be operating at every step rather than being violated.

  Worked by hand before writing this, for the three breach substitutions:

    VGSC T929  T->V   ACT(T) -> ATT(Ile = T929I) -> GTT(V)
    VGSC L1014 L->C   TTA(L) -> TTT(Phe = L1014F) -> TGT(C); also via TGG (L1014W)
    RDL  A301  A->N   no route through any documented allele; only Asp or Thr

  So four of the five events (T929V twice, L1014C twice) have such a route and
  one does not. L1014F is the canonical kdr mutation, which makes that route
  more than arithmetic.

  WHY THIS NEEDS A NULL.  The pattern was found by looking at the five
  breaches, which is exactly how spurious patterns are found. The question is
  not "do routes exist" - with 61 sense codons, routes usually exist - but
  whether OBSERVED inaccessible alleles have such routes MORE OFTEN than
  inaccessible alleles that were never observed. That comparison is what this
  script computes, over every species x site x allele combination in the atlas
  where the species carries the wild-type residue and the allele is out of
  one-step reach.

  TWO STRENGTHS OF CLAIM, both reported:
    BROAD   the intermediate residue is documented at that site in ANY species
    STRICT  the intermediate residue is documented at that site in THAT SAME
            species, so the intermediate demonstrably segregated in the
            population that made the second step. Strict is the falsifiable
            version and the one worth putting in a figure.

  WHAT WOULD REFUTE IT.  If observed breaches have documented-intermediate
  routes no more often than unobserved inaccessible alleles, the pattern is an
  artefact of the genetic code and the section should say so.

Reads   data/processed/codon_atlas.csv
        results/core_test_events.csv
Writes  results/rev09_two_step_routes.csv
        results/rev09_breach_detail.csv
        results/rev09_summary.md

stdlib + the project's codon_engine.
"""
from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PROJ = os.path.dirname(HERE)
RES = os.path.join(PROJ, "results")

from codon_engine import normalise, translate_codon               # noqa: E402

ATLAS = os.path.join(PROJ, "data", "processed", "codon_atlas.csv")


def resolve_results_file(stem, ext=".csv"):
    """
    Find a results file that may carry an analysis tag, and SAY WHICH ONE.

    step4_01_clean_results_tree.py moves superseded outputs to
    results/superseded/, and the project tags its canonical outputs "_primary".
    The old convention was

        for tag in ("_primary", "_t1", "_t5", "_final", ""):

    with a trailing "" that let a missing _primary fall through to an untagged
    file in silence - the C4 finding in step5_01_corrections.py. This resolver
    keeps the search order but refuses to be quiet about it: it prints the file
    it chose, and if the choice was not the preferred tag it says so loudly,
    because reading a superseded table without noticing is exactly how a
    corrected number gets un-corrected.
    """
    prefer = ("_primary", "_t1", "_t5", "_final", "")
    places = (RES, os.path.join(RES, "superseded"))
    found = []
    for d in places:
        for tag in prefer:
            cand = os.path.join(d, stem + tag + ext)
            if os.path.exists(cand):
                found.append((d == RES, prefer.index(tag), cand))
    if not found:
        sys.exit(
            "FATAL: no file matching %s*%s under results/ or results/superseded/.\n"
            "Run the core test first, or say which file to use." % (stem, ext)
        )
    found.sort(key=lambda t: (not t[0], t[1]))
    chosen = found[0][2]
    print("  input: %s" % os.path.relpath(chosen, PROJ))
    if os.path.basename(chosen) != stem + "_primary" + ext:
        print("  NOTE the preferred '%s_primary%s' was not found; using the above."
              % (stem, ext))
    if "superseded" in chosen:
        print("  WARNING this file lives in results/superseded/. It was classified as")
        print("          SUPERSEDED by step4_01_clean_results_tree.py. Check that it is")
        print("          really the table you want before trusting this output.")
    if len(found) > 1:
        print("  (other candidates present: %s)"
              % ", ".join(os.path.relpath(c, PROJ) for _, _, c in found[1:]))
    return chosen

BASES = "ACGT"
SENSE = [c for c in ("".join(p) for p in product(BASES, repeat=3))
         if translate_codon(c) != "*"]
AA_OF = {c: translate_codon(c) for c in SENSE}
CODONS_OF = defaultdict(list)
for c in SENSE:
    CODONS_OF[AA_OF[c]].append(c)


def one_step(a, b):
    return sum(1 for i in range(3) if a[i] != b[i]) == 1


NEIGHBOURS = {c: [d for d in SENSE if one_step(c, d)] for c in SENSE}


def routes(codon, target_aa):
    """Every codon -> mid -> target_codon path of exactly two substitutions."""
    out = []
    for mid in NEIGHBOURS.get(codon, ()):
        for tgt in CODONS_OF.get(target_aa, ()):
            if one_step(mid, tgt):
                out.append((mid, AA_OF[mid], tgt))
    return out


def fisher_exact(a, b, c, d):
    """Two-sided Fisher exact on [[a,b],[c,d]] by summing tables no more likely."""
    n = a + b + c + d
    r1, c1 = a + b, a + c

    def p_of(x):
        y, z, w = r1 - x, c1 - x, n - r1 - c1 + x
        if min(y, z, w) < 0:
            return 0.0
        return (math.comb(r1, x) * math.comb(n - r1, c1 - x)) / math.comb(n, c1)

    obs = p_of(a)
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    return min(1.0, sum(p_of(x) for x in range(lo, hi + 1) if p_of(x) <= obs * (1 + 1e-9)))


def main():
    if not os.path.exists(ATLAS):
        sys.exit("missing input: %s" % ATLAS)
    EVENTS = resolve_results_file("core_test_events")

    # ---------------------------------------------------- documented alleles
    observed = set()          # (species, gene, site, allele) actually used
    site_alleles = defaultdict(set)      # (gene, site) -> alleles documented anywhere
    sp_site_alleles = defaultdict(set)   # (species, gene, site) -> alleles in THAT species
    with open(EVENTS) as fh:
        for r in csv.DictReader(fh):
            key = (r["species"], r["gene"], r["site_label"], r["resistance_allele"])
            observed.add(key)
            site_alleles[(r["gene"], r["site_label"])].add(r["resistance_allele"])
            sp_site_alleles[(r["species"], r["gene"], r["site_label"])].add(r["resistance_allele"])

    print("=" * 78)
    print("STEP 9.  Two-step routes: are the breaches second steps?")
    print("=" * 78)
    print()
    print("  observed events                     : %d" % len(observed))
    print("  gene x site combinations documented : %d" % len(site_alleles))
    print()

    # ---------------------------------------------------- scan the atlas
    rows = []
    with open(ATLAS) as fh:
        for r in csv.DictReader(fh):
            if r["is_wildtype"] != "True":
                continue
            if r["already_resistant"] == "True":
                continue
            if r["one_step_accessible"] == "True":
                continue
            codon = normalise(r["codon"])
            if codon not in AA_OF:
                continue
            gene, site = r["gene"], r["site_label"]
            allele = r["resistance_allele"]
            sp = r["species"]
            rr = routes(codon, allele)
            mids = {m_aa for _, m_aa, _ in rr}
            doc_broad = sorted(mids & site_alleles.get((gene, site), set()))
            doc_strict = sorted(mids & sp_site_alleles.get((sp, gene, site), set()))
            syn = sorted(m for m in mids if m == r["wild_type_aa"])
            rows.append({
                "species": sp, "gene": gene, "site_label": site,
                "wild_type_aa": r["wild_type_aa"], "codon": codon,
                "resistance_allele": allele,
                "n_2step_routes": len(rr),
                "intermediates": "|".join(sorted(mids)),
                "route_via_documented_any": bool(doc_broad),
                "documented_any_mids": "|".join(doc_broad),
                "route_via_documented_same_species": bool(doc_strict),
                "documented_same_species_mids": "|".join(doc_strict),
                "route_via_synonymous": bool(syn),
                "observed": (sp, gene, site, allele) in observed,
            })

    n_inacc = len(rows)
    obs_rows = [r for r in rows if r["observed"]]
    unobs = [r for r in rows if not r["observed"]]
    print("-" * 78)
    print("THE COMPARISON SET")
    print("-" * 78)
    print()
    print("  species x site x allele combinations where the species carries the")
    print("  wild-type residue and the allele is NOT one substitution away: %d" % n_inacc)
    print("    of these, actually observed in the record (the breaches) : %d" % len(obs_rows))
    print("    never observed                                           : %d" % len(unobs))

    if not obs_rows:
        print()
        print("  No inaccessible combination is observed. Nothing to test.")
        return

    # ---------------------------------------------------- the tests
    def two_by_two(flag, label):
        a = sum(1 for r in obs_rows if r[flag])
        b = len(obs_rows) - a
        c = sum(1 for r in unobs if r[flag])
        d = len(unobs) - c
        p = fisher_exact(a, b, c, d)
        p_obs = a / (a + b) if a + b else float("nan")
        p_un = c / (c + d) if c + d else float("nan")
        or_ = float("inf") if (b == 0 or c == 0) else (a * d) / (b * c)
        print()
        print("  %s" % label)
        print("    observed breaches with such a route   : %d/%d = %.1f%%" % (a, a + b, 100 * p_obs))
        print("    unobserved inaccessible with a route  : %d/%d = %.1f%%" % (c, c + d, 100 * p_un))
        print("    odds ratio %s    Fisher exact P = %.4g"
              % ("inf" if or_ == float("inf") else "%.2f" % or_, p))
        return {"test": label, "obs_with": a, "obs_total": a + b,
                "unobs_with": c, "unobs_total": c + d,
                "pct_obs": round(100 * p_obs, 1), "pct_unobs": round(100 * p_un, 1),
                "odds_ratio": "inf" if or_ == float("inf") else round(or_, 3),
                "fisher_p": "%.4g" % p}

    print()
    print("-" * 78)
    print("DOES AN OBSERVED BREACH HAVE A ROUTE THROUGH A DOCUMENTED ALLELE")
    print("MORE OFTEN THAN AN INACCESSIBLE ALLELE THAT WAS NEVER USED?")
    print("-" * 78)
    tests = [
        two_by_two("route_via_documented_any",
                   "BROAD  intermediate documented at this site in ANY species"),
        two_by_two("route_via_documented_same_species",
                   "STRICT intermediate documented at this site in THAT species"),
        two_by_two("route_via_synonymous",
                   "CONTROL intermediate is synonymous (same residue)"),
    ]

    print()
    print("  Read the STRICT row as the real result. BROAD can be inflated because")
    print("  well-studied sites accumulate many documented alleles, which makes a")
    print("  documented intermediate easier to find for reasons unrelated to the")
    print("  species in question.")

    # ---------------------------------------------------- the five, in detail
    print()
    print("-" * 78)
    print("THE BREACHES, ONE BY ONE")
    print("-" * 78)
    detail = []
    for r in sorted(obs_rows, key=lambda x: (x["gene"], x["site_label"], x["species"])):
        print()
        print("  %s   %s %s  %s -> %s   codon %s"
              % (r["species"], r["gene"], r["site_label"], r["wild_type_aa"],
                 r["resistance_allele"], r["codon"]))
        print("    two-step routes: %d   intermediate residues: %s"
              % (r["n_2step_routes"], r["intermediates"] or "none"))
        if r["documented_same_species_mids"]:
            print("    route through an allele documented IN THIS SPECIES: %s"
                  % r["documented_same_species_mids"])
            print("    -> consistent with a second step from a codon this population already had")
        elif r["documented_any_mids"]:
            print("    route through an allele documented at this site elsewhere: %s"
                  % r["documented_any_mids"])
            print("    -> plausible, but the intermediate is not documented in this species;")
            print("       that is the check to make in the literature before claiming it")
        else:
            print("    NO route through any documented resistance allele.")
            print("    -> this one is a genuine exception and should be reported as such")
        for mid, mid_aa, tgt in routes(r["codon"], r["resistance_allele"])[:6]:
            print("      %s(%s) -> %s(%s) -> %s(%s)"
                  % (r["codon"], r["wild_type_aa"], mid, mid_aa, tgt, r["resistance_allele"]))
        detail.append(r)

    # ---------------------------------------------------- write
    os.makedirs(RES, exist_ok=True)
    keys = list(rows[0].keys())
    with open(os.path.join(RES, "rev09_two_step_routes.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(RES, "rev09_breach_detail.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in detail:
            w.writerow(r)
    with open(os.path.join(RES, "rev09_tests.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tests[0].keys()))
        w.writeheader()
        for t in tests:
            w.writerow(t)

    print()
    print("=" * 78)
    print("WROTE")
    print("  results/rev09_two_step_routes.csv   every inaccessible combination")
    print("  results/rev09_breach_detail.csv     the observed ones")
    print("  results/rev09_tests.csv             the 2x2 tests")
    print("=" * 78)


if __name__ == "__main__":
    main()
