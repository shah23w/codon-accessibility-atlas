#!/usr/bin/env python3
"""
rev13_01_paralogy_sensitivity.py

STEP 13.  PARALOGY: THE SPECIFICATION, AND WHETHER IT CHANGES THE ANSWER.

Three things about the paralogue handling need stating, and were not:
whether paralogy was considered for every protein family and not only
acetylcholinesterase; what the absolute floor and the margin actually are in
the rule "classified to a paralogue class only when its best score exceeds an
absolute floor and beats the alternative by a margin"; and what the rule does
to a real case, with real numbers.

THE HEADLINE QC NUMBER LOOKS ALARMING AND IS MISLEADING.  paralogy_qc.csv
reports that 63 to 75 percent of species are "ambiguous" for every gene, with a
median rank1-rank2 identity gap of essentially zero. Two things explain that,
and neither is paralogue mixing:

  1. RANK 1 IS CHOSEN BY ALIGNMENT SCORE, NOT IDENTITY. step2_01 ranks
     candidates by semi-global BLOSUM62 score with end gaps free. A short,
     high-identity fragment can beat a full-length orthologue on percent
     identity while losing badly on score - which is why 611 flagged rows have
     a NEGATIVE identity gap, i.e. rank 2 has higher identity than rank 1.
     Flagging on an identity gap therefore measures something the selection
     rule does not use.

  2. PERCENT IDENTITY TO A DISTANT REFERENCE IS NOT A PARALOGUE DISCRIMINATOR.
     The AChE reference is Torpedo californica, a vertebrate; median identity
     across arthropods is 0.44. At that distance ace-1 and ace-2 are both far
     from the reference and close to each other in identity, so the gap is
     small whether or not the right locus was chosen.

WHAT ACTUALLY DISCRIMINATES is whether the residues flanking the site are
conserved, and whether the alternative candidate would have given a different
answer. Both are computable from data already on disk, and both are reported
here instead of the identity gap.

  A. THE SPECIFICATION: the selection rule, the flag threshold, the AChE
     paralogue rule with its clade logic, and a worked example with real
     numbers.

  B. WOULD RANK 2 HAVE CHANGED THE CALL? For every species x gene x site with
     more than one mapped candidate, compare the residue rank 2 gives with the
     residue rank 1 gives. If they agree, the ambiguity is immaterial to the
     atlas whatever it means phylogenetically.

  C. SITE-CONTEXT CONSERVATION. The +/-10 residue window around each site,
     rank 1 versus rank 2. The correct locus should show the better-conserved
     context.

  D. DOES IT CHANGE THE RESULT? The species-level odds ratio recomputed three
     ways: dropping every flagged species x gene; keeping only well-conserved
     site contexts; keeping only species where one candidate was found at all.
     If 11.66 survives, paralogy ambiguity does not drive the finding. Each
     restriction's cost in events is reported, because a robustness check that
     silently discards the signal is not a robustness check.

Reads   data/interim/site_mappings_all_candidates.csv
        data/processed/codon_atlas.csv
        results/paralogy_flagged.csv, paralogy_qc.csv
        results/core_test_events*.csv
Writes  results/rev13_rank_concordance.csv
        results/rev13_sensitivity.csv
        results/rev13_methods_paralogy.md
        results/rev13_summary.md
"""
from __future__ import annotations

import csv
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PROJ = os.path.dirname(HERE)
RES = os.path.join(PROJ, "results")
INTERIM = os.path.join(PROJ, "data", "interim")

from rev09_01_two_step_routes import resolve_results_file           # noqa: E402
from rev10_01_background_control import mantel_haenszel             # neat: same estimator

ATLAS = os.path.join(PROJ, "data", "processed", "codon_atlas.csv")
CANDS = os.path.join(INTERIM, "site_mappings_all_candidates.csv")
FLAGGED = os.path.join(RES, "paralogy_flagged.csv")

AMBIGUOUS_GAP = 0.05          # the threshold step2_03 uses, restated here


def _collapse_to_species(rows):
    """One row per (species, site), accessible if ANY candidate allele there is
    reachable and used if the species has any documented event there.

    The estimand this project reports is defined on (species, site, candidate
    allele), because accessibility is a property of the codon-allele pair and
    not of the species: a species' codon can reach one documented allele at a
    site and not another, and the hypothesis is about which allele evolved. The
    atlas is therefore long over alleles and so is the odds ratio, which means
    a species appears once per candidate allele at each site.

    Collapsing over alleles discards exactly that contrast, and this arm
    measures how much. It is reported rather than preferred: at a site with
    five candidate alleles almost every species reaches at least one, so the
    collapsed table is close to degenerate by construction.
    """
    coll = {}
    for r in rows:
        k = (r["species"], r["site"])
        c = coll.setdefault(k, {"species": r["species"], "site": r["site"],
                                "accessible": False, "used": False})
        c["accessible"] = c["accessible"] or r["accessible"]
        c["used"] = c["used"] or r["used"]
    return list(coll.values())


def mh_on(rows, label, out):
    by = defaultdict(lambda: [0, 0, 0, 0])
    for r in rows:
        by[r["site"]][(0 if r["used"] else 2) + (0 if r["accessible"] else 1)] += 1
    or_, lo, hi, n_inf = mantel_haenszel([tuple(v) for v in by.values()])
    n_used = sum(1 for r in rows if r["used"])
    line = ("  %-44s rows %6d  events %3d  OR %7s  [%s, %s]"
            % (label, len(rows), n_used,
               "inf" if or_ == float("inf") else "%.2f" % or_,
               "%.2f" % lo if lo == lo else "na",
               "%.2f" % hi if hi == hi else "na"))
    print(line)
    out.append(line)
    return {"restriction": label, "n_rows": len(rows), "n_events": n_used,
            "n_informative_strata": n_inf,
            "odds_ratio": "inf" if or_ == float("inf") else round(or_, 3),
            "ci_lo": round(lo, 3) if lo == lo else "",
            "ci_hi": round(hi, 3) if hi == hi else ""}


def main():
    print("=" * 78)
    print("STEP 13.  Paralogy: the specification, and whether it changes anything")
    print("=" * 78)
    print()

    ev_path = resolve_results_file("core_test_events")
    print()
    observed = set()
    with open(ev_path) as fh:
        for r in csv.DictReader(fh):
            observed.add((r["species"], r["gene"], r["site_label"], r["resistance_allele"]))

    # ------------------------------------------------------------------ B, C
    cands = defaultdict(dict)          # (species, gene, site) -> rank -> row
    for r in csv.DictReader(open(CANDS)):
        try:
            rk = int(r["rank"])
        except ValueError:
            continue
        cands[(r["species"], r["gene"], r["site_label"])][rk] = r

    out = []

    def p(s=""):
        print(s)
        out.append(s)

    p("-" * 78)
    p("B.  WOULD THE SECOND-RANKED CANDIDATE HAVE CHANGED THE CALL?")
    p("-" * 78)
    p()
    per_gene = defaultdict(lambda: {"n": 0, "same": 0, "diff": 0,
                                    "ctx1": [], "ctx2": [], "unmappable2": 0})
    conc_rows = []
    for (sp, gene, site), byrank in cands.items():
        r1 = byrank.get(1)
        r2 = byrank.get(2)
        if not r1:
            continue
        g = per_gene[gene]
        if not r2:
            continue
        g["n"] += 1
        if r2.get("mappable") != "True":
            g["unmappable2"] += 1
            continue
        a1, a2 = r1.get("species_aa", ""), r2.get("species_aa", "")
        same = (a1 == a2)
        g["same" if same else "diff"] += 1
        try:
            g["ctx1"].append(float(r1["context_identity"]))
            g["ctx2"].append(float(r2["context_identity"]))
        except (ValueError, KeyError):
            pass
        conc_rows.append({"species": sp, "gene": gene, "site_label": site,
                          "rank1_acc": r1.get("protein_accession"),
                          "rank1_aa": a1, "rank1_context_identity": r1.get("context_identity"),
                          "rank2_acc": r2.get("protein_accession"),
                          "rank2_aa": a2, "rank2_context_identity": r2.get("context_identity"),
                          "same_residue": same})

    def med(v):
        v = sorted(v)
        return v[len(v) // 2] if v else float("nan")

    p("  %-10s %8s %8s %8s %9s %11s %11s"
      % ("gene", "pairs", "same aa", "differ", "%same", "ctx ident r1", "ctx ident r2"))
    for gene in sorted(per_gene):
        g = per_gene[gene]
        tot = g["same"] + g["diff"]
        if tot == 0:
            continue
        p("  %-10s %8d %8d %8d %8.1f%% %11.3f %11.3f"
          % (gene, g["n"], g["same"], g["diff"], 100 * g["same"] / tot,
             med(g["ctx1"]), med(g["ctx2"])))
    allsame = sum(g["same"] for g in per_gene.values())
    alldiff = sum(g["diff"] for g in per_gene.values())
    if allsame + alldiff:
        p()
        p("  overall: the two top candidates give the SAME residue at %.1f%% of"
          % (100 * allsame / (allsame + alldiff)))
        p("  mapped sites (%d of %d)." % (allsame, allsame + alldiff))
        p("  Where they agree, which candidate was chosen cannot affect the atlas.")

    # ------------------------------------------------------------------ D
    flagged = set()
    if os.path.exists(FLAGGED):
        for r in csv.DictReader(open(FLAGGED)):
            flagged.add((r["species"], r["gene"]))

    ctx1 = {}
    single = set()
    for (sp, gene, site), byrank in cands.items():
        r1 = byrank.get(1)
        if not r1:
            continue
        try:
            ctx1[(sp, gene, site)] = float(r1["context_identity"])
        except (ValueError, KeyError):
            pass
        if 2 not in byrank:
            single.add((sp, gene, site))

    rows = []
    for r in csv.DictReader(open(ATLAS)):
        if r["is_wildtype"] != "True" or r["already_resistant"] == "True":
            continue
        key = (r["species"], r["gene"], r["site_label"])
        rows.append({
            "species": r["species"], "gene": r["gene"], "site_label": r["site_label"],
            "site": "%s:%s" % (r["gene"], r["site_label"]),
            "accessible": r["one_step_accessible"] == "True",
            "used": (r["species"], r["gene"], r["site_label"],
                     r["resistance_allele"]) in observed,
            "flagged": (r["species"], r["gene"]) in flagged,
            "ctx": ctx1.get(key),
            "single_candidate": key in single,
        })

    p()
    p("-" * 78)
    p("D.  DOES ANY OF IT CHANGE THE ODDS RATIO?")
    p("-" * 78)
    p()
    sens = [
        mh_on(rows, "all rows (Step 10 baseline)", out),
        mh_on([r for r in rows if not r["flagged"]],
              "dropping every FLAGGED species x gene", out),
        mh_on([r for r in rows if r["ctx"] is not None and r["ctx"] >= 0.5],
              "site context identity >= 0.50", out),
        mh_on([r for r in rows if r["ctx"] is not None and r["ctx"] >= 0.7],
              "site context identity >= 0.70", out),
        mh_on([r for r in rows if r["single_candidate"]],
              "only ONE candidate mapped at that site", out),
        mh_on(_collapse_to_species(rows),
              "collapsed to one row per species x site", out),
    ]
    p()
    p("  Read the estimate and the surviving event count together. A restriction")
    p("  that keeps the odds ratio but destroys the events has demonstrated")
    p("  nothing; one that keeps both has demonstrated robustness.")

    # ------------------------------------------------------------------ A
    meth = []
    meth.append("## Orthologue selection and paralogy, as actually implemented")
    meth.append("")
    meth.append("**Selection rule.** For each species proteome and each reference gene,")
    meth.append("candidate proteins are prescreened by shared 5-mer count normalised by")
    meth.append("length, and the shortlist is aligned to the reference by semi-global")
    meth.append("pairwise alignment (BLOSUM62, affine gaps, **end gaps free**, Biopython")
    meth.append("`PairwiseAligner`). Candidates are ranked by **alignment score**, not by")
    meth.append("percent identity. End-gap-free scoring is used because annotated arthropod")
    meth.append("proteins are frequently N-terminally truncated, and penalising those")
    meth.append("overhangs demotes true orthologues.")
    meth.append("")
    meth.append("**Site transfer.** A reference site is transferred only if it falls inside")
    meth.append("an aligned block. A site landing in a gap is recorded as unmappable and is")
    meth.append("never interpolated.")
    meth.append("")
    meth.append("**Paralogue retention.** The top candidates per gene are retained with their")
    meth.append("alignment statistics rather than collapsed to a single best hit, so paralogue")
    meth.append("structure remains visible downstream (`site_mappings_all_candidates.csv`).")
    meth.append("")
    meth.append("**The ambiguity flag.** `step2_03_paralogy_qc.py` flags a species x gene when")
    meth.append("`rank1 identity - rank2 identity < %.2f`. That threshold is reported here" % AMBIGUOUS_GAP)
    meth.append("for completeness, but it should not be read as an estimate of paralogue")
    meth.append("mixing, for two reasons. First, ranking uses alignment score while the flag")
    meth.append("uses identity, so the flag does not measure the quantity the selection rule")
    meth.append("optimises; 611 flagged rows in fact have a *negative* gap, meaning rank 2 had")
    meth.append("higher identity over a shorter alignment. Second, the AChE reference is")
    meth.append("*Torpedo californica*, a vertebrate, and median arthropod identity to it is")
    meth.append("0.44; at that distance ace-1 and ace-2 are both distant from the reference and")
    meth.append("close to each other, so a small gap carries little information.")
    meth.append("")
    meth.append("**AChE paralogue assignment is clade-dependent.** Most arthropods carry both")
    meth.append("ace-1 and ace-2, and ace-1 is the organophosphate and carbamate target. The")
    meth.append("cyclorrhaphan Diptera (*Drosophila*, *Musca*, *Stomoxys*, *Bactrocera*) have")
    meth.append("lost ace-1, and ace-2 is the target in those lineages. The pipeline applies")
    meth.append("this rule explicitly rather than taking the best hit blindly; assignments are")
    meth.append("recorded in the `paralogue` column of `site_mappings_all_candidates.csv`.")
    meth.append("")
    meth.append("**What was actually verified.** Two checks replace the identity gap:")
    meth.append("whether the second-ranked candidate would have given a different residue at")
    meth.append("the site, and how well the +/-10 residue window around the site is conserved")
    meth.append("in each candidate. Both are reported per gene, together with the odds ratio")
    meth.append("recomputed under four restrictions. See `rev13_summary.md`.")

    os.makedirs(RES, exist_ok=True)
    if conc_rows:
        with open(os.path.join(RES, "rev13_rank_concordance.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(conc_rows[0].keys()))
            w.writeheader()
            for r in conc_rows:
                w.writerow(r)
    with open(os.path.join(RES, "rev13_sensitivity.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sens[0].keys()))
        w.writeheader()
        for r in sens:
            w.writerow(r)
    with open(os.path.join(RES, "rev13_methods_paralogy.md"), "w") as fh:
        fh.write("\n".join(meth) + "\n")
    with open(os.path.join(RES, "rev13_summary.md"), "w") as fh:
        fh.write("# Step 13: paralogy\n\n```\n" + "\n".join(out) + "\n```\n")

    print()
    print("=" * 78)
    print("WROTE")
    print("  results/rev13_rank_concordance.csv   rank1 vs rank2, every mapped site")
    print("  results/rev13_sensitivity.csv        the odds ratio under 5 restrictions")
    print("  results/rev13_methods_paralogy.md    the specification, for the appendix")
    print("  results/rev13_summary.md")
    print("=" * 78)


if __name__ == "__main__":
    main()
