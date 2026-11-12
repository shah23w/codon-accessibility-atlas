#!/usr/bin/env python3
"""
aprd02_sampling_fraction.py - STEP APRD-2.  The ascertainment denominator.

THE OBJECTION THIS ANSWERS
--------------------------
Ascertainment bias has several sources here, and the evolutionary histories of
the surveyed resistance genes differ. Stated sharply: the 108 events are not a
sample of evolution, they are a sample of what somebody chose to sequence. If species whose codon happens
to sit one step from a resistance allele are ALSO the species people sequence,
the odds ratio is manufactured by the sampling, not by mutation.

Until now that has been a stated limitation. APRD makes it a measured quantity,
because APRD records PHENOTYPIC resistance - who evolved resistance in the field
- with no reference to any mutation. This project's species panel came from NCBI
Datasets (594 annotated arthropod genomes), so the two lists are independent.

THREE ANALYSES, IN ASCENDING ORDER OF WHAT THEY SETTLE
------------------------------------------------------
A.  THE SAMPLING FRACTION (descriptive).  Of the species with documented
    field-evolved resistance to a gene's mode-of-action class, how many have a
    genome in the panel, and how many have a characterised target-site genotype
    in our record? Reported per gene and per taxonomic order.

B.  THE ODDS RATIO AMONG THE PHENOTYPICALLY RESISTANT (the decisive one).
    Re-run the site-stratified Mantel-Haenszel estimate restricted to species
    that APRD already records as resistant to that gene's class. Every species
    in that set had a documented reason to be studied, so "we only see what
    someone looked at" is held roughly constant across it. If the odds ratio
    survives there, ascertainment is not manufacturing it. If it collapses, it
    was.

C.  IS SAMPLING ITSELF ACCESSIBILITY-BIASED?  Among phenotypically resistant
    species that have a genome, compare the accessible fraction of those that
    were genotyped against those that were not, with a label permutation.
    Under the ascertainment story the genotyped ones should be the accessible
    ones; if the two groups look alike, they are not.

    READ C WITH ITS CAVEAT. "Genotyped" here means "appears in our event
    record", which conflates two things: somebody looked, AND they found a
    target-site mutation. A species studied and found to carry no target-site
    change is absent for the second reason, not the first. So C bounds the
    bias rather than isolating it, and it is reported as such. B does not have
    this problem, which is why B is the analysis the paper should lead with.

NAMES
-----
APRD names are matched to the atlas by exact lowercase binomial first, then by
the documented fallbacks: a trinomial is retried as its first two tokens, and
synonym notation ("cosmophila=anomis flava") is retried under each genus.
Genus-level records ("anopheles spp") are EXCLUDED from every denominator and
counted separately - a denominator inflated by records that name no species is
a wrong denominator that still produces a number. Every unmatched name is
written out for review rather than dropped in silence.

Reads   data/processed/aprd_species_by_moa.csv     (from aprd01_parse.py)
        data/processed/codon_atlas.csv
        results/core_test_events*.csv               (via the loud resolver)
Writes  results/aprd02_sampling_fraction.csv
        results/aprd02_restricted_odds_ratio.csv
        results/aprd02_unmatched_names.csv
        results/aprd02_summary.md
"""
from __future__ import annotations

import csv
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PROJ = os.path.dirname(HERE)
RES = os.path.join(PROJ, "results")
PROC = os.path.join(PROJ, "data", "processed")

from rev09_01_two_step_routes import resolve_results_file            # noqa: E402
from rev10_01_background_control import mantel_haenszel              # noqa: E402

ATLAS = os.path.join(PROJ, "data", "processed", "codon_atlas.csv")
APRD = os.path.join(PROC, "aprd_species_by_moa.csv")
PANEL = os.path.join(PROC, "species_panel.csv")
SYNONYMS = os.path.join(PROJ, "data", "curated", "aprd_species_synonyms.csv")

PERMS = 10000
SEED = 20260928

# An APRD gene label may cover more than one atlas gene.
ATLAS_GENES_FOR = {
    "AChE": ["AChE"],
    "RDL": ["RDL"],
    "VGSC": ["VGSC"],
    "GluCl": ["GluCl_Px", "GluCl_Tu"],
    "CHS1": ["CHS1"],
    "RyR": ["RyR"],
}


def norm(s):
    return " ".join((s or "").split()).lower()


def name_variants(sp):
    """
    Documented fallbacks, tried in order. Nothing here invents a name: each
    variant is a literal re-reading of the string APRD supplied.
    """
    out = [sp]
    toks = sp.split()
    if "=" in sp:
        # "cosmophila=anomis flava" -> "cosmophila flava", "anomis flava"
        head, rest = toks[0], toks[1:]
        for alt in head.split("="):
            if alt:
                out.append(" ".join([alt] + rest))
    if len(toks) > 2:
        out.append(" ".join(toks[:2]))
        if "=" in toks[0]:
            for alt in toks[0].split("="):
                if alt:
                    out.append(" ".join([alt, toks[1]]))
    seen, uniq = set(), []
    for v in out:
        v = norm(v)
        if v and v not in seen:
            seen.add(v)
            uniq.append(v)
    return uniq


def main():
    print("=" * 78)
    print("STEP APRD-2.  The ascertainment denominator")
    print("=" * 78)
    print()

    if not os.path.exists(APRD):
        sys.exit("FATAL: %s not found. Run src/aprd01_parse.py first." % APRD)

    # ---- atlas ------------------------------------------------------------
    atlas_rows = []
    species_order, species_by_lower = {}, {}
    with open(ATLAS) as fh:
        for r in csv.DictReader(fh):
            species_by_lower.setdefault(norm(r["species"]), r["species"])
            species_order[r["species"]] = r.get("order") or "UNPLACED"
            if r["is_wildtype"] != "True" or r["already_resistant"] == "True":
                continue
            atlas_rows.append({
                "species": r["species"], "gene": r["gene"],
                "site": "%s:%s" % (r["gene"], r["site_label"]),
                "site_label": r["site_label"],
                "resistance_allele": r["resistance_allele"],
                "accessible": r["one_step_accessible"] == "True",
            })
    print("  atlas: %d wild-type rows, %d species with codon data"
          % (len(atlas_rows), len(species_by_lower)))
    # The genome panel and the codon atlas are NOT the same set. Matching is
    # done against the atlas, because a species with no codon data cannot
    # contribute accessibility; naming the column "in panel" hid a real
    # difference and made one identical name look like a spelling variant.
    if os.path.exists(PANEL):
        with open(PANEL) as fh:
            panel = {norm(r["species"]) for r in csv.DictReader(fh)}
        gap = sorted(panel - set(species_by_lower))
        print("  panel: %d species; %d in the panel have NO codon data%s"
              % (len(panel), len(gap),
                 (" (%s)" % ", ".join(gap[:4])) if gap else ""))

    # ---- events -----------------------------------------------------------
    ev_path = resolve_results_file("core_test_events")
    observed, genotyped = set(), defaultdict(set)
    with open(ev_path) as fh:
        for r in csv.DictReader(fh):
            observed.add((r["species"], r["gene"], r["site_label"],
                          r["resistance_allele"]))
            genotyped[r["gene"]].add(r["species"])
    print("  events: %d, across %d genes\n" % (len(observed), len(genotyped)))

    # ---- APRD denominator, per APRD gene label ----------------------------
    aprd_raw = defaultdict(set)
    genus_level = defaultdict(set)
    with open(APRD) as fh:
        for r in csv.DictReader(fh):
            gene = r["gene"]
            if not gene:
                continue
            flags = (r.get("name_flags") or "").split(";")
            if "genus-level-record" in flags or "genus-only" in flags:
                genus_level[gene].add(r["species"])
                continue
            aprd_raw[gene].add(r["species"])

    # ---- reviewed synonyms -------------------------------------------------
    # An explicit, auditable list beats a fuzzy match: every line records which
    # APRD string maps to which atlas string and why. See
    # data/curated/README_aprd_synonyms.md.
    syn, syn_status, syn_kind = {}, {}, {}
    if os.path.exists(SYNONYMS):
        with open(SYNONYMS) as fh:
            for r in csv.DictReader(fh):
                a, t = norm(r["aprd_name"]), (r["atlas_name"] or "").strip()
                if not a or not t:
                    continue
                if norm(t) not in species_by_lower:
                    print("  SYNONYM IGNORED: %r -> %r, but %r is not in the atlas"
                          % (a, t, t))
                    continue
                syn[a] = species_by_lower[norm(t)]
                syn_status[a] = (r.get("status") or "proposed").strip().lower()
                syn_kind[a] = (r.get("kind") or "").strip()
        n_prop = sum(1 for v in syn_status.values() if v != "confirmed")
        print("  synonyms: %d reviewed mappings loaded from data/curated/" % len(syn))
        if n_prop:
            # A taxonomic judgement that has not been confirmed must not enter a
            # denominator silently. It is applied, and it is named.
            print("  %d of them are still marked PROPOSED and need confirming:" % n_prop)
            for a in sorted(syn):
                if syn_status.get(a) != "confirmed":
                    print("      %-32s -> %-32s (%s)" % (a, syn[a], syn_kind.get(a, "")))
    print()

    # ---- match APRD names to the atlas ------------------------------------
    matched, unmatched = defaultdict(dict), []
    n_by_synonym = 0
    for gene, names in aprd_raw.items():
        for sp in sorted(names):
            hit = None
            for v in name_variants(sp):
                if v in species_by_lower:
                    hit = species_by_lower[v]
                    break
            if hit is None and sp in syn:
                hit = syn[sp]
                n_by_synonym += 1
            if hit:
                matched[gene][sp] = hit
            else:
                unmatched.append({"gene": gene, "aprd_species": sp,
                                  "variants_tried": "|".join(name_variants(sp))})

    with open(os.path.join(RES, "aprd02_unmatched_names.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "aprd_species", "variants_tried"])
        w.writeheader()
        for r in sorted(unmatched, key=lambda r: (r["gene"], r["aprd_species"])):
            w.writerow(r)

    # =====================================================================
    # A.  THE SAMPLING FRACTION
    # =====================================================================
    print("-" * 78)
    print("A.  SAMPLING FRACTION - of species with documented field-evolved")
    print("    resistance, how many have a genome, and how many are genotyped?")
    print("-" * 78)
    print("  %-8s %10s %10s %10s %12s %12s"
          % ("gene", "APRD", "in atlas", "genotyped", "frac|APRD", "frac|atlas"))
    rows_A = []
    for gene in sorted(aprd_raw, key=lambda g: -len(aprd_raw[g])):
        atlas_gs = ATLAS_GENES_FOR.get(gene, [gene])
        # Count ORGANISMS, not APRD name-strings. APRD lists some species twice
        # under an old and a current genus - heliothis zea and helicoverpa zea
        # are both present - so a denominator built from raw strings counts one
        # animal as two, inflating the denominator and understating the
        # sampling fraction. Resolving through the synonym table first is what
        # makes the merge worth doing; without this, merging would have made
        # the count worse rather than better.
        resolved = {sp: matched[gene].get(sp, sp) for sp in aprd_raw[gene]}
        denom = set(resolved.values())
        n_collapsed = len(aprd_raw[gene]) - len(denom)
        in_panel = {v for sp, v in resolved.items() if sp in matched[gene]}
        gt = set()
        for ag in atlas_gs:
            for v in in_panel:
                if v in genotyped.get(ag, ()):
                    gt.add(v)
        f_all = len(gt) / len(denom) if denom else float("nan")
        f_pan = len(gt) / len(in_panel) if in_panel else float("nan")
        print("  %-8s %10d %10d %10d %11.1f%% %11.1f%%"
              % (gene, len(denom), len(in_panel), len(gt), 100 * f_all, 100 * f_pan))
        rows_A.append({"gene": gene,
                       "aprd_classes": "|".join(atlas_gs),
                       "n_aprd_species": len(denom),
                       "n_genus_level_excluded": len(genus_level.get(gene, ())),
                       "n_aprd_strings": len(aprd_raw[gene]),
                       "n_collapsed_by_synonym": n_collapsed,
                       "n_in_codon_atlas": len(in_panel),
                       "n_genotyped": len(gt),
                       "sampling_fraction_of_aprd": round(f_all, 5),
                       "sampling_fraction_of_atlas": round(f_pan, 5)})
    print()
    n_coll = sum(r["n_collapsed_by_synonym"] for r in rows_A)
    if n_coll:
        print("  APRD entries collapsed as the same organism under two genus "
              "names: %d (the denominator counts organisms, not name-strings)"
              % n_coll)
    print("  genus-level records excluded from every denominator: %d"
          % sum(len(v) for v in genus_level.values()))
    if n_by_synonym:
        print("  matches added by the reviewed synonym table: %d" % n_by_synonym)
    print("  APRD names with no codon data and no match: %d "
          "(results/aprd02_unmatched_names.csv)" % len(unmatched))
    print()

    with open(os.path.join(RES, "aprd02_sampling_fraction.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_A[0].keys()))
        w.writeheader()
        for r in rows_A:
            w.writerow(r)

    # =====================================================================
    # B.  THE ODDS RATIO AMONG THE PHENOTYPICALLY RESISTANT
    # =====================================================================
    print("-" * 78)
    print("B.  ODDS RATIO RESTRICTED TO SPECIES APRD ALREADY CALLS RESISTANT")
    print("    Every species in this set had a documented reason to be studied.")
    print("-" * 78)

    keep = set()                      # (atlas species, atlas gene)
    for gene, m in matched.items():
        for ag in ATLAS_GENES_FOR.get(gene, [gene]):
            for sp_atlas in m.values():
                keep.add((sp_atlas, ag))

    def strata(rows):
        by = defaultdict(lambda: [0, 0, 0, 0])
        for r in rows:
            used = (r["species"], r["gene"], r["site_label"],
                    r["resistance_allele"]) in observed
            by[r["site"]][(0 if used else 2) + (0 if r["accessible"] else 1)] += 1
        return [tuple(v) for v in by.values()]

    full = strata(atlas_rows)
    rest_rows = [r for r in atlas_rows if (r["species"], r["gene"]) in keep]
    rest = strata(rest_rows)

    o_f = mantel_haenszel(full)
    o_r = mantel_haenszel(rest)
    n_ev_rest = sum(1 for r in rest_rows
                    if (r["species"], r["gene"], r["site_label"],
                        r["resistance_allele"]) in observed)
    print("  %-34s %8s %9s %9s %8s %8s"
          % ("set", "rows", "events", "OR", "lo", "hi"))
    print("  %-34s %8d %9d %9.2f %8.2f %8.2f"
          % ("all species with a genome", len(atlas_rows), len(observed),
             o_f[0], o_f[1], o_f[2]))
    print("  %-34s %8d %9d %9.2f %8.2f %8.2f"
          % ("APRD phenotypically resistant", len(rest_rows), n_ev_rest,
             o_r[0], o_r[1], o_r[2]))
    print()

    # Which genes actually survive the restriction? A restricted estimate
    # carried almost entirely by one gene is a different claim from one spread
    # across seven, so it is shown rather than left implicit.
    print("  what the restriction keeps, per gene")
    print("  %-10s %9s %9s %9s %9s"
          % ("gene", "species", "kept", "events", "kept"))
    gene_rows = defaultdict(lambda: [set(), set(), 0, 0])
    for r in atlas_rows:
        g = gene_rows[r["gene"]]
        g[0].add(r["species"])
        used = (r["species"], r["gene"], r["site_label"],
                r["resistance_allele"]) in observed
        if used:
            g[2] += 1
        if (r["species"], r["gene"]) in keep:
            g[1].add(r["species"])
            if used:
                g[3] += 1
    for g in sorted(gene_rows, key=lambda g: -len(gene_rows[g][0])):
        allsp, kept, ev, evk = gene_rows[g]
        print("  %-10s %9d %9d %9d %9d" % (g, len(allsp), len(kept), ev, evk))
    print()

    with open(os.path.join(RES, "aprd02_restricted_odds_ratio.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["set", "n_rows", "n_events", "n_informative_strata",
                    "odds_ratio", "ci_lo", "ci_hi"])
        w.writerow(["all species with a genome", len(atlas_rows), len(observed),
                    o_f[3], round(o_f[0], 4), round(o_f[1], 4), round(o_f[2], 4)])
        w.writerow(["APRD phenotypically resistant", len(rest_rows), n_ev_rest,
                    o_r[3], round(o_r[0], 4), round(o_r[1], 4), round(o_r[2], 4)])

    # =====================================================================
    # C.  IS SAMPLING ITSELF ACCESSIBILITY-BIASED?
    # =====================================================================
    print("-" * 78)
    print("C.  AMONG PHENOTYPICALLY RESISTANT SPECIES WITH A GENOME, ARE THE")
    print("    GENOTYPED ONES THE ACCESSIBLE ONES?")
    print("-" * 78)

    frac = defaultdict(dict)
    tally = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in atlas_rows:
        t = tally[r["gene"]][r["species"]]
        t[0] += 1 if r["accessible"] else 0
        t[1] += 1
    for g, d in tally.items():
        for sp, (a, n) in d.items():
            frac[g][sp] = a / n if n else float("nan")

    rng = random.Random(SEED)
    rows_C = []
    print("  %-10s %7s %7s %11s %11s %10s %9s"
          % ("gene", "n gt", "n not", "acc|gt", "acc|not", "diff", "P"))
    for gene in sorted(aprd_raw, key=lambda g: -len(aprd_raw[g])):
        for ag in ATLAS_GENES_FOR.get(gene, [gene]):
            pool = sorted({matched[gene][sp] for sp in matched[gene]}
                          & set(frac.get(ag, {})))
            if len(pool) < 8:
                continue
            gt = set(genotyped.get(ag, ()))
            x = [frac[ag][sp] for sp in pool if sp in gt]
            y = [frac[ag][sp] for sp in pool if sp not in gt]
            if len(x) < 3 or len(y) < 3:
                continue
            obs = sum(x) / len(x) - sum(y) / len(y)
            vals = [frac[ag][sp] for sp in pool]
            k, hits = len(x), 0
            for _ in range(PERMS):
                rng.shuffle(vals)
                d = sum(vals[:k]) / k - sum(vals[k:]) / (len(vals) - k)
                if abs(d) >= abs(obs) - 1e-12:
                    hits += 1
            p = (hits + 1) / (PERMS + 1)
            print("  %-10s %7d %7d %10.3f %10.3f %+10.3f %9.4f"
                  % (ag, len(x), len(y), sum(x) / len(x), sum(y) / len(y), obs, p))
            rows_C.append({"atlas_gene": ag, "aprd_gene": gene,
                           "n_genotyped": len(x), "n_not_genotyped": len(y),
                           "mean_accessible_genotyped": round(sum(x) / len(x), 5),
                           "mean_accessible_not_genotyped": round(sum(y) / len(y), 5),
                           "difference": round(obs, 5), "p_permutation": round(p, 5),
                           "n_permutations": PERMS})
    if rows_C:
        print("  %d genes tested separately; no multiplicity adjustment applied."
              % len(rows_C))
        with open(os.path.join(RES, "aprd02_accessibility_of_sampling.csv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows_C[0].keys()))
            w.writeheader()
            for r in rows_C:
                w.writerow(r)
    else:
        print("  (no gene had enough species on both sides to test)")
    print()

    # ---- summary ----------------------------------------------------------
    md = ["# APRD-2. The ascertainment denominator", "",
          "APRD records phenotypic resistance with no reference to any mutation; "
          "the species panel came from NCBI Datasets. The two lists are "
          "independent, so the ratio below is not circular.", "",
          "## A. Sampling fraction", "",
          "| gene | APRD species | with codon data | genotyped | % of APRD | % with data |",
          "|---|---|---|---|---|---|"]
    for r in rows_A:
        md.append("| %s | %d | %d | %d | %.1f%% | %.1f%% |"
                  % (r["gene"], r["n_aprd_species"], r["n_in_codon_atlas"],
                     r["n_genotyped"], 100 * r["sampling_fraction_of_aprd"],
                     100 * r["sampling_fraction_of_atlas"]))
    md += ["", "Genus-level records (e.g. \"anopheles spp\") excluded from every "
           "denominator: %d. APRD names with no match in the atlas: %d, listed in "
           "`results/aprd02_unmatched_names.csv`."
           % (sum(len(v) for v in genus_level.values()), len(unmatched)),
           "", "## B. Odds ratio among the phenotypically resistant", "",
           "| set | rows | events | OR [95% CI] |", "|---|---|---|---|",
           "| all species with a genome | %d | %d | %.2f [%.2f, %.2f] |"
           % (len(atlas_rows), len(observed), o_f[0], o_f[1], o_f[2]),
           "| APRD phenotypically resistant | %d | %d | %.2f [%.2f, %.2f] |"
           % (len(rest_rows), n_ev_rest, o_r[0], o_r[1], o_r[2]), "",
           "Every species in the restricted set is one APRD already records as "
           "having evolved resistance in the field, so each had a documented "
           "reason to be studied. An odds ratio that holds here is not an "
           "artefact of who gets sequenced.", "",
           "| gene | species in panel | kept by the restriction | events | events kept |",
           "|---|---|---|---|---|"]
    for g in sorted(gene_rows, key=lambda g: -len(gene_rows[g][0])):
        allsp, kept, ev, evk = gene_rows[g]
        md.append("| %s | %d | %d | %d | %d |" % (g, len(allsp), len(kept), ev, evk))
    md += [""]
    if rows_C:
        md += ["## C. Is sampling itself accessibility-biased?", "",
               "| gene | n genotyped | n not | accessible, genotyped | accessible, not | difference | P |",
               "|---|---|---|---|---|---|---|"]
        for r in rows_C:
            md.append("| %s | %d | %d | %.3f | %.3f | %+.3f | %.4f |"
                      % (r["atlas_gene"], r["n_genotyped"], r["n_not_genotyped"],
                         r["mean_accessible_genotyped"],
                         r["mean_accessible_not_genotyped"],
                         r["difference"], r["p_permutation"]))
        md += ["", "%d genes were tested separately and no multiplicity "
               "adjustment is applied; read the individual P values with that in "
               "mind." % len(rows_C),
               "", "**Caveat.** \"Genotyped\" means \"appears in the event record\", "
               "which conflates *somebody looked* with *they found a target-site "
               "mutation*. A species studied and found to carry no target-site "
               "change is absent for the second reason, not the first. C therefore "
               "bounds the bias rather than isolating it; B does not have this "
               "problem, which is why B is the analysis to lead with.", ""]
    with open(os.path.join(RES, "aprd02_summary.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")

    print("  wrote results/aprd02_sampling_fraction.csv")
    print("  wrote results/aprd02_restricted_odds_ratio.csv")
    if rows_C:
        print("  wrote results/aprd02_accessibility_of_sampling.csv")
    print("  wrote results/aprd02_unmatched_names.csv")
    print("  wrote results/aprd02_summary.md")
    print("\nDone.")


if __name__ == "__main__":
    main()
