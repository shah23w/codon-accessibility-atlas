#!/usr/bin/env python3
"""
step9_04_codon_selection.py: has insecticide selection shaped the SYNONYMOUS codon at
target sites?

Why this question is worth asking
---------------------------------
The codon-constraint hypothesis itself is not new: Weill et al. (2004, Curr Biol 14:R552)
established it for AChE1 G119 across 44 mosquito species, and this project reproduces
their Table 1 at 13/13. What nobody has asked at scale is the inverse question.

Weill noticed, and did not pursue, that the bias "does not reflect global codon usage, as
GGR is used in 64.3% of Gly codons in An. gambiae". That observation implies a testable
claim: if which codon sits at a target site determines whether resistance can evolve, then
the codon at that site is itself a trait exposed to selection.

  * enriched for ACCESSIBLE codons in species under insecticide pressure -> the standing
    codon has been shaped in a way that permits resistance (selection for evolvability, or
    hitchhiking with resistance sweeps)
  * enriched for INACCESSIBLE codons -> selection against evolvability
  * no difference -> the codon at a target site is just the codon that happened to be there

The design
----------
Statistic: for a species' target site with wild-type residue X and codon C_site, compare
the accessibility of C_site to the accessibility of EVERY OTHER codon encoding X in the
SAME PROTEIN of the SAME SPECIES.

That is a within-protein null. It holds constant, by construction:
  - the species' genome-wide codon usage
  - the gene
  - the amino acid
  - GC content, expression level, and any other whole-protein property

None of the earlier nulls in this project do that. The atlas-pool null compares across
species; the within-record null compares across the studied subset. This one asks whether
the codon AT the site is unusual for that amino acid IN THAT PROTEIN.

THE ASCERTAINMENT PROBLEM, AND THE CONTRAST THAT SOLVES IT
----------------------------------------------------------
Resistance sites are known BECAUSE resistance evolved there, so sites whose codon happened
to be accessible are over-represented in the literature. A raw excess of accessibility at
target sites is therefore expected under pure ascertainment and proves nothing.

The discriminating comparison is between species, at the SAME sites:

    species with documented resistance   (decades of insecticide selection)
    versus
    species with no documented resistance at that site

Both groups are scored at the same sites with the same statistic. Ascertainment acts on
the SITE; selection acts on the SPECIES. If the codon has been shaped by insecticide
selection, the two groups must differ. If ascertainment alone is at work, they will not.

Outputs
-------
results/codon_selection_per_site.csv
results/codon_selection_summary.csv
logs/step9_04_codon_selection.log
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import sys
from collections import defaultdict

SEED = 20260731
N_BOOT = 2000

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import (  # noqa: E402
    DIR_INTERIM, DIR_PROC, DIR_REFS, DIR_RESULTS, Tee, ensure_dirs, parse_fasta,
)
from codon_engine import CODON_TABLE, codon_reaches, normalise  # noqa: E402


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def poisson_binomial_sf(ps, obs):
    d = [1.0]
    for p in ps:
        n = [0.0] * (len(d) + 1)
        for i, v in enumerate(d):
            if v:
                n[i] += v * (1 - p)
                n[i + 1] += v * p
        d = n
    return sum(d[obs:]) if obs < len(d) else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default="codon_atlas.csv")
    ap.add_argument("--record", default="resistance_record_corrected.csv")
    ap.add_argument("--reconciliation", default="species_reconciliation_final.csv")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    ensure_dirs()
    log = Tee("step9_04_codon_selection" + args.tag)
    log("=== Step 9.4  is the synonymous codon at target sites under selection? ===\n")

    # ---- CDS by protein accession -----------------------------------------
    cds = {}
    for h, s in parse_fasta(open(os.path.join(DIR_INTERIM, "target_cds.fasta"))):
        pid = h.split("protein_id=")[1].rstrip("]") if "protein_id=" in h else h.split()[0]
        cds[pid] = s.upper()
    log(f"CDS sequences available: {len(cds)}")

    # ---- documented alleles per site ---------------------------------------
    alleles = defaultdict(set)
    for r in read_csv(os.path.join(DIR_REFS, "reference_sites.csv")):
        alleles[(r["gene"], r["site_label"])].add(r["resistance_allele"])

    # ---- which species have a documented record at each site ---------------
    recon = {r["record_species"]: r["resolved_species"]
             for r in read_csv(os.path.join(DIR_PROC, args.reconciliation))
             if r["matched"] == "True"}
    documented_at = defaultdict(set)      # (gene, site) -> {species}
    any_record = set()
    for r in read_csv(os.path.join(DIR_PROC, args.record)):
        sp = recon.get(r["species"], r["species"])
        documented_at[(r["gene"], r["site_label"])].add(sp)
        any_record.add(sp)
    log(f"species appearing anywhere in the record: {len(any_record)}")

    # ---- per species x site: site codon vs within-protein background -------
    rng = random.Random(SEED)
    dropped_bg = []
    rows = []
    atlas = read_csv(os.path.join(DIR_PROC, args.atlas))
    seen = set()
    # Every filter is counted. The denominator of this analysis is not the
    # 40,312 wild-type rows of the atlas and not the 23,978 unique wild-type
    # calls either: it is the deduplicated calls that survive the filters
    # below, re-expanded over the curated candidate alleles at each site. The
    # ladder is written to codon_selection_reconciliation.csv so that the
    # Supplementary Information can state it rather than assert it.
    lad = defaultdict(int)
    n_nocds = 0
    for r in atlas:
        lad["atlas_rows"] += 1
        if r["is_wildtype"] != "True":
            continue
        lad["wild_type_rows"] += 1
        key = (r["species"], r["gene"], r["site_label"])
        if key in seen:
            continue
        seen.add(key)
        lad["unique_wild_type_calls"] += 1
        al = alleles.get((r["gene"], r["site_label"]))
        if not al:
            lad["dropped_no_curated_allele"] += 1
            continue
        lad["site_has_curated_allele"] += 1
        seq = cds.get(r["protein_accession"])
        if not seq:
            n_nocds += 1
            lad["dropped_no_cds"] += 1
            continue
        lad["cds_available"] += 1
        wt = r["wild_type_aa"]
        c_site = normalise(r["codon"])
        if CODON_TABLE.get(c_site) != wt:
            lad["dropped_codon_does_not_translate_to_wt"] += 1
            continue
        lad["codon_translates_to_wt"] += 1
        # every OTHER codon encoding the same residue in the same protein
        bg = []
        for i in range(0, len(seq) - 2, 3):
            c = seq[i:i + 3]
            if len(c) == 3 and CODON_TABLE.get(c) == wt:
                bg.append(c)
        # drop one instance of the site codon itself
        if c_site in bg:
            bg.remove(c_site)
        if len(bg) < 5:
            lad["dropped_fewer_than_5_background_codons"] += 1
            # Recorded, not merely counted. The Table S3 caption names the
            # comparisons this filter removes, and a caption that names them
            # must read them from here rather than restate them.
            dropped_bg.append({"species": r["species"], "gene": r["gene"],
                               "site_label": r["site_label"],
                               "wild_type_aa": wt, "site_codon": c_site,
                               "n_background_codons": len(bg)})
            continue
        lad["calls_entering_the_null"] += 1
        lad["allele_expansions"] += len(al)
        for a in sorted(al):
            obs = 1 if codon_reaches(c_site, a) else 0
            exp = sum(1 for c in bg if codon_reaches(c, a)) / len(bg)
            rows.append({
                "species": r["species"], "gene": r["gene"],
                "site_label": r["site_label"], "wild_type_aa": wt,
                "resistance_allele": a, "site_codon": c_site,
                "n_background_codons": len(bg),
                "site_accessible": obs,
                "background_accessible_fraction": round(exp, 6),
                "documented_here": r["species"] in documented_at[(r["gene"],
                                                                  r["site_label"])],
                "in_record_anywhere": r["species"] in any_record,
            })
    log(f"species x site x allele comparisons: {len(rows)}")
    log(f"  skipped, no CDS for the selected protein: {n_nocds}")
    lad["comparisons"] = len(rows)
    _ladder = [
        ("atlas_rows", "rows in codon_atlas.csv (long format, one per species x "
         "gene x site x curated allele)"),
        ("wild_type_rows", "of those, rows where the species carries the wild-type "
         "residue"),
        ("unique_wild_type_calls", "deduplicated to one codon call per (species, "
         "gene, site)"),
        ("dropped_no_curated_allele", "dropped: the site has no curated candidate "
         "allele"),
        ("dropped_no_cds", "dropped: no coding sequence for the selected protein"),
        ("dropped_codon_does_not_translate_to_wt", "dropped: the extracted codon "
         "does not translate to the wild-type residue"),
        ("dropped_fewer_than_5_background_codons", "dropped: fewer than five other "
         "codons for the same residue in the same protein, so no within-protein "
         "background exists"),
        ("calls_entering_the_null", "codon calls that survive every filter"),
        ("comparisons", "those calls re-expanded over the curated candidate alleles "
         "at each site; the denominator of Table S3"),
    ]
    _lp = os.path.join(DIR_RESULTS, "codon_selection_reconciliation" + args.tag + ".csv")
    with open(_lp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["stage", "count", "definition"])
        w.writeheader()
        for k, d in _ladder:
            w.writerow({"stage": k, "count": lad[k], "definition": d})
    log(f"  wrote {_lp}")
    for k, d in _ladder:
        log(f"    {k:46s} {lad[k]:>8,}")

    def block(sel, name):
        sub = [x for x in rows if sel(x)]
        if len(sub) < 10:
            log(f"\n  {name}: too few ({len(sub)})")
            return None
        o = sum(x["site_accessible"] for x in sub)
        ps = [x["background_accessible_fraction"] for x in sub]
        e = sum(ps)
        p_hi = poisson_binomial_sf(ps, o)
        p_lo = 1.0 - poisson_binomial_sf(ps, o + 1)
        log(f"\n  {name}")
        log(f"    comparisons              : {len(sub)}")
        log(f"    site codon accessible    : {o} ({o/len(sub):.1%})")
        log(f"    within-protein expectation: {e:.2f} ({e/len(sub):.1%})")
        log(f"    ratio                    : {o/e:.4f}x" if e else "")
        log(f"    exact p (excess)         : {p_hi:.4g}")
        log(f"    exact p (deficit)        : {p_lo:.4g}")
        # Resample SPECIES, not comparisons: a species contributes one row per
        # site per candidate allele, so resampling rows would treat those as
        # independent replicates. Supplementary Fig. 4c draws this interval.
        #
        # Each species' observed and expected counts are summed once, so a
        # resample costs one pass over the species rather than one pass over
        # every comparison. Written the other way round, the whole-atlas block
        # alone took longer than the rest of the script.
        lo = hi = ""
        if e:
            agg = defaultdict(lambda: [0.0, 0.0])
            for x in sub:
                v = agg[x["species"]]
                v[0] += x["site_accessible"]
                v[1] += x["background_accessible_fraction"]
            tot = [agg[k] for k in agg]
            n_sp = len(tot)
            rs = []
            for _ in range(N_BOOT):
                oo = tt = 0.0
                for _ in range(n_sp):
                    v = tot[rng.randrange(n_sp)]
                    oo += v[0]
                    tt += v[1]
                if tt > 0:
                    rs.append(oo / tt)
            rs.sort()
            lo = round(rs[int(.025 * len(rs))], 4)
            hi = round(rs[int(.975 * len(rs))], 4)
            log(f"    95% CI (species bootstrap over {n_sp} species): "
                f"[{lo}, {hi}]")
        return dict(name=name, n=len(sub), obs=o, exp=round(e, 3),
                    ratio=round(o / e, 4) if e else "", p_excess=p_hi, p_deficit=p_lo,
                    ci_lo=lo, ci_hi=hi)

    log("\n" + "=" * 74)
    log("1  ALL species: is the target-site codon unusual for its amino acid?")
    log("=" * 74)
    log("  NOTE: an excess here is EXPECTED under pure ascertainment - resistance sites")
    log("  are known because resistance evolved there. This block is a baseline, not a")
    log("  test of selection.")
    out = []
    a = block(lambda x: True, "all species x site x allele")
    if a:
        out.append(a)

    log("\n" + "=" * 74)
    log("2  THE DISCRIMINATING CONTRAST: same sites, species with vs without a record")
    log("=" * 74)
    log("  Ascertainment acts on the SITE and is shared by both groups. Selection acts on")
    log("  the SPECIES. Only selection can make these two differ.")
    b1 = block(lambda x: x["documented_here"], "species WITH a documented record at this site")
    b2 = block(lambda x: not x["documented_here"], "species WITHOUT a record at this site")
    for z in (b1, b2):
        if z:
            out.append(z)

    if b1 and b2:
        log("\n  --- comparison ---")
        log(f"    with a record   : {b1['ratio']}x  ({b1['obs']}/{b1['n']})")
        log(f"    without a record: {b2['ratio']}x  ({b2['obs']}/{b2['n']})")
        from scipy.stats import fisher_exact
        tab = [[b1["obs"], b1["n"] - b1["obs"]], [b2["obs"], b2["n"] - b2["obs"]]]
        orr, pf = fisher_exact(tab)
        log(f"    Fisher on raw accessibility: OR = {orr:.3f}, p = {pf:.4g}")
        log("    (raw comparison; the within-protein expectation differs between groups,")
        log("     so the RATIOS above are the quantities to compare, not the raw counts)")
        out.append({"name": "with_vs_without_record_fisher", "n": b1["n"] + b2["n"],
                    "obs": "", "exp": "", "ratio": round(orr, 4),
                    "p_excess": pf, "p_deficit": "", "ci_lo": "", "ci_hi": ""})

    log("\n" + "=" * 74)
    log("3  IS THE ENRICHMENT JUST THE HYPOTHESIS RESTATED? own allele vs other alleles")
    log("=" * 74)
    log("  Among species WITH a record at a site, split the comparisons by whether the")
    log("  allele being scored is the one THAT species actually carries.")
    log("    (a) the species' OWN allele  -> this IS the core hypothesis under test;")
    log("        the within-protein null is a legitimate control, but the EVENT SET is")
    log("        ascertained on resistance having evolved.")
    log("    (b) OTHER alleles at the same site -> the species did not evolve these, so")
    log("        no ascertainment on accessibility applies. If the codon at target sites")
    log("        had been shaped by insecticide selection, (b) should ALSO be enriched.")
    own = set()
    for r in read_csv(os.path.join(DIR_PROC, args.record)):
        sp = recon.get(r["species"], r["species"])
        own.add((sp, r["gene"], r["site_label"], r["resistance_allele"]))
    key4 = lambda x: (x["species"], x["gene"], x["site_label"], x["resistance_allele"])
    c1 = block(lambda x: x["documented_here"] and key4(x) in own,
               "(a) species' OWN documented allele")
    c2 = block(lambda x: x["documented_here"] and key4(x) not in own,
               "(b) OTHER alleles at the same site")
    # The same conditioned test, gene by gene. Table S3 breaks the
    # UNconditioned whole-atlas null down per gene and shows real
    # heterogeneity - RyR runs a significant deficit - so a reader cannot tell
    # whether the 1.246 above is carried by the genes whose target sites
    # already sit on more accessible codons. Only genes clearing block()'s
    # ten-comparison floor get a row; the rest are too thin to report and
    # their absence is the answer for them.
    for z in (c1, c2):
        if z:
            out.append(z)
    # Named so they cannot be taken for the pooled arm, and appended after it.
    # Figure 3 selects the pooled row from this file by exact name: a prefix
    # match with first match winning would let a per-gene row be plotted while
    # the text quoted the pool.
    for g in sorted({x["gene"] for x in rows}):
        z = block(lambda x, g=g: x["documented_here"] and key4(x) in own
                  and x["gene"] == g,
                  # No comma: the name is a CSV field, and a comma inside it
                  # forces quoting that breaks a plain cut or grep on the sheet.
                  "own documented allele per gene: %s" % g)
        if z:
            out.append(z)
    if c1 and c2:
        log("\n  READING")
        if c2["p_excess"] > 0.05:
            log(f"    (b) is NOT enriched ({c2['ratio']}x, p = {c2['p_excess']:.3g}).")
            log("    The excess in the documented group is confined to the alleles those")
            log("    species actually evolved. There is therefore NO evidence that")
            log("    insecticide selection has shaped the synonymous codon at target")
            log("    sites: the codon is simply the codon that happened to be there.")
            log(f"    (a) remains the cleanest statement of the core hypothesis")
            log(f"    ({c1['ratio']}x, p = {c1['p_excess']:.3g}), with species, gene and")
            log("    amino acid all held constant by the within-protein null.")
        else:
            log(f"    (b) IS enriched ({c2['ratio']}x, p = {c2['p_excess']:.3g}) - the codon")
            log("    at target sites is unusual even for alleles the species never evolved,")
            log("    which would be evidence of selection on the codon itself.")

    log("\n" + "=" * 74)
    log("4  BY GENE")
    log("=" * 74)
    for g in sorted({x["gene"] for x in rows}):
        z = block(lambda x, g=g: x["gene"] == g, f"gene {g}")
        if z:
            out.append(z)

    p1 = os.path.join(DIR_RESULTS, "codon_selection_per_site" + args.tag + ".csv")
    with open(p1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    p3 = os.path.join(DIR_RESULTS,
                      "codon_selection_excluded" + args.tag + ".csv")
    with open(p3, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["species", "gene", "site_label",
                                           "wild_type_aa", "site_codon",
                                           "n_background_codons"])
        w.writeheader(); w.writerows(dropped_bg)
    log("  wrote %s  (%d comparison(s) dropped for fewer than five background "
        "codons)" % (os.path.basename(p3), len(dropped_bg)))
    p2 = os.path.join(DIR_RESULTS, "codon_selection_summary" + args.tag + ".csv")
    with open(p2, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    log(f"\nwrote {p1}")
    log(f"wrote {p2}")
    log("\nStep 9.4 complete.")


if __name__ == "__main__":
    main()
