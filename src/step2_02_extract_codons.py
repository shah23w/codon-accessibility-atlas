#!/usr/bin/env python3
"""
step2_02_extract_codons.py: retrieve CDS for identified orthologs and read the codon
at every mapped resistance site, then compute mutational accessibility.

Pipeline position: consumes `data/interim/site_mappings.csv` from step2_01 and produces
the codon-accessibility atlas.

Method
------
1. Collect the unique protein accessions carrying at least one mapped site and fetch
   their coding sequences with `efetch(db=protein, rettype=fasta_cds_na)`, batched.
   This returns the exact CDS backing each protein, keyed by protein accession, so no
   genome coordinate arithmetic or exon stitching is needed.

2. For a residue at protein position p, the codon is CDS[3(p-1) : 3p]. This is VERIFIED,
   not assumed: the extracted codon must translate to the residue observed in the
   protein. Mismatches (partial CDS, annotated translational exceptions, selenocysteine,
   ribosomal frameshift) are recorded and the row is dropped, never patched.

3. Accessibility is computed with `codon_engine` for every candidate resistance allele
   at that site: one-step reachability, path counts split by transition/transversion,
   a Ts/Tv-weighted rate, and two-step path counts for the inaccessible cases.

Wild-type conditioning
----------------------
A species enters a site's analysis only if it carries the reference WILD-TYPE residue
there. Where a species already carries the resistance residue, asking whether its codon
"can reach" that residue is meaningless - it is not at risk at that site, and scoring it
as inaccessible would be a category error. Rows are retained with an explicit
`is_wildtype` flag so both the conditioned and unconditioned views are available.

Outputs
-------
data/interim/target_cds.fasta        retrieved CDS (cached; re-runs are offline)
data/processed/codon_atlas.csv       the atlas: species x gene x site x allele
data/processed/site_codon_summary.csv codon spectrum per site
logs/step2_02_extract_codons.log
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import (  # noqa: E402
    DIR_INTERIM, DIR_PROC, DIR_REFS, Tee, efetch, ensure_dirs, parse_fasta,
    write_fasta,
)
from codon_engine import (  # noqa: E402
    codon_path_counts, codon_reaches, codon_reaches_two_step, codon_weight_to,
    is_valid_codon, normalise, translate_codon,
)

KAPPA = 4.0   # transition:transversion weight; sensitivity-tested in step3


def fetch_cds(accessions, log, batch=150):
    """Fetch CDS nucleotide for protein accessions; cache to disk."""
    cache = os.path.join(DIR_INTERIM, "target_cds.fasta")
    have = {}
    if os.path.exists(cache):
        for h, s in parse_fasta(open(cache)):
            pid = None
            for tok in h.split("[") :
                if tok.startswith("protein_id="):
                    pid = tok.split("=", 1)[1].rstrip("] ")
            if pid:
                have[pid] = s
        log(f"  cache: {len(have)} CDS already retrieved")

    need = [a for a in accessions if a not in have]
    log(f"  need {len(need)} more CDS")
    for i in range(0, len(need), batch):
        chunk = need[i:i + batch]
        try:
            txt = efetch("protein", chunk, rettype="fasta_cds_na", retmode="text")
        except Exception as e:  # noqa: BLE001
            log(f"    batch {i//batch}: FAILED ({type(e).__name__}) - recorded as missing")
            continue
        n_new = 0
        for h, s in parse_fasta(txt):
            pid = None
            for tok in h.split("["):
                if tok.startswith("protein_id="):
                    pid = tok.split("=", 1)[1].rstrip("] ")
            if pid and s:
                have[pid] = s
                n_new += 1
        log(f"    batch {i//batch + 1}/{(len(need)+batch-1)//batch}: +{n_new} "
            f"(total {len(have)})")

    write_fasta(cache, [(f"[protein_id={p}]", s) for p, s in sorted(have.items())])
    log(f"  cached {len(have)} CDS to {cache}")
    return have


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=KAPPA)
    ap.add_argument("--context-floor", type=float, default=0.25,
                    help="absolute floor on site-context identity; removes broken alignments only")
    ap.add_argument("--context-drop", type=float, default=0.30,
                    help=("relative guard: drop rows whose site-context identity falls more "
                          "than this below the median for that gene x site"))
    args = ap.parse_args()

    ensure_dirs()
    log = Tee("step2_02_extract_codons")

    # ---- site -> candidate resistance alleles ------------------------------
    alleles = defaultdict(list)
    with open(os.path.join(DIR_REFS, "reference_sites.csv")) as fh:
        for r in csv.DictReader(fh):
            alleles[(r["gene"], r["site_label"])].append(r["resistance_allele"])
    log(f"reference sites: {len(alleles)}  "
        f"({sum(len(v) for v in alleles.values())} site x allele combinations)")

    with open(os.path.join(DIR_INTERIM, "site_mappings.csv")) as fh:
        all_maps = [r for r in csv.DictReader(fh) if r["mappable"] == "True"]
    log(f"mapped site records: {len(all_maps)}")

    # Context-identity filtering is RELATIVE to each gene x site, not absolute.
    # An absolute cutoff is invalid across genes because reference divergence differs:
    # the AChE reference is Torpedo californica (a fish), so a CORRECT insect ortholog
    # still shows ~0.48 median site-context identity, while VGSC (Musca reference)
    # shows ~1.00. Empirically an absolute 0.60 cutoff removed 553/565 AChE F331 rows
    # whose wild-type match rate (0.924) was indistinguishable from those retained
    # (0.917), and removed CHS1 rows that were MORE accurate than those kept
    # (0.963 vs 0.784) - i.e. it discarded correct data rather than errors.
    med = {}
    grp = defaultdict(list)
    for r in all_maps:
        grp[(r["gene"], r["site_label"])].append(float(r["context_identity"]))
    for k, v in grp.items():
        v.sort()
        med[k] = v[len(v) // 2]

    maps, rej = [], Counter()
    for r in all_maps:
        k = (r["gene"], r["site_label"])
        ci = float(r["context_identity"])
        if ci >= max(args.context_floor, med[k] - args.context_drop):
            maps.append(r)
        else:
            rej[k] += 1
    log(f"  context filter: floor {args.context_floor}, "
        f"or > {args.context_drop} below the per-site median")
    log(f"  retained {len(maps)}/{len(all_maps)} ({len(maps)/len(all_maps):.1%})")
    if rej:
        log("  excluded by low site-context identity:")
        for (g, s), n in rej.most_common(12):
            log(f"      {g:10s} {s:8s} {n:5d}  (site median {med[(g,s)]:.3f})")

    accs = sorted({m["protein_accession"] for m in maps})
    log(f"\nretrieving CDS for {len(accs)} distinct proteins...")
    cds = fetch_cds(accs, log)

    # ---- extract codons ----------------------------------------------------
    rows, drops = [], Counter()
    for m in maps:
        pid = m["protein_accession"]
        seq = cds.get(pid)
        if not seq:
            drops["no_CDS_retrieved"] += 1
            continue
        p = int(m["species_position"])
        codon = normalise(seq[3 * (p - 1): 3 * p])
        if len(codon) != 3 or not is_valid_codon(codon):
            drops["codon_out_of_range_or_ambiguous"] += 1
            continue
        aa_from_codon = translate_codon(codon)
        if aa_from_codon != m["species_aa"]:
            # CDS and protein disagree: partial CDS, translational exception, or
            # frameshift annotation. Dropped, never patched.
            drops["codon_translation_mismatch"] += 1
            continue

        wt = m["wild_type_aa"]
        is_wt = m["species_aa"] == wt
        for allele in alleles.get((m["gene"], m["site_label"]), []):
            n_paths, n_ts, n_tv = codon_path_counts(codon, allele)
            rows.append({
                "species": m["species"], "tax_id": m["tax_id"], "tier": m["tier"],
                "order": m["order"], "gene": m["gene"], "site_label": m["site_label"],
                "reference_position": m["reference_position"],
                "wild_type_aa": wt,
                "species_aa": m["species_aa"],
                "is_wildtype": is_wt,
                "species_position": p,
                "codon": codon,
                "resistance_allele": allele,
                "already_resistant": m["species_aa"] == allele,
                "one_step_accessible": codon_reaches(codon, allele),
                "n_1step_paths": n_paths,
                "n_transitions": n_ts,
                "n_transversions": n_tv,
                "ts_tv_weighted_rate": codon_weight_to(codon, allele, kappa=args.kappa),
                "n_2step_paths": codon_reaches_two_step(codon, allele),
                "protein_accession": pid,
                "align_identity": m["identity"],
                "align_coverage": m["coverage"],
                "context_identity": m["context_identity"],
            })

    log(f"\nextracted {len(rows)} atlas rows")
    if drops:
        log("dropped (recorded, never patched):")
        for k, v in drops.most_common():
            log(f"    {k:38s} {v}")

    out = os.path.join(DIR_PROC, "codon_atlas.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["gene"], r["site_label"],
                                                r["species"], r["resistance_allele"])))
    log(f"wrote {out}")

    # ---- codon spectrum per site ------------------------------------------
    log("\n--- codon spectrum at each site (wild-type species only) ---")
    summ = []
    seen = set()
    for r in rows:
        key = (r["gene"], r["site_label"], r["species"])
        if key in seen or not r["is_wildtype"]:
            continue
        seen.add(key)
    for (gene, site) in sorted(alleles):
        sub = {}
        for r in rows:
            if r["gene"] == gene and r["site_label"] == site and r["is_wildtype"]:
                sub[r["species"]] = r["codon"]
        if not sub:
            continue
        cnt = Counter(sub.values())
        spectrum = ", ".join(f"{c}:{n}" for c, n in cnt.most_common())
        log(f"  {gene:10s} {site:7s} n={len(sub):4d}  {len(cnt)} codons  {spectrum}")
        summ.append({"gene": gene, "site_label": site, "n_wildtype_species": len(sub),
                     "n_distinct_codons": len(cnt), "codon_spectrum": spectrum})
    sp = os.path.join(DIR_PROC, "site_codon_summary.csv")
    with open(sp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0].keys()))
        w.writeheader(); w.writerows(summ)
    log(f"\nwrote {sp}")
    log.close()


if __name__ == "__main__":
    main()
