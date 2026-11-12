#!/usr/bin/env python3
"""
step2_01_identify_orthologs.py: find target-gene orthologs and map resistance sites.

For every species proteome, for every reference gene, this identifies the best-matching
protein(s) and transfers each verified reference site position onto that protein by
pairwise alignment. No external aligner is required: Biopython's PairwiseAligner (C
implementation) does the work, which also avoids the MSA-quality confound that a
594-taxon multiple alignment would introduce.

Method
------
1. k-mer prescreen. Index the reference's 5-mers; score every protein in the proteome by
   shared-5-mer count normalised by length. Cheap (set intersections) and reduces ~30,000
   candidates to a handful before any alignment is attempted.

2. Semi-global pairwise alignment (BLOSUM62, affine gaps, end gaps free) of the reference
   against each shortlisted candidate. End-gap-free scoring matters because annotated
   proteins are frequently N-terminally truncated or carry extra first exons; penalising
   those overhangs would demote true orthologs.

3. Coordinate transfer. The alignment's aligned blocks give an exact reference-index ->
   species-index map. A site is reported ONLY if it falls inside an aligned block; if it
   lands in a gap the site is recorded as unmappable for that species. It is never
   approximated by interpolation.

4. Paralogy is not assumed away. The top `--keep` candidates are retained per gene with
   their alignment statistics, so paralog structure (e.g. insect ace-1 vs ace-2, or the
   several GABA/GluCl subunits) is visible downstream rather than silently collapsed to
   a single best hit.

Outputs
-------
data/interim/ortholog_hits.csv     every retained candidate with alignment statistics
data/interim/site_mappings.csv     species x gene x site -> protein accession, position, residue
logs/step2_01_identify_orthologs.log
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import sys
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import (  # noqa: E402
    DIR_GENOMES, DIR_INTERIM, DIR_PROC, DIR_REFS, Tee, ensure_dirs, parse_fasta,
)

from Bio.Align import PairwiseAligner, substitution_matrices  # noqa: E402
# The orthologue-selection rule, shared with step4_02's replicate arm so both arms
# choose a species' protein by the same criteria. Verified against this atlas by
# data09_verify_selection_refactor.py.
from ortholog_selection import (  # noqa: E402
    IDXCOV_TOL, MIN_MEAN_CONTEXT, select_one_protein_per_species_gene,
)

KMER = 5
CONTEXT = 10        # +/- residues around a site used for local context identity
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# How far down the k-mer ranking to look before building the shortlist. The old code
# stopped at keep*2 = 6, which is fewer slots than a well-annotated locus has isoforms.
PRESCREEN_K = 40
# Two candidates whose 5-mer Jaccard is at or above this are treated as isoforms of one
# locus and collapsed to the higher-scoring representative. Isoforms of a gene share
# nearly all their k-mers; genuine paralogs (insect ace-1 vs ace-2, ~40% identity) do not
# come close.
ISOFORM_JACCARD = 0.80
# MIN_MEAN_CONTEXT (a candidate whose mean local site-context identity falls below it is
# not treated as the gene at all) and IDXCOV_TOL (the band inside which context identity
# decides) are imported from ortholog_selection above, NOT re-declared here: the replicate
# arm of step4_02 must use the same values, and a second copy of a threshold is exactly how
# the two arms drifted apart.
# A paralogue call needs both an absolute and a relative floor: taxa that predate the
# insect ace-1/ace-2 duplication must come back unclassified rather than forced.
PARA_MIN_SCORE = 0.25
PARA_MIN_MARGIN = 0.03


# --------------------------------------------------------------------------- #
def make_aligner():
    a = PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11.0
    a.extend_gap_score = -1.0
    a.mode = "global"
    # End gaps free: tolerate N/C-terminal truncation and alternative first exons.
    # Biopython >=1.85 renamed these; set the new names when available.
    try:
        a.end_insertion_score = 0.0
        a.end_deletion_score = 0.0
    except AttributeError:                      # older Biopython
        a.target_end_gap_score = 0.0
        a.query_end_gap_score = 0.0
    return a


def sanitise(seq: str) -> str:
    """BLOSUM62 has no entry for U/O/J/B/Z/*; map anything non-standard to X."""
    return "".join(c if c in STANDARD_AA else "X" for c in seq.upper())


def kmers(seq: str, k: int = KMER) -> set:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def _jaccard(a: set, b: set) -> float:
    """Set similarity, used only to recognise isoforms of one locus."""
    if not a or not b:
        return 0.0
    inter = len(a & b)
    return inter / (len(a) + len(b) - inter)


def load_references(refs_fasta: str, sites_csv: str):
    refs = {}
    for h, s in parse_fasta(open(refs_fasta)):
        gene = h.split("|")[0]
        refs[gene] = sanitise(s)
    sites = defaultdict(list)
    with open(sites_csv) as fh:
        for r in csv.DictReader(fh):
            key = (r["gene"], r["site_label"])
            if key not in {(g, sl) for g, sl, *_ in
                           [(x["gene"], x["site_label"]) for x in []]}:
                pass
            sites[r["gene"]].append(r)
    # collapse to unique sites per gene (the CSV has one row per allele)
    uniq = {}
    for gene, rows in sites.items():
        seen = {}
        for r in rows:
            seen[r["site_label"]] = {
                "site_label": r["site_label"],
                "reference_position": int(r["reference_position"]),
                "published_position": int(r["published_position"]),
                "wild_type_aa": r["wild_type_aa"],
            }
        uniq[gene] = sorted(seen.values(), key=lambda d: d["reference_position"])
    return refs, uniq


def read_proteome(zip_path: str):
    with zipfile.ZipFile(zip_path) as zf:
        names = [n for n in zf.namelist() if n.endswith("protein.faa")]
        if not names:
            return None
        out = {}
        with zf.open(names[0]) as fh:
            txt = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
            for h, s in parse_fasta(txt):
                pid = h.split()[0]
                desc = h.split(" ", 1)[1] if " " in h else ""
                out[pid] = (desc, s)
    return out


def build_index(ref_positions):
    """Map reference index -> query index from an alignment's aligned blocks."""
    ref_blocks, qry_blocks = ref_positions
    mapping = {}
    for (rs, re_), (qs, qe) in zip(ref_blocks, qry_blocks):
        for off in range(re_ - rs):
            mapping[rs + off] = qs + off
    return mapping


def process_species(args):
    (species, tax_id, tier, order, acc, zip_path, refs, sites, keep, min_identity,
     para_refs) = args
    try:
        prot = read_proteome(zip_path)
    except Exception as e:  # noqa: BLE001
        return species, [], [], f"proteome read error: {type(e).__name__}: {e}"
    if not prot:
        return species, [], [], "no protein.faa in archive"

    aligner = make_aligner()
    clean = {pid: (desc, sanitise(s)) for pid, (desc, s) in prot.items() if len(s) >= 80}

    # ---- 1. k-mer prescreen, single pass over the proteome ----------------
    # Scored in ONE pass with each protein's k-mer set held only transiently. Caching a
    # k-mer set per protein would cost ~1.3 GB per species (30k proteins x ~400 5-mers),
    # which across parallel workers would exhaust memory on an 8-core desktop.
    ref_kmers = {gene: kmers(ref) for gene, ref in refs.items()}
    # The AChE paralogue references are prescreened TOO, not only used for classification.
    # Without this the shortlist is built from Torpedo similarity alone, so an ace-1
    # ortholog that happens to score poorly against a fish AChE never reaches the
    # shortlist and the paralogue step never gets to see it. Measured: Plutella xylostella
    # "acetylcholinesterase type 1" scores 0.4300 against Torpedo (below its own type 2 at
    # 0.4198 on the per-site criterion) but 0.5748 against insect ace-1. Before this fix,
    # 25 of 74 Lepidoptera were assigned ace-2 for exactly this reason.
    para_kmers = {f"__para_{lbl}": kmers(seq) for lbl, seq in (para_refs or {}).items()}
    scan = {**ref_kmers, **para_kmers}
    scored_by_gene = {gene: [] for gene in scan}
    for pid, (desc, s) in clean.items():
        pk = kmers(s)
        for gene, rk in scan.items():
            shared = len(rk & pk)
            if shared:
                # Normalise by the REFERENCE k-mer count, not by the shorter of the two
                # sequences. The old form, min(len(ref), len(s)), divides a short fragment
                # by its own tiny length and inflates it: a 90-aa esterase fragment could
                # outrank a full-length ortholog. Measured consequence before this fix -
                # for Palaemon carinicauda the real AChE ranked 319th while the shortlist
                # filled with cocaine esterase and juvenile-hormone esterase fragments.
                scored_by_gene[gene].append((shared / len(rk), shared, pid))
        del pk

    # ---- 1b. re-derive k-mer sets for the top candidates only --------------
    # Needed for the isoform-collapsing step below. Held for at most
    # PRESCREEN_K x n_genes proteins (~200), not the whole proteome, so the memory
    # argument above still holds.
    top_pids = set()
    for gene in scored_by_gene:
        scored_by_gene[gene].sort(reverse=True)
        top_pids.update(pid for _sc, _sh, pid in scored_by_gene[gene][:PRESCREEN_K])
    top_kmers = {pid: kmers(clean[pid][1]) for pid in top_pids}

    hits_out, maps_out = [], []
    for gene, ref in refs.items():
        scored = scored_by_gene[gene]
        if not scored:
            continue
        # Shortlist by DIVERSITY, not by rank alone.
        #
        # The previous version took scored[:keep*2] - six slots - on the stated
        # assumption that "the true ortholog is essentially always rank 1 by shared-k-mer
        # density". That assumption fails exactly when a paralog is annotated with more
        # isoforms than the target: for Anopheles gambiae AChE, all six slots were ace-2
        # isoforms and the true ace-1 ranked 7th, 8th and 10th, so it was never aligned.
        # Had it been aligned it would have won outright (identity 0.459 vs 0.415).
        #
        # Fix: walk further down the ranked list and admit a candidate only if it is not
        # a near-duplicate of one already admitted. Isoforms of one locus collapse to a
        # single representative, so the shortlist spans distinct loci and a genuine
        # paralog can no longer be crowded out by isoform redundancy.
        # For AChE the candidate pool is the UNION over the Torpedo reference and the two
        # paralogue references, so an ace-1 that only a paralogue reference recognises is
        # still available to the classification step. Each source contributes its own
        # best candidate before the list is filled by rank, so no source can be crowded
        # out by another.
        sources = [scored]
        if gene == "AChE":
            sources += [scored_by_gene[k] for k in scored_by_gene
                        if k.startswith("__para_")]
        pool, seen_pid = [], set()
        for src in sources:                       # seed: best from each source
            for sc, sh, pid in src[:1]:
                if pid not in seen_pid:
                    seen_pid.add(pid)
                    pool.append((sc, sh, pid))
        best_score = scored[0][0]
        for sc, sh, pid in scored[:PRESCREEN_K]:  # then fill by primary-reference rank
            if sc < 0.25 * best_score:
                break
            if pid not in seen_pid:
                seen_pid.add(pid)
                pool.append((sc, sh, pid))

        shortlist = []
        for _sc, _sh, pid in pool:
            pk = top_kmers.get(pid)
            if pk is not None and any(
                    _jaccard(pk, top_kmers[q]) >= ISOFORM_JACCARD
                    for q in shortlist if q in top_kmers):
                continue                      # isoform of an already-admitted locus
            shortlist.append(pid)
            if len(shortlist) >= keep * 2 + len(sources) - 1:
                break

        # ---- 2. align shortlist -------------------------------------------
        aligned = []
        for pid in shortlist:
            desc, s = clean[pid]
            try:
                aln = aligner.align(ref, s)[0]
            except Exception:  # noqa: BLE001
                continue
            idx = build_index(aln.aligned)
            n_aln = len(idx)
            n_id = sum(1 for ri, qi in idx.items() if ref[ri] == s[qi])
            ident = n_id / n_aln if n_aln else 0.0
            cov = n_aln / len(ref)
            aligned.append((ident * cov, ident, cov, pid, desc, s, idx, aln.score))
        if not aligned:
            continue
        aligned.sort(reverse=True, key=lambda t: t[0])

        # ---- 2b. AChE paralogue classification (ace-1 vs ace-2) ------------
        # Which AChE gene is the insecticide target is CLADE-DEPENDENT: most arthropods
        # carry both ace-1 and ace-2 and ace-1 is the OP/carbamate target, but
        # cyclorrhaphan Diptera (Drosophila, Musca, Stomoxys, Bactrocera) lost ace-1 and
        # their single Ace IS the target. Ranking against the Torpedo reference cannot
        # make that distinction - measured: for Plutella xylostella, Torpedo scores
        # type-2 0.4198 and type-1 0.4300 (nearly tied, and the wrong one won on the
        # per-site criterion), whereas an insect ace-1 reference scores them 0.3553 and
        # 0.5748 - a decisive separation.
        #
        # So candidates are classified against insect ace-1 and ace-2 references, and
        # selection prefers ace-1 where the species has one. Torpedo remains the
        # numbering frame, because "G119" IS Torpedo mature numbering in the literature.
        para_of, para_score = {}, {}
        if gene == "AChE" and para_refs:
            for _c, _i, _cv, pid, _d, s, _ix, _sc in aligned:
                best, best_lbl, scores = -1.0, "", {}
                for lbl, pref in para_refs.items():
                    try:
                        a2 = aligner.align(pref, s)[0]
                    except Exception:  # noqa: BLE001
                        continue
                    i2 = build_index(a2.aligned)
                    if not i2:
                        continue
                    sc = (sum(1 for ri, qi in i2.items() if pref[ri] == s[qi])
                          / len(i2)) * (len(i2) / len(pref))
                    scores[lbl] = sc
                    if sc > best:
                        best, best_lbl = sc, lbl
                # Only call a paralogue when the evidence supports it. Non-insect taxa
                # (Acari, Crustacea, Chelicerata) predate the insect ace-1/ace-2
                # duplication, so forcing them into one bin would be an invented fact.
                # Tetranychus urticae, for instance, is ~equally distant from both.
                second = max([v for k, v in scores.items() if k != best_lbl] or [0.0])
                if best < PARA_MIN_SCORE or (best - second) < PARA_MIN_MARGIN:
                    best_lbl = ""
                para_of[pid] = best_lbl
                para_score[pid] = round(best, 4)

        # ---- 3. record top candidates + map sites --------------------------
        for rank, (comb, ident, cov, pid, desc, s, idx, score) in enumerate(aligned[:keep], 1):
            hits_out.append({
                "species": species, "tax_id": tax_id, "tier": tier, "order": order,
                "assembly_accession": acc, "gene": gene, "rank": rank,
                "protein_accession": pid, "protein_length": len(s),
                "description": desc[:120], "align_score": round(score, 1),
                "identity": round(ident, 4), "coverage": round(cov, 4),
                "identity_x_coverage": round(comb, 4),
                "passes_threshold": ident >= min_identity,
                "paralogue": para_of.get(pid, ""),
                "paralogue_score": para_score.get(pid, ""),
            })
            # Map sites for EVERY retained candidate, not just rank 1. The final
            # per-site choice is then made downstream on local context identity, which
            # is the criterion that actually bears on site correctness. Ranking by
            # global identity x coverage alone can favour an adjacent gene fragment or
            # a close paralog.
            if ident < min_identity:
                continue
            for site in sites.get(gene, []):
                ri = site["reference_position"] - 1
                qi = idx.get(ri)
                # Local site-context identity over a +/-10 residue window.
                # Global identity is the wrong criterion for SITE correctness: a
                # close paralog (glycine receptor vs GluCl) or an adjacent gene
                # fragment can score well globally while placing the site wrongly.
                # The flanking context around a genuine target site is near-invariant,
                # so local identity discriminates where global identity cannot.
                lo_r, hi_r = max(0, ri - CONTEXT), min(len(ref), ri + CONTEXT + 1)
                ctx_pairs = [(ref[r], s[idx[r]]) for r in range(lo_r, hi_r) if r in idx]
                ctx_ident = (sum(1 for a, b in ctx_pairs if a == b) / len(ctx_pairs)
                             if ctx_pairs else 0.0)
                maps_out.append({
                    "species": species, "tax_id": tax_id, "tier": tier, "order": order,
                    "gene": gene, "rank": rank, "site_label": site["site_label"],
                    "paralogue": para_of.get(pid, ""),
                    "reference_position": site["reference_position"],
                    "wild_type_aa": site["wild_type_aa"],
                    "protein_accession": pid,
                    "species_position": (qi + 1) if qi is not None else "",
                    "species_aa": s[qi] if qi is not None else "",
                    "matches_reference_wt": (s[qi] == site["wild_type_aa"]) if qi is not None else False,
                    "mappable": qi is not None,
                    "identity": round(ident, 4), "coverage": round(cov, 4),
                    "context_identity": round(ctx_ident, 4),
                    "context_aligned_len": len(ctx_pairs),
                    "site_context_ref": ref[lo_r:hi_r],
                    "site_context_species": "".join(
                        s[idx[r]] if r in idx else "-" for r in range(lo_r, hi_r)),
                })
    return species, hits_out, maps_out, ""


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--keep", type=int, default=3, help="candidates retained per gene")
    ap.add_argument("--min-identity", type=float, default=0.35,
                    help="minimum aligned identity for a hit to be used for site mapping")
    ap.add_argument("--limit", type=int, default=0, help="process only N species (testing)")
    args = ap.parse_args()

    ensure_dirs()
    log = Tee("step2_01_identify_orthologs")

    refs, sites = load_references(
        os.path.join(DIR_REFS, "reference_proteins.fasta"),
        os.path.join(DIR_REFS, "reference_sites.csv"))
    log(f"references: {len(refs)} genes, "
        f"{sum(len(v) for v in sites.values())} sites")
    for g in refs:
        log(f"    {g:10s} {len(refs[g]):6d} aa   sites: "
            f"{', '.join(s['site_label'] for s in sites.get(g, []))}")

    # AChE paralogue references: insect ace-1 and ace-2, used ONLY to decide which
    # candidate is the target gene. Torpedo stays the numbering frame.
    para_refs = {}
    para_path = os.path.join(DIR_REFS, "ache_paralogue_refs.fasta")
    if os.path.exists(para_path):
        for h, seq in parse_fasta(open(para_path)):
            para_refs[h.split("|")[0]] = sanitise(seq)
        log(f"AChE paralogue references: "
            f"{ {k: len(v) for k, v in sorted(para_refs.items())} }")
    else:
        log("WARNING: no ache_paralogue_refs.fasta; AChE paralogue preference disabled")

    with open(os.path.join(DIR_PROC, "species_panel.csv")) as fh:
        panel = list(csv.DictReader(fh))
    tasks = []
    for r in panel:
        zp = os.path.join(DIR_GENOMES, f"{r['assembly_accession']}.zip")
        if os.path.exists(zp):
            tasks.append((r["species"], r["tax_id"], r["tier"], r["order"],
                          r["assembly_accession"], zp, refs, sites,
                          args.keep, args.min_identity, para_refs))
    if args.limit:
        tasks = tasks[:args.limit]
    log(f"\nproteomes available: {len(tasks)} / {len(panel)} panel species")
    log(f"workers: {args.workers}\n")

    all_hits, all_maps, errors = [], [], []
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(process_species, t) for t in tasks]
        for fut in as_completed(futs):
            species, hits, maps, err = fut.result()
            done += 1
            if err:
                errors.append({"species": species, "error": err})
                log(f"  [{done}/{len(tasks)}] {species}: {err}")
            all_hits.extend(hits)
            all_maps.extend(maps)
            if done % 25 == 0 or done == len(tasks):
                log(f"  [{done}/{len(tasks)}] hits={len(all_hits)} mappings={len(all_maps)}")

    hp = os.path.join(DIR_INTERIM, "ortholog_hits.csv")
    with open(hp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_hits[0].keys()))
        w.writeheader(); w.writerows(sorted(all_hits, key=lambda r: (r["species"], r["gene"], r["rank"])))
    log(f"\nwrote {hp}  ({len(all_hits)} rows)")

    mp = os.path.join(DIR_INTERIM, "site_mappings_all_candidates.csv")
    with open(mp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_maps[0].keys()))
        w.writeheader(); w.writerows(sorted(all_maps, key=lambda r: (r["species"], r["gene"], r["site_label"], r["rank"])))
    log(f"wrote {mp}  ({len(all_maps)} candidate rows)")

    # ---- select ONE candidate per species x gene x site --------------------
    # Chosen on LOCAL site-context identity first. The flanking residues of a genuine
    # target site are near-invariant across arthropods, so context identity separates
    # the true site from a close paralog (e.g. glycine receptor vs GluCl) or an
    # adjacent gene fragment, which global identity x coverage cannot.
    # ONE PROTEIN per species x gene, then all of that protein's sites.
    #
    # The previous version keyed on (species, gene, SITE), choosing independently at every
    # site. That let one species read its five AChE sites off five different proteins:
    # measured before this fix, 325 of 586 species used more than one protein across the
    # AChE sites, 317 had a pair less than 70% identical to each other, and 123 site-calls
    # in 65 species came from non-cholinesterases (juvenile-hormone esterase, cocaine
    # esterase, and in one species an ABC transporter). A site set read off different loci
    # is not a gene.
    #
    # Ranking: global identity x coverage first, mean site-context identity second,
    # accession last as a deterministic tie-break.
    #
    # Ranking by NUMBER OF SITES PLACED first is wrong: for Tetranychus urticae
    # it selects
    # XP_025017031.1 purely because it placed 5 sites, over two full-length
    # acetylcholinesterases that placed 4 but beat it on BOTH quality measures
    # (mean context 0.536 / 0.508 vs 0.447; identity x coverage 0.3362 / 0.3379 vs
    # 0.3158). Counting sites rewards a protein for spanning the reference, which a
    # long paralog does trivially.
    #
    # A candidate whose site contexts are poor is excluded outright rather than
    # out-ranked, so a globally-similar but site-wrong paralog cannot win on identity
    # alone. Disagreements between the two orderings are logged, not hidden.
    # The rule itself is in ortholog_selection.py so that the replicate arm of
    # step4_02_assembly_replication.py applies the SAME one, so the primary and
    # replicate arms cannot select proteins by different rules.
    all_maps, _sel = select_one_protein_per_species_gene(all_maps, log=log)

    mp = os.path.join(DIR_INTERIM, "site_mappings.csv")
    with open(mp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_maps[0].keys()))
        w.writeheader(); w.writerows(all_maps)
    log(f"wrote {mp}  ({len(all_maps)} rows)")

    ctx = sorted(m["context_identity"] for m in all_maps)
    log(f"context identity: median {ctx[len(ctx)//2]:.3f}, "
        f"10th pct {ctx[len(ctx)//10]:.3f}, "
        f">=0.7 in {sum(1 for c in ctx if c >= 0.7)}/{len(ctx)} "
        f"({sum(1 for c in ctx if c >= 0.7)/len(ctx):.1%})")

    # ---- QC summary --------------------------------------------------------
    log("\n--- per-gene recovery ---")
    log(f"  {'gene':10s} {'species':>8s} {'sites':>7s} {'mapped':>7s} {'wt_match':>9s} {'med_ident':>10s}")
    import statistics
    for gene in refs:
        gh = [h for h in all_hits if h["gene"] == gene and h["rank"] == 1 and h["passes_threshold"]]
        gm = [m for m in all_maps if m["gene"] == gene]
        mapped = [m for m in gm if m["mappable"]]
        wt = [m for m in mapped if m["matches_reference_wt"]]
        med = statistics.median([h["identity"] for h in gh]) if gh else 0
        log(f"  {gene:10s} {len(gh):8d} {len(gm):7d} {len(mapped):7d} "
            f"{len(wt):9d} {med:10.3f}")
    # ---- per-site conservation --------------------------------------------
    # A site is only comparable across species if the WILD-TYPE residue is conserved.
    # Where a species already carries the resistance residue, asking whether it "can
    # mutate to" that residue is meaningless: such species are not at risk at that site
    # and must be excluded from it, NOT scored as inaccessible. This filter is applied
    # in step2_02; here we quantify it so the threshold is chosen from data.
    log("\n--- per-site wild-type conservation (drives the inclusion filter) ---")
    log(f"  {'gene':10s} {'site':7s} {'wt':>3s} {'mapped':>7s} {'wt_frac':>8s}  most common residues")
    site_cons = []
    for gene in refs:
        for site in sites.get(gene, []):
            sl = site["site_label"]
            rows = [m for m in all_maps if m["gene"] == gene and m["site_label"] == sl and m["mappable"]]
            if not rows:
                continue
            cnt = Counter(m["species_aa"] for m in rows)
            wt = site["wild_type_aa"]
            frac = cnt.get(wt, 0) / len(rows)
            top = ", ".join(f"{a}:{n}" for a, n in cnt.most_common(4))
            log(f"  {gene:10s} {sl:7s} {wt:>3s} {len(rows):7d} {frac:8.3f}  {top}")
            site_cons.append({"gene": gene, "site_label": sl, "wild_type_aa": wt,
                              "n_mapped": len(rows), "wt_fraction": round(frac, 4),
                              "residue_counts": ";".join(f"{a}:{n}" for a, n in cnt.most_common())})
    cp = os.path.join(DIR_INTERIM, "site_conservation.csv")
    with open(cp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(site_cons[0].keys()))
        w.writeheader(); w.writerows(site_cons)
    log(f"\nwrote {cp}")

    if errors:
        log(f"\n{len(errors)} species errored (recorded, not substituted)")
    log.close()


if __name__ == "__main__":
    main()
