#!/usr/bin/env python3
"""
step4_02_assembly_replication.py: replicate every codon call in an INDEPENDENT assembly.

The atlas reads one codon per species from one reference assembly. Two distinct
objections follow, and this step addresses both with data rather than argument:

  (1) TECHNICAL: is the codon call reproducible, or an artefact of one assembly's
      sequencing, assembly, or annotation pipeline?
  (2) BIOLOGICAL: synonymous sites are polymorphic, so a single reference genome may
      not represent the species. This is the "reference-genome codons are not
      field-population codons" limitation stated in the proposal.

122 species in the panel have >= 2 independently annotated assemblies (148 extra
assemblies in total), frequently from different isolates, laboratories, and submitters.
Re-running the entire ortholog -> site -> codon pipeline on the SECOND assembly and
comparing gives a direct, quantitative measure of both error rates.

Two outcomes are reported separately, because they mean different things:
  * amino-acid disagreement  -> a site-mapping or annotation problem (technical);
  * synonymous codon disagreement at an identical amino acid -> genuine within-species
    codon polymorphism (biological), which is exactly the quantity the accessibility
    argument depends on.

Usage
-----
    python3 step4_02_assembly_replication.py --download     # fetch replicate assemblies
    python3 step4_02_assembly_replication.py                # compare (after download)

Outputs
-------
data/interim/replicate_assemblies.csv    the replicate assembly chosen per species
data/raw/genomes_replicate/*.zip         downloaded replicate proteomes
results/assembly_replication.csv         per site x species agreement
results/assembly_replication_summary.csv headline agreement rates
logs/step4_02_assembly_replication.log
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
import time
import urllib.request
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import (  # noqa: E402
    DATASETS, DIR_INTERIM, DIR_PROC, DIR_RAW, DIR_REFS, DIR_RESULTS, Tee, USER_AGENT,
    efetch, ensure_dirs, parse_fasta, species_key,
)
from step2_01_identify_orthologs import (  # noqa: E402
    CONTEXT, KMER, build_index, kmers, load_references, make_aligner, sanitise,
)
from codon_engine import is_valid_codon, normalise, translate_codon  # noqa: E402

DIR_REP = os.path.join(DIR_RAW, "genomes_replicate")
LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}


def choose_replicates(log):
    """For each panel species, pick the best assembly that is NOT the one already used."""
    with open(os.path.join(DIR_PROC, "species_panel.csv")) as fh:
        used = {r["species"]: r["assembly_accession"] for r in csv.DictReader(fh)}
    with open(os.path.join(DIR_INTERIM, "assemblies_all.json")) as fh:
        allasm = json.load(fh)

    cand = defaultdict(list)
    for src in ("refseq", "genbank"):
        for rep in allasm.get(src, []) or []:
            sp = species_key((rep.get("organism") or {}).get("organism_name", ""))
            if sp not in used or rep.get("accession") == used[sp]:
                continue
            ai = rep.get("assembly_info", {}) or {}
            st = rep.get("assembly_stats", {}) or {}
            try:
                n50 = int(st.get("contig_n50") or 0)
            except (TypeError, ValueError):
                n50 = 0
            cand[sp].append((LEVEL_RANK.get(ai.get("assembly_level", ""), 0), n50,
                             rep.get("accession"), src, ai.get("assembly_level", ""),
                             ai.get("submitter", "")))
    out = []
    for sp, lst in sorted(cand.items()):
        lst.sort(reverse=True)
        lvl, n50, acc, src, level, sub = lst[0]
        out.append({"species": sp, "primary_accession": used[sp],
                    "replicate_accession": acc, "source": src,
                    "assembly_level": level, "contig_n50": n50,
                    "submitter": (sub or "")[:60],
                    "n_alternatives": len(lst)})
    log(f"species with an independent replicate assembly: {len(out)}")
    path = os.path.join(DIR_INTERIM, "replicate_assemblies.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    log(f"wrote {path}")
    return out


def dl_one(acc, tries=4):
    path = os.path.join(DIR_REP, f"{acc}.zip")
    if os.path.exists(path) and os.path.getsize(path) > 10000:
        try:
            with zipfile.ZipFile(path) as zf:
                if any(n.endswith("protein.faa") for n in zf.namelist()):
                    return acc, "cached", os.path.getsize(path)
        except zipfile.BadZipFile:
            pass
    url = (f"{DATASETS}/genome/accession/{acc}/download"
           f"?include_annotation_type=PROT_FASTA&filename={acc}.zip")
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            tmp = path + ".part"
            with urllib.request.urlopen(req, timeout=900) as r, open(tmp, "wb") as fh:
                while True:
                    ch = r.read(1 << 20)
                    if not ch:
                        break
                    fh.write(ch)
            os.replace(tmp, path)
            with zipfile.ZipFile(path) as zf:
                if any(n.endswith("protein.faa") for n in zf.namelist()):
                    return acc, "downloaded", os.path.getsize(path)
            os.remove(path)
        except Exception:  # noqa: BLE001
            time.sleep(3 * (a + 1))
    return acc, "FAILED", 0


def map_sites_in_assembly(args):
    """Run the step2_01 mapping logic on one replicate assembly."""
    species, acc, refs, sites = args
    path = os.path.join(DIR_REP, f"{acc}.zip")
    if not os.path.exists(path):
        return species, acc, [], "archive missing"
    try:
        with zipfile.ZipFile(path) as zf:
            names = [n for n in zf.namelist() if n.endswith("protein.faa")]
            if not names:
                return species, acc, [], "no protein.faa"
            prot = {}
            with zf.open(names[0]) as fh:
                txt = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
                for h, s in parse_fasta(txt):
                    prot[h.split()[0]] = sanitise(s)
    except Exception as e:  # noqa: BLE001
        return species, acc, [], f"{type(e).__name__}: {e}"

    clean = {p: s for p, s in prot.items() if len(s) >= 80}
    if not clean:
        return species, acc, [], "empty proteome"
    aligner = make_aligner()
    ref_kmers = {g: kmers(r) for g, r in refs.items()}
    scored = {g: [] for g in refs}
    for pid, s in clean.items():
        pk = kmers(s)
        for g, rk in ref_kmers.items():
            sh = len(rk & pk)
            if sh:
                scored[g].append((sh / max(1, min(len(refs[g]), len(s)) - KMER + 1), pid))
        del pk

    out = []
    for gene, ref in refs.items():
        sc = scored[gene]
        if not sc:
            continue
        sc.sort(reverse=True)
        best = sc[0][0]
        shortlist = [p for v, p in sc[:6] if v >= 0.25 * best]
        cands = []
        for pid in shortlist:
            s = clean[pid]
            try:
                aln = aligner.align(ref, s)[0]
            except Exception:  # noqa: BLE001
                continue
            idx = build_index(aln.aligned)
            if not idx:
                continue
            nid = sum(1 for ri, qi in idx.items() if ref[ri] == s[qi])
            cands.append((nid / len(idx), len(idx) / len(ref), pid, s, idx))
        if not cands:
            continue
        for site in sites.get(gene, []):
            ri = site["reference_position"] - 1
            bestrow = None
            for ident, cov, pid, s, idx in cands:
                qi = idx.get(ri)
                if qi is None:
                    continue
                lo, hi = max(0, ri - CONTEXT), min(len(ref), ri + CONTEXT + 1)
                pairs = [(ref[r], s[idx[r]]) for r in range(lo, hi) if r in idx]
                ctx = sum(1 for a, b in pairs if a == b) / len(pairs) if pairs else 0.0
                key = (ctx, ident * cov)
                if bestrow is None or key > bestrow[0]:
                    bestrow = (key, pid, qi + 1, s[qi], ident, cov, ctx)
            if bestrow:
                _k, pid, pos, aa, ident, cov, ctx = bestrow
                out.append({"species": species, "assembly": acc, "gene": gene,
                            "site_label": site["site_label"], "protein_accession": pid,
                            "position": pos, "aa": aa, "identity": round(ident, 4),
                            "coverage": round(cov, 4), "context_identity": round(ctx, 4)})
    return species, acc, out, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()

    ensure_dirs()
    os.makedirs(DIR_REP, exist_ok=True)
    log = Tee("step4_02_assembly_replication")
    log("=== Step 4.2  independent-assembly replication ===\n")

    reps = choose_replicates(log)

    if args.download:
        log(f"\ndownloading {len(reps)} replicate proteomes...")
        t0, done, sizes = time.time(), 0, []
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(dl_one, r["replicate_accession"]): r for r in reps}
            for f in as_completed(futs):
                acc, status, n = f.result()
                done += 1
                sizes.append(n)
                if status == "FAILED":
                    log(f"  FAILED {acc}")
                elif done % 20 == 0:
                    log(f"  [{done}/{len(reps)}] {sum(sizes)/1e9:.2f} GB "
                        f"({(time.time()-t0)/60:.1f} min)")
        log(f"download complete: {sum(sizes)/1e9:.2f} GB in {(time.time()-t0)/60:.1f} min")
        return

    # ---- compare -----------------------------------------------------------
    refs, sites = load_references(os.path.join(DIR_REFS, "reference_proteins.fasta"),
                                 os.path.join(DIR_REFS, "reference_sites.csv"))
    have = {r["species"]: r["replicate_accession"] for r in reps
            if os.path.exists(os.path.join(DIR_REP, f"{r['replicate_accession']}.zip"))}
    log(f"\nreplicate proteomes on disk: {len(have)}")
    if not have:
        log("none downloaded - run with --download first")
        return

    tasks = [(sp, acc, refs, sites) for sp, acc in sorted(have.items())]
    rep_calls, errs = [], []
    done = 0
    with ProcessPoolExecutor(max_workers=max(1, (os.cpu_count() or 4) - 2)) as ex:
        futs = [ex.submit(map_sites_in_assembly, t) for t in tasks]
        for f in as_completed(futs):
            sp, acc, rows, err = f.result()
            done += 1
            if err:
                errs.append((sp, acc, err))
            rep_calls.extend(rows)
            if done % 20 == 0 or done == len(tasks):
                log(f"  [{done}/{len(tasks)}] replicate site calls: {len(rep_calls)}")

    # primary calls, with the CODON (from the atlas)
    primary = {}
    with open(os.path.join(DIR_PROC, "codon_atlas.csv")) as fh:
        for r in csv.DictReader(fh):
            primary[(r["species"], r["gene"], r["site_label"])] = (
                r["species_aa"], r["codon"], r["protein_accession"])
    log(f"primary atlas calls: {len(primary)}")

    # replicate codons: fetch CDS for the replicate proteins
    accs = sorted({r["protein_accession"] for r in rep_calls})
    log(f"\nretrieving CDS for {len(accs)} replicate proteins...")
    cds = {}
    B = 150
    for i in range(0, len(accs), B):
        chunk = accs[i:i + B]
        try:
            txt = efetch("protein", chunk, rettype="fasta_cds_na", retmode="text")
        except Exception as e:  # noqa: BLE001
            log(f"  batch {i//B+1}: FAILED ({type(e).__name__}) - those rows are dropped")
            continue
        for h, s in parse_fasta(txt):
            pid = None
            for tok in h.split("["):
                if tok.startswith("protein_id="):
                    pid = tok.split("=", 1)[1].rstrip("] ")
            if pid and s:
                cds[pid] = s
        if (i // B) % 4 == 0:
            log(f"  CDS {min(i+B, len(accs))}/{len(accs)}")
    log(f"  retrieved {len(cds)} CDS")

    rows = []
    for c in rep_calls:
        key = (c["species"], c["gene"], c["site_label"])
        if key not in primary:
            continue
        p_aa, p_codon, p_acc = primary[key]
        seq = cds.get(c["protein_accession"])
        r_codon = ""
        if seq:
            cand = normalise(seq[3 * (c["position"] - 1): 3 * c["position"]])
            if len(cand) == 3 and is_valid_codon(cand) and translate_codon(cand) == c["aa"]:
                r_codon = cand
        rows.append({
            "species": c["species"], "gene": c["gene"], "site_label": c["site_label"],
            "primary_assembly_protein": p_acc, "replicate_assembly": c["assembly"],
            "replicate_protein": c["protein_accession"],
            "primary_aa": p_aa, "replicate_aa": c["aa"],
            "aa_agree": p_aa == c["aa"],
            "primary_codon": p_codon, "replicate_codon": r_codon,
            "codon_agree": bool(r_codon) and r_codon == p_codon,
            "codon_comparable": bool(r_codon),
            "synonymous_difference": bool(r_codon) and r_codon != p_codon and p_aa == c["aa"],
            "replicate_context_identity": c["context_identity"],
        })

    out = os.path.join(DIR_RESULTS, "assembly_replication.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: (r["species"], r["gene"], r["site_label"])))
    log(f"\nwrote {out}  ({len(rows)} paired site calls)")

    n = len(rows)
    aa_ok = sum(1 for r in rows if r["aa_agree"])
    comp = [r for r in rows if r["codon_comparable"]]
    cd_ok = sum(1 for r in comp if r["codon_agree"])
    syn = sum(1 for r in comp if r["synonymous_difference"])
    log("\n=== REPLICATION ===")
    log(f"  species compared            : {len({r['species'] for r in rows})}")
    log(f"  paired site calls           : {n}")
    log(f"  amino-acid agreement        : {aa_ok}/{n} ({aa_ok/n:.2%})")
    log(f"  codon calls comparable      : {len(comp)}")
    log(f"  codon agreement             : {cd_ok}/{len(comp)} ({cd_ok/len(comp):.2%})")
    log(f"  synonymous differences      : {syn}/{len(comp)} ({syn/len(comp):.2%})"
        f"   <- within-species codon polymorphism")

    per_gene = defaultdict(lambda: [0, 0, 0, 0])
    for r in rows:
        g = per_gene[r["gene"]]
        g[0] += 1
        g[1] += r["aa_agree"]
        if r["codon_comparable"]:
            g[2] += 1
            g[3] += r["codon_agree"]
    log(f"\n  {'gene':10s} {'n':>5s} {'aa_agree':>9s} {'codon_agree':>12s}")
    for g, (nn, aa, cc, ca) in sorted(per_gene.items()):
        log(f"  {g:10s} {nn:5d} {aa/nn:9.2%} {(ca/cc if cc else float('nan')):12.2%}")

    summ = [{"metric": "species_compared", "value": len({r["species"] for r in rows})},
            {"metric": "paired_site_calls", "value": n},
            {"metric": "aa_agreement", "value": round(aa_ok / n, 4)},
            {"metric": "codon_calls_comparable", "value": len(comp)},
            {"metric": "codon_agreement", "value": round(cd_ok / len(comp), 4)},
            {"metric": "synonymous_difference_rate", "value": round(syn / len(comp), 4)}]
    sp = os.path.join(DIR_RESULTS, "assembly_replication_summary.csv")
    with open(sp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["metric", "value"])
        w.writeheader(); w.writerows(summ)
    log(f"\nwrote {sp}")
    if errs:
        log(f"\n{len(errs)} assemblies could not be processed (logged, not substituted)")
    log.close()


if __name__ == "__main__":
    main()
