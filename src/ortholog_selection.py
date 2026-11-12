#!/usr/bin/env python3
"""
ortholog_selection.py - ONE copy of the "which protein is this species' gene?" rule.

WHY THIS MODULE EXISTS
----------------------
The rule was written twice, and the two copies diverged. The correct form is:

  * a 40-deep k-mer prescreen, not top-6, because for Anopheles gambiae AChE all
    six slots were ace-2 isoforms and the true ace-1 ranked 7/8/10;
  * isoforms of one locus collapsed to one representative before truncation;
  * ONE protein per species x gene, not one per species x gene x SITE;
  * AChE candidates classified ace-1 / ace-2, ace-1 preferred where present;
  * identity x coverage decides, except inside a 5% band where mean site-context
    identity decides; a mean context below 0.30 is excluded outright.

The replicate side held an older copy (top-6 shortlist, per-site protein
choice, no paralogue classification), so the primary and replicate arms of the
replication analysis selected proteins by different rules.

The selection logic therefore lives here, once. Both arms import it.

WHAT IS AND IS NOT HERE
-----------------------
Here: the choice among already-aligned candidates (the code that lived inline in
step2_01's main(), moved verbatim).
Not here: the prescreen, alignment and site mapping - that is
step2_01_identify_orthologs.process_species(), which both arms also call.

VERIFIED
--------
A verification step runs this module over the committed
data/interim/site_mappings_all_candidates.csv and requires it to reproduce
data/interim/site_mappings.csv exactly, including the 483 ace-1 / 106 ace-2
split. A refactor that cannot reproduce the atlas it came from is not a
refactor.
"""
from __future__ import annotations

from collections import defaultdict

# Thresholds. These are the values step2_01 used for the committed atlas; they are
# defined here so the two arms cannot drift apart, and imported by step2_01 so
# there is still only one definition.
MIN_MEAN_CONTEXT = 0.30   # a candidate whose site contexts are poor is excluded
IDXCOV_TOL = 0.05         # tolerance band inside which context identity decides

_NUMERIC = ("identity", "coverage", "context_identity")


def coerce_map_rows(rows):
    """
    Make CSV-loaded mapping rows type-identical to in-memory ones.

    csv.DictReader yields strings, and the string "False" is TRUTHY. The original
    selection tests `if not m["mappable"]`, so feeding it raw CSV rows would keep
    every unmappable row and silently change the answer. That is the same class of
    error as the rest of this project's history - a value that looks like the thing
    standing in for the thing - so the coercion is explicit and it is loud.
    """
    out = []
    for r in rows:
        r = dict(r)
        m = r.get("mappable")
        if isinstance(m, str):
            if m not in ("True", "False"):
                raise ValueError(f"mappable is {m!r}; expected 'True' or 'False'")
            r["mappable"] = (m == "True")
        wt = r.get("matches_reference_wt")
        if isinstance(wt, str):
            r["matches_reference_wt"] = (wt == "True")
        for k in _NUMERIC:
            if k in r and isinstance(r[k], str):
                r[k] = float(r[k]) if r[k] != "" else 0.0
        out.append(r)
    return out


def _check_types(rows):
    """Refuse to run on rows that were not coerced. Fail loudly, not quietly."""
    for r in rows[:50]:
        if isinstance(r.get("mappable"), str):
            raise TypeError(
                "mappable is a string - call coerce_map_rows() first. "
                "'False' is truthy, so running on raw CSV rows would keep "
                "unmappable rows and change the selection."
            )
        for k in _NUMERIC:
            if k in r and isinstance(r[k], str):
                raise TypeError(f"{k} is a string - call coerce_map_rows() first.")


def select_one_protein_per_species_gene(all_maps, log=print):
    """
    Choose ONE protein per (species, gene) and return that protein's site rows.

    all_maps: candidate site-mapping rows (step2_01 process_species maps_out, or
              data/interim/site_mappings_all_candidates.csv through
              coerce_map_rows).

    Returns (selected_rows, stats). selected_rows is sorted by
    (species, gene, site_label), which is what step2_01 wrote.

    Ranking, unchanged from step2_01:
      * mean site-context identity below MIN_MEAN_CONTEXT -> excluded outright,
        so a globally similar but site-wrong paralogue cannot win on identity;
      * for AChE, ace-1 candidates only, where the species has one (cyclorrhaphan
        Diptera lost ace-1, so they fall through to ace-2, their real target);
      * identity x coverage decides, EXCEPT inside IDXCOV_TOL, where mean context
        identity decides. Global identity says which protein is the orthologue;
        local context says whether the sites are placed correctly; when the first
        is effectively tied the second must break it.
      * accession last, as a deterministic tie-break.
    """
    all_maps = list(all_maps)
    _check_types(all_maps)

    by_protein = defaultdict(list)
    for m in all_maps:
        if not m["mappable"]:
            continue
        by_protein[(m["species"], m["gene"], m["protein_accession"])].append(m)

    best_protein, ctx_rejected, n_ace1, n_ace2 = {}, 0, 0, 0
    ranked = defaultdict(list)
    for (sp, gene, pid), rows in by_protein.items():
        mean_ctx = sum(r["context_identity"] for r in rows) / len(rows)
        if mean_ctx < MIN_MEAN_CONTEXT:
            ctx_rejected += 1
            continue
        gscore = rows[0]["identity"] * rows[0]["coverage"]
        ranked[(sp, gene)].append((round(gscore, 6), round(mean_ctx, 6), pid,
                                   rows[0].get("paralogue", ""), rows))

    for key, cands in ranked.items():
        sp, gene = key
        if gene == "AChE":
            ace1 = [c for c in cands if c[3] == "ace1"]
            if ace1:
                cands, n_ace1 = ace1, n_ace1 + 1
            else:
                n_ace2 += 1
        best_id = max(c[0] for c in cands)
        eligible = [c for c in cands if c[0] >= best_id * (1.0 - IDXCOV_TOL)]
        chosen = max(eligible, key=lambda t: (t[1], t[0], t[2]))
        best_protein[key] = chosen[4]

    selected = sorted((r for rows in best_protein.values() for r in rows),
                      key=lambda r: (r["species"], r["gene"], r["site_label"]))

    stats = {
        "candidate_proteins": len(by_protein),
        "species_x_gene": len(best_protein),
        "alternatives_discarded": len(by_protein) - len(best_protein),
        "context_rejected": ctx_rejected,
        "ache_ace1": n_ace1,
        "ache_ace2": n_ace2,
        "selected_rows": len(selected),
        "non_rank1": sum(1 for m in selected if str(m.get("rank")) != "1"),
    }

    if log:
        log(f"pinned one protein per species x gene: {stats['species_x_gene']} species x gene; "
            f"{stats['alternatives_discarded']} alternative proteins discarded; "
            f"{stats['context_rejected']} candidates rejected on mean context "
            f"< {MIN_MEAN_CONTEXT}")
        log(f"  AChE paralogue assignment: {stats['ache_ace1']} species resolved to ace-1, "
            f"{stats['ache_ace2']} to ace-2 (no ace-1 present)")
        log(f"selected {stats['selected_rows']} site mappings "
            f"({stats['non_rank1']} took a non-rank-1 candidate on context identity)")

    return selected, stats
