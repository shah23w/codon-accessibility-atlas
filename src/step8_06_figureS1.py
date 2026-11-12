#!/usr/bin/env python3
"""
step8_06_figureS1.py: Supplementary Figure S1: what the atlas actually contains,
and how far its codon calls can be trusted.

Panel A  taxonomic breadth. Species per order, so a reader can see the panel is not
         four model insects.
Panel B  what is covered per gene: species with a codon call, and how many sites.
Panel C  the honesty check. Every call that could be replicated was re-derived from a
         SECOND, independently sequenced assembly. Agreement is shown per gene, at the
         amino-acid and at the codon level - including AChE, where paralogy makes the
         codon call least reliable.

Every number is read from a result file at run time; nothing is typed in.
"""
from __future__ import annotations

import csv
import os
import re
import sys
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atlas_common import DIR_FIGS, DIR_PROC, DIR_RESULTS, Tee, ensure_dirs  # noqa: E402
import fig00_style as st  # noqa: E402

TOP_ORDERS = 12


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def main():
    ensure_dirs()
    st.apply()
    log = Tee("step8_06_figureS1")
    log("=== Figure S1 ===")

    atlas = [r for r in read_csv(os.path.join(DIR_PROC, "codon_atlas.csv"))
             if r["is_wildtype"] == "True"]

    # ---- Panel A: species per order ---------------------------------------
    order_of = {r["species"]: r["order"] for r in atlas}
    orders = Counter(order_of.values())
    top = orders.most_common(TOP_ORDERS)
    rest_orders = len(orders) - len(top)
    rest_species = sum(orders.values()) - sum(n for _, n in top)
    n_species = len(order_of)
    log(f"  panel A: {n_species} species across {len(orders)} orders")

    # ---- Panel B: coverage per gene ---------------------------------------
    # GluCl_Px and GluCl_Tu are two SITES on one GluCl gene: ortho01 finds them
    # on the same protein accession in 236 of 569 species, at a near-constant
    # offset of +6 residues in 352. Drawing them as two bars would make this
    # panel contradict the paper's own gene count, so the labels are merged and
    # their site counts added, exactly as Figure 4 does.
    def gene_of(label):
        return "GluCl" if label.startswith("GluCl") else label

    gene_spp = defaultdict(set)
    gene_sites = defaultdict(set)
    gene_labels = set()
    for r in atlas:
        gene_labels.add(r["gene"])
        g = gene_of(r["gene"])
        gene_spp[g].add(r["species"])
        gene_sites[g].add(r["site_label"])
    genes = sorted(gene_spp, key=lambda g: -len(gene_spp[g]))

    # len(atlas) is the LONG-FORMAT ROW COUNT, which is a function of how many
    # alleles were curated at each site -- a property of the literature, not of
    # the atlas -- and canonical_quantities.csv warns explicitly against
    # quoting it as atlas size. It is 40,312; the atlas holds ONE codon per
    # species, gene and site, and there are 23,978 wild-type such calls, which
    # is what this panel reports.
    n_calls = len({(r["species"], r["gene"], r["site_label"]) for r in atlas})
    n_sites = sum(len(v) for v in gene_sites.values())
    log(f"  panel B: {n_calls} wild-type codon calls, {n_sites} sites, "
        f"{len(genes)} genes ({len(gene_labels)} gene labels)")

    # The panel asserts its own claims against the declared quantities and
    # aborts rather than drawing a number that has drifted from
    # canonical_quantities.csv.
    canon = {r["quantity"]: r["value"] for r in
             read_csv(os.path.join(DIR_RESULTS, "canonical_quantities.csv"))}

    def _declared(key):
        return int(float(canon[key]))

    for got, key, what in ((n_calls, "n_codon_calls_wildtype", "wild-type codon calls"),
                           (n_sites, "n_sites_atlas", "sites"),
                           (len(genes), "n_genes_atlas", "genes"),
                           (len(gene_labels), "n_gene_labels_atlas", "gene labels"),
                           (n_species, "n_species_atlas", "species")):
        assert got == _declared(key), (
            "panel B would misreport %s: computed %d, canonical_quantities.csv "
            "declares %d (%s)" % (what, got, _declared(key), key))

    # ---- Panel C: independent-assembly replication ------------------------
    # Pooled by GENE, using the same merge as panel B. This row is already an
    # aggregate over sites -- VGSC's pools 34 of them -- so leaving GluCl_Px and
    # GluCl_Tu as two rows would show six genes on the left of this figure and
    # seven on the right. Agreement is recomputed from the pooled counts, never
    # averaged across the two site rows.
    rep = defaultdict(lambda: [0, 0, 0])          # gene -> [n, aa_agree, codon_agree]
    for r in read_csv(os.path.join(DIR_RESULTS, "assembly_replication.csv")):
        v = rep[gene_of(r["gene"])]
        v[0] += 1
        v[1] += r["aa_agree"] == "True"
        v[2] += r["codon_agree"] == "True"
    summ = {r["metric"]: r["value"] for r in
            read_csv(os.path.join(DIR_RESULTS, "assembly_replication_summary.csv"))}
    rep_genes = sorted(rep, key=lambda g: -rep[g][0])
    log(f"  panel C: {summ['paired_site_calls']} paired calls in "
        f"{summ['species_compared']} species, aa {summ['aa_agreement']}, "
        f"codon {summ['codon_agreement']}")

    # -------------------------------------------------------------- figure
    fig = plt.figure(figsize=(st.W2, 3.15))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.00, 0.92, 1.06], wspace=0.62,
                          # 0.055 clipped the long order names once the
                          # figure was pinned to the 183 mm column:
                          # 'Hymenoptera' rendered as 'enoptera'.
                          left=0.105, right=0.985, top=0.86, bottom=0.17)

    # ---- A: taxonomic breadth ---------------------------------------------
    axA = fig.add_subplot(gs[0])
    labs = [o for o, _ in top] + [f"{rest_orders} further orders"]
    vals = [n for _, n in top] + [rest_species]
    cols = [st.C_OBS] * len(top) + [st.MUTED]
    yy = range(len(labs))
    axA.barh(list(yy), vals, height=0.68, color=cols, linewidth=0)
    axA.set_yticks(list(yy))
    axA.set_yticklabels(labs, fontsize=6.2)
    axA.invert_yaxis()
    axA.set_xlabel("Species", fontsize=6.8)
    axA.set_xlim(0, max(vals) * 1.20)
    axA.set_xticks([0, 50, 100, 150])
    st.recessive_grid(axA, "x")
    for i, v in enumerate(vals):
        axA.text(v + max(vals) * 0.022, i, str(v), va="center", fontsize=6.0,
                 color=st.INK if i < len(top) else st.INK2)
    st.panel_label(axA, "A", dx=-0.44, dy=1.10)

    # ---- B: coverage per gene, as a dot plot ------------------------------
    axB = fig.add_subplot(gs[1])
    vb = [len(gene_spp[g]) for g in genes]
    yy = list(range(len(genes)))
    for i, g in enumerate(genes):
        axB.plot([0, vb[i]], [i, i], color=st.GRID, lw=1.0, zorder=1)
        axB.plot(vb[i], i, "o", ms=5.4, color=st.C_OBS, mec=st.SURFACE, mew=0.9,
                 zorder=3)
    axB.set_yticks(yy)
    ns = [len(gene_sites[g]) for g in genes]
    axB.set_yticklabels([f"{g}\n{n} site{'s' if n > 1 else ''}"
                         for g, n in zip(genes, ns)], fontsize=6.0)
    axB.invert_yaxis()
    axB.set_xlabel("Species with a codon call", fontsize=6.8)
    axB.set_xlim(0, max(vb) * 1.22)
    st.recessive_grid(axB, "x")
    for i in yy:
        axB.text(vb[i] + max(vb) * 0.03, i, str(vb[i]), va="center", fontsize=6.0,
                 color=st.INK)
    st.panel_label(axB, "B", dx=-0.44, dy=1.10)

    # ---- C: replication, as a dumbbell ------------------------------------
    axC = fig.add_subplot(gs[2])
    order_g = sorted(rep, key=lambda g: 100 * rep[g][2] / rep[g][0])
    yy = list(range(len(order_g)))
    for i, g in enumerate(order_g):
        aa = 100 * rep[g][1] / rep[g][0]
        cd = 100 * rep[g][2] / rep[g][0]
        axC.plot([cd, aa], [i, i], color=st.GRID, lw=2.6, solid_capstyle="round",
                 zorder=1)
        # When the two agreements coincide -- RyR, where codon and amino acid
        # both replicate at 99.0% -- one marker lands exactly on the other and
        # the reader sees a single dot in a row of dumbbells, which looks like
        # a missing measurement rather than a perfect one. Draw the hidden one
        # as a ring around the visible one so both are legible at 5 pt.
        if abs(aa - cd) < 0.25:
            axC.plot(cd, i, "o", ms=8.4, mfc="none", mec=st.C_ALT, mew=1.1,
                     zorder=2)
            axC.plot(aa, i, "o", ms=5.0, color=st.C_OBS, mec=st.SURFACE,
                     mew=0.8, zorder=3)
        else:
            axC.plot(cd, i, "o", ms=5.0, color=st.C_ALT, mec=st.SURFACE, mew=0.8, zorder=3)
            axC.plot(aa, i, "o", ms=5.0, color=st.C_OBS, mec=st.SURFACE, mew=0.8, zorder=3)
    axC.set_yticks(yy)
    axC.set_yticklabels([f"{g}  (n = {rep[g][0]:,})" for g in order_g], fontsize=6.0)
    axC.invert_yaxis()
    axC.set_xlabel("Agreement with second assembly (%)", fontsize=6.8)
    # Panels a and b print their value beside every mark; this one did not,
    # and at VGSC and RyR the two dots are close enough that the gap cannot be
    # read off the axis. The pair is printed in a right-hand column instead,
    # the same idiom Fig. 2 uses, with the axis still ending at 100.
    for i, g in enumerate(order_g):
        axC.text(106.5, i, "%.1f / %.1f"
                 % (100 * rep[g][1] / rep[g][0], 100 * rep[g][2] / rep[g][0]),
                 va="center", ha="left", fontsize=5.8, color=st.INK,
                 clip_on=False)
    axC.text(106.5, -0.95, "AA / codon", va="center", ha="left",
             fontsize=5.8, color=st.MUTED, clip_on=False)
    # Room for that column inside the axes, since there is none outside it:
    # the limit runs past 100 but the ticks and the drawn spine stop there, so
    # the extra space reads as a margin rather than as axis.
    axC.set_xlim(60, 126)
    axC.set_xticks([60, 70, 80, 90, 100])
    axC.spines["bottom"].set_bounds(60, 100)
    st.recessive_grid(axC, "x")
    axC.legend(handles=[
        Line2D([], [], marker="o", ls="none", ms=5.0, color=st.C_OBS,
               label="Amino acid"),
        Line2D([], [], marker="o", ls="none", ms=5.0, color=st.C_ALT, label="Codon")],
        loc="lower left", bbox_to_anchor=(0.0, 1.00), ncol=2, handletextpad=0.3,
        borderpad=0.2, columnspacing=1.2, fontsize=6.2)
    st.panel_label(axC, "C", dx=-0.40, dy=1.10)

    _sd = [("a", _o, int(_n), "", "", "") for _o, _n in top]
    _sd.append(("a", "further orders pooled", int(sum(orders.values())
                                                  - sum(n for _, n in top)),
                "", "", "%d orders" % rest_orders))
    for _g in genes:
        _sd.append(("b", _g, len(gene_spp[_g]), "", "",
                    "%d sites" % len(gene_sites[_g])))
    for _g in rep:
        _sd.append(("c", "%s amino acid" % _g,
                    100.0 * rep[_g][1] / rep[_g][0], "", "",
                    "n = %d" % rep[_g][0]))
        _sd.append(("c", "%s codon" % _g,
                    100.0 * rep[_g][2] / rep[_g][0], "", "",
                    "n = %d" % rep[_g][0]))
    st.write_source_data("FigS1", _sd)

    st.write_legend_data("figureS1_atlas_composition", {
        "n_species": n_species,
        "n_codon_calls": n_calls,
        "n_sites": n_sites,
        "orders_shown": {o: int(n) for o, n in top},
        "n_orders_total": len(orders),
        "species_with_call_by_gene": {g: len(gene_spp[g]) for g in genes},
        "sites_by_gene": {g: len(gene_sites[g]) for g in genes},
        "replication": {g: {"n": rep[g][0],
                            "amino_acid_pct": 100.0 * rep[g][1] / rep[g][0],
                            "codon_pct": 100.0 * rep[g][2] / rep[g][0]}
                        for g in rep},
    })

    for ext in ("pdf", "png"):
        p = os.path.join(DIR_FIGS, f"figureS1_atlas_composition.{ext}")
        fig.savefig(p)
        if ext == "pdf":
            st.assert_no_type3(p)
            st.audit(fig, p)
        log(f"  wrote {p}")
    plt.close(fig)
    log("Figure S1 complete.")


if __name__ == "__main__":
    main()
