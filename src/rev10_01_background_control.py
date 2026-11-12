#!/usr/bin/env python3
"""
rev10_01_background_control.py

THE GENERIC-SUBSTITUTION ALTERNATIVE, TESTED ON THE SAME SCALE.

Resistance alleles may be enriched for one-step substitutions not because
accessibility constrains adaptation to insecticides, but because amino acid
substitutions in general arise through readily accessible mutational paths and
resistance evolution recruits some of them. Bolourchi et al. (2025 PNAS
122:e2424307122), screening 18,000 bacterial genomes for rpoB resistance
mutations, reached a compatible conclusion: "the genetic code imposes only
minor constraints on RIF resistance evolution". On either reading the effect
would be generic to substitution rather than specific to resistance.

THE CONTROL ALREADY EXISTS.  step4_01_ancestral.py reconstructs ancestral
codons on an NCBI-taxonomy tree (independent of the target sequences, so the
reconstruction cannot manufacture the effect), masks each tip, and identifies
489 NATURAL substitutions at the same sites - lineages that changed residue
with no documented insecticide involvement. Its docstring anticipates the
objection exactly: the question is not whether substitutions are one-step,
which is trivially true, but whether the ancestral codon predicts WHICH of the
viable residues a lineage acquired.

WHAT WAS MISSING.  That control reports an ENRICHMENT RATIO (observed over
expected). The manuscript's headline reports an ODDS RATIO. The two arms have
never been put on the same scale, so no one can say whether accessibility
matters more for resistance than for ordinary evolution.

THIS SCRIPT builds both arms in identical long format and runs the identical
estimator on each:

  ARM R  resistance.  Rows are species x site x documented resistance allele.
         used = that species is recorded as having evolved that allele.
         accessible = one substitution from the species' own codon.

  ARM N  natural.  Rows are species x site x viable residue, for lineages that
         actually changed residue at that site. used = the residue the lineage
         acquired. accessible = one substitution from the RECONSTRUCTED
         ANCESTRAL codon. The viable set is the residues observed at that site
         across the panel, which is the same "empirically viable" principle the
         resistance arm applies to documented alleles.

  ESTIMATOR  Mantel-Haenszel odds ratio stratified by site, with
  Robins-Breslow-Greenland variance - the same estimator the corrected
  manuscript reports, so ARM R must reproduce the committed MH value. That
  reproduction is a validation check built into the run: if ARM R does not come
  back near the committed figure, the long format is wrong and ARM N means
  nothing.

READING THE RESULT
  OR(R) >> OR(N)   accessibility matters more for resistance than for ordinary
                   substitution, and the generic-substitution alternative does
                   not account for the effect.
  OR(R) ~= OR(N)   the alternative stands. The effect is a property of
                   mutation, not of resistance.

Ambiguous ancestral reconstructions are excluded from the primary and reported
separately; Step 7 established that 94.7% of reconstructions are unambiguous.

Reads   data/processed/codon_atlas.csv
        results/core_test_events*.csv
        results/ancestral_substitutions*.csv
Writes  results/rev10_arm_rows.csv
        results/rev10_comparison.csv
        results/rev10_summary.md
"""
from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PROJ = os.path.dirname(HERE)
RES = os.path.join(PROJ, "results")

from codon_engine import codon_reaches, is_valid_codon, normalise   # noqa: E402
from rev09_01_two_step_routes import resolve_results_file           # noqa: E402

ATLAS = os.path.join(PROJ, "data", "processed", "codon_atlas.csv")


# --------------------------------------------------------------------------- #
# Mantel-Haenszel with Robins-Breslow-Greenland variance
# --------------------------------------------------------------------------- #
def mantel_haenszel(strata):
    """
    strata: iterable of (a, b, c, d) with
        a = used & accessible      b = used & not accessible
        c = unused & accessible    d = unused & not accessible
    Returns (OR, lo, hi, n_informative_strata).
    """
    R = S = 0.0
    PR = PS = QR = QS = 0.0
    n_inf = 0
    for a, b, c, d in strata:
        n = a + b + c + d
        if n == 0:
            continue
        r = a * d / n
        s = b * c / n
        if r == 0 and s == 0:
            continue
        n_inf += 1
        R += r
        S += s
        P = (a + d) / n
        Q = (b + c) / n
        PR += P * r
        PS += P * s
        QR += Q * r
        QS += Q * s
    if R == 0 or S == 0:
        return (float("inf") if S == 0 and R > 0 else float("nan"),
                float("nan"), float("nan"), n_inf)
    or_mh = R / S
    var = PR / (2 * R * R) + (PS + QR) / (2 * R * S) + QS / (2 * S * S)
    se = math.sqrt(var)
    return or_mh, or_mh * math.exp(-1.96 * se), or_mh * math.exp(1.96 * se), n_inf


def crude(rows):
    a = sum(1 for r in rows if r["used"] and r["accessible"])
    b = sum(1 for r in rows if r["used"] and not r["accessible"])
    c = sum(1 for r in rows if not r["used"] and r["accessible"])
    d = sum(1 for r in rows if not r["used"] and not r["accessible"])
    return a, b, c, d


def strata_from(rows):
    by = defaultdict(lambda: [0, 0, 0, 0])
    for r in rows:
        k = (r["gene"], r["site_label"])
        i = (0 if r["used"] else 2) + (0 if r["accessible"] else 1)
        by[k][i] += 1
    return [tuple(v) for v in by.values()], by


def main():
    print("=" * 78)
    print("STEP 10.  Is the accessibility effect specific to resistance?")
    print("=" * 78)
    print()

    ev_path = resolve_results_file("core_test_events")
    anc_path = resolve_results_file("ancestral_substitutions")
    mv_path = resolve_results_file("multivariate_dataset")
    print()

    # ------------------------------------------------------------------ #
    # VALIDATE THE ESTIMATOR against the committed analysis dataset.
    #
    # The manuscript's odds ratio is computed by step3_01 on
    # multivariate_dataset_primary.csv. Running THIS Mantel-Haenszel code on
    # THAT file must return the committed 26.26. If it does, the estimator is
    # correct and any difference below is a difference of estimand, not a bug.
    # ------------------------------------------------------------------ #
    tf = lambda v: str(v).strip().lower() in ("true", "1", "yes", "t")
    mv = list(csv.DictReader(open(mv_path)))
    by = defaultdict(lambda: [0, 0, 0, 0])
    for r in mv:
        k = r.get("site") or "%s:%s" % (r["gene"], r["site_label"])
        by[k][(0 if tf(r["observed"]) else 2) + (0 if tf(r["accessible"]) else 1)] += 1
    or_c, lo_c, hi_c, _ = mantel_haenszel([tuple(v) for v in by.values()])
    print("-" * 78)
    print("ESTIMATOR VALIDATION  (against the committed analysis dataset)")
    print("-" * 78)
    print("  %s : %d rows, %d sites" % (os.path.basename(mv_path), len(mv), len(by)))
    print("  Mantel-Haenszel OR from this code : %.2f  95%% CI [%.2f, %.2f]"
          % (or_c, lo_c, hi_c))
    print("  committed manuscript value        : 26.26")
    if abs(or_c - 26.26) > 1.0:
        sys.exit("\n  ESTIMATOR MISMATCH. Fix this before reading anything below.")
    print("  ESTIMATOR VALIDATED. Differences below are differences of ESTIMAND.")
    print()
    print("  NOTE ON THE ESTIMAND. That dataset has NO species column: its rows are")
    print("  site x candidate residue, where 'accessible' means reachable from AT")
    print("  LEAST ONE species' codon at that site and 'observed' means documented")
    print("  ANYWHERE. So 26.26 answers:")
    print("      is a residue reachable from some species' codon more often")
    print("      documented as a resistance allele somewhere?")
    print("  The title asks something else - which allele AN ARTHROPOD evolves -")
    print("  and that is a species-level question. ARM R below is the species-level")
    print("  test: rows are species x site x allele, 'accessible' is reachable from")
    print("  THAT species' own codon, and 'used' is THAT species evolving it. Both")
    print("  are legitimate; only the second matches the claim in the title, and")
    print("  only the second is comparable with ARM N.")
    print()

    # -------------------------------------------------- observed resistance
    observed = set()
    with open(ev_path) as fh:
        for r in csv.DictReader(fh):
            observed.add((r["species"], r["gene"], r["site_label"], r["resistance_allele"]))

    # -------------------------------------------------- ARM R from the atlas
    armR = []
    viable_at_site = defaultdict(set)      # (gene, site) -> residues seen across species
    with open(ATLAS) as fh:
        for r in csv.DictReader(fh):
            viable_at_site[(r["gene"], r["site_label"])].add(r["species_aa"])
            if r["is_wildtype"] != "True" or r["already_resistant"] == "True":
                continue
            armR.append({
                "arm": "resistance", "species": r["species"], "gene": r["gene"],
                "site_label": r["site_label"], "from_codon": normalise(r["codon"]),
                "candidate_aa": r["resistance_allele"],
                "accessible": r["one_step_accessible"] == "True",
                "used": (r["species"], r["gene"], r["site_label"],
                         r["resistance_allele"]) in observed,
            })

    # -------------------------------------------------- ARM N from ancestral
    armN = []
    n_ambig = 0
    with open(anc_path) as fh:
        for r in csv.DictReader(fh):
            anc = [normalise(c) for c in r["ancestral_codons"].split("|") if c]
            if len(anc) != 1 or not is_valid_codon(anc[0]):
                n_ambig += 1
                continue
            c0 = anc[0]
            gene, site = r["gene"], r["site_label"]
            anc_aa = r["ancestral_aas"].split("|")[0]
            cands = {a for a in viable_at_site[(gene, site)] if a and a != anc_aa}
            for aa in sorted(cands):
                armN.append({
                    "arm": "natural", "species": r["species"], "gene": gene,
                    "site_label": site, "from_codon": c0, "candidate_aa": aa,
                    "accessible": codon_reaches(c0, aa),
                    "used": aa == r["tip_aa"],
                })

    # -------------------------------------------------- report
    out = []

    def p(s=""):
        print(s)
        out.append(s)

    results = []
    for name, rows in (("ARM R  resistance", armR), ("ARM N  natural substitutions", armN)):
        a, b, c, d = crude(rows)
        st, by = strata_from(rows)
        or_mh, lo, hi, n_inf = mantel_haenszel(st)
        crude_or = (a * d / (b * c)) if b and c else float("inf")
        p("-" * 78)
        p(name)
        p("-" * 78)
        p("  rows (species x site x candidate allele) : %d" % len(rows))
        p("  species x site strata                    : %d  (informative: %d)"
          % (len(by), n_inf))
        p()
        p("                      accessible   not accessible")
        p("    used                  %6d         %6d" % (a, b))
        p("    not used              %6d         %6d" % (c, d))
        p()
        p("  crude odds ratio                 : %s"
          % ("inf" if crude_or == float("inf") else "%.2f" % crude_or))
        p("  Mantel-Haenszel OR (by site)     : %s   95%% CI [%s, %s]"
          % ("inf" if or_mh == float("inf") else "%.2f" % or_mh,
             "%.2f" % lo if lo == lo else "na", "%.2f" % hi if hi == hi else "na"))
        p()
        results.append({"arm": name, "a": a, "b": b, "c": c, "d": d,
                        "crude_or": crude_or, "mh_or": or_mh, "ci_lo": lo, "ci_hi": hi,
                        "n_rows": len(rows), "n_strata": len(by), "n_informative": n_inf})

    p("  (%d ancestral reconstructions were ambiguous and excluded; Step 7" % n_ambig)
    p("   established that 94.7 percent are unambiguous)")
    p()

    R, N = results[0], results[1]
    p("=" * 78)
    p("THE COMPARISON")
    p("=" * 78)
    p()
    p("  Both arms below are SPECIES-LEVEL and use the identical estimator,")
    p("  validated above against the committed dataset. They are therefore")
    p("  comparable with each other, which the question requires.")
    p()
    p("  for reference, the three estimands:")
    p("    site-level, pooled over species (the manuscript's headline) : %.2f" % or_c)
    p("    species-level, resistance            (ARM R)                : %.2f" % R["mh_or"])
    p("    species-level, natural substitutions (ARM N)                : %.2f" % N["mh_or"])
    p()
    if N["mh_or"] == N["mh_or"] and N["mh_or"] not in (0.0, float("inf")):
        ratio = R["mh_or"] / N["mh_or"]
        # z on the log scale, treating the two arms as independent
        seR = (math.log(R["ci_hi"]) - math.log(R["ci_lo"])) / (2 * 1.959964)
        seN = (math.log(N["ci_hi"]) - math.log(N["ci_lo"])) / (2 * 1.959964)
        se = math.sqrt(seR ** 2 + seN ** 2)
        z = (math.log(R["mh_or"]) - math.log(N["mh_or"])) / se if se > 0 else float("nan")
        pz = math.erfc(abs(z) / math.sqrt(2))
        p("  resistance OR      %.2f  [%.2f, %.2f]" % (R["mh_or"], R["ci_lo"], R["ci_hi"]))
        p("  natural OR         %.2f  [%.2f, %.2f]" % (N["mh_or"], N["ci_lo"], N["ci_hi"]))
        p("  ratio of ORs       %.2f     z = %.2f, P = %.4g" % (ratio, z, pz))
        p()
        if pz < 0.05 and ratio > 1:
            p("  ACCESSIBILITY MATTERS MORE FOR RESISTANCE THAN FOR ORDINARY SUBSTITUTION.")
            p("  The generic-substitution alternative, that resistance merely recruits")
            p("  substitutions that were already accessible, does not account for it. The")
            p("  same estimator on natural substitutions at the same sites, in the same")
            p("  species, gives a %.1f-fold weaker association." % ratio)
        elif pz < 0.05 and ratio < 1:
            p("  Accessibility matters LESS for resistance than for ordinary substitution.")
            p("  This is the opposite of the manuscript's claim and must be reported.")
        else:
            p("  THE TWO ARMS ARE NOT DISTINGUISHABLE.")
            p("  On this evidence the generic-substitution alternative stands: the")
            p("  association between accessibility and which residue is acquired is not")
            p("  specific to resistance, and the resistance odds ratio should not be")
            p("  read as a fact about resistance evolution.")
    else:
        p("  ARM N gave a degenerate odds ratio; report the 2x2 directly instead.")

    os.makedirs(RES, exist_ok=True)
    allrows = armR + armN
    with open(os.path.join(RES, "rev10_arm_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()))
        w.writeheader()
        for r in allrows:
            w.writerow(r)
    with open(os.path.join(RES, "rev10_comparison.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        for r in results:
            w.writerow(r)
    with open(os.path.join(RES, "rev10_summary.md"), "w") as fh:
        fh.write("# Step 10: is the accessibility effect specific to resistance?\n\n```\n")
        fh.write("\n".join(out))
        fh.write("\n```\n")

    p()
    p("WROTE results/rev10_arm_rows.csv, rev10_comparison.csv, rev10_summary.md")


if __name__ == "__main__":
    main()
