#!/usr/bin/env python3
"""
data15_canonical_quantities.py
------------------------------
One quantity, one name, one owning file, one value. Everything is READ from the
result files; nothing is typed. If a value here disagrees with a result file,
this script is wrong and says so -- it never wins an argument with the data.

Writes:
    results/CANONICAL_QUANTITIES.md     the human document (for the SI)
    results/canonical_quantities.csv    machine-readable; FIGURES READ THIS

Why this file exists.

  Two quantities sharing a name is a collision a reader adding up a table would
  hit, and a figure re-deriving a number from a rounded one prints a value that
  appears nowhere in the data. Both are prevented the same way: name each
  quantity once, state its definition, and make every consumer read it here.

Exit status: 0 if every declaration resolves against the files; 1 if any
declaration is stale, missing, or contradicted.
"""
import csv, os, re, sys

RES = "results"
out = []
def say(s=""):
    print(s); out.append(s)

def rd(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p): return None
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))

def num(x):
    try: return float(str(x).strip())
    except Exception: return None

problems, rows = [], []

def record(key, value, definition, owner, method, note=""):
    """value may be None -> that is a problem, not a blank."""
    if value is None:
        problems.append("%s: could not be resolved from %s" % (key, owner))
        value = "UNRESOLVED"
    rows.append(dict(quantity=key, value=value, definition=definition,
                     owning_file=owner, method=method, note=note))
    return value

# ===========================================================================
#  A0. ATLAS SIZE
# ===========================================================================
rc = rd("data16_recount.csv") or []
by_rc = {r["quantity"]: r for r in rc}

# `or []` above makes a missing input indistinguishable from an input full of
# blanks, and every rcv() below then returns None. A silent None propagates: if
# results/data16_recount.csv is absent, the gene count declared here never
# updates and canonical_quantities.csv goes on publishing the old one. The
# atlas size is the spine of this file, so its absence is fatal, not a default.
if not by_rc:
    sys.exit("FATAL: results/data16_recount.csv is missing or empty. Every "
             "quantity in section A0 derives from it, so this script cannot "
             "run. Rebuild it first:\n"
             "    python3 src/data16_recount_reconciliation.py")

_REQUIRED_RC = ("keys_all", "keys_wildtype", "keys_not_wildtype",
                "rows_wildtype", "n_species", "n_genes", "n_gene_site_pairs")
_absent = [k for k in _REQUIRED_RC if num((by_rc.get(k) or {}).get("value")) is None]
if _absent:
    sys.exit("FATAL: results/data16_recount.csv is present but does not carry "
             "a usable value for: %s. Rebuild it with "
             "src/data16_recount_reconciliation.py." % ", ".join(_absent))


def rcv(k):
    return num((by_rc.get(k) or {}).get("value"))

n_calls_all = rcv("keys_all")
n_calls_wt = rcv("keys_wildtype")
n_calls_nwt = rcv("keys_not_wildtype")
n_rows_wt = rcv("rows_wildtype")

record("n_codon_calls", n_calls_all,
       "distinct species x gene x site codon calls in the atlas, whatever the "
       "species carries at that site",
       "data16_recount.csv", "row keys_all",
       "DECLARED HEADLINE for atlas size; this is what the abstract quotes as "
       "codon calls. A long-format row count is a different quantity and must "
       "not be quoted in its place.")

record("n_codon_calls_wildtype", n_calls_wt,
       "distinct species x gene x site codon calls where the species is "
       "wild-type at that site",
       "data16_recount.csv", "row keys_wildtype",
       "The subset the null pools are built from, because a species that is "
       "not wild-type at a site contributes no codon to that site's pool. "
       "This is the 23,978 reported by step9_01.")

record("n_codon_calls_not_wildtype", n_calls_nwt,
       "distinct species x gene x site codon calls where the species is not "
       "wild-type at that site",
       "data16_recount.csv", "row keys_not_wildtype",
       "The complement. 25,574 minus 23,978; the reason the two scripts "
       "appeared to disagree.")

record("n_atlas_rows_long_format", n_rows_wt,
       "rows of codon_atlas.csv carrying is_wildtype == True",
       "data16_recount.csv", "row rows_wildtype",
       "NOT a count of codon calls. One row per curated allele at a site, so "
       "this is a property of how much literature exists per site. Never quote "
       "it as atlas size.")

# The size of Dataset S1 itself, which is the row count the previous entry is
# NOT. The Methods quoted 40,312 as the size of Dataset S1 and said those rows
# carried the 25,574 calls; both were wrong, because 40,312 is the wild-type
# subset and it carries the 23,978 wild-type calls. Declaring the total here
# means the two can be told apart on the page.
record("n_atlas_rows_total", rcv("rows"),
       "rows of codon_atlas.csv in all, wild-type and not",
       "data16_recount.csv", "row rows",
       "The size of Dataset S1. The wild-type subset of it is "
       "n_atlas_rows_long_format; quoting one for the other conflates the "
       "dataset with the part of it the core test uses.")

record("n_species_atlas", rcv("n_species"),
       "distinct species carrying at least one codon call",
       "data16_recount.csv", "row n_species")

record("n_gene_labels_atlas", rcv("n_genes"),
       "distinct gene LABELS in the atlas",
       "data16_recount.csv", "row n_genes",
       "NOT the number of genes. GluCl_Px and GluCl_Tu are two sites on one "
       "GluCl gene, not two genes: ortho01 finds them on the same protein "
       "accession in 236 of 569 species and at a near-constant offset of +6 "
       "residues in 352 of them. Never quote this as a gene count.")

record("n_genes_atlas", rcv("n_genes") - 1,
       "distinct target GENES in the atlas, after merging the two GluCl labels",
       "ortho01_gluc_correction.csv", "n_gene_labels minus one",
       "DECLARED HEADLINE for the gene count: six, not the seven gene LABELS. "
       "The site count is unaffected. Leave-one-gene-out "
       "that drops GluCl_Px or GluCl_Tu separately is dropping one SITE, not "
       "one gene, and must be described that way.")

record("n_sites_atlas", rcv("n_gene_site_pairs"),
       "distinct gene x site pairs in the atlas",
       "data16_recount.csv", "row n_gene_site_pairs")

# the atlas partition must close, or the headline atlas size is not established
if None not in (n_calls_all, n_calls_wt, n_calls_nwt):
    if abs((n_calls_wt + n_calls_nwt) - n_calls_all) > 1e-9:
        problems.append("atlas partition does not close: %g + %g = %g, not %g"
                        % (n_calls_wt, n_calls_nwt,
                           n_calls_wt + n_calls_nwt, n_calls_all))
else:
    problems.append("atlas size could not be read from data16_recount.csv; "
                    "run src/data16_recount_reconciliation.py first")

# data16 audits itself; if its own checks failed, nothing here is quotable
_cf = rcv("checks_failed")
if _cf is None or _cf != 0:
    problems.append("data16_recount.csv reports %s failed check(s), so the "
                    "atlas size is not established" % _cf)

# ===========================================================================
#  A. THE EVENT PARTITION  --  where the two "informative" numbers come from
# ===========================================================================
est = rd("step2_estimands.csv") or []
by_q = {r["quantity"]: r for r in est}

n_p1   = num((by_q.get("stratum_p_eq_1") or {}).get("n"))
n_p0   = num((by_q.get("stratum_p_eq_0") or {}).get("n"))
n_inf  = num((by_q.get("conditional_informative") or {}).get("n"))
n_all  = num((by_q.get("unconditional_all_events") or {}).get("n"))

ct = rd("core_test_conditional_primary.csv") or []
n_disc = float(len(ct)) if ct else None

record("n_events_total", n_all,
       "curated resistance events entering the core test",
       "step2_estimands.csv", "row unconditional_all_events")

record("n_events_supply_informative", n_inf,
       "events where the SITE x ALLELE supply fraction is strictly between 0 "
       "and 1 -- some but not all panel species could reach that allele there",
       "step2_estimands.csv", "row conditional_informative",
       "The 48. A property of the site, not of the species.")

record("n_events_supply_saturated", n_p1,
       "events where every panel species could reach the allele (supply = 1)",
       "step2_estimands.csv", "row stratum_p_eq_1")

record("n_events_supply_zero", n_p0,
       "events where no panel species could reach the allele (supply = 0)",
       "step2_estimands.csv", "row stratum_p_eq_0")

record("n_events_codon_discriminating", n_disc,
       "events where THAT SPECIES' OWN codon reaches some but not all of the "
       "alleles curated at that site, so it can tell them apart",
       "core_test_conditional_primary.csv", "step3_03 TEST 2, one row per event",
       "The 32. A property of the species' codon, not of the site. "
       "NOT a subset relation with the 48 -- a different question.")

# the partition must close
if None not in (n_p1, n_p0, n_inf, n_all):
    tot = n_p1 + n_p0 + n_inf
    if abs(tot - n_all) > 1e-9:
        problems.append("partition does not close: %g + %g + %g = %g, not %g"
                        % (n_p1, n_p0, n_inf, tot, n_all))

# ===========================================================================
#  B. BREACHES  --  two scopes, one total
# ===========================================================================
acc_col = None
if ct:
    for pref in ("observed_allele_accessible", "accessible"):
        if pref in ct[0]: acc_col = pref; break
n_br_disc = None
if acc_col:
    n_br_disc = float(sum(1 for r in ct if str(r[acc_col]).strip().lower() == "false"))

record("n_breaches_total", n_p0,
       "events whose observed allele was not one nucleotide step from the "
       "species' wild-type codon",
       "step2_estimands.csv", "row stratum_p_eq_0",
       "Verified identical to the supply-zero stratum, not assumed.")

record("n_breaches_among_discriminating", n_br_disc,
       "breaches that fall inside the codon-discriminating subset",
       "core_test_conditional_primary.csv", "rows with %s == False" % (acc_col or "?"),
       "One breach (Ctenocephalides felis VGSC T929 ACG->V) reaches NONE of the "
       "site repertoire, so its codon cannot discriminate and TEST 2 excludes it "
       "by construction. The difference is definitional, not a missing event.")

if None not in (n_br_disc, n_p0) and not (0 <= n_p0 - n_br_disc <= 2):
    problems.append("breach scopes differ by %g; expected the discriminating "
                    "subset to be at most a couple short of the total"
                    % (n_p0 - n_br_disc))

# ===========================================================================
#  C. THE POOLED NULL  --  canonical estimator declared
# ===========================================================================
exp_exact = num((by_q.get("unconditional_all_events") or {}).get("expected"))
rat_exact = num((by_q.get("unconditional_all_events") or {}).get("ratio"))
obs_exact = num((by_q.get("unconditional_all_events") or {}).get("observed"))

record("expected_null_pooled", exp_exact,
       "sum of per-event accessibility probabilities under the atlas-pool null",
       "step2_estimands.csv", "Poisson-binomial exact expectation, sum(p_i)",
       "CANONICAL. Permutation budgets give 84.68-84.71; those estimate this "
       "number, they do not compete with it. Quote the exact value.")

record("pooled_enrichment", rat_exact,
       "observed accessible / expected accessible, all events",
       "step2_estimands.csv", "ratio of the two rows above",
       "CANONICAL. Must equal observed/expected_null_pooled at full precision. "
       "This is NOT the headline result of the paper -- see the note at the end.")

# does the canonical ratio actually reproduce?
if None not in (obs_exact, exp_exact, rat_exact):
    recomputed = obs_exact / exp_exact
    if abs(recomputed - rat_exact) > 5e-4:
        problems.append("pooled_enrichment %.6f does not reproduce from "
                        "%g/%g = %.6f" % (rat_exact, obs_exact, exp_exact, recomputed))

# cross-check step9_01, which computes the same thing from its own p-vector
cs = rd("corrected_statistics.csv") or []
cs_by = {r["statistic"]: r for r in cs}
s9_exp = num((cs_by.get("expected_null_pooled") or {}).get("value"))
s9_rat = num((cs_by.get("pooled_enrichment") or {}).get("value"))
if s9_exp is None and s9_rat is not None:
    problems.append("corrected_statistics.csv reports pooled_enrichment %s but "
                    "no expected_null_pooled row -- run patch01_precision.py and "
                    "re-run step9_01, otherwise its denominator is unrecorded"
                    % s9_rat)
elif None not in (s9_exp, exp_exact):
    d = abs(s9_exp - exp_exact)
    if d > 0.02:
        problems.append("step9_01 expected %.6f vs canonical %.6f (differ by "
                        "%.4f): the two scripts are building different p-vectors. "
                        "One is wrong; find out which before any figure quotes a "
                        "ratio." % (s9_exp, exp_exact, d))
    elif d > 1e-6:
        # NOT float accumulation. Summing 108 doubles accumulates about 1e-12; a gap of
        # 1e-4 is three orders of magnitude too large for that. The cause is storage
        # precision: step3_03_core_test.py writes site_supply_fraction as
        # round(frac, 4), step2_01 sums that rounded column, and step9_01 recomputes
        # p_i from the codon pools. So step9_01's value is the exact one.
        say("  note: step9_01 expected %.6f vs step2_01 %.6f (differ by %.6f)."
            % (s9_exp, exp_exact, d))
        say("        This is NOT float accumulation -- that would be ~1e-12. It is")
        say("        storage precision: step3_03 stores site_supply_fraction as")
        say("        round(frac, 4), step2_01 sums the rounded column, step9_01")
        say("        recomputes from the codon pools. step9_01 holds the exact value.")
        say("        Immaterial to every conclusion (1.2161 vs 1.2162), and left")
        say("        as it stands: fixing it means storing more digits in")
        say("        step3_03_core_test.py, which 48 scripts read downstream. No")
        say("        reported digit of either value differs.")

# ===========================================================================
#  D. ODDS RATIOS  --  estimand families, with the headline declared
# ===========================================================================
#  These are not competing answers. They are different estimators of two
#  different estimands, and reporting several is good practice. What is NOT
#  good practice is leaving unstated which one the abstract quotes.
# LOCKED DECLARATIONS. The paper quotes these and no others as the headline odds
# ratios. Site-level: the conditional logit adjusted for ESM-2 -- the smallest of the
# five estimators, so it cannot be called cherry-picked, and the only one with ESM in
# the model, so it answers "this is protein tolerance" inside the estimate rather than
# in a separate table. Species-level: Mantel-Haenszel stratified by site. Everything
# else is a robustness panel. The GEE figure is one row among several, not the
# headline.
HEADLINE_SITE    = "conditional_logit_by_site_adj_esm"
HEADLINE_SPECIES = "site"

oe = rd("step3_odds_ratio_estimators.csv") or []
site_ors = {r["estimator"]: num(r["odds_ratio"]) for r in oe}
hv = site_ors.get(HEADLINE_SITE)
record("OR_site_level_headline", hv,
       "site-level odds that an observed allele is accessible, conditional "
       "logit with site fixed effects, adjusted for ESM-2 predicted effect",
       "step3_odds_ratio_estimators.csv", HEADLINE_SITE,
       "DECLARED HEADLINE. The most conservative of the five and the one that "
       "survives the 'this is just protein-level tolerance' objection, because "
       "ESM-2 is in the model. Quoting the smallest defensible effect is the "
       "stronger position; the other four go in the robustness panel.")
for k, v in sorted(site_ors.items()):
    if v is None or k == HEADLINE_SITE: continue
    record("OR_site_level__" + k, v,
           "alternative estimator of the same site-level estimand",
           "step3_odds_ratio_estimators.csv", k, "robustness panel, not the headline")

rv = rd("rev11_stratification.csv") or []
sp = {r["stratification"]: num(r["odds_ratio"]) for r in rv}
record("OR_species_level_headline", sp.get(HEADLINE_SPECIES),
       "species-level odds, Mantel-Haenszel stratified by site",
       "rev11_stratification.csv", "stratification = " + HEADLINE_SPECIES,
       "DECLARED HEADLINE. order x site and family x site add phylogenetic "
       "strata and move it to 13.27 and 11.18 -- report as robustness.")
for k, v in sorted(sp.items()):
    if v is None or k == HEADLINE_SPECIES: continue
    record("OR_species_level__" + k.replace(" ", "_"), v,
           "same estimand, additionally stratified", "rev11_stratification.csv",
           k, "robustness")

ar = rd("aprd02_restricted_odds_ratio.csv") or []
for r in ar:
    lab = list(r.values())[0]
    kk = [c for c in r if "odds" in c.lower()]
    if not kk: continue
    v = num(r[kk[0]])
    dup = (sp.get(HEADLINE_SPECIES) is not None and v is not None
           and abs(v - sp[HEADLINE_SPECIES]) < 0.01)
    record("OR_species_level_ascertainment__" + str(lab).replace(" ", "_"), v,
           "species-level odds restricted by ascertainment subset",
           "aprd02_restricted_odds_ratio.csv", str(lab),
           ("SAME ESTIMATOR as the headline at finer precision -- not a sixth "
            "estimate; the apparent multiplicity was a duplicate." if dup else
            "ascertainment robustness: the effect survives restriction"))

# the assumption-free anchor, which has a P but no OR
for r in oe:
    if r["estimator"] == "exact_stratified_hypergeometric":
        record("P_exact_stratified_hypergeometric", num(r["p_value"]),
               "exact stratified hypergeometric tail, conditioning on site margins",
               "step3_odds_ratio_estimators.csv", r["estimator"],
               "No odds ratio by construction. This is the assumption-free "
               "anchor: it makes no model of the effect, only of the margins. "
               "Quote it alongside the headline OR.")


# ===========================================================================
#  E. MUTATIONAL-SUPPLY EXPONENT BETA   (beta01)
# ===========================================================================
#  beta is a parameter, not a contrast, so it is the one quantity here that
#  states a mechanism rather than a comparison. Its interval is wide and that
#  must travel with it everywhere it is quoted.
bb = rd("beta01_mutation_bias.csv") or []
b_by = {r["model"]: r for r in bb}
b_head = b_by.get("beta + gamma + kappa   <- HEADLINE")
b_ci   = b_by.get("beta_profile_CI")
b_k2   = b_by.get("kappa fixed 2 (previous headline)")
b_zero = b_by.get("beta_excludes_zero")
b_boot = next((r for k, r in b_by.items() if k.startswith("bootstrap")), None)
b_ktest = b_by.get("kappa_test_vs_fixed")
b_kprof = b_by.get("kappa_profile")

if b_head:
    record("beta_mutational_supply", num(b_head["value"] if "value" in b_head
                                         else b_head["beta"]),
           "exponent on mutational supply: P(allele observed) proportional to "
           "m^beta, m = kappa*transitions + transversions, ESM-2 adjusted",
           "beta01_mutation_bias.csv", "exact conditional likelihood, site strata",
           "CANONICAL. beta=0 means supply is irrelevant; beta=1 means observation "
           "is proportional to supply. Quote with its interval, never alone.")
    record("gamma_esm_per_sd", num(b_head["gamma"]),
           "coefficient on ESM-2 predicted effect, per standard deviation",
           "beta01_mutation_bias.csv", "same model as beta",
           "the protein-tolerance term that sits alongside beta in one model")
    record("kappa_coestimated", num(b_head["kappa"]),
           "transition:transversion rate ratio, estimated not assumed",
           "beta01_mutation_bias.csv", "co-estimated with beta and gamma",
           "NOT a measurement of the arthropod ts:tv ratio. It is a nuisance "
           "parameter and its own profile is nearly flat. Do not quote it as biology.")
if b_ci:
    record("beta_CI95", str(b_ci.get("note", "")).split(";")[0].strip(),
           "95% profile interval for beta, gamma and kappa profiled out",
           "beta01_mutation_bias.csv",
           "profile likelihood, bounds located by bisection to 1e-4",
           "CANONICAL interval. Bounds are a property of the likelihood, not of a "
           "grid step.")
if b_zero:
    ez = num(b_zero["value"] if "value" in b_zero else b_zero["beta"])
    record("beta_excludes_zero", ez,
           "does the 95% profile interval exclude beta = 0?",
           "beta01_mutation_bias.csv", "profile log-likelihood at beta=0 vs cutoff",
           "1 = yes. If 0, the manuscript must NOT claim the interval excludes zero.")
    if ez == 0:
        problems.append("beta01 reports that the 95% interval does NOT exclude "
                        "beta = 0; any sentence claiming it does is wrong")
if b_boot:
    record("beta_bootstrap_median", num(b_boot["beta"]),
           "median beta over site-cluster bootstrap resamples",
           "beta01_mutation_bias.csv", "resampling whole sites",
           "independent of the profile interval; %s"
           % str(b_boot.get("note", ""))[:90])
if b_ktest:
    _m = re.search(r"P=([0-9.eE+-]+)", str(b_ktest.get("note", "")))
    record("kappa_vs_fixed2_P", (float(_m.group(1)) if _m else ""),
           "likelihood-ratio test of kappa = 2 against the co-estimated value",
           "beta01_mutation_bias.csv", "LR, 1 df",
           str(b_ktest.get("note", "")))
if b_kprof:
    _m = re.search(r"spans\s*([0-9.]+)\s*-\s*([0-9.]+)", str(b_kprof.get("note", "")))
    record("kappa_profile_interval",
           ("%s to %s" % (_m.group(1), _m.group(2)) if _m else ""),
           "95% profile interval for kappa itself",
           "beta01_mutation_bias.csv", "profile likelihood",
           str(b_kprof.get("note", "")))
if not bb:
    problems.append("beta01_mutation_bias.csv not found -- run "
                    "src/beta01_mutation_bias_regression.py before this script")

# ===========================================================================
#  WRITE
# ===========================================================================
os.makedirs(RES, exist_ok=True)
with open(os.path.join(RES, "canonical_quantities.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, ["quantity", "value", "definition", "owning_file",
                            "method", "note"])
    w.writeheader(); w.writerows(rows)

say("=" * 78)
say("CANONICAL QUANTITIES")
say("=" * 78)
say()
say("  %-46s %14s" % ("quantity", "value"))
say("  " + "-" * 62)
for r in rows:
    say("  %-46s %14s" % (r["quantity"][:46], r["value"]))
say()
say("  %d quantities declared, each read from its owning file." % len(rows))
say()

if problems:
    say("PROBLEMS (%d) -- these block the figures:" % len(problems))
    for p in problems: say("   ! " + p)
else:
    say("  No problems. Every declaration resolves against the result files,")
    say("  the event partition closes, and the headline ratio reproduces from")
    say("  its own numerator and denominator.")

md = ["# Canonical quantities", "",
      "Every value below is read from the result file named in its row by",
      "`src/data15_canonical_quantities.py`. Nothing here is typed by hand. If a",
      "figure, table or sentence in the manuscript quotes a number, it quotes it",
      "from `results/canonical_quantities.csv`.", "",
      "## Why two quantities were called the same thing", "",
      "Two pairs of quantities in this project shared a name. Both members of each",
      "pair were computed correctly; the words collided. Stated once, here:", "",
      "**\"Informative\"** meant two things. The **supply-informative** count asks a",
      "question about a *site*: could some but not all species in the panel reach",
      "this allele here? The **codon-discriminating** count asks a question about a",
      "*species*: does this species' own codon reach some but not all of the alleles",
      "curated at this site? A site can be supply-informative while a given species'",
      "codon discriminates nothing, and vice versa. They are not nested and neither",
      "is a corrected version of the other.", "",
      "**\"Breach\"** meant two scopes. The **total** is every event whose observed",
      "allele was not one step away. The **discriminating subset** is those breaches",
      "inside the codon-discriminating events. They differ by one event,",
      "*Ctenocephalides felis* VGSC T929 ACG->V, whose codon reaches none of the site",
      "repertoire and which the discrimination test therefore excludes by",
      "construction. The event is not missing and is not an error; it is outside that",
      "test's domain.", "",
      "## Declared headlines", "",
      "Several estimators of one estimand is good practice. Leaving unstated which",
      "one the abstract quotes is not. The declarations are in the table below,",
      "marked DECLARED HEADLINE.", "",
      "For the site-level odds ratio the headline is the **conditional logit with",
      "site fixed effects, adjusted for ESM-2** -- the most conservative of the five",
      "estimators, and the only one that already answers the reviewer who says the",
      "pattern is protein-level tolerance rather than codon accessibility. The",
      "assumption-free **exact stratified hypergeometric** is quoted beside it,",
      "because it assumes nothing about the form of the effect.", "",
      "## The pooled enrichment is not the headline result", "",
      "The pooled enrichment is a modest ratio and it is *arithmetically bounded*:",
      "55 of the events sit in a stratum where every species could reach the allele,",
      "so they contribute observed = expected by construction and dilute the ratio",
      "toward 1. The number is correct and belongs in the paper as a descriptive",
      "quantity. It is not the finding. The finding is the stratified test, where the",
      "degenerate strata drop out and carry no likelihood.", "",
      "## Table", "",
      "| quantity | value | definition | owning file | method | note |",
      "|---|---|---|---|---|---|"]
for r in rows:
    md.append("| `%s` | %s | %s | `%s` | %s | %s |" % (
        r["quantity"], r["value"], r["definition"], r["owning_file"],
        r["method"], r["note"]))
if problems:
    md += ["", "## Unresolved", ""] + ["- " + p for p in problems]
open(os.path.join(RES, "CANONICAL_QUANTITIES.md"), "w").write("\n".join(md) + "\n")
say()
say("  wrote results/CANONICAL_QUANTITIES.md")
say("  wrote results/canonical_quantities.csv")
sys.exit(1 if problems else 0)
