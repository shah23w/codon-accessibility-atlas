#!/usr/bin/env python3
"""
step9_02  The mutational-supply exponent beta, placed in context.

Three questions.

1  IS BETA CARRIED BY ONE GENE?
   `results/beta01_mutation_bias.csv` reports that dropping VGSC moves beta to
   -0.0909, which has to be reported or resolved
   it. Both of those treat -0.0909 as an estimate. It is not: that fit returns
   kappa = 100.0, which is exactly beta01's KAPPA_MAX, so the optimiser ran to
   its boundary on a subset that cannot identify the nuisance parameter.

   The right question is not "what is beta without VGSC" but "is there evidence
   that beta differs between VGSC and the rest". That is a likelihood-ratio
   test on an interaction, computed here, and it does not depend on either
   subset being individually well conditioned.

2  HOW DOES BETA COMPARE WITH THE PUBLISHED VALUES?
   Cano et al. 2022 PNAS 119:e2119720119 fit the same estimand: their
   codon-frequency offset cancels when you condition on the starting codon,
   leaving the model used here exactly.
   Their three species are clonal or non-recombining. Ours are sexual, which is
   the gap their own Discussion names as open.

   Two theories predict different values for a sexual pest and we test against
   both rather than choosing the flattering one.

3  WHAT WOULD IT TAKE TO SEPARATE THEM?
   Our interval contains both predictions. The honest follow-up is to say how
   much more data would be needed to tell them apart. Answered by parametric
   bootstrap under the fitted model.

Everything is recomputed here from the same source file beta01 reads, and
section 0 checks that this module reproduces beta01's published headline before
any new claim is made. If that check fails nothing below it is quotable.

writes
    results/step9_02_beta_context.csv       one row per reported quantity
    results/step9_02_gene_information.csv   the per-gene information census
    results/step9_02_power.csv              the power curve
    logs/step9_02_beta_in_context.log
"""
import csv, math, os, sys, time

T0 = time.time()
try:
    import numpy as np
except ImportError:
    sys.exit("needs numpy")
try:
    from scipy.optimize import minimize
    from scipy.stats import chi2
except ImportError:
    sys.exit("needs scipy")

DIR_RES = "results"
DIR_LOG = "logs"
TAG = "_primary"
SRC = os.path.join(DIR_RES, "mutation_spectrum_substitutions%s.csv" % TAG)
BETA01 = os.path.join(DIR_RES, "beta01_mutation_bias.csv")

os.makedirs(DIR_LOG, exist_ok=True)
_lines = []
def log(s=""):
    print(("[%6.1fs] %s" % (time.time() - T0, s)) if s else "")
    _lines.append(s)

FAILURES = []
def chk(label, condition, detail=""):
    """A check that can fail. Never assert x == x."""
    if condition:
        log("  [PASS] %s%s" % (label, ("   " + detail) if detail else ""))
    else:
        log("  [FAIL] %s%s" % (label, ("   " + detail) if detail else ""))
        FAILURES.append(label)
    return condition

RNG = np.random.default_rng(20261101)

# --------------------------------------------------------------------------- #
#  data and likelihood -- deliberately the same construction as beta01
# --------------------------------------------------------------------------- #
def load():
    if not os.path.exists(SRC):
        sys.exit("FATAL: missing %s -- run step6_04_mutation_spectrum.py first" % SRC)
    rows = list(csv.DictReader(open(SRC)))
    strata, order, dropped = {}, [], []
    for r in rows:
        key = (r["gene"], r["site_label"])
        try:
            ts = float(r["n_transition_paths"]); tv = float(r["n_transversion_paths"])
            esm = float(r["esm_score"]); obs = int(r["observed"])
        except (ValueError, KeyError):
            dropped.append(r.get("site_label", "?")); continue
        if ts + tv <= 0:
            dropped.append("%s:%s zero-supply" % (r["gene"], r["candidate_aa"]))
            continue
        if key not in strata:
            strata[key] = {"ts": [], "tv": [], "esm": [], "obs": [], "aa": []}
            order.append(key)
        d = strata[key]
        d["ts"].append(ts); d["tv"].append(tv); d["esm"].append(esm)
        d["obs"].append(obs); d["aa"].append(r["candidate_aa"])
    usable, degenerate = [], []
    for key in order:
        d = strata[key]; n = len(d["obs"]); k = sum(d["obs"])
        (usable if 0 < k < n else degenerate).append((key, n, k))
    return strata, usable, degenerate, dropped

strata, usable, degenerate, dropped = load()

log("=" * 74)
log("step9_02  beta in context: gene dependence, published values, power")
log("=" * 74)
log("  source            : %s" % SRC)
log("  sites read        : %d" % len(strata))
log("  usable strata     : %d" % len(usable))
log("  degenerate strata : %d" % len(degenerate))
if dropped:
    log("  rows set aside    : %d  %s" % (len(dropped), dropped[:4]))

all_esm = np.array([e for k, _, _ in usable for e in strata[k]["esm"]], float)
ESM_MU, ESM_SD = all_esm.mean(), all_esm.std() or 1.0
log("  ESM standardised  : mean %.3f, sd %.3f" % (ESM_MU, ESM_SD))

STRATA = []
for key, n, k in usable:
    d = strata[key]
    STRATA.append(dict(
        key=key, gene=key[0], n=n, k=k,
        ts=np.array(d["ts"], float), tv=np.array(d["tv"], float),
        esm=(np.array(d["esm"], float) - ESM_MU) / ESM_SD,
        obs=np.array(d["obs"], int), aa=d["aa"]))

KAPPA_MIN, KAPPA_MAX = 0.05, 100.0

def log_ek(w, k):
    e = np.zeros(k + 1); e[0] = 1.0
    for x in w:
        for j in range(min(k, len(e) - 1), 0, -1):
            e[j] += x * e[j - 1]
    return math.log(e[k]) if e[k] > 0 else -math.inf

def _ll_stratum(S, beta, gamma, kappa, obs=None):
    m = kappa * S["ts"] + S["tv"]
    if np.any(m <= 0):
        return -math.inf
    eta = beta * np.log(m) + gamma * S["esm"]
    eta = eta - eta.max()
    o = S["obs"] if obs is None else obs
    num = float(eta[o == 1].sum())
    den = log_ek(np.exp(eta), int(o.sum()))
    if not np.isfinite(den):
        return -math.inf
    return num - den

def neg_ll(beta, gamma, kappa, subset=None, obsmap=None):
    if not (KAPPA_MIN <= kappa <= KAPPA_MAX):
        return math.inf
    tot = 0.0
    for S in (STRATA if subset is None else subset):
        v = _ll_stratum(S, beta, gamma, kappa,
                        None if obsmap is None else obsmap[id(S)])
        if not np.isfinite(v):
            return math.inf
        tot += v
    return -tot

def fit(free_names, fixed, x0, subset=None, obsmap=None):
    """Nelder-Mead with the box enforced inside the objective, because
    Nelder-Mead ignores `bounds`."""
    def obj(x):
        p = dict(fixed)
        for nm, v in zip(free_names, x):
            p[nm] = v
        if not (KAPPA_MIN <= p.get("kappa", 2.0) <= KAPPA_MAX):
            return 1e12
        return neg_ll(p["beta"], p["gamma"], p["kappa"], subset, obsmap)
    r = minimize(obj, np.array(x0, float), method="Nelder-Mead",
                 options=dict(xatol=1e-6, fatol=1e-8, maxiter=8000, maxfev=8000))
    p = dict(fixed)
    for nm, v in zip(free_names, r.x):
        p[nm] = v
    return p, -float(r.fun)

def fit_best(free_names, fixed, starts, subset=None, obsmap=None):
    best = None
    for x0 in starts:
        p, ll = fit(free_names, fixed, x0, subset, obsmap)
        if best is None or ll > best[1]:
            best = (p, ll)
    return best

STARTS3 = [(0.5, 0.8, 2.0), (1.0, 0.5, 2.0), (0.0, 1.0, 4.0), (0.3, 0.8, 6.0)]
STARTS2 = [(0.5, 0.8), (1.0, 0.5), (0.0, 1.0)]

def simulate_obs_h0(S, beta, gamma, kappa=2.0):
    """Draw this stratum's k observed alleles without replacement, with
    probability proportional to exp(eta). That is exactly the sampling scheme
    the conditional likelihood assumes, so data simulated this way are drawn
    from the fitted model rather than from an approximation to it."""
    m = kappa * S["ts"] + S["tv"]
    eta = beta * np.log(m) + gamma * S["esm"]
    w = np.exp(eta - eta.max())
    idx = RNG.choice(len(w), size=S["k"], replace=False, p=w / w.sum())
    o = np.zeros(len(w), dtype=int); o[idx] = 1
    return o

ROWS = []
def emit(quantity, value, note=""):
    ROWS.append(dict(quantity=quantity, value=value, note=note))

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("0  DOES THIS MODULE REPRODUCE beta01's PUBLISHED HEADLINE?")
log("-" * 74)
log("   Recomputed from the same source file, not read from beta01's output.")

p, ll_full = fit_best(("beta", "gamma", "kappa"), {}, STARTS3)
BETA, GAMMA, KAPPA = p["beta"], p["gamma"], p["kappa"]
log("   refit here        : beta %.4f  gamma %.4f  kappa %.3f  logLik %.5f"
    % (BETA, GAMMA, KAPPA, ll_full))

pub_beta = pub_ll = None
if os.path.exists(BETA01):
    for r in csv.DictReader(open(BETA01)):
        if r["model"].startswith("beta + gamma + kappa"):
            pub_beta = float(r["beta"]); pub_ll = float(r["logLik"]); break
if pub_beta is None:
    log("   [FAIL] could not read beta01's headline row from %s" % BETA01)
    FAILURES.append("beta01 headline unreadable")
else:
    log("   beta01 published  : beta %.4f                    logLik %.5f"
        % (pub_beta, pub_ll))
    chk("this module reproduces beta01's beta",
        abs(BETA - pub_beta) < 2e-3, "|d| = %.2e" % abs(BETA - pub_beta))
    chk("this module reproduces beta01's logLik",
        abs(ll_full - pub_ll) < 1e-2, "|d| = %.2e" % abs(ll_full - pub_ll))

emit("beta_headline_refit", round(BETA, 4),
     "recomputed here; beta01 publishes %s" % (("%.4f" % pub_beta) if pub_beta else "n/a"))

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("1  IS BETA CARRIED BY ONE GENE?")
log("-" * 74)

genes = sorted({S["gene"] for S in STRATA})
log()
log("  1a  per-gene information census")
log("      a stratum informs beta only if 0 < observed < candidates; a gene's")
log("      weight is the number of such strata and the alleles inside them.")
log()
log("      %-12s %8s %12s %11s" % ("gene", "strata", "candidates", "observed"))
GROWS = []
for g in genes:
    sub = [S for S in STRATA if S["gene"] == g]
    nstr = len(sub); ncand = sum(S["n"] for S in sub); nobs = sum(S["k"] for S in sub)
    log("      %-12s %8d %12d %11d" % (g, nstr, ncand, nobs))
    GROWS.append(dict(gene=g, strata=nstr, candidates=ncand, observed=nobs))
tot_str = len(STRATA); tot_cand = sum(S["n"] for S in STRATA)
log("      %-12s %8d %12d %11d" % ("TOTAL", tot_str, tot_cand,
                                   sum(S["k"] for S in STRATA)))

vg = [S for S in STRATA if S["gene"] == "VGSC"]
share = len(vg) / tot_str if tot_str else 0.0
log()
log("      VGSC holds %d of %d informative strata (%.1f%%)"
    % (len(vg), tot_str, 100 * share))
emit("vgsc_share_of_informative_strata", round(share, 4),
     "%d of %d" % (len(vg), tot_str))

# The partitioned fits below hold kappa at 2, because co-estimating it without
# VGSC runs it to the boundary. The pooled fit at that SAME fixed kappa is the
# only like-for-like comparator for them; the published headline co-estimates
# kappa and is therefore not comparable with them number for number.
_p_k2, _ = fit_best(("beta", "gamma"), {"kappa": 2.0}, STARTS2)
log("      pooled beta at kappa = 2 (the partition's comparator): %.4f"
    % _p_k2["beta"])
emit("beta_pooled_kappa_fixed2", round(_p_k2["beta"], 4),
     "all %d strata, kappa fixed 2 and gamma free: the like-for-like "
     "comparator for the VGSC and non-VGSC fits, which fix kappa the same way"
     % tot_str)

# Information per stratum, by gene. The gene contrast invites the reading that
# sodium-channel sites tolerate more distinct replacements than catalytic or
# pore-lining residues elsewhere, so that supply has more alternatives to
# order there. In this record they tolerate FEWER, and that reading is
# therefore not available from these data.
_dps_vg = sum(S["k"] for S in vg) / len(vg)
_other = [S for S in STRATA if S["gene"] != "VGSC"]
_dps_ot = sum(S["k"] for S in _other) / len(_other)
log("      documented alleles per stratum: VGSC %.2f, other genes %.2f"
    % (_dps_vg, _dps_ot))
emit("docs_per_stratum_vgsc", round(_dps_vg, 2),
     "%d documented alleles over %d sodium-channel strata"
     % (sum(S["k"] for S in vg), len(vg)))
emit("docs_per_stratum_other", round(_dps_ot, 2),
     "%d documented alleles over %d strata in the other genes; per gene %s"
     % (sum(S["k"] for S in _other), len(_other),
        ", ".join("%s %.2f" % (d["gene"], d["observed"] / d["strata"])
                  for d in GROWS if d["gene"] != "VGSC")))
emit("n_strata_non_vgsc", len(_other),
     "informative strata outside the sodium channel, against %d within it"
     % len(vg))

log()
log("  1b  leave-one-gene-out, with the boundary made visible")
log("      kappa free, then kappa fixed at 2. beta01's drop-VGSC row reports")
log("      kappa = 100.0, which is exactly KAPPA_MAX; a fit sitting on its")
log("      boundary is a failure to identify, not an estimate.")
log()
log("      %-12s %7s %9s %9s %9s %9s" %
    ("dropped", "strata", "beta", "kappa", "at bound", "beta(k=2)"))
for g in genes:
    sub = [S for S in STRATA if S["gene"] != g]
    if len(sub) < 3:
        log("      %-12s %7d   too few strata to fit" % (g, len(sub)))
        GROWS = [d for d in GROWS]
        continue
    pf, _ = fit_best(("beta", "gamma", "kappa"), {}, STARTS3, subset=sub)
    at_bound = (pf["kappa"] >= KAPPA_MAX - 1e-6) or (pf["kappa"] <= KAPPA_MIN + 1e-6)
    pk, _ = fit_best(("beta", "gamma"), {"kappa": 2.0}, STARTS2, subset=sub)
    log("      %-12s %7d %9.4f %9.3f %9s %9.4f"
        % (g, len(sub), pf["beta"], pf["kappa"], "YES" if at_bound else "no",
           pk["beta"]))
    for d in GROWS:
        if d["gene"] == g:
            d["drop_beta_kappa_free"] = round(pf["beta"], 4)
            d["drop_kappa"] = round(pf["kappa"], 3)
            d["drop_kappa_at_bound"] = int(at_bound)
            d["drop_beta_kappa_fixed2"] = round(pk["beta"], 4)

sub_novg = [S for S in STRATA if S["gene"] != "VGSC"]
pf_novg, ll_novg = fit_best(("beta", "gamma", "kappa"), {}, STARTS3, subset=sub_novg)
bound_hit = pf_novg["kappa"] >= KAPPA_MAX - 1e-6
chk("the drop-VGSC fit sits on the kappa boundary (so -0.0909 is not an estimate)",
    bound_hit, "kappa = %.3f, KAPPA_MAX = %.1f" % (pf_novg["kappa"], KAPPA_MAX))
emit("drop_vgsc_kappa_at_bound", int(bound_hit),
     "kappa = %.3f against KAPPA_MAX = %.1f" % (pf_novg["kappa"], KAPPA_MAX))

log()
log("  1c  is beta estimable at all without VGSC?  profile interval, kappa = 2")
log("      kappa is fixed because the free-kappa fit above does not identify it.")

def profile_ci(subset, kappa_fix=2.0, lo=-4.0, hi=4.0, step=0.05):
    """Profile beta over a grid with gamma profiled out; return (mle, lo, hi,
    ll_at_mle). Bounds are the outermost grid crossings of the 95% cutoff, so a
    bound reported at the grid edge means 'not bounded by the data'."""
    best_b, best_ll = None, -math.inf
    grid = np.arange(lo, hi + 1e-9, step)
    lls = []
    for b in grid:
        pg, llb = fit_best(("gamma",), {"beta": float(b), "kappa": kappa_fix},
                           [(0.8,), (0.0,)], subset=subset)
        lls.append(llb)
        if llb > best_ll:
            best_ll, best_b = llb, float(b)
    cut = best_ll - chi2.ppf(0.95, 1) / 2.0
    inside = [float(b) for b, l in zip(grid, lls) if l >= cut]
    if not inside:
        return best_b, float("nan"), float("nan"), best_ll, cut
    return best_b, min(inside), max(inside), best_ll, cut

sub_vg = [S for S in STRATA if S["gene"] == "VGSC"]
b_nv, lo_nv, hi_nv, ll_nv, cut_nv = profile_ci(sub_novg)
width_nv = hi_nv - lo_nv
b_vg, lo_vg, hi_vg, ll_vg, cut_vg = profile_ci(sub_vg)
width_vg = hi_vg - lo_vg
b_all, lo_all, hi_all, ll_all2, cut_all = profile_ci(STRATA)
width_all = hi_all - lo_all
log("      all genes    : beta = %.3f   95%% profile [%.3f, %.3f]   width %.3f"
    % (b_all, lo_all, hi_all, width_all))
log("      VGSC only    : beta = %.3f   95%% profile [%.3f, %.3f]   width %.3f"
    % (b_vg, lo_vg, hi_vg, width_vg))
log("      without VGSC : beta = %.3f   95%% profile [%.3f, %.3f]   width %.3f"
    % (b_nv, lo_nv, hi_nv, width_nv))
log()
log("      Width is the diagnostic. If the non-VGSC interval is much wider AND")
log("      contains the pooled estimate, the subset is uninformative rather than")
log("      contradictory, and a point estimate from it means nothing.")
contains_pooled = lo_nv <= b_all <= hi_nv
log("      non-VGSC interval contains the pooled estimate %.3f : %s"
    % (b_all, "yes" if contains_pooled else "NO"))
log("      non-VGSC interval is %.2fx the width of the pooled one"
    % (width_nv / width_all if width_all else float("nan")))
emit("beta_vgsc_only", round(b_vg, 4),
     "95%% profile [%.3f, %.3f], kappa fixed 2, %d strata" % (lo_vg, hi_vg, len(sub_vg)))
emit("beta_non_vgsc_contains_pooled", int(contains_pooled),
     "non-VGSC 95%% [%.3f, %.3f] against pooled %.3f" % (lo_nv, hi_nv, b_all))
emit("beta_without_vgsc", round(b_nv, 4),
     "95%% profile [%.3f, %.3f], kappa fixed 2" % (lo_nv, hi_nv))
emit("beta_without_vgsc_ci_width", round(width_nv, 4),
     "against %.4f with all genes; a wider interval means less information, "
     "not a different beta" % width_all)

log()
log("  1d  THE DECISIVE TEST: does beta actually DIFFER between VGSC and the rest?")
log("      H0: one beta for all genes.   H1: a separate beta for VGSC.")
log("      Both fits use the whole dataset, so neither depends on a subset being")
log("      well conditioned. kappa fixed at 2 in both (beta01's own LR test finds")
log("      kappa = 2 and the co-estimated 4.703 indistinguishable, P = 0.298).")

def neg_ll_split(beta_vgsc, beta_rest, gamma, kappa=2.0):
    tot = 0.0
    for S in STRATA:
        b = beta_vgsc if S["gene"] == "VGSC" else beta_rest
        v = _ll_stratum(S, b, gamma, kappa)
        if not np.isfinite(v):
            return math.inf
        tot += v
    return -tot

p0, ll0 = fit_best(("beta", "gamma"), {"kappa": 2.0}, STARTS2)

def obj_split(x):
    return neg_ll_split(x[0], x[1], x[2])
best1 = None
for x0 in [(0.5, 0.5, 0.8), (1.0, 0.0, 0.8), (0.3, 0.3, 1.0), (0.7, -0.2, 0.6)]:
    r = minimize(obj_split, np.array(x0, float), method="Nelder-Mead",
                 options=dict(xatol=1e-6, fatol=1e-8, maxiter=8000, maxfev=8000))
    if best1 is None or -r.fun > best1[1]:
        best1 = (r.x, -float(r.fun))
x1, ll1 = best1
LR = 2.0 * (ll1 - ll0)
P_split = float(chi2.sf(LR, 1)) if LR > 0 else 1.0
log()
log("      H0  one beta          : beta %.4f                 logLik %.5f"
    % (p0["beta"], ll0))
log("      H1  VGSC / rest split : beta_VGSC %.4f  beta_rest %.4f  logLik %.5f"
    % (x1[0], x1[1], ll1))
log("      LR = %.4f  df = 1  P = %.4f" % (LR, P_split))
emit("beta_gene_interaction_LR", round(LR, 4), "df = 1")
emit("beta_gene_interaction_P", round(P_split, 5),
     "H0: one beta for all genes vs H1: separate beta for VGSC; kappa fixed 2")

log()
log("  1e  IS THAT P VALUE TRUSTWORTHY AT THIS SAMPLE SIZE?")
log("      The split above puts 34 strata on one side and 13 on the other, and a")
log("      1 df asymptotic chi-square is not credible there. The null")
log("      distribution of the LR statistic is therefore simulated directly:")
log("      generate data under H0 (one beta for every gene), refit both models,")
log("      and read the empirical P as the fraction of simulated LRs at least")
log("      as large as the observed one.")

NBOOT_LR = int(os.environ.get("STEP9_NBOOT_LR", "300"))
log()
log("      simulations : %d   (STEP9_NBOOT_LR)" % NBOOT_LR)

def _fit_h0(sim):
    def obj(x):
        tot = 0.0
        for S, o in sim:
            v = _ll_stratum(S, x[0], x[1], 2.0, obs=o)
            if not np.isfinite(v):
                return 1e12
            tot += v
        return -tot
    best = None
    for s in [(0.65, 0.75), (0.0, 1.0)]:
        r = minimize(obj, np.array(s, float), method="Nelder-Mead",
                     options=dict(xatol=1e-5, fatol=1e-7, maxiter=4000, maxfev=4000))
        if best is None or -r.fun > best:
            best = -float(r.fun)
    return best

def _fit_h1(sim):
    def obj(x):
        tot = 0.0
        for S, o in sim:
            b = x[0] if S["gene"] == "VGSC" else x[1]
            v = _ll_stratum(S, b, x[2], 2.0, obs=o)
            if not np.isfinite(v):
                return 1e12
            tot += v
        return -tot
    best = None
    for s in [(1.0, 0.0, 0.75), (0.65, 0.65, 0.75)]:
        r = minimize(obj, np.array(s, float), method="Nelder-Mead",
                     options=dict(xatol=1e-5, fatol=1e-7, maxiter=4000, maxfev=4000))
        if best is None or -r.fun > best:
            best = -float(r.fun)
    return best

B0, G0 = p0["beta"], p0["gamma"]
lr_null = []
for _ in range(NBOOT_LR):
    sim = [(S, simulate_obs_h0(S, B0, G0)) for S in STRATA]
    lr_null.append(2.0 * (_fit_h1(sim) - _fit_h0(sim)))
lr_null = np.array([x for x in lr_null if np.isfinite(x)])
P_boot = float((lr_null >= LR).sum() + 1) / (len(lr_null) + 1)
log()
log("      observed LR            : %.4f" % LR)
log("      asymptotic chi2(1) P   : %.4f" % P_split)
log("      bootstrap P (%d sims)  : %.4f" % (len(lr_null), P_boot))
log("      simulated LR quantiles : 50%% %.2f   90%% %.2f   95%% %.2f   99%% %.2f"
    % tuple(np.percentile(lr_null, [50, 90, 95, 99])))
emit("beta_gene_interaction_P_bootstrap", round(P_boot, 4),
     "%d simulations under H0; the asymptotic P is %.4f" % (len(lr_null), P_split))

# --------------------------------------------------------------------------- #
#  Does the exponent rest on sites with no documented resistance at all?
#
#  A stratum here is a site and its outcome is the site's curated repertoire,
#  so a site enters the fit whether or not any species is recorded as having
#  evolved resistance there. Ten of the 47 carry no species-level event, and
#  at those the supply covariate is the panel's codon composition in taxa that
#  have never been reported resistant. Refitting on the 37 sites that do carry
#  an event is the sensitivity that observation asks for. kappa is held at 2
#  so the estimate is comparable with the gene-partitioned fits above.
# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("4  SENSITIVITY: sites that carry a species-level event")
log("-" * 74)
_evp = os.path.join(DIR_RES, "core_test_events_primary.csv")
if not os.path.exists(_evp):
    sys.exit("FATAL: %s is missing, so the event-bearing-site sensitivity "
             "cannot be fitted." % _evp)
_ev_sites = {(r["gene"], r["site_label"])
             for r in csv.DictReader(open(_evp))}
sub_ev = [S for S in STRATA if S["key"] in _ev_sites]
if not 0 < len(sub_ev) < len(STRATA):
    sys.exit("FATAL: %d of %d strata carry an event, so the sensitivity is "
             "either empty or the whole fit." % (len(sub_ev), len(STRATA)))
p_all2, _ = fit_best(("beta", "gamma"), {"kappa": 2.0}, STARTS2)
p_ev2, _ = fit_best(("beta", "gamma"), {"kappa": 2.0}, STARTS2, subset=sub_ev)
log("      strata carrying an event : %d of %d" % (len(sub_ev), len(STRATA)))
log("      beta, all strata (k=2)   : %.4f" % p_all2["beta"])
log("      beta, event-bearing only : %.4f" % p_ev2["beta"])
log("      difference               : %.4f" % (p_ev2["beta"] - p_all2["beta"]))
emit("n_strata_with_event", len(sub_ev),
     "of %d; the rest carry a curated repertoire but no species-level event"
     % len(STRATA))
emit("beta_event_bearing_sites_kappa2", round(p_ev2["beta"], 4),
     "kappa fixed 2, %d strata that carry at least one of the curated events; "
     "the comparable all-strata value is %.4f"
     % (len(sub_ev), p_all2["beta"]))

log()
log("      VERDICT")
if P_boot >= 0.05:
    log("      -> The bootstrap does NOT reject a common beta (P = %.4f), although"
        % P_boot)
    log("         the asymptotic test did (P = %.4f). The asymptotic P was" % P_split)
    log("         optimistic at this sample size. Contradiction 6 is about")
    log("         INFORMATION, not heterogeneity: six genes supply %d informative"
        % len(sub_novg))
    log("         strata between them, which cannot pin beta down, and a point")
    log("         estimate from them means nothing. Report the interaction test")
    log("         and the interval widths; never report -0.0909 as a beta.")
else:
    log("      -> beta DIFFERS by gene and it survives the bootstrap (P = %.4f)."
        % P_boot)
    log("         beta_VGSC = %.3f [%.3f, %.3f] on %d strata; the remaining six"
        % (b_vg, lo_vg, hi_vg, len(sub_vg)))
    log("         genes give %.3f [%.3f, %.3f] on %d strata."
        % (b_nv, lo_nv, hi_nv, len(sub_novg)))
    log("         The pooled 0.4462 then averages unlike things and MUST NOT be")
    log("         presented as one property of the system. Two defensible ways")
    log("         forward, to be decided on evidence and not on convenience:")
    log("           (i)  report beta for VGSC, where it is measurable, and state")
    log("                that the other genes are too sparse to test;")
    log("           (ii) drop beta from the main text and keep accessibility,")
    log("                which survives every gene-drop, as the only claim.")
    log("         Note the non-VGSC interval still contains the pooled estimate,")
    log("         so this is a difference in what is MEASURABLE per gene, not a")
    log("         demonstration that mutational supply is irrelevant elsewhere.")

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("2  BETA AGAINST THE PUBLISHED VALUES AND THE TWO RIVAL PREDICTIONS")
log("-" * 74)
log("   Literature values are typed from the cited papers, not computed here,")
log("   and are labelled as such. Cano's beta and the one fitted here are the")
log("   same estimand; step9_03 recovers their values to confirm it.")
log()
log("   The interval used here is beta01's published profile interval, with")
log("   kappa profiled out, because the point estimate quoted alongside it is")
log("   the kappa-free one. The kappa = 2 intervals in section 1 are for")
log("   comparing subsets with each other, not with the literature; mixing the")
log("   two would pair a point estimate with an interval from a different model.")

LO_PUB = HI_PUB = None
for r in csv.DictReader(open(BETA01)):
    if r["model"] == "beta_profile_CI":
        import re as _re
        mm = _re.search(r"\[([-0-9.]+),\s*([-0-9.]+)\]", r["note"])
        if mm:
            LO_PUB, HI_PUB = float(mm.group(1)), float(mm.group(2))
        break
chk("beta01's published profile interval was read",
    LO_PUB is not None,
    "[%.4f, %.4f]" % (LO_PUB, HI_PUB) if LO_PUB is not None else "not found")
if LO_PUB is None:
    LO_PUB, HI_PUB = lo_all, hi_all

# (label, beta, ci_lo, ci_hi, reproduction, source)
LIT = [
    ("S. cerevisiae",   1.05, 0.89, 1.21, "clonal, laboratory",
     "Cano 2022 PNAS Table 1"),
    ("E. coli",         0.98, 0.71, 1.25, "clonal, laboratory",
     "Cano 2022 PNAS Table 1"),
    ("M. tuberculosis", 0.85, 0.31, 1.37, "little or no recombination",
     "Cano 2022 PNAS Table 1"),
]
log()
log("   %-18s %7s %18s  %-28s" % ("system", "beta", "95% CI", "reproduction"))
for lab, b, l, h, repro, src in LIT:
    log("   %-18s %7.2f  [%5.2f, %5.2f]  %-28s" % (lab, b, l, h, repro))
log("   %-18s %7.4f  [%5.4f, %5.4f]  %-28s"
    % ("THIS STUDY pooled", BETA, LO_PUB, HI_PUB, "sexual, outcrossing, natural"))
log("   %-18s %7.4f  [%5.4f, %5.4f]  %-28s"
    % ("  of which VGSC", b_vg, lo_vg, hi_vg, "kappa fixed 2, %d strata" % len(sub_vg)))

log()
log("   Overlap with each published interval (an overlap means our measurement")
log("   is compatible with theirs, not that the systems are the same):")
for lab, b, l, h, repro, src in LIT:
    ov = not (HI_PUB < l or LO_PUB > h)
    log("     %-18s %s" % (lab, "compatible" if ov else "INCOMPATIBLE"))
    emit("compatible_with_%s" % lab.replace(". ", "_").replace(" ", "_"),
         int(ov), "%s beta %.2f [%.2f, %.2f]" % (src, b, l, h))

log()
log("   The two rival predictions for a sexual pest:")
log()
log("     A  asexual mutational-supply curve.  Cano Fig. 4A read at the pest")
log("        per-site supply N*mu = 5.8e-9 * 1e7..1e8 = 0.058..0.58, which their")
log("        own code confirms is the x-axis quantity (MutRate = 1e-7 per")
log("        nucleotide, sizes.txt N = 1e3..1e7).          predicts beta ~ 0.1-0.6")
log("     B  sexual origin-fixation.  Gomez 2020 Proc R Soc B: 'in the absence of")
log("        strong linkage disequilibrium in more sexual populations, the")
log("        origin-fixation regime is likely to be a good approximation', where")
log("        mutation enters proportionally.                predicts beta ~ 1")
log()
log("     Cano's simulation sets initializeRecombinationRate(0) and")
log("     setCloningRate(1.0), so prediction A is an asexual extrapolation.")
log("     Our species are sexual. We report both and claim neither.")

A_LO, A_HI, B_POINT = 0.1, 0.6, 1.0
inA = not (HI_PUB < A_LO or LO_PUB > A_HI)
inB = (LO_PUB <= B_POINT <= HI_PUB)
log()
log("     our interval [%0.4f, %0.4f] overlaps prediction A (0.1-0.6) : %s"
    % (LO_PUB, HI_PUB, "yes" if inA else "no"))
log("     our interval contains prediction B (beta = 1)               : %s"
    % ("yes" if inB else "no"))
emit("overlaps_prediction_A_asexual_supply", int(inA), "A: beta 0.1-0.6")
emit("contains_prediction_B_sexual_origin_fixation", int(inB), "B: beta = 1")

discriminates = not (inA and inB)
chk("we correctly report that the data do NOT discriminate the two predictions",
    not discriminates,
    "both inside the interval" if not discriminates else
    "one prediction is excluded -- this is now a positive result, rewrite section 2")
emit("discriminates_A_vs_B", int(discriminates),
     "0 means the interval contains both predictions and no choice between "
     "them is supported")

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("3  HOW MUCH MORE DATA WOULD SEPARATE THE TWO PREDICTIONS?")
log("-" * 74)
log("   Parametric bootstrap. Simulate allele choices under the fitted model at")
log("   beta = 1 (prediction B), on a dataset built by replicating the observed")
log("   strata M times, then test H0: beta = 0.6 (the top of prediction A) by")
log("   profile likelihood ratio. Power is the rejection rate at alpha = 0.05.")
log("   kappa is fixed at 2 throughout so each simulated fit is two-parameter.")

NSIM = int(os.environ.get("STEP9_NSIM", "200"))
MULTS = [int(x) for x in os.environ.get("STEP9_MULTS", "1,2,4,6,8").split(",")]
BETA_TRUE, BETA_NULL = 1.0, 0.6
log()
log("   simulations per multiplier : %d   (STEP9_NSIM)" % NSIM)
log("   multipliers                : %s   (STEP9_MULTS)" % MULTS)
log("   truth beta = %.1f, tested against H0 beta = %.1f" % (BETA_TRUE, BETA_NULL))

simulate_obs = simulate_obs_h0

def fit_beta_on(sim_strata, kappa=2.0, beta_fixed=None):
    def obj(x):
        b = beta_fixed if beta_fixed is not None else x[0]
        g = x[-1]
        tot = 0.0
        for S, o in sim_strata:
            v = _ll_stratum(S, b, g, kappa, obs=o)
            if not np.isfinite(v):
                return 1e12
            tot += v
        return -tot
    x0 = [0.8] if beta_fixed is not None else [1.0, 0.8]
    best = None
    for s in ([(0.8,), (0.0,)] if beta_fixed is not None
              else [(1.0, 0.8), (0.4, 0.5)]):
        r = minimize(obj, np.array(s, float), method="Nelder-Mead",
                     options=dict(xatol=1e-5, fatol=1e-7, maxiter=4000, maxfev=4000))
        if best is None or -r.fun > best:
            best = -float(r.fun)
    return best

CRIT = chi2.ppf(0.95, 1)
PROWS = []
log()
log("   %6s %9s %9s %9s" % ("mult", "strata", "power", "mean beta"))
for M in MULTS:
    rej = 0; bs = []
    base = STRATA * M
    for _ in range(NSIM):
        sim = [(S, simulate_obs(S, BETA_TRUE, GAMMA)) for S in base]
        ll_free = fit_beta_on(sim)
        ll_null = fit_beta_on(sim, beta_fixed=BETA_NULL)
        if 2.0 * (ll_free - ll_null) > CRIT:
            rej += 1
    power = rej / NSIM
    log("   %6d %9d %9.3f %9s" % (M, len(base), power, "-"))
    PROWS.append(dict(multiplier=M, strata=len(base), power=round(power, 3),
                      nsim=NSIM, beta_true=BETA_TRUE, beta_null=BETA_NULL))

reached = [r for r in PROWS if r["power"] >= 0.80]
if reached:
    need = reached[0]
    log()
    log("   -> %dx the present informative strata (%d, from %d now) reaches 80%% power"
        % (need["multiplier"], need["strata"], len(STRATA)))
    emit("multiplier_for_80pc_power", need["multiplier"],
         "%d informative strata, to separate beta = 1.0 from beta = 0.6"
         % need["strata"])
else:
    log()
    log("   -> 80%% power not reached within the multipliers tested (max %dx)."
        % max(MULTS))
    log("      Report the power actually achieved rather than extrapolating.")
    emit("multiplier_for_80pc_power", "not reached within %dx" % max(MULTS),
         "highest power %.3f" % max(r["power"] for r in PROWS))

# --------------------------------------------------------------------------- #
log()
log("-" * 74)
log("WRITING")
log("-" * 74)
with open(os.path.join(DIR_RES, "step9_02_beta_context.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["quantity", "value", "note"],
                       lineterminator="\n")
    w.writeheader()
    for r in ROWS:
        w.writerow(r)
log("  wrote results/step9_02_beta_context.csv")

gf = ["gene", "strata", "candidates", "observed", "drop_beta_kappa_free",
      "drop_kappa", "drop_kappa_at_bound", "drop_beta_kappa_fixed2"]
with open(os.path.join(DIR_RES, "step9_02_gene_information.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=gf, lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    for r in GROWS:
        w.writerow(r)
log("  wrote results/step9_02_gene_information.csv")

with open(os.path.join(DIR_RES, "step9_02_power.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["multiplier", "strata", "power", "nsim",
                                       "beta_true", "beta_null"],
                       lineterminator="\n")
    w.writeheader()
    for r in PROWS:
        w.writerow(r)
log("  wrote results/step9_02_power.csv")

log()
log("=" * 74)
if FAILURES:
    log("RESULT: %d CHECK(S) FAILED" % len(FAILURES))
    for f in FAILURES:
        log("   - %s" % f)
else:
    log("RESULT: every check passed")
log("=" * 74)

with open(os.path.join(DIR_LOG, "step9_02_beta_in_context.log"), "w") as fh:
    fh.write("\n".join(_lines) + "\n")

sys.exit(1 if FAILURES else 0)
