#!/usr/bin/env python3
"""
step9_03  Does our estimator recover a published beta?

Cano et al.'s coefficient of mutational influence and the one fitted here are
the same estimand: their Eq. 2 carries codon frequency as an offset, so
conditioning on the starting codon cancels it and leaves P(a | c) proportional
to mu(c,a)^beta, which is this model. That is an argument. This module turns it
into a measurement.

We fit OUR conditional likelihood to THEIR published data and ask whether it
returns THEIR published beta. If it does, the equivalence is demonstrated and
not merely asserted, so anyone who doubts the algebra can check the number.

  Cano AV, Rozhonova H, Stoltzfus A, McCandlish DM, Payne JL (2022)
  Mutation bias shapes the spectrum of adaptive substitutions.
  PNAS 119:e2119720119.  doi:10.1073/pnas.2119720119
  Data: github.com/alejvcano/Mutbias2022 -> mutbias.zip -> files/data{MTB,Yeast,Ecoli}.txt

Published values to recover (their Table 1):
    S. cerevisiae    1.05  (95% CI 0.89-1.21)
    E. coli          0.98  (95% CI 0.71-1.25)
    M. tuberculosis  0.85  (95% CI 0.31-1.37)

ONE DIFFERENCE, STATED PLAINLY. Their response is a count of adaptive events
per codon-to-amino-acid type, fitted by negative binomial regression. Ours is
which allele a species carries at a site, so our atlas likelihood counts each
candidate at most once per stratum and uses an elementary symmetric polynomial.
Their data have many events per type, so the matching conditional here is the
multinomial-with-counts obtained by conditioning a Poisson count model on the
total events from each starting codon:

    log L = sum_c [ sum_a n(c,a) * eta(c,a)  -  N_c * log sum_a exp(eta(c,a)) ]
    eta(c,a) = beta * log mu(c,a)

So what is being benchmarked is the linear predictor and beta's meaning, which
is the thing in dispute, not the particular conditioning machinery the atlas
needs. Overdispersion is handled by their negative binomial and not by this
conditional, so exact agreement is not expected; agreement within their own
confidence interval is.

THEIR DATA ARE NOT REDISTRIBUTED. Point CANO_DIR at a local copy; the deposit
ships this script and the citation, never their files.

writes
    results/step9_03_benchmark_cano.csv
    logs/step9_03_benchmark_cano.log
"""
import csv, math, os, sys, time

T0 = time.time()
try:
    import numpy as np
except ImportError:
    sys.exit("needs numpy")
try:
    from scipy.optimize import minimize_scalar
    from scipy.stats import chi2
except ImportError:
    sys.exit("needs scipy")

DIR_RES = "results"
DIR_LOG = "logs"
CANO_DIR = os.environ.get("CANO_DIR", "external/cano2022")

os.makedirs(DIR_LOG, exist_ok=True)
_lines = []
def log(s=""):
    print(("[%6.1fs] %s" % (time.time() - T0, s)) if s else "")
    _lines.append(s)

FAILURES = []
def chk(label, condition, detail=""):
    if condition:
        log("  [PASS] %s%s" % (label, ("   " + detail) if detail else ""))
    else:
        log("  [FAIL] %s%s" % (label, ("   " + detail) if detail else ""))
        FAILURES.append(label)
    return condition

# (file, label, published beta, CI lo, CI hi)
SPECIES = [
    ("dataYeast.txt", "S. cerevisiae",   1.05, 0.89, 1.21),
    ("dataEcoli.txt", "E. coli",         0.98, 0.71, 1.25),
    ("dataMTB.txt",   "M. tuberculosis", 0.85, 0.31, 1.37),
]

log("=" * 74)
log("step9_03  benchmark: does our conditional likelihood recover Cano's beta?")
log("=" * 74)
log("  data directory : %s   (CANO_DIR)" % CANO_DIR)

if not os.path.isdir(CANO_DIR):
    log()
    log("  FATAL: %s not found." % CANO_DIR)
    log("  Obtain it once with:")
    log("    mkdir -p %s && cd %s" % (CANO_DIR, CANO_DIR))
    log("    curl -sSL -o mutbias.zip \\")
    log("      https://github.com/alejvcano/Mutbias2022/raw/main/mutbias.zip")
    log("    unzip -j mutbias.zip 'mutbias/files/data*.txt'")
    sys.exit(2)


def load_species(path):
    """Group their table into strata keyed by starting codon.

    mu(c,a) is taken as exp(log_mutAll), which is the summed codon-to-amino-acid
    mutation rate THEY computed, so the benchmark tests our estimator and not
    our reconstruction of their covariate."""
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    strata, order = {}, []
    bad = 0
    for r in rows:
        try:
            c = r["codon"]
            ev = float(r["events"])
            lm = float(r["log_mutAll"])
        except (KeyError, ValueError):
            bad += 1
            continue
        if not np.isfinite(lm):
            bad += 1
            continue
        if c not in strata:
            strata[c] = {"logmu": [], "n": []}
            order.append(c)
        strata[c]["logmu"].append(lm)
        strata[c]["n"].append(ev)
    out = []
    for c in order:
        d = strata[c]
        n = np.array(d["n"], float)
        if n.sum() <= 0 or len(n) < 2:
            continue          # a codon with no events, or no choice to make,
            # carries no conditional information -- excluded by the
            # conditioning itself, not by a decision of ours
        out.append(dict(codon=c, logmu=np.array(d["logmu"], float), n=n))
    return out, len(rows), bad


def neg_ll(beta, strata):
    tot = 0.0
    for S in strata:
        eta = beta * S["logmu"]
        m = eta.max()
        den = m + math.log(float(np.exp(eta - m).sum()))
        tot += float((S["n"] * eta).sum()) - float(S["n"].sum()) * den
    return -tot


def fit(strata):
    r = minimize_scalar(neg_ll, args=(strata,), bounds=(-3.0, 4.0),
                        method="bounded", options=dict(xatol=1e-6))
    return float(r.x), -float(r.fun)


def profile_ci(strata, bhat, llhat, step=0.002, span=3.0):
    """Walk outward from the MLE until the profile drops below the 95% cutoff.
    Reports a bound at the search edge as non-finite rather than silently
    returning the edge as if the data had bounded it."""
    cut = llhat - chi2.ppf(0.95, 1) / 2.0
    out = []
    for direction in (-1, +1):
        b = bhat
        bound = None
        while abs(b - bhat) <= span:
            b += direction * step
            if -neg_ll(b, strata) < cut:
                bound = b
                break
        out.append(bound if bound is not None else float("nan"))
    return out[0], out[1]


ROWS = []
log()
log("  %-18s %8s %9s %18s %20s" %
    ("species", "codons", "events", "beta (ours)", "beta (published)"))

for fn, label, pub, lo_pub, hi_pub in SPECIES:
    path = os.path.join(CANO_DIR, fn)
    if not os.path.exists(path):
        log("  %-18s MISSING %s" % (label, path))
        FAILURES.append("missing %s" % fn)
        continue
    strata, nrow, bad = load_species(path)
    nev = sum(S["n"].sum() for S in strata)
    b, ll = fit(strata)
    lo, hi = profile_ci(strata, b, ll)
    log("  %-18s %8d %9d %8.3f [%5.3f,%5.3f] %8.2f [%5.2f,%5.2f]"
        % (label, len(strata), int(nev), b, lo, hi, pub, lo_pub, hi_pub))
    ROWS.append(dict(species=label, source_file=fn, informative_codons=len(strata),
                     events=int(nev), beta_ours=round(b, 4),
                     ci_lo_ours=round(lo, 4), ci_hi_ours=round(hi, 4),
                     beta_published=pub, ci_lo_published=lo_pub,
                     ci_hi_published=hi_pub,
                     rows_read=nrow, rows_set_aside=bad))

log()
log("  WHAT IS BEING TESTED, AND WHAT IS NOT.")
log()
log("  The claim is that the two coefficients are the same ESTIMAND. The")
log("  evidence for that is agreement between the point estimates. It is not")
log("  agreement between the intervals, because the two models make different")
log("  variance assumptions: their negative binomial absorbs overdispersion,")
log("  this conditional multinomial assumes Poisson variance and will report an")
log("  interval that is too narrow wherever the events are concentrated on a")
log("  few paths. So the primary check is on the point estimate, and the")
log("  interval comparison is reported as a diagnostic of that difference.")
log()
TOL = 0.15
for r in ROWS:
    d = abs(r["beta_ours"] - r["beta_published"])
    chk("%-18s point estimate agrees within %.2f" % (r["species"], TOL), d <= TOL,
        "ours %.3f vs theirs %.2f, difference %.3f"
        % (r["beta_ours"], r["beta_published"], d))
    r["abs_difference"] = round(d, 4)

log()
log("  Diagnostic: interval containment and relative width.")
for r in ROWS:
    ours_has_theirs = (r["ci_lo_ours"] <= r["beta_published"] <= r["ci_hi_ours"])
    theirs_has_ours = (r["ci_lo_published"] <= r["beta_ours"] <= r["ci_hi_published"])
    w_ours = r["ci_hi_ours"] - r["ci_lo_ours"]
    w_pub = r["ci_hi_published"] - r["ci_lo_published"]
    ratio = w_ours / w_pub if w_pub else float("nan")
    ev_per_codon = r["events"] / r["informative_codons"]
    log("    %-18s ours contains theirs %-3s   theirs contains ours %-3s"
        % (r["species"], "yes" if ours_has_theirs else "NO",
           "yes" if theirs_has_ours else "NO"))
    log("    %-18s width ratio %.2f   events per informative codon %.0f"
        % ("", ratio, ev_per_codon))
    r["ours_contains_theirs"] = int(ours_has_theirs)
    r["theirs_contains_ours"] = int(theirs_has_ours)
    r["ci_width_ratio"] = round(ratio, 3)
    r["events_per_codon"] = round(ev_per_codon, 1)

chk("their interval contains our estimate in every species",
    all(r["theirs_contains_ours"] for r in ROWS),
    "%d of %d" % (sum(r["theirs_contains_ours"] for r in ROWS), len(ROWS)))

narrow = [r for r in ROWS if r["ci_width_ratio"] < 0.5]
if narrow:
    log()
    log("  As predicted by the variance argument above, our interval is much too")
    log("  narrow in %s. That is the species with by far the most"
        % ", ".join(r["species"] for r in narrow))
    log("  concentrated data (%s events per informative codon), which is exactly"
        % ", ".join("%.0f" % r["events_per_codon"] for r in narrow))
    log("  where ignoring overdispersion bites. CONSEQUENCE FOR THE PAPER: quote")
    log("  the agreement of the point estimates as evidence of a shared")
    log("  estimand, and never quote an interval this conditional produces on")
    log("  count data. Our own atlas is not affected, because there each")
    log("  candidate allele is observed at most once per stratum, so there are")
    log("  no repeated events for overdispersion to act on.")

log()
log("  Also worth stating in the paper: in every one of their species beta is")
log("  near 1, and all three are clonal or non-recombining. Our sexual")
log("  arthropods are the comparison their Discussion asks for.")

log()
log("-" * 74)
log("WRITING")
log("-" * 74)
fields = ["species", "source_file", "informative_codons", "events",
          "events_per_codon", "beta_ours", "ci_lo_ours", "ci_hi_ours",
          "beta_published", "ci_lo_published", "ci_hi_published",
          "abs_difference", "ci_width_ratio",
          "ours_contains_theirs", "theirs_contains_ours",
          "rows_read", "rows_set_aside"]
with open(os.path.join(DIR_RES, "step9_03_benchmark_cano.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    for r in ROWS:
        w.writerow(r)
log("  wrote results/step9_03_benchmark_cano.csv")

log()
log("=" * 74)
if FAILURES:
    log("RESULT: %d CHECK(S) FAILED" % len(FAILURES))
    for f in FAILURES:
        log("   - %s" % f)
else:
    log("RESULT: every check passed")
log("=" * 74)

with open(os.path.join(DIR_LOG, "step9_03_benchmark_cano.log"), "w") as fh:
    fh.write("\n".join(_lines) + "\n")

sys.exit(1 if FAILURES else 0)
