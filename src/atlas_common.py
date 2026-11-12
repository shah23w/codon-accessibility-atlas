#!/usr/bin/env python3
"""
atlas_common.py: shared infrastructure for the codon-accessibility atlas.

Everything that touches the network, the filesystem layout, or logging lives here
so that each numbered pipeline step stays short and standalone-runnable.

Design rules enforced here:
  * every remote call retries with exponential backoff and raises on final failure
    (a failed retrieval must never silently become an empty/placeholder record);
  * every raw payload is cached to data/raw so a re-run is reproducible offline
    and the exact bytes used in the analysis can be shipped with the paper;
  * NCBI E-utilities rate limits are respected (3 req/s anonymous, 10 req/s with
    an API key supplied via the NCBI_API_KEY environment variable).
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import datetime
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

# --------------------------------------------------------------------------- #
# Project layout
# --------------------------------------------------------------------------- #
# The repository root, inferred from this file: src/atlas_common.py -> src -> root.
# It must NOT default to a particular machine's path. This file is shipped in
# the public deposit, where that default is wrong and silently writes outputs
# somewhere the reader did not ask for.
#
# CODON_ATLAS_ROOT overrides the inference, but ONLY when the directory above
# this file is not itself a project root. A tree that has data/ and results/
# beside src/ is self-contained and wins, because an environment variable left
# set in the shell otherwise reaches across into it: run inside the public
# deposit with CODON_ATLAS_ROOT pointing at the working repository, the five
# figure scripts that take their paths from this module wrote their PDFs back
# into the working repository, the deposit's own verifier found nothing, and
# the deposit reported five figures missing while having drawn them correctly.
# A reader with any such variable set would have seen the same.
_HERE = os.path.dirname(os.path.abspath(__file__))
_INFERRED = os.path.dirname(_HERE)
_SELF_CONTAINED = any(os.path.isdir(os.path.join(_INFERRED, d))
                      for d in ("data", "results"))
PROJ = (_INFERRED if _SELF_CONTAINED
        else (os.environ.get("CODON_ATLAS_ROOT") or _INFERRED))

DIR_SRC = os.path.join(PROJ, "src")
DIR_RAW = os.path.join(PROJ, "data", "raw")
DIR_INTERIM = os.path.join(PROJ, "data", "interim")
DIR_PROC = os.path.join(PROJ, "data", "processed")
DIR_REFS = os.path.join(PROJ, "data", "refs")
DIR_RESULTS = os.path.join(PROJ, "results")
DIR_FIGS = os.path.join(PROJ, "results", "figures")
DIR_LOGS = os.path.join(PROJ, "logs")
DIR_DOCS = os.path.join(PROJ, "docs")
DIR_MS = os.path.join(PROJ, "manuscript")
DIR_GENOMES = os.path.join(DIR_RAW, "genomes")
DIR_ENTREZ = os.path.join(DIR_RAW, "entrez")

ALL_DIRS = [DIR_SRC, DIR_RAW, DIR_INTERIM, DIR_PROC, DIR_REFS, DIR_RESULTS,
            DIR_FIGS, DIR_LOGS, DIR_DOCS, DIR_MS, DIR_GENOMES, DIR_ENTREZ]


def ensure_dirs() -> None:
    for d in ALL_DIRS:
        os.makedirs(d, exist_ok=True)


# --------------------------------------------------------------------------- #
# Contact / rate-limit configuration
# --------------------------------------------------------------------------- #
CONTACT_EMAIL = os.environ.get("NCBI_EMAIL", "sakhawatshahkhaksar@gmail.com")
API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
TOOL_NAME = "codon-accessibility-atlas"

# NCBI: 3 requests/second anonymous, 10/second with a key. Stay just under.
_MIN_INTERVAL = 0.11 if API_KEY else 0.36
_last_call = [0.0]

USER_AGENT = f"{TOOL_NAME}/1.0 ({CONTACT_EMAIL})"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"


def _throttle() -> None:
    dt = time.time() - _last_call[0]
    if dt < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - dt)
    _last_call[0] = time.time()


# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #
class Tee:
    """Minimal logger writing to stdout and a per-step logfile."""

    def __init__(self, step: str):
        ensure_dirs()
        self.path = os.path.join(DIR_LOGS, f"{step}.log")
        self.fh = open(self.path, "a")
        self.t0 = time.time()
        self("=" * 72)
        self(f"RUN {step}  started {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
        self("=" * 72)

    def __call__(self, msg: str = "") -> None:
        line = f"[{time.time()-self.t0:7.1f}s] {msg}" if msg else ""
        print(line, flush=True)
        self.fh.write(line + "\n")
        self.fh.flush()

    def close(self) -> None:
        self(f"finished in {time.time()-self.t0:.1f}s")
        self.fh.close()


# --------------------------------------------------------------------------- #
# HTTP with retry + on-disk cache
# --------------------------------------------------------------------------- #
class RetrievalError(RuntimeError):
    """Raised when a remote resource could not be retrieved after all retries.

    Callers must let this propagate or explicitly record the item as MISSING.
    Substituting a placeholder value is forbidden by the project data policy.
    """


def http_get(url: str, timeout: int = 180, retries: int = 5,
             binary: bool = False, headers: Optional[dict] = None):
    """GET with exponential backoff. Raises RetrievalError on final failure."""
    hdrs = {"User-Agent": USER_AGENT}
    if headers:
        hdrs.update(headers)
    last = None
    for attempt in range(retries):
        try:
            _throttle()
            req = urllib.request.Request(url, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
            return data if binary else data.decode("utf-8", errors="replace")
        except Exception as e:  # noqa: BLE001 - we genuinely retry anything
            last = e
            code = getattr(e, "code", None)
            # 400-level errors other than 429 will not fix themselves
            if code is not None and 400 <= code < 500 and code != 429:
                raise RetrievalError(f"HTTP {code} for {url}: {e}") from e
            sleep = min(60, 2 ** attempt + 1)
            time.sleep(sleep)
    raise RetrievalError(f"failed after {retries} attempts: {url}\n  last error: {last}")


def cached_get(url: str, cache_name: str, subdir: str = "misc",
               binary: bool = False, timeout: int = 180, refresh: bool = False):
    """http_get with an on-disk cache under data/raw/<subdir>/."""
    d = os.path.join(DIR_RAW, subdir)
    os.makedirs(d, exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in cache_name)[:180]
    path = os.path.join(d, safe)
    if os.path.exists(path) and not refresh and os.path.getsize(path) > 0:
        mode = "rb" if binary else "r"
        with open(path, mode) as fh:
            return fh.read()
    data = http_get(url, timeout=timeout, binary=binary)
    mode = "wb" if binary else "w"
    with open(path, mode) as fh:
        fh.write(data)
    return data


# --------------------------------------------------------------------------- #
# NCBI E-utilities
# --------------------------------------------------------------------------- #
def _eutils_params(**kw) -> str:
    kw.setdefault("tool", TOOL_NAME)
    kw.setdefault("email", CONTACT_EMAIL)
    if API_KEY:
        kw.setdefault("api_key", API_KEY)
    return urllib.parse.urlencode({k: v for k, v in kw.items() if v is not None})


def esearch(db: str, term: str, retmax: int = 0, retstart: int = 0,
            usehistory: bool = False) -> dict:
    q = _eutils_params(db=db, term=term, retmax=retmax, retstart=retstart,
                       retmode="json", usehistory="y" if usehistory else None)
    return json.loads(http_get(EUTILS + "esearch.fcgi?" + q))["esearchresult"]


def esearch_all_ids(db: str, term: str, log=None, cap: Optional[int] = None,
                    page: int = 5000) -> List[str]:
    """Page through esearch and return every UID (up to `cap`)."""
    first = esearch(db, term, retmax=0)
    total = int(first.get("count", 0))
    if cap:
        total = min(total, cap)
    ids: List[str] = []
    start = 0
    while start < total:
        r = esearch(db, term, retmax=min(page, total - start), retstart=start)
        got = r.get("idlist", []) or []
        if not got:
            break
        ids.extend(got)
        start += len(got)
        if log:
            log(f"    esearch {db}: {len(ids)}/{total}")
        if len(got) < min(page, total - start + len(got)):
            pass
    return ids[:cap] if cap else ids


def efetch(db: str, ids: Iterable[str], rettype: str = "fasta",
           retmode: str = "text", batch: int = 200, log=None) -> str:
    """Batched efetch. Returns the concatenated payload."""
    ids = list(ids)
    out = []
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        q = _eutils_params(db=db, id=",".join(chunk), rettype=rettype, retmode=retmode)
        out.append(http_get(EUTILS + "efetch.fcgi?" + q, timeout=300))
        if log and (i // batch) % 5 == 0:
            log(f"    efetch {db}: {min(i+batch, len(ids))}/{len(ids)}")
    return "".join(out)


def esummary(db: str, ids: Iterable[str], batch: int = 200, log=None) -> Dict[str, dict]:
    ids = list(ids)
    res: Dict[str, dict] = {}
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        q = _eutils_params(db=db, id=",".join(chunk), retmode="json")
        d = json.loads(http_get(EUTILS + "esummary.fcgi?" + q, timeout=300))
        result = d.get("result", {})
        for uid in result.get("uids", []):
            res[uid] = result[uid]
        if log and (i // batch) % 5 == 0:
            log(f"    esummary {db}: {min(i+batch, len(ids))}/{len(ids)}")
    return res


# --------------------------------------------------------------------------- #
# NCBI Datasets
# --------------------------------------------------------------------------- #
def datasets_json(path: str, **params) -> dict:
    url = f"{DATASETS}/{path.lstrip('/')}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    return json.loads(http_get(url, timeout=300))


def datasets_paged(path: str, key: str = "reports", page_size: int = 1000,
                   log=None, **params) -> List[dict]:
    """Page through any Datasets endpoint that uses next_page_token."""
    out: List[dict] = []
    token = None
    page = 0
    while True:
        p = dict(params)
        p["page_size"] = page_size
        if token:
            p["page_token"] = token
        d = datasets_json(path, **p)
        rows = d.get(key, []) or []
        out.extend(rows)
        page += 1
        if log:
            log(f"    page {page}: +{len(rows)} (total {len(out)})")
        token = d.get("next_page_token")
        if not token or not rows:
            break
    return out


# --------------------------------------------------------------------------- #
# FASTA
# --------------------------------------------------------------------------- #
def parse_fasta(handle) -> Iterator[Tuple[str, str]]:
    """Stream (header, sequence) pairs from a text handle or string."""
    if isinstance(handle, (str, bytes)):
        if isinstance(handle, bytes):
            handle = handle.decode("utf-8", errors="replace")
        handle = io.StringIO(handle)
    name, seq = None, []
    for line in handle:
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:], []
        elif line:
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def write_fasta(path: str, records: Iterable[Tuple[str, str]], width: int = 60) -> int:
    n = 0
    with open(path, "w") as fh:
        for h, s in records:
            fh.write(f">{h}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")
            n += 1
    return n


def read_fasta_file(path: str) -> Dict[str, str]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return {h: s for h, s in parse_fasta(fh)}


# --------------------------------------------------------------------------- #
# Misc
# --------------------------------------------------------------------------- #
def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save_json(obj, path: str, indent: int = 1) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=indent, default=str)
    os.replace(tmp, path)


def load_json(path: str):
    with open(path) as fh:
        return json.load(fh)


def species_key(name: str) -> str:
    """Normalise an organism string to 'Genus species'."""
    if not name:
        return ""
    n = name.strip().replace("[", "").replace("]", "")
    for cut in (" str. ", " strain ", " isolate ", " subsp. ", " var. ",
                " cv. ", " serotype ", " biotype "):
        if cut in n:
            n = n.split(cut)[0].strip()
            break
    toks = n.split()
    if len(toks) >= 2:
        return f"{toks[0]} {toks[1]}"
    return n


if __name__ == "__main__":
    ensure_dirs()
    print(f"project root : {PROJ}")
    print(f"contact email: {CONTACT_EMAIL}")
    print(f"NCBI API key : {'set (10 req/s)' if API_KEY else 'not set (3 req/s)'}")
    print("directories ready:")
    for d in ALL_DIRS:
        print(f"  {d}")


def cluster_bootstrap_ci(labels, accessible, probs, n_boot: int = 5000,
                         seed: int = 20260731, conf: float = 0.95):
    """Percentile interval for observed/expected, resampling whole clusters.

    The three arguments are aligned, one entry per event: `labels` is the
    cluster the event belongs to, `accessible` is 1 where the evolved allele
    was one step from the species' codon, and `probs` is that event's exact
    null probability of being accessible.

    Clusters, not events, are drawn with replacement. A species contributing
    many events is therefore not treated as that many independent draws, which
    is the whole point: the enrichment's nominal exact P value conditions on
    independence between events and the robustness strata inherit that
    assumption. For the collapsed strata the caller passes the collapse unit
    (genus, family, order) as the label instead of the species, so that the
    resampled unit matches the unit the stratum is defined on.

    Returns (lo, hi) rounded to four decimals, or ("", "") where the interval
    is undefined - fewer than two clusters, or an expected count of zero in
    every resample.
    """
    import random as _random

    by = {}
    for lab, acc, p in zip(labels, accessible, probs):
        by.setdefault(lab, []).append((1.0 if acc else 0.0, float(p)))
    keys = sorted(by)
    if len(keys) < 2:
        return "", ""
    rng = _random.Random(seed)
    ratios = []
    for _ in range(n_boot):
        o = e = 0.0
        for k in (keys[rng.randrange(len(keys))] for _ in keys):
            for acc, p in by[k]:
                o += acc
                e += p
        if e > 0:
            ratios.append(o / e)
    if not ratios:
        return "", ""
    ratios.sort()
    a = (1.0 - conf) / 2.0
    lo = ratios[max(0, min(len(ratios) - 1, int(a * len(ratios))))]
    hi = ratios[max(0, min(len(ratios) - 1, int((1.0 - a) * len(ratios))))]
    return round(lo, 4), round(hi, 4)
