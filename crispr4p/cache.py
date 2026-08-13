"""Versioned result-cache operations independent of design computation."""

import os
import pickle


def ensure_cache_directory(cache_directory):
    """Create the configured cache directory using legacy semantics."""
    if not os.path.isdir(cache_directory):
        os.makedirs(cache_directory)


def build_cache_path(
    cache_directory,
    version,
    n_mismatch,
    chromosome,
    start,
    end,
    systematic_name=None,
):
    """Return the historical gene or coordinate cache filename."""
    if systematic_name:
        basename = "%s_v%s_n%s.pickle" % (
            systematic_name,
            version,
            n_mismatch,
        )
    else:
        basename = "%s_%s_%s_v%s_n%s.pickle" % (
            chromosome,
            start,
            end,
            version,
            n_mismatch,
        )
    return os.path.join(cache_directory, basename)


def cache_exists(cache_path):
    return os.path.isfile(cache_path)


def load_cached_result(cache_path):
    """Load and normalize the historical four-item cache payload."""
    with open(cache_path, "rb") as cache_file:
        payload = pickle.load(cache_file)
    return _normalize_result(payload)


def save_cached_result(cache_path, result):
    """Best-effort cache write with silent cleanup on failure."""
    try:
        data = pickle.dumps(_normalize_result(result), protocol=-1)
        with open(cache_path, "wb") as cache_file:
            cache_file.write(data)
    except Exception:
        discard_cached_result(cache_path)


def discard_cached_result(cache_path):
    """Remove a cache file if one exists."""
    if os.path.isfile(cache_path):
        os.remove(cache_path)


def load_or_compute(cache_path, compute, is_cached=cache_exists):
    """Return a cached result or compute it with legacy failure behavior."""
    if not is_cached(cache_path):
        result = _normalize_result(compute())
        save_cached_result(cache_path, result)
        return result

    try:
        return load_cached_result(cache_path)
    except Exception:
        # A bad warm cache is removed and recomputed but intentionally is not
        # rewritten until the next request, matching the original run path.
        discard_cached_result(cache_path)
        return _normalize_result(compute())


def _normalize_result(result):
    table, hr_dna, checking_primers, guide_matches = result
    return table, hr_dna, checking_primers, guide_matches
