"""Result cache."""

import os
import pickle


def ensure_cache_dir(cache_directory):
    """Create the cache directory if needed."""
    if not os.path.isdir(cache_directory):
        os.makedirs(cache_directory)


def cache_path(
    cache_directory,
    version,
    n_mismatch,
    chromosome,
    start,
    end,
    systematic_name=None,
):
    """Return a gene or coordinate cache path."""
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


def load_cache(cache_path):
    """Load a cached design result."""
    with open(cache_path, "rb") as cache_file:
        payload = pickle.load(cache_file)
    return _normalize_result(payload)


def save_cache(cache_path, result):
    """Write a result, removing a partial file on failure."""
    try:
        data = pickle.dumps(_normalize_result(result), protocol=-1)
        with open(cache_path, "wb") as cache_file:
            cache_file.write(data)
    except Exception:
        drop_cache(cache_path)


def drop_cache(cache_path):
    """Remove a cached result."""
    if os.path.isfile(cache_path):
        os.remove(cache_path)


def get_cached(cache_path, compute, is_cached=cache_exists):
    """Load a result or compute it."""
    if not is_cached(cache_path):
        result = _normalize_result(compute())
        save_cache(cache_path, result)
        return result

    try:
        return load_cache(cache_path)
    except Exception:
        # Match the old behavior: recompute a bad cache without rewriting it.
        drop_cache(cache_path)
        return _normalize_result(compute())


def _normalize_result(result):
    table, hr_dna, checking_primers, guide_matches = result
    return table, hr_dna, checking_primers, guide_matches
