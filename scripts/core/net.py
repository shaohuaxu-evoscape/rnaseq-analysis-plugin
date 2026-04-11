"""Shared network utilities with retry logic."""

import time
import urllib.request

from .logger import setup_logger

log = setup_logger("net")


def fetch_with_retry(url, timeout=120, retries=3, backoff=2.0):
    """Fetch a URL with retry and exponential backoff.

    Parameters
    ----------
    url : str
    timeout : int
        Request timeout in seconds.
    retries : int
        Number of retry attempts.
    backoff : float
        Base backoff multiplier in seconds.

    Returns
    -------
    str
        Response body as text.
    """
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Python/rnaseq-analysis"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8")
        except Exception as e:
            last_err = e
            if attempt < retries - 1:
                wait = backoff * (2 ** attempt)
                log.warning(f"  Fetch attempt {attempt + 1}/{retries} failed: {e}; retrying in {wait:.0f}s")
                time.sleep(wait)
    raise RuntimeError(f"Failed to fetch {url[:80]} after {retries} attempts: {last_err}")


def fetch_cached(url, cache_dir, cache_file, *, timeout=120):
    """Fetch a URL, caching the result to a local file.

    Returns the response body as text. If the cache file already exists,
    returns its contents without making a network request.
    """
    import os
    path = os.path.join(cache_dir, cache_file)
    if os.path.exists(path):
        log.info(f"  Using cached: {cache_file}")
        with open(path, "r") as f:
            return f.read()
    log.info(f"  Fetching: {url[:80]}...")
    text = fetch_with_retry(url, timeout=timeout)
    with open(path, "w") as f:
        f.write(text)
    return text
