"""
Caching utilities for libephemeris hot path optimization.

This module provides LRU caches for expensive calculations that are called
repeatedly with the same inputs. Key optimizations include:

1. Nutation caching - IAU 2006/2000A nutation calculations (via pyerfa)
2. Obliquity caching - True obliquity of the ecliptic (IAU 2006)
3. Timescale object caching - Skyfield timescale objects

These caches significantly improve performance when calculating multiple
planetary positions for the same Julian Day (common in chart calculations).

Precision:
    - Nutation: IAU 2006/2000A model via pyerfa (~0.01-0.05 mas precision)
    - Obliquity: IAU 2006 model via pyerfa (consistent across all code paths)

Thread Safety:
    The caches use functools.lru_cache which is thread-safe for reads but
    not for simultaneous writes with the same key. In practice, this is
    acceptable since we're caching pure functions with deterministic outputs.
"""

import math
from functools import lru_cache
from typing import Tuple

import erfa


# Cache size limits
# - NUTATION_CACHE_SIZE: Caches nutation angles (dpsi, deps) for Julian Days
#   Most chart calculations use 1-2 JDs, but time ranges may use more
# - OBLIQUITY_CACHE_SIZE: Caches obliquity values, same usage pattern
_NUTATION_CACHE_SIZE = 128
_OBLIQUITY_CACHE_SIZE = 128

# J2000.0 epoch in Julian Days
_J2000_JD = 2451545.0


@lru_cache(maxsize=_NUTATION_CACHE_SIZE)
def get_cached_nutation(jd_tt: float) -> Tuple[float, float]:
    """
    Get cached nutation angles (dpsi, deps) for a Julian Day.

    Uses the IAU 2006/2000A nutation model via pyerfa (erfa.nut06a), which
    is the most precise nutation model currently adopted by the IAU.
    This includes the IAU 2000A lunisolar+planetary nutation (1365 terms)
    with IAU 2006 J2 secular variation correction.

    Precision: ~0.01-0.05 milliarcsecond

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (dpsi_radians, deps_radians):
            - dpsi: Nutation in longitude (radians)
            - deps: Nutation in obliquity (radians)

    Performance:
        - Uncached: ~0.02ms per call
        - Cached hit: ~0.0001ms per call (200x faster)
    """
    # erfa.nut06a expects (jd1, jd2) where jd1 + jd2 = JD in TT
    # Using 2-part JD for maximum numerical precision (SOFA convention)
    dpsi, deps = erfa.nut06a(_J2000_JD, jd_tt - _J2000_JD)
    return float(dpsi), float(deps)


def get_nutation_degrees(jd_tt: float) -> Tuple[float, float]:
    """
    Get nutation in degrees (convenience wrapper).

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (dpsi_degrees, deps_degrees)
    """
    dpsi_rad, deps_rad = get_cached_nutation(jd_tt)
    return math.degrees(dpsi_rad), math.degrees(deps_rad)


@lru_cache(maxsize=_OBLIQUITY_CACHE_SIZE)
def get_cached_obliquity(jd_tt: float) -> Tuple[float, float]:
    """
    Get cached mean and true obliquity for a Julian Day.

    Uses IAU 2006 mean obliquity (via erfa.obl06) and IAU 2006/2000A
    nutation (via erfa.nut06a) for true obliquity, providing consistency
    with all other code paths in libephemeris.

    Mean obliquity formula (IAU 2006, Capitaine et al. 2003):
        eps0 = 84381.406″ - 46.836769″T - 0.0001831″T² + 0.00200340″T³
               - 0.000000576″T⁴ - 0.0000000434″T⁵

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (mean_obliquity, true_obliquity) in degrees
    """
    # Mean obliquity via IAU 2006 (erfa.obl06)
    # Returns mean obliquity in radians
    eps0_rad = erfa.obl06(_J2000_JD, jd_tt - _J2000_JD)
    eps0 = math.degrees(eps0_rad)

    # Nutation in obliquity from IAU 2006/2000A
    _, deps_rad = get_cached_nutation(jd_tt)
    deps = math.degrees(deps_rad)

    # True obliquity = mean + nutation in obliquity
    eps = eps0 + deps

    return eps0, eps


def get_mean_obliquity(jd_tt: float) -> float:
    """Get mean obliquity in degrees (IAU 2006)."""
    eps0, _ = get_cached_obliquity(jd_tt)
    return eps0


def get_true_obliquity(jd_tt: float) -> float:
    """Get true obliquity in degrees (IAU 2006 + nutation)."""
    _, eps = get_cached_obliquity(jd_tt)
    return eps


def clear_caches() -> None:
    """
    Clear all computation caches.

    Call this when:
    - Changing ephemeris files
    - After calling swe_close()
    - When memory needs to be freed
    """
    get_cached_nutation.cache_clear()
    get_cached_obliquity.cache_clear()


def get_cache_info() -> dict:
    """
    Get cache statistics for debugging and optimization.

    Returns:
        Dictionary with cache statistics for each cached function
    """
    return {
        "nutation": get_cached_nutation.cache_info()._asdict(),
        "obliquity": get_cached_obliquity.cache_info()._asdict(),
    }
