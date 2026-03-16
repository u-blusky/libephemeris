"""
Fixtures for medium tier (de440.bsp, 1550-2650) LEB comparison tests.

Provides date fixtures and tier-specific tolerances for medium tier validation.

Note: The parent conftest's ``compare`` fixture already targets medium tier,
so we do NOT override it here (unlike base/ which redirects to compare_base).
"""

from __future__ import annotations

import pytest

from tests.test_leb.compare.conftest import (
    TierTolerances,
    generate_test_dates,
    year_to_jd,
)


# =============================================================================
# TIER CONFIG
# =============================================================================

TIER = "medium"

_MEDIUM_START = year_to_jd(1560)
_MEDIUM_END = year_to_jd(2640)

TOLS_MEDIUM = TierTolerances.for_tier(TIER)


# =============================================================================
# ASTEROID DATE FILTERING
# =============================================================================

# Asteroid SPK coverage ranges for the medium tier (de440.bsp).
#
# The LEB file stores the full tier range (1550-2650) for all bodies, but
# asteroid SPK files only cover ~1900-2100. Outside SPK coverage, the LEB
# generator used Keplerian orbital element fallback which produces
# catastrophically wrong data (7,000-97,000" errors).
#
# The safe range below (1920-2080) provides ~20 years of margin from the
# ~1900/~2100 SPK boundaries, which is sufficient to avoid all contaminated
# segments.
_ASTEROID_SPK_COVERAGE: dict[int, tuple[float, float]] = {
    15: (year_to_jd(1920), year_to_jd(2080)),  # Chiron
    17: (year_to_jd(1920), year_to_jd(2080)),  # Ceres
    18: (year_to_jd(1920), year_to_jd(2080)),  # Pallas
    19: (year_to_jd(1920), year_to_jd(2080)),  # Juno
    20: (year_to_jd(1920), year_to_jd(2080)),  # Vesta
}


def filter_dates_for_body(
    dates: list[float], leb_path: str, body_id: int, margin: float = 30.0
) -> list[float]:
    """Filter test dates to only include those within a body's SPK coverage.

    For non-asteroid bodies, returns dates unchanged (no filtering needed).
    For asteroids, uses hardcoded SPK coverage ranges to exclude dates where
    the LEB file contains Keplerian-fallback data.

    Args:
        dates: list of JD values to filter.
        leb_path: path to the LEB file (unused, kept for API compatibility).
        body_id: SE body ID.
        margin: days of margin inside the boundaries.

    Returns:
        Filtered list of JDs within [jd_start + margin, jd_end - margin].
    """
    if body_id not in _ASTEROID_SPK_COVERAGE:
        return dates
    jd_start, jd_end = _ASTEROID_SPK_COVERAGE[body_id]
    return [jd for jd in dates if jd_start + margin <= jd <= jd_end - margin]


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture(scope="session")
def medium_dates_300() -> list[float]:
    """300 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(300, _MEDIUM_START, _MEDIUM_END)


@pytest.fixture(scope="session")
def medium_dates_200() -> list[float]:
    """200 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(200, _MEDIUM_START, _MEDIUM_END)


@pytest.fixture(scope="session")
def medium_dates_150() -> list[float]:
    """150 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(150, _MEDIUM_START, _MEDIUM_END)


@pytest.fixture(scope="session")
def medium_dates_100() -> list[float]:
    """100 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(100, _MEDIUM_START, _MEDIUM_END)


@pytest.fixture(scope="session")
def medium_dates_50() -> list[float]:
    """50 uniformly-spaced JDs across 1560-2640."""
    return generate_test_dates(50, _MEDIUM_START, _MEDIUM_END)
