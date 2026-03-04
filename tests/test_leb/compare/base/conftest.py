"""
Fixtures for base tier (de440s.bsp, 1850-2150) LEB comparison tests.

Provides date fixtures and tier-specific tolerances for base tier validation.
"""

from __future__ import annotations

from typing import Generator

import pytest

from tests.test_leb.compare.conftest import (
    CompareHelper,
    TierTolerances,
    generate_test_dates,
    year_to_jd,
)


# =============================================================================
# TIER CONFIG
# =============================================================================

TIER = "base"

_BASE_START = year_to_jd(1860)
_BASE_END = year_to_jd(2140)

TOLS_BASE = TierTolerances.for_tier(TIER)


# =============================================================================
# FIXTURES
# =============================================================================


# =============================================================================
# ASTEROID DATE FILTERING
# =============================================================================

# Asteroid SPK coverage ranges for the base tier (de440s.bsp).
#
# The LEB file stores the full tier range (1850-2150) for all bodies, but
# asteroid SPK files only cover ~1900-2100. Outside SPK coverage, the LEB
# generator used Keplerian orbital element fallback which produces
# catastrophically wrong data (7,000-14,000" errors). The LEB BodyEntry
# jd_start/jd_end reflect the full tier range, NOT the actual SPK coverage.
#
# Additionally, Chebyshev segments near the SPK boundary are contaminated:
# the fitting window may include dates outside SPK coverage, producing bad
# Chebyshev coefficients for the entire segment. This contamination extends
# well beyond the 32-day segment interval because the underlying Keplerian
# orbit diverges from the true orbit.
#
# The safe range below (1920-2080) provides ~20 years of margin from the
# ~1900/~2100 SPK boundaries, which is sufficient to avoid all contaminated
# segments. After LEB regeneration with proper SPK data (Phase 4), this
# range can be widened to match the full SPK coverage.
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
def base_dates_300() -> list[float]:
    """300 uniformly-spaced JDs across 1860-2140."""
    return generate_test_dates(300, _BASE_START, _BASE_END)


@pytest.fixture(scope="session")
def base_dates_150() -> list[float]:
    """150 uniformly-spaced JDs across 1860-2140."""
    return generate_test_dates(150, _BASE_START, _BASE_END)


@pytest.fixture(scope="session")
def base_dates_100() -> list[float]:
    """100 uniformly-spaced JDs across 1860-2140."""
    return generate_test_dates(100, _BASE_START, _BASE_END)


@pytest.fixture(scope="session")
def base_dates_50() -> list[float]:
    """50 uniformly-spaced JDs across 1860-2140."""
    return generate_test_dates(50, _BASE_START, _BASE_END)


@pytest.fixture
def compare(compare_base: CompareHelper) -> Generator[CompareHelper, None, None]:
    """Alias compare_base as compare for base tier tests."""
    yield compare_base
