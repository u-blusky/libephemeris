"""
Fixtures for extended tier (de441.bsp, -5000 to 5000) LEB comparison tests.

Provides date fixtures split into sub-ranges and tier-specific tolerances
for extended tier validation across 10,000 years of coverage.
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

TIER = "extended"

_EXT_START = year_to_jd(-4990)
_EXT_END = year_to_jd(4990)

# Sub-range boundaries
_ANCIENT_START = year_to_jd(-4990)
_ANCIENT_END = year_to_jd(-1000)

_MODERN_START = year_to_jd(-1000)
_MODERN_END = year_to_jd(3000)

_FUTURE_START = year_to_jd(3000)
_FUTURE_END = year_to_jd(4990)

TOLS_EXT = TierTolerances.for_tier(TIER)


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture(scope="session")
def ext_dates_500() -> list[float]:
    """500 uniformly-spaced JDs across -4990 to 4990."""
    return generate_test_dates(500, _EXT_START, _EXT_END)


@pytest.fixture(scope="session")
def ext_dates_200() -> list[float]:
    """200 uniformly-spaced JDs across -4990 to 4990."""
    return generate_test_dates(200, _EXT_START, _EXT_END)


@pytest.fixture(scope="session")
def ext_dates_100() -> list[float]:
    """100 uniformly-spaced JDs across -4990 to 4990."""
    return generate_test_dates(100, _EXT_START, _EXT_END)


@pytest.fixture(scope="session")
def ext_dates_50() -> list[float]:
    """50 uniformly-spaced JDs across -4990 to 4990."""
    return generate_test_dates(50, _EXT_START, _EXT_END)


@pytest.fixture(scope="session")
def ext_ancient_dates() -> list[float]:
    """150 dates in the ancient sub-range (-4990 to -1000)."""
    return generate_test_dates(150, _ANCIENT_START, _ANCIENT_END)


@pytest.fixture(scope="session")
def ext_modern_dates() -> list[float]:
    """200 dates in the modern sub-range (-1000 to 3000)."""
    return generate_test_dates(200, _MODERN_START, _MODERN_END)


@pytest.fixture(scope="session")
def ext_future_dates() -> list[float]:
    """150 dates in the future sub-range (3000 to 4990)."""
    return generate_test_dates(150, _FUTURE_START, _FUTURE_END)


@pytest.fixture(scope="session")
def ext_boundary_dates() -> list[float]:
    """40 dates concentrated at the boundaries of the extended range."""
    start_boundary = generate_test_dates(20, _EXT_START, year_to_jd(-4900), margin=5.0)
    end_boundary = generate_test_dates(20, year_to_jd(4900), _EXT_END, margin=5.0)
    return start_boundary + end_boundary


@pytest.fixture
def compare(compare_extended: CompareHelper) -> Generator[CompareHelper, None, None]:
    """Alias compare_extended as compare for extended tier tests."""
    yield compare_extended
