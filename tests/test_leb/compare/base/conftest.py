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
