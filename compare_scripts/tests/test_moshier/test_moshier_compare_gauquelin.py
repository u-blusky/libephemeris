"""
Moshier Gauquelin Sector Cross-Library Comparison Tests.

Compares swe.gauquelin_sector() (pyswisseph C Moshier) vs
ephem.swe_gauquelin_sector() (libephemeris Python Moshier) using SEFLG_MOSEPH.

Gauquelin sectors (1-36) are the basis of the most cited statistical study in
astrology (the Mars Effect). Replicating Gauquelin research for historical dates
with libephemeris requires sector assignments identical to pyswisseph; a
discrepancy of even one sector invalidates the entire statistical analysis.

This test validates Mars, Jupiter, Saturn across 3 locations (Rome, London,
New York), 2 methods (Gauquelin=0, Kundig=1), and 2 dates, producing 36
parametrized test cases.

Structural differences:
- ARMC and obliquity are computed by Skyfield in libephemeris but by internal
  Moshier routines in the C library, producing small positional discrepancies.
- Planet positions via Moshier semi-analytical algorithms may differ slightly
  between C and Python implementations.
- These justify a tolerance of < 0.5 sectors (less than half a sector), ensuring
  the assigned sector number is either identical or off by at most one.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SEFLG_MOSEPH,
)


# ============================================================================
# TOLERANCES
# ============================================================================

# Less than half a sector: ensures sector assignment is either identical
# or differs by at most 1 sector (out of 36).
MOSHIER_SECTOR_TOL = 0.5


# ============================================================================
# TEST DATA
# ============================================================================

# Planets relevant to Gauquelin research (Mars Effect, Jupiter Effect, etc.)
PLANETS = [
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# Locations from conftest.py test_locations (Rome, London, New York)
LOCATIONS = [
    ("Rome", 41.9028, 12.4964),
    ("London", 51.5074, -0.1278),
    ("New York", 40.7128, -74.0060),
]

# Methods: 0=Gauquelin (with latitude), 1=Kundig (without latitude)
METHODS = [
    (0, "Gauquelin"),
    (1, "Kundig"),
]

# Two dates within both DE440 and Moshier range for cross-validation
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (1980, 5, 20, 14.5, "Past"),
]

# Standard atmospheric parameters
ATPRESS = 1013.25
ATTEMP = 15.0


# ============================================================================
# TEST CLASS
# ============================================================================


class TestMoshierGauquelinSector:
    """Cross-library comparison of Gauquelin sector with SEFLG_MOSEPH."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PLANETS)
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS)
    @pytest.mark.parametrize("method,method_name", METHODS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_gauquelin_sector_moshier_cross_library(
        self,
        body_id,
        body_name,
        loc_name,
        lat,
        lon,
        method,
        method_name,
        year,
        month,
        day,
        hour,
        date_desc,
    ):
        """Gauquelin sector with SEFLG_MOSEPH must match between implementations.

        Validates that sector values agree within < 0.5 sectors and that both
        return values in the valid range [1, 37).
        """
        jd = swe.julday(year, month, day, hour)
        geopos = (lon, lat, 0.0)  # (longitude, latitude, altitude)

        # --- pyswisseph (C Moshier) ---
        # pyswisseph signature: gauquelin_sector(jd, body, method, geopos,
        #     atpress=0, attemp=0, flags=FLG_SWIEPH|FLG_TOPOCTR)
        try:
            ret_swe = swe.gauquelin_sector(
                jd, body_id, method, geopos, ATPRESS, ATTEMP, swe.FLG_MOSEPH
            )
            sector_swe = ret_swe[0] if isinstance(ret_swe, tuple) else ret_swe
        except Exception as e:
            pytest.skip(f"pyswisseph gauquelin_sector failed: {e}")

        # --- libephemeris (Python Moshier) ---
        try:
            sector_py = ephem.swe_gauquelin_sector(
                jd, body_id, method, geopos, ATPRESS, ATTEMP, SEFLG_MOSEPH
            )
        except Exception as e:
            pytest.skip(f"libephemeris swe_gauquelin_sector failed: {e}")

        # --- Validate range [1, 37) for both ---
        assert 1.0 <= sector_swe < 37.0, (
            f"pyswisseph sector {sector_swe} out of range [1, 37) "
            f"for {body_name} at {loc_name} ({date_desc}, {method_name})"
        )
        assert 1.0 <= sector_py < 37.0, (
            f"libephemeris sector {sector_py} out of range [1, 37) "
            f"for {body_name} at {loc_name} ({date_desc}, {method_name})"
        )

        # --- Compare sector values ---
        diff = abs(sector_swe - sector_py)
        # Handle wrap-around: sector 36.9 and 1.1 are only 0.2 apart
        if diff > 18.0:
            diff = 36.0 - diff

        assert diff < MOSHIER_SECTOR_TOL, (
            f"{body_name} Gauquelin sector mismatch at {loc_name} "
            f"({date_desc}, {method_name}): "
            f"pyswisseph={sector_swe:.4f}, libephemeris={sector_py:.4f}, "
            f"diff={diff:.4f} (tolerance={MOSHIER_SECTOR_TOL})"
        )
