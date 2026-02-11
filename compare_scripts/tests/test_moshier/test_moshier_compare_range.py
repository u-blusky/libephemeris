"""
Moshier Ephemeris Range Comparison: pyswisseph (C) vs libephemeris.

Compares Moshier (SEFLG_MOSEPH) range limits between pyswisseph (C Swiss
Ephemeris) and libephemeris (Python reimplementation) to map the complete
range discrepancy.

Empirically determined ranges:
  - pyswisseph (C):  JD 625000.50 to JD 2818000.50 (~-2999 to ~+3003 CE)
  - libephemeris:    JD 625673.50 to JD 3182395.50 (~-2999 to ~+4000 CE)

The discrepancy is the *reverse* of the initial hypothesis: libephemeris has
a WIDER range than pyswisseph, not narrower. Dates ~+3004 to ~+4000 CE work
in libephemeris but are rejected by pyswisseph. This means libephemeris does
NOT regress on archaeological astronomy dates (Sumeria, Egypt) - both
libraries reject dates before ~-2999 CE equally.

Test matrix (15 dates):
  Category 1 (5): Inside both ranges (~-2999 to ~+3003) - both accept,
                   positions agree within tolerance.
  Category 2 (5): Inside libephemeris only (~+3004 to ~+4000) - libephemeris
                   accepts, pyswisseph rejects. Documents where libephemeris
                   EXCEEDS C capability.
  Category 3 (5): Outside both ranges - both reject.
"""

from __future__ import annotations

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MARS,
    SEFLG_MOSEPH,
)
from libephemeris.exceptions import (
    EphemerisRangeError,
    MOSHIER_JD_START,
    MOSHIER_JD_END,
)


# =============================================================================
# RANGE CONSTANTS (empirically determined from pyswisseph error messages)
# =============================================================================

# pyswisseph reports: "outside Moshier planet range 625000.50 .. 2818000.50"
C_MOSHIER_JD_START = 625000.50
C_MOSHIER_JD_END = 2818000.50


# =============================================================================
# DATE CONFIGURATIONS (15 dates total)
# =============================================================================

# Category 1: Inside both C and libephemeris Moshier ranges.
# Both libraries should accept these and return concordant positions.
# The overlap region is approximately JD 625673.5 to JD 2818000.5
# (~-2999 to ~+3003 CE).
DATES_INSIDE_BOTH = [
    (-2500, 6, 15, 12.0, "-2500 CE (early Bronze Age)"),
    (-1000, 3, 21, 12.0, "-1000 CE (Iron Age)"),
    (0, 1, 1, 12.0, "0 CE (1 BCE)"),
    (1000, 6, 15, 12.0, "1000 CE (Medieval)"),
    (2500, 1, 1, 12.0, "2500 CE (near C range end)"),
]

# Category 2: Inside libephemeris range but outside C range.
# libephemeris accepts these; pyswisseph raises swisseph.Error.
# This documents where libephemeris EXCEEDS C library capability.
# libephemeris range extends to JD 3182395.5 (~+4000 CE) while C stops
# at JD 2818000.5 (~+3003 CE).
DATES_LIBEPHEMERIS_ONLY = [
    (3100, 1, 1, 12.0, "3100 CE"),
    (3200, 6, 15, 12.0, "3200 CE"),
    (3500, 3, 21, 12.0, "3500 CE"),
    (3800, 9, 1, 12.0, "3800 CE"),
    (3999, 12, 31, 12.0, "3999 CE (near libephemeris end)"),
]

# Category 3: Outside both ranges - both reject.
# Includes archaeological dates (Sumeria, Egypt) AND far future dates.
# Documents that NEITHER library supports very ancient dates - there is
# no regression for archaeological astronomy migration.
DATES_OUTSIDE_BOTH = [
    (-5000, 1, 1, 12.0, "-5000 CE (Sumerian)"),
    (-4000, 6, 15, 12.0, "-4000 CE (pre-Dynastic Egypt)"),
    (5000, 1, 1, 12.0, "+5000 CE"),
    (7000, 6, 15, 12.0, "+7000 CE"),
    (10000, 1, 1, 12.0, "+10000 CE"),
]

# Tolerance for Moshier C-vs-Python position comparison.
# For dates far from J2000 (up to ~4500 years away), the VSOP87 C-vs-Python
# implementation differences grow beyond the 0.03° used for near-J2000 dates.
# Observed at -2000 CE: Sun ~0.033°, Mars ~0.056°.
# 0.1° (~6 arcmin) accommodates dates up to ~4500 years from J2000.
LONGITUDE_TOL = 0.1  # degrees (~360 arcsec)
LATITUDE_TOL = 0.1  # degrees


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# =============================================================================
# CATEGORY 1: Inside both ranges - both accept, positions agree
# =============================================================================


class TestDatesInsideBothRanges:
    """Dates inside both libephemeris and C Moshier ranges.

    Both libraries should accept these dates. Sun and Mars positions from both
    implementations should agree within tolerance, validating that the Python
    Moshier reimplementation is consistent with the C original across the full
    overlapping range.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_INSIDE_BOTH)
    def test_both_accept_sun(self, year, month, day, hour, date_desc):
        """Both accept and Sun positions agree within tolerance."""
        jd = swe.julday(year, month, day, hour)

        pos_c, _ = swe.calc_ut(jd, swe.SUN, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)

        assert len(pos_c) == 6
        assert len(pos_py) == 6

        diff_lon = angular_diff(pos_c[0], pos_py[0])
        diff_lat = abs(pos_c[1] - pos_py[1])

        assert diff_lon < LONGITUDE_TOL, (
            f"Sun lon at {date_desc}: "
            f"C={pos_c[0]:.6f}, Py={pos_py[0]:.6f}, "
            f"diff={diff_lon:.6f} > tol={LONGITUDE_TOL}"
        )
        assert diff_lat < LATITUDE_TOL, (
            f"Sun lat at {date_desc}: "
            f"C={pos_c[1]:.6f}, Py={pos_py[1]:.6f}, "
            f"diff={diff_lat:.6f} > tol={LATITUDE_TOL}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_INSIDE_BOTH)
    def test_both_accept_mars(self, year, month, day, hour, date_desc):
        """Both accept and Mars positions agree within tolerance."""
        jd = swe.julday(year, month, day, hour)

        pos_c, _ = swe.calc_ut(jd, swe.MARS, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MARS, SEFLG_MOSEPH)

        diff_lon = angular_diff(pos_c[0], pos_py[0])

        assert diff_lon < LONGITUDE_TOL, (
            f"Mars lon at {date_desc}: "
            f"C={pos_c[0]:.6f}, Py={pos_py[0]:.6f}, "
            f"diff={diff_lon:.6f} > tol={LONGITUDE_TOL}"
        )


# =============================================================================
# CATEGORY 2: Inside libephemeris only - libephemeris exceeds C range
# =============================================================================


class TestDatesLibephemerisOnly:
    """Dates inside libephemeris range but outside C Moshier range.

    libephemeris extends Moshier coverage to ~+4000 CE (JD 3182395.5) while
    the C implementation stops at ~+3003 CE (JD 2818000.5). This is a region
    where libephemeris EXCEEDS pyswisseph capability.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_LIBEPHEMERIS_ONLY)
    def test_libephemeris_accepts(self, year, month, day, hour, date_desc):
        """libephemeris accepts and returns valid positions."""
        jd = ephem.swe_julday(year, month, day, hour)

        pos, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)

        assert len(pos) == 6, (
            f"libephemeris at {date_desc}: expected 6 elements, got {len(pos)}"
        )
        assert 0 <= pos[0] < 360, (
            f"libephemeris Sun lon at {date_desc}: {pos[0]:.6f} not in [0, 360)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_LIBEPHEMERIS_ONLY)
    def test_pyswisseph_rejects(self, year, month, day, hour, date_desc):
        """pyswisseph raises error for dates beyond its Moshier range.

        Documents that the C Swiss Ephemeris has a NARROWER Moshier range
        than libephemeris for future dates.
        """
        jd = swe.julday(year, month, day, hour)

        with pytest.raises(swe.Error, match=r"outside Moshier planet range"):
            swe.calc_ut(jd, swe.SUN, swe.FLG_MOSEPH)


# =============================================================================
# CATEGORY 3: Outside both ranges - both reject
# =============================================================================


class TestDatesOutsideBothRanges:
    """Dates outside both libephemeris and C Moshier ranges.

    Both libraries should reject these. This includes the archaeological dates
    (Sumeria ~-5000, Egypt ~-4000) that were hypothesized to be a migration
    regression - empirical testing shows NEITHER library supports them.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_OUTSIDE_BOTH)
    def test_libephemeris_rejects(self, year, month, day, hour, date_desc):
        """libephemeris raises EphemerisRangeError."""
        jd = swe.julday(year, month, day, hour)

        with pytest.raises(EphemerisRangeError) as exc_info:
            ephem.swe_calc_ut(jd, SE_SUN, SEFLG_MOSEPH)

        err = exc_info.value
        assert err.requested_jd == jd
        assert "Moshier" in err.message
        assert err.start_jd == MOSHIER_JD_START
        assert err.end_jd == MOSHIER_JD_END

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", DATES_OUTSIDE_BOTH)
    def test_pyswisseph_rejects(self, year, month, day, hour, date_desc):
        """pyswisseph also rejects - no migration regression."""
        jd = swe.julday(year, month, day, hour)

        with pytest.raises(swe.Error, match=r"outside Moshier planet range"):
            swe.calc_ut(jd, swe.SUN, swe.FLG_MOSEPH)


# =============================================================================
# Range boundary documentation
# =============================================================================


class TestRangeConstants:
    """Document the empirically determined range boundaries."""

    @pytest.mark.comparison
    def test_c_range_boundaries(self):
        """C Moshier range is JD 625000.50 to 2818000.50.

        Note: The C range check is applied to the TT Julian Day (after
        delta-T correction), so UT values near 625000.50 may still pass.
        We use large offsets (+-1000) to avoid boundary ambiguity.
        """
        # Well inside C range - should succeed
        pos, _ = swe.calc_ut(C_MOSHIER_JD_START + 1000, swe.SUN, swe.FLG_MOSEPH)
        assert len(pos) == 6

        # Well outside C range - should fail (1000 days offset exceeds any
        # delta-T uncertainty at the range boundaries)
        with pytest.raises(swe.Error):
            swe.calc_ut(C_MOSHIER_JD_START - 1000, swe.SUN, swe.FLG_MOSEPH)
        with pytest.raises(swe.Error):
            swe.calc_ut(C_MOSHIER_JD_END + 1000, swe.SUN, swe.FLG_MOSEPH)

    @pytest.mark.comparison
    def test_libephemeris_range_boundaries(self):
        """libephemeris Moshier range is JD 625673.5 to 3182395.5."""
        assert MOSHIER_JD_START == 625673.5
        assert MOSHIER_JD_END == 3182395.5

        # Just inside - should succeed
        pos, _ = ephem.swe_calc_ut(MOSHIER_JD_START, SE_SUN, SEFLG_MOSEPH)
        assert len(pos) == 6
        pos, _ = ephem.swe_calc_ut(MOSHIER_JD_END, SE_SUN, SEFLG_MOSEPH)
        assert len(pos) == 6

        # Just outside - should fail
        with pytest.raises(EphemerisRangeError):
            ephem.swe_calc_ut(MOSHIER_JD_START - 1.0, SE_SUN, SEFLG_MOSEPH)
        with pytest.raises(EphemerisRangeError):
            ephem.swe_calc_ut(MOSHIER_JD_END + 1.0, SE_SUN, SEFLG_MOSEPH)

    @pytest.mark.comparison
    def test_libephemeris_range_is_wider(self):
        """libephemeris extends ~1000 years beyond C range at the upper end.

        C range ends at JD 2818000.50 (~+3003 CE).
        libephemeris extends to JD 3182395.5 (~+4000 CE).
        The lower bound is nearly identical (~-2999 CE in both).
        """
        assert MOSHIER_JD_END > C_MOSHIER_JD_END
        range_extension_days = MOSHIER_JD_END - C_MOSHIER_JD_END
        range_extension_years = range_extension_days / 365.25
        # libephemeris extends ~997 years beyond C range
        assert range_extension_years > 900
