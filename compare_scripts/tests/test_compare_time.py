"""
Time Functions Comparison Tests.

Compares all time-related functions between pyswisseph and libephemeris:
- Julian Day conversions (julday, revjul)
- Delta T calculations
- UTC conversions
- Sidereal time
- Equation of time
"""

import pytest
import swisseph as swe
import libephemeris as ephem


# ============================================================================
# TOLERANCES
# ============================================================================

JULIAN_DAY_TOL = 1e-10
DELTA_T_TOL = 0.01  # seconds
SIDEREAL_TIME_TOL = 0.0001  # hours
TIME_EQU_TOL = 1e-6  # days


# ============================================================================
# TEST DATA
# ============================================================================

JULDAY_TEST_CASES = [
    # (year, month, day, hour, calendar_type, description)
    (2000, 1, 1, 12.0, 1, "J2000.0 Epoch"),
    (1900, 1, 1, 0.0, 1, "1900 Start"),
    (2024, 2, 29, 12.0, 1, "Leap Year 2024"),
    (1999, 12, 31, 23.999, 1, "End of 1999"),
    (2050, 6, 15, 6.5, 1, "Future Date"),
    (1582, 10, 15, 12.0, 1, "Gregorian Start"),
    (1582, 10, 4, 12.0, 0, "Julian Calendar Last Day"),
    (-4713, 1, 1, 12.0, 0, "JD Zero Point"),
    (1, 1, 1, 0.0, 0, "Year 1 AD Julian"),
]

REVJUL_TEST_CASES = [
    (2451545.0, "J2000.0"),
    (2460000.0, "JD 2460000"),
    (2415020.5, "1900-01-01"),
    (2488070.0, "2100-01-01"),
    (0.0, "JD Zero"),
]

DELTAT_TEST_DATES = [
    (1900, 1, 1, 0.0, "1900"),
    (1950, 1, 1, 0.0, "1950"),
    (2000, 1, 1, 12.0, "2000"),
    (2020, 1, 1, 0.0, "2020"),
    (2024, 6, 15, 12.0, "2024"),
]

UTC_TEST_CASES = [
    (2000, 1, 1, 12, 0, 0.0, "J2000.0"),
    (2024, 6, 15, 14, 30, 45.5, "Modern date"),
    (1999, 12, 31, 23, 59, 59.0, "End of 1999"),
    (2024, 2, 29, 0, 0, 0.0, "Leap day 2024"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestJulday:
    """Compare julday function between swisseph and libephemeris."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,cal_type,desc", JULDAY_TEST_CASES)
    def test_julday(self, year, month, day, hour, cal_type, desc):
        """Test julday returns identical values."""
        jd_swe = swe.julday(year, month, day, hour, cal_type)
        jd_py = ephem.julday(year, month, day, hour, cal_type)

        diff = abs(jd_swe - jd_py)

        assert diff < JULIAN_DAY_TOL, (
            f"{desc}: julday diff {diff:.15f} exceeds tolerance"
        )


class TestRevjul:
    """Compare revjul function between swisseph and libephemeris."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("jd,desc", REVJUL_TEST_CASES)
    def test_revjul(self, jd, desc):
        """Test revjul returns identical values."""
        y_swe, m_swe, d_swe, h_swe = swe.revjul(jd)
        y_py, m_py, d_py, h_py = ephem.revjul(jd)

        assert int(y_swe) == int(y_py), f"{desc}: year mismatch"
        assert int(m_swe) == int(m_py), f"{desc}: month mismatch"
        assert int(d_swe) == int(d_py), f"{desc}: day mismatch"

        hour_diff = abs(h_swe - h_py)
        assert hour_diff < 1e-8, f"{desc}: hour diff {hour_diff:.12f} exceeds tolerance"


class TestDeltat:
    """Compare Delta T calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", DELTAT_TEST_DATES)
    def test_deltat(self, year, month, day, hour, desc):
        """Test deltat returns values within tolerance."""
        jd = swe.julday(year, month, day, hour)

        dt_swe = swe.deltat(jd)
        dt_py = ephem.deltat(jd)

        diff_days = abs(dt_swe - dt_py)
        diff_seconds = diff_days * 86400

        assert diff_seconds < DELTA_T_TOL, (
            f"{desc}: Delta T diff {diff_seconds:.6f}s exceeds tolerance"
        )


class TestSiderealTime:
    """Compare sidereal time calculations."""

    TEST_DATES = [
        (2000, 1, 1, 12.0, "J2000"),
        (2024, 6, 15, 0.0, "Summer 2024"),
        (1950, 3, 21, 12.0, "Equinox 1950"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_sidtime(self, year, month, day, hour, desc):
        """Test sidtime returns identical values."""
        jd = swe.julday(year, month, day, hour)

        st_swe = swe.sidtime(jd)
        st_py = ephem.sidtime(jd)

        diff = abs(st_swe - st_py)
        # Handle wrap-around at 24 hours
        if diff > 12:
            diff = 24 - diff

        assert diff < SIDEREAL_TIME_TOL, (
            f"{desc}: sidtime diff {diff:.10f}h exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_sidtime0_with_obliquity(self):
        """Test sidtime0 with specified obliquity and nutation."""
        jd = 2451545.0  # J2000
        eps = 23.439291  # Mean obliquity at J2000
        nut = 0.00256  # Typical nutation in longitude

        st_swe = swe.sidtime0(jd, eps, nut)
        st_py = ephem.sidtime0(jd, eps, nut)

        diff = abs(st_swe - st_py)
        if diff > 12:
            diff = 24 - diff

        assert diff < SIDEREAL_TIME_TOL, f"sidtime0 diff {diff:.10f}h exceeds tolerance"


class TestUtcToJd:
    """Compare UTC to JD conversions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,minute,second,desc", UTC_TEST_CASES)
    def test_utc_to_jd(self, year, month, day, hour, minute, second, desc):
        """Test utc_to_jd returns identical values."""
        jd_et_swe, jd_ut_swe = swe.utc_to_jd(year, month, day, hour, minute, second, 1)
        jd_et_py, jd_ut_py = ephem.utc_to_jd(year, month, day, hour, minute, second, 1)

        diff_et = abs(jd_et_swe - jd_et_py)
        diff_ut = abs(jd_ut_swe - jd_ut_py)

        assert diff_et < JULIAN_DAY_TOL, (
            f"{desc}: JD_ET diff {diff_et:.15f} exceeds tolerance"
        )
        assert diff_ut < JULIAN_DAY_TOL, (
            f"{desc}: JD_UT diff {diff_ut:.15f} exceeds tolerance"
        )


class TestEquationOfTime:
    """Compare equation of time calculations."""

    EQU_TEST_CASES = [
        (2024, 2, 11, 12.0, "Feb 11 (max negative)"),
        (2024, 5, 14, 12.0, "May 14 (zero crossing)"),
        (2024, 7, 26, 12.0, "Jul 26 (zero crossing)"),
        (2024, 11, 3, 12.0, "Nov 3 (max positive)"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", EQU_TEST_CASES)
    def test_time_equ(self, year, month, day, hour, desc):
        """Test time_equ returns values within tolerance."""
        jd = swe.julday(year, month, day, hour)

        teq_swe = swe.time_equ(jd)
        teq_py = ephem.time_equ(jd)

        # time_equ returns tuple in some implementations
        e_swe = teq_swe[0] if isinstance(teq_swe, tuple) else teq_swe
        e_py = teq_py[0] if isinstance(teq_py, tuple) else teq_py

        diff = abs(e_swe - e_py)

        assert diff < TIME_EQU_TOL, (
            f"{desc}: time_equ diff {diff:.12f} days exceeds tolerance"
        )


class TestRoundtrip:
    """Test roundtrip conversions (julday -> revjul -> julday)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,cal_type,desc",
        JULDAY_TEST_CASES[:5],  # Use first 5 test cases
    )
    def test_julday_revjul_roundtrip(self, year, month, day, hour, cal_type, desc):
        """Test that julday -> revjul gives back original values."""
        # Forward
        jd = ephem.julday(year, month, day, hour, cal_type)

        # Reverse
        y2, m2, d2, h2 = ephem.revjul(jd, cal_type)

        # Compare
        assert int(y2) == year, f"{desc}: year mismatch in roundtrip"
        assert int(m2) == month, f"{desc}: month mismatch in roundtrip"
        assert int(d2) == day, f"{desc}: day mismatch in roundtrip"

        hour_diff = abs(h2 - hour)
        assert hour_diff < 1e-8, f"{desc}: hour diff {hour_diff:.12f} in roundtrip"
