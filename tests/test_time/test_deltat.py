"""
Comprehensive tests for Delta T calculation (swe_deltat).

Delta T = TT - UT (Terrestrial Time minus Universal Time)
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestDeltatBasicValues:
    """Test Delta T returns reasonable values for known epochs."""

    @pytest.mark.unit
    def test_deltat_j2000(self):
        """Delta T at J2000 should be around 63-64 seconds."""
        jd = 2451545.0  # J2000
        dt = ephem.swe_deltat(jd)
        # Convert to seconds (dt is in days)
        dt_seconds = dt * 86400
        assert 60 < dt_seconds < 70, f"Delta T at J2000 = {dt_seconds}s, expected ~63s"

    @pytest.mark.unit
    def test_deltat_1950(self):
        """Delta T around 1950 should be ~29 seconds."""
        jd = ephem.swe_julday(1950, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd)
        dt_seconds = dt * 86400
        assert 25 < dt_seconds < 35, f"Delta T at 1950 = {dt_seconds}s, expected ~29s"

    @pytest.mark.unit
    def test_deltat_1900(self):
        """Delta T around 1900 should be small negative or near zero."""
        jd = ephem.swe_julday(1900, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd)
        dt_seconds = dt * 86400
        # Around 1900, Delta T was about -3 to +5 seconds
        assert -10 < dt_seconds < 10, f"Delta T at 1900 = {dt_seconds}s"

    @pytest.mark.unit
    def test_deltat_2020(self):
        """Delta T around 2020 should be ~69 seconds."""
        jd = ephem.swe_julday(2020, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd)
        dt_seconds = dt * 86400
        assert 65 < dt_seconds < 75, f"Delta T at 2020 = {dt_seconds}s, expected ~69s"


class TestDeltatProperties:
    """Test properties of Delta T values."""

    @pytest.mark.unit
    def test_deltat_positive_modern(self):
        """Delta T should be positive for modern dates (after ~1902)."""
        jd = ephem.swe_julday(2000, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd)
        assert dt > 0, "Delta T should be positive for modern dates"

    @pytest.mark.unit
    def test_deltat_increases_over_time_modern(self):
        """Delta T should generally increase over modern era."""
        jd_1980 = ephem.swe_julday(1980, 1, 1, 12.0)
        jd_2000 = ephem.swe_julday(2000, 1, 1, 12.0)
        jd_2020 = ephem.swe_julday(2020, 1, 1, 12.0)

        dt_1980 = ephem.swe_deltat(jd_1980)
        dt_2000 = ephem.swe_deltat(jd_2000)
        dt_2020 = ephem.swe_deltat(jd_2020)

        assert dt_1980 < dt_2000 < dt_2020, "Delta T should increase over time"

    @pytest.mark.unit
    def test_deltat_continuous(self):
        """Delta T should change smoothly (no large jumps)."""
        jd_start = ephem.swe_julday(2000, 1, 1, 12.0)
        prev_dt = ephem.swe_deltat(jd_start)

        for i in range(1, 365):
            jd = jd_start + i
            dt = ephem.swe_deltat(jd)
            # Change should be small (less than 0.1 seconds per day)
            change = abs(dt - prev_dt) * 86400
            assert change < 0.1, f"Delta T changed by {change}s in one day"
            prev_dt = dt


class TestDeltatVsPyswisseph:
    """Compare Delta T with pyswisseph."""

    @pytest.mark.comparison
    def test_deltat_j2000_matches_swe(self):
        """Delta T at J2000 should match pyswisseph."""
        jd = 2451545.0
        dt_lib = ephem.swe_deltat(jd)
        dt_swe = swe.deltat(jd)
        # Allow 0.5 second tolerance (in days: 0.5/86400)
        assert dt_lib == pytest.approx(dt_swe, abs=0.5 / 86400)

    @pytest.mark.comparison
    @pytest.mark.parametrize("year", [1900, 1920, 1950, 1980, 2000, 2010, 2020, 2040])
    def test_deltat_various_years_match_swe(self, year):
        """Delta T should match pyswisseph for various years."""
        jd = ephem.swe_julday(year, 6, 15, 12.0)
        dt_lib = ephem.swe_deltat(jd)
        dt_swe = swe.deltat(jd)
        # Allow 1 second tolerance
        assert dt_lib == pytest.approx(dt_swe, abs=1.0 / 86400), (
            f"Year {year}: lib={dt_lib * 86400:.2f}s, swe={dt_swe * 86400:.2f}s"
        )

    @pytest.mark.comparison
    def test_deltat_100_dates_match_swe(self, random_dates_in_de421_range):
        """100 random dates should match pyswisseph."""
        dates = random_dates_in_de421_range(100)
        max_diff = 0
        for _, _, _, _, jd in dates:
            dt_lib = ephem.swe_deltat(jd)
            dt_swe = swe.deltat(jd)
            diff = abs(dt_lib - dt_swe) * 86400  # in seconds
            max_diff = max(max_diff, diff)
            # Allow 2 seconds tolerance
            assert diff < 2.0, f"Delta T diff = {diff}s at JD {jd}"
        print(f"Max Delta T difference: {max_diff:.3f} seconds")


class TestDeltatUTtoTT:
    """Test using Delta T for UT to TT conversion."""

    @pytest.mark.unit
    def test_tt_greater_than_ut(self):
        """TT should be greater than UT for modern dates."""
        jd_ut = ephem.swe_julday(2000, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd_ut)
        jd_tt = jd_ut + dt
        assert jd_tt > jd_ut

    @pytest.mark.unit
    def test_tt_ut_difference_seconds(self):
        """TT - UT should equal Delta T in seconds."""
        jd_ut = ephem.swe_julday(2000, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd_ut)
        # Delta T in seconds
        dt_seconds = dt * 86400
        # This should be approximately 63 seconds at J2000
        assert 60 < dt_seconds < 70


class TestDeltatEdgeCases:
    """Test edge cases for Delta T."""

    @pytest.mark.edge_case
    def test_deltat_de421_range_start(self):
        """Delta T at DE421 range start (1900)."""
        jd = ephem.swe_julday(1900, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd)
        # Should return a valid number
        assert isinstance(dt, float)
        assert -0.001 < dt < 0.001  # Within ±86 seconds

    @pytest.mark.edge_case
    def test_deltat_de421_range_end(self):
        """Delta T at DE421 range end (2050)."""
        jd = ephem.swe_julday(2050, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd)
        # Should return a valid number (probably ~80+ seconds)
        assert isinstance(dt, float)
        dt_seconds = dt * 86400
        assert 70 < dt_seconds < 100

    @pytest.mark.unit
    def test_deltat_returns_float(self):
        """Delta T should always return a float."""
        jd = 2451545.0
        dt = ephem.swe_deltat(jd)
        assert isinstance(dt, float)

    @pytest.mark.unit
    def test_deltat_reasonable_magnitude(self):
        """Delta T should have reasonable magnitude for any date in range."""
        for year in range(1900, 2051, 10):
            jd = ephem.swe_julday(year, 1, 1, 12.0)
            dt = ephem.swe_deltat(jd)
            dt_seconds = abs(dt * 86400)
            # Delta T should be less than 200 seconds for dates 1900-2050
            assert dt_seconds < 200, (
                f"Delta T at {year} = {dt_seconds}s seems too large"
            )
