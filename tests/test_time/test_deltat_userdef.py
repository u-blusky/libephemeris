"""
Tests for user-defined Delta T functionality (set_delta_t_userdef).

This feature allows overriding the computed Delta T value with a fixed
user-specified value, useful for testing and for ancient/future dates
where Delta T is uncertain.
"""

import pytest
import libephemeris as ephem


class TestSetDeltaTUserdef:
    """Test the set_delta_t_userdef function."""

    @pytest.fixture(autouse=True)
    def reset_delta_t(self):
        """Reset user-defined Delta T before and after each test."""
        ephem.set_delta_t_userdef(None)
        yield
        ephem.set_delta_t_userdef(None)

    @pytest.mark.unit
    def test_set_delta_t_userdef_sets_value(self):
        """Setting a user-defined Delta T should store the value."""
        dt_value = 0.001  # ~86 seconds in days
        ephem.set_delta_t_userdef(dt_value)
        assert ephem.get_delta_t_userdef() == dt_value

    @pytest.mark.unit
    def test_set_delta_t_userdef_none_clears_value(self):
        """Setting None should clear the user-defined value."""
        ephem.set_delta_t_userdef(0.001)
        assert ephem.get_delta_t_userdef() is not None
        ephem.set_delta_t_userdef(None)
        assert ephem.get_delta_t_userdef() is None

    @pytest.mark.unit
    def test_get_delta_t_userdef_default_is_none(self):
        """get_delta_t_userdef should return None by default."""
        assert ephem.get_delta_t_userdef() is None

    @pytest.mark.unit
    def test_set_delta_t_userdef_zero_is_valid(self):
        """Zero should be a valid user-defined Delta T value."""
        ephem.set_delta_t_userdef(0.0)
        assert ephem.get_delta_t_userdef() == 0.0

    @pytest.mark.unit
    def test_set_delta_t_userdef_negative_is_valid(self):
        """Negative values should be valid (for dates before ~1902)."""
        dt_value = -0.0001  # ~-8.6 seconds
        ephem.set_delta_t_userdef(dt_value)
        assert ephem.get_delta_t_userdef() == dt_value


class TestDeltatWithUserdef:
    """Test swe_deltat when user-defined value is set."""

    @pytest.fixture(autouse=True)
    def reset_delta_t(self):
        """Reset user-defined Delta T before and after each test."""
        ephem.set_delta_t_userdef(None)
        yield
        ephem.set_delta_t_userdef(None)

    @pytest.mark.unit
    def test_swe_deltat_returns_userdef_value(self):
        """swe_deltat should return user-defined value when set."""
        dt_value = 0.00075  # ~65 seconds in days
        ephem.set_delta_t_userdef(dt_value)

        jd = 2451545.0  # J2000
        result = ephem.swe_deltat(jd)
        assert result == dt_value

    @pytest.mark.unit
    def test_swe_deltat_ignores_jd_when_userdef_set(self):
        """swe_deltat should return same value regardless of JD when userdef is set."""
        dt_value = 0.00075
        ephem.set_delta_t_userdef(dt_value)

        jd1 = 2451545.0  # J2000
        jd2 = 2460000.0  # A future date
        jd3 = 2400000.0  # A past date

        assert ephem.swe_deltat(jd1) == dt_value
        assert ephem.swe_deltat(jd2) == dt_value
        assert ephem.swe_deltat(jd3) == dt_value

    @pytest.mark.unit
    def test_swe_deltat_uses_computed_when_cleared(self):
        """swe_deltat should use computed value after clearing userdef."""
        dt_userdef = 0.00075
        ephem.set_delta_t_userdef(dt_userdef)

        jd = 2451545.0  # J2000
        assert ephem.swe_deltat(jd) == dt_userdef

        # Clear the user-defined value
        ephem.set_delta_t_userdef(None)

        # Now should get computed value (not equal to our arbitrary value)
        computed_dt = ephem.swe_deltat(jd)
        assert computed_dt != dt_userdef
        # Should be approximately 63-64 seconds at J2000
        dt_seconds = computed_dt * 86400
        assert 60 < dt_seconds < 70

    @pytest.mark.unit
    def test_swe_deltat_userdef_zero_returns_zero(self):
        """swe_deltat should return exactly 0.0 when userdef is set to 0."""
        ephem.set_delta_t_userdef(0.0)
        jd = 2451545.0  # J2000
        assert ephem.swe_deltat(jd) == 0.0


class TestDeltatExWithUserdef:
    """Test swe_deltat_ex when user-defined value is set."""

    @pytest.fixture(autouse=True)
    def reset_delta_t(self):
        """Reset user-defined Delta T before and after each test."""
        ephem.set_delta_t_userdef(None)
        yield
        ephem.set_delta_t_userdef(None)

    @pytest.mark.unit
    def test_swe_deltat_ex_returns_userdef_value(self):
        """swe_deltat_ex should return user-defined value when set."""
        dt_value = 0.00075  # ~65 seconds in days
        ephem.set_delta_t_userdef(dt_value)

        jd = 2451545.0  # J2000
        result, serr = ephem.swe_deltat_ex(jd, ephem.SEFLG_SWIEPH)
        assert result == dt_value
        assert serr == ""

    @pytest.mark.unit
    def test_swe_deltat_ex_userdef_ignores_flag(self):
        """swe_deltat_ex should return userdef regardless of ephemeris flag."""
        dt_value = 0.00075
        ephem.set_delta_t_userdef(dt_value)

        jd = 2451545.0

        result_swieph, _ = ephem.swe_deltat_ex(jd, ephem.SEFLG_SWIEPH)
        result_jpleph, _ = ephem.swe_deltat_ex(jd, ephem.SEFLG_JPLEPH)
        result_moseph, serr_moseph = ephem.swe_deltat_ex(jd, ephem.SEFLG_MOSEPH)

        # All should return the user-defined value
        assert result_swieph == dt_value
        assert result_jpleph == dt_value
        assert result_moseph == dt_value
        # When userdef is set, no warning should be generated
        assert serr_moseph == ""

    @pytest.mark.unit
    def test_swe_deltat_ex_matches_swe_deltat_with_userdef(self):
        """swe_deltat and swe_deltat_ex should match when userdef is set."""
        dt_value = 0.00075
        ephem.set_delta_t_userdef(dt_value)

        jd = 2451545.0
        dt = ephem.swe_deltat(jd)
        dt_ex, _ = ephem.swe_deltat_ex(jd, ephem.SEFLG_SWIEPH)

        assert dt == dt_ex == dt_value


class TestDeltaTUserdefAliases:
    """Test that both swe_ and non-prefixed aliases work."""

    @pytest.fixture(autouse=True)
    def reset_delta_t(self):
        """Reset user-defined Delta T before and after each test."""
        ephem.set_delta_t_userdef(None)
        yield
        ephem.set_delta_t_userdef(None)

    @pytest.mark.unit
    def test_swe_set_delta_t_userdef_exists(self):
        """swe_set_delta_t_userdef alias should exist."""
        assert hasattr(ephem, "swe_set_delta_t_userdef")

    @pytest.mark.unit
    def test_swe_get_delta_t_userdef_exists(self):
        """swe_get_delta_t_userdef alias should exist."""
        assert hasattr(ephem, "swe_get_delta_t_userdef")

    @pytest.mark.unit
    def test_set_delta_t_userdef_exists(self):
        """set_delta_t_userdef alias should exist."""
        assert hasattr(ephem, "set_delta_t_userdef")

    @pytest.mark.unit
    def test_get_delta_t_userdef_exists(self):
        """get_delta_t_userdef alias should exist."""
        assert hasattr(ephem, "get_delta_t_userdef")

    @pytest.mark.unit
    def test_aliases_are_same_function(self):
        """swe_ and non-prefixed aliases should be the same function."""
        assert ephem.swe_set_delta_t_userdef is ephem.set_delta_t_userdef
        assert ephem.swe_get_delta_t_userdef is ephem.get_delta_t_userdef


class TestDeltaTUserdefUseCases:
    """Test practical use cases for user-defined Delta T."""

    @pytest.fixture(autouse=True)
    def reset_delta_t(self):
        """Reset user-defined Delta T before and after each test."""
        ephem.set_delta_t_userdef(None)
        yield
        ephem.set_delta_t_userdef(None)

    @pytest.mark.unit
    def test_ancient_date_with_large_delta_t(self):
        """Test using a large Delta T for ancient dates."""
        # For 2000 BCE, Delta T could be several hours
        # Set 3 hours = 0.125 days
        large_dt = 3.0 / 24.0  # 3 hours in days
        ephem.set_delta_t_userdef(large_dt)

        jd_ancient = ephem.swe_julday(-2000, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd_ancient)

        assert dt == large_dt
        assert dt * 24 == 3.0  # Verify it's 3 hours

    @pytest.mark.unit
    def test_testing_reproducibility(self):
        """Test using fixed Delta T for reproducible calculations."""
        # Set a fixed Delta T for testing
        fixed_dt = 64.0 / 86400.0  # exactly 64 seconds

        ephem.set_delta_t_userdef(fixed_dt)

        # Multiple calls should always return the same value
        results = [ephem.swe_deltat(2451545.0 + i) for i in range(10)]
        assert all(r == fixed_dt for r in results)

    @pytest.mark.unit
    def test_toggle_between_computed_and_fixed(self):
        """Test toggling between computed and user-defined Delta T."""
        jd = 2451545.0  # J2000

        # Get computed value
        computed = ephem.swe_deltat(jd)

        # Set user-defined value
        fixed_dt = 0.001
        ephem.set_delta_t_userdef(fixed_dt)
        assert ephem.swe_deltat(jd) == fixed_dt

        # Toggle back to computed
        ephem.set_delta_t_userdef(None)
        assert ephem.swe_deltat(jd) == computed

        # Toggle to user-defined again
        ephem.set_delta_t_userdef(fixed_dt)
        assert ephem.swe_deltat(jd) == fixed_dt

    @pytest.mark.unit
    def test_future_date_uncertainty(self):
        """Test using estimated Delta T for far future dates."""
        # For year 3000, Delta T is unknown - we might estimate ~500 seconds
        future_dt = 500.0 / 86400.0  # 500 seconds in days
        ephem.set_delta_t_userdef(future_dt)

        jd_future = ephem.swe_julday(3000, 1, 1, 12.0)
        dt = ephem.swe_deltat(jd_future)

        assert dt == future_dt
        assert dt * 86400 == pytest.approx(500.0)
