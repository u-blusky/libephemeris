"""
Tests for SEFLG_MOSEPH documentation and API visibility.

This module verifies that:
1. SEFLG_MOSEPH is exported in the public API
2. Documentation correctly describes Moshier ephemeris capabilities
3. The flag is usable and correctly routes calculations
"""

import pytest

import libephemeris as ephem


class TestSeflgMosephExport:
    """Test that SEFLG_MOSEPH is properly exported in the public API."""

    def test_seflg_moseph_exported(self):
        """SEFLG_MOSEPH should be exported from libephemeris."""
        assert hasattr(ephem, "SEFLG_MOSEPH")
        assert ephem.SEFLG_MOSEPH == 4

    def test_flg_moseph_exported(self):
        """FLG_MOSEPH alias should be exported from libephemeris."""
        assert hasattr(ephem, "FLG_MOSEPH")
        assert ephem.FLG_MOSEPH == 4

    def test_seflg_moseph_in_all(self):
        """SEFLG_MOSEPH should be in __all__ list."""
        assert "SEFLG_MOSEPH" in ephem.__all__

    def test_flg_moseph_in_all(self):
        """FLG_MOSEPH should be in __all__ list."""
        assert "FLG_MOSEPH" in ephem.__all__

    def test_all_ephemeris_flags_exported(self):
        """All three ephemeris selection flags should be exported."""
        # SEFLG_* prefix versions
        assert hasattr(ephem, "SEFLG_JPLEPH")
        assert hasattr(ephem, "SEFLG_SWIEPH")
        assert hasattr(ephem, "SEFLG_MOSEPH")
        # FLG_* prefix versions
        assert hasattr(ephem, "FLG_JPLEPH")
        assert hasattr(ephem, "FLG_SWIEPH")
        assert hasattr(ephem, "FLG_MOSEPH")

    def test_ephemeris_flag_values(self):
        """Ephemeris flags should have correct values (1, 2, 4)."""
        assert ephem.SEFLG_JPLEPH == 1
        assert ephem.SEFLG_SWIEPH == 2
        assert ephem.SEFLG_MOSEPH == 4


class TestSeflgMosephFunctionality:
    """Test that SEFLG_MOSEPH works correctly in calculations."""

    @pytest.mark.unit
    def test_moshier_calc_sun(self):
        """SEFLG_MOSEPH should work for Sun calculations."""
        jd = 2451545.0  # J2000
        pos, retflag = ephem.swe_calc_ut(jd, ephem.SE_SUN, ephem.SEFLG_MOSEPH)
        # Sun should be near 280 degrees longitude at J2000
        assert 279.0 < pos[0] < 282.0

    @pytest.mark.unit
    def test_moshier_calc_moon(self):
        """SEFLG_MOSEPH should work for Moon calculations."""
        jd = 2451545.0  # J2000
        pos, retflag = ephem.swe_calc_ut(jd, ephem.SE_MOON, ephem.SEFLG_MOSEPH)
        # Result should be valid ecliptic coordinates
        assert 0.0 <= pos[0] < 360.0
        assert -90.0 <= pos[1] <= 90.0

    @pytest.mark.unit
    def test_moshier_extended_range(self):
        """SEFLG_MOSEPH should work for dates outside DE440 range."""
        # Ancient date: 1000 BCE (outside DE440 range of 1550-2650 CE)
        jd = 1355808.0  # ~1000 BCE
        pos, retflag = ephem.swe_calc_ut(jd, ephem.SE_MARS, ephem.SEFLG_MOSEPH)
        # Result should be valid ecliptic coordinates
        assert 0.0 <= pos[0] < 360.0

    @pytest.mark.unit
    def test_moshier_deltat_no_warning(self):
        """swe_deltat_ex with SEFLG_MOSEPH should not produce a warning."""
        jd = 2451545.0  # J2000
        dt, serr = ephem.swe_deltat_ex(jd, ephem.SEFLG_MOSEPH)
        # No warning - Moshier uses the same Skyfield Delta T model
        assert serr == ""
        # Should return valid Delta T
        assert isinstance(dt, float)


class TestSeflgMosephDocumentation:
    """Test that SEFLG_MOSEPH documentation is accurate."""

    def test_swe_calc_ut_docstring_mentions_moseph(self):
        """swe_calc_ut docstring should document SEFLG_MOSEPH."""
        docstring = ephem.swe_calc_ut.__doc__
        assert "SEFLG_MOSEPH" in docstring
        assert "Moshier" in docstring

    def test_swe_calc_ut_docstring_describes_range(self):
        """swe_calc_ut docstring should mention Moshier range."""
        docstring = ephem.swe_calc_ut.__doc__
        # Should mention the extended date range
        assert "-3000" in docstring or "3000" in docstring

    def test_swe_calc_ut_docstring_describes_limitations(self):
        """swe_calc_ut docstring should mention Moshier limitations."""
        docstring = ephem.swe_calc_ut.__doc__
        # Should mention asteroids are not supported
        assert "asteroid" in docstring.lower() or "SE_AST_OFFSET" in docstring


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
