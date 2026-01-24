"""
Comprehensive tests for all 43 ayanamsha modes.

Tests sidereal calculations and ayanamsha values.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestAyanamshaBasicValues:
    """Test basic ayanamsha value calculations."""

    @pytest.mark.unit
    def test_lahiri_at_j2000(self):
        """Lahiri ayanamsha at J2000 should be ~23.9 degrees."""
        jd = 2451545.0
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)
        ayan = ephem.swe_get_ayanamsa_ut(jd)

        assert 23.0 < ayan < 25.0, f"Lahiri ayanamsha {ayan} unexpected"

    @pytest.mark.unit
    def test_fagan_bradley_at_j2000(self):
        """Fagan-Bradley ayanamsha at J2000 should be ~24.7 degrees."""
        jd = 2451545.0
        ephem.swe_set_sid_mode(SE_SIDM_FAGAN_BRADLEY)
        ayan = ephem.swe_get_ayanamsa_ut(jd)

        assert 24.0 < ayan < 26.0, f"Fagan-Bradley ayanamsha {ayan} unexpected"

    @pytest.mark.unit
    def test_ayanamsha_positive(self):
        """Ayanamsha should be positive for modern dates."""
        jd = 2451545.0
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)
        ayan = ephem.swe_get_ayanamsa_ut(jd)

        assert ayan > 0


class TestAyanamshaVsPyswisseph:
    """Compare ayanamsha values with pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "sid_mode,name",
        [
            (SE_SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
            (SE_SIDM_LAHIRI, "Lahiri"),
            (SE_SIDM_RAMAN, "Raman"),
            (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
            (SE_SIDM_J2000, "J2000"),
            (SE_SIDM_B1950, "B1950"),
        ],
    )
    def test_standard_ayanamsha_matches_swe(self, sid_mode, name):
        """Standard ayanamshas should match pyswisseph closely."""
        jd = 2451545.0
        tolerance = 0.06  # degrees

        ephem.swe_set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        ayan_lib = ephem.swe_get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        diff = abs(ayan_lib - ayan_swe)
        assert diff < tolerance, f"{name} diff {diff} >= {tolerance}"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "sid_mode,name",
        [
            (SE_SIDM_TRUE_CITRA, "True Citra"),
            (SE_SIDM_TRUE_REVATI, "True Revati"),
            (SE_SIDM_GALCENT_0SAG, "Galcent 0 Sag"),
        ],
    )
    def test_star_based_ayanamsha_relaxed(self, sid_mode, name):
        """Star-based ayanamshas may need relaxed tolerance."""
        jd = 2451545.0
        tolerance = 1.0  # degrees (relaxed)

        ephem.swe_set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        ayan_lib = ephem.swe_get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        diff = abs(ayan_lib - ayan_swe)
        assert diff < tolerance, f"{name} diff {diff} >= {tolerance}"


class TestSiderealPositions:
    """Test sidereal planet positions."""

    @pytest.mark.unit
    def test_sidereal_sun_differs_from_tropical(self):
        """Sidereal Sun should differ from tropical by ayanamsha."""
        jd = 2451545.0

        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        pos_trop, _ = ephem.swe_calc_ut(jd, SE_SUN, 0)
        pos_sid, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)

        diff = pos_trop[0] - pos_sid[0]
        if diff < 0:
            diff += 360

        # Difference should be close to ayanamsha
        ayan = ephem.swe_get_ayanamsa_ut(jd)
        assert abs(diff - ayan) < 0.01, f"Diff {diff} != ayanamsha {ayan}"

    @pytest.mark.comparison
    def test_sidereal_sun_matches_pyswisseph(self):
        """Sidereal Sun should match pyswisseph."""
        jd = 2451545.0

        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)
        swe.set_sid_mode(SE_SIDM_LAHIRI)

        pos_lib, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)
        pos_swe, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)

        diff = abs(pos_lib[0] - pos_swe[0])
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Sidereal Sun diff {diff}"


class TestAyanamshaProgression:
    """Test that ayanamsha progresses correctly over time."""

    @pytest.mark.unit
    def test_ayanamsha_increases_over_time(self):
        """Ayanamsha should increase with time (precession)."""
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        jd_2000 = 2451545.0
        jd_2020 = ephem.swe_julday(2020, 1, 1, 12.0)

        ayan_2000 = ephem.swe_get_ayanamsa_ut(jd_2000)
        ayan_2020 = ephem.swe_get_ayanamsa_ut(jd_2020)

        assert ayan_2020 > ayan_2000, "Ayanamsha should increase over time"

    @pytest.mark.unit
    def test_ayanamsha_precession_rate(self):
        """Ayanamsha should increase ~50 arcsec/year."""
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        jd_2000 = 2451545.0
        jd_2020 = ephem.swe_julday(2020, 1, 1, 12.0)

        ayan_2000 = ephem.swe_get_ayanamsa_ut(jd_2000)
        ayan_2020 = ephem.swe_get_ayanamsa_ut(jd_2020)

        # ~20 years × 50"/year = 1000" = ~0.28°
        diff = ayan_2020 - ayan_2000
        assert 0.2 < diff < 0.4, f"20-year ayanamsha change {diff} unexpected"


class TestAyanamshaName:
    """Test swe_get_ayanamsa_name function."""

    @pytest.mark.unit
    def test_lahiri_name(self):
        """Should return 'Lahiri' for SE_SIDM_LAHIRI."""
        name = ephem.swe_get_ayanamsa_name(SE_SIDM_LAHIRI)
        assert "Lahiri" in name or "lahiri" in name.lower()

    @pytest.mark.unit
    def test_fagan_bradley_name(self):
        """Should return 'Fagan' or 'Bradley' for SE_SIDM_FAGAN_BRADLEY."""
        name = ephem.swe_get_ayanamsa_name(SE_SIDM_FAGAN_BRADLEY)
        assert "Fagan" in name or "Bradley" in name


class TestAllAyanamshas:
    """Test all 43 ayanamsha modes work."""

    @pytest.mark.unit
    @pytest.mark.parametrize("sid_mode", range(43))
    def test_ayanamsha_mode_valid(self, sid_mode):
        """Each ayanamsha mode should return valid value."""
        jd = 2451545.0

        try:
            ephem.swe_set_sid_mode(sid_mode)
            ayan = ephem.swe_get_ayanamsa_ut(jd)

            # Should be reasonable (some modes might give very different values)
            assert isinstance(ayan, float)
            assert -360 < ayan < 360  # Very permissive
        except Exception as e:
            pytest.fail(f"Mode {sid_mode} failed: {e}")
