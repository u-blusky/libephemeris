"""
Tests for SEFLG_EQUATORIAL conversion of lunar nodes and Lilith.

Regression tests ensuring that SEFLG_EQUATORIAL correctly converts
ecliptic coordinates (lon, lat) to equatorial coordinates (RA, Dec)
for lunar nodes, Lilith, and interpolated apogee/perigee.

Previously, these bodies had early returns in _calc_body() that bypassed
the equatorial conversion entirely, returning ecliptic coordinates
regardless of the SEFLG_EQUATORIAL flag.
"""

import math
import pytest

from libephemeris import swe_calc_ut
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
    SE_INTP_APOG,
    SE_INTP_PERG,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
)

# J2000.0 epoch
JD_J2000 = 2451545.0
# Summer solstice 2024 - Sun near max declination, good test date
JD_2024 = 2460479.5


class TestEquatorialConversionLunarNodes:
    """Tests that SEFLG_EQUATORIAL produces different results from ecliptic for nodes."""

    @pytest.mark.unit
    def test_mean_node_equatorial_differs_from_ecliptic(self):
        """Mean node ecliptic and equatorial coordinates must differ."""
        ecl, _ = swe_calc_ut(JD_J2000, SE_MEAN_NODE, SEFLG_SPEED)
        equ, _ = swe_calc_ut(JD_J2000, SE_MEAN_NODE, SEFLG_SPEED | SEFLG_EQUATORIAL)

        # Ecliptic lat for mean node is always 0 (on ecliptic by definition)
        # But declination should NOT be 0 (unless node is at 0° or 180° lon)
        assert ecl[1] == 0.0, "Mean node ecliptic latitude should be 0"
        assert equ[1] != 0.0, (
            f"Mean node declination should NOT be 0 "
            f"(was {equ[1]}, same as ecliptic lat - conversion not applied)"
        )

        # RA and lon should differ (obliquity rotation changes longitude)
        assert abs(ecl[0] - equ[0]) > 0.01, (
            f"Mean node RA ({equ[0]:.4f}) should differ from lon ({ecl[0]:.4f})"
        )

    @pytest.mark.unit
    def test_true_node_equatorial_differs_from_ecliptic(self):
        """True node ecliptic and equatorial coordinates must differ."""
        ecl, _ = swe_calc_ut(JD_J2000, SE_TRUE_NODE, SEFLG_SPEED)
        equ, _ = swe_calc_ut(JD_J2000, SE_TRUE_NODE, SEFLG_SPEED | SEFLG_EQUATORIAL)

        # For true node, ecliptic lat is near 0 but not exactly 0
        # Declination should be significantly different
        assert abs(ecl[1] - equ[1]) > 0.1 or abs(ecl[0] - equ[0]) > 0.01, (
            f"True node equatorial coords ({equ[0]:.4f}, {equ[1]:.4f}) "
            f"should differ from ecliptic ({ecl[0]:.4f}, {ecl[1]:.4f})"
        )

    @pytest.mark.unit
    def test_south_node_equatorial_differs_from_ecliptic(self):
        """South node (negative ID) should also convert to equatorial."""
        ecl, _ = swe_calc_ut(JD_J2000, -SE_MEAN_NODE, SEFLG_SPEED)
        equ, _ = swe_calc_ut(JD_J2000, -SE_MEAN_NODE, SEFLG_SPEED | SEFLG_EQUATORIAL)

        # South node is 180° from north - its declination should also be non-zero
        assert equ[1] != 0.0, (
            f"South node declination should NOT be 0 "
            f"(was {equ[1]}, conversion not applied)"
        )

    @pytest.mark.unit
    def test_mean_node_declination_range(self):
        """Mean node declination should be within a reasonable range."""
        equ, _ = swe_calc_ut(JD_J2000, SE_MEAN_NODE, SEFLG_SPEED | SEFLG_EQUATORIAL)
        dec = equ[1]

        # Declination must be within [-90, +90]
        assert -90.0 <= dec <= 90.0, f"Declination {dec} out of range"

        # Mean node declination should typically be within obliquity range (~23.4°)
        # since it lies on the ecliptic (lat=0)
        assert abs(dec) <= 25.0, (
            f"Mean node declination {dec:.4f} seems too large "
            f"(should be within obliquity range)"
        )


class TestEquatorialConversionLilith:
    """Tests that SEFLG_EQUATORIAL produces correct results for Lilith bodies."""

    @pytest.mark.unit
    def test_mean_lilith_equatorial_differs_from_ecliptic(self):
        """Mean Lilith ecliptic and equatorial coordinates must differ."""
        ecl, _ = swe_calc_ut(JD_J2000, SE_MEAN_APOG, SEFLG_SPEED)
        equ, _ = swe_calc_ut(JD_J2000, SE_MEAN_APOG, SEFLG_SPEED | SEFLG_EQUATORIAL)

        # The coordinates should be different due to obliquity rotation
        coords_differ = abs(ecl[0] - equ[0]) > 0.01 or abs(ecl[1] - equ[1]) > 0.01
        assert coords_differ, (
            f"Mean Lilith equatorial ({equ[0]:.4f}, {equ[1]:.4f}) "
            f"should differ from ecliptic ({ecl[0]:.4f}, {ecl[1]:.4f}). "
            f"SEFLG_EQUATORIAL not applied."
        )

    @pytest.mark.unit
    def test_oscu_lilith_equatorial_differs_from_ecliptic(self):
        """Osculating Lilith (True Lilith) equatorial coords must differ from ecliptic."""
        ecl, _ = swe_calc_ut(JD_J2000, SE_OSCU_APOG, SEFLG_SPEED)
        equ, _ = swe_calc_ut(JD_J2000, SE_OSCU_APOG, SEFLG_SPEED | SEFLG_EQUATORIAL)

        coords_differ = abs(ecl[0] - equ[0]) > 0.01 or abs(ecl[1] - equ[1]) > 0.01
        assert coords_differ, (
            f"Oscu Lilith equatorial ({equ[0]:.4f}, {equ[1]:.4f}) "
            f"should differ from ecliptic ({ecl[0]:.4f}, {ecl[1]:.4f}). "
            f"SEFLG_EQUATORIAL not applied."
        )

    @pytest.mark.unit
    def test_mean_lilith_declination_range(self):
        """Mean Lilith declination should be within physically reasonable range."""
        equ, _ = swe_calc_ut(JD_J2000, SE_MEAN_APOG, SEFLG_SPEED | SEFLG_EQUATORIAL)
        dec = equ[1]

        # Declination must be within [-90, +90]
        assert -90.0 <= dec <= 90.0, f"Declination {dec} out of range"

        # Mean Lilith has small ecliptic latitude (~5°), so declination
        # should be within obliquity + lat range
        assert abs(dec) <= 30.0, f"Mean Lilith declination {dec:.4f} seems too large"


class TestEquatorialConversionInterpolatedApogeePerigee:
    """Tests for SEFLG_EQUATORIAL with interpolated apogee/perigee."""

    @pytest.mark.unit
    def test_intp_apogee_equatorial_differs_from_ecliptic(self):
        """Interpolated apogee equatorial coords must differ from ecliptic."""
        ecl, _ = swe_calc_ut(JD_J2000, SE_INTP_APOG, SEFLG_SPEED)
        equ, _ = swe_calc_ut(JD_J2000, SE_INTP_APOG, SEFLG_SPEED | SEFLG_EQUATORIAL)

        coords_differ = abs(ecl[0] - equ[0]) > 0.01 or abs(ecl[1] - equ[1]) > 0.01
        assert coords_differ, (
            f"Intp apogee equatorial ({equ[0]:.4f}, {equ[1]:.4f}) "
            f"should differ from ecliptic ({ecl[0]:.4f}, {ecl[1]:.4f}). "
            f"SEFLG_EQUATORIAL not applied."
        )

    @pytest.mark.unit
    def test_intp_perigee_equatorial_differs_from_ecliptic(self):
        """Interpolated perigee equatorial coords must differ from ecliptic."""
        ecl, _ = swe_calc_ut(JD_J2000, SE_INTP_PERG, SEFLG_SPEED)
        equ, _ = swe_calc_ut(JD_J2000, SE_INTP_PERG, SEFLG_SPEED | SEFLG_EQUATORIAL)

        coords_differ = abs(ecl[0] - equ[0]) > 0.01 or abs(ecl[1] - equ[1]) > 0.01
        assert coords_differ, (
            f"Intp perigee equatorial ({equ[0]:.4f}, {equ[1]:.4f}) "
            f"should differ from ecliptic ({ecl[0]:.4f}, {ecl[1]:.4f}). "
            f"SEFLG_EQUATORIAL not applied."
        )


class TestEquatorialConsistencyWithPlanets:
    """Cross-check: equatorial conversion for lunar bodies should behave like planets."""

    @pytest.mark.unit
    def test_sun_equatorial_as_reference(self):
        """Verify Sun equatorial conversion works (reference for comparison)."""
        ecl, _ = swe_calc_ut(JD_J2000, SE_SUN, SEFLG_SPEED)
        equ, _ = swe_calc_ut(JD_J2000, SE_SUN, SEFLG_SPEED | SEFLG_EQUATORIAL)

        # Sun should have significantly different lat vs dec
        assert abs(ecl[1] - equ[1]) > 0.01 or abs(ecl[0] - equ[0]) > 0.01, (
            "Sun equatorial should differ from ecliptic (reference check)"
        )

    @pytest.mark.unit
    def test_all_lunar_bodies_equatorial_ra_in_range(self):
        """All lunar body RA values should be in [0, 360) range."""
        bodies = [
            (SE_MEAN_NODE, "Mean Node"),
            (SE_TRUE_NODE, "True Node"),
            (SE_MEAN_APOG, "Mean Lilith"),
            (SE_OSCU_APOG, "Oscu Lilith"),
            (SE_INTP_APOG, "Intp Apogee"),
            (SE_INTP_PERG, "Intp Perigee"),
        ]
        for body_id, name in bodies:
            equ, _ = swe_calc_ut(JD_J2000, body_id, SEFLG_SPEED | SEFLG_EQUATORIAL)
            ra = equ[0]
            dec = equ[1]
            assert 0.0 <= ra < 360.0, f"{name} RA {ra} not in [0, 360)"
            assert -90.0 <= dec <= 90.0, f"{name} Dec {dec} not in [-90, 90]"

    @pytest.mark.unit
    def test_equatorial_without_speed_flag(self):
        """SEFLG_EQUATORIAL should work without SEFLG_SPEED too."""
        ecl, _ = swe_calc_ut(JD_J2000, SE_MEAN_NODE, 0)
        equ, _ = swe_calc_ut(JD_J2000, SE_MEAN_NODE, SEFLG_EQUATORIAL)

        # Should still convert position even without speed
        assert ecl[1] == 0.0, "Mean node ecliptic latitude should be 0"
        assert equ[1] != 0.0, (
            "Mean node declination should NOT be 0 even without SEFLG_SPEED"
        )

    @pytest.mark.unit
    def test_equatorial_speed_transformed(self):
        """Equatorial velocity should also be transformed, not just position."""
        ecl, _ = swe_calc_ut(JD_J2000, SE_MEAN_NODE, SEFLG_SPEED)
        equ, _ = swe_calc_ut(JD_J2000, SE_MEAN_NODE, SEFLG_SPEED | SEFLG_EQUATORIAL)

        # Ecliptic: dlon != 0, dlat = 0 (mean node moves only in longitude)
        assert ecl[3] != 0.0, "Mean node ecliptic dlon should be non-zero"
        assert ecl[4] == 0.0, "Mean node ecliptic dlat should be 0"

        # Equatorial: both dRA and dDec should potentially be non-zero
        # since the obliquity rotation mixes lon/lat into RA/Dec
        assert equ[3] != 0.0, "Mean node dRA should be non-zero"
        # dDec might be very small but the velocity should be transformed
        # (not identical to ecliptic speed)
        assert abs(ecl[3] - equ[3]) > 0.001, (
            f"Mean node dRA ({equ[3]:.6f}) should differ from dlon ({ecl[3]:.6f})"
        )
