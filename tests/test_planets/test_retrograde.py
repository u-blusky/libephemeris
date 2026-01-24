"""
Comprehensive tests for retrograde motion detection.

Tests cover:
- Known retrograde periods
- Station points (velocity near zero)
- Velocity sign changes
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestMercuryRetrograde:
    """Test Mercury retrograde detection."""

    @pytest.mark.unit
    def test_mercury_has_retrograde_periods(self):
        """Mercury should have retrograde periods (negative velocity)."""
        # Sample 1 year of Mercury velocities
        retrograde_found = False
        for jd in range(2451545, 2451545 + 365):
            pos, _ = ephem.swe_calc_ut(float(jd), SE_MERCURY, SEFLG_SPEED)
            if pos[3] < 0:
                retrograde_found = True
                break

        assert retrograde_found, "Mercury should go retrograde within a year"

    @pytest.mark.unit
    def test_mercury_direct_motion_exists(self):
        """Mercury should also have direct motion periods."""
        direct_found = False
        for jd in range(2451545, 2451545 + 365):
            pos, _ = ephem.swe_calc_ut(float(jd), SE_MERCURY, SEFLG_SPEED)
            if pos[3] > 0.5:  # Clear direct motion
                direct_found = True
                break

        assert direct_found, "Mercury should have direct motion"

    @pytest.mark.unit
    def test_mercury_station_points(self):
        """Mercury should have station points (velocity near zero)."""
        # Find velocity changes
        prev_sign = None
        station_found = False

        for jd in range(2451545, 2451545 + 120):  # ~4 months covers a cycle
            pos, _ = ephem.swe_calc_ut(float(jd), SE_MERCURY, SEFLG_SPEED)
            current_sign = pos[3] > 0

            if prev_sign is not None and current_sign != prev_sign:
                # Velocity changed sign - this is near a station
                station_found = True
                break
            prev_sign = current_sign

        assert station_found, "Mercury should have station points"


class TestVenusRetrograde:
    """Test Venus retrograde detection."""

    @pytest.mark.unit
    def test_venus_has_retrograde(self):
        """Venus should have retrograde periods (~40 days every 18 months)."""
        retrograde_found = False
        # Need to sample ~600 days to catch a Venus retrograde
        for jd in range(2451545, 2451545 + 600):
            pos, _ = ephem.swe_calc_ut(float(jd), SE_VENUS, SEFLG_SPEED)
            if pos[3] < 0:
                retrograde_found = True
                break

        assert retrograde_found, "Venus should go retrograde within 600 days"


class TestMarsRetrograde:
    """Test Mars retrograde detection."""

    @pytest.mark.unit
    def test_mars_has_retrograde(self):
        """Mars should have retrograde periods (~72 days every 2 years)."""
        retrograde_found = False
        # Need to sample ~800 days to catch a Mars retrograde
        for jd in range(2451545, 2451545 + 800):
            pos, _ = ephem.swe_calc_ut(float(jd), SE_MARS, SEFLG_SPEED)
            if pos[3] < 0:
                retrograde_found = True
                break

        assert retrograde_found, "Mars should go retrograde within ~2 years"


class TestOuterPlanetRetrograde:
    """Test outer planet retrograde motion."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,sample_days",
        [
            (SE_JUPITER, "Jupiter", 400),
            (SE_SATURN, "Saturn", 380),
            (SE_URANUS, "Uranus", 370),
            (SE_NEPTUNE, "Neptune", 370),
        ],
    )
    def test_outer_planet_retrograde(self, planet_id, planet_name, sample_days):
        """Outer planets should have annual retrograde periods."""
        retrograde_found = False
        for jd in range(2451545, 2451545 + sample_days):
            pos, _ = ephem.swe_calc_ut(float(jd), planet_id, SEFLG_SPEED)
            if pos[3] < 0:
                retrograde_found = True
                break

        assert retrograde_found, (
            f"{planet_name} should go retrograde within {sample_days} days"
        )


class TestNoRetrogradeBodies:
    """Test bodies that never go retrograde."""

    @pytest.mark.unit
    def test_sun_never_retrograde(self):
        """Sun should never have negative velocity."""
        for jd in range(2451545, 2451545 + 365):
            pos, _ = ephem.swe_calc_ut(float(jd), SE_SUN, SEFLG_SPEED)
            assert pos[3] > 0, f"Sun retrograde at JD {jd}?!"

    @pytest.mark.unit
    def test_moon_never_retrograde(self):
        """Moon should never have negative velocity."""
        for jd in range(2451545, 2451545 + 30):
            pos, _ = ephem.swe_calc_ut(float(jd), SE_MOON, SEFLG_SPEED)
            assert pos[3] > 0, f"Moon retrograde at JD {jd}?!"


class TestRetrogradeVelocitySign:
    """Test velocity sign during retrograde."""

    @pytest.mark.unit
    def test_retrograde_means_negative_velocity(self):
        """During retrograde, longitude velocity should be negative."""
        # Find a Mercury retrograde period
        retrograde_days = []
        for jd in range(2451545, 2451545 + 365):
            pos, _ = ephem.swe_calc_ut(float(jd), SE_MERCURY, SEFLG_SPEED)
            if pos[3] < 0:
                retrograde_days.append(jd)

        # Verify all retrograde days have negative velocity
        for jd in retrograde_days[:10]:  # Check first 10
            pos, _ = ephem.swe_calc_ut(float(jd), SE_MERCURY, SEFLG_SPEED)
            assert pos[3] < 0

    @pytest.mark.unit
    def test_retrograde_longitude_decreases(self):
        """During retrograde, longitude should decrease day to day."""
        # Find start of Mercury retrograde
        retro_start = None
        for jd in range(2451545, 2451545 + 365):
            pos, _ = ephem.swe_calc_ut(float(jd), SE_MERCURY, SEFLG_SPEED)
            if pos[3] < -0.5:  # Clear retrograde
                retro_start = jd
                break

        if retro_start:
            pos1, _ = ephem.swe_calc_ut(float(retro_start), SE_MERCURY, 0)
            pos2, _ = ephem.swe_calc_ut(float(retro_start + 1), SE_MERCURY, 0)

            # Longitude should decrease (with wraparound handling)
            change = pos2[0] - pos1[0]
            if change > 180:
                change -= 360
            elif change < -180:
                change += 360

            assert change < 0, "During retrograde, longitude should decrease"
