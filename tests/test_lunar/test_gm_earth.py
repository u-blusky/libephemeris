"""
Tests for Earth's gravitational parameter (GM) in lunar calculations.

The GM value used in orbital calculations must match the IAU 2015 Resolution B3
TDB-compatible value used by the DE440 ephemeris for consistency.

Reference:
    IAU 2015 Resolution B3: Recommended Nominal Values for Selected Solar and
    Planetary Quantities
    GM_Earth (TDB) = 398600.435436 km^3/s^2
"""

import re
import math
from pathlib import Path


# The IAU 2015 Resolution B3 TDB-compatible GM value for Earth in km^3/s^2
# This is the value used by DE440/DE441 ephemeris files
GM_EARTH_IAU2015_KM3_S2 = 398600.435436

# Conversion factors
AU_KM = 149597870.7  # km per AU (IAU 2012 definition)
SECONDS_PER_DAY = 86400.0


class TestGMEarthConstant:
    """Tests for the GM_Earth constant value used in lunar.py."""

    def test_gm_value_in_source_code(self):
        """Verify that lunar.py uses the IAU 2015 GM_Earth value."""
        lunar_py = Path(__file__).parent.parent.parent / "libephemeris" / "lunar.py"
        content = lunar_py.read_text()

        # Find the mu calculation line
        # Pattern: mu = 398600.xxx... / (149597870.7**3) * (86400**2)
        pattern = r"mu\s*=\s*(\d+\.?\d*)\s*/\s*\(149597870\.7\*\*3\)"
        match = re.search(pattern, content)

        assert match is not None, (
            "Could not find GM_Earth constant (mu = ...) in lunar.py"
        )

        gm_value = float(match.group(1))

        # Verify it matches IAU 2015 value within a small tolerance
        assert abs(gm_value - GM_EARTH_IAU2015_KM3_S2) < 0.001, (
            f"GM_Earth value {gm_value} does not match IAU 2015 value "
            f"{GM_EARTH_IAU2015_KM3_S2} km^3/s^2"
        )

    def test_gm_value_documented_in_docstring(self):
        """Verify that the GM value is documented in the docstring."""
        lunar_py = Path(__file__).parent.parent.parent / "libephemeris" / "lunar.py"
        content = lunar_py.read_text()

        # Check that the IAU 2015 value is mentioned in docstring
        assert "398600.435436" in content, (
            "GM_Earth value 398600.435436 should be documented in lunar.py"
        )

        # Check that IAU 2015 is referenced
        assert "IAU 2015" in content, (
            "IAU 2015 Resolution B3 should be referenced in lunar.py"
        )

    def test_gm_conversion_to_au_day(self):
        """Test the conversion of GM from km^3/s^2 to AU^3/day^2."""
        # Convert from km^3/s^2 to AU^3/day^2
        # Formula: GM (AU^3/day^2) = GM (km^3/s^2) / (AU_km^3) * (sec_per_day^2)
        gm_au3_day2 = GM_EARTH_IAU2015_KM3_S2 / (AU_KM**3) * (SECONDS_PER_DAY**2)

        # Expected value approximately 8.887e-10 AU^3/day^2
        # This is used in orbital mechanics for Moon's orbit
        assert gm_au3_day2 > 0, "GM in AU^3/day^2 must be positive"
        assert 8.8e-10 < gm_au3_day2 < 9.0e-10, (
            f"GM conversion result {gm_au3_day2} is outside expected range"
        )

    def test_gm_relative_precision(self):
        """
        Test that the IAU 2015 value provides sufficient precision.

        The relative difference between IAU 2015 and older values (e.g., WGS-84)
        is small but significant for precise orbital calculations.
        """
        # Old WGS-84/EGM96 value (circa 1996)
        gm_old = 398600.4418

        # Relative difference
        rel_diff = abs(GM_EARTH_IAU2015_KM3_S2 - gm_old) / GM_EARTH_IAU2015_KM3_S2

        # The difference should be on the order of 1e-8 (about 16 ppb)
        assert rel_diff < 1e-6, "Relative difference is unexpectedly large"
        assert rel_diff > 1e-10, "Values should not be exactly equal"

        # The new value is slightly smaller (better satellite tracking data)
        assert GM_EARTH_IAU2015_KM3_S2 < gm_old, (
            "IAU 2015 value should be smaller than old WGS-84 value"
        )


class TestGMImpactOnOrbitalCalculations:
    """Tests for the impact of GM value on orbital calculations."""

    def test_orbital_velocity_formula(self):
        """
        Test that circular orbital velocity formula works with the GM value.

        v = sqrt(GM/r) where v is velocity, GM is gravitational parameter, r is radius
        """
        # Mean Moon distance in km
        moon_distance_km = 384400.0

        # Calculate circular orbital velocity
        v_km_s = math.sqrt(GM_EARTH_IAU2015_KM3_S2 / moon_distance_km)

        # Moon's actual mean orbital velocity is about 1.022 km/s
        assert 1.0 < v_km_s < 1.1, (
            f"Calculated velocity {v_km_s} km/s is outside expected range"
        )

    def test_orbital_period_formula(self):
        """
        Test that orbital period formula works with the GM value.

        T = 2*pi*sqrt(a^3/GM) where T is period, a is semi-major axis
        """
        # Moon's mean semi-major axis in km
        moon_sma_km = 384748.0

        # Calculate orbital period in seconds
        T_seconds = 2 * math.pi * math.sqrt(moon_sma_km**3 / GM_EARTH_IAU2015_KM3_S2)

        # Convert to days
        T_days = T_seconds / SECONDS_PER_DAY

        # Moon's sidereal orbital period is about 27.32 days
        assert 27.0 < T_days < 28.0, (
            f"Calculated period {T_days} days is outside expected range"
        )
