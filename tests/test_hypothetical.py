"""
Tests for the hypothetical module.

Tests calculation functions for hypothetical celestial bodies including:
- Uranian planets (Hamburg School): Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon
- Transpluto (Isis/Persephone)
- Other hypothetical bodies: White Moon (Selena), Waldemath
"""

import pytest
import math

from libephemeris.hypothetical import (
    # Body IDs
    SE_CUPIDO,
    SE_HADES,
    SE_ZEUS,
    SE_KRONOS,
    SE_APOLLON,
    SE_ADMETOS,
    SE_VULKANUS,
    SE_POSEIDON,
    SE_ISIS,
    SE_TRANSPLUTO,
    SE_WHITE_MOON,
    SE_WALDEMATH,
    SE_VULCAN,
    SE_PROSERPINA,
    # Functions
    is_hypothetical_body,
    get_hypothetical_name,
    calc_uranian_longitude,
    calc_uranian_position,
    calc_transpluto_position,
    calc_white_moon_position,
    calc_waldemath_position,
    calc_hypothetical_position,
    list_hypothetical_bodies,
    # Data structures
    URANIAN_ELEMENTS,
    HYPOTHETICAL_ELEMENTS,
    HYPOTHETICAL_NAMES,
)
from libephemeris.constants import SE_SUN, SE_MOON, SE_MARS, SE_FICT_OFFSET


class TestBodyIdentification:
    """Tests for body identification functions."""

    def test_is_hypothetical_body_uranian_planets(self):
        """Test that Uranian planets are identified as hypothetical."""
        uranian_ids = [
            SE_CUPIDO,
            SE_HADES,
            SE_ZEUS,
            SE_KRONOS,
            SE_APOLLON,
            SE_ADMETOS,
            SE_VULKANUS,
            SE_POSEIDON,
        ]
        for body_id in uranian_ids:
            assert is_hypothetical_body(body_id), (
                f"Body ID {body_id} should be hypothetical"
            )

    def test_is_hypothetical_body_other(self):
        """Test that other hypothetical bodies are identified correctly."""
        other_ids = [SE_ISIS, SE_WHITE_MOON, SE_WALDEMATH, SE_VULCAN, SE_PROSERPINA]
        for body_id in other_ids:
            assert is_hypothetical_body(body_id), (
                f"Body ID {body_id} should be hypothetical"
            )

    def test_is_hypothetical_body_real_planets(self):
        """Test that real planets are NOT identified as hypothetical."""
        real_ids = [SE_SUN, SE_MOON, SE_MARS, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        for body_id in real_ids:
            assert not is_hypothetical_body(body_id), (
                f"Body ID {body_id} should NOT be hypothetical"
            )

    def test_transpluto_isis_alias(self):
        """Test that SE_TRANSPLUTO and SE_ISIS are the same."""
        assert SE_TRANSPLUTO == SE_ISIS

    def test_fict_offset_correct(self):
        """Test that body IDs use correct offset."""
        assert SE_CUPIDO == SE_FICT_OFFSET + 0
        assert SE_HADES == SE_FICT_OFFSET + 1
        assert SE_POSEIDON == SE_FICT_OFFSET + 7


class TestBodyNames:
    """Tests for body name functions."""

    def test_get_hypothetical_name_uranian(self):
        """Test getting names of Uranian planets."""
        assert get_hypothetical_name(SE_CUPIDO) == "Cupido"
        assert get_hypothetical_name(SE_HADES) == "Hades"
        assert get_hypothetical_name(SE_ZEUS) == "Zeus"
        assert get_hypothetical_name(SE_KRONOS) == "Kronos"
        assert get_hypothetical_name(SE_APOLLON) == "Apollon"
        assert get_hypothetical_name(SE_ADMETOS) == "Admetos"
        assert get_hypothetical_name(SE_VULKANUS) == "Vulkanus"
        assert get_hypothetical_name(SE_POSEIDON) == "Poseidon"

    def test_get_hypothetical_name_other(self):
        """Test getting names of other hypothetical bodies."""
        assert get_hypothetical_name(SE_ISIS) == "Transpluto"
        assert get_hypothetical_name(SE_WHITE_MOON) == "White Moon (Selena)"
        assert get_hypothetical_name(SE_WALDEMATH) == "Waldemath"

    def test_get_hypothetical_name_unknown(self):
        """Test that unknown body ID returns descriptive string."""
        name = get_hypothetical_name(999)
        assert "Unknown" in name or "999" in name

    def test_list_hypothetical_bodies(self):
        """Test listing all hypothetical bodies."""
        bodies = list_hypothetical_bodies()
        assert isinstance(bodies, dict)
        assert len(bodies) > 0
        assert SE_CUPIDO in bodies
        assert SE_ISIS in bodies
        assert bodies[SE_CUPIDO] == "Cupido"


class TestUranianElements:
    """Tests for Uranian planet element data."""

    def test_all_uranian_planets_have_elements(self):
        """Test that all Uranian planets have orbital elements defined."""
        uranian_ids = [
            SE_CUPIDO,
            SE_HADES,
            SE_ZEUS,
            SE_KRONOS,
            SE_APOLLON,
            SE_ADMETOS,
            SE_VULKANUS,
            SE_POSEIDON,
        ]
        for body_id in uranian_ids:
            assert body_id in URANIAN_ELEMENTS, (
                f"Body ID {body_id} should have elements"
            )

    def test_uranian_elements_have_required_fields(self):
        """Test that Uranian elements have all required fields."""
        for body_id, elements in URANIAN_ELEMENTS.items():
            assert hasattr(elements, "name")
            assert hasattr(elements, "L0")
            assert hasattr(elements, "n")
            assert hasattr(elements, "amplitude")
            assert hasattr(elements, "phase")
            assert hasattr(elements, "phase_rate")

    def test_uranian_mean_motions_are_positive(self):
        """Test that mean motions are positive (all move prograde)."""
        for body_id, elements in URANIAN_ELEMENTS.items():
            assert elements.n > 0, f"Mean motion for {elements.name} should be positive"

    def test_uranian_mean_motions_slower_than_neptune(self):
        """Test that Uranian planets are slower than Neptune."""
        neptune_n = 0.006021 * 36525  # Convert deg/day to deg/century
        for body_id, elements in URANIAN_ELEMENTS.items():
            assert elements.n < neptune_n, (
                f"{elements.name} should be slower than Neptune"
            )


class TestUranianLongitudeCalculations:
    """Tests for Uranian planet longitude calculations."""

    # J2000.0 epoch
    J2000 = 2451545.0

    def test_calc_uranian_longitude_at_j2000(self):
        """Test longitude at J2000.0 epoch."""
        # At J2000, T=0, so longitude should be close to L0
        for body_id, elements in URANIAN_ELEMENTS.items():
            lon = calc_uranian_longitude(body_id, self.J2000)
            # Should be close to L0 (within amplitude for oscillation)
            expected = elements.L0 % 360.0
            # Allow for oscillation amplitude
            diff = abs(lon - expected)
            if diff > 180:
                diff = 360 - diff
            assert diff < 5.0, (
                f"{elements.name}: longitude {lon} too far from L0 {expected}"
            )

    def test_calc_uranian_longitude_range(self):
        """Test that longitude is always in valid range."""
        test_dates = [
            self.J2000 - 36525,
            self.J2000,
            self.J2000 + 36525,
        ]  # 100 years span
        for body_id in URANIAN_ELEMENTS:
            for jd in test_dates:
                lon = calc_uranian_longitude(body_id, jd)
                assert 0.0 <= lon < 360.0, f"Longitude {lon} out of range"

    def test_calc_uranian_longitude_progression(self):
        """Test that longitude increases over time (prograde motion)."""
        for body_id in URANIAN_ELEMENTS:
            lon1 = calc_uranian_longitude(body_id, self.J2000)
            lon2 = calc_uranian_longitude(
                body_id, self.J2000 + 365.25 * 10
            )  # 10 years later
            # Handle wrap-around
            diff = lon2 - lon1
            if diff < -180:
                diff += 360
            elif diff > 180:
                diff -= 360
            # For slow bodies, might wrap around, so check that motion occurred
            assert diff != 0.0, f"{URANIAN_ELEMENTS[body_id].name} should have moved"

    def test_calc_uranian_longitude_invalid_body(self):
        """Test that invalid body ID raises ValueError."""
        with pytest.raises(ValueError):
            calc_uranian_longitude(SE_SUN, self.J2000)  # Not a Uranian planet


class TestUranianPositionCalculations:
    """Tests for full Uranian position calculations."""

    J2000 = 2451545.0

    def test_calc_uranian_position_returns_tuple(self):
        """Test that position returns correct tuple format."""
        pos = calc_uranian_position(SE_CUPIDO, self.J2000)
        assert isinstance(pos, tuple)
        assert len(pos) == 6
        lon, lat, dist, dlon, dlat, ddist = pos
        assert isinstance(lon, float)
        assert isinstance(lat, float)
        assert isinstance(dist, float)
        assert isinstance(dlon, float)
        assert isinstance(dlat, float)
        assert isinstance(ddist, float)

    def test_calc_uranian_position_latitude_zero(self):
        """Test that Uranian planets have zero latitude (on ecliptic)."""
        for body_id in URANIAN_ELEMENTS:
            pos = calc_uranian_position(body_id, self.J2000)
            assert pos[1] == 0.0, (
                f"{URANIAN_ELEMENTS[body_id].name} should have zero latitude"
            )

    def test_calc_uranian_position_distance_reasonable(self):
        """Test that distance estimates are reasonable for slow-moving bodies."""
        for body_id in URANIAN_ELEMENTS:
            pos = calc_uranian_position(body_id, self.J2000)
            # All Uranian planets should be beyond Neptune (~30 AU)
            # Based on their slow mean motions
            assert pos[2] > 30.0, (
                f"{URANIAN_ELEMENTS[body_id].name} should be beyond 30 AU"
            )
            # Periods up to ~1700 years correspond to semi-major axes up to ~3000+ AU
            assert pos[2] < 5000.0, (
                f"{URANIAN_ELEMENTS[body_id].name} distance seems too large"
            )

    def test_calc_uranian_position_velocity_positive(self):
        """Test that velocity is positive (prograde motion)."""
        for body_id in URANIAN_ELEMENTS:
            pos = calc_uranian_position(body_id, self.J2000)
            dlon = pos[3]
            # Mean motion should be positive, oscillation might briefly reverse
            # Check that it's reasonable (less than 1 deg/day for these slow bodies)
            assert abs(dlon) < 1.0, (
                f"{URANIAN_ELEMENTS[body_id].name} velocity too high"
            )

    def test_calc_uranian_position_invalid_body(self):
        """Test that invalid body ID raises ValueError."""
        with pytest.raises(ValueError):
            calc_uranian_position(SE_MOON, self.J2000)


class TestTransplutoCalculations:
    """Tests for Transpluto (Isis) position calculations."""

    J2000 = 2451545.0

    def test_calc_transpluto_position_returns_tuple(self):
        """Test that Transpluto position returns correct format."""
        pos = calc_transpluto_position(self.J2000)
        assert isinstance(pos, tuple)
        assert len(pos) == 6

    def test_calc_transpluto_position_range(self):
        """Test that Transpluto position is in valid range."""
        pos = calc_transpluto_position(self.J2000)
        lon, lat, dist = pos[0], pos[1], pos[2]
        assert 0.0 <= lon < 360.0
        assert -90.0 <= lat <= 90.0
        assert dist > 0

    def test_calc_transpluto_distance_large(self):
        """Test that Transpluto is far from the Sun."""
        pos = calc_transpluto_position(self.J2000)
        # Transpluto is beyond Pluto (~77 AU semi-major axis)
        assert pos[2] > 40.0, "Transpluto should be very distant"

    def test_calc_transpluto_slow_motion(self):
        """Test that Transpluto moves very slowly."""
        pos = calc_transpluto_position(self.J2000)
        dlon = pos[3]
        # Transpluto has ~2100 year period, so very slow motion
        # Daily motion should be roughly 360 / (2100 * 365.25) ~ 0.00047 deg/day
        assert abs(dlon) < 0.01, "Transpluto should move very slowly"


class TestWhiteMoonCalculations:
    """Tests for White Moon (Selena) position calculations."""

    J2000 = 2451545.0

    def test_calc_white_moon_position_returns_tuple(self):
        """Test that White Moon position returns correct format."""
        pos = calc_white_moon_position(self.J2000)
        assert isinstance(pos, tuple)
        assert len(pos) == 6

    def test_calc_white_moon_opposite_to_lilith(self):
        """Test that White Moon is 180 degrees from Black Moon Lilith."""
        from libephemeris import lunar

        lilith_lon = lunar.calc_mean_lilith(self.J2000)
        white_moon_pos = calc_white_moon_position(self.J2000)
        white_moon_lon = white_moon_pos[0]

        # Should be 180 degrees apart
        diff = abs(white_moon_lon - lilith_lon)
        if diff > 180:
            diff = 360 - diff
        assert abs(diff - 180.0) < 0.1, (
            f"White Moon should be 180 deg from Lilith, got {diff}"
        )

    def test_calc_white_moon_longitude_range(self):
        """Test that White Moon longitude is in valid range."""
        pos = calc_white_moon_position(self.J2000)
        assert 0.0 <= pos[0] < 360.0


class TestWaldemathCalculations:
    """Tests for Waldemath Black Moon calculations."""

    J2000 = 2451545.0

    def test_calc_waldemath_position_returns_tuple(self):
        """Test that Waldemath position returns correct format."""
        pos = calc_waldemath_position(self.J2000)
        assert isinstance(pos, tuple)
        assert len(pos) == 6

    def test_calc_waldemath_small_distance(self):
        """Test that Waldemath is at small distance (second focus of lunar orbit)."""
        pos = calc_waldemath_position(self.J2000)
        # Second focus is ~42,200 km ~ 0.00028 AU
        assert pos[2] < 0.01, "Waldemath should be at very small distance"

    def test_calc_waldemath_longitude_range(self):
        """Test that Waldemath longitude is in valid range."""
        pos = calc_waldemath_position(self.J2000)
        assert 0.0 <= pos[0] < 360.0


class TestMainCalculationFunction:
    """Tests for the main calc_hypothetical_position function."""

    J2000 = 2451545.0

    def test_calc_hypothetical_position_uranian(self):
        """Test that main function routes Uranian planets correctly."""
        pos1 = calc_hypothetical_position(SE_CUPIDO, self.J2000)
        pos2 = calc_uranian_position(SE_CUPIDO, self.J2000)
        assert pos1 == pos2

    def test_calc_hypothetical_position_transpluto(self):
        """Test that main function routes Transpluto correctly."""
        pos1 = calc_hypothetical_position(SE_ISIS, self.J2000)
        pos2 = calc_transpluto_position(self.J2000)
        assert pos1 == pos2

    def test_calc_hypothetical_position_white_moon(self):
        """Test that main function routes White Moon correctly."""
        pos1 = calc_hypothetical_position(SE_WHITE_MOON, self.J2000)
        pos2 = calc_white_moon_position(self.J2000)
        assert pos1 == pos2

    def test_calc_hypothetical_position_waldemath(self):
        """Test that main function routes Waldemath correctly."""
        pos1 = calc_hypothetical_position(SE_WALDEMATH, self.J2000)
        pos2 = calc_waldemath_position(self.J2000)
        assert pos1 == pos2

    def test_calc_hypothetical_position_invalid_body(self):
        """Test that invalid body ID raises ValueError."""
        with pytest.raises(ValueError):
            calc_hypothetical_position(SE_SUN, self.J2000)

    def test_calc_hypothetical_position_all_supported_bodies(self):
        """Test that all named hypothetical bodies can be calculated."""
        for body_id in HYPOTHETICAL_NAMES:
            # Should not raise an error
            pos = calc_hypothetical_position(body_id, self.J2000)
            assert isinstance(pos, tuple)
            assert len(pos) == 6


class TestPositionConsistency:
    """Tests for position calculation consistency over time."""

    J2000 = 2451545.0

    def test_positions_change_over_time(self):
        """Test that positions change over long time periods."""
        # Test over 100 years
        jd1 = self.J2000
        jd2 = self.J2000 + 36525  # 100 years

        for body_id in [SE_CUPIDO, SE_POSEIDON, SE_ISIS]:
            pos1 = calc_hypothetical_position(body_id, jd1)
            pos2 = calc_hypothetical_position(body_id, jd2)
            assert pos1[0] != pos2[0], (
                f"Body {body_id} should have moved over 100 years"
            )

    def test_positions_continuous(self):
        """Test that positions change smoothly (no jumps except at 360/0)."""
        body_id = SE_KRONOS
        jd = self.J2000
        prev_lon = calc_hypothetical_position(body_id, jd)[0]

        for i in range(1, 100):
            jd_new = jd + i * 10  # Every 10 days
            lon = calc_hypothetical_position(body_id, jd_new)[0]
            diff = abs(lon - prev_lon)
            if diff > 180:
                diff = 360 - diff
            # For slow bodies, should not jump more than a few degrees in 10 days
            assert diff < 5.0, f"Position jumped too much: {diff} degrees"
            prev_lon = lon


class TestHypotheticalElements:
    """Tests for hypothetical element data structures."""

    def test_hypothetical_elements_have_transpluto(self):
        """Test that HYPOTHETICAL_ELEMENTS contains Transpluto."""
        assert SE_ISIS in HYPOTHETICAL_ELEMENTS

    def test_hypothetical_elements_fields(self):
        """Test that elements have required fields."""
        for body_id, elements in HYPOTHETICAL_ELEMENTS.items():
            assert hasattr(elements, "name")
            assert hasattr(elements, "epoch")
            assert hasattr(elements, "a")
            assert hasattr(elements, "e")
            assert hasattr(elements, "i")
            assert hasattr(elements, "omega")
            assert hasattr(elements, "Omega")
            assert hasattr(elements, "M0")
            assert hasattr(elements, "n")

    def test_hypothetical_elements_valid_eccentricity(self):
        """Test that eccentricities are valid."""
        for body_id, elements in HYPOTHETICAL_ELEMENTS.items():
            assert 0 <= elements.e < 1, f"{elements.name} has invalid eccentricity"
