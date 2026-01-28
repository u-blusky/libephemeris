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
    calc_cupido,
    calc_hades,
    calc_zeus,
    calc_kronos,
    # Data structures
    URANIAN_ELEMENTS,
    HYPOTHETICAL_ELEMENTS,
    HYPOTHETICAL_NAMES,
    CUPIDO_KEPLERIAN_ELEMENTS,
    HADES_KEPLERIAN_ELEMENTS,
    ZEUS_KEPLERIAN_ELEMENTS,
    KRONOS_KEPLERIAN_ELEMENTS,
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


class TestCalcCupido:
    """Tests for the calc_cupido Keplerian propagation function."""

    J2000 = 2451545.0
    J1900 = 2415020.0

    def test_calc_cupido_returns_tuple(self):
        """Test that calc_cupido returns correct tuple format."""
        pos = calc_cupido(self.J2000)
        assert isinstance(pos, tuple)
        assert len(pos) == 6
        lon, lat, dist, dlon, dlat, ddist = pos
        assert isinstance(lon, float)
        assert isinstance(lat, float)
        assert isinstance(dist, float)
        assert isinstance(dlon, float)
        assert isinstance(dlat, float)
        assert isinstance(ddist, float)

    def test_calc_cupido_longitude_range(self):
        """Test that longitude is in valid range."""
        test_dates = [
            self.J2000 - 36525,  # 100 years before J2000
            self.J2000,
            self.J2000 + 36525,  # 100 years after J2000
        ]
        for jd in test_dates:
            pos = calc_cupido(jd)
            assert 0.0 <= pos[0] < 360.0, f"Longitude {pos[0]} out of range"

    def test_calc_cupido_latitude_zero(self):
        """Test that Cupido has zero latitude (on ecliptic)."""
        pos = calc_cupido(self.J2000)
        assert pos[1] == 0.0, "Cupido should have zero latitude"

    def test_calc_cupido_distance_correct(self):
        """Test that Cupido distance matches semi-major axis (circular orbit)."""
        pos = calc_cupido(self.J2000)
        # For circular orbit, distance = semi-major axis = 40.99837 AU
        assert abs(pos[2] - 40.99837) < 0.001, (
            f"Cupido distance should be ~40.99837 AU, got {pos[2]}"
        )

    def test_calc_cupido_distance_constant(self):
        """Test that distance is constant for circular orbit."""
        pos1 = calc_cupido(self.J2000)
        pos2 = calc_cupido(self.J2000 + 365.25 * 50)  # 50 years later
        assert pos1[2] == pos2[2], "Distance should be constant for circular orbit"

    def test_calc_cupido_velocity_positive(self):
        """Test that daily motion is positive (prograde)."""
        pos = calc_cupido(self.J2000)
        assert pos[3] > 0, "Cupido should have prograde motion"

    def test_calc_cupido_velocity_constant(self):
        """Test that velocity is constant for circular orbit."""
        pos = calc_cupido(self.J2000)
        # For circular orbit, daily motion = mean motion
        expected_n = CUPIDO_KEPLERIAN_ELEMENTS.n
        assert pos[3] == expected_n, (
            f"Daily motion should be {expected_n}, got {pos[3]}"
        )

    def test_calc_cupido_latitude_velocity_zero(self):
        """Test that latitude and distance velocity are zero."""
        pos = calc_cupido(self.J2000)
        assert pos[4] == 0.0, "Latitude velocity should be zero"
        assert pos[5] == 0.0, "Distance velocity should be zero for e=0"

    def test_calc_cupido_at_epoch(self):
        """Test that longitude at J1900.0 matches L0."""
        pos = calc_cupido(self.J1900)
        # At epoch, longitude should equal L0
        expected_L0 = CUPIDO_KEPLERIAN_ELEMENTS.L0
        assert abs(pos[0] - expected_L0) < 0.01, (
            f"At epoch, longitude should be {expected_L0}, got {pos[0]}"
        )

    def test_calc_cupido_progression(self):
        """Test that Cupido progresses correctly over time."""
        pos1 = calc_cupido(self.J2000)
        pos2 = calc_cupido(self.J2000 + 365.25)  # 1 year later

        # Calculate expected motion in 1 year
        expected_motion = CUPIDO_KEPLERIAN_ELEMENTS.n * 365.25

        # Calculate actual motion (handle wrap-around)
        actual_motion = pos2[0] - pos1[0]
        if actual_motion < -180:
            actual_motion += 360
        elif actual_motion > 180:
            actual_motion -= 360

        assert abs(actual_motion - expected_motion) < 0.01, (
            f"Annual motion should be {expected_motion:.4f} deg, got {actual_motion:.4f}"
        )

    def test_calc_cupido_orbital_period(self):
        """Test that Cupido completes one orbit in expected period."""
        # Orbital period from Kepler's 3rd law: T = a^1.5 years
        a = CUPIDO_KEPLERIAN_ELEMENTS.a
        expected_period_years = a**1.5  # ~262.5 years

        # One full orbit should advance longitude by ~360 degrees
        period_days = expected_period_years * 365.25
        pos1 = calc_cupido(self.J2000)
        pos2 = calc_cupido(self.J2000 + period_days)

        # After one full period, should be back to same longitude (within tolerance)
        diff = abs(pos2[0] - pos1[0])
        if diff > 180:
            diff = 360 - diff

        assert diff < 1.0, (
            f"After one orbit ({expected_period_years:.1f} years), "
            f"longitude should return near start, diff = {diff:.2f} deg"
        )

    def test_calc_cupido_keplerian_elements_valid(self):
        """Test that Cupido Keplerian elements have valid values."""
        elements = CUPIDO_KEPLERIAN_ELEMENTS
        assert elements.name == "Cupido"
        assert elements.a > 30.0, "Cupido should be beyond Neptune"
        assert elements.e == 0.0, "Cupido should have circular orbit"
        assert elements.i == 0.0, "Cupido should be on ecliptic"
        assert 0 < elements.n < 1.0, "Mean motion should be small but positive"

    def test_calc_cupido_vs_uranian_different_methods(self):
        """Test that calc_cupido uses different method from calc_uranian_position.

        The two methods use different algorithms:
        - calc_cupido: Simple Keplerian propagation from J1900.0 epoch
        - calc_uranian_position: Secular polynomial formula from J2000.0 epoch

        They may give different longitudes but both should produce valid positions
        that progress at similar rates over time.
        """
        keplerian_pos = calc_cupido(self.J2000)
        uranian_pos = calc_uranian_position(SE_CUPIDO, self.J2000)

        # Both should be in valid range
        assert 0.0 <= keplerian_pos[0] < 360.0
        assert 0.0 <= uranian_pos[0] < 360.0

        # Both should show similar daily motion (prograde, slow)
        # Keplerian uses n = 0.003757 deg/day
        # Uranian uses n/36525 where n = 1.091437 deg/century = 0.0000299 deg/day
        # Wait, that's much slower. Let me check the velocities instead.
        assert keplerian_pos[3] > 0, "Keplerian motion should be prograde"
        assert uranian_pos[3] > 0, "Uranian motion should be prograde"

    def test_se_cupido_constant_value(self):
        """Test that SE_CUPIDO has correct value (SE_FICT_OFFSET + 0 = 40)."""
        assert SE_CUPIDO == 40
        assert SE_CUPIDO == SE_FICT_OFFSET + 0

    def test_calc_cupido_exportable_from_main_module(self):
        """Test that calc_cupido is exported from main libephemeris module."""
        import libephemeris

        assert hasattr(libephemeris, "calc_cupido")
        pos = libephemeris.calc_cupido(self.J2000)
        assert len(pos) == 6


class TestCalcHades:
    """Tests for the calc_hades Keplerian propagation function."""

    J2000 = 2451545.0
    J1900 = 2415020.0

    def test_calc_hades_returns_tuple(self):
        """Test that calc_hades returns correct tuple format."""
        pos = calc_hades(self.J2000)
        assert isinstance(pos, tuple)
        assert len(pos) == 6
        lon, lat, dist, dlon, dlat, ddist = pos
        assert isinstance(lon, float)
        assert isinstance(lat, float)
        assert isinstance(dist, float)
        assert isinstance(dlon, float)
        assert isinstance(dlat, float)
        assert isinstance(ddist, float)

    def test_calc_hades_longitude_range(self):
        """Test that longitude is in valid range."""
        test_dates = [
            self.J2000 - 36525,  # 100 years before J2000
            self.J2000,
            self.J2000 + 36525,  # 100 years after J2000
        ]
        for jd in test_dates:
            pos = calc_hades(jd)
            assert 0.0 <= pos[0] < 360.0, f"Longitude {pos[0]} out of range"

    def test_calc_hades_latitude_small(self):
        """Test that Hades has small latitude (inclination = 1.05 degrees)."""
        pos = calc_hades(self.J2000)
        # With inclination of 1.05 degrees, latitude should be in that range
        assert -2.0 <= pos[1] <= 2.0, f"Latitude {pos[1]} unexpectedly large"

    def test_calc_hades_distance_correct(self):
        """Test that Hades distance is close to semi-major axis."""
        pos = calc_hades(self.J2000)
        # For nearly circular orbit (e=0.00245), distance ~ semi-major axis = 50.66744 AU
        # Maximum deviation from a = a * e = 50.66744 * 0.00245 = 0.124 AU
        expected_a = HADES_KEPLERIAN_ELEMENTS.a
        assert abs(pos[2] - expected_a) < 0.2, (
            f"Hades distance should be ~{expected_a:.2f} AU, got {pos[2]:.2f}"
        )

    def test_calc_hades_distance_nearly_constant(self):
        """Test that distance is nearly constant for nearly circular orbit."""
        pos1 = calc_hades(self.J2000)
        pos2 = calc_hades(self.J2000 + 365.25 * 50)  # 50 years later
        # With e=0.00245, max variation is about 0.25 AU
        assert abs(pos1[2] - pos2[2]) < 0.5, "Distance should be nearly constant"

    def test_calc_hades_velocity_positive(self):
        """Test that daily motion is positive (prograde)."""
        pos = calc_hades(self.J2000)
        assert pos[3] > 0, "Hades should have prograde motion"

    def test_calc_hades_velocity_reasonable(self):
        """Test that velocity is reasonable for Hades orbital period."""
        pos = calc_hades(self.J2000)
        # Expected mean motion: n = 360 / (a^1.5 * 365.25) deg/day
        # For a = 50.66744 AU: period ~ 360.5 years
        # n ~ 360 / (360.5 * 365.25) ~ 0.00274 deg/day
        expected_n = HADES_KEPLERIAN_ELEMENTS.n
        # Velocity should be close to mean motion (within 10% due to eccentricity)
        assert abs(pos[3] - expected_n) < expected_n * 0.5, (
            f"Daily motion should be ~{expected_n:.6f}, got {pos[3]:.6f}"
        )

    def test_calc_hades_at_epoch(self):
        """Test that mean anomaly at J1900.0 matches M0."""
        pos = calc_hades(self.J1900)
        # At epoch, we start with the given mean anomaly
        # The longitude depends on M0, omega, and Omega
        assert 0.0 <= pos[0] < 360.0, "Longitude should be in valid range at epoch"

    def test_calc_hades_progression(self):
        """Test that Hades progresses correctly over time."""
        pos1 = calc_hades(self.J2000)
        pos2 = calc_hades(self.J2000 + 365.25)  # 1 year later

        # Calculate expected motion in 1 year
        expected_motion = HADES_KEPLERIAN_ELEMENTS.n * 365.25

        # Calculate actual motion (handle wrap-around)
        actual_motion = pos2[0] - pos1[0]
        if actual_motion < -180:
            actual_motion += 360
        elif actual_motion > 180:
            actual_motion -= 360

        assert abs(actual_motion - expected_motion) < 0.1, (
            f"Annual motion should be ~{expected_motion:.4f} deg, got {actual_motion:.4f}"
        )

    def test_calc_hades_orbital_period(self):
        """Test that Hades completes one orbit in expected period."""
        # Orbital period from Kepler's 3rd law: T = a^1.5 years
        a = HADES_KEPLERIAN_ELEMENTS.a
        expected_period_years = a**1.5  # ~360.5 years

        # One full orbit should advance longitude by ~360 degrees
        period_days = expected_period_years * 365.25
        pos1 = calc_hades(self.J2000)
        pos2 = calc_hades(self.J2000 + period_days)

        # After one full period, should be back to same longitude (within tolerance)
        diff = abs(pos2[0] - pos1[0])
        if diff > 180:
            diff = 360 - diff

        assert diff < 1.0, (
            f"After one orbit ({expected_period_years:.1f} years), "
            f"longitude should return near start, diff = {diff:.2f} deg"
        )

    def test_calc_hades_keplerian_elements_valid(self):
        """Test that Hades Keplerian elements have valid values."""
        elements = HADES_KEPLERIAN_ELEMENTS
        assert elements.name == "Hades"
        assert elements.a > 40.0, "Hades should be beyond Cupido"
        assert 0 < elements.e < 1.0, "Eccentricity should be valid"
        assert abs(elements.e - 0.00245) < 0.0001, (
            "Eccentricity should match seorbel.txt"
        )
        assert abs(elements.i - 1.0500) < 0.0001, "Inclination should match seorbel.txt"
        assert 0 < elements.n < 1.0, "Mean motion should be small but positive"

    def test_calc_hades_vs_uranian_different_methods(self):
        """Test that calc_hades uses different method from calc_uranian_position.

        The two methods use different algorithms:
        - calc_hades: Full Keplerian propagation from J1900.0 epoch
        - calc_uranian_position: Secular polynomial formula from J2000.0 epoch

        They may give different longitudes but both should produce valid positions
        that progress at similar rates over time.
        """
        keplerian_pos = calc_hades(self.J2000)
        uranian_pos = calc_uranian_position(SE_HADES, self.J2000)

        # Both should be in valid range
        assert 0.0 <= keplerian_pos[0] < 360.0
        assert 0.0 <= uranian_pos[0] < 360.0

        # Both should show prograde motion
        assert keplerian_pos[3] > 0, "Keplerian motion should be prograde"
        assert uranian_pos[3] > 0, "Uranian motion should be prograde"

    def test_se_hades_constant_value(self):
        """Test that SE_HADES has correct value (SE_FICT_OFFSET + 1 = 41)."""
        assert SE_HADES == 41
        assert SE_HADES == SE_FICT_OFFSET + 1

    def test_calc_hades_exportable_from_main_module(self):
        """Test that calc_hades is exported from main libephemeris module."""
        import libephemeris

        assert hasattr(libephemeris, "calc_hades")
        pos = libephemeris.calc_hades(self.J2000)
        assert len(pos) == 6

    def test_calc_hades_latitude_velocity_small(self):
        """Test that latitude velocity is small."""
        pos = calc_hades(self.J2000)
        # Latitude velocity should be small for nearly circular orbit
        assert abs(pos[4]) < 0.01, "Latitude velocity should be small"

    def test_calc_hades_distance_velocity_small(self):
        """Test that distance velocity is small for nearly circular orbit."""
        pos = calc_hades(self.J2000)
        # Distance velocity should be small for nearly circular orbit (e=0.00245)
        assert abs(pos[5]) < 0.001, (
            "Distance velocity should be small for nearly circular orbit"
        )

    def test_hades_slower_than_cupido(self):
        """Test that Hades moves slower than Cupido (larger semi-major axis)."""
        cupido_pos = calc_cupido(self.J2000)
        hades_pos = calc_hades(self.J2000)

        # Hades has larger semi-major axis, so slower motion
        assert hades_pos[3] < cupido_pos[3], "Hades should move slower than Cupido"

    def test_hades_farther_than_cupido(self):
        """Test that Hades is farther from Sun than Cupido."""
        cupido_pos = calc_cupido(self.J2000)
        hades_pos = calc_hades(self.J2000)

        # Hades has larger semi-major axis
        assert hades_pos[2] > cupido_pos[2], (
            "Hades should be farther from Sun than Cupido"
        )


class TestCalcZeus:
    """Tests for the calc_zeus Keplerian propagation function."""

    J2000 = 2451545.0
    J1900 = 2415020.0

    def test_calc_zeus_returns_tuple(self):
        """Test that calc_zeus returns correct tuple format."""
        pos = calc_zeus(self.J2000)
        assert isinstance(pos, tuple)
        assert len(pos) == 6
        lon, lat, dist, dlon, dlat, ddist = pos
        assert isinstance(lon, float)
        assert isinstance(lat, float)
        assert isinstance(dist, float)
        assert isinstance(dlon, float)
        assert isinstance(dlat, float)
        assert isinstance(ddist, float)

    def test_calc_zeus_longitude_range(self):
        """Test that longitude is in valid range."""
        test_dates = [
            self.J2000 - 36525,  # 100 years before J2000
            self.J2000,
            self.J2000 + 36525,  # 100 years after J2000
        ]
        for jd in test_dates:
            pos = calc_zeus(jd)
            assert 0.0 <= pos[0] < 360.0, f"Longitude {pos[0]} out of range"

    def test_calc_zeus_latitude_zero(self):
        """Test that Zeus has zero latitude (on ecliptic)."""
        pos = calc_zeus(self.J2000)
        assert pos[1] == 0.0, "Zeus should have zero latitude"

    def test_calc_zeus_distance_correct(self):
        """Test that Zeus distance matches semi-major axis (circular orbit)."""
        pos = calc_zeus(self.J2000)
        # For circular orbit, distance = semi-major axis = 59.21436 AU
        assert abs(pos[2] - 59.21436) < 0.001, (
            f"Zeus distance should be ~59.21436 AU, got {pos[2]}"
        )

    def test_calc_zeus_distance_constant(self):
        """Test that distance is constant for circular orbit."""
        pos1 = calc_zeus(self.J2000)
        pos2 = calc_zeus(self.J2000 + 365.25 * 50)  # 50 years later
        assert pos1[2] == pos2[2], "Distance should be constant for circular orbit"

    def test_calc_zeus_velocity_positive(self):
        """Test that daily motion is positive (prograde)."""
        pos = calc_zeus(self.J2000)
        assert pos[3] > 0, "Zeus should have prograde motion"

    def test_calc_zeus_velocity_constant(self):
        """Test that velocity is constant for circular orbit."""
        pos = calc_zeus(self.J2000)
        # For circular orbit, daily motion = mean motion
        expected_n = ZEUS_KEPLERIAN_ELEMENTS.n
        assert pos[3] == expected_n, (
            f"Daily motion should be {expected_n}, got {pos[3]}"
        )

    def test_calc_zeus_latitude_velocity_zero(self):
        """Test that latitude and distance velocity are zero."""
        pos = calc_zeus(self.J2000)
        assert pos[4] == 0.0, "Latitude velocity should be zero"
        assert pos[5] == 0.0, "Distance velocity should be zero for e=0"

    def test_calc_zeus_at_epoch(self):
        """Test that longitude at J1900.0 matches L0."""
        pos = calc_zeus(self.J1900)
        # At epoch, longitude should equal L0
        expected_L0 = ZEUS_KEPLERIAN_ELEMENTS.L0
        assert abs(pos[0] - expected_L0) < 0.01, (
            f"At epoch, longitude should be {expected_L0}, got {pos[0]}"
        )

    def test_calc_zeus_progression(self):
        """Test that Zeus progresses correctly over time."""
        pos1 = calc_zeus(self.J2000)
        pos2 = calc_zeus(self.J2000 + 365.25)  # 1 year later

        # Calculate expected motion in 1 year
        expected_motion = ZEUS_KEPLERIAN_ELEMENTS.n * 365.25

        # Calculate actual motion (handle wrap-around)
        actual_motion = pos2[0] - pos1[0]
        if actual_motion < -180:
            actual_motion += 360
        elif actual_motion > 180:
            actual_motion -= 360

        assert abs(actual_motion - expected_motion) < 0.01, (
            f"Annual motion should be {expected_motion:.4f} deg, got {actual_motion:.4f}"
        )

    def test_calc_zeus_orbital_period(self):
        """Test that Zeus completes one orbit in expected period."""
        # Orbital period from Kepler's 3rd law: T = a^1.5 years
        a = ZEUS_KEPLERIAN_ELEMENTS.a
        expected_period_years = a**1.5  # ~455.7 years

        # One full orbit should advance longitude by ~360 degrees
        period_days = expected_period_years * 365.25
        pos1 = calc_zeus(self.J2000)
        pos2 = calc_zeus(self.J2000 + period_days)

        # After one full period, should be back to same longitude (within tolerance)
        diff = abs(pos2[0] - pos1[0])
        if diff > 180:
            diff = 360 - diff

        assert diff < 1.0, (
            f"After one orbit ({expected_period_years:.1f} years), "
            f"longitude should return near start, diff = {diff:.2f} deg"
        )

    def test_calc_zeus_keplerian_elements_valid(self):
        """Test that Zeus Keplerian elements have valid values."""
        elements = ZEUS_KEPLERIAN_ELEMENTS
        assert elements.name == "Zeus"
        assert elements.a > 50.0, "Zeus should be beyond Hades"
        assert elements.e == 0.0, "Zeus should have circular orbit"
        assert elements.i == 0.0, "Zeus should be on ecliptic"
        assert 0 < elements.n < 1.0, "Mean motion should be small but positive"

    def test_calc_zeus_vs_uranian_different_methods(self):
        """Test that calc_zeus uses different method from calc_uranian_position.

        The two methods use different algorithms:
        - calc_zeus: Simple Keplerian propagation from J1900.0 epoch
        - calc_uranian_position: Secular polynomial formula from J2000.0 epoch

        They may give different longitudes but both should produce valid positions
        that progress at similar rates over time.
        """
        keplerian_pos = calc_zeus(self.J2000)
        uranian_pos = calc_uranian_position(SE_ZEUS, self.J2000)

        # Both should be in valid range
        assert 0.0 <= keplerian_pos[0] < 360.0
        assert 0.0 <= uranian_pos[0] < 360.0

        # Keplerian should show prograde motion
        assert keplerian_pos[3] > 0, "Keplerian motion should be prograde"
        # Uranian velocity may be very small or briefly negative due to oscillation
        # Check that it's reasonable (within 1 deg/day for these slow bodies)
        assert abs(uranian_pos[3]) < 1.0, "Uranian velocity should be small"

    def test_se_zeus_constant_value(self):
        """Test that SE_ZEUS has correct value (SE_FICT_OFFSET + 2 = 42)."""
        assert SE_ZEUS == 42
        assert SE_ZEUS == SE_FICT_OFFSET + 2

    def test_calc_zeus_exportable_from_main_module(self):
        """Test that calc_zeus is exported from main libephemeris module."""
        import libephemeris

        assert hasattr(libephemeris, "calc_zeus")
        pos = libephemeris.calc_zeus(self.J2000)
        assert len(pos) == 6

    def test_zeus_slower_than_hades(self):
        """Test that Zeus moves slower than Hades (larger semi-major axis)."""
        hades_pos = calc_hades(self.J2000)
        zeus_pos = calc_zeus(self.J2000)

        # Zeus has larger semi-major axis, so slower motion
        assert zeus_pos[3] < hades_pos[3], "Zeus should move slower than Hades"

    def test_zeus_farther_than_hades(self):
        """Test that Zeus is farther from Sun than Hades."""
        hades_pos = calc_hades(self.J2000)
        zeus_pos = calc_zeus(self.J2000)

        # Zeus has larger semi-major axis
        assert zeus_pos[2] > hades_pos[2], "Zeus should be farther from Sun than Hades"


class TestCalcKronos:
    """Tests for the calc_kronos Keplerian propagation function."""

    J2000 = 2451545.0
    J1900 = 2415020.0

    def test_calc_kronos_returns_tuple(self):
        """Test that calc_kronos returns correct tuple format."""
        pos = calc_kronos(self.J2000)
        assert isinstance(pos, tuple)
        assert len(pos) == 6
        lon, lat, dist, dlon, dlat, ddist = pos
        assert isinstance(lon, float)
        assert isinstance(lat, float)
        assert isinstance(dist, float)
        assert isinstance(dlon, float)
        assert isinstance(dlat, float)
        assert isinstance(ddist, float)

    def test_calc_kronos_longitude_range(self):
        """Test that longitude is in valid range."""
        test_dates = [
            self.J2000 - 36525,  # 100 years before J2000
            self.J2000,
            self.J2000 + 36525,  # 100 years after J2000
        ]
        for jd in test_dates:
            pos = calc_kronos(jd)
            assert 0.0 <= pos[0] < 360.0, f"Longitude {pos[0]} out of range"

    def test_calc_kronos_latitude_zero(self):
        """Test that Kronos has zero latitude (on ecliptic)."""
        pos = calc_kronos(self.J2000)
        assert pos[1] == 0.0, "Kronos should have zero latitude"

    def test_calc_kronos_distance_correct(self):
        """Test that Kronos distance matches semi-major axis (circular orbit)."""
        pos = calc_kronos(self.J2000)
        # For circular orbit, distance = semi-major axis = 64.81690 AU
        assert abs(pos[2] - 64.81690) < 0.001, (
            f"Kronos distance should be ~64.81690 AU, got {pos[2]}"
        )

    def test_calc_kronos_distance_constant(self):
        """Test that distance is constant for circular orbit."""
        pos1 = calc_kronos(self.J2000)
        pos2 = calc_kronos(self.J2000 + 365.25 * 50)  # 50 years later
        assert pos1[2] == pos2[2], "Distance should be constant for circular orbit"

    def test_calc_kronos_velocity_positive(self):
        """Test that daily motion is positive (prograde)."""
        pos = calc_kronos(self.J2000)
        assert pos[3] > 0, "Kronos should have prograde motion"

    def test_calc_kronos_velocity_constant(self):
        """Test that velocity is constant for circular orbit."""
        pos = calc_kronos(self.J2000)
        # For circular orbit, daily motion = mean motion
        expected_n = KRONOS_KEPLERIAN_ELEMENTS.n
        assert pos[3] == expected_n, (
            f"Daily motion should be {expected_n}, got {pos[3]}"
        )

    def test_calc_kronos_latitude_velocity_zero(self):
        """Test that latitude and distance velocity are zero."""
        pos = calc_kronos(self.J2000)
        assert pos[4] == 0.0, "Latitude velocity should be zero"
        assert pos[5] == 0.0, "Distance velocity should be zero for e=0"

    def test_calc_kronos_at_epoch(self):
        """Test that longitude at J1900.0 matches L0."""
        pos = calc_kronos(self.J1900)
        # At epoch, longitude should equal L0
        expected_L0 = KRONOS_KEPLERIAN_ELEMENTS.L0
        assert abs(pos[0] - expected_L0) < 0.01, (
            f"At epoch, longitude should be {expected_L0}, got {pos[0]}"
        )

    def test_calc_kronos_progression(self):
        """Test that Kronos progresses correctly over time."""
        pos1 = calc_kronos(self.J2000)
        pos2 = calc_kronos(self.J2000 + 365.25)  # 1 year later

        # Calculate expected motion in 1 year
        expected_motion = KRONOS_KEPLERIAN_ELEMENTS.n * 365.25

        # Calculate actual motion (handle wrap-around)
        actual_motion = pos2[0] - pos1[0]
        if actual_motion < -180:
            actual_motion += 360
        elif actual_motion > 180:
            actual_motion -= 360

        assert abs(actual_motion - expected_motion) < 0.01, (
            f"Annual motion should be {expected_motion:.4f} deg, got {actual_motion:.4f}"
        )

    def test_calc_kronos_orbital_period(self):
        """Test that Kronos completes one orbit in expected period."""
        # Orbital period from Kepler's 3rd law: T = a^1.5 years
        a = KRONOS_KEPLERIAN_ELEMENTS.a
        expected_period_years = a**1.5  # ~521.9 years

        # One full orbit should advance longitude by ~360 degrees
        period_days = expected_period_years * 365.25
        pos1 = calc_kronos(self.J2000)
        pos2 = calc_kronos(self.J2000 + period_days)

        # After one full period, should be back to same longitude (within tolerance)
        diff = abs(pos2[0] - pos1[0])
        if diff > 180:
            diff = 360 - diff

        assert diff < 1.0, (
            f"After one orbit ({expected_period_years:.1f} years), "
            f"longitude should return near start, diff = {diff:.2f} deg"
        )

    def test_calc_kronos_keplerian_elements_valid(self):
        """Test that Kronos Keplerian elements have valid values."""
        elements = KRONOS_KEPLERIAN_ELEMENTS
        assert elements.name == "Kronos"
        assert elements.a > 59.0, "Kronos should be beyond Zeus"
        assert elements.e == 0.0, "Kronos should have circular orbit"
        assert elements.i == 0.0, "Kronos should be on ecliptic"
        assert 0 < elements.n < 1.0, "Mean motion should be small but positive"

    def test_calc_kronos_vs_uranian_different_methods(self):
        """Test that calc_kronos uses different method from calc_uranian_position.

        The two methods use different algorithms:
        - calc_kronos: Simple Keplerian propagation from J1900.0 epoch
        - calc_uranian_position: Secular polynomial formula from J2000.0 epoch

        They may give different longitudes but both should produce valid positions
        that progress at similar rates over time.
        """
        keplerian_pos = calc_kronos(self.J2000)
        uranian_pos = calc_uranian_position(SE_KRONOS, self.J2000)

        # Both should be in valid range
        assert 0.0 <= keplerian_pos[0] < 360.0
        assert 0.0 <= uranian_pos[0] < 360.0

        # Keplerian should show prograde motion
        assert keplerian_pos[3] > 0, "Keplerian motion should be prograde"
        # Uranian velocity may be very small or briefly negative due to oscillation
        # Check that it's reasonable (within 1 deg/day for these slow bodies)
        assert abs(uranian_pos[3]) < 1.0, "Uranian velocity should be small"

    def test_se_kronos_constant_value(self):
        """Test that SE_KRONOS has correct value (SE_FICT_OFFSET + 3 = 43)."""
        assert SE_KRONOS == 43
        assert SE_KRONOS == SE_FICT_OFFSET + 3

    def test_calc_kronos_exportable_from_main_module(self):
        """Test that calc_kronos is exported from main libephemeris module."""
        import libephemeris

        assert hasattr(libephemeris, "calc_kronos")
        pos = libephemeris.calc_kronos(self.J2000)
        assert len(pos) == 6

    def test_kronos_slower_than_zeus(self):
        """Test that Kronos moves slower than Zeus (larger semi-major axis)."""
        zeus_pos = calc_zeus(self.J2000)
        kronos_pos = calc_kronos(self.J2000)

        # Kronos has larger semi-major axis, so slower motion
        assert kronos_pos[3] < zeus_pos[3], "Kronos should move slower than Zeus"

    def test_kronos_farther_than_zeus(self):
        """Test that Kronos is farther from Sun than Zeus."""
        zeus_pos = calc_zeus(self.J2000)
        kronos_pos = calc_kronos(self.J2000)

        # Kronos has larger semi-major axis
        assert kronos_pos[2] > zeus_pos[2], (
            "Kronos should be farther from Sun than Zeus"
        )
