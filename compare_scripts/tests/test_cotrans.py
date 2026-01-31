"""
Tests for the cotrans coordinate transformation function.

Tests verify that cotrans correctly transforms coordinates between
ecliptic and equatorial systems, matching pyswisseph's swe.cotrans() behavior.
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestCotransBasic:
    """Basic functionality tests for cotrans."""

    def test_cotrans_exported(self):
        """Test that cotrans is exported from the package."""
        assert hasattr(ephem, "cotrans")
        assert callable(ephem.cotrans)

    def test_cotrans_origin_point(self):
        """Test transformation at the origin (0, 0)."""
        # At origin, both systems should give (0, 0)
        obliquity = 23.4392911
        coord = (0.0, 0.0, 1.0)

        result = ephem.cotrans(coord, obliquity)

        assert abs(result[0]) < 1e-10 or abs(result[0] - 360.0) < 1e-10
        assert abs(result[1]) < 1e-10
        assert result[2] == 1.0  # Distance unchanged

    def test_cotrans_preserves_distance(self):
        """Test that distance is preserved in transformation."""
        obliquity = 23.4
        coord = (45.0, 30.0, 2.5)

        result = ephem.cotrans(coord, obliquity)

        assert result[2] == coord[2]

    def test_cotrans_returns_tuple(self):
        """Test that cotrans returns a tuple of three elements."""
        coord = (90.0, 0.0, 1.0)
        obliquity = 23.4

        result = ephem.cotrans(coord, obliquity)

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_cotrans_longitude_in_range(self):
        """Test that output longitude is in [0, 360) range."""
        test_cases = [
            ((0.0, 0.0, 1.0), 23.4),
            ((90.0, 0.0, 1.0), 23.4),
            ((180.0, 45.0, 1.0), 23.4),
            ((270.0, -30.0, 1.0), 23.4),
            ((359.9, 10.0, 1.0), 23.4),
        ]

        for coord, obl in test_cases:
            result = ephem.cotrans(coord, obl)
            assert 0.0 <= result[0] < 360.0, (
                f"Longitude {result[0]} out of range for input {coord}"
            )


class TestCotransVsSwisseph:
    """Comparison tests with pyswisseph's swe.cotrans()."""

    @pytest.mark.parametrize(
        "coord,obliquity",
        [
            ((0.0, 0.0, 1.0), 23.4392911),
            ((90.0, 0.0, 1.0), 23.4392911),
            ((180.0, 0.0, 1.0), 23.4392911),
            ((270.0, 0.0, 1.0), 23.4392911),
            ((45.0, 30.0, 1.0), 23.4392911),
            ((120.0, -15.0, 1.0), 23.4392911),
            ((200.0, 5.0, 2.5), 23.4392911),
            ((359.0, 0.0, 1.0), 23.4392911),
        ],
    )
    def test_cotrans_ecl_to_eq_matches_swisseph(self, coord, obliquity):
        """Test ecliptic to equatorial transformation matches pyswisseph."""
        result_lib = ephem.cotrans(coord, obliquity)
        result_swe = swe.cotrans(coord, obliquity)

        # Compare longitude (handle wrap-around)
        lon_diff = abs(result_lib[0] - result_swe[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 0.0001, (
            f"Longitude mismatch: lib={result_lib[0]}, swe={result_swe[0]}"
        )

        # Compare latitude
        assert abs(result_lib[1] - result_swe[1]) < 0.0001, (
            f"Latitude mismatch: lib={result_lib[1]}, swe={result_swe[1]}"
        )

        # Compare distance
        assert abs(result_lib[2] - result_swe[2]) < 1e-10, (
            f"Distance mismatch: lib={result_lib[2]}, swe={result_swe[2]}"
        )

    @pytest.mark.parametrize(
        "coord,obliquity",
        [
            ((0.0, 0.0, 1.0), -23.4392911),
            ((90.0, 0.0, 1.0), -23.4392911),
            ((180.0, 0.0, 1.0), -23.4392911),
            ((270.0, 0.0, 1.0), -23.4392911),
            ((45.0, 30.0, 1.0), -23.4392911),
            ((120.0, -15.0, 1.0), -23.4392911),
            ((200.0, 5.0, 2.5), -23.4392911),
            ((359.0, 0.0, 1.0), -23.4392911),
        ],
    )
    def test_cotrans_eq_to_ecl_matches_swisseph(self, coord, obliquity):
        """Test equatorial to ecliptic transformation matches pyswisseph."""
        result_lib = ephem.cotrans(coord, obliquity)
        result_swe = swe.cotrans(coord, obliquity)

        # Compare longitude (handle wrap-around)
        lon_diff = abs(result_lib[0] - result_swe[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 0.0001, (
            f"Longitude mismatch: lib={result_lib[0]}, swe={result_swe[0]}"
        )

        # Compare latitude
        assert abs(result_lib[1] - result_swe[1]) < 0.0001, (
            f"Latitude mismatch: lib={result_lib[1]}, swe={result_swe[1]}"
        )

        # Compare distance
        assert abs(result_lib[2] - result_swe[2]) < 1e-10, (
            f"Distance mismatch: lib={result_lib[2]}, swe={result_swe[2]}"
        )


class TestCotransRoundTrip:
    """Test that round-trip transformations return to the original coordinates."""

    @pytest.mark.parametrize(
        "coord",
        [
            (0.0, 0.0, 1.0),
            (45.0, 23.0, 1.0),
            (90.0, 0.0, 1.0),
            (120.0, -10.0, 1.5),
            (180.0, 45.0, 1.0),
            (270.0, -30.0, 2.0),
            (359.9, 5.0, 1.0),
        ],
    )
    def test_round_trip_ecl_eq_ecl(self, coord):
        """Test ecliptic -> equatorial -> ecliptic round trip."""
        obliquity = 23.4392911

        # Ecliptic to equatorial
        equatorial = ephem.cotrans(coord, obliquity)

        # Equatorial back to ecliptic
        result = ephem.cotrans(equatorial, -obliquity)

        # Compare with original
        lon_diff = abs(result[0] - coord[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 0.0001, (
            f"Longitude not preserved in round trip: {coord[0]} -> {result[0]}"
        )
        assert abs(result[1] - coord[1]) < 0.0001, (
            f"Latitude not preserved in round trip: {coord[1]} -> {result[1]}"
        )
        assert abs(result[2] - coord[2]) < 1e-10, (
            f"Distance not preserved: {coord[2]} -> {result[2]}"
        )


class TestCotransEdgeCases:
    """Edge case tests for cotrans."""

    def test_cotrans_at_poles(self):
        """Test transformation at polar coordinates."""
        obliquity = 23.4392911

        # Near north pole
        coord_north = (0.0, 89.0, 1.0)
        result = ephem.cotrans(coord_north, obliquity)
        swe_result = swe.cotrans(coord_north, obliquity)
        assert abs(result[1] - swe_result[1]) < 0.001

        # Near south pole
        coord_south = (0.0, -89.0, 1.0)
        result = ephem.cotrans(coord_south, obliquity)
        swe_result = swe.cotrans(coord_south, obliquity)
        assert abs(result[1] - swe_result[1]) < 0.001

    def test_cotrans_zero_obliquity(self):
        """Test with zero obliquity (no transformation)."""
        coord = (45.0, 30.0, 1.5)
        result = ephem.cotrans(coord, 0.0)

        # With zero obliquity, coordinates should be unchanged
        lon_diff = abs(result[0] - coord[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 1e-10
        assert abs(result[1] - coord[1]) < 1e-10
        assert result[2] == coord[2]

    def test_cotrans_various_obliquities(self):
        """Test with various obliquity values."""
        coord = (60.0, 15.0, 1.0)
        obliquities = [20.0, 23.0, 23.4392911, 25.0, 30.0]

        for obl in obliquities:
            result_lib = ephem.cotrans(coord, obl)
            result_swe = swe.cotrans(coord, obl)

            lon_diff = abs(result_lib[0] - result_swe[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            assert lon_diff < 0.0001, f"Failed for obliquity {obl}"
            assert abs(result_lib[1] - result_swe[1]) < 0.0001, (
                f"Failed for obliquity {obl}"
            )


class TestCotransRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_coordinates(self, random_longitudes):
        """Test with random longitude values."""
        obliquity = 23.4392911
        lons = random_longitudes(50)

        for lon in lons:
            lat = (lon % 60) - 30  # Generate a latitude based on longitude
            coord = (lon, lat, 1.0)

            result_lib = ephem.cotrans(coord, obliquity)
            result_swe = swe.cotrans(coord, obliquity)

            lon_diff = abs(result_lib[0] - result_swe[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            assert lon_diff < 0.0001, f"Longitude mismatch for input {coord}"
            assert abs(result_lib[1] - result_swe[1]) < 0.0001, (
                f"Latitude mismatch for input {coord}"
            )


class TestCotransSpBasic:
    """Basic functionality tests for cotrans_sp."""

    def test_cotrans_sp_exported(self):
        """Test that cotrans_sp is exported from the package."""
        assert hasattr(ephem, "cotrans_sp")
        assert callable(ephem.cotrans_sp)

    def test_cotrans_sp_origin_point(self):
        """Test transformation at the origin (0, 0)."""
        obliquity = 23.4392911
        coord = (0.0, 0.0, 1.0)
        speed = (1.0, 0.0, 0.0)

        coord_result, speed_result = ephem.cotrans_sp(coord, speed, obliquity)

        assert abs(coord_result[0]) < 1e-10 or abs(coord_result[0] - 360.0) < 1e-10
        assert abs(coord_result[1]) < 1e-10
        assert coord_result[2] == 1.0  # Distance unchanged

    def test_cotrans_sp_preserves_distance(self):
        """Test that distance and distance speed are preserved."""
        obliquity = 23.4
        coord = (45.0, 30.0, 2.5)
        speed = (1.0, 0.5, 0.1)

        coord_result, speed_result = ephem.cotrans_sp(coord, speed, obliquity)

        assert coord_result[2] == coord[2]
        assert speed_result[2] == speed[2]

    def test_cotrans_sp_returns_tuples(self):
        """Test that cotrans_sp returns two tuples of three elements each."""
        coord = (90.0, 0.0, 1.0)
        speed = (1.0, 0.0, 0.0)
        obliquity = 23.4

        result = ephem.cotrans_sp(coord, speed, obliquity)

        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], tuple)
        assert isinstance(result[1], tuple)
        assert len(result[0]) == 3
        assert len(result[1]) == 3

    def test_cotrans_sp_zero_speed(self):
        """Test with zero speed returns consistent position."""
        obliquity = 23.4392911
        coord = (90.0, 23.4, 1.0)
        speed = (0.0, 0.0, 0.0)

        coord_result, speed_result = ephem.cotrans_sp(coord, speed, obliquity)

        # Position should match cotrans
        cotrans_result = ephem.cotrans(coord, obliquity)
        lon_diff = abs(coord_result[0] - cotrans_result[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 1e-10
        assert abs(coord_result[1] - cotrans_result[1]) < 1e-10

        # Speed should be zero when input speed is zero
        assert abs(speed_result[0]) < 1e-10
        assert abs(speed_result[1]) < 1e-10
        assert abs(speed_result[2]) < 1e-10


class TestCotransSpVsSwisseph:
    """Comparison tests with pyswisseph's swe.cotrans_sp()."""

    @pytest.mark.parametrize(
        "coord,speed,obliquity",
        [
            ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 23.4392911),
            ((90.0, 0.0, 1.0), (0.5, 0.1, 0.0), 23.4392911),
            ((180.0, 0.0, 1.0), (1.0, 0.5, 0.0), 23.4392911),
            ((270.0, 0.0, 1.0), (0.0, 1.0, 0.0), 23.4392911),
            ((45.0, 30.0, 1.0), (0.8, 0.2, 0.1), 23.4392911),
            ((120.0, -15.0, 1.0), (1.2, -0.3, 0.0), 23.4392911),
            ((200.0, 5.0, 2.5), (0.5, 0.1, 0.05), 23.4392911),
            ((359.0, 0.0, 1.0), (1.0, 0.0, 0.0), 23.4392911),
        ],
    )
    def test_cotrans_sp_ecl_to_eq_matches_swisseph(self, coord, speed, obliquity):
        """Test ecliptic to equatorial transformation matches pyswisseph."""
        coord_result, speed_result = ephem.cotrans_sp(coord, speed, obliquity)

        # pyswisseph takes a 6-tuple and returns 6 values
        input_6 = coord + speed
        swe_result = swe.cotrans_sp(input_6, obliquity)

        # Compare position
        lon_diff = abs(coord_result[0] - swe_result[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 0.0001, (
            f"Longitude mismatch: lib={coord_result[0]}, swe={swe_result[0]}"
        )
        assert abs(coord_result[1] - swe_result[1]) < 0.0001, (
            f"Latitude mismatch: lib={coord_result[1]}, swe={swe_result[1]}"
        )
        assert abs(coord_result[2] - swe_result[2]) < 1e-10, (
            f"Distance mismatch: lib={coord_result[2]}, swe={swe_result[2]}"
        )

        # Compare speed
        assert abs(speed_result[0] - swe_result[3]) < 0.0001, (
            f"Lon speed mismatch: lib={speed_result[0]}, swe={swe_result[3]}"
        )
        assert abs(speed_result[1] - swe_result[4]) < 0.0001, (
            f"Lat speed mismatch: lib={speed_result[1]}, swe={swe_result[4]}"
        )
        assert abs(speed_result[2] - swe_result[5]) < 1e-10, (
            f"Dist speed mismatch: lib={speed_result[2]}, swe={swe_result[5]}"
        )

    @pytest.mark.parametrize(
        "coord,speed,obliquity",
        [
            ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), -23.4392911),
            ((90.0, 0.0, 1.0), (0.5, 0.1, 0.0), -23.4392911),
            ((180.0, 0.0, 1.0), (1.0, 0.5, 0.0), -23.4392911),
            ((270.0, 0.0, 1.0), (0.0, 1.0, 0.0), -23.4392911),
            ((45.0, 30.0, 1.0), (0.8, 0.2, 0.1), -23.4392911),
            ((120.0, -15.0, 1.0), (1.2, -0.3, 0.0), -23.4392911),
            ((200.0, 5.0, 2.5), (0.5, 0.1, 0.05), -23.4392911),
            ((359.0, 0.0, 1.0), (1.0, 0.0, 0.0), -23.4392911),
        ],
    )
    def test_cotrans_sp_eq_to_ecl_matches_swisseph(self, coord, speed, obliquity):
        """Test equatorial to ecliptic transformation matches pyswisseph."""
        coord_result, speed_result = ephem.cotrans_sp(coord, speed, obliquity)

        # pyswisseph takes a 6-tuple and returns 6 values
        input_6 = coord + speed
        swe_result = swe.cotrans_sp(input_6, obliquity)

        # Compare position
        lon_diff = abs(coord_result[0] - swe_result[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 0.0001, (
            f"Longitude mismatch: lib={coord_result[0]}, swe={swe_result[0]}"
        )
        assert abs(coord_result[1] - swe_result[1]) < 0.0001, (
            f"Latitude mismatch: lib={coord_result[1]}, swe={swe_result[1]}"
        )
        assert abs(coord_result[2] - swe_result[2]) < 1e-10, (
            f"Distance mismatch: lib={coord_result[2]}, swe={swe_result[2]}"
        )

        # Compare speed
        assert abs(speed_result[0] - swe_result[3]) < 0.0001, (
            f"Lon speed mismatch: lib={speed_result[0]}, swe={swe_result[3]}"
        )
        assert abs(speed_result[1] - swe_result[4]) < 0.0001, (
            f"Lat speed mismatch: lib={speed_result[1]}, swe={swe_result[4]}"
        )
        assert abs(speed_result[2] - swe_result[5]) < 1e-10, (
            f"Dist speed mismatch: lib={speed_result[2]}, swe={swe_result[5]}"
        )


class TestCotransSpRoundTrip:
    """Test that round-trip transformations return to the original values."""

    @pytest.mark.parametrize(
        "coord,speed",
        [
            ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
            ((45.0, 23.0, 1.0), (0.5, 0.2, 0.0)),
            ((90.0, 0.0, 1.0), (0.8, 0.1, 0.0)),
            ((120.0, -10.0, 1.5), (1.0, -0.5, 0.1)),
            ((180.0, 45.0, 1.0), (0.3, 0.3, 0.0)),
            ((270.0, -30.0, 2.0), (0.7, 0.0, 0.0)),
        ],
    )
    def test_round_trip_ecl_eq_ecl(self, coord, speed):
        """Test ecliptic -> equatorial -> ecliptic round trip."""
        obliquity = 23.4392911

        # Ecliptic to equatorial
        eq_coord, eq_speed = ephem.cotrans_sp(coord, speed, obliquity)

        # Equatorial back to ecliptic
        result_coord, result_speed = ephem.cotrans_sp(eq_coord, eq_speed, -obliquity)

        # Compare with original position
        lon_diff = abs(result_coord[0] - coord[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 0.0001, (
            f"Longitude not preserved in round trip: {coord[0]} -> {result_coord[0]}"
        )
        assert abs(result_coord[1] - coord[1]) < 0.0001, (
            f"Latitude not preserved in round trip: {coord[1]} -> {result_coord[1]}"
        )
        assert abs(result_coord[2] - coord[2]) < 1e-10, (
            f"Distance not preserved: {coord[2]} -> {result_coord[2]}"
        )

        # Compare with original speed
        assert abs(result_speed[0] - speed[0]) < 0.0001, (
            f"Lon speed not preserved: {speed[0]} -> {result_speed[0]}"
        )
        assert abs(result_speed[1] - speed[1]) < 0.0001, (
            f"Lat speed not preserved: {speed[1]} -> {result_speed[1]}"
        )
        assert abs(result_speed[2] - speed[2]) < 1e-10, (
            f"Dist speed not preserved: {speed[2]} -> {result_speed[2]}"
        )


class TestCotransSpPositionConsistency:
    """Test that cotrans_sp position matches cotrans."""

    @pytest.mark.parametrize(
        "coord,obliquity",
        [
            ((0.0, 0.0, 1.0), 23.4392911),
            ((90.0, 23.4, 1.0), 23.4392911),
            ((180.0, -15.0, 1.0), 23.4392911),
            ((270.0, 45.0, 1.0), 23.4392911),
            ((45.0, 30.0, 1.0), -23.4392911),
            ((120.0, -5.0, 1.0), -23.4392911),
        ],
    )
    def test_position_matches_cotrans(self, coord, obliquity):
        """Test that position from cotrans_sp matches cotrans."""
        speed = (1.0, 0.5, 0.1)

        coord_sp, _ = ephem.cotrans_sp(coord, speed, obliquity)
        coord_plain = ephem.cotrans(coord, obliquity)

        lon_diff = abs(coord_sp[0] - coord_plain[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 1e-10, (
            f"cotrans_sp lon {coord_sp[0]} != cotrans lon {coord_plain[0]}"
        )
        assert abs(coord_sp[1] - coord_plain[1]) < 1e-10, (
            f"cotrans_sp lat {coord_sp[1]} != cotrans lat {coord_plain[1]}"
        )
        assert abs(coord_sp[2] - coord_plain[2]) < 1e-10
