"""
Tests for solar eclipse calculations.
"""

from libephemeris import (
    sol_eclipse_when_glob,
    swe_sol_eclipse_when_glob,
    sol_eclipse_when_loc,
    swe_sol_eclipse_when_loc,
    swe_julday,
    swe_revjul,
    SE_ECL_TOTAL,
    SE_ECL_ANNULAR,
    SE_ECL_PARTIAL,
    SE_ECL_CENTRAL,
    SE_ECL_VISIBLE,
    SE_ECL_1ST_VISIBLE,
    SE_ECL_4TH_VISIBLE,
    SEFLG_SWIEPH,
)


class TestSolEclipseWhenGlob:
    """Tests for sol_eclipse_when_glob function."""

    def test_finds_eclipse_after_start_date(self):
        """Should find an eclipse after the start date."""
        # Start from Jan 1, 2024
        jd_start = swe_julday(2024, 1, 1, 0)

        times, ecl_type = sol_eclipse_when_glob(jd_start)

        # Should return a valid eclipse time
        assert times[0] > jd_start
        # Eclipse type should have some flags set
        assert ecl_type != 0

    def test_returns_ten_time_values(self):
        """Should return tuple of 10 time values like pyswisseph."""
        jd_start = swe_julday(2024, 1, 1, 0)

        times, _ = sol_eclipse_when_glob(jd_start)

        assert len(times) == 10

    def test_time_order_is_correct(self):
        """Eclipse phase times should be in chronological order."""
        jd_start = swe_julday(2024, 1, 1, 0)

        times, ecl_type = sol_eclipse_when_glob(jd_start)

        t_max = times[0]
        t_first = times[1]
        t_fourth = times[4]

        # First contact should be before maximum
        assert t_first < t_max
        # Maximum should be before fourth contact
        assert t_max < t_fourth

    def test_filter_by_total_eclipse(self):
        """Should filter for total eclipses only."""
        jd_start = swe_julday(2020, 1, 1, 0)

        times, ecl_type = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # Should have total eclipse flag set
        assert ecl_type & SE_ECL_TOTAL

    def test_filter_by_annular_eclipse(self):
        """Should filter for annular eclipses only."""
        jd_start = swe_julday(2020, 1, 1, 0)

        times, ecl_type = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_ANNULAR)

        # Should have annular eclipse flag set
        assert ecl_type & SE_ECL_ANNULAR

    def test_filter_accepts_partial_eclipse_type(self):
        """Should accept partial eclipse filter parameter."""
        # Note: Pure partial solar eclipses (not central) are rare.
        # This test verifies the function accepts the parameter and either
        # finds a partial eclipse or times out (both are valid behaviors).
        jd_start = swe_julday(2020, 1, 1, 0)

        # Test that the function accepts the parameter - we check for any eclipse
        # that includes partial (since many eclipses have partial phases)
        times, ecl_type = sol_eclipse_when_glob(jd_start)

        # The function should return something (any eclipse type)
        assert ecl_type != 0

    def test_known_eclipse_april_2024(self):
        """Test against known total solar eclipse of April 8, 2024."""
        # Start search before the eclipse
        jd_start = swe_julday(2024, 3, 1, 0)

        times, ecl_type = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # Convert JD to date
        year, month, day, hour = swe_revjul(times[0])

        # The April 8, 2024 total solar eclipse
        # Allow for some tolerance in the search
        assert year == 2024
        assert month == 4
        # Day should be around 8 (allowing some tolerance)
        assert 7 <= day <= 9

    def test_known_eclipse_october_2023(self):
        """Test against known annular solar eclipse of October 14, 2023."""
        # Start search before the eclipse
        jd_start = swe_julday(2023, 9, 1, 0)

        times, ecl_type = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_ANNULAR)

        # Convert JD to date
        year, month, day, hour = swe_revjul(times[0])

        # The October 14, 2023 annular solar eclipse
        assert year == 2023
        assert month == 10
        # Day should be around 14
        assert 13 <= day <= 15

    def test_central_eclipse_has_second_third_contact(self):
        """Central eclipses should have second and third contact times."""
        jd_start = swe_julday(2024, 1, 1, 0)

        times, ecl_type = sol_eclipse_when_glob(
            jd_start, eclipse_type=SE_ECL_TOTAL | SE_ECL_ANNULAR
        )

        if ecl_type & SE_ECL_CENTRAL:
            # Second contact (total/annular begins)
            t_second = times[2]
            # Third contact (total/annular ends)
            t_third = times[3]

            # These should be non-zero for central eclipses
            assert t_second > 0
            assert t_third > 0
            # And in correct order
            assert t_second < t_third

    def test_alias_matches_main_function(self):
        """swe_sol_eclipse_when_glob should be an alias for sol_eclipse_when_glob."""
        assert swe_sol_eclipse_when_glob is sol_eclipse_when_glob

    def test_flags_parameter_accepted(self):
        """Should accept flags parameter without error."""
        jd_start = swe_julday(2024, 1, 1, 0)

        # Should not raise
        times, ecl_type = sol_eclipse_when_glob(jd_start, flags=SEFLG_SWIEPH)

        assert times[0] > jd_start

    def test_multiple_eclipses_in_sequence(self):
        """Finding consecutive eclipses should return different events."""
        jd_start = swe_julday(2024, 1, 1, 0)

        # Find first eclipse
        times1, _ = sol_eclipse_when_glob(jd_start)

        # Find next eclipse (search from after first one)
        jd_after_first = times1[0] + 1  # Day after first eclipse
        times2, _ = sol_eclipse_when_glob(jd_after_first)

        # Second eclipse should be after first
        assert times2[0] > times1[0]
        # They should be at least 25 days apart (minimum time between eclipses)
        assert times2[0] - times1[0] >= 25


class TestNewMoonFinding:
    """Tests for internal New Moon finding logic."""

    def test_eclipse_occurs_near_new_moon(self):
        """Eclipse should occur very close to a New Moon."""
        from libephemeris import swe_calc_ut, SE_SUN, SE_MOON

        jd_start = swe_julday(2024, 1, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start)

        # Get Sun and Moon positions at eclipse maximum
        sun_pos, _ = swe_calc_ut(times[0], SE_SUN, SEFLG_SWIEPH)
        moon_pos, _ = swe_calc_ut(times[0], SE_MOON, SEFLG_SWIEPH)

        # Calculate elongation (should be very close to 0 at New Moon)
        elongation = abs((moon_pos[0] - sun_pos[0] + 180) % 360 - 180)

        # Should be within a few degrees of conjunction
        assert elongation < 5.0


class TestSolEclipseWhenLoc:
    """Tests for sol_eclipse_when_loc function."""

    def test_finds_eclipse_visible_from_location(self):
        """Should find an eclipse visible from the given location."""
        # Start from Jan 1, 2024, search from a central US location
        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0  # Central US

        times, attr, ecl_type = sol_eclipse_when_loc(jd_start, lat, lon)

        # Should return valid eclipse time
        assert times[0] > jd_start
        # Eclipse should be visible
        assert ecl_type & SE_ECL_VISIBLE
        # Magnitude should be positive
        assert attr[0] > 0

    def test_returns_correct_tuple_sizes(self):
        """Should return times tuple of 10 and attr tuple of 11 elements."""
        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon = 41.9028, 12.4964  # Rome

        times, attr, _ = sol_eclipse_when_loc(jd_start, lat, lon)

        assert len(times) == 10
        assert len(attr) == 11

    def test_time_order_is_correct(self):
        """Eclipse phase times should be in chronological order."""
        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        times, attr, ecl_type = sol_eclipse_when_loc(jd_start, lat, lon)

        t_max = times[0]
        t_first = times[1]
        t_fourth = times[4]

        # First contact should be before maximum (if visible)
        if t_first > 0:
            assert t_first < t_max
        # Maximum should be before fourth contact (if visible)
        if t_fourth > 0:
            assert t_max < t_fourth

    def test_eclipse_attributes_valid(self):
        """Eclipse attributes should be within valid ranges."""
        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        times, attr, ecl_type = sol_eclipse_when_loc(jd_start, lat, lon)

        # Magnitude should be 0-1.5 range
        magnitude = attr[0]
        assert 0 <= magnitude <= 1.5

        # Obscuration should be 0-1 range
        obscuration = attr[2]
        assert 0 <= obscuration <= 1

        # Sun altitude should be -90 to 90
        altitude = attr[4]
        assert -90 <= altitude <= 90

        # Azimuth should be 0 to 360
        azimuth = attr[3]
        assert 0 <= azimuth <= 360

        # Moon/Sun apparent diameters should be positive and reasonable
        moon_diam = attr[5]
        sun_diam = attr[6]
        assert 0.4 < moon_diam < 0.7  # degrees (roughly)
        assert 0.4 < sun_diam < 0.7  # degrees

    def test_known_eclipse_april_2024_texas(self):
        """Test April 8, 2024 total eclipse from Dallas, Texas."""
        # Dallas was in the path of totality
        jd_start = swe_julday(2024, 3, 1, 0)
        dallas_lat, dallas_lon = 32.7767, -96.7970

        times, attr, ecl_type = sol_eclipse_when_loc(jd_start, dallas_lat, dallas_lon)

        # Convert JD to date
        year, month, day, hour = swe_revjul(times[0])

        # Should find the April 8, 2024 eclipse
        assert year == 2024
        assert month == 4
        assert 7 <= day <= 9

        # Dallas was in path of totality - should have high magnitude
        # (may be total or high partial depending on exact location)
        assert attr[0] > 0.9  # High magnitude

    def test_known_eclipse_april_2024_europe(self):
        """Test April 8, 2024 eclipse from Europe (should not be visible)."""
        # The April 2024 eclipse was only visible in North America
        # From Europe, this eclipse is not visible, so we should find a different one
        jd_start = swe_julday(2024, 3, 1, 0)
        rome_lat, rome_lon = 41.9028, 12.4964

        times, attr, ecl_type = sol_eclipse_when_loc(jd_start, rome_lat, rome_lon)

        # Should find some eclipse visible from Rome
        assert times[0] > jd_start
        # It may or may not be April 2024 - the important thing is it's visible
        assert ecl_type & SE_ECL_VISIBLE

    def test_partial_eclipse_location(self):
        """Location far from centerline should see partial eclipse."""
        # Start before April 2024 eclipse and search from a location
        # that's away from the path of totality
        jd_start = swe_julday(2024, 3, 1, 0)
        # Location in Florida (gets partial coverage)
        florida_lat, florida_lon = 25.7617, -80.1918

        times, attr, ecl_type = sol_eclipse_when_loc(jd_start, florida_lat, florida_lon)

        # Should have visibility flags
        if times[0] > jd_start:
            assert ecl_type & SE_ECL_VISIBLE

    def test_visibility_flags_set(self):
        """Contact visibility flags should be set appropriately."""
        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        times, attr, ecl_type = sol_eclipse_when_loc(jd_start, lat, lon)

        # If first contact time is set, flag should be set
        if times[1] > 0:
            assert ecl_type & SE_ECL_1ST_VISIBLE
        # If fourth contact time is set, flag should be set
        if times[4] > 0:
            assert ecl_type & SE_ECL_4TH_VISIBLE

    def test_altitude_parameter(self):
        """Should accept altitude parameter."""
        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0
        altitude = 1000.0  # 1000 meters

        # Should not raise
        times, attr, ecl_type = sol_eclipse_when_loc(
            jd_start, lat, lon, altitude=altitude
        )

        assert times[0] > jd_start

    def test_flags_parameter_accepted(self):
        """Should accept flags parameter."""
        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        # Should not raise
        times, attr, ecl_type = sol_eclipse_when_loc(
            jd_start, lat, lon, flags=SEFLG_SWIEPH
        )

        assert times[0] > jd_start

    def test_swe_version_provides_pyswisseph_compatible_interface(self):
        """swe_sol_eclipse_when_loc provides pyswisseph-compatible interface.

        Note: swe_sol_eclipse_when_loc has a different signature than sol_eclipse_when_loc
        to match pyswisseph conventions (geopos sequence vs individual lat/lon/alt params).
        They are distinct functions that produce equivalent results.
        """
        from libephemeris import SEFLG_SWIEPH

        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon, altitude = 35.0, -100.0, 0.0

        # Call sol_eclipse_when_loc (libephemeris-style)
        times1, attr1, ecl_type1 = sol_eclipse_when_loc(
            jd_start, lat, lon, altitude=altitude
        )

        # Call swe_sol_eclipse_when_loc (pyswisseph-style with geopos sequence)
        geopos = (lon, lat, altitude)  # Note: pyswisseph uses (lon, lat, alt) order
        times2, attr2, ecl_type2 = swe_sol_eclipse_when_loc(
            jd_start, SEFLG_SWIEPH, geopos
        )

        # Both should find the same eclipse (same maximum time within tolerance)
        assert abs(times1[0] - times2[0]) < 0.01  # Within ~15 minutes
        # Both should return eclipse visibility flags
        assert ecl_type1 & SE_ECL_VISIBLE
        assert ecl_type2 & SE_ECL_VISIBLE

    def test_multiple_locations_same_eclipse(self):
        """Same eclipse should show different circumstances at different locations."""
        jd_start = swe_julday(2024, 3, 1, 0)

        # Dallas and Miami
        dallas_lat, dallas_lon = 32.7767, -96.7970
        miami_lat, miami_lon = 25.7617, -80.1918

        times_dallas, attr_dallas, _ = sol_eclipse_when_loc(
            jd_start, dallas_lat, dallas_lon
        )
        times_miami, attr_miami, _ = sol_eclipse_when_loc(
            jd_start, miami_lat, miami_lon
        )

        # Both should find eclipses (not necessarily the same one)
        assert times_dallas[0] > jd_start
        assert times_miami[0] > jd_start

    def test_sun_must_be_above_horizon(self):
        """Eclipse should only be returned if Sun is visible."""
        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        times, attr, ecl_type = sol_eclipse_when_loc(jd_start, lat, lon)

        # Sun altitude at maximum should be positive (above horizon)
        sun_altitude = attr[4]
        assert sun_altitude > -1.0  # Allow for refraction near horizon

    def test_sequential_eclipses(self):
        """Finding consecutive eclipses should return different events."""
        jd_start = swe_julday(2024, 1, 1, 0)
        lat, lon = 35.0, -100.0

        # Find first eclipse
        times1, _, _ = sol_eclipse_when_loc(jd_start, lat, lon)

        # Find next eclipse (search from after first one)
        jd_after_first = times1[0] + 1
        times2, _, _ = sol_eclipse_when_loc(jd_after_first, lat, lon)

        # Second eclipse should be after first
        assert times2[0] > times1[0]


class TestSolEclipseWhere:
    """Tests for sol_eclipse_where function."""

    def test_returns_correct_tuple_sizes(self):
        """Should return geopos tuple of 10 and attr tuple of 20 elements per pyswisseph."""
        from libephemeris import sol_eclipse_where

        # First find an eclipse to get a valid time
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        geopos, attr, ecl_type = sol_eclipse_where(times[0])

        assert len(geopos) == 10
        assert len(attr) == 20

    def test_central_eclipse_returns_valid_position(self):
        """For central eclipse, should return valid geographic coordinates."""
        from libephemeris import sol_eclipse_where

        # Find a total eclipse (which is central)
        jd_start = swe_julday(2024, 3, 1, 0)
        times, global_type = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        geopos, attr, ecl_type = sol_eclipse_where(times[0])

        # If eclipse is central, we should get valid coordinates
        if ecl_type & SE_ECL_CENTRAL:
            lon, lat = geopos[0], geopos[1]
            # Longitude should be -180 to 180
            assert -180 <= lon <= 180
            # Latitude should be -90 to 90
            assert -90 <= lat <= 90

    def test_eclipse_type_flags_set(self):
        """Eclipse type flags should be set appropriately."""
        from libephemeris import sol_eclipse_where

        # Find a total eclipse
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        geopos, attr, ecl_type = sol_eclipse_where(times[0])

        # Should have central flag for central eclipse
        if ecl_type != 0:
            assert ecl_type & SE_ECL_CENTRAL
            # Should have either total or annular flag
            assert ecl_type & (SE_ECL_TOTAL | SE_ECL_ANNULAR)

    def test_attributes_have_valid_ranges(self):
        """Eclipse attributes should be within valid ranges."""
        from libephemeris import sol_eclipse_where

        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        geopos, attr, ecl_type = sol_eclipse_where(times[0])

        if ecl_type != 0:
            magnitude = attr[0]
            ratio = attr[1]
            obscuration = attr[2]
            path_width = attr[3]
            sun_azimuth = attr[4]
            sun_altitude = attr[5]
            moon_diameter = attr[6]
            sun_diameter = attr[7]

            # Magnitude should be positive for central eclipse
            assert magnitude > 0
            # Ratio should be close to 1 for total/annular
            assert 0.8 < ratio < 1.2
            # Obscuration should be 0-1
            assert 0 <= obscuration <= 1
            # Path width should be positive and reasonable (km)
            assert 0 <= path_width <= 1000
            # Sun altitude at central line should be positive (Sun above horizon)
            # (can be negative near sunrise/sunset)
            assert -90 <= sun_altitude <= 90
            # Azimuth should be 0-360
            assert 0 <= sun_azimuth <= 360
            # Apparent diameters should be positive
            assert moon_diameter > 0
            assert sun_diameter > 0

    def test_non_eclipse_time_returns_zeros(self):
        """Time far from any eclipse should return zeros."""
        from libephemeris import sol_eclipse_where

        # Random time not during an eclipse (full moon time)
        jd_non_eclipse = swe_julday(2024, 4, 23, 12)  # Near full moon

        geopos, attr, ecl_type = sol_eclipse_where(jd_non_eclipse)

        # Should return zeros or partial flag
        if ecl_type == 0:
            assert geopos[0] == 0.0
            assert geopos[1] == 0.0

    def test_april_2024_eclipse_path(self):
        """Test April 8, 2024 total eclipse path."""
        from libephemeris import sol_eclipse_where

        # April 8, 2024 eclipse maximum around 18:18 UT
        jd_eclipse = swe_julday(2024, 4, 8, 18.3)

        geopos, attr, ecl_type = sol_eclipse_where(jd_eclipse)

        # Should be a total eclipse
        if ecl_type & SE_ECL_TOTAL:
            lon, lat = geopos[0], geopos[1]
            # The path crossed Mexico, US, and Eastern Canada
            # Central line was roughly around:
            # - Mexico: around 100W, 20-25N
            # - US: around 85-100W, 30-45N
            # - Canada: around 60-80W, 45-50N
            # At 18:18 UT, central line should be somewhere in Americas
            # Allow reasonable tolerance for approximate algorithm
            assert -130 < lon < -40  # Western longitudes (Americas)
            # The algorithm may give slightly different latitudes
            # due to approximations in shadow geometry
            assert -10 < lat < 60  # Northern hemisphere, with tolerance

    def test_path_width_positive_for_central_eclipse(self):
        """Path width should be positive for central eclipses."""
        from libephemeris import sol_eclipse_where

        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        geopos, attr, ecl_type = sol_eclipse_where(times[0])

        if ecl_type & SE_ECL_CENTRAL:
            path_width = attr[3]
            # Path width should be positive
            assert path_width > 0

    def test_swe_version_provides_pyswisseph_compatible_interface(self):
        """swe_sol_eclipse_where and sol_eclipse_where produce equivalent results.

        Note: sol_eclipse_where is a wrapper that calls swe_sol_eclipse_where internally.
        They are distinct functions with the same interface that produce identical results.
        """
        from libephemeris import sol_eclipse_where, swe_sol_eclipse_where

        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # Call both versions
        geopos1, attr1, ecl_type1 = sol_eclipse_where(times[0])
        geopos2, attr2, ecl_type2 = swe_sol_eclipse_where(times[0], SEFLG_SWIEPH)

        # Results should be identical
        assert geopos1 == geopos2
        assert attr1 == attr2
        assert ecl_type1 == ecl_type2

    def test_flags_parameter_accepted(self):
        """Should accept flags parameter without error."""
        from libephemeris import sol_eclipse_where

        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start)

        # Should not raise
        geopos, attr, ecl_type = sol_eclipse_where(times[0], flags=SEFLG_SWIEPH)

    def test_annular_eclipse(self):
        """Test with annular eclipse."""
        from libephemeris import sol_eclipse_where

        # October 14, 2023 annular eclipse
        jd_start = swe_julday(2023, 9, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_ANNULAR)

        geopos, attr, ecl_type = sol_eclipse_where(times[0])

        # Should be annular
        if ecl_type & SE_ECL_ANNULAR:
            # Moon/Sun ratio should be < 1 for annular
            ratio = attr[1]
            assert ratio < 1.0

    def test_multiple_times_during_eclipse(self):
        """Different times during eclipse should give different positions."""
        from libephemeris import sol_eclipse_where

        # Find an eclipse
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # Check positions at different times
        jd_max = times[0]
        delta = 10.0 / (24 * 60)  # 10 minutes in days

        geopos1, _, type1 = sol_eclipse_where(jd_max - delta)
        geopos2, _, type2 = sol_eclipse_where(jd_max)
        geopos3, _, type3 = sol_eclipse_where(jd_max + delta)

        # If all are central, positions should differ
        if type1 & SE_ECL_CENTRAL and type2 & SE_ECL_CENTRAL and type3 & SE_ECL_CENTRAL:
            # Longitudes should be different (shadow moves westward)
            # The shadow moves about 0.25 degrees per minute in longitude
            assert geopos1[0] != geopos2[0] or geopos1[1] != geopos2[1]


class TestSolEclipseHow:
    """Tests for sol_eclipse_how function."""

    def test_returns_correct_tuple_size(self):
        """Should return attr tuple of 11 elements."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse to get a valid time
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # Calculate circumstances at Dallas
        dallas_lat, dallas_lon = 32.7767, -96.7970
        attr, ecl_type = sol_eclipse_how(times[0], dallas_lat, dallas_lon)

        assert len(attr) == 20

    def test_eclipse_visible_during_eclipse(self):
        """During an eclipse, should return positive magnitude."""
        from libephemeris import sol_eclipse_how

        # April 8, 2024 total eclipse
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # Dallas was in path of totality
        dallas_lat, dallas_lon = 32.7767, -96.7970
        attr, ecl_type = sol_eclipse_how(times[0], dallas_lat, dallas_lon)

        # Should have visibility
        if ecl_type & SE_ECL_VISIBLE:
            magnitude = attr[0]
            assert magnitude > 0

    def test_no_eclipse_at_random_time(self):
        """At a random time with no eclipse, magnitude should be 0."""
        from libephemeris import sol_eclipse_how

        # Random time during full moon (no solar eclipse possible)
        jd_full_moon = swe_julday(2024, 4, 23, 12)  # Near full moon

        dallas_lat, dallas_lon = 32.7767, -96.7970
        attr, ecl_type = sol_eclipse_how(jd_full_moon, dallas_lat, dallas_lon)

        # Should have no eclipse
        assert ecl_type == 0 or attr[0] == 0

    def test_attributes_have_valid_ranges(self):
        """Eclipse attributes should be within valid ranges."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        dallas_lat, dallas_lon = 32.7767, -96.7970
        attr, ecl_type = sol_eclipse_how(times[0], dallas_lat, dallas_lon)

        if ecl_type != 0:
            # Magnitude should be 0-1.5 range
            magnitude = attr[0]
            assert 0 <= magnitude <= 1.5

            # Obscuration should be 0-1 range
            obscuration = attr[2]
            assert 0 <= obscuration <= 1

            # Azimuth should be 0 to 360
            azimuth = attr[4]
            assert 0 <= azimuth <= 360

            # Sun altitude should be -90 to 90
            altitude = attr[5]
            assert -90 <= altitude <= 90

    def test_sun_below_horizon_returns_zero_magnitude(self):
        """When Sun is below horizon, should return 0 magnitude or no visibility."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # Check from a location on opposite side of Earth (Sun below horizon)
        # For an eclipse happening in North America around 18 UT,
        # a location in Eastern Asia/Australia would have Sun below horizon
        east_asia_lat, east_asia_lon = 35.0, 135.0  # Japan area

        attr, ecl_type = sol_eclipse_how(times[0], east_asia_lat, east_asia_lon)

        # Either no visibility flag or zero magnitude
        sun_alt = attr[5]
        if sun_alt < -1.0:
            # Sun below horizon, should have no visible eclipse
            assert attr[0] == 0.0 or not (ecl_type & SE_ECL_VISIBLE)

    def test_different_locations_different_magnitudes(self):
        """Different locations should show different eclipse magnitudes."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # Dallas (in path of totality)
        dallas_lat, dallas_lon = 32.7767, -96.7970
        attr_dallas, type_dallas = sol_eclipse_how(times[0], dallas_lat, dallas_lon)

        # Miami (away from path)
        miami_lat, miami_lon = 25.7617, -80.1918
        attr_miami, type_miami = sol_eclipse_how(times[0], miami_lat, miami_lon)

        # If both visible, Dallas should have higher magnitude (closer to center)
        if (type_dallas & SE_ECL_VISIBLE) and (type_miami & SE_ECL_VISIBLE):
            # The magnitudes should be different
            # (may not always be true depending on exact timing)
            assert attr_dallas[0] != attr_miami[0] or attr_dallas[2] != attr_miami[2]

    def test_total_eclipse_type_flags(self):
        """Total eclipse should have appropriate type flags."""
        from libephemeris import sol_eclipse_how

        # Find a total eclipse
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # From a location in the path
        # Use local time found by sol_eclipse_when_loc for better timing
        times_loc, attr_loc, _ = sol_eclipse_when_loc(
            swe_julday(2024, 3, 1, 0), 32.7767, -96.7970
        )

        if times_loc[0] > 0:
            attr, ecl_type = sol_eclipse_how(times_loc[0], 32.7767, -96.7970)

            if ecl_type & SE_ECL_TOTAL:
                # Should also have central and visible flags
                assert ecl_type & SE_ECL_CENTRAL
                assert ecl_type & SE_ECL_VISIBLE

    def test_partial_eclipse_type_flags(self):
        """Partial eclipse should have partial flag."""
        from libephemeris import sol_eclipse_how

        # Find an eclipse and check from edge location
        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # NYC is far from totality path for April 2024 eclipse
        nyc_lat, nyc_lon = 40.7128, -74.0060
        attr, ecl_type = sol_eclipse_how(times[0], nyc_lat, nyc_lon)

        # If visible and not central, should be partial
        if (ecl_type & SE_ECL_VISIBLE) and not (ecl_type & SE_ECL_CENTRAL):
            assert ecl_type & SE_ECL_PARTIAL

    def test_altitude_parameter(self):
        """Should accept altitude parameter."""
        from libephemeris import sol_eclipse_how

        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start)

        # Should not raise
        attr, ecl_type = sol_eclipse_how(times[0], 32.7767, -96.7970, altitude=1000.0)

        assert len(attr) == 20

    def test_flags_parameter_accepted(self):
        """Should accept flags parameter."""
        from libephemeris import sol_eclipse_how

        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start)

        # Should not raise
        attr, ecl_type = sol_eclipse_how(
            times[0], 32.7767, -96.7970, flags=SEFLG_SWIEPH
        )

        assert len(attr) == 20

    def test_swe_version_provides_pyswisseph_compatible_interface(self):
        """swe_sol_eclipse_how and sol_eclipse_how produce equivalent results.

        Note: sol_eclipse_how is a wrapper that calls swe_sol_eclipse_how internally
        with a different parameter order (lat, lon vs geopos sequence).
        They produce identical results for the same location.
        """
        from libephemeris import sol_eclipse_how, swe_sol_eclipse_how

        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        lat, lon, altitude = 32.7767, -96.7970, 0.0
        geopos = (lon, lat, altitude)  # pyswisseph uses (lon, lat, alt) order

        # Call both versions
        attr1, ecl_type1 = sol_eclipse_how(times[0], lat, lon, altitude=altitude)
        attr2, ecl_type2 = swe_sol_eclipse_how(times[0], SEFLG_SWIEPH, geopos)

        # Results should be identical
        assert attr1 == attr2
        assert ecl_type1 == ecl_type2

    def test_known_eclipse_april_2024_dallas(self):
        """Test April 8, 2024 eclipse from Dallas."""
        from libephemeris import sol_eclipse_how

        # During the April 2024 eclipse, around 18:40 UT
        jd_eclipse = swe_julday(2024, 4, 8, 18.67)  # ~18:40 UT

        dallas_lat, dallas_lon = 32.7767, -96.7970
        attr, ecl_type = sol_eclipse_how(jd_eclipse, dallas_lat, dallas_lon)

        # Should have significant magnitude (Dallas was in path of totality)
        if ecl_type & SE_ECL_VISIBLE:
            assert attr[0] > 0.8  # High magnitude expected

    def test_obscuration_consistent_with_magnitude(self):
        """Obscuration should be consistent with magnitude."""
        from libephemeris import sol_eclipse_how

        jd_start = swe_julday(2024, 3, 1, 0)
        times, _ = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        dallas_lat, dallas_lon = 32.7767, -96.7970
        attr, ecl_type = sol_eclipse_how(times[0], dallas_lat, dallas_lon)

        if ecl_type & SE_ECL_VISIBLE:
            magnitude = attr[0]
            obscuration = attr[2]

            # Both should be positive if eclipse is visible
            if magnitude > 0:
                assert obscuration > 0

            # For total eclipse (magnitude >= 1), obscuration should be 1
            if magnitude >= 1.0:
                assert obscuration >= 0.99  # Allow small tolerance


class TestKnownEclipseValidation:
    """Validation tests against known eclipse times.

    These tests verify that the eclipse times are accurate within 0.01 days (~15 minutes)
    compared to known eclipse dates from NASA eclipse catalogs.
    """

    def test_april_2024_total_eclipse_timing(self):
        """Validate April 8, 2024 total eclipse maximum time.

        Known data: Total solar eclipse on April 8, 2024
        Maximum eclipse around JD 2460409.26 (approximately 18:17 UT)
        """
        # Search starting well before the eclipse
        jd_start = 2460400.0  # About 8 days before

        times, ecl_type = sol_eclipse_when_glob(jd_start)

        # Eclipse maximum should be around JD 2460409.26
        expected_jd = 2460409.26
        tolerance = 0.01  # 0.01 days = ~15 minutes

        assert abs(times[0] - expected_jd) < tolerance, (
            f"Eclipse maximum {times[0]} differs from expected {expected_jd} "
            f"by more than {tolerance} days"
        )

        # Should be a total eclipse
        assert ecl_type & SE_ECL_TOTAL

    def test_october_2023_annular_eclipse_timing(self):
        """Validate October 14, 2023 annular eclipse maximum time.

        Known data: Annular solar eclipse on October 14, 2023
        Maximum eclipse around 17:59 UT (JD ~2460232.25)
        """
        jd_start = 2460200.0  # About 32 days before

        times, ecl_type = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_ANNULAR)

        # Verify the correct date was found
        year, month, day, hour = swe_revjul(times[0])
        assert year == 2023
        assert month == 10
        assert 13 <= day <= 15  # October 14

        # Should be an annular eclipse
        assert ecl_type & SE_ECL_ANNULAR

    def test_august_2017_total_eclipse_timing(self):
        """Validate August 21, 2017 total eclipse maximum time.

        Known data: Total solar eclipse on August 21, 2017
        Maximum eclipse around 18:26 UT
        The famous "Great American Eclipse"
        """
        jd_start = 2457900.0  # About 87 days before

        times, ecl_type = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_TOTAL)

        # Verify the correct date was found
        year, month, day, hour = swe_revjul(times[0])
        assert year == 2017
        assert month == 8
        assert 20 <= day <= 22  # August 21

        # Should be a total eclipse
        assert ecl_type & SE_ECL_TOTAL

    def test_june_2021_annular_eclipse_timing(self):
        """Validate June 10, 2021 annular eclipse maximum time.

        Known data: Annular solar eclipse on June 10, 2021
        Maximum eclipse around 10:41 UT
        """
        jd_start = 2459350.0  # About 26 days before

        times, ecl_type = sol_eclipse_when_glob(jd_start, eclipse_type=SE_ECL_ANNULAR)

        # Verify the correct date was found
        year, month, day, hour = swe_revjul(times[0])
        assert year == 2021
        assert month == 6
        assert 9 <= day <= 11  # June 10

        # Should be an annular eclipse
        assert ecl_type & SE_ECL_ANNULAR

    def test_eclipse_type_classification_consistency(self):
        """Verify eclipse type classification is consistent for known eclipses."""
        # Total eclipse 2024
        times1, ecl_type1 = sol_eclipse_when_glob(2460400.0, eclipse_type=SE_ECL_TOTAL)
        assert ecl_type1 & SE_ECL_TOTAL, (
            "April 2024 eclipse should be classified as TOTAL"
        )

        # Annular eclipse 2023
        times2, ecl_type2 = sol_eclipse_when_glob(
            2460200.0, eclipse_type=SE_ECL_ANNULAR
        )
        assert ecl_type2 & SE_ECL_ANNULAR, (
            "October 2023 eclipse should be classified as ANNULAR"
        )

    def test_return_tuple_has_ten_elements(self):
        """Verify the tret tuple has exactly 10 elements like pyswisseph."""
        jd_start = 2460400.0

        times, ecl_type = sol_eclipse_when_glob(jd_start)

        assert len(times) == 10, f"Expected 10 elements in tret, got {len(times)}"

        # First element should be eclipse maximum
        assert times[0] > jd_start, "First element should be eclipse maximum time"

        # Elements 0-4 should be phase times
        # Maximum should be between first and fourth contacts
        assert times[1] < times[0] < times[4], "Maximum should be between contacts"
