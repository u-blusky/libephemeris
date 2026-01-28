"""
Tests for the seorbel.txt parser.

This module tests the parse_seorbel() function and related utilities
for parsing Swiss Ephemeris seorbel.txt orbital elements files.
"""

import math
import tempfile
from pathlib import Path

import pytest

from libephemeris.hypothetical import (
    parse_seorbel,
    SeorbelElements,
    TPolynomial,
    get_seorbel_body_by_name,
    calc_seorbel_position,
    get_bundled_seorbel_path,
    load_bundled_seorbel,
    _parse_t_polynomial,
    _parse_epoch_or_equinox,
)


class TestTPolynomial:
    """Tests for the TPolynomial dataclass."""

    def test_constant_only(self):
        """Test polynomial with only constant term."""
        poly = TPolynomial(constant=42.5, linear=0.0)
        assert poly.evaluate(0.0) == 42.5
        assert poly.evaluate(1.0) == 42.5
        assert poly.evaluate(-1.0) == 42.5

    def test_with_linear_term(self):
        """Test polynomial with linear term."""
        poly = TPolynomial(constant=100.0, linear=50.0)
        assert poly.evaluate(0.0) == 100.0
        assert poly.evaluate(1.0) == 150.0
        assert poly.evaluate(2.0) == 200.0
        assert poly.evaluate(-1.0) == 50.0

    def test_negative_linear(self):
        """Test polynomial with negative linear term."""
        poly = TPolynomial(constant=100.0, linear=-10.0)
        assert poly.evaluate(0.0) == 100.0
        assert poly.evaluate(1.0) == 90.0
        assert poly.evaluate(5.0) == 50.0


class TestParseEpochOrEquinox:
    """Tests for _parse_epoch_or_equinox helper."""

    def test_j1900(self):
        """Test parsing J1900 epoch."""
        jd, is_jdate = _parse_epoch_or_equinox("J1900")
        assert jd == 2415020.0
        assert is_jdate is False

    def test_j2000(self):
        """Test parsing J2000 epoch."""
        jd, is_jdate = _parse_epoch_or_equinox("J2000")
        assert jd == 2451545.0
        assert is_jdate is False

    def test_b1950(self):
        """Test parsing B1950 epoch."""
        jd, is_jdate = _parse_epoch_or_equinox("B1950")
        assert jd == pytest.approx(2433282.42345905, abs=0.0001)
        assert is_jdate is False

    def test_jdate(self):
        """Test parsing JDATE equinox."""
        jd, is_jdate = _parse_epoch_or_equinox("JDATE")
        assert jd is None
        assert is_jdate is True

    def test_numeric_jd(self):
        """Test parsing numeric Julian Day."""
        jd, is_jdate = _parse_epoch_or_equinox("2451545.0")
        assert jd == 2451545.0
        assert is_jdate is False

    def test_case_insensitive(self):
        """Test that parsing is case-insensitive."""
        jd1, _ = _parse_epoch_or_equinox("j2000")
        jd2, _ = _parse_epoch_or_equinox("J2000")
        assert jd1 == jd2


class TestParseTPolynomial:
    """Tests for _parse_t_polynomial helper."""

    def test_simple_number(self):
        """Test parsing a simple number."""
        poly = _parse_t_polynomial("163.7409")
        assert poly.constant == pytest.approx(163.7409)
        assert poly.linear == 0.0

    def test_negative_number(self):
        """Test parsing a negative number."""
        poly = _parse_t_polynomial("-45.5")
        assert poly.constant == pytest.approx(-45.5)
        assert poly.linear == 0.0

    def test_polynomial_plus(self):
        """Test parsing polynomial with + operator."""
        poly = _parse_t_polynomial("252.8987988 + 707550.7341 * T")
        assert poly.constant == pytest.approx(252.8987988)
        assert poly.linear == pytest.approx(707550.7341)

    def test_polynomial_minus(self):
        """Test parsing polynomial with - operator."""
        poly = _parse_t_polynomial("47.787931 - 1670.056 * T")
        assert poly.constant == pytest.approx(47.787931)
        assert poly.linear == pytest.approx(-1670.056)

    def test_polynomial_compact(self):
        """Test parsing polynomial without spaces."""
        poly = _parse_t_polynomial("322.212069+1670.056*T")
        assert poly.constant == pytest.approx(322.212069)
        assert poly.linear == pytest.approx(1670.056)

    def test_polynomial_compact_minus(self):
        """Test parsing compact polynomial with minus."""
        poly = _parse_t_polynomial("47.787931-1670.056*T")
        assert poly.constant == pytest.approx(47.787931)
        assert poly.linear == pytest.approx(-1670.056)

    def test_zero_constant(self):
        """Test parsing with zero constant."""
        poly = _parse_t_polynomial("0.0")
        assert poly.constant == 0.0
        assert poly.linear == 0.0


class TestParseSeorbel:
    """Tests for parse_seorbel function."""

    def test_parse_real_seorbel_file(self):
        """Test parsing the actual Swiss Ephemeris seorbel.txt file."""
        seorbel_path = (
            Path(__file__).parent.parent / "swisseph" / "ephe" / "seorbel.txt"
        )
        if not seorbel_path.exists():
            pytest.skip("seorbel.txt not found in swisseph/ephe/")

        elements = parse_seorbel(seorbel_path)

        # Should have parsed multiple elements
        assert len(elements) > 0

        # Check that Cupido is first (it's line #1 in the file)
        cupido = get_seorbel_body_by_name(elements, "Cupido")
        assert cupido is not None
        assert cupido.semi_axis == pytest.approx(40.99837)
        assert cupido.eccentricity.constant == pytest.approx(0.00460)
        assert cupido.epoch_jd == 2415020.0  # J1900

    def test_parse_simple_file(self):
        """Test parsing a simple custom seorbel file."""
        content = """\
# Custom test file
J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, TestPlanet
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_seorbel(filepath)
            assert len(elements) == 1

            elem = elements[0]
            assert elem.name == "TestPlanet"
            assert elem.epoch_jd == 2451545.0  # J2000
            assert elem.semi_axis == 100.0
            assert elem.eccentricity.constant == 0.1
            assert elem.arg_perihelion.constant == 45.0
            assert elem.asc_node.constant == 30.0
            assert elem.inclination.constant == 5.0
            assert elem.mean_anomaly.constant == 0.0
        finally:
            Path(filepath).unlink()

    def test_parse_geocentric_body(self):
        """Test parsing a geocentric body (with ', geo' suffix)."""
        content = """\
J2000, JDATE, 248.8833, 0.0029833, 0.0, 0.0, 0.0, 0.0, Waldemath, geo
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_seorbel(filepath)
            assert len(elements) == 1

            elem = elements[0]
            assert elem.name == "Waldemath"
            assert elem.is_geocentric is True
            assert elem.equinox_is_jdate is True
        finally:
            Path(filepath).unlink()

    def test_parse_t_polynomial_elements(self):
        """Test parsing elements with T-polynomial expressions."""
        # From the actual seorbel.txt Vulcan entry
        content = """\
J1900,JDATE, 252.8987988 + 707550.7341 * T, 0.13744, 0.019, 322.212069+1670.056*T, 47.787931-1670.056*T, 7.5, Vulcan
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_seorbel(filepath)
            assert len(elements) == 1

            elem = elements[0]
            assert elem.name == "Vulcan"
            assert elem.epoch_jd == 2415020.0  # J1900
            assert elem.equinox_is_jdate is True

            # Check T-polynomial parsing
            assert elem.mean_anomaly.constant == pytest.approx(252.8987988)
            assert elem.mean_anomaly.linear == pytest.approx(707550.7341)

            assert elem.arg_perihelion.constant == pytest.approx(322.212069)
            assert elem.arg_perihelion.linear == pytest.approx(1670.056)

            assert elem.asc_node.constant == pytest.approx(47.787931)
            assert elem.asc_node.linear == pytest.approx(-1670.056)

            assert elem.semi_axis == pytest.approx(0.13744)
            assert elem.eccentricity.constant == pytest.approx(0.019)
            assert elem.inclination.constant == pytest.approx(7.5)
        finally:
            Path(filepath).unlink()

    def test_parse_inline_comments(self):
        """Test that inline comments are handled correctly."""
        content = """\
J1900, J1900, 163.7409, 40.99837, 0.00460, 171.4333, 129.8325, 1.0833, Cupido   # 1
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_seorbel(filepath)
            assert len(elements) == 1
            assert elements[0].name == "Cupido"
        finally:
            Path(filepath).unlink()

    def test_skip_comment_lines(self):
        """Test that comment-only lines are skipped."""
        content = """\
# This is a comment
    # Indented comment
J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, Planet1
# Another comment
J2000, J2000, 90.0, 60.0, 0.0, 0.0, 0.0, 0.0, Planet2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_seorbel(filepath)
            assert len(elements) == 2
            assert elements[0].name == "Planet1"
            assert elements[1].name == "Planet2"
        finally:
            Path(filepath).unlink()

    def test_skip_empty_lines(self):
        """Test that empty lines are skipped."""
        content = """\

J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, Planet1

J2000, J2000, 90.0, 60.0, 0.0, 0.0, 0.0, 0.0, Planet2

"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_seorbel(filepath)
            assert len(elements) == 2
        finally:
            Path(filepath).unlink()

    def test_file_not_found(self):
        """Test that FileNotFoundError is raised for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_seorbel("/nonexistent/path/seorbel.txt")

    def test_parse_numeric_epoch(self):
        """Test parsing with numeric Julian Day epoch."""
        content = """\
2368547.66, 2431456.5, 0.0, 77.775, 0.3, 0.7, 0, 0, Isis-Transpluto
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_seorbel(filepath)
            assert len(elements) == 1

            elem = elements[0]
            assert elem.name == "Isis-Transpluto"
            assert elem.epoch_jd == pytest.approx(2368547.66)
            assert elem.equinox_jd == pytest.approx(2431456.5)
            assert elem.semi_axis == pytest.approx(77.775)
            assert elem.eccentricity.constant == pytest.approx(0.3)
        finally:
            Path(filepath).unlink()

    def test_preserves_order(self):
        """Test that elements are returned in file order."""
        content = """\
J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, First
J2000, J2000, 0.0, 60.0, 0.0, 0.0, 0.0, 0.0, Second
J2000, J2000, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, Third
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_seorbel(filepath)
            assert len(elements) == 3
            assert elements[0].name == "First"
            assert elements[1].name == "Second"
            assert elements[2].name == "Third"
        finally:
            Path(filepath).unlink()


class TestGetSeorbelBodyByName:
    """Tests for get_seorbel_body_by_name function."""

    def test_find_existing_body(self):
        """Test finding an existing body by name."""
        elements = [
            SeorbelElements(
                name="TestBody",
                epoch_jd=2451545.0,
                equinox_jd=2451545.0,
                equinox_is_jdate=False,
                mean_anomaly=TPolynomial(0.0, 0.0),
                semi_axis=50.0,
                eccentricity=TPolynomial(0.1, 0.0),
                arg_perihelion=TPolynomial(0.0, 0.0),
                asc_node=TPolynomial(0.0, 0.0),
                inclination=TPolynomial(0.0, 0.0),
            )
        ]
        result = get_seorbel_body_by_name(elements, "TestBody")
        assert result is not None
        assert result.name == "TestBody"

    def test_case_insensitive_search(self):
        """Test that name search is case-insensitive."""
        elements = [
            SeorbelElements(
                name="TestBody",
                epoch_jd=2451545.0,
                equinox_jd=2451545.0,
                equinox_is_jdate=False,
                mean_anomaly=TPolynomial(0.0, 0.0),
                semi_axis=50.0,
                eccentricity=TPolynomial(0.1, 0.0),
                arg_perihelion=TPolynomial(0.0, 0.0),
                asc_node=TPolynomial(0.0, 0.0),
                inclination=TPolynomial(0.0, 0.0),
            )
        ]
        result = get_seorbel_body_by_name(elements, "testbody")
        assert result is not None
        assert result.name == "TestBody"

    def test_not_found_returns_none(self):
        """Test that searching for non-existent body returns None."""
        elements = [
            SeorbelElements(
                name="TestBody",
                epoch_jd=2451545.0,
                equinox_jd=2451545.0,
                equinox_is_jdate=False,
                mean_anomaly=TPolynomial(0.0, 0.0),
                semi_axis=50.0,
                eccentricity=TPolynomial(0.1, 0.0),
                arg_perihelion=TPolynomial(0.0, 0.0),
                asc_node=TPolynomial(0.0, 0.0),
                inclination=TPolynomial(0.0, 0.0),
            )
        ]
        result = get_seorbel_body_by_name(elements, "NonExistent")
        assert result is None


class TestCalcSeorbelPosition:
    """Tests for calc_seorbel_position function."""

    def test_circular_orbit(self):
        """Test position calculation for circular orbit."""
        elem = SeorbelElements(
            name="CircularPlanet",
            epoch_jd=2451545.0,  # J2000.0
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            mean_anomaly=TPolynomial(0.0, 0.0),  # At 0 deg at epoch
            semi_axis=1.0,  # 1 AU
            eccentricity=TPolynomial(0.0, 0.0),  # Circular
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(0.0, 0.0),
        )

        pos = calc_seorbel_position(elem, 2451545.0)  # At epoch
        lon, lat, dist, dlon, dlat, ddist = pos

        # At epoch with M=0, omega=0, Omega=0, longitude should be ~0
        assert lon == pytest.approx(0.0, abs=0.1)
        assert lat == pytest.approx(0.0, abs=0.01)
        assert dist == pytest.approx(1.0, abs=0.001)

        # Mean motion should be ~360/(365.25) deg/day = ~0.986 deg/day
        assert dlon == pytest.approx(360.0 / 365.25, rel=0.01)
        assert dlat == pytest.approx(0.0, abs=0.001)
        assert ddist == pytest.approx(0.0, abs=0.001)

    def test_elliptic_orbit(self):
        """Test position calculation for elliptic orbit."""
        elem = SeorbelElements(
            name="EllipticPlanet",
            epoch_jd=2451545.0,  # J2000.0
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            mean_anomaly=TPolynomial(0.0, 0.0),  # At perihelion
            semi_axis=2.0,  # 2 AU
            eccentricity=TPolynomial(0.5, 0.0),  # High eccentricity
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(0.0, 0.0),
        )

        pos = calc_seorbel_position(elem, 2451545.0)  # At epoch
        lon, lat, dist, dlon, dlat, ddist = pos

        # At perihelion, distance = a(1-e) = 2(1-0.5) = 1 AU
        assert dist == pytest.approx(1.0, abs=0.001)

    def test_inclined_orbit(self):
        """Test position calculation for inclined orbit."""
        elem = SeorbelElements(
            name="InclinedPlanet",
            epoch_jd=2451545.0,  # J2000.0
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            mean_anomaly=TPolynomial(90.0, 0.0),  # 90 deg from perihelion
            semi_axis=1.0,
            eccentricity=TPolynomial(0.0, 0.0),  # Circular
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(30.0, 0.0),  # 30 deg inclination
        )

        pos = calc_seorbel_position(elem, 2451545.0)  # At epoch
        lon, lat, dist, dlon, dlat, ddist = pos

        # With i=30 deg and body at 90 deg from node, latitude should be significant
        # At u = 90 deg (ascending node + 90), latitude = arcsin(sin(i) * sin(u))
        # = arcsin(sin(30) * sin(90)) = arcsin(0.5 * 1) = 30 deg
        assert lat == pytest.approx(30.0, abs=1.0)

    def test_time_dependent_elements(self):
        """Test position calculation with time-dependent elements."""
        # Create an element with a high mean motion polynomial
        elem = SeorbelElements(
            name="FastMover",
            epoch_jd=2451545.0,  # J2000.0
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            # 100 + 36000*T means 36000 deg/century motion = ~0.986 deg/day secular
            mean_anomaly=TPolynomial(100.0, 36000.0),
            semi_axis=1.0,
            eccentricity=TPolynomial(0.0, 0.0),
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(0.0, 0.0),
        )

        # At epoch
        pos1 = calc_seorbel_position(elem, 2451545.0)
        # 1 century later
        pos2 = calc_seorbel_position(elem, 2451545.0 + 36525.0)

        # The longitude should have advanced by 36000 deg (modulo 360)
        # 36000 % 360 = 0, so longitudes should be the same
        expected_advance = 36000.0 % 360.0
        actual_advance = (pos2[0] - pos1[0]) % 360.0
        assert actual_advance == pytest.approx(expected_advance, abs=1.0)

    def test_position_matches_cupido(self):
        """Test that calculated position is reasonable for Cupido."""
        seorbel_path = (
            Path(__file__).parent.parent / "swisseph" / "ephe" / "seorbel.txt"
        )
        if not seorbel_path.exists():
            pytest.skip("seorbel.txt not found in swisseph/ephe/")

        elements = parse_seorbel(seorbel_path)
        cupido = get_seorbel_body_by_name(elements, "Cupido")
        assert cupido is not None

        # Calculate position at J2000.0
        pos = calc_seorbel_position(cupido, 2451545.0)
        lon, lat, dist, dlon, dlat, ddist = pos

        # Check reasonable values
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0
        assert dist == pytest.approx(cupido.semi_axis, rel=0.1)

        # Mean motion should be very slow (Cupido has ~262 year period)
        # n = 360 / (a^1.5 * 365.25) deg/day
        expected_n = 360.0 / (cupido.semi_axis**1.5 * 365.25)
        assert dlon == pytest.approx(expected_n, rel=0.1)


class TestSeorbelElementsDataclass:
    """Tests for the SeorbelElements dataclass."""

    def test_get_mean_motion(self):
        """Test mean motion calculation from semi-major axis."""
        elem = SeorbelElements(
            name="Test",
            epoch_jd=2451545.0,
            equinox_jd=2451545.0,
            equinox_is_jdate=False,
            mean_anomaly=TPolynomial(0.0, 0.0),
            semi_axis=1.0,  # 1 AU = 1 year period
            eccentricity=TPolynomial(0.0, 0.0),
            arg_perihelion=TPolynomial(0.0, 0.0),
            asc_node=TPolynomial(0.0, 0.0),
            inclination=TPolynomial(0.0, 0.0),
        )

        n = elem.get_mean_motion()
        # For a = 1 AU, period = 1 year, n = 360/365.25 deg/day
        expected = 360.0 / 365.25
        assert n == pytest.approx(expected, rel=0.001)

    def test_line_number_tracking(self):
        """Test that line numbers are tracked correctly."""
        content = """\
# Comment
J2000, J2000, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, Planet1
# Another comment
J2000, J2000, 0.0, 60.0, 0.0, 0.0, 0.0, 0.0, Planet2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(content)
            f.flush()
            filepath = f.name

        try:
            elements = parse_seorbel(filepath)
            assert elements[0].line_number == 2  # After the first comment
            assert elements[1].line_number == 4  # After another comment
        finally:
            Path(filepath).unlink()


class TestBundledSeorbel:
    """Tests for bundled seorbel.txt functionality."""

    def test_get_bundled_seorbel_path_exists(self):
        """Test that the bundled seorbel.txt path exists."""
        path = get_bundled_seorbel_path()
        assert path.exists(), "Bundled seorbel.txt should exist"
        assert path.is_file(), "Bundled seorbel.txt should be a file"
        assert path.name == "seorbel.txt", "File should be named seorbel.txt"

    def test_get_bundled_seorbel_path_in_package(self):
        """Test that the bundled file is inside the libephemeris package."""
        path = get_bundled_seorbel_path()
        # The file should be in the libephemeris package directory
        assert "libephemeris" in str(path), "Path should contain 'libephemeris'"

    def test_load_bundled_seorbel_returns_list(self):
        """Test that load_bundled_seorbel returns a list of elements."""
        elements = load_bundled_seorbel()
        assert isinstance(elements, list), "Should return a list"
        assert len(elements) > 0, "Should have at least one element"
        for elem in elements:
            assert isinstance(elem, SeorbelElements), (
                "Each element should be SeorbelElements"
            )

    def test_load_bundled_seorbel_contains_expected_bodies(self):
        """Test that load_bundled_seorbel contains expected hypothetical bodies."""
        elements = load_bundled_seorbel()
        names = [elem.name for elem in elements]

        # Check for Uranian planets (Hamburg School)
        expected_bodies = [
            "Cupido",
            "Hades",
            "Zeus",
            "Kronos",
            "Apollon",
            "Admetos",
            # Note: The file has "Vulcanus" not "Vulkanus"
            "Poseidon",
        ]

        for body in expected_bodies:
            assert any(body.lower() in name.lower() for name in names), (
                f"Expected body '{body}' not found in bundled seorbel.txt"
            )

    def test_load_bundled_seorbel_cupido_elements(self):
        """Test that Cupido has correct orbital elements from bundled file."""
        elements = load_bundled_seorbel()
        cupido = get_seorbel_body_by_name(elements, "Cupido")

        assert cupido is not None, "Cupido should be found"
        assert abs(cupido.semi_axis - 40.99837) < 0.001, "Cupido semi-axis should match"
        assert abs(cupido.eccentricity.constant - 0.00460) < 0.0001, (
            "Cupido eccentricity should match"
        )

    def test_load_bundled_seorbel_transpluto_elements(self):
        """Test that Transpluto/Isis has correct orbital elements."""
        elements = load_bundled_seorbel()
        transpluto = get_seorbel_body_by_name(elements, "Isis-Transpluto")

        assert transpluto is not None, "Transpluto should be found"
        assert abs(transpluto.semi_axis - 77.775) < 0.001, (
            "Transpluto semi-axis should match"
        )
        assert abs(transpluto.eccentricity.constant - 0.3) < 0.01, (
            "Transpluto eccentricity should match"
        )

    def test_load_bundled_seorbel_geocentric_bodies(self):
        """Test that geocentric bodies are correctly identified."""
        elements = load_bundled_seorbel()

        # Waldemath is geocentric (orbits Earth, not Sun)
        waldemath = get_seorbel_body_by_name(elements, "Waldemath")
        assert waldemath is not None, "Waldemath should be found"
        assert waldemath.is_geocentric, "Waldemath should be geocentric"

        # White Moon (Selena) is geocentric
        selena = get_seorbel_body_by_name(elements, "Selena/White Moon")
        assert selena is not None, "Selena/White Moon should be found"
        assert selena.is_geocentric, "Selena should be geocentric"

    def test_load_bundled_seorbel_can_calculate_position(self):
        """Test that we can calculate positions from bundled elements."""
        elements = load_bundled_seorbel()
        cupido = get_seorbel_body_by_name(elements, "Cupido")

        assert cupido is not None, "Cupido should be found"

        # Calculate position at J2000.0
        pos = calc_seorbel_position(cupido, 2451545.0)

        # Basic sanity checks
        assert 0 <= pos[0] < 360, "Longitude should be in range [0, 360)"
        assert -90 <= pos[1] <= 90, "Latitude should be in range [-90, 90]"
        assert pos[2] > 0, "Distance should be positive"

    def test_bundled_seorbel_matches_swisseph_source(self):
        """Test that bundled seorbel.txt matches the original Swiss Ephemeris source."""
        bundled_path = get_bundled_seorbel_path()
        bundled_elements = parse_seorbel(bundled_path)

        # If the original file exists, compare
        original_path = (
            Path(__file__).parent.parent / "swisseph" / "ephe" / "seorbel.txt"
        )
        if original_path.exists():
            original_elements = parse_seorbel(original_path)

            # Should have same number of bodies
            assert len(bundled_elements) == len(original_elements), (
                "Bundled file should have same number of bodies as original"
            )

            # Check that first few bodies match
            for i in range(min(5, len(bundled_elements))):
                assert bundled_elements[i].name == original_elements[i].name, (
                    f"Body {i} name should match"
                )
                assert (
                    abs(bundled_elements[i].semi_axis - original_elements[i].semi_axis)
                    < 0.0001
                ), f"Body {i} semi-axis should match"
        else:
            # Skip comparison if original file doesn't exist
            pytest.skip("Original swisseph/ephe/seorbel.txt not found for comparison")
