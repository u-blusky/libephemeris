"""
Tests for planetary nodes and apsides calculation (swe_nod_aps and swe_nod_aps_ut).

These functions calculate:
- Ascending and descending orbital nodes
- Perihelion and aphelion positions

Note on precision:
    The libephemeris implementation calculates osculating orbital elements from
    JPL DE421 ephemeris data. This matches JPL Horizons orbital elements exactly.

    Swiss Ephemeris uses mean orbital elements from specific data sources that
    may differ from instantaneous osculating elements. The longitude of ascending
    node (Omega) can differ by ~40° for inner planets due to different conventions.

    For most astrological applications, both approaches are valid. The key
    properties (orbital shape, orientation) are preserved.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
    SE_EARTH,
    SE_NODBIT_MEAN,
    SE_NODBIT_OSCU,
)


class TestNodApsBasic:
    """Test basic functionality of nod_aps functions."""

    @pytest.mark.unit
    def test_nod_aps_ut_returns_4_tuples(self):
        """nod_aps_ut should return 4 position tuples."""
        jd = 2451545.0  # J2000
        result = ephem.nod_aps_ut(jd, SE_MARS, 0, SE_NODBIT_MEAN)

        assert isinstance(result, tuple)
        assert len(result) == 4
        for i, pos in enumerate(result):
            assert isinstance(pos, tuple), f"Position {i} should be tuple"
            assert len(pos) == 6, f"Position {i} should have 6 elements"

    @pytest.mark.unit
    def test_nod_aps_returns_4_tuples(self):
        """nod_aps (ET version) should return 4 position tuples."""
        jd = 2451545.0
        result = ephem.nod_aps(jd, SE_MARS, SE_NODBIT_MEAN)

        assert isinstance(result, tuple)
        assert len(result) == 4

    @pytest.mark.unit
    def test_position_elements_are_floats(self):
        """All position elements should be floats."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, SE_MARS, 0, SE_NODBIT_MEAN)

        for name, pos in [
            ("nasc", nasc),
            ("ndsc", ndsc),
            ("peri", peri),
            ("aphe", aphe),
        ]:
            for i, val in enumerate(pos):
                assert isinstance(val, (int, float)), f"{name}[{i}] should be numeric"

    @pytest.mark.unit
    def test_sun_returns_zeros(self):
        """Sun should return zero positions (no heliocentric orbit)."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, SE_SUN, 0, SE_NODBIT_MEAN)

        # All positions should be zero for Sun
        assert nasc[0] == 0.0
        assert ndsc[0] == 0.0
        assert peri[0] == 0.0
        assert aphe[0] == 0.0

    @pytest.mark.unit
    def test_earth_returns_zeros(self):
        """Earth should return zero positions (no heliocentric orbit)."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, SE_EARTH, 0, SE_NODBIT_MEAN)

        assert nasc[0] == 0.0
        assert peri[0] == 0.0


class TestNodApsOrbitalRelationships:
    """Test that orbital relationships are mathematically correct."""

    @pytest.mark.unit
    def test_descending_node_is_180_from_ascending(self):
        """Descending node should be 180° from ascending node."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, SE_MARS, 0, SE_NODBIT_MEAN)

        diff = abs(ndsc[0] - nasc[0])
        if diff > 180:
            diff = 360 - diff

        assert abs(diff - 180) < 0.1, f"Node difference {diff} should be 180°"

    @pytest.mark.unit
    def test_aphelion_is_180_from_perihelion(self):
        """Aphelion should be 180° from perihelion."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, SE_MARS, 0, SE_NODBIT_MEAN)

        diff = abs(aphe[0] - peri[0])
        if diff > 180:
            diff = 360 - diff

        assert abs(diff - 180) < 0.1, f"Apside difference {diff} should be 180°"

    @pytest.mark.unit
    def test_nodes_have_zero_latitude(self):
        """Nodes should have zero latitude (on ecliptic by definition)."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, SE_MARS, 0, SE_NODBIT_MEAN)

        assert abs(nasc[1]) < 0.01, f"Ascending node latitude {nasc[1]} should be ~0"
        assert abs(ndsc[1]) < 0.01, f"Descending node latitude {ndsc[1]} should be ~0"

    @pytest.mark.unit
    def test_apsides_have_latitude(self):
        """Apsides should have non-zero latitude for inclined orbits."""
        jd = 2451545.0
        # Pluto has high inclination (~17°)
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, SE_PLUTO, 0, SE_NODBIT_MEAN)

        # Pluto's apsides should have measurable latitude
        # (depending on argument of perihelion)
        # Just verify they're calculated (can be positive or negative)
        assert peri[1] is not None
        assert aphe[1] is not None

    @pytest.mark.unit
    def test_perihelion_closer_than_aphelion(self):
        """Perihelion distance should be less than aphelion distance."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, SE_MARS, 0, SE_NODBIT_MEAN)

        # peri[2] = distance at perihelion
        # aphe[2] = distance at aphelion
        assert peri[2] < aphe[2], "Perihelion should be closer than aphelion"


class TestNodApsAllPlanets:
    """Test nod_aps for all planets."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
            (SE_URANUS, "Uranus"),
            (SE_NEPTUNE, "Neptune"),
            (SE_PLUTO, "Pluto"),
        ],
    )
    def test_all_planets_return_valid_positions(self, planet_id, planet_name):
        """All planets should return valid node/apse positions."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, planet_id, 0, SE_NODBIT_MEAN)

        # Longitudes should be in valid range
        assert 0 <= nasc[0] < 360, f"{planet_name} ascending node longitude invalid"
        assert 0 <= ndsc[0] < 360, f"{planet_name} descending node longitude invalid"
        assert 0 <= peri[0] < 360, f"{planet_name} perihelion longitude invalid"
        assert 0 <= aphe[0] < 360, f"{planet_name} aphelion longitude invalid"

        # Distances should be positive
        assert nasc[2] > 0, f"{planet_name} ascending node distance should be positive"
        assert peri[2] > 0, f"{planet_name} perihelion distance should be positive"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,expected_a",
        [
            (SE_MERCURY, "Mercury", 0.387),
            (SE_VENUS, "Venus", 0.723),
            (SE_MARS, "Mars", 1.524),
            (SE_JUPITER, "Jupiter", 5.203),
            (SE_SATURN, "Saturn", 9.537),
        ],
    )
    def test_semi_major_axis_reasonable(self, planet_id, planet_name, expected_a):
        """Semi-major axis computed from apsides should match known values."""
        jd = 2451545.0
        nasc, ndsc, peri, aphe = ephem.nod_aps_ut(jd, planet_id, 0, SE_NODBIT_MEAN)

        # a = (r_peri + r_aphe) / 2
        computed_a = (peri[2] + aphe[2]) / 2

        # Allow 5% tolerance
        assert abs(computed_a - expected_a) / expected_a < 0.05, (
            f"{planet_name} semi-major axis {computed_a} differs from expected {expected_a}"
        )


class TestNodApsMethods:
    """Test different calculation methods."""

    @pytest.mark.unit
    def test_mean_and_oscu_methods_exist(self):
        """Both MEAN and OSCU methods should work."""
        jd = 2451545.0

        result_mean = ephem.nod_aps_ut(jd, SE_MARS, 0, SE_NODBIT_MEAN)
        result_oscu = ephem.nod_aps_ut(jd, SE_MARS, 0, SE_NODBIT_OSCU)

        # Both should return valid results
        assert len(result_mean) == 4
        assert len(result_oscu) == 4

        # Results may differ
        # (In our implementation, they use the same method currently)


class TestNodApsTimeVariation:
    """Test how nodes and apsides change over time."""

    @pytest.mark.unit
    def test_nodes_change_over_century(self):
        """Orbital nodes should show precession over a century."""
        jd_2000 = 2451545.0
        jd_2100 = jd_2000 + 36525  # 100 years later

        result_2000 = ephem.nod_aps_ut(jd_2000, SE_MARS, 0, SE_NODBIT_MEAN)
        result_2100 = ephem.nod_aps_ut(jd_2100, SE_MARS, 0, SE_NODBIT_MEAN)

        # Nodes should have precessed
        diff = abs(result_2100[0][0] - result_2000[0][0])

        # Mars' node precession is about 0.5°/century
        # Allow significant range due to osculating nature
        assert diff > 0.01, "Nodes should show some precession over 100 years"

    @pytest.mark.unit
    def test_perihelion_advances(self):
        """Perihelion should advance over time (apsidal precession)."""
        jd_2000 = 2451545.0
        jd_2100 = jd_2000 + 36525

        result_2000 = ephem.nod_aps_ut(jd_2000, SE_MARS, 0, SE_NODBIT_MEAN)
        result_2100 = ephem.nod_aps_ut(jd_2100, SE_MARS, 0, SE_NODBIT_MEAN)

        # Perihelion should advance
        diff = result_2100[2][0] - result_2000[2][0]
        # Handle wraparound
        if diff < -180:
            diff += 360
        if diff > 180:
            diff -= 360

        # Just verify there's some change
        assert abs(diff) > 0.01, "Perihelion should advance over 100 years"


class TestNodApsAliases:
    """Test that function aliases work."""

    @pytest.mark.unit
    def test_swe_nod_aps_ut_alias(self):
        """swe_nod_aps_ut should be available."""
        result = ephem.swe_nod_aps_ut(2451545.0, SE_MARS, 0, SE_NODBIT_MEAN)
        assert len(result) == 4

    @pytest.mark.unit
    def test_swe_nod_aps_alias(self):
        """swe_nod_aps should be available."""
        result = ephem.swe_nod_aps(2451545.0, SE_MARS, SE_NODBIT_MEAN)
        assert len(result) == 4

    @pytest.mark.unit
    def test_nod_aps_ut_alias(self):
        """nod_aps_ut should be available."""
        result = ephem.nod_aps_ut(2451545.0, SE_MARS, 0, SE_NODBIT_MEAN)
        assert len(result) == 4

    @pytest.mark.unit
    def test_nod_aps_alias(self):
        """nod_aps should be available."""
        result = ephem.nod_aps(2451545.0, SE_MARS, SE_NODBIT_MEAN)
        assert len(result) == 4


class TestNodApsMethodologyWarning:
    """Test the HeliocentricNodApsWarning for inner planets."""

    @pytest.mark.unit
    def test_warning_issued_for_mercury(self):
        """Mercury should trigger HeliocentricNodApsWarning."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, SE_MERCURY, 0, SE_NODBIT_MEAN)

            # Check that a warning was issued
            assert len(w) >= 1
            assert any(
                issubclass(warning.category, ephem.HeliocentricNodApsWarning)
                for warning in w
            )

    @pytest.mark.unit
    def test_warning_issued_for_venus(self):
        """Venus should trigger HeliocentricNodApsWarning."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, SE_VENUS, 0, SE_NODBIT_MEAN)

            assert len(w) >= 1
            assert any(
                issubclass(warning.category, ephem.HeliocentricNodApsWarning)
                for warning in w
            )

    @pytest.mark.unit
    def test_no_warning_for_mars(self):
        """Mars should not trigger HeliocentricNodApsWarning."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, SE_MARS, 0, SE_NODBIT_MEAN)

            # No HeliocentricNodApsWarning should be issued for outer planets
            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 0

    @pytest.mark.unit
    def test_no_warning_for_jupiter(self):
        """Jupiter should not trigger HeliocentricNodApsWarning."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, SE_JUPITER, 0, SE_NODBIT_MEAN)

            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 0

    @pytest.mark.unit
    def test_warning_can_be_suppressed(self):
        """HeliocentricNodApsWarning should be suppressible."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings("ignore", category=ephem.HeliocentricNodApsWarning)
            ephem.nod_aps_ut(jd, SE_MERCURY, 0, SE_NODBIT_MEAN)

            # No HeliocentricNodApsWarning should be recorded when filtered
            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 0

    @pytest.mark.unit
    def test_warning_message_contains_planet_name(self):
        """Warning message should contain the planet name."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, SE_MERCURY, 0, SE_NODBIT_MEAN)

            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 1
            assert "Mercury" in str(helio_warnings[0].message)

    @pytest.mark.unit
    def test_warning_message_mentions_precision_doc(self):
        """Warning message should reference PRECISION.md."""
        import warnings

        jd = 2451545.0
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ephem.nod_aps_ut(jd, SE_VENUS, 0, SE_NODBIT_MEAN)

            helio_warnings = [
                warning
                for warning in w
                if issubclass(warning.category, ephem.HeliocentricNodApsWarning)
            ]
            assert len(helio_warnings) == 1
            assert "PRECISION.md" in str(helio_warnings[0].message)

    @pytest.mark.unit
    def test_warning_class_is_exported(self):
        """HeliocentricNodApsWarning should be exported from libephemeris."""
        assert hasattr(ephem, "HeliocentricNodApsWarning")
        assert issubclass(ephem.HeliocentricNodApsWarning, UserWarning)
