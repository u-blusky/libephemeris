"""
Tests for planetary phenomena calculations (swe_pheno, swe_pheno_ut).

Tests cover:
- Phase angle calculations
- Phase (illumination fraction)
- Elongation from Sun
- Apparent diameter
- Visual magnitude
- Comparison with pyswisseph
"""

import pytest
import math
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
)


class TestPhenoBasic:
    """Test basic pheno function behavior."""

    @pytest.mark.unit
    def test_pheno_ut_returns_tuple(self):
        """pheno_ut should return a tuple with attr and flag."""
        jd = 2451545.0  # J2000
        result = ephem.pheno_ut(jd, SE_MARS, 0)

        assert isinstance(result, tuple)
        assert len(result) == 2
        attr, flag = result
        assert isinstance(attr, tuple)
        assert len(attr) >= 5  # At least 5 phenomenon values

    @pytest.mark.unit
    def test_pheno_returns_tuple(self):
        """pheno should return a tuple with attr and flag."""
        jd = 2451545.0  # J2000 TT
        result = ephem.pheno(jd, SE_MARS, 0)

        assert isinstance(result, tuple)
        assert len(result) == 2
        attr, flag = result
        assert isinstance(attr, tuple)
        assert len(attr) >= 5

    @pytest.mark.unit
    def test_swe_prefix_aliases_exist(self):
        """Both swe_pheno and swe_pheno_ut should be accessible."""
        assert hasattr(ephem, "swe_pheno")
        assert hasattr(ephem, "swe_pheno_ut")
        assert hasattr(ephem, "pheno")
        assert hasattr(ephem, "pheno_ut")


class TestSunPhenomena:
    """Test Sun phenomena (edge case)."""

    @pytest.mark.unit
    def test_sun_phase_is_full(self):
        """Sun should always have phase = 1.0 (fully illuminated)."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, SE_SUN, 0)

        phase_angle = attr[0]
        phase = attr[1]
        elongation = attr[2]

        assert phase_angle == 0.0, "Sun phase angle should be 0"
        assert phase == 1.0, "Sun phase should be 1.0"
        assert elongation == 0.0, "Sun elongation from itself should be 0"

    @pytest.mark.unit
    def test_sun_diameter_reasonable(self):
        """Sun apparent diameter should be around 32 arcminutes."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, SE_SUN, 0)

        diameter = attr[3]
        # Sun diameter ~32 arcmin = ~1920 arcsec
        assert 1850 < diameter < 1980, f"Sun diameter {diameter} arcsec unexpected"

    @pytest.mark.unit
    def test_sun_magnitude_reasonable(self):
        """Sun magnitude should be around -26.7."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, SE_SUN, 0)

        magnitude = attr[4]
        assert -27.0 < magnitude < -26.0, f"Sun magnitude {magnitude} unexpected"


class TestMoonPhenomena:
    """Test Moon phenomena calculations."""

    @pytest.mark.unit
    def test_moon_phase_varies(self):
        """Moon phase should vary between 0 and 1."""
        # Test at different times to see variation
        jd_new = 2451550.1  # Near new moon (around Jan 6, 2000)
        jd_full = 2451565.0  # Near full moon (around Jan 21, 2000)

        attr_new, _ = ephem.pheno_ut(jd_new, SE_MOON, 0)
        attr_full, _ = ephem.pheno_ut(jd_full, SE_MOON, 0)

        phase_new = attr_new[1]
        phase_full = attr_full[1]

        # Phases should be in valid range
        assert 0 <= phase_new <= 1.0
        assert 0 <= phase_full <= 1.0

    @pytest.mark.unit
    def test_moon_diameter_reasonable(self):
        """Moon apparent diameter should be around 30-34 arcminutes."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, SE_MOON, 0)

        diameter = attr[3]
        # Moon diameter ~31 arcmin = ~1860 arcsec, varies with distance
        assert 1700 < diameter < 2100, f"Moon diameter {diameter} arcsec unexpected"

    @pytest.mark.unit
    def test_moon_elongation_varies(self):
        """Moon elongation should be between 0 and 180 degrees."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, SE_MOON, 0)

        elongation = attr[2]
        assert 0 <= elongation <= 180, f"Moon elongation {elongation} out of range"


class TestPlanetPhenomena:
    """Test planetary phenomena calculations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
        ],
    )
    def test_planet_phase_in_range(self, planet_id, planet_name):
        """Planet phase should be between 0 and 1."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        phase = attr[1]
        assert 0 <= phase <= 1.0, f"{planet_name} phase {phase} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
        ],
    )
    def test_planet_elongation_in_range(self, planet_id, planet_name):
        """Planet elongation should be between 0 and 180 degrees."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        elongation = attr[2]
        assert 0 <= elongation <= 180, (
            f"{planet_name} elongation {elongation} out of range"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
        ],
    )
    def test_planet_diameter_positive(self, planet_id, planet_name):
        """Planet apparent diameter should be positive."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        diameter = attr[3]
        assert diameter > 0, f"{planet_name} diameter should be positive"

    @pytest.mark.unit
    def test_jupiter_diameter_large(self):
        """Jupiter should have the largest apparent diameter of planets."""
        jd = 2451545.0

        diameters = {}
        for planet_id in [SE_MARS, SE_JUPITER, SE_SATURN]:
            attr, _ = ephem.pheno_ut(jd, planet_id, 0)
            diameters[planet_id] = attr[3]

        # Jupiter should have larger diameter than Mars
        assert diameters[SE_JUPITER] > diameters[SE_MARS]


class TestMagnitudes:
    """Test visual magnitude calculations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,min_mag,max_mag,planet_name",
        [
            (SE_VENUS, -5.0, -3.0, "Venus"),  # Venus is always bright
            (SE_JUPITER, -3.0, -1.0, "Jupiter"),  # Jupiter is bright
            (SE_SATURN, -0.5, 1.5, "Saturn"),  # Saturn varies more
            (SE_MARS, -3.0, 2.0, "Mars"),  # Mars varies significantly
        ],
    )
    def test_planet_magnitude_typical_range(
        self, planet_id, min_mag, max_mag, planet_name
    ):
        """Planet magnitudes should be in typical observable range."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        magnitude = attr[4]
        assert min_mag < magnitude < max_mag, (
            f"{planet_name} magnitude {magnitude:.2f} "
            f"outside typical range [{min_mag}, {max_mag}]"
        )


class TestComparisonWithSwissEph:
    """Compare pheno results with pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
        ],
    )
    def test_phase_angle_matches_swe(self, planet_id, planet_name):
        """Phase angle should match Swiss Ephemeris within tolerance."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        # pyswisseph returns just the tuple, not (tuple, flag)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr, lib_flag = ephem.pheno_ut(jd, planet_id, 0)

        # Phase angle comparison (allow 1 degree tolerance)
        swe_phase_angle = swe_attr[0]
        lib_phase_angle = lib_attr[0]

        diff = abs(swe_phase_angle - lib_phase_angle)
        assert diff < 1.0, (
            f"{planet_name} phase angle: SWE={swe_phase_angle:.4f}, "
            f"LIB={lib_phase_angle:.4f}, diff={diff:.4f}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
        ],
    )
    def test_phase_matches_swe(self, planet_id, planet_name):
        """Phase (illumination) should match Swiss Ephemeris within tolerance."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        swe_phase = swe_attr[1]
        lib_phase = lib_attr[1]

        diff = abs(swe_phase - lib_phase)
        assert diff < 0.02, (
            f"{planet_name} phase: SWE={swe_phase:.4f}, "
            f"LIB={lib_phase:.4f}, diff={diff:.4f}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
        ],
    )
    def test_elongation_matches_swe(self, planet_id, planet_name):
        """Elongation should match Swiss Ephemeris within tolerance."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        swe_elong = swe_attr[2]
        lib_elong = lib_attr[2]

        diff = abs(swe_elong - lib_elong)
        assert diff < 1.0, (
            f"{planet_name} elongation: SWE={swe_elong:.4f}, "
            f"LIB={lib_elong:.4f}, diff={diff:.4f}"
        )

    @pytest.mark.comparison
    def test_moon_phenomena_matches_swe(self):
        """Moon phenomena should match Swiss Ephemeris."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, SE_MOON, 0)
        lib_attr, _ = ephem.pheno_ut(jd, SE_MOON, 0)

        # Phase angle
        diff_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_angle < 2.0, f"Moon phase angle diff: {diff_angle:.4f}"

        # Phase
        diff_phase = abs(swe_attr[1] - lib_attr[1])
        assert diff_phase < 0.05, f"Moon phase diff: {diff_phase:.4f}"

        # Elongation
        diff_elong = abs(swe_attr[2] - lib_attr[2])
        assert diff_elong < 2.0, f"Moon elongation diff: {diff_elong:.4f}"


class TestOuterPlanets:
    """Test outer planet phenomena."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_URANUS, "Uranus"),
            (SE_NEPTUNE, "Neptune"),
            (SE_PLUTO, "Pluto"),
        ],
    )
    def test_outer_planet_nearly_full(self, planet_id, planet_name):
        """Outer planets should be nearly fully illuminated (phase ~1.0)."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        phase = attr[1]
        # Outer planets are almost fully illuminated from Earth
        assert phase > 0.95, f"{planet_name} phase {phase} should be > 0.95"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_URANUS, "Uranus"),
            (SE_NEPTUNE, "Neptune"),
            (SE_PLUTO, "Pluto"),
        ],
    )
    def test_outer_planet_small_phase_angle(self, planet_id, planet_name):
        """Outer planets should have small phase angles."""
        jd = 2451545.0
        attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        phase_angle = attr[0]
        # Outer planets have small phase angles (usually < 5 degrees)
        assert phase_angle < 10.0, (
            f"{planet_name} phase angle {phase_angle} should be < 10"
        )


class TestPhaseAngleFormula:
    """Test phase angle and phase relationship."""

    @pytest.mark.unit
    def test_phase_from_phase_angle(self):
        """Phase should be calculated correctly from phase angle."""
        jd = 2451545.0

        for planet_id in [SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER]:
            attr, _ = ephem.pheno_ut(jd, planet_id, 0)

            phase_angle = attr[0]
            phase = attr[1]

            # Phase = (1 + cos(phase_angle)) / 2
            expected_phase = (1 + math.cos(math.radians(phase_angle))) / 2

            diff = abs(phase - expected_phase)
            assert diff < 0.001, (
                f"Planet {planet_id}: Phase {phase:.6f} doesn't match "
                f"expected {expected_phase:.6f} from angle {phase_angle:.2f}"
            )


class TestEdgeCases:
    """Test edge cases and error handling."""

    @pytest.mark.unit
    def test_unsupported_body_returns_zeros(self):
        """Calculating phenomena for unsupported body should return zeros."""
        jd = 2451545.0
        # Use a body ID that's not in _PLANET_MAP (e.g., a fictitious planet)
        unsupported_id = 99  # Not a valid planet ID
        attr, _ = ephem.pheno_ut(jd, unsupported_id, 0)

        # Should return zeros for unsupported body
        assert all(v == 0.0 for v in attr[:5]), (
            "Unsupported body phenomena should return zeros"
        )

    @pytest.mark.unit
    def test_different_dates(self):
        """Phenomena should vary with date."""
        jd1 = 2451545.0  # J2000
        jd2 = 2451645.0  # 100 days later

        attr1, _ = ephem.pheno_ut(jd1, SE_MARS, 0)
        attr2, _ = ephem.pheno_ut(jd2, SE_MARS, 0)

        # At least one value should differ significantly
        differences = [abs(attr1[i] - attr2[i]) for i in range(5)]
        assert max(differences) > 0.01, "Phenomena should change over 100 days"
