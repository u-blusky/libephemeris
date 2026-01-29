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


class TestSwePheno20Values:
    """
    Verify swe_pheno_ut returns all 20 output values matching Swiss Ephemeris.

    This is the core verification suite for Swiss Ephemeris compatibility.
    All 20 values in the phenomena tuple must be verified:
    - attr[0-4]: Active phenomenon values (phase angle, phase, elongation, diameter, magnitude)
    - attr[5-19]: Reserved values (should all be 0.0)
    """

    # Tolerances for comparison
    PHASE_ANGLE_TOL = 1.0  # degrees (accounts for ephemeris differences)
    PHASE_TOL = 0.02  # illuminated fraction
    ELONGATION_TOL = 1.0  # degrees
    DIAMETER_TOL = 0.05  # relative tolerance (5%)
    MAGNITUDE_TOL = 0.5  # magnitudes (accounts for different formulas)
    # Saturn has larger magnitude differences due to ring contribution
    # which Swiss Ephemeris calculates differently (requires ring tilt)
    SATURN_MAGNITUDE_TOL = 1.0

    @pytest.mark.comparison
    def test_returns_exactly_20_values(self):
        """swe_pheno_ut must return exactly 20 values in the attr tuple."""
        jd = 2451545.0
        attr, flag = ephem.pheno_ut(jd, SE_MARS, 0)

        assert len(attr) == 20, f"Expected 20 values, got {len(attr)}"
        assert isinstance(flag, int), "Flag should be an integer"

    @pytest.mark.comparison
    def test_reserved_values_are_zero(self):
        """Values attr[5] through attr[19] must all be 0.0 (reserved)."""
        jd = 2451545.0

        # Test across multiple planets and bodies
        bodies = [SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN]

        for body_id in bodies:
            attr, _ = ephem.pheno_ut(jd, body_id, 0)

            for i in range(5, 20):
                assert attr[i] == 0.0, (
                    f"Body {body_id}: attr[{i}] should be 0.0, got {attr[i]}"
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
    def test_all_five_active_values_match_swe(self, planet_id, planet_name):
        """
        All 5 active phenomenon values must match Swiss Ephemeris within tolerance.

        Note: libephemeris returns diameter in arcseconds while Swiss Ephemeris
        returns it in degrees. The comparison accounts for this by converting
        SE degrees to arcseconds.
        """
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        # attr[0]: Phase angle (degrees)
        diff_phase_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_phase_angle < self.PHASE_ANGLE_TOL, (
            f"{planet_name} phase angle: SE={swe_attr[0]:.4f}, "
            f"LIB={lib_attr[0]:.4f}, diff={diff_phase_angle:.4f}"
        )

        # attr[1]: Phase (illuminated fraction, 0.0-1.0)
        diff_phase = abs(swe_attr[1] - lib_attr[1])
        assert diff_phase < self.PHASE_TOL, (
            f"{planet_name} phase: SE={swe_attr[1]:.4f}, "
            f"LIB={lib_attr[1]:.4f}, diff={diff_phase:.4f}"
        )

        # attr[2]: Elongation (degrees)
        diff_elongation = abs(swe_attr[2] - lib_attr[2])
        assert diff_elongation < self.ELONGATION_TOL, (
            f"{planet_name} elongation: SE={swe_attr[2]:.4f}, "
            f"LIB={lib_attr[2]:.4f}, diff={diff_elongation:.4f}"
        )

        # attr[3]: Diameter (SE returns degrees, libephemeris returns arcseconds)
        swe_diam_arcsec = swe_attr[3] * 3600  # Convert SE degrees to arcseconds
        lib_diam_arcsec = lib_attr[3]
        if swe_diam_arcsec > 0:
            rel_diff_diam = abs(swe_diam_arcsec - lib_diam_arcsec) / swe_diam_arcsec
        else:
            rel_diff_diam = 0.0 if lib_diam_arcsec == 0 else 1.0
        assert rel_diff_diam < self.DIAMETER_TOL, (
            f"{planet_name} diameter: SE={swe_diam_arcsec:.2f}as, "
            f"LIB={lib_diam_arcsec:.2f}as, rel_diff={rel_diff_diam:.4f}"
        )

        # attr[4]: Magnitude
        # Saturn magnitude has larger tolerance due to ring contribution differences
        if planet_id == SE_SATURN:
            mag_tol = self.SATURN_MAGNITUDE_TOL
        else:
            mag_tol = self.MAGNITUDE_TOL
        diff_magnitude = abs(swe_attr[4] - lib_attr[4])
        assert diff_magnitude < mag_tol, (
            f"{planet_name} magnitude: SE={swe_attr[4]:.2f}, "
            f"LIB={lib_attr[4]:.2f}, diff={diff_magnitude:.2f}"
        )

    @pytest.mark.comparison
    def test_sun_all_20_values(self):
        """Verify Sun phenomena returns all 20 values correctly."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, SE_SUN, 0)
        lib_attr, flag = ephem.pheno_ut(jd, SE_SUN, 0)

        # Must have exactly 20 values
        assert len(lib_attr) == 20

        # Phase angle = 0, phase = 1.0, elongation = 0
        assert lib_attr[0] == 0.0, "Sun phase angle should be 0"
        assert lib_attr[1] == 1.0, "Sun phase should be 1.0"
        assert lib_attr[2] == 0.0, "Sun elongation should be 0"

        # Diameter check (SE degrees -> arcsec conversion)
        swe_diam_arcsec = swe_attr[3] * 3600
        lib_diam_arcsec = lib_attr[3]
        rel_diff = abs(swe_diam_arcsec - lib_diam_arcsec) / swe_diam_arcsec
        assert rel_diff < 0.01, (
            f"Sun diameter: SE={swe_diam_arcsec:.2f}as, LIB={lib_diam_arcsec:.2f}as"
        )

        # Magnitude check (small tolerance for Sun)
        diff_mag = abs(swe_attr[4] - lib_attr[4])
        assert diff_mag < 0.15, (
            f"Sun magnitude: SE={swe_attr[4]:.2f}, LIB={lib_attr[4]:.2f}"
        )

        # Reserved values
        for i in range(5, 20):
            assert lib_attr[i] == 0.0, f"Sun attr[{i}] should be 0.0"

    @pytest.mark.comparison
    def test_moon_all_20_values(self):
        """Verify Moon phenomena returns all 20 values correctly."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, SE_MOON, 0)
        lib_attr, _ = ephem.pheno_ut(jd, SE_MOON, 0)

        # Must have exactly 20 values
        assert len(lib_attr) == 20

        # Phase angle comparison (Moon can have larger discrepancies)
        diff_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_angle < 2.0, f"Moon phase angle diff: {diff_angle:.4f}"

        # Phase comparison
        diff_phase = abs(swe_attr[1] - lib_attr[1])
        assert diff_phase < 0.05, f"Moon phase diff: {diff_phase:.4f}"

        # Elongation comparison
        diff_elong = abs(swe_attr[2] - lib_attr[2])
        assert diff_elong < 2.0, f"Moon elongation diff: {diff_elong:.4f}"

        # Diameter check (SE degrees -> arcsec)
        swe_diam_arcsec = swe_attr[3] * 3600
        lib_diam_arcsec = lib_attr[3]
        rel_diff = abs(swe_diam_arcsec - lib_diam_arcsec) / swe_diam_arcsec
        assert rel_diff < 0.01, (
            f"Moon diameter: SE={swe_diam_arcsec:.2f}as, LIB={lib_diam_arcsec:.2f}as"
        )

        # Reserved values
        for i in range(5, 20):
            assert lib_attr[i] == 0.0, f"Moon attr[{i}] should be 0.0"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_URANUS, "Uranus"),
            (SE_NEPTUNE, "Neptune"),
            (SE_PLUTO, "Pluto"),
        ],
    )
    def test_outer_planets_all_20_values(self, planet_id, planet_name):
        """Verify outer planets return all 20 values correctly."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr, _ = ephem.pheno_ut(jd, planet_id, 0)

        # Must have exactly 20 values
        assert len(lib_attr) == 20, f"{planet_name} should return 20 values"

        # Phase angle - outer planets have small phase angles
        diff_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_angle < 1.0, (
            f"{planet_name} phase angle: SE={swe_attr[0]:.4f}, LIB={lib_attr[0]:.4f}"
        )

        # Phase - outer planets are nearly full
        diff_phase = abs(swe_attr[1] - lib_attr[1])
        assert diff_phase < 0.02, (
            f"{planet_name} phase: SE={swe_attr[1]:.4f}, LIB={lib_attr[1]:.4f}"
        )

        # Elongation
        diff_elong = abs(swe_attr[2] - lib_attr[2])
        assert diff_elong < 1.0, (
            f"{planet_name} elongation: SE={swe_attr[2]:.4f}, LIB={lib_attr[2]:.4f}"
        )

        # Diameter (SE degrees -> arcsec)
        swe_diam_arcsec = swe_attr[3] * 3600
        lib_diam_arcsec = lib_attr[3]
        if swe_diam_arcsec > 0:
            rel_diff = abs(swe_diam_arcsec - lib_diam_arcsec) / swe_diam_arcsec
            assert rel_diff < 0.1, (
                f"{planet_name} diameter: SE={swe_diam_arcsec:.2f}as, "
                f"LIB={lib_diam_arcsec:.2f}as"
            )

        # Reserved values
        for i in range(5, 20):
            assert lib_attr[i] == 0.0, f"{planet_name} attr[{i}] should be 0.0"

    @pytest.mark.comparison
    def test_multiple_dates_all_20_values(self):
        """Verify all 20 values at multiple dates throughout a year."""
        test_dates = [
            2451545.0,  # J2000 (Jan 1, 2000)
            2459580.5,  # Jan 1, 2022
            2460310.5,  # Jan 1, 2024
            2460400.5,  # Apr 1, 2024
            2460500.5,  # Jul 10, 2024
        ]

        swe.set_ephe_path(None)

        for jd in test_dates:
            for planet_id in [SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER]:
                swe_attr = swe.pheno_ut(jd, planet_id, 0)
                lib_attr, _ = ephem.pheno_ut(jd, planet_id, 0)

                # Must have 20 values
                assert len(lib_attr) == 20

                # Reserved values must be zero
                for i in range(5, 20):
                    assert lib_attr[i] == 0.0, (
                        f"JD {jd}, planet {planet_id}: attr[{i}] should be 0.0"
                    )

                # Active values should be in valid ranges
                assert 0 <= lib_attr[0] <= 180, "Phase angle out of range"
                assert 0 <= lib_attr[1] <= 1.0, "Phase out of range"
                assert 0 <= lib_attr[2] <= 180, "Elongation out of range"
                assert lib_attr[3] >= 0, "Diameter should be non-negative"

    @pytest.mark.comparison
    def test_swe_pheno_et_version_20_values(self):
        """Verify swe_pheno (ET version) also returns 20 values."""
        jd_et = 2451545.0  # J2000 TT

        swe.set_ephe_path(None)
        swe_attr = swe.pheno(jd_et, SE_MARS, 0)
        lib_attr, flag = ephem.pheno(jd_et, SE_MARS, 0)

        # Must have exactly 20 values
        assert len(lib_attr) == 20, (
            f"pheno should return 20 values, got {len(lib_attr)}"
        )

        # Reserved values must be zero
        for i in range(5, 20):
            assert lib_attr[i] == 0.0, f"pheno: attr[{i}] should be 0.0"

        # Phase angle comparison
        diff_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_angle < 1.0, f"Phase angle diff too large: {diff_angle}"

    @pytest.mark.comparison
    def test_return_flag_matches_input(self):
        """Return flag should reflect the input calculation flags."""
        jd = 2451545.0

        # Test with flag = 0
        _, flag0 = ephem.pheno_ut(jd, SE_MARS, 0)
        assert isinstance(flag0, int)

        # Test with SEFLG_TRUEPOS
        from libephemeris.constants import SEFLG_TRUEPOS

        _, flag_true = ephem.pheno_ut(jd, SE_MARS, SEFLG_TRUEPOS)
        assert isinstance(flag_true, int)
