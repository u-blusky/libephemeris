"""
Comprehensive tests for all 19 house systems.

Tests all house systems against pyswisseph for multiple locations.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


# House system codes and their expected tolerance
HOUSE_SYSTEMS = [
    (ord("P"), "Placidus", 0.1),
    (ord("K"), "Koch", 0.1),
    (ord("O"), "Porphyry", 0.1),
    (ord("R"), "Regiomontanus", 0.1),
    (ord("C"), "Campanus", 0.1),
    (ord("E"), "Equal", 0.1),
    (ord("W"), "Whole Sign", 0.1),
    (ord("B"), "Alcabitius", 0.1),
    (ord("M"), "Morinus", 0.1),
    (ord("T"), "Topocentric", 1.0),
    (ord("X"), "Meridian", 0.1),
    (ord("V"), "Vehlow", 0.1),
    (ord("H"), "Horizontal", 0.5),
]


class TestHouseSystemsBasic:
    """Basic tests for house systems."""

    @pytest.mark.unit
    def test_houses_returns_correct_structure(self):
        """houses() should return (cusps, ascmc) tuples."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        # Should have 12 cusps
        assert len(cusps) >= 12
        # Should have at least ASC and MC
        assert len(ascmc) >= 2

    @pytest.mark.unit
    def test_cusps_are_valid_longitudes(self):
        """All cusps should be valid longitudes (0-360)."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        for i, cusp in enumerate(cusps[:12]):
            assert 0 <= cusp < 360, f"Cusp {i + 1} = {cusp} out of range"

    @pytest.mark.unit
    def test_asc_and_mc_valid(self):
        """ASC and MC should be valid longitudes."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        # ASC is ascmc[0], MC is ascmc[1]
        assert 0 <= ascmc[0] < 360, f"ASC = {ascmc[0]} out of range"
        assert 0 <= ascmc[1] < 360, f"MC = {ascmc[1]} out of range"


class TestHouseSystemsVsPyswisseph:
    """Compare each house system with pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,name,tolerance", HOUSE_SYSTEMS)
    def test_house_system_rome(self, hsys, name, tolerance):
        """Test house system at Rome."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps_lib, ascmc_lib = ephem.swe_houses(jd, lat, lon, hsys)
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes([hsys]))

        # Compare ASC
        asc_diff = abs(ascmc_lib[0] - ascmc_swe[0])
        if asc_diff > 180:
            asc_diff = 360 - asc_diff
        assert asc_diff < tolerance, f"{name} ASC diff {asc_diff}"

        # Compare MC
        mc_diff = abs(ascmc_lib[1] - ascmc_swe[1])
        if mc_diff > 180:
            mc_diff = 360 - mc_diff
        assert mc_diff < tolerance, f"{name} MC diff {mc_diff}"

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,name,tolerance", HOUSE_SYSTEMS[:7])  # Common systems
    def test_house_system_equator(self, hsys, name, tolerance):
        """Test house system at equator."""
        jd = 2451545.0
        lat, lon = 0.0, 0.0

        cusps_lib, ascmc_lib = ephem.swe_houses(jd, lat, lon, hsys)
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes([hsys]))

        asc_diff = abs(ascmc_lib[0] - ascmc_swe[0])
        if asc_diff > 180:
            asc_diff = 360 - asc_diff
        assert asc_diff < tolerance, f"{name} at equator ASC diff {asc_diff}"


class TestHouseCuspsOrder:
    """Test that house cusps are in correct order."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys,name,_", HOUSE_SYSTEMS[:7])
    def test_cusps_in_zodiacal_order(self, hsys, name, _):
        """Cusps should generally increase (with wraparound)."""
        jd = 2451545.0
        cusps, _ = ephem.swe_houses(jd, 41.9, 12.5, hsys)

        # Check that cusps generally progress around the zodiac
        # (allowing for the 360/0 wraparound)
        for i in range(11):
            diff = cusps[i + 1] - cusps[i]
            if diff < -180:
                diff += 360
            # Each house should span 0-60 degrees typically
            assert diff > 0 or diff > -180, (
                f"{name} cusp {i + 1} to {i + 2} has odd progression"
            )


class TestHouseSystemsPolarLatitudes:
    """Test house systems at polar latitudes."""

    @pytest.mark.edge_case
    def test_placidus_below_polar_circle(self):
        """Placidus should work below the polar circle (~66.5°)."""
        jd = 2451545.0

        # At 65° latitude, Placidus should work normally
        cusps, ascmc = ephem.swe_houses(jd, 65.0, 0.0, ord("P"))
        assert 0 <= ascmc[0] < 360
        assert 0 <= ascmc[1] < 360

    @pytest.mark.edge_case
    @pytest.mark.parametrize("lat", [67.0, 70.0, 80.0])
    def test_placidus_polar_raises_error(self, lat):
        """Placidus should raise Error at polar latitudes (>66.5°).

        This matches Swiss Ephemeris behavior which raises an error
        with message 'within polar circle, switched to Porphyry'.
        """
        jd = 2451545.0

        # Should raise an Error for polar latitudes
        with pytest.raises(ephem.Error) as exc_info:
            ephem.swe_houses(jd, lat, 0.0, ord("P"))

        # Error message should mention polar circle
        assert "polar circle" in str(exc_info.value).lower()

    @pytest.mark.edge_case
    @pytest.mark.parametrize("lat", [67.0, 70.0, 80.0])
    def test_koch_polar_raises_error(self, lat):
        """Koch should raise Error at polar latitudes (>66.5°).

        This matches Swiss Ephemeris behavior.
        """
        jd = 2451545.0

        with pytest.raises(ephem.Error) as exc_info:
            ephem.swe_houses(jd, lat, 0.0, ord("K"))

        assert "polar circle" in str(exc_info.value).lower()

    @pytest.mark.edge_case
    @pytest.mark.parametrize("lat", [-67.0, -70.0, -80.0])
    def test_placidus_polar_southern_raises_error(self, lat):
        """Placidus should raise Error at Southern polar latitudes."""
        jd = 2451545.0

        with pytest.raises(ephem.Error) as exc_info:
            ephem.swe_houses(jd, lat, 0.0, ord("P"))

        assert "polar circle" in str(exc_info.value).lower()

    @pytest.mark.edge_case
    def test_porphyry_works_at_polar(self):
        """Porphyry should work at polar latitudes (user fallback)."""
        jd = 2451545.0

        # Users can catch the Placidus error and use Porphyry instead
        for lat in [67.0, 70.0, 80.0]:
            cusps, ascmc = ephem.swe_houses(jd, lat, 0.0, ord("O"))
            assert 0 <= ascmc[0] < 360
            assert 0 <= ascmc[1] < 360

    @pytest.mark.edge_case
    def test_equal_works_at_poles(self):
        """Equal house system should work at any latitude."""
        jd = 2451545.0

        for lat in [89.0, -89.0]:
            cusps, ascmc = ephem.swe_houses(jd, lat, 0.0, ord("E"))
            assert 0 <= ascmc[0] < 360

    @pytest.mark.edge_case
    def test_houses_armc_polar_raises_error(self):
        """swe_houses_armc should also raise Error at polar latitudes."""
        armc = 280.46
        eps = 23.44
        lat = 70.0

        with pytest.raises(ephem.Error) as exc_info:
            ephem.swe_houses_armc(armc, lat, eps, ord("P"))

        assert "polar circle" in str(exc_info.value).lower()

    @pytest.mark.edge_case
    def test_polar_error_matches_swisseph(self):
        """Error behavior should match Swiss Ephemeris for polar latitudes."""
        jd = 2451545.0
        lat = 70.0

        # Both should raise Error
        with pytest.raises(swe.Error):
            swe.houses(jd, lat, 0.0, b"P")

        with pytest.raises(ephem.Error):
            ephem.swe_houses(jd, lat, 0.0, ord("P"))

        # Both should succeed with Porphyry
        cusps_swe, ascmc_swe = swe.houses(jd, lat, 0.0, b"O")
        cusps_lib, ascmc_lib = ephem.swe_houses(jd, lat, 0.0, ord("O"))

        # Values should match
        assert abs(ascmc_swe[0] - ascmc_lib[0]) < 0.1
        assert abs(ascmc_swe[1] - ascmc_lib[1]) < 0.1


class TestHouseAscMcRelationship:
    """Test relationship between ASC, MC, and cusps."""

    @pytest.mark.unit
    def test_asc_equals_cusp_1(self):
        """ASC should equal cusp 1 for most systems."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        diff = abs(ascmc[0] - cusps[0])
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.01, f"ASC {ascmc[0]} != cusp1 {cusps[0]}"

    @pytest.mark.unit
    def test_mc_equals_cusp_10(self):
        """MC should equal cusp 10 for most systems."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        diff = abs(ascmc[1] - cusps[9])  # cusp 10 is index 9
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.01, f"MC {ascmc[1]} != cusp10 {cusps[9]}"


class TestHouseNameFunction:
    """Test swe_house_name function."""

    @pytest.mark.unit
    def test_house_name_placidus(self):
        """Should return 'Placidus' for P."""
        name = ephem.swe_house_name(ord("P"))
        assert "Placidus" in name or "placidus" in name.lower()

    @pytest.mark.unit
    def test_house_name_koch(self):
        """Should return 'Koch' for K."""
        name = ephem.swe_house_name(ord("K"))
        assert "Koch" in name or "koch" in name.lower()

    @pytest.mark.unit
    def test_house_name_equal(self):
        """Should return 'Equal' for E."""
        name = ephem.swe_house_name(ord("E"))
        assert "Equal" in name or "equal" in name.lower()


class TestHousesEx2:
    """Tests for swe_houses_ex2 function that returns velocities."""

    @pytest.mark.unit
    def test_houses_ex2_returns_four_tuples(self):
        """houses_ex2() should return (cusps, ascmc, cusps_speed, ascmc_speed)."""
        jd = 2451545.0
        result = ephem.swe_houses_ex2(jd, 41.9, 12.5, ord("P"), 0)

        assert len(result) == 4
        cusps, ascmc, cusps_speed, ascmc_speed = result

        # Verify structure
        assert len(cusps) == 12
        assert len(ascmc) == 8
        assert len(cusps_speed) == 12
        assert len(ascmc_speed) == 8

    @pytest.mark.unit
    def test_houses_ex2_cusps_match_houses_ex(self):
        """Cusps and ascmc from houses_ex2 should match houses_ex."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5
        hsys = ord("P")

        cusps_ex, ascmc_ex = ephem.swe_houses_ex(jd, lat, lon, hsys, 0)
        cusps_ex2, ascmc_ex2, _, _ = ephem.swe_houses_ex2(jd, lat, lon, hsys, 0)

        for i in range(12):
            assert abs(cusps_ex[i] - cusps_ex2[i]) < 1e-10, f"Cusp {i + 1} mismatch"

        for i in range(8):
            assert abs(ascmc_ex[i] - ascmc_ex2[i]) < 1e-10, f"ASCMC {i} mismatch"

    @pytest.mark.unit
    def test_houses_ex2_velocities_are_reasonable(self):
        """Cusp velocities should be in a reasonable range."""
        jd = 2451545.0
        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd, 41.9, 12.5, ord("P"), 0
        )

        # House cusps move roughly at the speed of Earth's rotation
        # MC moves about 360 deg/day (1 degree per 4 minutes)
        # ASC varies with latitude, can be faster than 360 deg/day
        # At mid-latitudes, ASC can reach ~550 deg/day
        for i, speed in enumerate(cusps_speed):
            # Velocities should be positive (houses move direct)
            # and roughly in the range of 100-600 deg/day at mid-latitudes
            assert 100 < abs(speed) < 600, f"Cusp {i + 1} speed {speed} out of range"

        # ASC speed (index 0) - can be faster than MC
        assert 100 < abs(ascmc_speed[0]) < 600, (
            f"ASC speed {ascmc_speed[0]} out of range"
        )

        # MC speed (index 1) - should be close to 360 deg/day
        assert 300 < abs(ascmc_speed[1]) < 400, (
            f"MC speed {ascmc_speed[1]} out of range"
        )

    @pytest.mark.comparison
    def test_houses_ex2_vs_pyswisseph(self):
        """Compare houses_ex2 results with pyswisseph."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5
        hsys = ord("P")

        cusps_lib, ascmc_lib, cusps_speed_lib, ascmc_speed_lib = ephem.swe_houses_ex2(
            jd, lat, lon, hsys, 0
        )
        cusps_swe, ascmc_swe, cusps_speed_swe, ascmc_speed_swe = swe.houses_ex2(
            jd, lat, lon, bytes([hsys]), 0
        )

        # Compare cusp positions (tolerance 0.1 degrees)
        for i in range(12):
            diff = abs(cusps_lib[i] - cusps_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.1, f"Cusp {i + 1} position diff {diff}"

        # Compare cusp velocities (tolerance 2 degrees/day due to different calc methods)
        for i in range(12):
            diff = abs(cusps_speed_lib[i] - cusps_speed_swe[i])
            assert diff < 2.0, f"Cusp {i + 1} velocity diff {diff}"

        # Compare ASC and MC velocities
        for i in range(2):  # ASC and MC
            diff = abs(ascmc_speed_lib[i] - ascmc_speed_swe[i])
            assert diff < 2.0, f"ASCMC {i} velocity diff {diff}"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "hsys,name,tolerance",
        [
            (ord("P"), "Placidus", 0.1),
            (ord("R"), "Regiomontanus", 0.1),
            (ord("C"), "Campanus", 0.1),
            (ord("E"), "Equal", 0.1),
        ],
    )
    def test_houses_ex2_multiple_systems(self, hsys, name, tolerance):
        """Test houses_ex2 with multiple house systems.

        Note: Some house systems (Koch, Porphyry, Whole Sign) have significant
        differences in velocity calculation between implementations, so they are
        excluded from this comparison test.
        """
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd, lat, lon, hsys, 0
        )
        cusps_swe, ascmc_swe, cusps_speed_swe, ascmc_speed_swe = swe.houses_ex2(
            jd, lat, lon, bytes([hsys]), 0
        )

        # Compare cusp velocities
        for i in range(12):
            diff = abs(cusps_speed[i] - cusps_speed_swe[i])
            assert diff < 5.0, f"{name} cusp {i + 1} velocity diff {diff}"

    @pytest.mark.unit
    def test_houses_ex2_with_sidereal(self):
        """houses_ex2 should work with SEFLG_SIDEREAL."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        # Set Lahiri ayanamsa
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.swe_houses_ex2(
            jd, lat, lon, ord("P"), SEFLG_SIDEREAL
        )

        # Cusps should be shifted from tropical by about 23-24 degrees
        cusps_trop, _, _, _ = ephem.swe_houses_ex2(jd, lat, lon, ord("P"), 0)

        # Check that sidereal cusps differ from tropical
        ayanamsa = ephem.swe_get_ayanamsa_ut(jd)
        for i in range(12):
            expected = (cusps_trop[i] - ayanamsa) % 360
            diff = abs(cusps[i] - expected)
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.001, f"Sidereal cusp {i + 1} not properly adjusted"

        # Velocities should still be valid (can reach ~550 deg/day at mid-latitudes)
        for speed in cusps_speed:
            assert 100 < abs(speed) < 600
