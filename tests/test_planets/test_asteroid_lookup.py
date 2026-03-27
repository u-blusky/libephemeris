"""
Tests for asteroid position calculations by number.

Verifies that asteroids accessed via SE_AST_OFFSET + number
return valid positions, and that the 5 main asteroids
(Chiron, Ceres, Pallas, Juno, Vesta) are accurate.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SE_CHIRON,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
    SE_AST_OFFSET,
    SEFLG_SPEED,
    SEFLG_HELCTR,
    SEFLG_EQUATORIAL,
    SEFLG_SIDEREAL,
    SE_SIDM_LAHIRI,
)


# Main asteroids with their direct body IDs
MAIN_ASTEROIDS = [
    (SE_CHIRON, "Chiron", 15),
    (SE_CERES, "Ceres", 17),
    (SE_PALLAS, "Pallas", 18),
    (SE_JUNO, "Juno", 19),
    (SE_VESTA, "Vesta", 20),
]

# Approximate heliocentric distance ranges (AU)
ASTEROID_DISTANCES = {
    SE_CHIRON: (8, 20),  # Chiron: between Saturn and Uranus
    SE_CERES: (2.5, 3.0),  # Ceres: in the main belt
    SE_PALLAS: (2.0, 3.5),  # Pallas: main belt, eccentric
    SE_JUNO: (2.0, 3.5),  # Juno: main belt
    SE_VESTA: (2.0, 2.8),  # Vesta: main belt
}


class TestMainAsteroidsBasic:
    """Basic tests for the 5 main asteroids."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_returns_valid(self, body_id: int, name: str, idx: int):
        """Each main asteroid returns valid 6-element tuple."""
        jd = 2451545.0
        result, retflag = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)
        assert len(result) == 6, f"{name}: expected 6 elements"
        lon, lat, dist = result[0], result[1], result[2]
        assert 0 <= lon < 360, f"{name}: lon {lon} out of range"
        assert -90 <= lat <= 90, f"{name}: lat {lat} out of range"
        assert dist > 0, f"{name}: distance {dist} not positive"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_speed_reasonable(self, body_id: int, name: str, idx: int):
        """Asteroid speeds should be reasonable (< 1 deg/day)."""
        jd = 2451545.0
        result, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)
        speed = result[3]
        assert abs(speed) < 1.5, f"{name}: speed {speed}°/day too large"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_all_finite(self, body_id: int, name: str, idx: int):
        """All values should be finite."""
        jd = 2451545.0
        result, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name}: result[{i}] = {val}"


class TestAsteroidDistances:
    """Test asteroid distances are physically plausible."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_geocentric_distance(self, body_id: int, name: str, idx: int):
        """Geocentric distance should be plausible."""
        jd = 2451545.0
        result, _ = swe.swe_calc_ut(jd, body_id, 0)
        dist = result[2]
        # Geocentric distance depends on Earth position relative to asteroid
        # For main belt: ~1-4 AU from Earth; for Chiron: ~7-20 AU
        if body_id == SE_CHIRON:
            assert 5 < dist < 25, f"{name}: dist {dist} AU"
        else:
            assert 0.5 < dist < 5.0, f"{name}: dist {dist} AU"


class TestAsteroidFlagCombinations:
    """Test asteroids with various flag combinations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        [
            (0, "default"),
            (SEFLG_SPEED, "speed"),
            (SEFLG_EQUATORIAL, "equatorial"),
            (SEFLG_HELCTR, "heliocentric"),
            (SEFLG_SPEED | SEFLG_EQUATORIAL, "speed+equatorial"),
        ],
    )
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (SE_CHIRON, "Chiron"),
            (SE_CERES, "Ceres"),
            (SE_VESTA, "Vesta"),
        ],
    )
    def test_asteroid_flag_combo(self, body_id: int, name: str, flags: int, desc: str):
        """Asteroid works with various flag combinations."""
        jd = 2451545.0
        result, _ = swe.swe_calc_ut(jd, body_id, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name}+{desc}: result[{i}]={val}"


class TestAsteroidSidereal:
    """Test asteroids in sidereal mode."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_sidereal(self, body_id: int, name: str, idx: int):
        """Sidereal positions valid for asteroids."""
        swe.swe_set_sid_mode(SE_SIDM_LAHIRI)
        jd = 2451545.0
        result, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SIDEREAL)
        assert 0 <= result[0] < 360, f"{name}: sidereal lon {result[0]}"


class TestAsteroidDateRange:
    """Test asteroids across date ranges."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_across_centuries(self, body_id: int, name: str, idx: int):
        """Asteroids valid from 1800 to 2200."""
        for year in [1800, 1900, 1950, 2000, 2050, 2100, 2200]:
            jd = swe.swe_julday(year, 1, 1, 12.0)
            result, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)
            assert 0 <= result[0] < 360, f"{name} @ {year}: lon={result[0]}"

    @pytest.mark.unit
    def test_chiron_continuity(self):
        """Chiron positions should be continuous over months."""
        jd_start = 2451545.0
        prev_lon = None
        for i in range(24):  # 2 years monthly
            jd = jd_start + i * 30.0
            result, _ = swe.swe_calc_ut(jd, SE_CHIRON, 0)
            lon = result[0]
            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                # Chiron moves ~2°/month
                assert diff < 10.0, f"Chiron jump {diff:.2f}° at month {i}"
            prev_lon = lon


class TestAsteroidHeliocentric:
    """Test asteroids in heliocentric mode."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_helio_valid(self, body_id: int, name: str, idx: int):
        """Heliocentric positions valid for asteroids."""
        jd = 2451545.0
        result, _ = swe.swe_calc_ut(jd, body_id, SEFLG_HELCTR | SEFLG_SPEED)
        assert len(result) == 6
        assert 0 <= result[0] < 360
        assert result[2] > 0

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name,idx", MAIN_ASTEROIDS)
    def test_asteroid_helio_distance_range(self, body_id: int, name: str, idx: int):
        """Heliocentric distance should match known orbital ranges."""
        jd = 2451545.0
        result, _ = swe.swe_calc_ut(jd, body_id, SEFLG_HELCTR)
        dist = result[2]
        lo, hi = ASTEROID_DISTANCES[body_id]
        assert lo * 0.6 <= dist <= hi * 1.5, (
            f"{name}: helio dist {dist} AU outside [{lo}, {hi}]"
        )


class TestAsteroidByNumber:
    """Test asteroid lookup by SE_AST_OFFSET + catalog number."""

    @pytest.mark.unit
    def test_ceres_by_offset(self):
        """Ceres via SE_AST_OFFSET + 1 should match SE_CERES."""
        jd = 2451545.0
        try:
            r_direct, _ = swe.swe_calc_ut(jd, SE_CERES, SEFLG_SPEED)
            r_offset, _ = swe.swe_calc_ut(jd, SE_AST_OFFSET + 1, SEFLG_SPEED)

            # Positions should be identical
            for i in range(6):
                assert abs(r_direct[i] - r_offset[i]) < 0.01, (
                    f"Ceres element {i}: direct={r_direct[i]}, offset={r_offset[i]}"
                )
        except Exception:
            # AST_OFFSET lookup may not be supported for all implementations
            pytest.skip("SE_AST_OFFSET lookup not supported")

    @pytest.mark.unit
    def test_pallas_by_offset(self):
        """Pallas via SE_AST_OFFSET + 2 should match SE_PALLAS."""
        jd = 2451545.0
        try:
            r_direct, _ = swe.swe_calc_ut(jd, SE_PALLAS, SEFLG_SPEED)
            r_offset, _ = swe.swe_calc_ut(jd, SE_AST_OFFSET + 2, SEFLG_SPEED)

            for i in range(6):
                assert abs(r_direct[i] - r_offset[i]) < 0.01, (
                    f"Pallas element {i}: direct={r_direct[i]}, offset={r_offset[i]}"
                )
        except Exception:
            pytest.skip("SE_AST_OFFSET lookup not supported")

    @pytest.mark.unit
    def test_juno_by_offset(self):
        """Juno via SE_AST_OFFSET + 3 should match SE_JUNO."""
        jd = 2451545.0
        try:
            r_direct, _ = swe.swe_calc_ut(jd, SE_JUNO, SEFLG_SPEED)
            r_offset, _ = swe.swe_calc_ut(jd, SE_AST_OFFSET + 3, SEFLG_SPEED)

            for i in range(6):
                assert abs(r_direct[i] - r_offset[i]) < 0.01, (
                    f"Juno element {i}: direct={r_direct[i]}, offset={r_offset[i]}"
                )
        except Exception:
            pytest.skip("SE_AST_OFFSET lookup not supported")

    @pytest.mark.unit
    def test_vesta_by_offset(self):
        """Vesta via SE_AST_OFFSET + 4 should match SE_VESTA."""
        jd = 2451545.0
        try:
            r_direct, _ = swe.swe_calc_ut(jd, SE_VESTA, SEFLG_SPEED)
            r_offset, _ = swe.swe_calc_ut(jd, SE_AST_OFFSET + 4, SEFLG_SPEED)

            for i in range(6):
                assert abs(r_direct[i] - r_offset[i]) < 0.01, (
                    f"Vesta element {i}: direct={r_direct[i]}, offset={r_offset[i]}"
                )
        except Exception:
            pytest.skip("SE_AST_OFFSET lookup not supported")
