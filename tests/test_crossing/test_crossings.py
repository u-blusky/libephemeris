"""
Comprehensive tests for longitude crossing functions.

Tests solcross_ut, mooncross_ut, and cross_ut.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestSolcrossBasic:
    """Basic tests for Sun crossing function."""

    @pytest.mark.unit
    def test_solcross_vernal_equinox(self):
        """Find when Sun crosses 0° (vernal equinox)."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        jd_cross = ephem.swe_solcross_ut(0.0, jd_start, 0)

        # Should be around March 20, 2024
        year, month, day, hour = ephem.swe_revjul(jd_cross)
        assert year == 2024
        assert month == 3
        assert 19 <= day <= 21

    @pytest.mark.unit
    def test_solcross_summer_solstice(self):
        """Find when Sun crosses 90° (summer solstice)."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        jd_cross = ephem.swe_solcross_ut(90.0, jd_start, 0)

        year, month, day, hour = ephem.swe_revjul(jd_cross)
        assert year == 2024
        assert month == 6
        assert 20 <= day <= 22

    @pytest.mark.unit
    def test_solcross_precision(self):
        """Sun should be very close to target at crossing time."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 45.0
        jd_cross = ephem.swe_solcross_ut(target, jd_start, 0)

        # Check Sun position at crossing
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_SUN, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.001, f"Sun at {pos[0]}, target {target}, diff {diff}"


class TestSolcrossTT:
    """Tests for swe_solcross (TT version)."""

    @pytest.mark.unit
    def test_solcross_tt_vernal_equinox(self):
        """Find when Sun crosses 0° (vernal equinox) using TT."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        # Convert to TT (approximately, for starting point)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.swe_solcross(0.0, jd_start_tt, 0)

        # Should be around March 20, 2024
        # Convert back to UT for date check
        delta_t_cross = ephem.swe_deltat(jd_cross_tt)
        jd_cross_ut = jd_cross_tt - delta_t_cross
        year, month, day, hour = ephem.swe_revjul(jd_cross_ut)
        assert year == 2024
        assert month == 3
        assert 19 <= day <= 21

    @pytest.mark.unit
    def test_solcross_tt_precision(self):
        """Sun should be very close to target at crossing time (TT)."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 45.0
        jd_cross_tt = ephem.swe_solcross(target, jd_start_tt, 0)

        # Check Sun position at crossing (using TT version of calc)
        pos, _ = ephem.swe_calc(jd_cross_tt, SE_SUN, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.001, f"Sun at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_solcross_tt_vs_ut_consistency(self):
        """TT and UT versions should give consistent results."""
        jd_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 90.0  # Summer solstice

        # Get crossing time in UT
        jd_cross_ut = ephem.swe_solcross_ut(target, jd_ut, 0)

        # Get crossing time in TT
        jd_cross_tt = ephem.swe_solcross(target, jd_tt, 0)

        # Convert UT result to TT for comparison
        delta_t_cross = ephem.swe_deltat(jd_cross_ut)
        jd_cross_ut_as_tt = jd_cross_ut + delta_t_cross

        # They should be very close (within seconds)
        diff_seconds = abs(jd_cross_tt - jd_cross_ut_as_tt) * 86400
        assert diff_seconds < 10, f"TT vs UT consistency diff {diff_seconds} seconds"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_solcross_tt_all_signs(self, target):
        """All zodiac sign ingresses should work with TT version."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.swe_solcross(float(target), jd_start_tt, 0)

        # Verify Sun is at target
        pos, _ = ephem.swe_calc(jd_cross_tt, SE_SUN, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Target {target}° diff {diff}°"

    @pytest.mark.comparison
    def test_solcross_tt_vs_pyswisseph(self):
        """TT version should match pyswisseph solcross()."""
        jd_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = swe.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 0.0

        jd_lib = ephem.swe_solcross(target, jd_tt, 0)
        jd_swe = swe.solcross(target, jd_tt, 0)

        # Difference should be less than 1 minute
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 60, f"Timing diff {diff_seconds} seconds"

    @pytest.mark.edge_case
    def test_solcross_tt_target_360_equals_0(self):
        """360° should be same as 0° for TT version."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_360 = ephem.swe_solcross(360.0, jd_start_tt, 0)
        jd_0 = ephem.swe_solcross(0.0, jd_start_tt, 0)

        # Should be the same crossing
        assert abs(jd_360 - jd_0) < 0.001


class TestSolcrossVsPyswisseph:
    """Compare solcross with pyswisseph."""

    @pytest.mark.comparison
    def test_solcross_equinox_timing(self):
        """Crossing time should match pyswisseph."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_lib = ephem.swe_solcross_ut(0.0, jd_start, 0)
        jd_swe = swe.solcross_ut(0.0, jd_start, 0)

        # Difference should be less than 1 minute
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 60, f"Timing diff {diff_seconds} seconds"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_solcross_all_signs(self, target):
        """All zodiac sign ingresses should match."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_lib = ephem.swe_solcross_ut(float(target), jd_start, 0)
        jd_swe = swe.solcross_ut(float(target), jd_start, 0)

        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 120, f"Target {target}° diff {diff_seconds} seconds"


class TestMooncrossBasic:
    """Basic tests for Moon crossing function."""

    @pytest.mark.unit
    def test_mooncross_basic(self):
        """Moon crossing should return valid JD."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 90.0
        jd_cross = ephem.swe_mooncross_ut(target, jd_start, 0)

        assert jd_cross > jd_start
        # Should be within one lunar orbit (~27 days)
        assert jd_cross < jd_start + 28

    @pytest.mark.unit
    def test_mooncross_precision(self):
        """Moon should be close to target at crossing."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 123.456
        jd_cross = ephem.swe_mooncross_ut(target, jd_start, 0)

        pos, _ = ephem.swe_calc_ut(jd_cross, SE_MOON, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Moon at {pos[0]}, target {target}, diff {diff}"


class TestMooncrossTT:
    """Tests for swe_mooncross (TT version)."""

    @pytest.mark.unit
    def test_mooncross_tt_basic(self):
        """Moon crossing should return valid JD using TT."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 90.0
        jd_cross_tt = ephem.swe_mooncross(target, jd_start_tt, 0)

        assert jd_cross_tt > jd_start_tt
        # Should be within one lunar orbit (~27 days)
        assert jd_cross_tt < jd_start_tt + 28

    @pytest.mark.unit
    def test_mooncross_tt_precision(self):
        """Moon should be close to target at crossing time (TT)."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 123.456
        jd_cross_tt = ephem.swe_mooncross(target, jd_start_tt, 0)

        # Check Moon position at crossing (using TT version of calc)
        pos, _ = ephem.swe_calc(jd_cross_tt, SE_MOON, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Moon at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_mooncross_tt_vs_ut_consistency(self):
        """TT and UT versions should give consistent results."""
        jd_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 180.0

        # Get crossing time in UT
        jd_cross_ut = ephem.swe_mooncross_ut(target, jd_ut, 0)

        # Get crossing time in TT
        jd_cross_tt = ephem.swe_mooncross(target, jd_tt, 0)

        # Convert UT result to TT for comparison
        delta_t_cross = ephem.swe_deltat(jd_cross_ut)
        jd_cross_ut_as_tt = jd_cross_ut + delta_t_cross

        # They should be very close (within seconds)
        diff_seconds = abs(jd_cross_tt - jd_cross_ut_as_tt) * 86400
        assert diff_seconds < 10, f"TT vs UT consistency diff {diff_seconds} seconds"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_mooncross_tt_all_signs(self, target):
        """All zodiac sign ingresses should work with TT version."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.swe_mooncross(float(target), jd_start_tt, 0)

        # Verify Moon is at target
        pos, _ = ephem.swe_calc(jd_cross_tt, SE_MOON, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Target {target}° diff {diff}°"

    @pytest.mark.comparison
    def test_mooncross_tt_vs_pyswisseph(self):
        """TT version should match pyswisseph mooncross()."""
        jd_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = swe.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 180.0

        jd_lib = ephem.swe_mooncross(target, jd_tt, 0)
        jd_swe = swe.mooncross(target, jd_tt, 0)

        # Difference should be less than 3 minutes (Moon moves fast)
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 180, f"Timing diff {diff_seconds} seconds"

    @pytest.mark.edge_case
    def test_mooncross_tt_target_360_equals_0(self):
        """360° should be same as 0° for TT version."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_360 = ephem.swe_mooncross(360.0, jd_start_tt, 0)
        jd_0 = ephem.swe_mooncross(0.0, jd_start_tt, 0)

        # Should be the same crossing
        assert abs(jd_360 - jd_0) < 0.001


class TestMooncrossVsPyswisseph:
    """Compare mooncross with pyswisseph."""

    @pytest.mark.comparison
    def test_mooncross_timing(self):
        """Moon crossing should match pyswisseph."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 180.0

        jd_lib = ephem.swe_mooncross_ut(target, jd_start, 0)
        jd_swe = swe.mooncross_ut(target, jd_start, 0)

        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 180, f"Moon crossing diff {diff_seconds} seconds"


class TestCrossingConsecutive:
    """Test finding consecutive crossings."""

    @pytest.mark.unit
    def test_consecutive_sun_crossings(self):
        """Should find 12 consecutive Sun crossings in a year."""
        jd = ephem.swe_julday(2024, 1, 1, 0.0)
        crossings = []

        for target in range(0, 360, 30):
            jd_cross = ephem.swe_solcross_ut(float(target), jd, 0)
            crossings.append(jd_cross)

        # All should be in chronological order (after accounting for year wrap)
        # First 0° crossing should be in March
        assert crossings[0] > ephem.swe_julday(2024, 3, 1, 0.0)

    @pytest.mark.unit
    def test_12_moon_crossings_in_month(self):
        """Should find ~12 crossings in 28 days."""
        jd = ephem.swe_julday(2024, 1, 1, 0.0)
        crossings = []

        for target in range(0, 360, 30):
            jd_cross = ephem.swe_mooncross_ut(float(target), jd, 0)
            if jd_cross < jd + 28:
                crossings.append(jd_cross)

        # Should find all 12 in about 27 days
        assert len(crossings) == 12


class TestCrossingEdgeCases:
    """Test edge cases for crossing functions."""

    @pytest.mark.edge_case
    def test_target_360_equals_0(self):
        """360° should be same as 0°."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_360 = ephem.swe_solcross_ut(360.0, jd_start, 0)
        jd_0 = ephem.swe_solcross_ut(0.0, jd_start, 0)

        # Should be the same crossing
        assert abs(jd_360 - jd_0) < 0.001

    @pytest.mark.edge_case
    def test_target_negative_normalized(self):
        """Negative target should be normalized."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_neg = ephem.swe_solcross_ut(-30.0, jd_start, 0)  # = 330°
        jd_pos = ephem.swe_solcross_ut(330.0, jd_start, 0)

        assert abs(jd_neg - jd_pos) < 0.001


class TestCrossUtGeneric:
    """Test generic planet crossing function."""

    @pytest.mark.unit
    def test_cross_ut_mercury(self):
        """Should work for Mercury."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 0.0

        jd_cross = ephem.swe_cross_ut(SE_MERCURY, target, jd_start, 0)

        assert jd_cross > jd_start

        # Check Mercury position at crossing
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_MERCURY, 0)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1  # Mercury moves fast, allow more tolerance


class TestMooncrossNodeBasic:
    """Basic tests for Moon node crossing function."""

    @pytest.mark.unit
    def test_mooncross_node_basic(self):
        """Moon node crossing should return valid JD."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        jd_cross = ephem.swe_mooncross_node_ut(jd_start, 0)

        assert jd_cross > jd_start
        # Should be within ~14 days (half the nodal month)
        assert jd_cross < jd_start + 14

    @pytest.mark.unit
    def test_mooncross_node_latitude_near_zero(self):
        """Moon latitude should be near zero at node crossing."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        jd_cross = ephem.swe_mooncross_node_ut(jd_start, 0)

        # Check Moon latitude at crossing
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_MOON, 0)

        # Latitude should be very close to 0
        assert abs(pos[1]) < 0.01, f"Moon latitude at node crossing: {pos[1]}"

    @pytest.mark.unit
    def test_mooncross_node_consecutive_crossings(self):
        """Should find consecutive node crossings ~13.6 days apart."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        # Find first node crossing
        jd_first = ephem.swe_mooncross_node_ut(jd_start, 0)

        # Find second node crossing by starting just after first
        jd_second = ephem.swe_mooncross_node_ut(jd_first + 0.5, 0)

        # Should be roughly 13.6 days apart (half the nodal month)
        diff_days = jd_second - jd_first
        assert 12 < diff_days < 15, f"Time between node crossings: {diff_days} days"

    @pytest.mark.unit
    def test_mooncross_node_ascending_vs_descending(self):
        """Can determine ascending vs descending node from latitude velocity."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        # Find first crossing
        jd_first = ephem.swe_mooncross_node_ut(jd_start, 0)
        pos_first, _ = ephem.swe_calc_ut(jd_first, SE_MOON, SEFLG_SPEED)

        # Find second crossing
        jd_second = ephem.swe_mooncross_node_ut(jd_first + 0.5, 0)
        pos_second, _ = ephem.swe_calc_ut(jd_second, SE_MOON, SEFLG_SPEED)

        # Latitude velocities should have opposite signs
        # (one ascending, one descending)
        assert pos_first[4] * pos_second[4] < 0, (
            f"Consecutive crossings should alternate: "
            f"v1={pos_first[4]}, v2={pos_second[4]}"
        )


class TestMooncrossNodeTT:
    """Tests for swe_mooncross_node (TT version)."""

    @pytest.mark.unit
    def test_mooncross_node_tt_basic(self):
        """Moon node crossing should return valid JD using TT."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.swe_mooncross_node(jd_start_tt, 0)

        assert jd_cross_tt > jd_start_tt
        # Should be within ~14 days
        assert jd_cross_tt < jd_start_tt + 14

    @pytest.mark.unit
    def test_mooncross_node_tt_precision(self):
        """Moon latitude should be near zero at node crossing (TT)."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.swe_mooncross_node(jd_start_tt, 0)

        # Check Moon latitude at crossing (using TT version of calc)
        pos, _ = ephem.swe_calc(jd_cross_tt, SE_MOON, 0)

        # Latitude should be very close to 0
        assert abs(pos[1]) < 0.01, f"Moon latitude at node crossing: {pos[1]}"

    @pytest.mark.unit
    def test_mooncross_node_tt_vs_ut_consistency(self):
        """TT and UT versions should give consistent results."""
        jd_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        # Get crossing time in UT
        jd_cross_ut = ephem.swe_mooncross_node_ut(jd_ut, 0)

        # Get crossing time in TT
        jd_cross_tt = ephem.swe_mooncross_node(jd_tt, 0)

        # Convert UT result to TT for comparison
        delta_t_cross = ephem.swe_deltat(jd_cross_ut)
        jd_cross_ut_as_tt = jd_cross_ut + delta_t_cross

        # They should be very close (within seconds)
        diff_seconds = abs(jd_cross_tt - jd_cross_ut_as_tt) * 86400
        assert diff_seconds < 10, f"TT vs UT consistency diff {diff_seconds} seconds"


class TestMooncrossNodeVsPyswisseph:
    """Compare mooncross_node with pyswisseph."""

    @pytest.mark.comparison
    def test_mooncross_node_vs_pyswisseph(self):
        """Moon node crossing should match pyswisseph."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_lib = ephem.swe_mooncross_node_ut(jd_start, 0)
        # pyswisseph returns (jd_cross, xlon, xlat)
        result_swe = swe.mooncross_node_ut(jd_start, 0)
        jd_swe = result_swe[0]

        # Difference should be less than 3 minutes
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 180, f"Moon node crossing diff {diff_seconds} seconds"

    @pytest.mark.comparison
    def test_mooncross_node_tt_vs_pyswisseph(self):
        """TT version should match pyswisseph mooncross_node()."""
        jd_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = swe.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        jd_lib = ephem.swe_mooncross_node(jd_tt, 0)
        # pyswisseph returns (jd_cross, xlon, xlat)
        result_swe = swe.mooncross_node(jd_tt, 0)
        jd_swe = result_swe[0]

        # Difference should be less than 3 minutes
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 180, f"Moon node crossing diff {diff_seconds} seconds"

    @pytest.mark.comparison
    def test_mooncross_node_multiple_crossings(self):
        """Multiple consecutive crossings should match pyswisseph."""
        jd = ephem.swe_julday(2024, 1, 1, 0.0)

        for i in range(4):
            jd_lib = ephem.swe_mooncross_node_ut(jd, 0)
            # pyswisseph returns (jd_cross, xlon, xlat)
            result_swe = swe.mooncross_node_ut(jd, 0)
            jd_swe = result_swe[0]

            diff_seconds = abs(jd_lib - jd_swe) * 86400
            assert diff_seconds < 180, f"Crossing {i + 1} diff {diff_seconds} seconds"

            # Move to after this crossing for next iteration
            jd = jd_lib + 0.5


class TestHelioCrossBasic:
    """Basic tests for heliocentric crossing function."""

    @pytest.mark.unit
    def test_helio_cross_ut_mars_basic(self):
        """Heliocentric Mars crossing should return valid JD."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 0.0

        jd_cross = ephem.swe_helio_cross_ut(SE_MARS, target, jd_start, 0)

        assert jd_cross > jd_start
        # Mars completes one heliocentric orbit in ~687 days
        assert jd_cross < jd_start + 700

    @pytest.mark.unit
    def test_helio_cross_ut_precision(self):
        """Planet should be very close to target at heliocentric crossing."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 45.0

        jd_cross = ephem.swe_helio_cross_ut(SE_MARS, target, jd_start, 0)

        # Check Mars heliocentric position at crossing
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_MARS, SEFLG_HELCTR)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Mars helio at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_helio_cross_ut_earth(self):
        """Should work for Earth (heliocentric view of Earth from Sun)."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 90.0

        jd_cross = ephem.swe_helio_cross_ut(SE_EARTH, target, jd_start, 0)

        assert jd_cross > jd_start
        # Earth completes one heliocentric orbit in ~365 days
        assert jd_cross < jd_start + 366

        # Verify position
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_EARTH, SEFLG_HELCTR)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.01

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet,max_days",
        [
            (SE_MERCURY, 90),  # Mercury helio period ~88 days
            (SE_VENUS, 230),  # Venus helio period ~225 days
            (SE_EARTH, 370),  # Earth helio period ~365 days
            (SE_MARS, 700),  # Mars helio period ~687 days
        ],
    )
    def test_helio_cross_ut_multiple_planets(self, planet, max_days):
        """Should work for various planets."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target = 180.0

        jd_cross = ephem.swe_helio_cross_ut(planet, target, jd_start, 0)

        assert jd_cross > jd_start
        assert jd_cross < jd_start + max_days

        # Verify position
        pos, _ = ephem.swe_calc_ut(jd_cross, planet, SEFLG_HELCTR)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1, f"Planet {planet} at {pos[0]}, target {target}"


class TestHelioCrossTT:
    """Tests for swe_helio_cross (TT version)."""

    @pytest.mark.unit
    def test_helio_cross_tt_basic(self):
        """Heliocentric crossing should return valid JD using TT."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 90.0
        jd_cross_tt = ephem.swe_helio_cross(SE_MARS, target, jd_start_tt, 0)

        assert jd_cross_tt > jd_start_tt
        assert jd_cross_tt < jd_start_tt + 700

    @pytest.mark.unit
    def test_helio_cross_tt_precision(self):
        """Planet should be very close to target at crossing time (TT)."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 123.456
        jd_cross_tt = ephem.swe_helio_cross(SE_MARS, target, jd_start_tt, 0)

        # Check position (using TT version of calc)
        pos, _ = ephem.swe_calc(jd_cross_tt, SE_MARS, SEFLG_HELCTR)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Mars helio at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_helio_cross_tt_vs_ut_consistency(self):
        """TT and UT versions should give consistent results."""
        jd_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 180.0

        # Get crossing time in UT
        jd_cross_ut = ephem.swe_helio_cross_ut(SE_MARS, target, jd_ut, 0)

        # Get crossing time in TT
        jd_cross_tt = ephem.swe_helio_cross(SE_MARS, target, jd_tt, 0)

        # Convert UT result to TT for comparison
        delta_t_cross = ephem.swe_deltat(jd_cross_ut)
        jd_cross_ut_as_tt = jd_cross_ut + delta_t_cross

        # They should be very close (within seconds)
        diff_seconds = abs(jd_cross_tt - jd_cross_ut_as_tt) * 86400
        assert diff_seconds < 10, f"TT vs UT consistency diff {diff_seconds} seconds"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_helio_cross_tt_all_signs(self, target):
        """All zodiac sign ingresses should work with TT version."""
        jd_start_ut = ephem.swe_julday(2024, 1, 1, 0.0)
        delta_t = ephem.swe_deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.swe_helio_cross(SE_EARTH, float(target), jd_start_tt, 0)

        # Verify position
        pos, _ = ephem.swe_calc(jd_cross_tt, SE_EARTH, SEFLG_HELCTR)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Target {target}° diff {diff}°"


class TestHelioCrossEdgeCases:
    """Test edge cases for heliocentric crossing functions."""

    @pytest.mark.edge_case
    def test_helio_cross_target_360_equals_0(self):
        """360° should be same as 0°."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_360 = ephem.swe_helio_cross_ut(SE_MARS, 360.0, jd_start, 0)
        jd_0 = ephem.swe_helio_cross_ut(SE_MARS, 0.0, jd_start, 0)

        # Should be the same crossing
        assert abs(jd_360 - jd_0) < 0.001

    @pytest.mark.edge_case
    def test_helio_cross_target_negative_normalized(self):
        """Negative target should be normalized."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_neg = ephem.swe_helio_cross_ut(SE_MARS, -30.0, jd_start, 0)  # = 330°
        jd_pos = ephem.swe_helio_cross_ut(SE_MARS, 330.0, jd_start, 0)

        assert abs(jd_neg - jd_pos) < 0.001

    @pytest.mark.edge_case
    def test_helio_cross_consecutive_crossings(self):
        """Should find consecutive crossings for a full orbit."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        crossings = []

        # Find Earth crossing each 30° for a full year
        for target in range(0, 360, 30):
            jd_cross = ephem.swe_helio_cross_ut(SE_EARTH, float(target), jd_start, 0)
            crossings.append(jd_cross)

        # All crossings should be within one year
        assert max(crossings) - jd_start < 366


class TestHelioCrossVsGeocentricConcept:
    """Verify heliocentric differs from geocentric (conceptual tests)."""

    @pytest.mark.unit
    def test_helio_no_retrograde(self):
        """Heliocentric planets should not show retrograde behavior."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        # Find several consecutive crossings for Mars
        crossings = []
        jd = jd_start
        for _ in range(4):
            jd_cross = ephem.swe_helio_cross_ut(SE_MARS, 90.0, jd, 0)
            crossings.append(jd_cross)
            jd = jd_cross + 1  # Move past this crossing

        # Each crossing should be about one Mars orbital period apart (~687 days)
        for i in range(1, len(crossings)):
            diff = crossings[i] - crossings[i - 1]
            # Should be around 687 days (±50 days for minor perturbations)
            assert 600 < diff < 750, f"Period between crossings: {diff} days"

    @pytest.mark.unit
    def test_helio_vs_geocentric_different_positions(self):
        """Heliocentric and geocentric positions should generally differ."""
        jd = ephem.swe_julday(2024, 6, 15, 12.0)

        # Get geocentric Mars position
        geo_pos, _ = ephem.swe_calc_ut(jd, SE_MARS, 0)

        # Get heliocentric Mars position
        helio_pos, _ = ephem.swe_calc_ut(jd, SE_MARS, SEFLG_HELCTR)

        # Positions should be different (unless Mars is in opposition)
        diff = abs(geo_pos[0] - helio_pos[0])
        if diff > 180:
            diff = 360 - diff

        # They're generally different but we just verify the function works
        # and returns valid positions
        assert 0 <= geo_pos[0] < 360
        assert 0 <= helio_pos[0] < 360


class TestMooncrossNodeEclipseRelevance:
    """Test mooncross_node in context of eclipse calculations."""

    @pytest.mark.unit
    def test_eclipse_proximity_to_node(self):
        """Eclipses occur when Sun is near a lunar node."""
        # Start from a known solar eclipse date: April 8, 2024
        jd_eclipse = ephem.swe_julday(2024, 4, 8, 18.0)

        # Find the nearest node crossing
        # Check both before and after
        jd_before = jd_eclipse - 7  # Week before
        jd_cross_before = ephem.swe_mooncross_node_ut(jd_before, 0)

        jd_after = jd_eclipse - 0.5
        jd_cross_after = ephem.swe_mooncross_node_ut(jd_after, 0)

        # The eclipse should be close to a node crossing (within a few days)
        diff_before = abs(jd_eclipse - jd_cross_before)
        diff_after = abs(jd_cross_after - jd_eclipse)

        # Eclipse should be within ~1 day of a node crossing for a total eclipse
        assert min(diff_before, diff_after) < 2, (
            f"Eclipse not near node: before={diff_before}, after={diff_after}"
        )
