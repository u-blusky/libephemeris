"""
Tests for the migration guide examples.

This test module verifies that all code examples from the migration guide
(docs/migration-guide.md) work correctly.

Note: EphemerisContext calculation tests are in test_context_thread_safety.py.
The tests here focus on the module-level API which is documented in the migration guide.
"""

import libephemeris as swe
from libephemeris import EphemerisContext
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SEFLG_SPEED,
    SEFLG_SIDEREAL,
    SE_SIDM_LAHIRI,
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_TRUE_CITRA,
)


class TestMigrationGuideBasicUsage:
    """Test basic usage examples from the migration guide."""

    def test_drop_in_replacement_import(self):
        """Test that libephemeris works as a drop-in replacement for swisseph."""
        # This is the main claim: you can replace `import swisseph as swe`
        # with `import libephemeris as swe`
        jd = swe.swe_julday(2000, 1, 1, 12.0)
        pos, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)

        assert pos is not None
        assert len(pos) == 6
        assert 0 <= pos[0] < 360  # Valid longitude

    def test_constants_from_module(self):
        """Test accessing constants from the module."""
        assert swe.SE_SUN == 0
        assert swe.SE_MOON == 1
        assert swe.SEFLG_SPEED == 256
        assert swe.SEFLG_SIDEREAL == 64 * 1024

    def test_prefixed_and_unprefixed_functions(self):
        """Test that both swe_ prefixed and unprefixed functions work."""
        jd = 2451545.0

        # Prefixed version
        pos1, _ = swe.swe_calc_ut(jd, SE_SUN, 0)

        # Unprefixed version
        pos2, _ = swe.calc_ut(jd, SE_SUN, 0)

        # Both should give identical results
        assert pos1[0] == pos2[0]
        assert pos1[1] == pos2[1]
        assert pos1[2] == pos2[2]


class TestMigrationGuideHouseCusps:
    """Test house cusp indexing as documented in the migration guide."""

    def test_house_cusps_0_indexed(self):
        """Test that libephemeris house cusps are 0-indexed."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5  # Rome

        cusps, ascmc = swe.swe_houses(jd, lat, lon, ord("P"))

        # Verify we get 12 cusps
        assert len(cusps) == 12

        # First house cusp is at index 0
        cusp1 = cusps[0]
        assert 0 <= cusp1 < 360

        # ascmc contains angles
        assert len(ascmc) >= 4
        asc = ascmc[0]
        mc = ascmc[1]
        assert 0 <= asc < 360
        assert 0 <= mc < 360


class TestMigrationGuidePrecision:
    """Test precision claims from the migration guide."""

    def test_planetary_longitude_valid_range(self):
        """Test that planetary longitudes are in valid range."""
        jd = 2451545.0

        for planet in [SE_SUN, SE_MOON]:
            pos, _ = swe.swe_calc_ut(jd, planet, SEFLG_SPEED)

            assert 0 <= pos[0] < 360, f"Planet {planet} longitude out of range"
            assert -90 <= pos[1] <= 90, f"Planet {planet} latitude out of range"
            assert pos[2] > 0, f"Planet {planet} distance should be positive"

    def test_velocities_are_computed(self):
        """Test that velocities are computed when SEFLG_SPEED is set."""
        jd = 2451545.0

        pos, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)

        # Sun moves about 1 degree per day
        assert 0.9 < pos[3] < 1.1, "Sun speed should be ~1 deg/day"

        pos, _ = swe.swe_calc_ut(jd, SE_MOON, SEFLG_SPEED)

        # Moon moves about 12-15 degrees per day
        assert 11 < pos[3] < 16, "Moon speed should be ~12-15 deg/day"


class TestMigrationGuideAyanamshas:
    """Test ayanamsha functionality from the migration guide."""

    def test_sidereal_mode_setting(self):
        """Test setting sidereal mode and calculating sidereal positions."""
        jd = 2451545.0

        # Set Lahiri ayanamsha
        swe.swe_set_sid_mode(SE_SIDM_LAHIRI)

        # Calculate sidereal position
        pos_sidereal, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)

        # Calculate tropical position
        pos_tropical, _ = swe.swe_calc_ut(jd, SE_SUN, 0)

        # Sidereal should be less than tropical (ayanamsha is positive)
        # The difference should be approximately the ayanamsha value
        ayanamsha = swe.swe_get_ayanamsa_ut(jd)
        expected_sidereal = (pos_tropical[0] - ayanamsha) % 360

        diff = abs(pos_sidereal[0] - expected_sidereal)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.001, f"Sidereal calculation mismatch: {diff}"

    def test_multiple_ayanamshas_give_different_results(self):
        """Test that different ayanamshas produce different sidereal positions."""
        jd = 2451545.0

        results = {}
        for mode in [SE_SIDM_FAGAN_BRADLEY, SE_SIDM_LAHIRI, SE_SIDM_TRUE_CITRA]:
            swe.swe_set_sid_mode(mode)
            pos, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)
            results[mode] = pos[0]

        # All three should be different
        assert len(set(round(v, 4) for v in results.values())) == 3


class TestMigrationGuideEphemerisContextState:
    """Test EphemerisContext state management from the migration guide."""

    def test_ephemeris_context_creation(self):
        """Test creating an EphemerisContext and setting state."""
        ctx = EphemerisContext()

        # Should be able to set state
        ctx.set_topo(12.5, 41.9, 0)  # Rome
        ctx.set_sid_mode(SE_SIDM_LAHIRI)

        # Verify state was set
        assert ctx.get_sid_mode() == SE_SIDM_LAHIRI
        topo = ctx.get_topo()
        assert topo is not None
        assert abs(topo.latitude.degrees - 41.9) < 0.01

    def test_ephemeris_context_isolation(self):
        """Test that EphemerisContext instances have isolated state."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        # Set different sidereal modes
        ctx1.set_sid_mode(SE_SIDM_LAHIRI)
        ctx2.set_sid_mode(SE_SIDM_FAGAN_BRADLEY)

        # Verify isolation
        assert ctx1.get_sid_mode() == SE_SIDM_LAHIRI
        assert ctx2.get_sid_mode() == SE_SIDM_FAGAN_BRADLEY

        # Set different locations
        ctx1.set_topo(12.5, 41.9, 0)  # Rome
        ctx2.set_topo(-0.1, 51.5, 0)  # London

        # Verify isolation
        topo1 = ctx1.get_topo()
        topo2 = ctx2.get_topo()

        assert topo1 is not None
        assert topo2 is not None
        assert abs(topo1.latitude.degrees - 41.9) < 0.01
        assert abs(topo2.latitude.degrees - 51.5) < 0.01


class TestMigrationGuideNotImplemented:
    """Test features documented as not implemented in the migration guide."""

    def test_fixed_star_velocities_are_small(self):
        """Test that fixed star velocities are small but non-zero due to precession.

        Fixed stars have small velocities (~3.8e-05 deg/day) due to precession
        of equinoxes, consistent with pyswisseph behavior.
        """
        jd = 2451545.0

        result = swe.fixstar_ut("Aldebaran", jd, SEFLG_SPEED)
        # fixstar_ut returns (pos, flag, starname)
        pos = result[0]

        # Velocities are very small but non-zero due to precession of equinoxes
        assert abs(pos[3]) < 0.001, "Fixed star lon velocity should be very small"
        assert abs(pos[4]) < 0.001, "Fixed star lat velocity should be very small"
        assert abs(pos[5]) < 0.001, "Fixed star dist velocity should be very small"
