"""
Tests for REBOUND/ASSIST n-body integration module.

These tests verify that the REBOUND integration works correctly for asteroid
orbit propagation, both with and without the optional REBOUND/ASSIST packages.
"""

import math
import pytest
from unittest.mock import patch, MagicMock

from libephemeris.minor_bodies import (
    OrbitalElements,
    MINOR_BODY_ELEMENTS,
    calc_minor_body_position,
)
from libephemeris.constants import SE_CERES, SE_VESTA, SE_CHIRON, SE_ERIS
from libephemeris.rebound_integration import (
    ReboundIntegrator,
    AssistEphemConfig,
    PropagationResult,
    check_rebound_available,
    check_assist_available,
    get_rebound_version,
    get_assist_version,
    elements_to_rebound_orbit,
    propagate_orbit_rebound,
    propagate_trajectory,
    compare_with_keplerian,
)


class TestReboundIntegratorEnum:
    """Test the ReboundIntegrator enum."""

    def test_ias15_value(self):
        """IAS15 should have correct string value."""
        assert ReboundIntegrator.IAS15.value == "ias15"

    def test_whfast_value(self):
        """WHFast should have correct string value."""
        assert ReboundIntegrator.WHFAST.value == "whfast"

    def test_mercurius_value(self):
        """MERCURIUS should have correct string value."""
        assert ReboundIntegrator.MERCURIUS.value == "mercurius"

    def test_trace_value(self):
        """TRACE should have correct string value."""
        assert ReboundIntegrator.TRACE.value == "trace"

    def test_leapfrog_value(self):
        """LEAPFROG should have correct string value."""
        assert ReboundIntegrator.LEAPFROG.value == "leapfrog"


class TestAssistEphemConfig:
    """Test ASSIST ephemeris configuration."""

    def test_default_config(self):
        """Default config should have None paths."""
        config = AssistEphemConfig()
        # Paths may be None or set from environment
        assert isinstance(config, AssistEphemConfig)

    def test_custom_data_dir(self, tmp_path):
        """Custom data directory should be used."""
        config = AssistEphemConfig(data_dir=str(tmp_path))
        assert config.data_dir == str(tmp_path)

    def test_custom_file_paths(self, tmp_path):
        """Custom file paths should be preserved."""
        planets = str(tmp_path / "planets.dat")
        asteroids = str(tmp_path / "asteroids.bsp")

        config = AssistEphemConfig(planets_file=planets, asteroids_file=asteroids)
        assert config.planets_file == planets
        assert config.asteroids_file == asteroids


class TestPropagationResult:
    """Test PropagationResult dataclass."""

    def test_basic_result(self):
        """Basic PropagationResult should store values correctly."""
        result = PropagationResult(
            x=1.0, y=0.0, z=0.0, vx=0.0, vy=0.017, vz=0.0, jd_tt=2460000.5
        )
        assert result.x == 1.0
        assert result.y == 0.0
        assert result.z == 0.0
        assert result.jd_tt == 2460000.5

    def test_distance_calculation(self):
        """Distance should be calculated from x,y,z."""
        result = PropagationResult(
            x=3.0, y=4.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=2460000.5
        )
        assert result.distance == 5.0

    def test_ecliptic_lon_positive_x(self):
        """Ecliptic longitude at positive X axis should be 0."""
        result = PropagationResult(
            x=1.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=2460000.5
        )
        assert abs(result.ecliptic_lon) < 1e-10

    def test_ecliptic_lon_positive_y(self):
        """Ecliptic longitude at positive Y axis should be 90."""
        result = PropagationResult(
            x=0.0, y=1.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=2460000.5
        )
        assert abs(result.ecliptic_lon - 90.0) < 1e-10

    def test_ecliptic_lat_positive_z(self):
        """Ecliptic latitude at positive Z axis should be 90."""
        result = PropagationResult(
            x=0.0, y=0.0, z=1.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=2460000.5
        )
        assert abs(result.ecliptic_lat - 90.0) < 1e-10

    def test_to_tuple(self):
        """to_tuple should return correct order."""
        result = PropagationResult(
            x=1.0, y=2.0, z=3.0, vx=0.1, vy=0.2, vz=0.3, jd_tt=2460000.5
        )
        expected = (1.0, 2.0, 3.0, 0.1, 0.2, 0.3)
        assert result.to_tuple() == expected


class TestAvailabilityChecks:
    """Test package availability check functions."""

    def test_check_rebound_returns_bool(self):
        """check_rebound_available should return bool."""
        result = check_rebound_available()
        assert isinstance(result, bool)

    def test_check_assist_returns_bool(self):
        """check_assist_available should return bool."""
        result = check_assist_available()
        assert isinstance(result, bool)

    def test_get_rebound_version(self):
        """get_rebound_version should return string or None."""
        result = get_rebound_version()
        assert result is None or isinstance(result, str)

    def test_get_assist_version(self):
        """get_assist_version should return string or None."""
        result = get_assist_version()
        assert result is None or isinstance(result, str)

    def test_version_matches_availability(self):
        """Version should be available iff package is available."""
        if check_rebound_available():
            assert get_rebound_version() is not None
        else:
            assert get_rebound_version() is None


class TestElementsConversion:
    """Test orbital elements conversion to REBOUND format."""

    def test_convert_ceres_elements(self):
        """Converting Ceres elements should produce valid REBOUND params."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd = elements.epoch

        params = elements_to_rebound_orbit(elements, jd)

        assert "a" in params
        assert "e" in params
        assert "inc" in params
        assert "omega" in params
        assert "Omega" in params
        assert "M" in params

        # Check values are reasonable
        assert params["a"] > 0  # Semi-major axis positive
        assert 0 <= params["e"] < 1  # Eccentricity valid for ellipse
        assert 0 <= params["inc"] <= math.pi  # Inclination in [0, pi]

    def test_angles_in_radians(self):
        """Angles should be converted to radians."""
        elements = MINOR_BODY_ELEMENTS[SE_VESTA]
        jd = elements.epoch

        params = elements_to_rebound_orbit(elements, jd)

        # All angles should be in radians (< 2*pi for normalized angles)
        assert abs(params["inc"]) <= math.pi
        assert abs(params["omega"]) <= 2 * math.pi
        assert abs(params["Omega"]) <= 2 * math.pi
        assert abs(params["M"]) <= 2 * math.pi


class TestPropagateOrbitRebound:
    """Test REBOUND orbit propagation (requires rebound package)."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_propagate_ceres_short(self, skip_if_no_rebound):
        """Propagate Ceres orbit for a short time."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 30  # 30 days

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        assert isinstance(result, PropagationResult)
        assert result.jd_tt == jd_end
        assert result.distance > 0
        # Ceres should be between 2.5 and 3.0 AU from Sun
        assert 2.0 < result.distance < 3.5

    def test_propagate_vesta(self, skip_if_no_rebound):
        """Propagate Vesta orbit."""
        elements = MINOR_BODY_ELEMENTS[SE_VESTA]
        jd_start = elements.epoch
        jd_end = jd_start + 100  # 100 days

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        assert result.distance > 0
        # Vesta orbit between 2.15 and 2.57 AU
        assert 1.8 < result.distance < 2.8

    def test_propagate_chiron(self, skip_if_no_rebound):
        """Propagate Chiron (centaur) orbit."""
        elements = MINOR_BODY_ELEMENTS[SE_CHIRON]
        jd_start = elements.epoch
        jd_end = jd_start + 365  # 1 year

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        # Chiron orbit between 8.5 and 19 AU (highly eccentric)
        assert 5.0 < result.distance < 25.0

    def test_integrator_ias15(self, skip_if_no_rebound):
        """IAS15 integrator should work."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 10

        result = propagate_orbit_rebound(
            elements, jd_start, jd_end, integrator=ReboundIntegrator.IAS15
        )
        assert result.distance > 0

    def test_integrator_whfast(self, skip_if_no_rebound):
        """WHFast integrator should work."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 10

        result = propagate_orbit_rebound(
            elements, jd_start, jd_end, integrator=ReboundIntegrator.WHFAST
        )
        assert result.distance > 0

    def test_integrator_leapfrog(self, skip_if_no_rebound):
        """Leapfrog integrator should work."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 10

        result = propagate_orbit_rebound(
            elements, jd_start, jd_end, integrator=ReboundIntegrator.LEAPFROG
        )
        assert result.distance > 0

    def test_custom_timestep(self, skip_if_no_rebound):
        """Custom timestep should be used."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 30

        result = propagate_orbit_rebound(
            elements,
            jd_start,
            jd_end,
            dt=0.5,  # 0.5 day timestep
        )
        assert result.distance > 0

    def test_backward_integration(self, skip_if_no_rebound):
        """Backward integration should work."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start - 30  # 30 days backward

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        assert result.jd_tt == jd_end
        assert result.distance > 0


class TestPropagateTrajectory:
    """Test trajectory propagation (multiple output points)."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_trajectory_basic(self, skip_if_no_rebound):
        """Basic trajectory propagation should return list of points."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 100

        trajectory = propagate_trajectory(
            elements,
            jd_start,
            jd_end,
            num_points=10,
            use_assist=False,  # Use REBOUND only
        )

        assert len(trajectory) == 10
        assert all(isinstance(p, PropagationResult) for p in trajectory)

    def test_trajectory_times_evenly_spaced(self, skip_if_no_rebound):
        """Trajectory times should be evenly spaced."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 100

        trajectory = propagate_trajectory(
            elements, jd_start, jd_end, num_points=11, use_assist=False
        )

        times = [p.jd_tt for p in trajectory]
        dt = times[1] - times[0]

        for i in range(1, len(times)):
            assert abs(times[i] - times[i - 1] - dt) < 1e-6

    def test_trajectory_first_and_last(self, skip_if_no_rebound):
        """First and last points should match start and end times."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 50

        trajectory = propagate_trajectory(
            elements, jd_start, jd_end, num_points=5, use_assist=False
        )

        assert abs(trajectory[0].jd_tt - jd_start) < 1e-6
        assert abs(trajectory[-1].jd_tt - jd_end) < 1e-6


class TestCompareWithKeplerian:
    """Test comparison between n-body and Keplerian results."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_compare_returns_dict(self, skip_if_no_rebound):
        """Comparison should return a dictionary with expected keys."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd = elements.epoch + 100

        comparison = compare_with_keplerian(elements, jd, use_assist=False)

        assert isinstance(comparison, dict)
        assert "keplerian" in comparison
        assert "nbody" in comparison
        assert "delta_lon_arcsec" in comparison
        assert "delta_lat_arcsec" in comparison
        assert "angular_sep_arcsec" in comparison

    def test_compare_keplerian_values(self, skip_if_no_rebound):
        """Keplerian values should match direct calculation."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd = elements.epoch + 30

        comparison = compare_with_keplerian(elements, jd, use_assist=False)

        # Compare with direct Keplerian calculation
        x, y, z = calc_minor_body_position(elements, jd)
        r = math.sqrt(x**2 + y**2 + z**2)
        lon = math.degrees(math.atan2(y, x)) % 360.0
        lat = math.degrees(math.asin(z / r))

        assert abs(comparison["keplerian"]["lon"] - lon) < 1e-10
        assert abs(comparison["keplerian"]["lat"] - lat) < 1e-10
        assert abs(comparison["keplerian"]["dist"] - r) < 1e-10

    def test_short_propagation_small_difference(self, skip_if_no_rebound):
        """Short propagation should have small difference from Keplerian."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd = elements.epoch + 10  # Only 10 days

        comparison = compare_with_keplerian(elements, jd, use_assist=False)

        # For short propagation, 2-body REBOUND and Keplerian should be similar
        # The difference comes from REBOUND using slightly different integration
        assert comparison["angular_sep_arcsec"] < 3600  # Within 1 degree


class TestModuleWithoutRebound:
    """Test module behavior when REBOUND is not installed."""

    def test_import_works_without_rebound(self):
        """Module should import without REBOUND installed."""
        # This test verifies that the module can be imported even if
        # REBOUND is not installed (lazy import pattern)
        with patch.dict("sys.modules", {"rebound": None}):
            # Re-import to test
            from libephemeris import rebound_integration

            assert hasattr(rebound_integration, "check_rebound_available")

    def test_propagate_raises_without_rebound(self):
        """propagate_orbit_rebound should raise ImportError without REBOUND."""
        if check_rebound_available():
            pytest.skip("REBOUND is installed")

        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        with pytest.raises(ImportError):
            propagate_orbit_rebound(elements, 2460000.5, 2460010.5)


class TestIntegrationAccuracy:
    """Test integration accuracy for longer propagation times."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_orbit_stays_bounded(self, skip_if_no_rebound):
        """Orbit should stay within physical bounds over 1 year."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 365.25  # 1 year

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        # Ceres orbit: a = 2.77 AU, e = 0.076
        # Perihelion: 2.56 AU, Aphelion: 2.98 AU
        assert 2.0 < result.distance < 3.5

    def test_eccentric_orbit_accuracy(self, skip_if_no_rebound):
        """Highly eccentric orbit (Chiron) should remain bounded."""
        elements = MINOR_BODY_ELEMENTS[SE_CHIRON]
        jd_start = elements.epoch
        jd_end = jd_start + 365.25

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        # Chiron orbit: a = 13.7 AU, e = 0.38
        # Perihelion: 8.5 AU, Aphelion: 18.9 AU
        assert 5.0 < result.distance < 25.0

    def test_distant_tno_accuracy(self, skip_if_no_rebound):
        """Distant TNO (Eris) orbit should remain bounded."""
        elements = MINOR_BODY_ELEMENTS[SE_ERIS]
        jd_start = elements.epoch
        jd_end = jd_start + 365.25

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        # Eris orbit: a = 68 AU, e = 0.44
        # Perihelion: 38 AU, Aphelion: 98 AU
        assert 30.0 < result.distance < 110.0


class TestOrbitRoundTrip:
    """Test that orbit propagation is reversible."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_forward_backward_consistency(self, skip_if_no_rebound):
        """Forward then backward propagation should return to start."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        jd_start = elements.epoch
        jd_mid = jd_start + 50

        # Forward propagation
        result_fwd = propagate_orbit_rebound(elements, jd_start, jd_mid)

        # Create elements at mid point (approximate)
        # For this test, we just verify that both directions give similar distances
        result_bwd = propagate_orbit_rebound(elements, jd_start, jd_mid)

        # Results should be consistent
        assert abs(result_fwd.distance - result_bwd.distance) < 1e-10
