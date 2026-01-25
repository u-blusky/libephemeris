"""
Tests for parabolic and hyperbolic orbit handling in Kepler equation solvers.

Tests verify:
- Elliptic orbits (e < 1) work correctly with the standard Kepler equation
- Parabolic orbits (e = 1) use Barker's equation correctly
- Hyperbolic orbits (e > 1) use the hyperbolic Kepler equation correctly
- Position calculations work for all orbit types
"""

import pytest
import math
from libephemeris.minor_bodies import (
    solve_kepler_equation,
    solve_kepler_equation_elliptic,
    solve_kepler_equation_hyperbolic,
    solve_barker_equation,
    calc_minor_body_position,
    OrbitalElements,
)


class TestEllipticKeplerEquation:
    """Test the standard Kepler equation for elliptic orbits (e < 1)."""

    @pytest.mark.parametrize(
        "e,M_deg",
        [
            (0.0, 0.0),
            (0.0, 90.0),
            (0.0, 180.0),
            (0.1, 45.0),
            (0.3, 120.0),
            (0.5, 60.0),
            (0.7, 150.0),
            (0.9, 30.0),
        ],
    )
    def test_elliptic_kepler_consistency(self, e, M_deg):
        """Verify E - e*sin(E) = M for elliptic orbits."""
        M = math.radians(M_deg)
        E = solve_kepler_equation_elliptic(M, e)

        # Check the Kepler equation: M = E - e*sin(E)
        M_check = E - e * math.sin(E)
        assert abs(M_check - M) < 1e-8, (
            f"Kepler equation not satisfied: e={e}, M={M_deg}°, "
            f"E={math.degrees(E):.6f}°, residual={M_check - M}"
        )

    def test_elliptic_circular_orbit(self):
        """For circular orbit (e=0), E should equal M."""
        for M_deg in [0, 45, 90, 135, 180, 225, 270, 315]:
            M = math.radians(M_deg)
            E = solve_kepler_equation_elliptic(M, 0.0)
            assert abs(E - M) < 1e-10, (
                f"E should equal M for circular orbit at M={M_deg}°"
            )


class TestHyperbolicKeplerEquation:
    """Test the hyperbolic Kepler equation for hyperbolic orbits (e > 1)."""

    @pytest.mark.parametrize(
        "e,M",
        [
            (1.5, 0.0),
            (1.5, 1.0),
            (1.5, -1.0),
            (2.0, 2.0),
            (2.0, -2.0),
            (3.0, 5.0),
            (1.1, 0.5),
            (5.0, 10.0),
        ],
    )
    def test_hyperbolic_kepler_consistency(self, e, M):
        """Verify e*sinh(H) - H = M for hyperbolic orbits."""
        H = solve_kepler_equation_hyperbolic(M, e)

        # Check the hyperbolic Kepler equation: M = e*sinh(H) - H
        M_check = e * math.sinh(H) - H
        assert abs(M_check - M) < 1e-7, (
            f"Hyperbolic Kepler equation not satisfied: e={e}, M={M}, "
            f"H={H:.6f}, residual={M_check - M}"
        )

    def test_hyperbolic_zero_mean_anomaly(self):
        """For M=0, H should be 0 for any eccentricity."""
        for e in [1.1, 1.5, 2.0, 3.0, 5.0]:
            H = solve_kepler_equation_hyperbolic(0.0, e)
            assert abs(H) < 1e-10, f"H should be 0 for M=0 at e={e}"


class TestBarkerEquation:
    """Test Barker's equation for parabolic orbits (e = 1)."""

    def test_barker_zero_mean_anomaly(self):
        """For M=0, true anomaly should be 0."""
        nu = solve_barker_equation(0.0)
        assert abs(nu) < 1e-10, "True anomaly should be 0 for M=0"

    @pytest.mark.parametrize(
        "M",
        [0.1, 0.5, 1.0, 2.0, -0.1, -0.5, -1.0, -2.0],
    )
    def test_barker_consistency(self, M):
        """Verify Barker's equation: M = tan(ν/2)/2 + tan³(ν/2)/6."""
        nu = solve_barker_equation(M)

        # Check Barker's equation: M = tan(ν/2)/2 + tan³(ν/2)/6
        tan_half_nu = math.tan(nu / 2)
        M_check = tan_half_nu / 2 + tan_half_nu**3 / 6
        assert abs(M_check - M) < 1e-8, (
            f"Barker's equation not satisfied: M={M}, ν={math.degrees(nu):.4f}°, "
            f"residual={M_check - M}"
        )

    def test_barker_symmetry(self):
        """Barker's equation should be antisymmetric: ν(-M) = -ν(M)."""
        for M in [0.1, 0.5, 1.0, 2.0, 5.0]:
            nu_pos = solve_barker_equation(M)
            nu_neg = solve_barker_equation(-M)
            assert abs(nu_pos + nu_neg) < 1e-10, f"Symmetry failed for M={M}"


class TestUnifiedKeplerSolver:
    """Test the unified solve_kepler_equation function dispatches correctly."""

    def test_dispatches_to_elliptic(self):
        """For e < 1, should use elliptic solver."""
        M = math.radians(60.0)
        e = 0.5

        result_unified = solve_kepler_equation(M, e)
        result_elliptic = solve_kepler_equation_elliptic(M, e)

        assert abs(result_unified - result_elliptic) < 1e-12

    def test_dispatches_to_hyperbolic(self):
        """For e > 1, should use hyperbolic solver."""
        M = 1.5
        e = 2.0

        result_unified = solve_kepler_equation(M, e)
        result_hyperbolic = solve_kepler_equation_hyperbolic(M, e)

        assert abs(result_unified - result_hyperbolic) < 1e-12

    def test_dispatches_to_parabolic(self):
        """For e = 1, should use Barker's equation."""
        M = 1.0
        e = 1.0

        result_unified = solve_kepler_equation(M, e)
        result_barker = solve_barker_equation(M)

        assert abs(result_unified - result_barker) < 1e-12

    def test_near_parabolic_uses_barker(self):
        """For e very close to 1, should use Barker's equation."""
        M = 1.0
        e = 1.0 + 1e-12  # Just above 1 but within tolerance

        # This should use Barker's equation due to tolerance
        result = solve_kepler_equation(M, e)

        # Should match Barker's result
        result_barker = solve_barker_equation(M)
        assert abs(result - result_barker) < 1e-10


class TestParabolicPositionCalculation:
    """Test position calculations for parabolic orbits."""

    @pytest.fixture
    def parabolic_comet(self):
        """Create orbital elements for a parabolic comet."""
        return OrbitalElements(
            name="Test Parabolic Comet",
            epoch=2460000.5,  # Arbitrary epoch
            a=1.0,  # Perihelion distance q = 1 AU for parabolic orbits
            e=1.0,  # Exactly parabolic
            i=30.0,
            omega=90.0,
            Omega=45.0,
            M0=0.0,  # At perihelion at epoch
            n=0.5,  # Arbitrary mean motion (deg/day)
        )

    def test_parabolic_at_perihelion(self, parabolic_comet):
        """At perihelion (M=0), distance should equal perihelion distance."""
        x, y, z = calc_minor_body_position(parabolic_comet, parabolic_comet.epoch)

        r = math.sqrt(x**2 + y**2 + z**2)

        # At perihelion, r = q = 1 AU (elements.a stores q for parabolic orbits)
        assert abs(r - parabolic_comet.a) < 1e-6, (
            f"At perihelion, distance should be {parabolic_comet.a} AU, got {r} AU"
        )

    def test_parabolic_position_changes(self, parabolic_comet):
        """Position should change with time for parabolic orbit."""
        pos1 = calc_minor_body_position(parabolic_comet, parabolic_comet.epoch)
        pos2 = calc_minor_body_position(parabolic_comet, parabolic_comet.epoch + 10.0)

        # Positions should be different
        dist = math.sqrt(
            (pos2[0] - pos1[0]) ** 2
            + (pos2[1] - pos1[1]) ** 2
            + (pos2[2] - pos1[2]) ** 2
        )
        assert dist > 0.01, "Position should change over 10 days"


class TestHyperbolicPositionCalculation:
    """Test position calculations for hyperbolic orbits."""

    @pytest.fixture
    def hyperbolic_comet(self):
        """Create orbital elements for a hyperbolic comet (like 'Oumuamua)."""
        return OrbitalElements(
            name="Test Hyperbolic Object",
            epoch=2460000.5,  # Arbitrary epoch
            a=-1.0,  # Negative semi-major axis for hyperbolic (|a| = 1 AU)
            e=1.5,  # Hyperbolic
            i=45.0,
            omega=120.0,
            Omega=60.0,
            M0=0.0,  # At perihelion at epoch
            n=1.0,  # Arbitrary mean motion (deg/day)
        )

    def test_hyperbolic_at_perihelion(self, hyperbolic_comet):
        """At perihelion (M=0), should compute valid position."""
        x, y, z = calc_minor_body_position(hyperbolic_comet, hyperbolic_comet.epoch)

        r = math.sqrt(x**2 + y**2 + z**2)

        # For hyperbolic orbit at perihelion: r_perihelion = a(e - 1)
        # But with negative a: r_perihelion = |a|(e - 1)
        e = hyperbolic_comet.e
        a = hyperbolic_comet.a
        expected_r = abs(a) * (e - 1)

        assert abs(r - expected_r) < 1e-6, (
            f"At perihelion, distance should be {expected_r} AU, got {r} AU"
        )

    def test_hyperbolic_position_changes(self, hyperbolic_comet):
        """Position should change with time for hyperbolic orbit."""
        pos1 = calc_minor_body_position(hyperbolic_comet, hyperbolic_comet.epoch)
        pos2 = calc_minor_body_position(hyperbolic_comet, hyperbolic_comet.epoch + 5.0)

        # Positions should be different
        dist = math.sqrt(
            (pos2[0] - pos1[0]) ** 2
            + (pos2[1] - pos1[1]) ** 2
            + (pos2[2] - pos1[2]) ** 2
        )
        assert dist > 0.01, "Position should change over 5 days"


class TestEllipticPositionCalculation:
    """Test that elliptic orbits still work correctly after changes."""

    @pytest.fixture
    def elliptic_asteroid(self):
        """Create orbital elements for a typical elliptic asteroid."""
        return OrbitalElements(
            name="Test Asteroid",
            epoch=2460000.5,
            a=2.5,  # Semi-major axis
            e=0.2,  # Moderate eccentricity
            i=15.0,
            omega=45.0,
            Omega=90.0,
            M0=180.0,  # At aphelion at epoch
            n=0.25,  # Mean motion
        )

    def test_elliptic_at_aphelion(self, elliptic_asteroid):
        """At aphelion (M=180°), distance should be a(1+e)."""
        x, y, z = calc_minor_body_position(elliptic_asteroid, elliptic_asteroid.epoch)

        r = math.sqrt(x**2 + y**2 + z**2)

        expected_r = elliptic_asteroid.a * (1 + elliptic_asteroid.e)

        assert abs(r - expected_r) < 1e-4, (
            f"At aphelion, distance should be {expected_r} AU, got {r} AU"
        )

    def test_elliptic_orbital_period(self, elliptic_asteroid):
        """After one orbital period, position should return to start (pure Keplerian)."""
        # Period in days
        period_days = 360.0 / elliptic_asteroid.n

        # Use include_perturbations=False to test pure Keplerian behavior
        # With perturbations, omega and Omega precess, so position drifts slightly
        pos1 = calc_minor_body_position(
            elliptic_asteroid, elliptic_asteroid.epoch, include_perturbations=False
        )
        pos2 = calc_minor_body_position(
            elliptic_asteroid,
            elliptic_asteroid.epoch + period_days,
            include_perturbations=False,
        )

        # Positions should be nearly identical
        dist = math.sqrt(
            (pos2[0] - pos1[0]) ** 2
            + (pos2[1] - pos1[1]) ** 2
            + (pos2[2] - pos1[2]) ** 2
        )
        assert dist < 1e-6, f"Position should return after one period, diff={dist}"
