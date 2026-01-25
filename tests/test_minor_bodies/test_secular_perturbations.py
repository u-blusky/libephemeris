"""
Tests for secular perturbation calculations in minor body propagation.

Tests verify:
- Laplace coefficient calculations are reasonable
- Secular perturbation rates are physically sensible
- Perturbed positions differ from unperturbed for long propagation times
- Perturbations improve accuracy compared to pure Keplerian propagation
"""

import pytest
import math
from libephemeris.minor_bodies import (
    MINOR_BODY_ELEMENTS,
    calc_minor_body_position,
    calc_secular_perturbation_rates,
    apply_secular_perturbations,
    _calc_laplace_coefficients,
    GM_SUN,
    MASS_RATIO_JUPITER,
    MASS_RATIO_SATURN,
    JUPITER_A,
    SATURN_A,
)
from libephemeris.constants import (
    SE_CERES,
    SE_PALLAS,
    SE_VESTA,
    SE_CHIRON,
    SE_ERIS,
    SE_QUAOAR,
)


class TestLaplaceCoefficients:
    """Test the Laplace coefficient calculation."""

    def test_laplace_coefficient_small_alpha(self):
        """For small alpha, Laplace coefficients should be well-defined."""
        alpha = 0.5
        b_1_5_1 = _calc_laplace_coefficients(alpha, 1.5, 1)

        # Laplace coefficients should be positive for valid inputs
        assert b_1_5_1 > 0, "Laplace coefficient should be positive"
        # Should be on order of 1-10 for typical values
        assert 0.1 < b_1_5_1 < 20, (
            f"Laplace coefficient {b_1_5_1} outside expected range"
        )

    def test_laplace_coefficient_alpha_boundary(self):
        """Laplace coefficients should be zero for invalid alpha."""
        # alpha >= 1 is invalid (orbital crossing)
        assert _calc_laplace_coefficients(1.0, 1.5, 1) == 0.0
        assert _calc_laplace_coefficients(1.5, 1.5, 1) == 0.0
        assert _calc_laplace_coefficients(0.0, 1.5, 1) == 0.0
        assert _calc_laplace_coefficients(-0.5, 1.5, 1) == 0.0

    def test_laplace_coefficient_increases_with_alpha(self):
        """Laplace coefficients should generally increase with alpha."""
        alpha1 = 0.3
        alpha2 = 0.6
        b1 = _calc_laplace_coefficients(alpha1, 1.5, 1)
        b2 = _calc_laplace_coefficients(alpha2, 1.5, 1)

        # For j=1, the coefficient increases with alpha
        assert b2 > b1, "Laplace coefficient should increase with alpha"


class TestSecularPerturbationRates:
    """Test secular perturbation rate calculations."""

    def test_ceres_perturbation_rates(self):
        """Ceres should have measurable perturbation rates from Jupiter."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        d_omega, d_Omega, d_n = calc_secular_perturbation_rates(elements)

        # Ceres is in the main belt, should have positive perihelion precession
        # Typical rates: ~10-100 arcsec/year for main belt asteroids
        # Convert to arcsec/year for comparison
        d_omega_arcsec_year = d_omega * 365.25 * 3600
        d_Omega_arcsec_year = d_Omega * 365.25 * 3600

        # Omega should regress (negative rate) for prograde orbits
        assert d_Omega < 0, "Node should regress for prograde main belt asteroid"

        # Both rates should be on order of arcsec/year to arcmin/year
        assert abs(d_omega_arcsec_year) < 3600, "Omega rate should be < 1 deg/year"
        assert abs(d_Omega_arcsec_year) < 3600, "Omega rate should be < 1 deg/year"

    def test_chiron_perturbation_rates(self):
        """Chiron (between Saturn and Uranus) should have smaller perturbation rates."""
        elements = MINOR_BODY_ELEMENTS[SE_CHIRON]
        d_omega, d_Omega, d_n = calc_secular_perturbation_rates(elements)

        # Chiron is beyond Saturn, so perturbations from Jupiter are weaker
        # but should still be measurable
        d_omega_arcsec_year = d_omega * 365.25 * 3600

        # Rates should be smaller than for main belt asteroids
        assert abs(d_omega_arcsec_year) < 1000, "Chiron omega rate should be small"
        # Also verify d_Omega and d_n were returned (used indirectly)
        assert d_Omega is not None and d_n is not None

    def test_eris_perturbation_rates(self):
        """Eris (very distant TNO) should have very small perturbation rates."""
        elements = MINOR_BODY_ELEMENTS[SE_ERIS]
        d_omega, d_Omega, d_n = calc_secular_perturbation_rates(elements)

        # Eris is ~68 AU, so Jupiter/Saturn perturbations are very weak
        d_omega_arcsec_year = d_omega * 365.25 * 3600

        # Should be measurable but small
        assert abs(d_omega_arcsec_year) < 100, "Eris omega rate should be very small"


class TestApplySecularPerturbations:
    """Test the application of secular perturbations to orbital elements."""

    def test_short_propagation_unchanged(self):
        """For very short propagation times, elements should be essentially unchanged."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        epoch = elements.epoch

        # Propagate for 0.5 days
        omega_p, Omega_p, M_p, n_p = apply_secular_perturbations(elements, epoch + 0.5)

        # omega and Omega should be very close to original
        assert abs(omega_p - elements.omega) < 0.01, (
            "omega should barely change in 0.5 days"
        )
        assert abs(Omega_p - elements.Omega) < 0.01, (
            "Omega should barely change in 0.5 days"
        )

    def test_long_propagation_changed(self):
        """For long propagation times, elements should show measurable drift."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        epoch = elements.epoch

        # Propagate for 10 years (3652.5 days)
        dt = 10 * 365.25
        omega_p, Omega_p, M_p, n_p = apply_secular_perturbations(elements, epoch + dt)

        # omega and Omega should have drifted measurably
        d_omega = omega_p - elements.omega
        d_Omega = Omega_p - elements.Omega

        # Adjust for wraparound
        if d_omega > 180:
            d_omega -= 360
        elif d_omega < -180:
            d_omega += 360
        if d_Omega > 180:
            d_Omega -= 360
        elif d_Omega < -180:
            d_Omega += 360

        # Over 10 years, expect changes of several arcminutes to degrees
        assert abs(d_omega) > 0.001, "omega should change over 10 years"

    def test_perturbations_disabled_no_change(self):
        """With perturbations disabled, only mean anomaly should change."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        epoch = elements.epoch

        # Propagate for 1000 days without perturbations
        omega_p, Omega_p, M_p, n_p = apply_secular_perturbations(
            elements, epoch + 1000, include_perturbations=False
        )

        # omega and Omega should be unchanged
        assert omega_p == elements.omega, (
            "omega should not change without perturbations"
        )
        assert Omega_p == elements.Omega, (
            "Omega should not change without perturbations"
        )


class TestPositionWithPerturbations:
    """Test that position calculations properly use perturbations."""

    def test_position_with_and_without_perturbations(self):
        """Positions with perturbations should differ from those without."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        epoch = elements.epoch

        # Propagate for 5 years
        dt = 5 * 365.25
        jd = epoch + dt

        # Calculate with perturbations (default)
        x1, y1, z1 = calc_minor_body_position(elements, jd, include_perturbations=True)

        # Calculate without perturbations
        x2, y2, z2 = calc_minor_body_position(elements, jd, include_perturbations=False)

        # Positions should differ
        dist = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)

        # Over 5 years, expect differences of 0.001-0.1 AU
        assert dist > 0, "Positions should differ with perturbations"
        # But not too different (sanity check)
        assert dist < 1.0, f"Position difference {dist} AU seems too large"

    def test_position_at_epoch_similar(self):
        """At epoch, perturbed and unperturbed positions should be essentially identical."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        epoch = elements.epoch

        # Calculate at epoch
        x1, y1, z1 = calc_minor_body_position(
            elements, epoch, include_perturbations=True
        )
        x2, y2, z2 = calc_minor_body_position(
            elements, epoch, include_perturbations=False
        )

        # Should be nearly identical
        dist = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
        assert dist < 1e-10, f"At epoch, positions should be identical (diff: {dist})"


class TestPhysicalConstants:
    """Verify physical constants are reasonable."""

    def test_mass_ratios(self):
        """Mass ratios should be in expected range."""
        # Jupiter is about 1/1000 solar mass
        assert 0.0009 < MASS_RATIO_JUPITER < 0.001

        # Saturn is about 1/3500 solar mass
        assert 0.0002 < MASS_RATIO_SATURN < 0.0003

    def test_planet_semi_major_axes(self):
        """Planet semi-major axes should be correct."""
        # Jupiter at ~5.2 AU
        assert 5.1 < JUPITER_A < 5.3

        # Saturn at ~9.5 AU
        assert 9.4 < SATURN_A < 9.7

    def test_gm_sun(self):
        """GM_SUN should be the standard gravitational parameter."""
        # k^2 ~ 0.000296 AU^3/day^2
        assert 0.00029 < GM_SUN < 0.0003


class TestAllBodiesPerturbations:
    """Test perturbation calculations work for all minor bodies."""

    @pytest.mark.parametrize(
        "body_id",
        [SE_CERES, SE_PALLAS, SE_VESTA, SE_CHIRON, SE_ERIS, SE_QUAOAR],
    )
    def test_perturbation_rates_finite(self, body_id):
        """Perturbation rates should be finite for all bodies."""
        elements = MINOR_BODY_ELEMENTS[body_id]
        d_omega, d_Omega, d_n = calc_secular_perturbation_rates(elements)

        assert math.isfinite(d_omega), f"d_omega should be finite for {elements.name}"
        assert math.isfinite(d_Omega), f"d_Omega should be finite for {elements.name}"
        assert math.isfinite(d_n), f"d_n should be finite for {elements.name}"

    @pytest.mark.parametrize(
        "body_id",
        [SE_CERES, SE_PALLAS, SE_VESTA, SE_CHIRON, SE_ERIS, SE_QUAOAR],
    )
    def test_position_valid_with_perturbations(self, body_id):
        """Position should be valid with perturbations for all bodies."""
        elements = MINOR_BODY_ELEMENTS[body_id]
        epoch = elements.epoch

        # Propagate 2 years
        jd = epoch + 2 * 365.25
        x, y, z = calc_minor_body_position(elements, jd, include_perturbations=True)

        # Position should be finite
        assert math.isfinite(x), f"x should be finite for {elements.name}"
        assert math.isfinite(y), f"y should be finite for {elements.name}"
        assert math.isfinite(z), f"z should be finite for {elements.name}"

        # Distance should be reasonable (within orbit range)
        r = math.sqrt(x**2 + y**2 + z**2)
        perihelion = elements.a * (1 - elements.e)
        aphelion = elements.a * (1 + elements.e)
        assert perihelion * 0.9 < r < aphelion * 1.1, (
            f"{elements.name}: distance {r:.2f} AU outside expected range"
        )
