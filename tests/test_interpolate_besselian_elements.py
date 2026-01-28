"""
Tests for interpolate_besselian_elements function in libephemeris.

Tests the interpolation of Besselian elements to any time during an eclipse
using the elements and derivatives at the reference time.

The interpolation uses first-order Taylor series expansion:
    element(t) = element(t0) + d_element_dt * (t - t0) * 24

References:
- Meeus "Astronomical Algorithms" Chapter 54
- Explanatory Supplement to the Astronomical Almanac
"""

import pytest
import math
from libephemeris import (
    julday,
    BesselianElements,
    interpolate_besselian_elements,
    calc_besselian_x,
    calc_besselian_y,
    calc_besselian_d,
    calc_besselian_l1,
    calc_besselian_l2,
    calc_besselian_mu,
    calc_besselian_dx_dt,
    calc_besselian_dy_dt,
    calc_besselian_dd_dt,
    calc_besselian_dl1_dt,
    calc_besselian_dl2_dt,
    calc_besselian_dmu_dt,
)


class TestInterpolateBesselianElementsBasic:
    """Test basic functionality of interpolate_besselian_elements."""

    def test_function_exists_and_callable(self):
        """Test that interpolate_besselian_elements function exists and is callable."""
        from libephemeris.eclipse import interpolate_besselian_elements

        assert callable(interpolate_besselian_elements)

    def test_function_exported_from_main_module(self):
        """Test that interpolate_besselian_elements is exported from main module."""
        from libephemeris import interpolate_besselian_elements

        assert callable(interpolate_besselian_elements)

    def test_returns_besselian_elements(self):
        """Test that function returns a BesselianElements object."""
        elements = BesselianElements(
            t0=julday(2024, 4, 8, 18.0),
            x=0.3145,
            y=0.2731,
            d=7.5821,
            l1=0.5436,
            l2=-0.0047,
            mu=89.1234,
            dx_dt=0.5123,
            dy_dt=0.1456,
            dd_dt=0.0012,
            dl1_dt=-0.0001,
            dl2_dt=-0.0001,
            dmu_dt=15.0041,
        )
        target_jd = julday(2024, 4, 8, 18.5)
        result = interpolate_besselian_elements(elements, target_jd)
        assert isinstance(result, BesselianElements)

    def test_at_reference_time_returns_same_values(self):
        """Test that interpolating at t0 returns the same element values."""
        t0 = julday(2024, 4, 8, 18.0)
        elements = BesselianElements(
            t0=t0,
            x=0.3145,
            y=0.2731,
            d=7.5821,
            l1=0.5436,
            l2=-0.0047,
            mu=89.1234,
            dx_dt=0.5123,
            dy_dt=0.1456,
            dd_dt=0.0012,
            dl1_dt=-0.0001,
            dl2_dt=-0.0001,
            dmu_dt=15.0041,
        )
        result = interpolate_besselian_elements(elements, t0)

        assert result.t0 == t0
        assert result.x == pytest.approx(elements.x, abs=1e-10)
        assert result.y == pytest.approx(elements.y, abs=1e-10)
        assert result.d == pytest.approx(elements.d, abs=1e-10)
        assert result.l1 == pytest.approx(elements.l1, abs=1e-10)
        assert result.l2 == pytest.approx(elements.l2, abs=1e-10)
        assert result.mu == pytest.approx(elements.mu, abs=1e-10)

    def test_derivatives_are_preserved(self):
        """Test that derivatives are copied to the new BesselianElements."""
        elements = BesselianElements(
            t0=julday(2024, 4, 8, 18.0),
            x=0.3145,
            y=0.2731,
            d=7.5821,
            l1=0.5436,
            l2=-0.0047,
            mu=89.1234,
            dx_dt=0.5123,
            dy_dt=0.1456,
            dd_dt=0.0012,
            dl1_dt=-0.0001,
            dl2_dt=-0.0001,
            dmu_dt=15.0041,
        )
        target_jd = julday(2024, 4, 8, 18.5)
        result = interpolate_besselian_elements(elements, target_jd)

        assert result.dx_dt == elements.dx_dt
        assert result.dy_dt == elements.dy_dt
        assert result.dd_dt == elements.dd_dt
        assert result.dl1_dt == elements.dl1_dt
        assert result.dl2_dt == elements.dl2_dt
        assert result.dmu_dt == elements.dmu_dt

    def test_t0_updated_to_target_time(self):
        """Test that the result's t0 is set to the target time."""
        t0 = julday(2024, 4, 8, 18.0)
        target_jd = julday(2024, 4, 8, 18.5)
        elements = BesselianElements(
            t0=t0,
            x=0.3145,
            y=0.2731,
            d=7.5821,
            l1=0.5436,
            l2=-0.0047,
            mu=89.1234,
            dx_dt=0.5123,
            dy_dt=0.1456,
            dd_dt=0.0012,
            dl1_dt=-0.0001,
            dl2_dt=-0.0001,
            dmu_dt=15.0041,
        )
        result = interpolate_besselian_elements(elements, target_jd)

        assert result.t0 == target_jd


class TestInterpolateBesselianElementsMath:
    """Test mathematical correctness of interpolation."""

    def test_x_interpolation_forward(self):
        """Test x coordinate interpolation forward in time."""
        t0 = julday(2024, 4, 8, 18.0)
        dx_dt = 0.5  # Earth radii per hour
        elements = BesselianElements(
            t0=t0,
            x=0.0,
            y=0.0,
            d=0.0,
            l1=0.5,
            l2=-0.005,
            mu=90.0,
            dx_dt=dx_dt,
            dy_dt=0.0,
            dd_dt=0.0,
            dl1_dt=0.0,
            dl2_dt=0.0,
            dmu_dt=15.0,
        )

        # 1 hour forward
        target_jd = t0 + 1.0 / 24.0
        result = interpolate_besselian_elements(elements, target_jd)

        expected_x = 0.0 + dx_dt * 1.0
        assert result.x == pytest.approx(expected_x, abs=1e-8)

    def test_x_interpolation_backward(self):
        """Test x coordinate interpolation backward in time."""
        t0 = julday(2024, 4, 8, 18.0)
        dx_dt = 0.5  # Earth radii per hour
        elements = BesselianElements(
            t0=t0,
            x=1.0,
            y=0.0,
            d=0.0,
            l1=0.5,
            l2=-0.005,
            mu=90.0,
            dx_dt=dx_dt,
            dy_dt=0.0,
            dd_dt=0.0,
            dl1_dt=0.0,
            dl2_dt=0.0,
            dmu_dt=15.0,
        )

        # 1 hour backward
        target_jd = t0 - 1.0 / 24.0
        result = interpolate_besselian_elements(elements, target_jd)

        expected_x = 1.0 + dx_dt * (-1.0)
        assert result.x == pytest.approx(expected_x, abs=1e-8)

    def test_all_elements_interpolation(self):
        """Test that all elements are interpolated correctly."""
        t0 = julday(2024, 4, 8, 18.0)
        elements = BesselianElements(
            t0=t0,
            x=0.1,
            y=0.2,
            d=7.5,
            l1=0.54,
            l2=-0.005,
            mu=90.0,
            dx_dt=0.5,
            dy_dt=0.15,
            dd_dt=0.001,
            dl1_dt=-0.0001,
            dl2_dt=-0.0001,
            dmu_dt=15.0,
        )

        dt_hours = 0.5  # 30 minutes
        target_jd = t0 + dt_hours / 24.0
        result = interpolate_besselian_elements(elements, target_jd)

        assert result.x == pytest.approx(0.1 + 0.5 * 0.5, abs=1e-8)
        assert result.y == pytest.approx(0.2 + 0.15 * 0.5, abs=1e-8)
        assert result.d == pytest.approx(7.5 + 0.001 * 0.5, abs=1e-8)
        assert result.l1 == pytest.approx(0.54 + (-0.0001) * 0.5, abs=1e-8)
        assert result.l2 == pytest.approx(-0.005 + (-0.0001) * 0.5, abs=1e-8)
        assert result.mu == pytest.approx((90.0 + 15.0 * 0.5) % 360, abs=1e-6)

    def test_mu_normalization_wraparound(self):
        """Test that mu is normalized to [0, 360) after interpolation."""
        t0 = julday(2024, 4, 8, 18.0)
        elements = BesselianElements(
            t0=t0,
            x=0.0,
            y=0.0,
            d=0.0,
            l1=0.5,
            l2=-0.005,
            mu=350.0,
            dx_dt=0.0,
            dy_dt=0.0,
            dd_dt=0.0,
            dl1_dt=0.0,
            dl2_dt=0.0,
            dmu_dt=15.0,
        )

        # 2 hours forward: 350 + 15 * 2 = 380 -> should be 20
        target_jd = t0 + 2.0 / 24.0
        result = interpolate_besselian_elements(elements, target_jd)

        assert result.mu == pytest.approx(20.0, abs=1e-6)

    def test_mu_normalization_large_time(self):
        """Test mu normalization for larger time spans."""
        t0 = julday(2024, 4, 8, 18.0)
        elements = BesselianElements(
            t0=t0,
            x=0.0,
            y=0.0,
            d=0.0,
            l1=0.5,
            l2=-0.005,
            mu=0.0,
            dx_dt=0.0,
            dy_dt=0.0,
            dd_dt=0.0,
            dl1_dt=0.0,
            dl2_dt=0.0,
            dmu_dt=15.0,
        )

        # 25 hours forward: 0 + 15 * 25 = 375 -> 15 degrees
        target_jd = t0 + 25.0 / 24.0
        result = interpolate_besselian_elements(elements, target_jd)

        expected_mu = (0.0 + 15.0 * 25.0) % 360.0
        assert result.mu == pytest.approx(expected_mu, abs=1e-6)


class TestInterpolateBesselianElementsAccuracy:
    """Test interpolation accuracy against actual calculated values."""

    def test_interpolation_accuracy_15_minutes(self):
        """Test accuracy for 15-minute interpolation."""
        jd0 = julday(2024, 4, 8, 18.0)

        # Build BesselianElements from calculated values
        elements = BesselianElements(
            t0=jd0,
            x=calc_besselian_x(jd0),
            y=calc_besselian_y(jd0),
            d=calc_besselian_d(jd0),
            l1=calc_besselian_l1(jd0),
            l2=calc_besselian_l2(jd0),
            mu=calc_besselian_mu(jd0),
            dx_dt=calc_besselian_dx_dt(jd0),
            dy_dt=calc_besselian_dy_dt(jd0),
            dd_dt=calc_besselian_dd_dt(jd0),
            dl1_dt=calc_besselian_dl1_dt(jd0),
            dl2_dt=calc_besselian_dl2_dt(jd0),
            dmu_dt=calc_besselian_dmu_dt(jd0),
        )

        # Interpolate 15 minutes forward
        dt_hours = 0.25
        jd1 = jd0 + dt_hours / 24.0
        interpolated = interpolate_besselian_elements(elements, jd1)

        # Compare with actual calculated values
        x_actual = calc_besselian_x(jd1)
        y_actual = calc_besselian_y(jd1)
        d_actual = calc_besselian_d(jd1)
        l1_actual = calc_besselian_l1(jd1)
        l2_actual = calc_besselian_l2(jd1)
        mu_actual = calc_besselian_mu(jd1)

        # x, y, l1, l2 should be accurate to 0.01 Earth radii
        assert abs(interpolated.x - x_actual) < 0.01, (
            f"x error: {abs(interpolated.x - x_actual):.6f}"
        )
        assert abs(interpolated.y - y_actual) < 0.01, (
            f"y error: {abs(interpolated.y - y_actual):.6f}"
        )
        assert abs(interpolated.l1 - l1_actual) < 0.001, (
            f"l1 error: {abs(interpolated.l1 - l1_actual):.6f}"
        )
        assert abs(interpolated.l2 - l2_actual) < 0.001, (
            f"l2 error: {abs(interpolated.l2 - l2_actual):.6f}"
        )

        # d should be accurate to 0.01 degrees
        assert abs(interpolated.d - d_actual) < 0.01, (
            f"d error: {abs(interpolated.d - d_actual):.6f}"
        )

        # mu needs special handling for wraparound
        mu_error = abs(interpolated.mu - mu_actual)
        if mu_error > 180:
            mu_error = 360 - mu_error
        assert mu_error < 0.1, f"mu error: {mu_error:.6f}"

    def test_interpolation_accuracy_1_hour(self):
        """Test accuracy for 1-hour interpolation."""
        jd0 = julday(2024, 4, 8, 18.0)

        elements = BesselianElements(
            t0=jd0,
            x=calc_besselian_x(jd0),
            y=calc_besselian_y(jd0),
            d=calc_besselian_d(jd0),
            l1=calc_besselian_l1(jd0),
            l2=calc_besselian_l2(jd0),
            mu=calc_besselian_mu(jd0),
            dx_dt=calc_besselian_dx_dt(jd0),
            dy_dt=calc_besselian_dy_dt(jd0),
            dd_dt=calc_besselian_dd_dt(jd0),
            dl1_dt=calc_besselian_dl1_dt(jd0),
            dl2_dt=calc_besselian_dl2_dt(jd0),
            dmu_dt=calc_besselian_dmu_dt(jd0),
        )

        # Interpolate 1 hour forward
        dt_hours = 1.0
        jd1 = jd0 + dt_hours / 24.0
        interpolated = interpolate_besselian_elements(elements, jd1)

        x_actual = calc_besselian_x(jd1)
        y_actual = calc_besselian_y(jd1)

        # Allow slightly larger tolerance for 1 hour (second-order effects)
        assert abs(interpolated.x - x_actual) < 0.05, (
            f"x error: {abs(interpolated.x - x_actual):.6f}"
        )
        assert abs(interpolated.y - y_actual) < 0.05, (
            f"y error: {abs(interpolated.y - y_actual):.6f}"
        )

    def test_interpolation_forward_and_backward(self):
        """Test interpolation works in both time directions."""
        jd0 = julday(2024, 4, 8, 18.0)

        elements = BesselianElements(
            t0=jd0,
            x=calc_besselian_x(jd0),
            y=calc_besselian_y(jd0),
            d=calc_besselian_d(jd0),
            l1=calc_besselian_l1(jd0),
            l2=calc_besselian_l2(jd0),
            mu=calc_besselian_mu(jd0),
            dx_dt=calc_besselian_dx_dt(jd0),
            dy_dt=calc_besselian_dy_dt(jd0),
            dd_dt=calc_besselian_dd_dt(jd0),
            dl1_dt=calc_besselian_dl1_dt(jd0),
            dl2_dt=calc_besselian_dl2_dt(jd0),
            dmu_dt=calc_besselian_dmu_dt(jd0),
        )

        # Forward 15 min
        jd_forward = jd0 + 0.25 / 24.0
        forward = interpolate_besselian_elements(elements, jd_forward)
        x_actual_forward = calc_besselian_x(jd_forward)
        assert abs(forward.x - x_actual_forward) < 0.01

        # Backward 15 min
        jd_backward = jd0 - 0.25 / 24.0
        backward = interpolate_besselian_elements(elements, jd_backward)
        x_actual_backward = calc_besselian_x(jd_backward)
        assert abs(backward.x - x_actual_backward) < 0.01


class TestInterpolateBesselianElementsChaining:
    """Test chaining of interpolations."""

    def test_successive_interpolations(self):
        """Test that successive interpolations are consistent."""
        jd0 = julday(2024, 4, 8, 18.0)

        elements = BesselianElements(
            t0=jd0,
            x=calc_besselian_x(jd0),
            y=calc_besselian_y(jd0),
            d=calc_besselian_d(jd0),
            l1=calc_besselian_l1(jd0),
            l2=calc_besselian_l2(jd0),
            mu=calc_besselian_mu(jd0),
            dx_dt=calc_besselian_dx_dt(jd0),
            dy_dt=calc_besselian_dy_dt(jd0),
            dd_dt=calc_besselian_dd_dt(jd0),
            dl1_dt=calc_besselian_dl1_dt(jd0),
            dl2_dt=calc_besselian_dl2_dt(jd0),
            dmu_dt=calc_besselian_dmu_dt(jd0),
        )

        # Interpolate to +10 min, then to +20 min
        jd1 = jd0 + (10.0 / 60.0) / 24.0
        jd2 = jd0 + (20.0 / 60.0) / 24.0

        step1 = interpolate_besselian_elements(elements, jd1)
        step2 = interpolate_besselian_elements(step1, jd2)

        # Directly to +20 min
        direct = interpolate_besselian_elements(elements, jd2)

        # Both paths should give same result
        # (since derivatives are preserved and formula is linear)
        assert step2.x == pytest.approx(direct.x, abs=1e-10)
        assert step2.y == pytest.approx(direct.y, abs=1e-10)

    def test_interpolate_to_original_time(self):
        """Test interpolating back to original time gives original values."""
        jd0 = julday(2024, 4, 8, 18.0)

        elements = BesselianElements(
            t0=jd0,
            x=0.3145,
            y=0.2731,
            d=7.5821,
            l1=0.5436,
            l2=-0.0047,
            mu=89.1234,
            dx_dt=0.5123,
            dy_dt=0.1456,
            dd_dt=0.0012,
            dl1_dt=-0.0001,
            dl2_dt=-0.0001,
            dmu_dt=15.0041,
        )

        # Interpolate forward then back
        jd1 = jd0 + 0.5 / 24.0
        forward = interpolate_besselian_elements(elements, jd1)
        back = interpolate_besselian_elements(forward, jd0)

        assert back.x == pytest.approx(elements.x, abs=1e-10)
        assert back.y == pytest.approx(elements.y, abs=1e-10)
        assert back.d == pytest.approx(elements.d, abs=1e-10)
        assert back.l1 == pytest.approx(elements.l1, abs=1e-10)
        assert back.l2 == pytest.approx(elements.l2, abs=1e-10)
        assert back.mu == pytest.approx(elements.mu, abs=1e-10)


class TestInterpolateBesselianElementsDifferentEclipses:
    """Test interpolation for different eclipse dates."""

    @pytest.mark.parametrize(
        "year,month,day,hour",
        [
            (2024, 4, 8, 18.0),  # Total eclipse 2024
            (2023, 10, 14, 17.0),  # Annular eclipse 2023
            (2025, 3, 29, 10.0),  # Partial eclipse 2025
        ],
    )
    def test_interpolation_multiple_eclipses(self, year, month, day, hour):
        """Test interpolation works for different eclipse dates."""
        jd0 = julday(year, month, day, hour)

        elements = BesselianElements(
            t0=jd0,
            x=calc_besselian_x(jd0),
            y=calc_besselian_y(jd0),
            d=calc_besselian_d(jd0),
            l1=calc_besselian_l1(jd0),
            l2=calc_besselian_l2(jd0),
            mu=calc_besselian_mu(jd0),
            dx_dt=calc_besselian_dx_dt(jd0),
            dy_dt=calc_besselian_dy_dt(jd0),
            dd_dt=calc_besselian_dd_dt(jd0),
            dl1_dt=calc_besselian_dl1_dt(jd0),
            dl2_dt=calc_besselian_dl2_dt(jd0),
            dmu_dt=calc_besselian_dmu_dt(jd0),
        )

        # Interpolate 15 minutes forward
        jd1 = jd0 + 0.25 / 24.0
        interpolated = interpolate_besselian_elements(elements, jd1)

        # Basic sanity checks
        assert math.isfinite(interpolated.x)
        assert math.isfinite(interpolated.y)
        assert math.isfinite(interpolated.d)
        assert math.isfinite(interpolated.l1)
        assert math.isfinite(interpolated.l2)
        assert 0 <= interpolated.mu < 360
