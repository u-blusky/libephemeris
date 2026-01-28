"""
Tests for BesselianElements dataclass in libephemeris.

Tests the data structure that holds Besselian elements and their time
derivatives for solar eclipse calculations.
"""

import pytest
from dataclasses import is_dataclass, fields


class TestBesselianElementsBasicFunctionality:
    """Test basic functionality of BesselianElements dataclass."""

    def test_besselian_elements_exists_in_eclipse_module(self):
        """Test that BesselianElements exists in eclipse module."""
        from libephemeris.eclipse import BesselianElements

        assert BesselianElements is not None

    def test_besselian_elements_exported_from_main_module(self):
        """Test that BesselianElements is exported from main module."""
        from libephemeris import BesselianElements

        assert BesselianElements is not None

    def test_besselian_elements_is_dataclass(self):
        """Test that BesselianElements is a dataclass."""
        from libephemeris import BesselianElements

        assert is_dataclass(BesselianElements)

    def test_besselian_elements_has_required_fields(self):
        """Test that BesselianElements has all required fields."""
        from libephemeris import BesselianElements

        field_names = {f.name for f in fields(BesselianElements)}

        # Check core elements
        assert "t0" in field_names, "Missing reference time field t0"
        assert "x" in field_names, "Missing x coordinate field"
        assert "y" in field_names, "Missing y coordinate field"
        assert "d" in field_names, "Missing declination field d"
        assert "l1" in field_names, "Missing penumbral radius field l1"
        assert "l2" in field_names, "Missing umbral radius field l2"
        assert "mu" in field_names, "Missing hour angle field mu"

        # Check derivatives
        assert "dx_dt" in field_names, "Missing dx/dt derivative"
        assert "dy_dt" in field_names, "Missing dy/dt derivative"
        assert "dd_dt" in field_names, "Missing dd/dt derivative"
        assert "dl1_dt" in field_names, "Missing dl1/dt derivative"
        assert "dl2_dt" in field_names, "Missing dl2/dt derivative"
        assert "dmu_dt" in field_names, "Missing dmu/dt derivative"

    def test_besselian_elements_has_exactly_13_fields(self):
        """Test that BesselianElements has exactly 13 fields (t0 + 6 elements + 6 derivatives)."""
        from libephemeris import BesselianElements

        element_fields = fields(BesselianElements)
        assert len(element_fields) == 13


class TestBesselianElementsInstantiation:
    """Test instantiation of BesselianElements."""

    def test_can_create_instance_with_all_fields(self):
        """Test that we can create a BesselianElements instance."""
        from libephemeris import BesselianElements

        elements = BesselianElements(
            t0=2460408.5,  # Julian Day
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

        assert elements is not None

    def test_fields_are_accessible(self):
        """Test that all fields are accessible after creation."""
        from libephemeris import BesselianElements

        elements = BesselianElements(
            t0=2460408.5,
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

        assert elements.t0 == 2460408.5
        assert elements.x == 0.3145
        assert elements.y == 0.2731
        assert elements.d == 7.5821
        assert elements.l1 == 0.5436
        assert elements.l2 == -0.0047
        assert elements.mu == 89.1234
        assert elements.dx_dt == 0.5123
        assert elements.dy_dt == 0.1456
        assert elements.dd_dt == 0.0012
        assert elements.dl1_dt == -0.0001
        assert elements.dl2_dt == -0.0001
        assert elements.dmu_dt == 15.0041

    def test_missing_field_raises_error(self):
        """Test that missing fields raise TypeError."""
        from libephemeris import BesselianElements

        with pytest.raises(TypeError):
            # Missing required fields should raise TypeError
            BesselianElements(t0=2460408.5, x=0.3)  # type: ignore

    def test_can_store_negative_values(self):
        """Test that negative values (like l2 for total eclipses) are stored correctly."""
        from libephemeris import BesselianElements

        elements = BesselianElements(
            t0=2460408.5,
            x=-0.5,
            y=-0.3,
            d=-10.5,
            l1=0.54,
            l2=-0.01,  # Negative for total eclipse
            mu=-45.0,
            dx_dt=-0.5,
            dy_dt=-0.1,
            dd_dt=-0.001,
            dl1_dt=-0.0001,
            dl2_dt=0.0001,
            dmu_dt=15.0,
        )

        assert elements.x == -0.5
        assert elements.y == -0.3
        assert elements.d == -10.5
        assert elements.l2 == -0.01
        assert elements.mu == -45.0
        assert elements.dx_dt == -0.5


class TestBesselianElementsFieldTypes:
    """Test that field types are correct."""

    def test_all_fields_are_float_type(self):
        """Test that all fields are annotated as float."""
        from libephemeris import BesselianElements

        for field in fields(BesselianElements):
            assert field.type is float, (
                f"Field {field.name} should be float, got {field.type}"
            )


class TestBesselianElementsEquality:
    """Test equality comparison of BesselianElements."""

    def test_equal_elements_are_equal(self):
        """Test that two BesselianElements with same values are equal."""
        from libephemeris import BesselianElements

        elements1 = BesselianElements(
            t0=2460408.5,
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

        elements2 = BesselianElements(
            t0=2460408.5,
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

        assert elements1 == elements2

    def test_different_elements_are_not_equal(self):
        """Test that BesselianElements with different values are not equal."""
        from libephemeris import BesselianElements

        elements1 = BesselianElements(
            t0=2460408.5,
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

        elements2 = BesselianElements(
            t0=2460408.5,
            x=0.5000,  # Different x value
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

        assert elements1 != elements2


class TestBesselianElementsRepresentation:
    """Test string representation of BesselianElements."""

    def test_has_repr(self):
        """Test that BesselianElements has a repr."""
        from libephemeris import BesselianElements

        elements = BesselianElements(
            t0=2460408.5,
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

        repr_str = repr(elements)
        assert "BesselianElements" in repr_str
        assert "t0=" in repr_str
        assert "x=" in repr_str
        assert "dx_dt=" in repr_str


class TestBesselianElementsInterpolation:
    """Test that elements can be used for time interpolation."""

    def test_interpolate_x_coordinate(self):
        """Test that x can be interpolated using dx_dt."""
        from libephemeris import BesselianElements

        elements = BesselianElements(
            t0=2460408.5,
            x=0.3145,
            y=0.2731,
            d=7.5821,
            l1=0.5436,
            l2=-0.0047,
            mu=89.1234,
            dx_dt=0.5123,  # Earth radii per hour
            dy_dt=0.1456,
            dd_dt=0.0012,
            dl1_dt=-0.0001,
            dl2_dt=-0.0001,
            dmu_dt=15.0041,
        )

        # Interpolate 1 hour later
        dt_hours = 1.0
        dt_days = dt_hours / 24.0
        x_interpolated = elements.x + elements.dx_dt * dt_hours

        expected_x = 0.3145 + 0.5123 * 1.0  # = 0.8268
        assert abs(x_interpolated - expected_x) < 1e-10

    def test_interpolate_all_elements(self):
        """Test interpolation for all elements."""
        from libephemeris import BesselianElements

        elements = BesselianElements(
            t0=2460408.5,
            x=0.0,
            y=0.0,
            d=10.0,
            l1=0.5,
            l2=-0.01,
            mu=0.0,
            dx_dt=0.5,
            dy_dt=0.2,
            dd_dt=0.01,
            dl1_dt=-0.001,
            dl2_dt=0.0005,
            dmu_dt=15.0,
        )

        # After 2 hours
        dt = 2.0

        x_new = elements.x + elements.dx_dt * dt
        y_new = elements.y + elements.dy_dt * dt
        d_new = elements.d + elements.dd_dt * dt
        l1_new = elements.l1 + elements.dl1_dt * dt
        l2_new = elements.l2 + elements.dl2_dt * dt
        mu_new = elements.mu + elements.dmu_dt * dt

        assert abs(x_new - 1.0) < 1e-10
        assert abs(y_new - 0.4) < 1e-10
        assert abs(d_new - 10.02) < 1e-10
        assert abs(l1_new - 0.498) < 1e-10
        assert abs(l2_new - (-0.009)) < 1e-10
        assert abs(mu_new - 30.0) < 1e-10
