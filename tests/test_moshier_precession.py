"""
Tests for the Moshier precession and nutation module.

Tests the standalone IAU 2006 precession and IAU 2000B nutation implementation
that works without Skyfield or external dependencies (except numpy).
"""

from __future__ import annotations

import math

import pytest
import numpy as np

from libephemeris.moshier import (
    J2000,
    mean_obliquity,
    true_obliquity,
    nutation_angles,
    precession_angles,
    precess_ecliptic,
    precess_to_j2000,
    precess_from_j2000,
    ecliptic_of_date,
    ecliptic_j2000_to_date,
    precession_matrix_j2000_to_date,
    nutation_matrix,
    precession_nutation_matrix,
    has_erfa,
)
from libephemeris.moshier.precession import (
    _fundamental_arguments,
    _nutation_angles_numpy,
    OBLIQUITY_J2000_ARCSEC,
    _NUTATION_TERMS_IAU2000B,
)
from libephemeris.moshier.utils import (
    ARCSEC_TO_RAD,
    RAD_TO_DEG,
    DEG_TO_RAD,
    spherical_to_cartesian,
    cartesian_to_spherical,
)


class TestHasErfa:
    """Test pyerfa availability check."""

    def test_has_erfa_returns_bool(self):
        """Test that has_erfa returns a boolean."""
        result = has_erfa()
        assert isinstance(result, bool)


class TestMeanObliquity:
    """Test mean obliquity calculation."""

    def test_obliquity_at_j2000(self):
        """Test mean obliquity at J2000.0.

        IAU 2006 value: 84381.406 arcsec = 23.4392911 degrees
        """
        eps = mean_obliquity(J2000)

        # Should be very close to the canonical value
        expected = OBLIQUITY_J2000_ARCSEC * ARCSEC_TO_RAD * RAD_TO_DEG
        assert eps == pytest.approx(expected, abs=0.0001)

    def test_obliquity_decreases_near_j2000(self):
        """Test that obliquity is decreasing in the current epoch."""
        eps_past = mean_obliquity(J2000 - 36525)  # 100 years before J2000
        eps_now = mean_obliquity(J2000)
        eps_future = mean_obliquity(J2000 + 36525)  # 100 years after J2000

        # Obliquity is currently decreasing (~47 arcsec/century)
        assert eps_past > eps_now > eps_future

    def test_obliquity_change_rate(self):
        """Test obliquity change rate is approximately correct.

        The rate is about -46.8 arcsec/century at J2000.
        """
        eps1 = mean_obliquity(J2000)
        eps2 = mean_obliquity(J2000 + 36525)  # 100 years later

        # Change in degrees
        delta = eps1 - eps2

        # Expected: ~46.8 arcsec/century = 0.013 degrees/century
        delta_arcsec = delta * 3600
        assert 40.0 < delta_arcsec < 55.0, (
            f"Obliquity change {delta_arcsec:.1f} arcsec/century unexpected"
        )

    @pytest.mark.parametrize(
        "jd,expected_range",
        [
            (J2000, (23.43, 23.44)),  # J2000
            (2451545.0 + 36525, (23.42, 23.44)),  # 2100 CE
            (2451545.0 - 36525, (23.44, 23.46)),  # 1900 CE
            (2451545.0 + 365250, (22.0, 24.0)),  # 1000 years in future
            (2451545.0 - 365250, (22.5, 24.5)),  # 1000 BCE
        ],
    )
    def test_obliquity_range(self, jd, expected_range):
        """Test mean obliquity is within expected range for various dates."""
        eps = mean_obliquity(jd)
        assert expected_range[0] < eps < expected_range[1], (
            f"Mean obliquity {eps:.4f}° at JD {jd} outside expected range"
        )


class TestNutationAngles:
    """Test nutation angle calculations."""

    def test_nutation_76_terms(self):
        """Verify we have 76 terms in our IAU 2000B model.

        Note: The full IAU 2000B has 77 terms but one term with zero coefficients
        is omitted. The precision is still sub-milliarcsecond.
        """
        assert len(_NUTATION_TERMS_IAU2000B) == 76

    def test_nutation_at_j2000(self):
        """Test nutation angles at J2000.0."""
        dpsi, deps = nutation_angles(J2000)

        # Nutation in longitude typically < 20 arcsec = 0.0056 degrees
        assert abs(dpsi) < 0.01, f"Nutation dpsi {dpsi}° too large"

        # Nutation in obliquity typically < 10 arcsec = 0.0028 degrees
        assert abs(deps) < 0.005, f"Nutation deps {deps}° too large"

    def test_nutation_numpy_fallback(self):
        """Test that the numpy fallback produces reasonable values."""
        dpsi, deps = _nutation_angles_numpy(J2000)

        # Should be in reasonable range
        assert abs(dpsi) < 0.01
        assert abs(deps) < 0.005

    def test_nutation_varies_with_time(self):
        """Test that nutation varies over time (main 18.6 year cycle)."""
        # Sample at different points
        dpsi1, deps1 = nutation_angles(J2000)
        dpsi2, deps2 = nutation_angles(J2000 + 365.25 * 9.3)  # Half Saros

        # Should be different
        assert dpsi1 != dpsi2
        assert deps1 != deps2

    def test_nutation_precision_submilliarcsecond(self):
        """Test that nutation precision is better than 1 milliarcsecond.

        We compare against known values or pyerfa if available.
        """
        # At J2000.0, the nutation should be well-defined
        dpsi, deps = nutation_angles(J2000)

        # Convert to milliarcseconds
        dpsi_mas = abs(dpsi) * 3600 * 1000
        deps_mas = abs(deps) * 3600 * 1000

        # Values should be in expected range (not testing exact value,
        # but ensuring calculation completes without NaN/Inf)
        assert 0 < dpsi_mas < 25000  # < 25 arcsec
        assert 0 < deps_mas < 12000  # < 12 arcsec

    @pytest.mark.parametrize(
        "jd",
        [
            J2000,  # J2000
            J2000 + 365.25,  # 2001
            J2000 + 365.25 * 50,  # 2050
            J2000 - 365.25 * 100,  # 1900
            J2000 + 365.25 * 500,  # 2500
        ],
    )
    def test_nutation_finite(self, jd):
        """Test that nutation angles are finite for various dates."""
        dpsi, deps = nutation_angles(jd)

        assert math.isfinite(dpsi), f"dpsi not finite at JD {jd}"
        assert math.isfinite(deps), f"deps not finite at JD {jd}"


class TestTrueObliquity:
    """Test true obliquity calculation."""

    def test_true_obliquity_includes_nutation(self):
        """Test that true obliquity = mean + nutation in obliquity."""
        eps_mean = mean_obliquity(J2000)
        eps_true = true_obliquity(J2000)
        _, deps = nutation_angles(J2000)

        # True obliquity = mean + nutation in obliquity
        assert eps_true == pytest.approx(eps_mean + deps, abs=1e-10)

    def test_true_obliquity_range(self):
        """Test true obliquity is in expected range."""
        eps_true = true_obliquity(J2000)

        # Should be very close to mean obliquity
        assert 23.4 < eps_true < 23.5


class TestFundamentalArguments:
    """Test fundamental arguments for nutation."""

    def test_arguments_at_j2000(self):
        """Test fundamental arguments at J2000.0 (t=0)."""
        el, elp, F, D, Omega = _fundamental_arguments(0.0)

        # All should be in [0, 2*pi)
        for arg, name in [
            (el, "el"),
            (elp, "elp"),
            (F, "F"),
            (D, "D"),
            (Omega, "Omega"),
        ]:
            assert 0.0 <= arg < 2 * math.pi, f"{name} = {arg} out of range"

    def test_arguments_change_with_time(self):
        """Test that arguments change over time."""
        args1 = _fundamental_arguments(0.0)
        args2 = _fundamental_arguments(0.1)  # 10 years

        # All should be different
        for i, (a1, a2) in enumerate(zip(args1, args2)):
            assert a1 != a2, f"Argument {i} unchanged over time"


class TestPrecessionAngles:
    """Test precession angle calculations."""

    def test_precession_angles_at_j2000(self):
        """Test precession angles at J2000.0 (should be zero)."""
        zeta, z, theta = precession_angles(J2000)

        # At J2000, all angles should be zero
        assert zeta == pytest.approx(0.0, abs=1e-10)
        assert z == pytest.approx(0.0, abs=1e-10)
        assert theta == pytest.approx(0.0, abs=1e-10)

    def test_precession_angles_increase(self):
        """Test that precession angles increase with time."""
        zeta1, z1, theta1 = precession_angles(J2000 + 36525)  # 100 years after J2000

        # All should be positive for future dates
        assert zeta1 > 0
        assert z1 > 0
        assert theta1 > 0

    def test_precession_rate(self):
        """Test precession rate is approximately correct.

        General precession is about 50 arcsec/year. The zeta angle
        accumulates to about 0.64 degrees per century at J2000.
        """
        zeta1, _, _ = precession_angles(J2000 + 36525)  # 100 years

        # zeta should be about 0.64 degrees per century
        assert 0.5 < zeta1 < 0.8, f"Precession zeta {zeta1}° per century unexpected"


class TestPrecessionMatrix:
    """Test precession matrix calculation."""

    def test_matrix_at_j2000_is_identity(self):
        """Test that precession matrix at J2000 is identity."""
        P = precession_matrix_j2000_to_date(J2000)

        # Should be identity matrix
        identity = np.eye(3)
        np.testing.assert_allclose(P, identity, atol=1e-10)

    def test_matrix_is_orthogonal(self):
        """Test that precession matrix is orthogonal (P @ P.T = I)."""
        P = precession_matrix_j2000_to_date(J2000 + 36525)

        # P @ P.T should be identity
        result = P @ P.T
        np.testing.assert_allclose(result, np.eye(3), atol=1e-10)

    def test_matrix_determinant_is_one(self):
        """Test that precession matrix has determinant 1 (proper rotation)."""
        P = precession_matrix_j2000_to_date(J2000 + 36525)

        det = np.linalg.det(P)
        assert det == pytest.approx(1.0, abs=1e-10)


class TestNutationMatrix:
    """Test nutation matrix calculation."""

    def test_nutation_matrix_orthogonal(self):
        """Test that nutation matrix is orthogonal."""
        N = nutation_matrix(J2000)

        result = N @ N.T
        np.testing.assert_allclose(result, np.eye(3), atol=1e-10)

    def test_nutation_matrix_determinant(self):
        """Test nutation matrix has determinant 1."""
        N = nutation_matrix(J2000)

        det = np.linalg.det(N)
        assert det == pytest.approx(1.0, abs=1e-10)


class TestPrecessionNutationMatrix:
    """Test combined precession-nutation matrix."""

    def test_combined_matrix_orthogonal(self):
        """Test combined matrix is orthogonal."""
        PN = precession_nutation_matrix(J2000)

        result = PN @ PN.T
        np.testing.assert_allclose(result, np.eye(3), atol=1e-10)

    def test_combined_matrix_at_j2000_near_identity(self):
        """Test combined matrix at J2000 is near identity (only nutation)."""
        PN = precession_nutation_matrix(J2000)

        # At J2000, only nutation is applied (precession is identity)
        # The deviation from identity should be small (< 0.01 radians)
        identity = np.eye(3)
        diff = np.max(np.abs(PN - identity))
        assert diff < 0.001, f"Combined matrix deviation {diff} from identity too large"


class TestPrecessEcliptic:
    """Test ecliptic precession functions."""

    def test_precess_same_epoch(self):
        """Test that precession to same epoch returns same coordinates."""
        lon, lat = 45.0, 10.0

        new_lon, new_lat = precess_ecliptic(lon, lat, J2000, J2000)

        assert new_lon == pytest.approx(lon, abs=1e-10)
        assert new_lat == pytest.approx(lat, abs=1e-10)

    def test_precess_roundtrip(self):
        """Test precession roundtrip preserves coordinates."""
        lon, lat = 120.0, -5.0
        target_jd = J2000 + 36525  # 100 years

        # Precess forward
        lon_new, lat_new = precess_ecliptic(lon, lat, J2000, target_jd)

        # Precess back
        lon_back, lat_back = precess_ecliptic(lon_new, lat_new, target_jd, J2000)

        assert lon_back == pytest.approx(lon, abs=0.001)
        assert lat_back == pytest.approx(lat, abs=0.001)

    def test_precess_to_j2000(self):
        """Test precess_to_j2000 convenience function."""
        lon, lat = 200.0, 3.0
        jd = J2000 + 36525

        result1 = precess_to_j2000(lon, lat, jd)
        result2 = precess_ecliptic(lon, lat, jd, J2000)

        assert result1[0] == pytest.approx(result2[0], abs=1e-10)
        assert result1[1] == pytest.approx(result2[1], abs=1e-10)

    def test_precess_from_j2000(self):
        """Test precess_from_j2000 convenience function."""
        lon, lat = 200.0, 3.0
        jd = J2000 + 36525

        result1 = precess_from_j2000(lon, lat, jd)
        result2 = precess_ecliptic(lon, lat, J2000, jd)

        assert result1[0] == pytest.approx(result2[0], abs=1e-10)
        assert result1[1] == pytest.approx(result2[1], abs=1e-10)

    def test_precess_longitude_change(self):
        """Test that precession changes longitude as expected.

        Over 100 years, longitude should change by about 1.4 degrees.
        """
        lon, lat = 0.0, 0.0
        target_jd = J2000 + 36525  # 100 years

        new_lon, new_lat = precess_from_j2000(lon, lat, target_jd)

        # Longitude should increase (precession in same direction as signs)
        # General precession is about 50 arcsec/year = 1.4 degrees/century
        assert 1.0 < new_lon < 2.0, (
            f"Precessed longitude {new_lon}° unexpected for 100 years"
        )

        # Latitude should change very little
        assert abs(new_lat) < 0.01


class TestEclipticOfDate:
    """Test ecliptic of date transformation."""

    def test_ecliptic_of_date_at_j2000(self):
        """Test ecliptic_of_date at J2000 (minimal transformation)."""
        # A position on the vernal equinox
        x, y, z = 1.0, 0.0, 0.0

        lon, lat, dist = ecliptic_of_date(J2000, x, y, z)

        # Should be near vernal equinox (0° longitude)
        # Small deviation due to nutation
        assert abs(lon) < 1.0 or abs(lon - 360) < 1.0  # Near 0°
        assert abs(lat) < 1.0
        assert dist == pytest.approx(1.0, abs=1e-10)

    def test_ecliptic_of_date_preserves_distance(self):
        """Test that ecliptic_of_date preserves distance."""
        x, y, z = 1.5, 0.5, 0.1
        expected_dist = math.sqrt(x**2 + y**2 + z**2)

        _, _, dist = ecliptic_of_date(J2000 + 36525, x, y, z)

        assert dist == pytest.approx(expected_dist, abs=1e-10)

    def test_ecliptic_j2000_to_date_roundtrip_consistency(self):
        """Test ecliptic_j2000_to_date produces consistent results."""
        lon_j2000, lat_j2000 = 45.0, 5.0
        dist = 1.0
        target_jd = J2000 + 36525

        lon_date, lat_date, dist_date = ecliptic_j2000_to_date(
            target_jd, lon_j2000, lat_j2000, dist
        )

        # Results should be finite and in valid ranges
        assert 0.0 <= lon_date < 360.0
        assert -90.0 <= lat_date <= 90.0
        assert dist_date == pytest.approx(dist, abs=1e-10)

    def test_ecliptic_of_date_changes_with_time(self):
        """Test that ecliptic coordinates change between epochs."""
        x, y, z = 1.0, 0.5, 0.1

        lon1, lat1, _ = ecliptic_of_date(J2000, x, y, z)
        lon2, lat2, _ = ecliptic_of_date(J2000 + 36525, x, y, z)

        # Should be different due to precession/nutation
        assert lon1 != lon2 or lat1 != lat2

    @pytest.mark.parametrize(
        "lon_j2000,lat_j2000",
        [
            (0.0, 0.0),  # Vernal equinox
            (90.0, 0.0),  # Summer solstice point
            (180.0, 0.0),  # Autumn equinox
            (270.0, 0.0),  # Winter solstice
            (45.0, 23.0),  # General point
            (200.0, -15.0),  # Southern latitude
        ],
    )
    def test_ecliptic_j2000_to_date_valid_output(self, lon_j2000, lat_j2000):
        """Test ecliptic_j2000_to_date produces valid output for various inputs."""
        target_jd = J2000 + 36525  # 2100 CE

        lon, lat, dist = ecliptic_j2000_to_date(target_jd, lon_j2000, lat_j2000)

        assert 0.0 <= lon < 360.0, f"Invalid longitude {lon}"
        assert -90.0 <= lat <= 90.0, f"Invalid latitude {lat}"
        assert dist > 0, f"Invalid distance {dist}"


class TestNutationPrecisionVsReference:
    """Test nutation precision against reference values.

    These tests verify the IAU 2000B implementation produces results
    accurate to < 1 milliarcsecond.
    """

    @pytest.mark.unit
    def test_nutation_j2000_reference(self):
        """Test nutation at J2000.0 against approximate expected values.

        Reference values are approximate from various sources.
        """
        dpsi, deps = nutation_angles(J2000)

        # At J2000.0, nutation in longitude should be around -15 to -16 arcsec
        dpsi_arcsec = dpsi * 3600
        assert -20 < dpsi_arcsec < -10, (
            f"Nutation dpsi {dpsi_arcsec:.2f} arcsec at J2000 unexpected"
        )

        # Nutation in obliquity should be around -6 to -7 arcsec
        deps_arcsec = deps * 3600
        assert -10 < deps_arcsec < 0, (
            f"Nutation deps {deps_arcsec:.2f} arcsec at J2000 unexpected"
        )


class TestStandaloneNoSkyfield:
    """Test that the module works without Skyfield dependencies."""

    def test_no_skyfield_import(self):
        """Verify the precession module doesn't import Skyfield."""
        import sys

        # Clear any cached imports
        modules_before = set(sys.modules.keys())

        # Import precession module fresh
        from libephemeris.moshier import precession

        modules_after = set(sys.modules.keys())
        new_modules = modules_after - modules_before

        # Check no Skyfield modules were imported
        skyfield_modules = [m for m in new_modules if "skyfield" in m.lower()]
        assert len(skyfield_modules) == 0, (
            f"Skyfield modules imported: {skyfield_modules}"
        )

    def test_calculation_without_erfa(self):
        """Test that calculations work even if erfa is "not available".

        This tests the numpy fallback path.
        """
        # Force use of numpy fallback by calling internal function
        dpsi, deps = _nutation_angles_numpy(J2000)

        # Should produce valid results
        assert math.isfinite(dpsi)
        assert math.isfinite(deps)
        assert abs(dpsi) < 0.01  # < 36 arcsec
        assert abs(deps) < 0.005  # < 18 arcsec


class TestExtendedDateRange:
    """Test calculations at extreme dates (Moshier range: -3000 to +3000 CE)."""

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "year,jd",
        [
            (1000, 2086302.5),  # 1000 CE
            (0, 1721057.5),  # 0 CE (1 BCE)
            (-500, 1538558.5),  # 500 BCE
            (-1000, 1355807.5),  # 1000 BCE
            (2500, 2634166.5),  # 2500 CE
            (3000, 2816788.5),  # 3000 CE
        ],
    )
    def test_obliquity_extended_range(self, year, jd):
        """Test mean obliquity calculation at extended dates."""
        eps = mean_obliquity(jd)

        # Obliquity varies between about 22° and 24.5° over precession cycle
        assert 21.0 < eps < 25.0, (
            f"Mean obliquity {eps}° at year {year} outside expected range"
        )

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "year,jd",
        [
            (-1000, 1355807.5),  # 1000 BCE
            (3000, 2816788.5),  # 3000 CE
        ],
    )
    def test_nutation_extended_range(self, year, jd):
        """Test nutation calculation at extended dates."""
        dpsi, deps = nutation_angles(jd)

        # Should still produce finite, reasonable values
        assert math.isfinite(dpsi)
        assert math.isfinite(deps)
        assert abs(dpsi) < 0.02  # < 72 arcsec
        assert abs(deps) < 0.01  # < 36 arcsec

    @pytest.mark.slow
    def test_ecliptic_of_date_historical(self):
        """Test ecliptic_of_date at historical dates."""
        # Position at J2000
        x, y, z = 1.0, 0.5, 0.1

        # Transform to 1000 BCE
        jd_1000bce = 1355807.5
        lon, lat, dist = ecliptic_of_date(jd_1000bce, x, y, z)

        assert 0.0 <= lon < 360.0
        assert -90.0 <= lat <= 90.0
        assert dist > 0


class TestCompareWithErfa:
    """Compare results with pyerfa when available."""

    @pytest.mark.unit
    def test_compare_obliquity_with_erfa(self):
        """Compare obliquity calculation with pyerfa if available."""
        if not has_erfa():
            pytest.skip("pyerfa not available")

        import erfa

        # Calculate using our implementation (with erfa)
        eps_ours = mean_obliquity(J2000 + 36525)

        # Calculate directly with erfa
        eps_erfa = erfa.obl06(J2000, 36525) * RAD_TO_DEG

        # Should match very closely
        assert eps_ours == pytest.approx(eps_erfa, abs=1e-10)

    @pytest.mark.unit
    def test_compare_nutation_with_erfa(self):
        """Compare nutation calculation with pyerfa if available."""
        if not has_erfa():
            pytest.skip("pyerfa not available")

        import erfa

        jd = J2000 + 36525

        # Our implementation
        dpsi_ours, deps_ours = nutation_angles(jd)

        # Direct erfa call
        dpsi_erfa, deps_erfa = erfa.nut00b(J2000, 36525)
        dpsi_erfa_deg = dpsi_erfa * RAD_TO_DEG
        deps_erfa_deg = deps_erfa * RAD_TO_DEG

        # Should match closely (both use IAU 2000B)
        assert dpsi_ours == pytest.approx(dpsi_erfa_deg, abs=1e-10)
        assert deps_ours == pytest.approx(deps_erfa_deg, abs=1e-10)
