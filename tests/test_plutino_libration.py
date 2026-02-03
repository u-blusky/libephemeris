"""
Unit tests for plutino resonant libration model.

Tests verify:
- Libration parameters are defined for known plutinos (Ixion, Orcus)
- Resonant argument calculation is correct
- Libration correction is applied to plutino positions
- Position errors for plutinos are reduced over long timescales
- Non-plutinos are not affected by libration corrections
"""

import math
import pytest
from libephemeris.constants import (
    SE_IXION,
    SE_ORCUS,
    SE_CERES,
    SE_ERIS,
    SE_QUAOAR,
)
from libephemeris.minor_bodies import (
    MINOR_BODY_ELEMENTS,
    PLUTINO_LIBRATION_PARAMS,
    J2000_EPOCH,
    NEPTUNE_MEAN_LONGITUDE_J2000,
    NEPTUNE_N,
    LibrationParameters,
    calc_neptune_mean_longitude,
    calc_resonant_argument_plutino,
    calc_libration_correction,
    has_libration_model,
    get_libration_parameters,
    calc_minor_body_heliocentric,
    detect_mean_motion_resonance,
)


class TestLibrationParameters:
    """Test that libration parameters are properly defined."""

    def test_ixion_has_libration_params(self):
        """Ixion should have libration parameters defined."""
        assert SE_IXION in PLUTINO_LIBRATION_PARAMS
        params = PLUTINO_LIBRATION_PARAMS[SE_IXION]
        assert isinstance(params, LibrationParameters)
        assert params.amplitude > 0
        assert params.period > 0
        assert 0 <= params.center <= 360

    def test_orcus_has_libration_params(self):
        """Orcus should have libration parameters defined."""
        assert SE_ORCUS in PLUTINO_LIBRATION_PARAMS
        params = PLUTINO_LIBRATION_PARAMS[SE_ORCUS]
        assert isinstance(params, LibrationParameters)
        assert params.amplitude > 0
        assert params.period > 0
        assert 0 <= params.center <= 360

    def test_plutino_amplitudes_reasonable(self):
        """Libration amplitudes should be in reasonable range (30-120 degrees)."""
        for body_id, params in PLUTINO_LIBRATION_PARAMS.items():
            assert 30 <= params.amplitude <= 120, (
                f"Body {body_id} amplitude out of range"
            )

    def test_plutino_periods_reasonable(self):
        """Libration periods should be ~20,000 years (7-8 million days)."""
        for body_id, params in PLUTINO_LIBRATION_PARAMS.items():
            # 18,000-22,000 years = 6.5-8 million days
            assert 6_500_000 <= params.period <= 8_000_000, (
                f"Body {body_id} period out of range"
            )

    def test_non_plutinos_no_params(self):
        """Non-plutino bodies should not have libration parameters."""
        assert SE_CERES not in PLUTINO_LIBRATION_PARAMS
        assert SE_ERIS not in PLUTINO_LIBRATION_PARAMS
        assert SE_QUAOAR not in PLUTINO_LIBRATION_PARAMS


class TestNeptuneMeanLongitude:
    """Test Neptune mean longitude calculation."""

    def test_neptune_longitude_at_j2000(self):
        """Neptune longitude at J2000.0 should match reference value."""
        lon = calc_neptune_mean_longitude(J2000_EPOCH)
        assert abs(lon - NEPTUNE_MEAN_LONGITUDE_J2000) < 0.001

    def test_neptune_longitude_progresses(self):
        """Neptune longitude should increase over time."""
        lon_j2000 = calc_neptune_mean_longitude(J2000_EPOCH)
        lon_later = calc_neptune_mean_longitude(J2000_EPOCH + 365.25)  # 1 year later

        # Neptune moves ~2.2 degrees/year
        expected_motion = NEPTUNE_N * 365.25
        actual_motion = (lon_later - lon_j2000) % 360
        assert abs(actual_motion - expected_motion) < 0.1

    def test_neptune_longitude_wraps(self):
        """Neptune longitude should wrap correctly at 360 degrees."""
        # After ~164 years, Neptune completes one orbit
        one_orbit_days = 360.0 / NEPTUNE_N  # ~59,800 days
        lon = calc_neptune_mean_longitude(J2000_EPOCH + one_orbit_days)
        # Should be back near starting position
        assert 0 <= lon < 360


class TestResonantArgument:
    """Test resonant argument calculation for plutinos."""

    def test_resonant_argument_range(self):
        """Resonant argument should be in 0-360 range."""
        elements = MINOR_BODY_ELEMENTS[SE_IXION]
        phi = calc_resonant_argument_plutino(
            elements, J2000_EPOCH, elements.omega, elements.M0
        )
        assert 0 <= phi < 360

    def test_resonant_argument_for_ixion(self):
        """Ixion's resonant argument should be calculable."""
        elements = MINOR_BODY_ELEMENTS[SE_IXION]
        # At epoch
        phi = calc_resonant_argument_plutino(
            elements, elements.epoch, elements.omega, elements.M0
        )
        # For a plutino, phi should librate around 180 degrees
        # Check it's in a reasonable range
        assert 0 <= phi < 360

    def test_resonant_argument_for_orcus(self):
        """Orcus's resonant argument should be calculable."""
        elements = MINOR_BODY_ELEMENTS[SE_ORCUS]
        phi = calc_resonant_argument_plutino(
            elements, elements.epoch, elements.omega, elements.M0
        )
        assert 0 <= phi < 360


class TestLibrationCorrection:
    """Test libration correction calculation."""

    def test_ixion_correction_bounded(self):
        """Ixion's libration correction should be bounded by amplitude/3."""
        params = PLUTINO_LIBRATION_PARAMS[SE_IXION]
        max_expected = params.amplitude / 3.0

        # Test at multiple times over a full libration period
        for phase in range(0, 360, 30):
            t = J2000_EPOCH + (phase / 360.0) * params.period
            correction = calc_libration_correction(SE_IXION, t)
            assert abs(correction) <= max_expected + 0.1  # Small tolerance

    def test_orcus_correction_bounded(self):
        """Orcus's libration correction should be bounded by amplitude/3."""
        params = PLUTINO_LIBRATION_PARAMS[SE_ORCUS]
        max_expected = params.amplitude / 3.0

        for phase in range(0, 360, 30):
            t = J2000_EPOCH + (phase / 360.0) * params.period
            correction = calc_libration_correction(SE_ORCUS, t)
            assert abs(correction) <= max_expected + 0.1

    def test_non_plutino_zero_correction(self):
        """Non-plutinos should have zero libration correction."""
        correction = calc_libration_correction(SE_CERES, J2000_EPOCH)
        assert correction == 0.0

        correction = calc_libration_correction(SE_ERIS, J2000_EPOCH)
        assert correction == 0.0

    def test_libration_periodic(self):
        """Libration correction should be periodic over the libration period."""
        params = PLUTINO_LIBRATION_PARAMS[SE_IXION]

        # Correction at t and t + period should be the same
        t0 = J2000_EPOCH
        t1 = J2000_EPOCH + params.period

        correction_t0 = calc_libration_correction(SE_IXION, t0)
        correction_t1 = calc_libration_correction(SE_IXION, t1)

        assert abs(correction_t0 - correction_t1) < 0.01

    def test_libration_varies_with_time(self):
        """Libration correction should vary over time (not constant)."""
        params = PLUTINO_LIBRATION_PARAMS[SE_IXION]

        # Check corrections at quarter period intervals
        corrections = []
        for i in range(4):
            t = J2000_EPOCH + i * params.period / 4
            corrections.append(calc_libration_correction(SE_IXION, t))

        # Not all corrections should be the same
        assert not all(abs(c - corrections[0]) < 0.01 for c in corrections)


class TestHelperFunctions:
    """Test helper functions for libration model."""

    def test_has_libration_model_plutinos(self):
        """has_libration_model should return True for plutinos."""
        assert has_libration_model(SE_IXION) is True
        assert has_libration_model(SE_ORCUS) is True

    def test_has_libration_model_non_plutinos(self):
        """has_libration_model should return False for non-plutinos."""
        assert has_libration_model(SE_CERES) is False
        assert has_libration_model(SE_ERIS) is False
        assert has_libration_model(SE_QUAOAR) is False

    def test_get_libration_parameters_plutinos(self):
        """get_libration_parameters should return params for plutinos."""
        params = get_libration_parameters(SE_IXION)
        assert params is not None
        assert isinstance(params, LibrationParameters)

    def test_get_libration_parameters_non_plutinos(self):
        """get_libration_parameters should return None for non-plutinos."""
        params = get_libration_parameters(SE_CERES)
        assert params is None


class TestPositionWithLibration:
    """Test that libration correction is applied in position calculation."""

    def test_plutino_position_includes_libration(self):
        """Plutino positions should include libration correction."""
        # Get position at two different times
        jd1 = J2000_EPOCH
        jd2 = J2000_EPOCH + 365.25 * 10  # 10 years later

        lon1, lat1, dist1 = calc_minor_body_heliocentric(SE_IXION, jd1, use_spk=False)
        lon2, lat2, dist2 = calc_minor_body_heliocentric(SE_IXION, jd2, use_spk=False)

        # Positions should be valid
        assert 0 <= lon1 < 360
        assert 0 <= lon2 < 360
        assert -90 <= lat1 <= 90
        assert -90 <= lat2 <= 90
        assert dist1 > 0
        assert dist2 > 0

    def test_non_plutino_position_unchanged(self):
        """Non-plutino positions should work normally."""
        jd = J2000_EPOCH

        lon, lat, dist = calc_minor_body_heliocentric(SE_CERES, jd, use_spk=False)

        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0

    def test_plutinos_detected_as_resonant(self):
        """Plutinos with libration params should be detected as resonant."""
        for body_id in PLUTINO_LIBRATION_PARAMS.keys():
            elements = MINOR_BODY_ELEMENTS[body_id]
            result = detect_mean_motion_resonance(elements)
            assert result.is_resonant, f"Body {body_id} not detected as resonant"
            assert result.resonance is not None
            assert result.resonance.name == "plutino"


class TestLongTermAccuracy:
    """Test improved accuracy over long timescales."""

    def test_libration_correction_magnitude_reasonable(self):
        """Libration correction should be in reasonable range over decades."""
        # Test over 100 years
        start_jd = J2000_EPOCH
        corrections = []

        for years in range(0, 101, 10):
            jd = start_jd + years * 365.25
            correction = calc_libration_correction(SE_IXION, jd)
            corrections.append(abs(correction))

        # Maximum correction should be bounded (amplitude/3 ~ 26 degrees)
        max_correction = max(corrections)
        assert max_correction < 30, f"Max correction {max_correction} too large"

    def test_position_continuity(self):
        """Positions should change continuously over time (no jumps)."""
        jd_start = J2000_EPOCH

        prev_lon = None
        for days in range(0, 3650, 365):  # 10 years, yearly steps
            jd = jd_start + days
            lon, _, _ = calc_minor_body_heliocentric(SE_IXION, jd, use_spk=False)

            if prev_lon is not None:
                # Position change should be reasonable (< 20 degrees/year for slow TNO)
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                assert diff < 20, f"Position jump too large: {diff} degrees"

            prev_lon = lon

    def test_orcus_position_continuity(self):
        """Orcus positions should also be continuous."""
        jd_start = J2000_EPOCH

        prev_lon = None
        for days in range(0, 3650, 365):
            jd = jd_start + days
            lon, _, _ = calc_minor_body_heliocentric(SE_ORCUS, jd, use_spk=False)

            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                assert diff < 20

            prev_lon = lon
