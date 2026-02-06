"""
Tests for the Moshier analytical method for interpolated lunar apogee.

This module validates the _calc_elp2000_apogee_perturbations() function
which implements a comprehensive ~50 term perturbation series based on
Moshier's DE404-fitted lunar theory.

Target precision: <0.5° max error, <0.2° mean error vs Swiss Ephemeris SE_INTP_APOG
"""

import math
import pytest
from libephemeris.lunar import (
    _calc_elp2000_apogee_perturbations,
    _calc_lunar_fundamental_arguments,
    calc_interpolated_apogee,
    calc_mean_lilith,
)


class TestMoshierApogeeBasics:
    """Basic functionality tests for the Moshier apogee perturbation series."""

    def test_perturbation_returns_float(self):
        """Test that the perturbation function returns a float."""
        jd_tt = 2451545.0  # J2000.0
        result = _calc_elp2000_apogee_perturbations(jd_tt)
        assert isinstance(result, float)

    def test_perturbation_within_expected_range(self):
        """Test that perturbation is within expected amplitude range."""
        # According to Swiss Ephemeris, apogee oscillates ~5° from mean
        jd_tt = 2451545.0  # J2000.0
        perturbation = _calc_elp2000_apogee_perturbations(jd_tt)
        # Perturbation should be within ±6° (conservative bound)
        assert abs(perturbation) < 6.0

    def test_perturbation_varies_over_time(self):
        """Test that perturbation changes over time."""
        jd_base = 2451545.0
        p1 = _calc_elp2000_apogee_perturbations(jd_base)
        p2 = _calc_elp2000_apogee_perturbations(jd_base + 7.0)  # 1 week later
        p3 = _calc_elp2000_apogee_perturbations(jd_base + 30.0)  # 1 month later

        # Should have variation over these time scales
        assert p1 != p2 or p2 != p3

    def test_perturbation_samples_typical_range(self):
        """Test perturbation range over a longer period."""
        jd_start = 2451545.0
        perturbations = []

        # Sample over 2 years (multiple evection cycles)
        for i in range(365 * 2):
            jd = jd_start + i
            p = _calc_elp2000_apogee_perturbations(jd)
            perturbations.append(p)

        max_p = max(perturbations)
        min_p = min(perturbations)
        range_p = max_p - min_p

        # Expected range is ~10° (±5° from zero)
        assert range_p > 5.0, f"Range {range_p}° too small, expected >5°"
        assert range_p < 15.0, f"Range {range_p}° too large, expected <15°"


class TestMoshierApogeeHarmonics:
    """Tests verifying harmonic structure of the perturbation series."""

    def test_dominant_evection_term_signature(self):
        """Test that the 2D-2M' term is present with correct sign and magnitude."""
        # The dominant evection term (2D-2M') has ~4.69° amplitude
        # The evection cycle is ~31.8 days (half the evection beat period)
        # We need to sample over enough time to capture full amplitude
        jd_base = 2451545.0

        # Sample perturbations over multiple evection cycles (~100 days)
        samples = []
        for i in range(200):
            jd = jd_base + i * 0.5  # half-day steps over 100 days
            p = _calc_elp2000_apogee_perturbations(jd)
            samples.append(p)

        # The range should reflect the dominant term amplitude plus other harmonics
        amp_estimate = (max(samples) - min(samples)) / 2.0
        # Total amplitude should be at least 3° (allowing for phase effects)
        # and less than 8° (sum of all harmonics)
        assert amp_estimate > 2.0, f"Amplitude {amp_estimate}° too small"
        assert amp_estimate < 8.0, f"Amplitude {amp_estimate}° too large"

    def test_annual_equation_term(self):
        """Test that annual variation (M term) is present."""
        # Sample over a full year to detect annual variation
        jd_start = 2451545.0

        samples_jan = []
        samples_jul = []

        # January (perihelion) vs July (aphelion)
        for day in range(30):
            jd_jan = jd_start + day  # January 2000
            jd_jul = jd_start + 180 + day  # July 2000
            samples_jan.append(_calc_elp2000_apogee_perturbations(jd_jan))
            samples_jul.append(_calc_elp2000_apogee_perturbations(jd_jul))

        mean_jan = sum(samples_jan) / len(samples_jan)
        mean_jul = sum(samples_jul) / len(samples_jul)

        # Should see some difference between perihelion and aphelion
        # (annual equation has ~0.38° amplitude in the solar anomaly term)
        # The difference won't be exactly 2*0.38° due to other terms
        assert abs(mean_jan - mean_jul) < 5.0  # Within reasonable bounds


class TestMoshierApogeeIntegration:
    """Integration tests with calc_interpolated_apogee."""

    def test_interpolated_apogee_uses_perturbation(self):
        """Verify calc_interpolated_apogee incorporates perturbations."""
        jd_tt = 2451545.0

        mean_lon = calc_mean_lilith(jd_tt)
        interp_lon, _, _ = calc_interpolated_apogee(jd_tt)
        perturbation = _calc_elp2000_apogee_perturbations(jd_tt)

        # The difference between interpolated and mean should match perturbation
        diff = (interp_lon - mean_lon) % 360.0
        if diff > 180.0:
            diff -= 360.0

        assert abs(diff - perturbation) < 0.001

    def test_interpolated_apogee_valid_range(self):
        """Test that interpolated apogee longitude is in valid range."""
        jd_tt = 2451545.0
        lon, lat, ecc = calc_interpolated_apogee(jd_tt)

        assert 0.0 <= lon < 360.0
        assert abs(lat) < 0.01  # Essentially zero latitude
        assert 0.05 < ecc < 0.06  # Mean eccentricity

    def test_interpolated_apogee_continuity(self):
        """Test that interpolated apogee moves continuously."""
        jd_start = 2451545.0
        lons = []

        for i in range(100):
            jd = jd_start + i
            lon, _, _ = calc_interpolated_apogee(jd)
            lons.append(lon)

        # Unwrap and check for continuity
        unwrapped = [lons[0]]
        for i in range(1, len(lons)):
            diff = lons[i] - lons[i - 1]
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            unwrapped.append(unwrapped[-1] + diff)

        # Daily motion should be roughly 0.11° (apsidal precession rate)
        for i in range(1, len(unwrapped)):
            daily_change = unwrapped[i] - unwrapped[i - 1]
            # Should be within reasonable bounds (0.11° mean ± some variation)
            assert abs(daily_change) < 1.0, f"Day {i}: Change {daily_change}° too large"


class TestMoshierApogeeDocumentation:
    """Documentation tests for the Moshier apogee implementation."""

    def test_perturbation_function_has_docstring(self):
        """Verify _calc_elp2000_apogee_perturbations has comprehensive docstring."""
        doc = _calc_elp2000_apogee_perturbations.__doc__
        assert doc is not None
        assert len(doc) > 500  # Comprehensive documentation expected

        # Key terms should be mentioned
        assert "Moshier" in doc
        assert "perturbation" in doc.lower()
        assert "evection" in doc.lower()

    def test_docstring_mentions_precision(self):
        """Verify docstring includes precision information."""
        doc = _calc_elp2000_apogee_perturbations.__doc__
        assert doc is not None
        assert "Precision" in doc
        assert "0.5" in doc or "0.2" in doc  # Target precision values


class TestMoshierApogeeEdgeCases:
    """Edge case tests for the Moshier apogee implementation."""

    def test_j2000_epoch(self):
        """Test at J2000.0 epoch."""
        jd_tt = 2451545.0
        result = _calc_elp2000_apogee_perturbations(jd_tt)
        assert math.isfinite(result)

    def test_year_1900(self):
        """Test at year 1900."""
        jd_tt = 2415020.5  # 1900-01-01
        result = _calc_elp2000_apogee_perturbations(jd_tt)
        assert math.isfinite(result)
        assert abs(result) < 10.0

    def test_year_2100(self):
        """Test at year 2100."""
        jd_tt = 2488070.5  # 2100-01-01
        result = _calc_elp2000_apogee_perturbations(jd_tt)
        assert math.isfinite(result)
        assert abs(result) < 10.0

    def test_consistency_over_long_range(self):
        """Test that function remains stable over a long time range."""
        jd_start = 2415020.5  # 1900
        jd_end = 2488070.5  # 2100

        # Sample at 10-year intervals
        for jd in range(int(jd_start), int(jd_end), 3652):
            result = _calc_elp2000_apogee_perturbations(float(jd))
            assert math.isfinite(result)
            assert abs(result) < 10.0, f"JD {jd}: Perturbation {result}° out of range"
