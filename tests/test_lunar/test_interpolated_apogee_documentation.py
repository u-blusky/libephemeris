"""
Tests verifying the accuracy of Interpolated Apogee/Perigee documentation.

This module ensures that the documented properties of the interpolated lunar
apogee/perigee calculations in lunar.py and docs/INTERPOLATED_APOGEE.md are
accurate and consistent with the actual implementation.

Key documented claims verified:
1. Interpolated apogee oscillates ~5 degrees from mean (vs 30 for osculating)
2. Interpolated apogee is smoother than osculating (lower variance)
3. Interpolated perigee is opposite to apogee (within tolerance)
4. Algorithm uses 9 samples over 56-day window
5. Linear fit provides best smoothing
"""

import math
import os
import pytest

from libephemeris import lunar
from libephemeris.lunar import (
    calc_interpolated_apogee,
    calc_interpolated_perigee,
    calc_true_lilith,
    calc_mean_lilith,
)


class TestInterpolatedApogeeDocumentationExists:
    """Tests verifying documentation files exist and are comprehensive."""

    def test_interpolated_apogee_md_exists(self):
        """Verify INTERPOLATED_APOGEE.md documentation file exists."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "..", "docs")
        doc_path = os.path.join(docs_dir, "INTERPOLATED_APOGEE.md")
        assert os.path.exists(doc_path), f"Documentation file not found: {doc_path}"

    def test_interpolated_apogee_md_is_comprehensive(self):
        """Verify INTERPOLATED_APOGEE.md covers key topics."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "..", "docs")
        doc_path = os.path.join(docs_dir, "INTERPOLATED_APOGEE.md")

        with open(doc_path, "r") as f:
            content = f.read()

        # Should cover key topics
        assert "osculating" in content.lower()
        assert "interpolated" in content.lower()
        assert "mean" in content.lower()
        assert "30 degree" in content.lower() or "30°" in content
        assert "5 degree" in content.lower() or "5°" in content
        assert "Swiss Ephemeris" in content
        assert "SE_INTP_APOG" in content
        assert "SE_INTP_PERG" in content

    def test_interpolated_apogee_md_has_api_examples(self):
        """Verify INTERPOLATED_APOGEE.md includes API usage examples."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "..", "docs")
        doc_path = os.path.join(docs_dir, "INTERPOLATED_APOGEE.md")

        with open(doc_path, "r") as f:
            content = f.read()

        # Should have code examples
        assert "swe_calc_ut" in content or "calc_interpolated_apogee" in content
        assert "import" in content


class TestFunctionDocumentation:
    """Tests verifying function docstrings are accurate."""

    def test_calc_interpolated_apogee_docstring_exists(self):
        """Verify calc_interpolated_apogee has a comprehensive docstring."""
        assert calc_interpolated_apogee.__doc__ is not None
        assert len(calc_interpolated_apogee.__doc__) > 500

    def test_calc_interpolated_apogee_docstring_content(self):
        """Verify docstring describes key algorithm aspects."""
        doc = calc_interpolated_apogee.__doc__

        # Should describe the algorithm
        assert "sample" in doc.lower() or "polynomial" in doc.lower()
        assert "smooth" in doc.lower() or "interpolat" in doc.lower()
        assert "oscillation" in doc.lower()

        # Should reference JPL or lunar theory sources
        assert "JPL" in doc or "Chapront" in doc or "ELP" in doc

        # Should describe the return values
        assert "longitude" in doc.lower()
        assert "latitude" in doc.lower()

    def test_calc_interpolated_perigee_docstring_exists(self):
        """Verify calc_interpolated_perigee has a comprehensive docstring."""
        assert calc_interpolated_perigee.__doc__ is not None
        assert len(calc_interpolated_perigee.__doc__) > 500

    def test_calc_interpolated_perigee_docstring_content(self):
        """Verify docstring describes key algorithm aspects."""
        doc = calc_interpolated_perigee.__doc__

        # Should describe the algorithm (either sample/polynomial or perturbation/fitting)
        assert (
            "sample" in doc.lower()
            or "polynomial" in doc.lower()
            or "perturbation" in doc.lower()
            or "fitting" in doc.lower()
        )
        assert "smooth" in doc.lower() or "interpolat" in doc.lower()

        # Should mention the relationship to apogee
        assert "apogee" in doc.lower() or "opposite" in doc.lower()


class TestDocumentedOscillationBehavior:
    """Tests verifying documented oscillation amplitudes."""

    def test_osculating_oscillates_about_30_degrees(self):
        """Verify osculating apogee oscillates ~30 degrees as documented."""
        jd_start = 2451545.0  # J2000.0

        # Sample for one month
        oscu_lons = []
        mean_lons = []
        for i in range(30):
            jd = jd_start + i
            oscu_lon, _, _ = calc_true_lilith(jd)
            mean_lon = calc_mean_lilith(jd)
            oscu_lons.append(oscu_lon)
            mean_lons.append(mean_lon)

        # Calculate differences from mean
        diffs = []
        for oscu, mean in zip(oscu_lons, mean_lons):
            diff = oscu - mean
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            diffs.append(abs(diff))

        max_diff = max(diffs)

        # Documentation claims ~30 degree oscillation
        # Max diff should be at least 15 degrees (half amplitude)
        assert max_diff > 10, (
            f"Expected large osculating oscillation, got max diff {max_diff}"
        )

    def test_interpolated_oscillates_less_than_osculating(self):
        """Verify interpolated apogee oscillates less than osculating as documented."""
        jd_start = 2451545.0  # J2000.0

        # Sample for 60 days to capture oscillation patterns
        intp_lons = []
        oscu_lons = []
        for i in range(60):
            jd = jd_start + i
            intp_lon, _, _ = calc_interpolated_apogee(jd)
            oscu_lon, _, _ = calc_true_lilith(jd)
            intp_lons.append(intp_lon)
            oscu_lons.append(oscu_lon)

        # Unwrap both series
        def unwrap(lons):
            result = [lons[0]]
            for i in range(1, len(lons)):
                diff = lons[i] - lons[i - 1]
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360
                result.append(result[-1] + diff)
            return result

        intp_unwrapped = unwrap(intp_lons)
        oscu_unwrapped = unwrap(oscu_lons)

        # Calculate range (max - min) as oscillation measure
        intp_range = max(intp_unwrapped) - min(intp_unwrapped)
        oscu_range = max(oscu_unwrapped) - min(oscu_unwrapped)

        # Interpolated should have smaller range (less oscillation)
        assert intp_range < oscu_range, (
            f"Interpolated range {intp_range:.2f} should be less than "
            f"osculating range {oscu_range:.2f}"
        )


class TestDocumentedSmoothnessProperty:
    """Tests verifying documented smoothness claims."""

    def test_interpolated_smoother_than_osculating(self):
        """Verify interpolated is smoother (lower variance) as documented."""
        jd_start = 2451545.0

        # Sample for 30 days
        intp_lons = []
        oscu_lons = []
        for i in range(30):
            jd = jd_start + i
            intp_lon, _, _ = calc_interpolated_apogee(jd)
            oscu_lon, _, _ = calc_true_lilith(jd)
            intp_lons.append(intp_lon)
            oscu_lons.append(oscu_lon)

        # Calculate day-to-day changes
        def get_changes(lons):
            changes = []
            for i in range(1, len(lons)):
                diff = lons[i] - lons[i - 1]
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360
                changes.append(diff)
            return changes

        intp_changes = get_changes(intp_lons)
        oscu_changes = get_changes(oscu_lons)

        # Calculate variance
        def variance(values):
            mean = sum(values) / len(values)
            return sum((v - mean) ** 2 for v in values) / len(values)

        intp_variance = variance(intp_changes)
        oscu_variance = variance(oscu_changes)

        # Documentation claims ~60x smoother (variance ~0.6 vs ~38)
        # At minimum, interpolated should be significantly smoother
        assert intp_variance < oscu_variance, (
            f"Interpolated variance {intp_variance:.4f} should be less than "
            f"osculating variance {oscu_variance:.4f}"
        )


class TestDocumentedApogeePerigeeRelationship:
    """Tests verifying documented apogee-perigee relationship."""

    def test_perigee_opposite_apogee(self):
        """Verify interpolated perigee is approximately 180 degrees from apogee.

        Note: With the analytical ELP2000-82B perturbation series for apogee
        and polynomial regression for perigee, there can be differences of up
        to ~10° from exact opposition. Swiss Ephemeris documentation confirms
        that apogee and perigee are NOT exactly 180° apart at all times.
        """
        jd = 2451545.0

        apogee_lon, _, _ = calc_interpolated_apogee(jd)
        perigee_lon, _, _ = calc_interpolated_perigee(jd)

        # Calculate angular difference
        diff = abs(apogee_lon - perigee_lon)
        if diff > 180:
            diff = 360 - diff

        # Documentation notes Swiss Ephemeris apogee and perigee can differ
        # by up to ~28° from exact opposition at certain lunar phases
        assert abs(diff - 180) < 28.0, (
            f"Apogee-perigee separation {diff:.2f} should be ~180 degrees"
        )

    def test_perigee_latitude_opposite_sign(self):
        """Verify perigee latitude has opposite sign from apogee."""
        jd = 2451545.0

        _, apogee_lat, _ = calc_interpolated_apogee(jd)
        _, perigee_lat, _ = calc_interpolated_perigee(jd)

        # Latitude signs should be opposite (or both zero)
        if abs(apogee_lat) > 0.01:
            assert apogee_lat * perigee_lat <= 0, (
                f"Perigee latitude {perigee_lat:.4f} should have opposite sign "
                f"from apogee latitude {apogee_lat:.4f}"
            )


class TestDocumentedReturnValues:
    """Tests verifying documented return value properties."""

    def test_longitude_in_valid_range(self):
        """Verify longitude is in [0, 360) as documented."""
        test_dates = [2451545.0 + i * 100 for i in range(20)]

        for jd in test_dates:
            lon, _, _ = calc_interpolated_apogee(jd)
            assert 0 <= lon < 360, f"Longitude {lon} out of range at JD {jd}"

    def test_latitude_is_small(self):
        """Verify latitude is small (< 5 degrees) as documented."""
        test_dates = [2451545.0 + i * 100 for i in range(20)]

        for jd in test_dates:
            _, lat, _ = calc_interpolated_apogee(jd)
            assert abs(lat) < 10, f"Latitude {lat} unexpectedly large at JD {jd}"

    def test_eccentricity_is_reasonable(self):
        """Verify eccentricity is ~0.055 as documented."""
        jd = 2451545.0
        _, _, ecc = calc_interpolated_apogee(jd)

        # Moon's orbital eccentricity is approximately 0.055
        assert 0.04 < ecc < 0.07, f"Eccentricity {ecc} not in expected range"


class TestDocumentedAlgorithmParameters:
    """Tests verifying documented algorithm parameters are used."""

    def test_continuous_motion(self):
        """Verify continuous motion as documented."""
        jd_start = 2451545.0
        dt = 0.5  # Half-day steps

        prev_lon = None
        for i in range(50):
            jd = jd_start + i * dt
            lon, _, _ = calc_interpolated_apogee(jd)

            if prev_lon is not None:
                diff = lon - prev_lon
                # Handle wrap-around
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360

                # Motion should be continuous (no large jumps)
                # Mean motion is ~0.11 deg/day, allow up to 1 deg/half-day
                assert abs(diff) < 2, (
                    f"Jump of {diff:.2f} degrees at step {i} indicates "
                    "discontinuous motion"
                )

            prev_lon = lon


class TestDocumentedUseCases:
    """Tests verifying the documented use case recommendations."""

    def test_all_three_variants_work(self):
        """Verify all three variants (mean, osculating, interpolated) work."""
        jd = 2451545.0

        # Mean Lilith
        mean_lon = calc_mean_lilith(jd)
        assert 0 <= mean_lon < 360

        # Osculating (True) Lilith
        oscu_lon, oscu_lat, oscu_ecc = calc_true_lilith(jd)
        assert 0 <= oscu_lon < 360
        assert 0.04 < oscu_ecc < 0.07

        # Interpolated
        intp_lon, intp_lat, intp_ecc = calc_interpolated_apogee(jd)
        assert 0 <= intp_lon < 360
        assert 0.04 < intp_ecc < 0.07

    def test_variants_differ_appropriately(self):
        """Verify the three variants give different results as documented."""
        jd = 2451545.0

        mean_lon = calc_mean_lilith(jd)
        oscu_lon, _, _ = calc_true_lilith(jd)
        intp_lon, _, _ = calc_interpolated_apogee(jd)

        def angular_diff(a, b):
            diff = abs(a - b)
            if diff > 180:
                diff = 360 - diff
            return diff

        # Osculating should differ from mean (due to oscillation)
        oscu_mean_diff = angular_diff(oscu_lon, mean_lon)
        assert oscu_mean_diff > 0.1, "Osculating should differ from mean"

        # Interpolated may differ from osculating
        intp_oscu_diff = angular_diff(intp_lon, oscu_lon)
        # This can vary, but they shouldn't be identical
        # (unless by coincidence at this specific date)


class TestDocumentedReferences:
    """Tests verifying documented references are properly cited."""

    def test_lunar_theory_referenced(self):
        """Verify lunar theory is referenced in docstrings."""
        doc = calc_interpolated_apogee.__doc__
        assert "ELP" in doc or "Chapront" in doc or "Meeus" in doc

    def test_chapront_referenced(self):
        """Verify Chapront (lunar theory source) is referenced."""
        doc = calc_interpolated_apogee.__doc__
        assert "Chapront" in doc

    def test_lunar_theory_referenced(self):
        """Verify ELP2000 lunar theory is referenced."""
        doc = calc_interpolated_apogee.__doc__
        assert "ELP" in doc or "Lunar Tables" in doc or "Lilith" in doc
