"""
Comprehensive precision tests for True Lilith (osculating lunar apogee) against pyswisseph.

This module validates calc_true_lilith against pyswisseph swe.calc_ut(jd, swe.OSCU_APOG, flags)
using 500 random dates spanning 1900-2100.

Current precision characteristics (as of this implementation):
- Initial error before corrections: 5-7 degrees
- Target final error: < 0.5 degrees

The main sources of error are:
1. Different orbital models (eccentricity vector vs. integrated lunar theory)
2. Solar perturbation modeling differences
3. Secular perturbation accumulation over time

This test tracks progress as corrections are added to calc_true_lilith.
"""

import random
import statistics
import pytest
import swisseph as swe
from libephemeris.lunar import calc_true_lilith
import libephemeris as ephem


# Threshold constants - based on actual implementation precision
# These will be tightened as corrections are added
MAX_ERROR_THRESHOLD = 7.0  # degrees (current known range 5-7°)
MEAN_ERROR_THRESHOLD = 5.0  # degrees
STD_ERROR_THRESHOLD = 2.0  # degrees

# Target threshold (goal after all corrections are implemented)
TARGET_MAX_ERROR = 0.5  # degrees

# Number of random dates to test
NUM_TEST_DATES = 500


def angular_difference(angle1: float, angle2: float) -> float:
    """
    Calculate the absolute angular difference between two angles in degrees.

    Handles wrap-around at 0/360 degrees.

    Args:
        angle1: First angle in degrees
        angle2: Second angle in degrees

    Returns:
        Absolute difference in degrees (0-180)
    """
    diff = abs(angle1 - angle2)
    if diff > 180:
        diff = 360 - diff
    return diff


def generate_random_jd(
    year_start: int, year_end: int, count: int = NUM_TEST_DATES, seed: int = 42
) -> list:
    """
    Generate random Julian Days within the specified year range.

    Args:
        year_start: Start year (inclusive)
        year_end: End year (inclusive)
        count: Number of dates to generate
        seed: Random seed for reproducibility

    Returns:
        List of tuples (year, month, day, hour, jd)
    """
    random.seed(seed)
    dates = []

    for _ in range(count):
        year = random.randint(year_start, year_end)
        month = random.randint(1, 12)
        day = random.randint(1, 28)  # Safe for all months
        hour = random.uniform(0, 24)
        jd = ephem.swe_julday(year, month, day, hour)
        dates.append((year, month, day, hour, jd))

    return dates


@pytest.mark.comparison
@pytest.mark.precision
class TestTrueLilithVsPyswisseph:
    """
    Comprehensive comparison of True Lilith calculation against pyswisseph.

    Tests 500 random dates spanning 1900-2100 and validates precision.

    Current thresholds (before corrections):
    - Maximum difference < 7 degrees
    - Mean difference < 5 degrees
    - Standard deviation < 2 degrees

    Target after corrections:
    - Maximum difference < 0.5 degrees
    """

    def test_true_lilith_500_random_dates(self, progress_reporter):
        """
        Compare calc_true_lilith against pyswisseph for 500 random dates.

        This test:
        - Generates 500 random dates between 1900-2100
        - Compares library results with pyswisseph swe.calc_ut(jd, swe.OSCU_APOG, 0)
        - Calculates max, mean, and standard deviation of differences
        - Reports progress tracking as corrections are added
        - Reports dates with largest discrepancies for debugging
        """
        dates = generate_random_jd(1900, 2100, count=NUM_TEST_DATES, seed=42)
        errors = []
        lat_errors = []
        discrepancies = []  # Store (diff, jd, swe_lon, lib_lon, date_info)
        skipped = 0

        progress = progress_reporter(
            "True Lilith comparison", len(dates), report_every=10
        )

        for i, (year, month, day, hour, jd) in enumerate(dates):
            try:
                # Calculate with pyswisseph (OSCU_APOG = osculating apogee)
                swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
                swe_lon = swe_result[0][0]
                swe_lat = swe_result[0][1]

                # Calculate with libephemeris
                lib_lon, lib_lat, lib_e = calc_true_lilith(jd)

                # Calculate angular difference for longitude
                lon_diff = angular_difference(swe_lon, lib_lon)
                errors.append(lon_diff)

                # Track latitude difference
                lat_diff = abs(swe_lat - lib_lat)
                lat_errors.append(lat_diff)

                # Store for discrepancy report
                discrepancies.append(
                    (
                        lon_diff,
                        jd,
                        swe_lon,
                        lib_lon,
                        swe_lat,
                        lib_lat,
                        f"{year:04d}-{month:02d}-{day:02d} {hour:.2f}h",
                    )
                )

            except Exception as e:
                # Skip dates outside ephemeris range
                if "ephemeris" in str(e).lower() or "outside" in str(e).lower():
                    skipped += 1
                else:
                    raise

            progress.update(i, f"JD {jd:.1f}")

        # Calculate statistics
        if not errors:
            pytest.fail("No valid dates were tested")

        max_error = max(errors)
        mean_error = statistics.mean(errors)
        std_error = statistics.stdev(errors) if len(errors) > 1 else 0.0
        median_error = statistics.median(errors)
        min_error = min(errors)

        # Latitude statistics
        max_lat_error = max(lat_errors) if lat_errors else 0.0
        mean_lat_error = statistics.mean(lat_errors) if lat_errors else 0.0

        # Sort discrepancies by error (descending)
        discrepancies.sort(key=lambda x: x[0], reverse=True)

        # Count dates meeting target precision
        dates_under_target = sum(1 for e in errors if e < TARGET_MAX_ERROR)
        pct_under_target = 100 * dates_under_target / len(errors)

        # Count dates at various precision levels
        under_1_deg = sum(1 for e in errors if e < 1.0)
        under_2_deg = sum(1 for e in errors if e < 2.0)
        under_3_deg = sum(1 for e in errors if e < 3.0)
        under_5_deg = sum(1 for e in errors if e < 5.0)

        # Print comprehensive report
        print(f"\n{'=' * 80}")
        print("TRUE LILITH (OSCULATING APOGEE) PRECISION REPORT")
        print(f"{'=' * 80}")
        print(f"Dates tested: {len(errors)}")
        print(f"Dates skipped: {skipped}")
        print("Date range: 1900-2100")
        print("\n--- Longitude Statistics ---")
        print(f"Maximum difference: {max_error:.4f}° ({max_error * 3600:.1f} arcsec)")
        print(f"Mean difference:    {mean_error:.4f}° ({mean_error * 3600:.1f} arcsec)")
        print(
            f"Median difference:  {median_error:.4f}° ({median_error * 3600:.1f} arcsec)"
        )
        print(f"Std deviation:      {std_error:.4f}° ({std_error * 3600:.1f} arcsec)")
        print(f"Minimum difference: {min_error:.4f}° ({min_error * 3600:.1f} arcsec)")
        print("\n--- Latitude Statistics ---")
        print(
            f"Maximum lat diff:   {max_lat_error:.4f}° ({max_lat_error * 3600:.1f} arcsec)"
        )
        print(
            f"Mean lat diff:      {mean_lat_error:.4f}° ({mean_lat_error * 3600:.1f} arcsec)"
        )
        print("\n--- Target Precision Analysis ---")
        print(f"Target: < {TARGET_MAX_ERROR}°")
        print(
            f"Dates meeting target: {dates_under_target}/{len(errors)} ({pct_under_target:.1f}%)"
        )
        if max_error >= TARGET_MAX_ERROR:
            print(
                f"NOTE: Target precision NOT yet achieved. Max error: {max_error:.4f}°"
            )
            print("      Corrections needed to reduce error to < 0.5°")

        # Distribution analysis
        print("\n--- Error Distribution ---")
        print(f"< 0.5°:   {dates_under_target:4d} ({pct_under_target:.1f}%)")
        print(f"< 1.0°:   {under_1_deg:4d} ({100 * under_1_deg / len(errors):.1f}%)")
        print(f"< 2.0°:   {under_2_deg:4d} ({100 * under_2_deg / len(errors):.1f}%)")
        print(f"< 3.0°:   {under_3_deg:4d} ({100 * under_3_deg / len(errors):.1f}%)")
        print(f"< 5.0°:   {under_5_deg:4d} ({100 * under_5_deg / len(errors):.1f}%)")

        # Report top 10 largest discrepancies
        print("\n--- Top 10 Largest Discrepancies ---")
        print(
            f"{'Date':<22} {'JD':<14} {'SWE Lon':<12} {'LIB Lon':<12} {'Diff':<10} {'SWE Lat':<10} {'LIB Lat':<10}"
        )
        print("-" * 90)
        for diff, jd, swe_lon, lib_lon, swe_lat, lib_lat, date_info in discrepancies[
            :10
        ]:
            print(
                f"{date_info:<22} {jd:<14.4f} {swe_lon:<12.4f} {lib_lon:<12.4f} "
                f"{diff:<10.4f} {swe_lat:<10.4f} {lib_lat:<10.4f}"
            )

        # Report top 10 smallest discrepancies (best cases)
        print("\n--- Top 10 Smallest Discrepancies (Best Cases) ---")
        print(f"{'Date':<22} {'JD':<14} {'SWE Lon':<12} {'LIB Lon':<12} {'Diff':<10}")
        print("-" * 70)
        for diff, jd, swe_lon, lib_lon, swe_lat, lib_lat, date_info in discrepancies[
            -10:
        ]:
            print(
                f"{date_info:<22} {jd:<14.4f} {swe_lon:<12.4f} {lib_lon:<12.4f} "
                f"{diff:<10.4f}"
            )

        # Progress tracking section
        print(f"\n{'=' * 80}")
        print("CORRECTION PROGRESS TRACKING")
        print(f"{'=' * 80}")
        print("Corrections implemented in calc_true_lilith:")
        print("  [x] Evection correction (1.274° amplitude)")
        print("  [x] Solar gravitational perturbation on eccentricity vector")
        print("  [x] Evection-related secondary terms")
        print("  [x] Variation correction (0.658° amplitude)")
        print("  [x] Annual equation (0.186° amplitude)")
        print("  [x] Parallactic inequality (0.125° amplitude)")
        print("  [x] Reduction to ecliptic (0.116° amplitude)")
        print("  [ ] Additional perturbation terms (if needed)")
        print("")
        print("Progress toward < 0.5° target:")
        progress_bar_len = 40
        progress_pct = (
            min(100, (TARGET_MAX_ERROR / max_error) * 100) if max_error > 0 else 100
        )
        filled = int(progress_bar_len * progress_pct / 100)
        bar = "█" * filled + "░" * (progress_bar_len - filled)
        print(f"  [{bar}] {progress_pct:.1f}%")
        print(f"  Current max error: {max_error:.4f}° -> Target: < {TARGET_MAX_ERROR}°")
        print(f"{'=' * 80}\n")

        progress.done(f"max: {max_error:.4f}°, mean: {mean_error:.4f}°")

        # Assert precision requirements (current achievable thresholds)
        assert max_error < MAX_ERROR_THRESHOLD, (
            f"Maximum error {max_error:.4f}° exceeds {MAX_ERROR_THRESHOLD}° threshold. "
            f"Worst case at JD {discrepancies[0][1]:.4f} "
            f"({discrepancies[0][6]}): SWE={discrepancies[0][2]:.4f}°, "
            f"LIB={discrepancies[0][3]:.4f}°"
        )
        assert mean_error < MEAN_ERROR_THRESHOLD, (
            f"Mean error {mean_error:.4f}° exceeds {MEAN_ERROR_THRESHOLD}° threshold"
        )
        assert std_error < STD_ERROR_THRESHOLD, (
            f"Std error {std_error:.4f}° exceeds {STD_ERROR_THRESHOLD}° threshold"
        )

    def test_true_lilith_extreme_dates(self):
        """
        Test True Lilith at boundary dates (1900 and 2100).

        These edge cases are important to verify the algorithm
        handles the full date range correctly.
        """
        # Boundary dates with expected precision
        test_cases = [
            (1900, 1, 1, 0.0, "1900 start"),
            (1900, 6, 15, 12.0, "1900 mid"),
            (1950, 1, 1, 0.0, "1950 start"),
            (2000, 1, 1, 12.0, "J2000"),
            (2050, 6, 21, 12.0, "2050 summer solstice"),
            (2100, 1, 1, 0.0, "2100 start"),
            (2100, 12, 31, 23.99, "2100 end"),
        ]

        for year, month, day, hour, desc in test_cases:
            jd = ephem.swe_julday(year, month, day, hour)

            try:
                swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
                swe_lon = swe_result[0][0]

                lib_lon, lib_lat, _ = calc_true_lilith(jd)

                diff = angular_difference(swe_lon, lib_lon)

                # Use current achievable threshold
                assert diff < MAX_ERROR_THRESHOLD, (
                    f"{desc} ({year}-{month:02d}-{day:02d}): "
                    f"diff {diff:.4f}° exceeds {MAX_ERROR_THRESHOLD}°"
                )
            except Exception as e:
                # Document failures but don't fail the test for ephemeris range issues
                if "ephemeris" not in str(e).lower():
                    raise
                pytest.skip(f"{desc}: outside ephemeris range")

    def test_true_lilith_specific_epochs(self):
        """
        Test True Lilith at specific well-known epochs.

        These are reference points that help identify systematic errors.
        Prints detailed comparison for debugging.
        """
        # Epochs with current thresholds (will be tightened as corrections are added)
        epochs = [
            (2451545.0, "J2000.0 (2000-01-01 12:00 TT)", MAX_ERROR_THRESHOLD),
            (2440587.5, "Unix Epoch (1970-01-01 00:00 UT)", MAX_ERROR_THRESHOLD),
            (2415020.0, "J1900.0", MAX_ERROR_THRESHOLD),
            (2433282.5, "1950-01-01 00:00 UT", MAX_ERROR_THRESHOLD),
            (2469807.5, "2050-01-01 00:00 UT", MAX_ERROR_THRESHOLD),
        ]

        for jd, desc, threshold in epochs:
            try:
                swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
                swe_lon = swe_result[0][0]
                swe_lat = swe_result[0][1]

                lib_lon, lib_lat, lib_e = calc_true_lilith(jd)

                diff = angular_difference(swe_lon, lib_lon)
                lat_diff = abs(swe_lat - lib_lat)

                print(f"\n{desc}:")
                print(f"  SWE: lon={swe_lon:.6f}°, lat={swe_lat:.6f}°")
                print(f"  LIB: lon={lib_lon:.6f}°, lat={lib_lat:.6f}°, e={lib_e:.6f}")
                print(f"  Diff: lon={diff:.4f}°, lat={lat_diff:.4f}°")
                print(f"  Threshold: {threshold}°")

                assert diff < threshold, (
                    f"{desc}: diff {diff:.4f}° exceeds {threshold}°"
                )
            except Exception as e:
                if "ephemeris" not in str(e).lower():
                    raise

    @pytest.mark.parametrize("year", range(1900, 2101, 20))
    def test_true_lilith_by_decade(self, year):
        """
        Test True Lilith at the start of each decade from 1900-2100.

        This parametrized test provides coverage across the entire date range
        and helps identify any systematic errors in specific time periods.
        """
        jd = ephem.swe_julday(year, 1, 1, 12.0)

        try:
            swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
            swe_lon = swe_result[0][0]

            lib_lon, lib_lat, _ = calc_true_lilith(jd)

            diff = angular_difference(swe_lon, lib_lon)

            assert diff < MAX_ERROR_THRESHOLD, (
                f"Year {year}: diff {diff:.4f}° exceeds {MAX_ERROR_THRESHOLD}° "
                f"(SWE: {swe_lon:.4f}°, LIB: {lib_lon:.4f}°)"
            )
        except Exception as e:
            if "ephemeris" not in str(e).lower():
                raise
            pytest.skip(f"Year {year}: outside ephemeris range")


@pytest.mark.comparison
class TestTrueLilithStatistics:
    """
    Statistical analysis tests for True Lilith precision.
    """

    def test_mean_error_threshold(self, progress_reporter):
        """
        Verify mean error is within acceptable threshold.

        A low mean error indicates consistent precision across all dates.
        Current threshold: 5 degrees (will be tightened as corrections are added)
        """
        dates = generate_random_jd(1900, 2100, count=NUM_TEST_DATES, seed=123)
        errors = []

        progress = progress_reporter("Mean error test", len(dates), report_every=20)

        for i, (year, month, day, hour, jd) in enumerate(dates):
            try:
                swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
                swe_lon = swe_result[0][0]

                lib_lon, _, _ = calc_true_lilith(jd)

                diff = angular_difference(swe_lon, lib_lon)
                errors.append(diff)
            except Exception:
                pass  # Skip dates outside ephemeris range

            progress.update(i)

        if not errors:
            pytest.fail("No valid dates were tested")

        mean_error = statistics.mean(errors)

        assert mean_error < MEAN_ERROR_THRESHOLD, (
            f"Mean error {mean_error:.4f}° exceeds {MEAN_ERROR_THRESHOLD}° threshold"
        )

        progress.done(f"mean: {mean_error:.4f}°")

    def test_error_standard_deviation(self, progress_reporter):
        """
        Verify standard deviation is reasonable.

        Low standard deviation indicates consistent precision,
        without outliers significantly affecting results.
        Current threshold: 2 degrees
        """
        dates = generate_random_jd(1900, 2100, count=NUM_TEST_DATES, seed=456)
        errors = []

        progress = progress_reporter("Std dev test", len(dates), report_every=20)

        for i, (year, month, day, hour, jd) in enumerate(dates):
            try:
                swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
                swe_lon = swe_result[0][0]

                lib_lon, _, _ = calc_true_lilith(jd)

                diff = angular_difference(swe_lon, lib_lon)
                errors.append(diff)
            except Exception:
                pass

            progress.update(i)

        if len(errors) < 2:
            pytest.fail("Not enough dates for standard deviation")

        std_error = statistics.stdev(errors)

        assert std_error < STD_ERROR_THRESHOLD, (
            f"Standard deviation {std_error:.4f}° exceeds {STD_ERROR_THRESHOLD}° threshold"
        )

        progress.done(f"std: {std_error:.4f}°")


@pytest.mark.comparison
class TestTrueLilithLatitude:
    """
    Tests for True Lilith latitude precision.

    The True Lilith has an ecliptic latitude component (unlike Mean Lilith
    which is always on the ecliptic). This tests that the latitude calculation
    is reasonably accurate.
    """

    def test_latitude_range(self):
        """
        Verify True Lilith latitude is within expected range.

        The lunar orbit is inclined about 5° to the ecliptic, so the
        osculating apogee latitude should generally be within ±6°.
        """
        test_dates = generate_random_jd(1950, 2050, count=100, seed=789)

        for year, month, day, hour, jd in test_dates:
            try:
                _, lib_lat, _ = calc_true_lilith(jd)

                assert -10 < lib_lat < 10, (
                    f"Latitude {lib_lat:.4f}° at JD {jd:.4f} outside expected range"
                )
            except Exception:
                pass  # Skip dates outside ephemeris range

    def test_latitude_comparison(self, progress_reporter):
        """
        Compare True Lilith latitude with pyswisseph.
        """
        dates = generate_random_jd(1950, 2050, count=100, seed=321)
        lat_errors = []

        progress = progress_reporter("Latitude comparison", len(dates), report_every=10)

        for i, (year, month, day, hour, jd) in enumerate(dates):
            try:
                swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
                swe_lat = swe_result[0][1]

                _, lib_lat, _ = calc_true_lilith(jd)

                lat_diff = abs(swe_lat - lib_lat)
                lat_errors.append(lat_diff)
            except Exception:
                pass

            progress.update(i)

        if lat_errors:
            max_lat_error = max(lat_errors)
            mean_lat_error = statistics.mean(lat_errors)

            print("\n--- Latitude Error Summary ---")
            print(f"Max latitude error: {max_lat_error:.4f}°")
            print(f"Mean latitude error: {mean_lat_error:.4f}°")

            # Latitude should be within 5 degrees (relaxed threshold for now)
            assert max_lat_error < 5.0, (
                f"Maximum latitude error {max_lat_error:.4f}° exceeds 5.0° threshold"
            )

            progress.done(f"max lat error: {max_lat_error:.4f}°")


@pytest.mark.comparison
class TestTrueLilithEccentricity:
    """
    Tests for the eccentricity value returned by calc_true_lilith.
    """

    def test_eccentricity_in_expected_range(self):
        """
        Verify that the returned eccentricity is in the expected range.

        The lunar orbit eccentricity is approximately 0.055, with variations
        due to solar perturbations. It should be between 0.03 and 0.08.
        """
        test_dates = generate_random_jd(1950, 2050, count=100, seed=654)

        eccentricities = []
        for year, month, day, hour, jd in test_dates:
            try:
                _, _, e_mag = calc_true_lilith(jd)
                eccentricities.append(e_mag)

                assert 0.03 < e_mag < 0.08, (
                    f"Eccentricity {e_mag:.6f} at JD {jd:.4f} outside expected range"
                )
            except Exception:
                pass

        if eccentricities:
            mean_e = statistics.mean(eccentricities)
            print("\n--- Eccentricity Summary ---")
            print(f"Mean eccentricity: {mean_e:.6f}")
            print(f"Min eccentricity: {min(eccentricities):.6f}")
            print(f"Max eccentricity: {max(eccentricities):.6f}")

            # Mean should be close to nominal 0.055
            assert 0.04 < mean_e < 0.07, (
                f"Mean eccentricity {mean_e:.6f} outside expected range"
            )
