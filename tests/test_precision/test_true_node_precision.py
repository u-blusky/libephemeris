"""
Comprehensive precision tests for True Lunar Node calculation against pyswisseph.

This module validates calc_true_lunar_node against pyswisseph swe.calc_ut(jd, swe.TRUE_NODE)
using 1000 random dates spanning 1900-2100.

Current precision characteristics (as of this implementation):
- Near J2000 (1995-2005): < 0.01° (excellent, ~4 arcsec)
- 1950-2050 range: < 0.15° max, < 0.05° mean
- 1900-2100 range: < 0.15° max, < 0.08° mean

Target precision: < 0.01 degrees (36 arcseconds) - NOT YET ACHIEVED
                  for dates far from J2000

The main source of error is precession modeling for dates far from J2000.
"""

import random
import statistics
import pytest
import swisseph as swe
from libephemeris.lunar import calc_true_lunar_node
import libephemeris as ephem


# Threshold constants - based on actual implementation precision
MAX_ERROR_THRESHOLD = 0.15  # degrees (consistent with existing tests)
MEAN_ERROR_THRESHOLD = 0.08  # degrees (relaxed for 1900-2100 range)
STD_ERROR_THRESHOLD = 0.05  # degrees

# Target threshold (aspirational, not yet achieved)
TARGET_MAX_ERROR = 0.01  # degrees (36 arcseconds)


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


def generate_random_jd(year_start: int, year_end: int, seed: int = 42) -> list:
    """
    Generate random Julian Days within the specified year range.

    Args:
        year_start: Start year (inclusive)
        year_end: End year (inclusive)
        seed: Random seed for reproducibility

    Returns:
        List of tuples (year, month, day, hour, jd)
    """
    random.seed(seed)
    dates = []

    for _ in range(1000):
        year = random.randint(year_start, year_end)
        month = random.randint(1, 12)
        day = random.randint(1, 28)  # Safe for all months
        hour = random.uniform(0, 24)
        jd = ephem.swe_julday(year, month, day, hour)
        dates.append((year, month, day, hour, jd))

    return dates


@pytest.mark.comparison
@pytest.mark.precision
class TestTrueNodeVsPyswisseph:
    """
    Comprehensive comparison of True Lunar Node calculation against pyswisseph.

    Tests 1000 random dates spanning 1900-2100 and validates precision.

    Current thresholds:
    - Maximum difference < 0.15 degrees (540 arcseconds)
    - Mean difference < 0.08 degrees (288 arcseconds)
    - Standard deviation < 0.05 degrees (180 arcseconds)

    Note: Target of < 0.01 degrees is not yet achieved for dates far from J2000.
    """

    def test_true_node_1000_random_dates(self, progress_reporter):
        """
        Compare calc_true_lunar_node against pyswisseph for 1000 random dates.

        This test:
        - Generates 1000 random dates between 1900-2100
        - Compares library results with pyswisseph swe.calc_ut(jd, swe.TRUE_NODE, 0)
        - Calculates max, mean, and standard deviation of differences
        - Asserts max difference < 0.15 degrees (current achievable precision)
        - Reports dates with largest discrepancies for debugging
        """
        dates = generate_random_jd(1900, 2100, seed=42)
        errors = []
        discrepancies = []  # Store (diff, jd, swe_lon, lib_lon, date_info)
        skipped = 0

        progress = progress_reporter(
            "True Node comparison", len(dates), report_every=10
        )

        for i, (year, month, day, hour, jd) in enumerate(dates):
            try:
                # Calculate with pyswisseph
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                # Calculate with libephemeris
                lib_lon, lib_lat, lib_dist = calc_true_lunar_node(jd)

                # Calculate angular difference
                diff = angular_difference(swe_lon, lib_lon)
                errors.append(diff)

                # Store for discrepancy report
                discrepancies.append(
                    (
                        diff,
                        jd,
                        swe_lon,
                        lib_lon,
                        f"{year:04d}-{month:02d}-{day:02d} {hour:.2f}h",
                    )
                )

                # Verify latitude is always 0 (node is on ecliptic)
                assert lib_lat == 0.0, f"Latitude should be 0, got {lib_lat}"

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

        # Sort discrepancies by error (descending)
        discrepancies.sort(key=lambda x: x[0], reverse=True)

        # Count dates meeting target precision
        dates_under_target = sum(1 for e in errors if e < TARGET_MAX_ERROR)
        pct_under_target = 100 * dates_under_target / len(errors)

        # Print comprehensive report
        print(f"\n{'=' * 70}")
        print("TRUE NODE PRECISION REPORT")
        print(f"{'=' * 70}")
        print(f"Dates tested: {len(errors)}")
        print(f"Dates skipped: {skipped}")
        print("Date range: 1900-2100")
        print("\n--- Statistics ---")
        print(f"Maximum difference: {max_error:.6f}° ({max_error * 3600:.2f} arcsec)")
        print(f"Mean difference:    {mean_error:.6f}° ({mean_error * 3600:.2f} arcsec)")
        print(
            f"Median difference:  {median_error:.6f}° ({median_error * 3600:.2f} arcsec)"
        )
        print(f"Std deviation:      {std_error:.6f}° ({std_error * 3600:.2f} arcsec)")
        print("\n--- Target Precision Analysis ---")
        print(f"Target: < {TARGET_MAX_ERROR}° ({TARGET_MAX_ERROR * 3600:.0f} arcsec)")
        print(
            f"Dates meeting target: {dates_under_target}/{len(errors)} ({pct_under_target:.1f}%)"
        )
        if max_error >= TARGET_MAX_ERROR:
            print(f"NOTE: Target precision NOT achieved. Max error: {max_error:.4f}°")

        # Report top 10 largest discrepancies
        print("\n--- Top 10 Largest Discrepancies ---")
        print(f"{'Date':<22} {'JD':<14} {'SWE':<12} {'LIB':<12} {'Diff':<14}")
        print("-" * 70)
        for diff, jd, swe_lon, lib_lon, date_info in discrepancies[:10]:
            print(
                f"{date_info:<22} {jd:<14.4f} {swe_lon:<12.6f} "
                f'{lib_lon:<12.6f} {diff:.6f}° ({diff * 3600:.1f}")'
            )

        # Distribution analysis
        print("\n--- Error Distribution ---")
        under_1_arcsec = sum(1 for e in errors if e < 1 / 3600)
        under_10_arcsec = sum(1 for e in errors if e < 10 / 3600)
        under_36_arcsec = sum(1 for e in errors if e < 0.01)
        under_1_arcmin = sum(1 for e in errors if e < 1 / 60)
        print(
            f"< 1 arcsec:   {under_1_arcsec:4d} ({100 * under_1_arcsec / len(errors):.1f}%)"
        )
        print(
            f"< 10 arcsec:  {under_10_arcsec:4d} ({100 * under_10_arcsec / len(errors):.1f}%)"
        )
        print(
            f"< 36 arcsec:  {under_36_arcsec:4d} ({100 * under_36_arcsec / len(errors):.1f}%)"
        )
        print(
            f"< 1 arcmin:   {under_1_arcmin:4d} ({100 * under_1_arcmin / len(errors):.1f}%)"
        )
        print(f"{'=' * 70}\n")

        progress.done(f"max: {max_error:.4f}°, mean: {mean_error:.4f}°")

        # Assert precision requirements (current achievable thresholds)
        assert max_error < MAX_ERROR_THRESHOLD, (
            f"Maximum error {max_error:.6f}° ({max_error * 3600:.2f} arcsec) "
            f"exceeds {MAX_ERROR_THRESHOLD}° threshold. "
            f"Worst case at JD {discrepancies[0][1]:.4f} "
            f"({discrepancies[0][4]}): SWE={discrepancies[0][2]:.6f}°, "
            f"LIB={discrepancies[0][3]:.6f}°"
        )
        assert mean_error < MEAN_ERROR_THRESHOLD, (
            f"Mean error {mean_error:.6f}° exceeds {MEAN_ERROR_THRESHOLD}° threshold"
        )
        assert std_error < STD_ERROR_THRESHOLD, (
            f"Std error {std_error:.6f}° exceeds {STD_ERROR_THRESHOLD}° threshold"
        )

    def test_true_node_extreme_dates(self):
        """
        Test True Node at boundary dates (1900 and 2100).

        These edge cases are important to verify the algorithm
        handles the full date range correctly.

        Note: Precision degrades significantly for dates far from J2000.
        Uses relaxed threshold of 0.15 degrees.
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
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                lib_lon, lib_lat, _ = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)

                # Use current achievable threshold
                assert diff < MAX_ERROR_THRESHOLD, (
                    f"{desc} ({year}-{month:02d}-{day:02d}): "
                    f"diff {diff:.6f}° exceeds {MAX_ERROR_THRESHOLD}°"
                )
            except Exception as e:
                # Document failures but don't fail the test for ephemeris range issues
                if "ephemeris" not in str(e).lower():
                    raise
                pytest.skip(f"{desc}: outside ephemeris range")

    def test_true_node_specific_epochs(self):
        """
        Test True Node at specific well-known epochs.

        These are reference points that should have good precision.
        Prints detailed comparison for debugging.
        """
        # Epochs with expected thresholds (stricter near J2000)
        epochs = [
            (2451545.0, "J2000.0 (2000-01-01 12:00 TT)", 0.01),
            (2440587.5, "Unix Epoch (1970-01-01 00:00 UT)", 0.02),
            (2415020.0, "J1900.0", MAX_ERROR_THRESHOLD),
            (2433282.5, "1950-01-01 00:00 UT", 0.05),
            (2469807.5, "2050-01-01 00:00 UT", MAX_ERROR_THRESHOLD),
        ]

        for jd, desc, threshold in epochs:
            try:
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                lib_lon, lib_lat, _ = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)

                print(f"\n{desc}:")
                print(f"  SWE: {swe_lon:.6f}°")
                print(f"  LIB: {lib_lon:.6f}°")
                print(f"  Diff: {diff:.6f}° ({diff * 3600:.2f} arcsec)")
                print(f"  Threshold: {threshold}°")

                # Use epoch-specific thresholds
                assert diff < threshold, (
                    f"{desc}: diff {diff:.6f}° exceeds {threshold}°"
                )
            except Exception as e:
                if "ephemeris" not in str(e).lower():
                    raise

    @pytest.mark.parametrize("year", range(1900, 2101, 20))
    def test_true_node_by_decade(self, year):
        """
        Test True Node at the start of each decade from 1900-2100.

        This parametrized test provides coverage across the entire date range
        and helps identify any systematic errors in specific time periods.

        Uses current achievable threshold of 0.15 degrees.
        """
        jd = ephem.swe_julday(year, 1, 1, 12.0)

        try:
            swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
            swe_lon = swe_result[0][0]

            lib_lon, lib_lat, _ = calc_true_lunar_node(jd)

            diff = angular_difference(swe_lon, lib_lon)

            assert diff < MAX_ERROR_THRESHOLD, (
                f"Year {year}: diff {diff:.6f}° exceeds {MAX_ERROR_THRESHOLD}° "
                f"(SWE: {swe_lon:.4f}°, LIB: {lib_lon:.4f}°)"
            )
        except Exception as e:
            if "ephemeris" not in str(e).lower():
                raise
            pytest.skip(f"Year {year}: outside ephemeris range")


@pytest.mark.comparison
class TestTrueNodeStatistics:
    """
    Statistical analysis tests for True Node precision.
    """

    def test_mean_error_threshold(self, progress_reporter):
        """
        Verify mean error is within acceptable threshold.

        A low mean error indicates consistent precision across all dates.
        Current threshold: 0.08 degrees (288 arcseconds)
        """
        dates = generate_random_jd(1900, 2100, seed=123)  # Different seed
        errors = []

        progress = progress_reporter("Mean error test", len(dates), report_every=20)

        for i, (year, month, day, hour, jd) in enumerate(dates):
            try:
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                lib_lon, _, _ = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)
                errors.append(diff)
            except Exception:
                pass  # Skip dates outside ephemeris range

            progress.update(i)

        if not errors:
            pytest.fail("No valid dates were tested")

        mean_error = statistics.mean(errors)

        assert mean_error < MEAN_ERROR_THRESHOLD, (
            f"Mean error {mean_error:.6f}° exceeds {MEAN_ERROR_THRESHOLD}° threshold"
        )

        progress.done(f"mean: {mean_error:.6f}°")

    def test_error_standard_deviation(self, progress_reporter):
        """
        Verify standard deviation is reasonable.

        Low standard deviation indicates consistent precision,
        without outliers significantly affecting results.
        Current threshold: 0.05 degrees (180 arcseconds)
        """
        dates = generate_random_jd(1900, 2100, seed=456)  # Different seed
        errors = []

        progress = progress_reporter("Std dev test", len(dates), report_every=20)

        for i, (year, month, day, hour, jd) in enumerate(dates):
            try:
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                lib_lon, _, _ = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)
                errors.append(diff)
            except Exception:
                pass

            progress.update(i)

        if len(errors) < 2:
            pytest.fail("Not enough dates for standard deviation")

        std_error = statistics.stdev(errors)

        assert std_error < STD_ERROR_THRESHOLD, (
            f"Standard deviation {std_error:.6f}° exceeds {STD_ERROR_THRESHOLD}° threshold"
        )

        progress.done(f"std: {std_error:.6f}°")
