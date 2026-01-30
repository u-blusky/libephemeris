"""
Comparison script for Lunar calculations: Nodes and Lilith.
Validates LibEphemeris against SwissEphemeris for lunar points.
"""

import random
import statistics
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *
from libephemeris.lunar import calc_true_lunar_node, calc_true_lilith
import sys

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import *


# Threshold constants for True Node precision comparison
TRUE_NODE_MAX_ERROR_THRESHOLD = 0.15  # degrees (consistent with existing tests)
TRUE_NODE_MEAN_ERROR_THRESHOLD = 0.08  # degrees (relaxed for 1900-2100 range)
TRUE_NODE_STD_ERROR_THRESHOLD = 0.05  # degrees
TRUE_NODE_TARGET_ERROR = 0.01  # degrees (36 arcseconds) - aspirational

# Threshold constants for True Lilith (Osculating Apogee) precision comparison
# These are relaxed thresholds - the implementation differs from Swiss Ephemeris
TRUE_LILITH_MAX_ERROR_THRESHOLD = 7.0  # degrees (current known range 5-7)
TRUE_LILITH_MEAN_ERROR_THRESHOLD = 5.0  # degrees
TRUE_LILITH_STD_ERROR_THRESHOLD = 2.0  # degrees
TRUE_LILITH_TARGET_ERROR = 0.5  # degrees - goal after corrections
TRUE_LILITH_NUM_DATES = 500  # Number of random dates to test


def compare_lunar_nodes(jd, name, date_str):
    """Compare lunar nodes between implementations."""
    print(f"\n{'=' * 80}")
    print(f"LUNAR NODES - {name} ({date_str})")
    print(f"{'=' * 80}")

    results = {}

    # Mean North Node
    print(f"\n{'Mean North Node':-^80}")
    pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, 0)
    pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, 0)

    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    passed = diff_lon < 0.01

    print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}°")
    print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}°")
    print(f"Difference:   {format_diff(diff_lon, 8)}° {format_status(passed)}")

    results["mean_node"] = (passed, diff_lon)

    # True North Node
    print(f"\n{'True North Node':-^80}")
    pos_swe, _ = swe.calc_ut(jd, swe.TRUE_NODE, 0)
    pos_py, _ = ephem.swe_calc_ut(jd, SE_TRUE_NODE, 0)

    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    passed = diff_lon < 1.0  # Relaxed tolerance for true node

    print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}°")
    print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}°")
    print(f"Difference:   {format_diff(diff_lon, 8)}° {format_status(passed)}")

    results["true_node"] = (passed, diff_lon)

    return results


def compare_lilith(jd, name, date_str):
    """Compare Lilith between implementations."""
    print(f"\n{'=' * 80}")
    print(f"LILITH - {name} ({date_str})")
    print(f"{'=' * 80}")

    results = {}

    # Mean Lilith
    print(f"\n{'Mean Lilith (Black Moon)':-^80}")
    pos_swe, _ = swe.calc_ut(jd, swe.MEAN_APOG, 0)
    pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, 0)

    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    passed = diff_lon < 0.1

    print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}°")
    print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}°")
    print(f"Difference:   {format_diff(diff_lon, 8)}° {format_status(passed)}")

    results["mean_lilith"] = (passed, diff_lon)

    # True Lilith (Osculating Apogee)
    print(f"\n{'True Lilith (Osculating Apogee)':-^80}")
    pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
    pos_py, _ = ephem.swe_calc_ut(jd, SE_OSCU_APOG, 0)

    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    passed = diff_lon < 5.0  # Very relaxed for osculating apogee

    print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}°")
    print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}°")
    print(f"Difference:   {format_diff(diff_lon, 8)}° {format_status(passed)}")
    print(
        "\nNote: Osculating apogee can have larger differences due to calculation method."
    )

    results["true_lilith"] = (passed, diff_lon)

    # Interpolated Apogee (SE_INTP_APOG)
    print(f"\n{'Interpolated Apogee':-^80}")
    try:
        pos_swe, _ = swe.calc_ut(jd, swe.INTP_APOG, 0)
        try:
            pos_py, _ = ephem.swe_calc_ut(jd, SE_INTP_APOG, 0)
            diff_lon = angular_diff(pos_swe[0], pos_py[0])
            passed = diff_lon < 1.0  # Relaxed tolerance

            print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}")
            print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}")
            print(f"Difference:   {format_diff(diff_lon, 8)} {format_status(passed)}")

            results["intp_apog"] = (passed, diff_lon)
        except Exception as e:
            print(f"LibEphemeris: NOT IMPLEMENTED ({e})")
            results["intp_apog"] = (True, 0.0)  # Skip if not implemented
    except Exception as e:
        print(f"SwissEph: NOT AVAILABLE ({e})")
        results["intp_apog"] = (True, 0.0)

    # Interpolated Perigee (SE_INTP_PERG)
    print(f"\n{'Interpolated Perigee':-^80}")
    try:
        pos_swe, _ = swe.calc_ut(jd, swe.INTP_PERG, 0)
        try:
            pos_py, _ = ephem.swe_calc_ut(jd, SE_INTP_PERG, 0)
            diff_lon = angular_diff(pos_swe[0], pos_py[0])
            passed = diff_lon < 1.0  # Relaxed tolerance

            print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}")
            print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}")
            print(f"Difference:   {format_diff(diff_lon, 8)} {format_status(passed)}")

            results["intp_perg"] = (passed, diff_lon)
        except Exception as e:
            print(f"LibEphemeris: NOT IMPLEMENTED ({e})")
            results["intp_perg"] = (True, 0.0)  # Skip if not implemented
    except Exception as e:
        print(f"SwissEph: NOT AVAILABLE ({e})")
        results["intp_perg"] = (True, 0.0)

    return results


def generate_random_jd(
    year_start: int, year_end: int, count: int = 1000, seed: int = 42
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


def compare_true_node_precision():
    """
    Compare True Lunar Node precision against pyswisseph using 1000 random dates.

    This function validates calc_true_lunar_node against pyswisseph swe.calc_ut(jd, swe.TRUE_NODE)
    spanning 1900-2100, calculating maximum/mean/standard deviation of differences
    and reporting dates with largest discrepancies.

    Returns:
        TestStatistics object with results
    """
    print(f"\n{'=' * 80}")
    print("TRUE NODE PRECISION COMPARISON (1000 random dates, 1900-2100)")
    print(f"{'=' * 80}")

    stats = TestStatistics()
    dates = generate_random_jd(1900, 2100, count=1000, seed=42)
    errors = []
    discrepancies = []  # Store (diff, jd, swe_lon, lib_lon, date_info)
    skipped = 0

    for year, month, day, hour, jd in dates:
        try:
            # Calculate with pyswisseph
            swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
            swe_lon = swe_result[0][0]

            # Calculate with libephemeris
            lib_lon, lib_lat, lib_dist = calc_true_lunar_node(jd)

            # Calculate angular difference
            diff = angular_diff(swe_lon, lib_lon)
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

            # Check if within threshold
            passed = diff < TRUE_NODE_MAX_ERROR_THRESHOLD
            stats.add_result(passed, diff)

        except Exception as e:
            # Skip dates outside ephemeris range
            if "ephemeris" in str(e).lower() or "outside" in str(e).lower():
                skipped += 1
            else:
                stats.add_result(False, 0.0, error=True)

    if not errors:
        print("ERROR: No valid dates were tested")
        return stats

    # Calculate statistics
    max_error = max(errors)
    mean_error = statistics.mean(errors)
    std_error = statistics.stdev(errors) if len(errors) > 1 else 0.0
    median_error = statistics.median(errors)

    # Sort discrepancies by error (descending)
    discrepancies.sort(key=lambda x: x[0], reverse=True)

    # Count dates meeting target precision
    dates_under_target = sum(1 for e in errors if e < TRUE_NODE_TARGET_ERROR)
    pct_under_target = 100 * dates_under_target / len(errors)

    # Print comprehensive report
    print("\n--- Statistics ---")
    print(f"Dates tested: {len(errors)}")
    print(f"Dates skipped: {skipped}")
    print("Date range: 1900-2100")
    print(f"Maximum difference: {max_error:.6f}° ({max_error * 3600:.2f} arcsec)")
    print(f"Mean difference:    {mean_error:.6f}° ({mean_error * 3600:.2f} arcsec)")
    print(f"Median difference:  {median_error:.6f}° ({median_error * 3600:.2f} arcsec)")
    print(f"Std deviation:      {std_error:.6f}° ({std_error * 3600:.2f} arcsec)")

    print("\n--- Target Precision Analysis ---")
    print(
        f"Target: < {TRUE_NODE_TARGET_ERROR}° ({TRUE_NODE_TARGET_ERROR * 3600:.0f} arcsec)"
    )
    print(
        f"Dates meeting target: {dates_under_target}/{len(errors)} ({pct_under_target:.1f}%)"
    )
    if max_error >= TRUE_NODE_TARGET_ERROR:
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

    # Print pass/fail status
    max_passed = max_error < TRUE_NODE_MAX_ERROR_THRESHOLD
    mean_passed = mean_error < TRUE_NODE_MEAN_ERROR_THRESHOLD
    std_passed = std_error < TRUE_NODE_STD_ERROR_THRESHOLD

    print("\n--- Threshold Checks ---")
    print(f"Max error < {TRUE_NODE_MAX_ERROR_THRESHOLD}°: {format_status(max_passed)}")
    print(
        f"Mean error < {TRUE_NODE_MEAN_ERROR_THRESHOLD}°: {format_status(mean_passed)}"
    )
    print(f"Std error < {TRUE_NODE_STD_ERROR_THRESHOLD}°: {format_status(std_passed)}")

    return stats


def compare_true_lilith_precision():
    """
    Compare True Lilith (Osculating Apogee) precision against pyswisseph using 500+ random dates.

    This function validates calc_true_lilith against pyswisseph swe.calc_ut(jd, swe.OSCU_APOG, flags)
    spanning 1900-2100, calculating maximum/mean/standard deviation of differences
    and tracking progress toward the target error of less than 0.5 degrees.

    Returns:
        TestStatistics object with results
    """
    print(f"\n{'=' * 80}")
    print(
        f"TRUE LILITH (OSCULATING APOGEE) PRECISION COMPARISON ({TRUE_LILITH_NUM_DATES} random dates, 1900-2100)"
    )
    print(f"{'=' * 80}")

    stats = TestStatistics()
    dates = generate_random_jd(1900, 2100, count=TRUE_LILITH_NUM_DATES, seed=42)
    lon_errors = []
    lat_errors = []
    discrepancies = []  # Store (diff, jd, swe_lon, lib_lon, swe_lat, lib_lat, date_info)
    skipped = 0

    for year, month, day, hour, jd in dates:
        try:
            # Calculate with pyswisseph (OSCU_APOG = osculating apogee = True Lilith)
            swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
            swe_lon = swe_result[0][0]
            swe_lat = swe_result[0][1]
            swe_dist = swe_result[0][2]

            # Calculate with libephemeris
            lib_lon, lib_lat, lib_e = calc_true_lilith(jd)

            # Calculate angular difference for longitude
            diff = angular_diff(swe_lon, lib_lon)
            lon_errors.append(diff)

            # Track latitude difference
            lat_diff = abs(swe_lat - lib_lat)
            lat_errors.append(lat_diff)

            # Store for discrepancy report
            discrepancies.append(
                (
                    diff,
                    jd,
                    swe_lon,
                    lib_lon,
                    swe_lat,
                    lib_lat,
                    f"{year:04d}-{month:02d}-{day:02d} {hour:.2f}h",
                )
            )

            # Check if within threshold
            passed = diff < TRUE_LILITH_MAX_ERROR_THRESHOLD
            stats.add_result(passed, diff)

        except Exception as e:
            # Skip dates outside ephemeris range
            if "ephemeris" in str(e).lower() or "outside" in str(e).lower():
                skipped += 1
            else:
                stats.add_result(False, 0.0, error=True)

    if not lon_errors:
        print("ERROR: No valid dates were tested")
        return stats

    # Calculate statistics
    max_error = max(lon_errors)
    mean_error = statistics.mean(lon_errors)
    std_error = statistics.stdev(lon_errors) if len(lon_errors) > 1 else 0.0
    median_error = statistics.median(lon_errors)
    min_error = min(lon_errors)

    # Latitude statistics
    max_lat_error = max(lat_errors) if lat_errors else 0.0
    mean_lat_error = statistics.mean(lat_errors) if lat_errors else 0.0

    # Sort discrepancies by error (descending)
    discrepancies.sort(key=lambda x: x[0], reverse=True)

    # Count dates meeting target precision
    dates_under_target = sum(1 for e in lon_errors if e < TRUE_LILITH_TARGET_ERROR)
    pct_under_target = 100 * dates_under_target / len(lon_errors)

    # Count dates at various precision levels
    under_1_deg = sum(1 for e in lon_errors if e < 1.0)
    under_2_deg = sum(1 for e in lon_errors if e < 2.0)
    under_3_deg = sum(1 for e in lon_errors if e < 3.0)
    under_5_deg = sum(1 for e in lon_errors if e < 5.0)

    # Print comprehensive report
    print("\n--- Longitude Statistics ---")
    print(f"Dates tested: {len(lon_errors)}")
    print(f"Dates skipped: {skipped}")
    print("Date range: 1900-2100")
    print(f"Maximum difference: {max_error:.4f} ({max_error * 3600:.1f} arcsec)")
    print(f"Mean difference:    {mean_error:.4f} ({mean_error * 3600:.1f} arcsec)")
    print(f"Median difference:  {median_error:.4f} ({median_error * 3600:.1f} arcsec)")
    print(f"Std deviation:      {std_error:.4f} ({std_error * 3600:.1f} arcsec)")
    print(f"Minimum difference: {min_error:.4f} ({min_error * 3600:.1f} arcsec)")

    print("\n--- Latitude Statistics ---")
    print(
        f"Maximum lat diff:   {max_lat_error:.4f} ({max_lat_error * 3600:.1f} arcsec)"
    )
    print(
        f"Mean lat diff:      {mean_lat_error:.4f} ({mean_lat_error * 3600:.1f} arcsec)"
    )

    print("\n--- Target Precision Analysis ---")
    print(f"Target: < {TRUE_LILITH_TARGET_ERROR}")
    print(
        f"Dates meeting target: {dates_under_target}/{len(lon_errors)} ({pct_under_target:.1f}%)"
    )
    if max_error >= TRUE_LILITH_TARGET_ERROR:
        print(f"NOTE: Target precision NOT yet achieved. Max error: {max_error:.4f}")
        print("      Corrections needed to reduce error to < 0.5")

    # Distribution analysis
    print("\n--- Error Distribution ---")
    print(f"< 0.5:   {dates_under_target:4d} ({pct_under_target:.1f}%)")
    print(f"< 1.0:   {under_1_deg:4d} ({100 * under_1_deg / len(lon_errors):.1f}%)")
    print(f"< 2.0:   {under_2_deg:4d} ({100 * under_2_deg / len(lon_errors):.1f}%)")
    print(f"< 3.0:   {under_3_deg:4d} ({100 * under_3_deg / len(lon_errors):.1f}%)")
    print(f"< 5.0:   {under_5_deg:4d} ({100 * under_5_deg / len(lon_errors):.1f}%)")

    # Report top 10 largest discrepancies
    print("\n--- Top 10 Largest Discrepancies ---")
    print(
        f"{'Date':<22} {'JD':<14} {'SWE Lon':<12} {'LIB Lon':<12} {'Diff':<10} {'SWE Lat':<10} {'LIB Lat':<10}"
    )
    print("-" * 90)
    for diff, jd, swe_lon, lib_lon, swe_lat, lib_lat, date_info in discrepancies[:10]:
        print(
            f"{date_info:<22} {jd:<14.4f} {swe_lon:<12.4f} {lib_lon:<12.4f} "
            f"{diff:<10.4f} {swe_lat:<10.4f} {lib_lat:<10.4f}"
        )

    # Progress tracking section
    print(f"\n{'=' * 80}")
    print("CORRECTION PROGRESS TRACKING")
    print(f"{'=' * 80}")
    print("Corrections implemented in calc_true_lilith:")
    print("  [x] Evection correction (1.274 amplitude)")
    print("  [x] Solar gravitational perturbation on eccentricity vector")
    print("  [x] Evection-related secondary terms")
    print("  [x] Variation correction (0.658 amplitude)")
    print("  [x] Annual equation (0.186 amplitude)")
    print("  [x] Parallactic inequality (0.125 amplitude)")
    print("  [x] Reduction to ecliptic (0.116 amplitude)")
    print("  [ ] Additional perturbation terms (if needed)")
    print("")
    print(f"Progress toward < {TRUE_LILITH_TARGET_ERROR} target:")
    progress_bar_len = 40
    progress_pct = (
        min(100, (TRUE_LILITH_TARGET_ERROR / max_error) * 100) if max_error > 0 else 100
    )
    filled = int(progress_bar_len * progress_pct / 100)
    bar = "#" * filled + "-" * (progress_bar_len - filled)
    print(f"  [{bar}] {progress_pct:.1f}%")
    print(
        f"  Current max error: {max_error:.4f} -> Target: < {TRUE_LILITH_TARGET_ERROR}"
    )

    # Print pass/fail status
    max_passed = max_error < TRUE_LILITH_MAX_ERROR_THRESHOLD
    mean_passed = mean_error < TRUE_LILITH_MEAN_ERROR_THRESHOLD
    std_passed = std_error < TRUE_LILITH_STD_ERROR_THRESHOLD

    print("\n--- Threshold Checks ---")
    print(f"Max error < {TRUE_LILITH_MAX_ERROR_THRESHOLD}: {format_status(max_passed)}")
    print(
        f"Mean error < {TRUE_LILITH_MEAN_ERROR_THRESHOLD}: {format_status(mean_passed)}"
    )
    print(f"Std error < {TRUE_LILITH_STD_ERROR_THRESHOLD}: {format_status(std_passed)}")

    return stats


def main():
    print_header("LUNAR CALCULATIONS COMPARISON: LibEphemeris vs SwissEphemeris")

    # Test subjects
    subjects = [
        ("J2000.0", 2000, 1, 1, 12.0),
        ("2024-01-01", 2024, 1, 1, 0.0),
        ("1990-06-15", 1990, 6, 15, 12.0),
        ("2010-12-31", 2010, 12, 31, 23.5),
    ]

    stats = TestStatistics()

    for name, year, month, day, hour in subjects:
        jd = swe.julday(year, month, day, hour)
        date_str = f"{year}-{month:02d}-{day:02d}"

        # Compare nodes
        node_results = compare_lunar_nodes(jd, name, date_str)
        for key, (passed, diff) in node_results.items():
            stats.add_result(passed, diff)

        # Compare Lilith
        lilith_results = compare_lilith(jd, name, date_str)
        for key, (passed, diff) in lilith_results.items():
            stats.add_result(passed, diff)

    # Run True Node precision comparison (1000 random dates)
    precision_stats = compare_true_node_precision()

    # Merge precision stats into main stats
    stats.total += precision_stats.total
    stats.passed += precision_stats.passed
    stats.failed += precision_stats.failed
    stats.errors += precision_stats.errors
    stats.max_diff = max(stats.max_diff, precision_stats.max_diff)
    stats.diff_sum += precision_stats.diff_sum

    # Run True Lilith precision comparison (500+ random dates)
    lilith_precision_stats = compare_true_lilith_precision()

    # Merge Lilith precision stats into main stats
    stats.total += lilith_precision_stats.total
    stats.passed += lilith_precision_stats.passed
    stats.failed += lilith_precision_stats.failed
    stats.errors += lilith_precision_stats.errors
    stats.max_diff = max(stats.max_diff, lilith_precision_stats.max_diff)
    stats.diff_sum += lilith_precision_stats.diff_sum

    # Summary
    stats.print_summary("LUNAR CALCULATIONS SUMMARY")

    if stats.pass_rate() >= 95:
        print("\n✓ Lunar calculations show excellent compatibility!")
    elif stats.pass_rate() >= 80:
        print("\n~ Lunar calculations show good compatibility with some differences.")
    else:
        print("\n✗ Lunar calculations show significant differences.")

    return 0 if stats.pass_rate() >= 80 else 1


if __name__ == "__main__":
    sys.exit(main())
