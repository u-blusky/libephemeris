"""
Comprehensive Performance Benchmark: LibEphemeris vs pyswisseph

This script measures and compares the execution speed of various calculations
between libephemeris (pure Python) and pyswisseph (C extension wrapper).

Benchmarks cover:
- Planetary positions (Sun through Pluto)
- House calculations (all 19 systems)
- Lunar points (Mean/True Node, Lilith)
- Minor bodies (Chiron, Ceres, etc.)
- Utility functions (Julian Day, ayanamsa, coordinate transforms)

Each calculation type is run 10000+ times to get statistically meaningful results.
Reports mean, median, min, max execution times and performance ratios.
"""

import statistics
import sys
import time
from dataclasses import dataclass
from typing import Callable, Optional

import swisseph as swe

import libephemeris as ephem
from libephemeris.constants import (
    SE_CERES,
    SE_CHIRON,
    SE_JUNO,
    SE_JUPITER,
    SE_MARS,
    SE_MEAN_APOG,
    SE_MEAN_NODE,
    SE_MERCURY,
    SE_MOON,
    SE_NEPTUNE,
    SE_OSCU_APOG,
    SE_PALLAS,
    SE_PLUTO,
    SE_SATURN,
    SE_SUN,
    SE_TRUE_NODE,
    SE_URANUS,
    SE_VENUS,
    SE_VESTA,
    SEFLG_SIDEREAL,
    SEFLG_SPEED,
    SEFLG_SWIEPH,
)

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import print_header, print_section

# ============================================================================
# BENCHMARK CONFIGURATION
# ============================================================================

# Default number of iterations for benchmarks
DEFAULT_ITERATIONS = 10000

# Minimum iterations for slower operations
MIN_ITERATIONS = 1000

# Test dates (Julian Days) for benchmarks
TEST_JDS = [
    2451545.0,  # J2000.0 (2000-01-01 12:00)
    2460310.0,  # 2024-01-15 12:00
    2440587.5,  # 1970-01-01 00:00 (Unix epoch)
    2415020.5,  # 1900-01-01 00:00
]

# Test locations for house calculations
TEST_LOCATIONS = [
    (41.9028, 12.4964),  # Rome
    (40.7128, -74.0060),  # New York
    (35.6762, 139.6503),  # Tokyo
    (0.0, 0.0),  # Equator
]

# House systems to benchmark (subset for performance)
BENCHMARK_HOUSE_SYSTEMS = ["P", "K", "R", "C", "E", "W", "M", "O"]


# ============================================================================
# BENCHMARK RESULT CLASSES
# ============================================================================


@dataclass
class BenchmarkResult:
    """Stores results from a single benchmark run."""

    name: str
    iterations: int
    times_ns: list  # All individual timing measurements in nanoseconds
    category: str = "General"

    @property
    def mean_ns(self) -> float:
        """Mean execution time in nanoseconds."""
        return statistics.mean(self.times_ns) if self.times_ns else 0.0

    @property
    def median_ns(self) -> float:
        """Median execution time in nanoseconds."""
        return statistics.median(self.times_ns) if self.times_ns else 0.0

    @property
    def min_ns(self) -> float:
        """Minimum execution time in nanoseconds."""
        return min(self.times_ns) if self.times_ns else 0.0

    @property
    def max_ns(self) -> float:
        """Maximum execution time in nanoseconds."""
        return max(self.times_ns) if self.times_ns else 0.0

    @property
    def stdev_ns(self) -> float:
        """Standard deviation in nanoseconds."""
        if len(self.times_ns) > 1:
            return statistics.stdev(self.times_ns)
        return 0.0

    @property
    def total_ms(self) -> float:
        """Total execution time in milliseconds."""
        return sum(self.times_ns) / 1_000_000

    def format_time(self, ns: float) -> str:
        """Format nanoseconds to appropriate unit."""
        if ns >= 1_000_000:
            return f"{ns / 1_000_000:.3f} ms"
        elif ns >= 1_000:
            return f"{ns / 1_000:.3f} us"
        else:
            return f"{ns:.1f} ns"


@dataclass
class ComparisonResult:
    """Compares benchmark results between libephemeris and pyswisseph."""

    name: str
    libephem: BenchmarkResult
    pyswisseph: BenchmarkResult
    category: str = "General"

    @property
    def ratio(self) -> float:
        """
        Ratio of libephemeris to pyswisseph execution time.
        >1 means libephemeris is slower, <1 means libephemeris is faster.
        """
        if self.pyswisseph.mean_ns > 0:
            return self.libephem.mean_ns / self.pyswisseph.mean_ns
        return float("inf")

    @property
    def speedup(self) -> str:
        """Human-readable speedup description."""
        r = self.ratio
        if r < 1:
            return f"{1 / r:.2f}x faster"
        elif r > 1:
            return f"{r:.2f}x slower"
        else:
            return "same speed"


# ============================================================================
# BENCHMARK UTILITIES
# ============================================================================


def run_benchmark(
    name: str,
    func: Callable,
    iterations: int = DEFAULT_ITERATIONS,
    category: str = "General",
    warmup: int = 100,
) -> BenchmarkResult:
    """
    Run a benchmark for a given function.

    Args:
        name: Name of the benchmark
        func: Zero-argument callable to benchmark
        iterations: Number of iterations to run
        category: Category for grouping results
        warmup: Number of warmup iterations before timing

    Returns:
        BenchmarkResult with timing statistics
    """
    # Warmup phase
    for _ in range(warmup):
        func()

    # Actual benchmark
    times = []
    for _ in range(iterations):
        start = time.perf_counter_ns()
        func()
        end = time.perf_counter_ns()
        times.append(end - start)

    return BenchmarkResult(
        name=name, iterations=iterations, times_ns=times, category=category
    )


def compare_benchmark(
    name: str,
    libephem_func: Callable,
    pyswisseph_func: Callable,
    iterations: int = DEFAULT_ITERATIONS,
    category: str = "General",
) -> ComparisonResult:
    """
    Compare benchmark results between libephemeris and pyswisseph.

    Args:
        name: Name of the benchmark
        libephem_func: Zero-argument callable for libephemeris
        pyswisseph_func: Zero-argument callable for pyswisseph
        iterations: Number of iterations
        category: Category for grouping

    Returns:
        ComparisonResult with both benchmarks and ratio
    """
    libephem_result = run_benchmark(
        name=f"libephem:{name}",
        func=libephem_func,
        iterations=iterations,
        category=category,
    )

    pyswisseph_result = run_benchmark(
        name=f"pyswisseph:{name}",
        func=pyswisseph_func,
        iterations=iterations,
        category=category,
    )

    return ComparisonResult(
        name=name,
        libephem=libephem_result,
        pyswisseph=pyswisseph_result,
        category=category,
    )


# ============================================================================
# PLANETARY POSITION BENCHMARKS
# ============================================================================


def benchmark_planets(iterations: int = DEFAULT_ITERATIONS) -> list:
    """Benchmark planetary position calculations."""
    results = []
    jd = TEST_JDS[0]
    flags = SEFLG_SWIEPH | SEFLG_SPEED

    planets = [
        (SE_SUN, swe.SUN, "Sun"),
        (SE_MOON, swe.MOON, "Moon"),
        (SE_MERCURY, swe.MERCURY, "Mercury"),
        (SE_VENUS, swe.VENUS, "Venus"),
        (SE_MARS, swe.MARS, "Mars"),
        (SE_JUPITER, swe.JUPITER, "Jupiter"),
        (SE_SATURN, swe.SATURN, "Saturn"),
        (SE_URANUS, swe.URANUS, "Uranus"),
        (SE_NEPTUNE, swe.NEPTUNE, "Neptune"),
        (SE_PLUTO, swe.PLUTO, "Pluto"),
    ]

    for planet_py, planet_swe, name in planets:
        result = compare_benchmark(
            name=f"Planet: {name}",
            libephem_func=lambda p=planet_py: ephem.swe_calc_ut(jd, p, flags),
            pyswisseph_func=lambda p=planet_swe: swe.calc_ut(jd, p, flags),
            iterations=iterations,
            category="Planets",
        )
        results.append(result)

    # Also benchmark all planets in a single call sequence
    def calc_all_planets_libephem():
        for planet_py, _, _ in planets:
            ephem.swe_calc_ut(jd, planet_py, flags)

    def calc_all_planets_pyswisseph():
        for _, planet_swe, _ in planets:
            swe.calc_ut(jd, planet_swe, flags)

    results.append(
        compare_benchmark(
            name="All Planets (batch)",
            libephem_func=calc_all_planets_libephem,
            pyswisseph_func=calc_all_planets_pyswisseph,
            iterations=iterations // 10,
            category="Planets",
        )
    )

    return results


# ============================================================================
# HOUSE CALCULATION BENCHMARKS
# ============================================================================


def benchmark_houses(iterations: int = DEFAULT_ITERATIONS) -> list:
    """Benchmark house calculations."""
    results = []
    jd = TEST_JDS[0]
    lat, lon = TEST_LOCATIONS[0]

    for hsys in BENCHMARK_HOUSE_SYSTEMS:
        hsys_name = {
            "P": "Placidus",
            "K": "Koch",
            "R": "Regiomontanus",
            "C": "Campanus",
            "E": "Equal",
            "W": "Whole Sign",
            "M": "Morinus",
            "O": "Porphyry",
        }.get(hsys, hsys)

        result = compare_benchmark(
            name=f"Houses: {hsys_name}",
            libephem_func=lambda h=hsys: ephem.swe_houses(jd, lat, lon, h),
            pyswisseph_func=lambda h=hsys: swe.houses(jd, lat, lon, h.encode("ascii")),
            iterations=iterations,
            category="Houses",
        )
        results.append(result)

    return results


# ============================================================================
# LUNAR POINT BENCHMARKS
# ============================================================================


def benchmark_lunar_points(iterations: int = DEFAULT_ITERATIONS) -> list:
    """Benchmark lunar point calculations (nodes, Lilith)."""
    results = []
    jd = TEST_JDS[0]
    flags = SEFLG_SWIEPH | SEFLG_SPEED

    lunar_points = [
        (SE_MEAN_NODE, swe.MEAN_NODE, "Mean Node"),
        (SE_TRUE_NODE, swe.TRUE_NODE, "True Node"),
        (SE_MEAN_APOG, swe.MEAN_APOG, "Mean Lilith"),
        (SE_OSCU_APOG, swe.OSCU_APOG, "True Lilith"),
    ]

    for point_py, point_swe, name in lunar_points:
        result = compare_benchmark(
            name=f"Lunar: {name}",
            libephem_func=lambda p=point_py: ephem.swe_calc_ut(jd, p, flags),
            pyswisseph_func=lambda p=point_swe: swe.calc_ut(jd, p, flags),
            iterations=iterations,
            category="Lunar Points",
        )
        results.append(result)

    return results


# ============================================================================
# MINOR BODY BENCHMARKS
# ============================================================================


def benchmark_minor_bodies(iterations: int = MIN_ITERATIONS) -> list:
    """Benchmark minor body calculations."""
    results = []
    jd = TEST_JDS[0]
    flags = SEFLG_SWIEPH

    # Only include bodies that are commonly available
    minor_bodies = [
        (SE_CHIRON, swe.CHIRON, "Chiron"),
        (SE_CERES, swe.CERES, "Ceres"),
        (SE_PALLAS, swe.PALLAS, "Pallas"),
        (SE_JUNO, swe.JUNO, "Juno"),
        (SE_VESTA, swe.VESTA, "Vesta"),
    ]

    for body_py, body_swe, name in minor_bodies:
        # Wrap in try/except for pyswisseph since some asteroid files may be missing
        try:
            swe.calc_ut(jd, body_swe, flags)
            pyswisseph_available = True
        except Exception:
            pyswisseph_available = False

        if pyswisseph_available:
            result = compare_benchmark(
                name=f"Minor: {name}",
                libephem_func=lambda b=body_py: ephem.swe_calc_ut(jd, b, flags),
                pyswisseph_func=lambda b=body_swe: swe.calc_ut(jd, b, flags),
                iterations=iterations,
                category="Minor Bodies",
            )
            results.append(result)

    return results


# ============================================================================
# UTILITY FUNCTION BENCHMARKS
# ============================================================================


def benchmark_utilities(iterations: int = DEFAULT_ITERATIONS) -> list:
    """Benchmark utility functions."""
    results = []
    jd = TEST_JDS[0]

    # Julian Day calculation
    results.append(
        compare_benchmark(
            name="Julian Day (julday)",
            libephem_func=lambda: ephem.swe_julday(2000, 1, 1, 12.0),
            pyswisseph_func=lambda: swe.julday(2000, 1, 1, 12.0),
            iterations=iterations,
            category="Utilities",
        )
    )

    # Reverse Julian Day
    results.append(
        compare_benchmark(
            name="Reverse Julian Day (revjul)",
            libephem_func=lambda: ephem.swe_revjul(jd),
            pyswisseph_func=lambda: swe.revjul(jd),
            iterations=iterations,
            category="Utilities",
        )
    )

    # Delta T
    results.append(
        compare_benchmark(
            name="Delta T (deltat)",
            libephem_func=lambda: ephem.swe_deltat(jd),
            pyswisseph_func=lambda: swe.deltat(jd),
            iterations=iterations,
            category="Utilities",
        )
    )

    # Ayanamsa
    ephem.swe_set_sid_mode(1)  # Lahiri
    swe.set_sid_mode(1)
    results.append(
        compare_benchmark(
            name="Ayanamsa (Lahiri)",
            libephem_func=lambda: ephem.swe_get_ayanamsa_ut(jd),
            pyswisseph_func=lambda: swe.get_ayanamsa_ut(jd),
            iterations=iterations,
            category="Utilities",
        )
    )

    # Sidereal time - note: libephemeris sidtime0 requires obliquity and nutation
    # while pyswisseph sidtime uses internal nutation calculation
    # We use sidtime0 with typical values for comparison
    obliquity = 23.4393  # approximate obliquity
    nutation = 0.0  # approximate nutation in longitude
    results.append(
        compare_benchmark(
            name="Sidereal Time (sidtime0)",
            libephem_func=lambda: ephem.sidtime0(jd, obliquity, nutation),
            pyswisseph_func=lambda: swe.sidtime(jd),
            iterations=iterations,
            category="Utilities",
        )
    )

    return results


# ============================================================================
# SIDEREAL CALCULATION BENCHMARKS
# ============================================================================


def benchmark_sidereal(iterations: int = DEFAULT_ITERATIONS) -> list:
    """Benchmark sidereal (vedic) calculations."""
    results = []
    jd = TEST_JDS[0]
    flags = SEFLG_SWIEPH | SEFLG_SIDEREAL | SEFLG_SPEED

    # Set ayanamsa mode (Lahiri)
    ephem.swe_set_sid_mode(1)
    swe.set_sid_mode(1)

    planets = [
        (SE_SUN, swe.SUN, "Sun"),
        (SE_MOON, swe.MOON, "Moon"),
        (SE_MARS, swe.MARS, "Mars"),
    ]

    for planet_py, planet_swe, name in planets:
        result = compare_benchmark(
            name=f"Sidereal: {name}",
            libephem_func=lambda p=planet_py: ephem.swe_calc_ut(jd, p, flags),
            pyswisseph_func=lambda p=planet_swe: swe.calc_ut(jd, p, flags),
            iterations=iterations,
            category="Sidereal",
        )
        results.append(result)

    return results


# ============================================================================
# COORDINATE TRANSFORMATION BENCHMARKS
# ============================================================================


def benchmark_coordinates(iterations: int = DEFAULT_ITERATIONS) -> list:
    """Benchmark coordinate transformation functions."""
    results = []

    # azalt (ecliptic to horizontal)
    # Note: libephemeris azalt has different parameter order than pyswisseph
    # libephem: azalt(jd, calc_flag, lat, lon, altitude, pressure, temperature, coord)
    # pyswisseph: azalt(jd, calc_flag, (lon, lat, alt), pressure, temperature, coord)
    jd = TEST_JDS[0]
    lat, lon = TEST_LOCATIONS[0]
    xin = (123.45, 1.23, 1.0)  # ecliptic lon/lat/distance

    results.append(
        compare_benchmark(
            name="Coordinate: azalt",
            libephem_func=lambda: ephem.azalt(jd, 0, lat, lon, 0, 0.0, 0.0, xin),
            pyswisseph_func=lambda: swe.azalt(jd, 0, (lon, lat, 0), 0.0, 0.0, xin),
            iterations=iterations,
            category="Coordinates",
        )
    )

    # cotrans (coordinate transformation)
    eps = 23.4392911  # obliquity
    results.append(
        compare_benchmark(
            name="Coordinate: cotrans",
            libephem_func=lambda: ephem.cotrans((123.45, 1.23, 1.0), eps),
            pyswisseph_func=lambda: swe.cotrans((123.45, 1.23, 1.0), eps),
            iterations=iterations,
            category="Coordinates",
        )
    )

    return results


# ============================================================================
# REPORT GENERATION
# ============================================================================


def print_results_table(results: list, title: str):
    """Print benchmark results as a formatted table."""
    if not results:
        return

    print_section(title)
    print(
        f"{'Benchmark':<30} {'LibEphem':>12} {'pyswisseph':>12} {'Ratio':>10} {'Status':<15}"
    )
    print("-" * 85)

    for r in results:
        ratio_str = f"{r.ratio:.2f}x"
        status = r.speedup

        # Color coding for terminal (if supported)
        if r.ratio < 1:
            status_marker = "+"  # faster
        elif r.ratio > 2:
            status_marker = "!"  # significantly slower
        elif r.ratio > 1:
            status_marker = "-"  # slower
        else:
            status_marker = "="  # same

        print(
            f"{r.name:<30} "
            f"{r.libephem.format_time(r.libephem.mean_ns):>12} "
            f"{r.pyswisseph.format_time(r.pyswisseph.mean_ns):>12} "
            f"{ratio_str:>10} "
            f"{status_marker} {status:<14}"
        )


def print_summary_statistics(all_results: list):
    """Print summary statistics across all benchmarks."""
    print_section("SUMMARY STATISTICS")

    if not all_results:
        print("No benchmark results to summarize.")
        return

    ratios = [r.ratio for r in all_results if r.ratio != float("inf")]

    if not ratios:
        print("No valid ratio data available.")
        return

    mean_ratio = statistics.mean(ratios)
    median_ratio = statistics.median(ratios)
    min_ratio = min(ratios)
    max_ratio = max(ratios)

    faster_count = sum(1 for r in ratios if r < 1)
    same_count = sum(1 for r in ratios if 0.9 <= r <= 1.1)
    slower_count = sum(1 for r in ratios if r > 1.1)

    print(f"Total benchmarks:    {len(all_results)}")
    print()
    print("Performance Ratio (libephemeris / pyswisseph):")
    print(f"  Mean:              {mean_ratio:.2f}x")
    print(f"  Median:            {median_ratio:.2f}x")
    print(f"  Min (fastest):     {min_ratio:.2f}x")
    print(f"  Max (slowest):     {max_ratio:.2f}x")
    print()
    print("Distribution:")
    print(
        f"  Faster (<1x):      {faster_count} ({100 * faster_count / len(ratios):.1f}%)"
    )
    print(f"  Same (~1x):        {same_count} ({100 * same_count / len(ratios):.1f}%)")
    print(
        f"  Slower (>1x):      {slower_count} ({100 * slower_count / len(ratios):.1f}%)"
    )

    # Category breakdown
    print()
    print("By Category:")
    categories = {}
    for r in all_results:
        if r.category not in categories:
            categories[r.category] = []
        if r.ratio != float("inf"):
            categories[r.category].append(r.ratio)

    for cat, cat_ratios in sorted(categories.items()):
        if cat_ratios:
            cat_mean = statistics.mean(cat_ratios)
            print(f"  {cat:<20} Mean ratio: {cat_mean:.2f}x")


def print_detailed_stats(results: list, title: str):
    """Print detailed statistics for each benchmark."""
    if not results:
        return

    print_section(f"DETAILED STATISTICS: {title}")
    print(
        f"{'Benchmark':<25} {'Mean':>10} {'Median':>10} {'Min':>10} {'Max':>10} {'StdDev':>10}"
    )
    print("-" * 85)

    for r in results:
        # LibEphemeris stats
        print(
            f"{'libephem:':<25} "
            f"{r.libephem.format_time(r.libephem.mean_ns):>10} "
            f"{r.libephem.format_time(r.libephem.median_ns):>10} "
            f"{r.libephem.format_time(r.libephem.min_ns):>10} "
            f"{r.libephem.format_time(r.libephem.max_ns):>10} "
            f"{r.libephem.format_time(r.libephem.stdev_ns):>10}"
        )
        # pyswisseph stats
        print(
            f"{'  pyswisseph:':<25} "
            f"{r.pyswisseph.format_time(r.pyswisseph.mean_ns):>10} "
            f"{r.pyswisseph.format_time(r.pyswisseph.median_ns):>10} "
            f"{r.pyswisseph.format_time(r.pyswisseph.min_ns):>10} "
            f"{r.pyswisseph.format_time(r.pyswisseph.max_ns):>10} "
            f"{r.pyswisseph.format_time(r.pyswisseph.stdev_ns):>10}"
        )
        print()


def generate_performance_report(all_results: list, detailed: bool = False) -> str:
    """Generate a complete performance report as a string."""
    lines = []

    lines.append("=" * 80)
    lines.append("LIBEPHEMERIS vs PYSWISSEPH PERFORMANCE REPORT")
    lines.append("=" * 80)
    lines.append("")

    # Summary table
    lines.append("PERFORMANCE COMPARISON (Ratio = libephemeris / pyswisseph)")
    lines.append("-" * 80)
    lines.append(f"{'Category':<20} {'Benchmark':<25} {'Ratio':>10} {'Status':<15}")
    lines.append("-" * 80)

    for r in all_results:
        ratio_str = f"{r.ratio:.2f}x" if r.ratio != float("inf") else "N/A"
        lines.append(f"{r.category:<20} {r.name:<25} {ratio_str:>10} {r.speedup:<15}")

    lines.append("")

    # Summary statistics
    ratios = [r.ratio for r in all_results if r.ratio != float("inf")]
    if ratios:
        lines.append("SUMMARY:")
        lines.append(f"  Total benchmarks: {len(all_results)}")
        lines.append(f"  Mean ratio:       {statistics.mean(ratios):.2f}x")
        lines.append(f"  Median ratio:     {statistics.median(ratios):.2f}x")

    return "\n".join(lines)


# ============================================================================
# PUBLIC API FOR TESTS
# ============================================================================


def run_benchmarks(
    iterations: int = DEFAULT_ITERATIONS,
    categories: Optional[list] = None,
    verbose: bool = True,
    detailed: bool = False,
) -> dict:
    """
    Run all benchmarks and return results.

    Args:
        iterations: Number of iterations per benchmark
        categories: List of categories to run, or None for all
        verbose: Print results to stdout
        detailed: Include detailed statistics

    Returns:
        Dictionary with benchmark results and statistics
    """
    all_categories = {
        "planets": benchmark_planets,
        "houses": benchmark_houses,
        "lunar": benchmark_lunar_points,
        "minor": benchmark_minor_bodies,
        "utilities": benchmark_utilities,
        "sidereal": benchmark_sidereal,
        "coordinates": benchmark_coordinates,
    }

    if categories is None:
        categories = list(all_categories.keys())

    all_results = []

    if verbose:
        print_header("LIBEPHEMERIS vs PYSWISSEPH PERFORMANCE BENCHMARK")
        print(f"Iterations per benchmark: {iterations}")
        print()

    for category in categories:
        if category in all_categories:
            if verbose:
                print(f"Running {category} benchmarks...")
            results = all_categories[category](iterations)
            all_results.extend(results)

            if verbose:
                print_results_table(results, category.upper())
                if detailed:
                    print_detailed_stats(results, category.upper())

    if verbose:
        print_summary_statistics(all_results)

    # Calculate summary statistics
    ratios = [r.ratio for r in all_results if r.ratio != float("inf")]
    summary = {
        "total_benchmarks": len(all_results),
        "mean_ratio": statistics.mean(ratios) if ratios else 0.0,
        "median_ratio": statistics.median(ratios) if ratios else 0.0,
        "min_ratio": min(ratios) if ratios else 0.0,
        "max_ratio": max(ratios) if ratios else 0.0,
    }

    return {
        "results": all_results,
        "summary": summary,
        "report": generate_performance_report(all_results, detailed),
    }


def get_quick_benchmark_stats(iterations: int = 1000) -> dict:
    """
    Run a quick benchmark with fewer iterations for testing.

    Returns dict with mean_ratio, median_ratio, etc.
    """
    result = run_benchmarks(
        iterations=iterations, categories=["planets", "utilities"], verbose=False
    )
    return result["summary"]


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_benchmark.py [OPTIONS]")
    print()
    print("Options:")
    print("  -n, --iterations N   Number of iterations (default: 10000)")
    print("  -d, --detailed       Show detailed statistics")
    print("  -q, --quick          Quick mode (1000 iterations)")
    print("  --category CAT       Run only specific category")
    print("                       (planets, houses, lunar, minor, utilities,")
    print("                        sidereal, coordinates)")
    print("  -h, --help           Show this help message")
    print()
    print("Examples:")
    print("  python compare_benchmark.py              # Full benchmark")
    print("  python compare_benchmark.py -q           # Quick benchmark")
    print("  python compare_benchmark.py --category planets")
    print()


def main():
    """Main entry point."""
    # Parse arguments
    args = sys.argv[1:]

    if "-h" in args or "--help" in args:
        print_help()
        sys.exit(0)

    iterations = DEFAULT_ITERATIONS
    detailed = False
    categories = None

    if "-q" in args or "--quick" in args:
        iterations = 1000

    if "-n" in args:
        try:
            idx = args.index("-n")
            iterations = int(args[idx + 1])
        except (IndexError, ValueError):
            print("Error: -n requires an integer argument")
            sys.exit(1)

    if "--iterations" in args:
        try:
            idx = args.index("--iterations")
            iterations = int(args[idx + 1])
        except (IndexError, ValueError):
            print("Error: --iterations requires an integer argument")
            sys.exit(1)

    if "-d" in args or "--detailed" in args:
        detailed = True

    if "--category" in args:
        try:
            idx = args.index("--category")
            categories = [args[idx + 1]]
        except IndexError:
            print("Error: --category requires an argument")
            sys.exit(1)

    # Run benchmarks
    result = run_benchmarks(
        iterations=iterations, categories=categories, verbose=True, detailed=detailed
    )

    # Print final report
    print()
    print("=" * 80)
    print("BENCHMARK COMPLETE")
    print("=" * 80)
    summary = result["summary"]
    print(f"Total benchmarks run: {summary['total_benchmarks']}")
    print(f"Average performance ratio: {summary['mean_ratio']:.2f}x")
    print()

    if summary["mean_ratio"] > 1:
        print(
            f"LibEphemeris is on average {summary['mean_ratio']:.1f}x slower than pyswisseph."
        )
        print(
            "This is expected as libephemeris is pure Python while pyswisseph is a C extension."
        )
    else:
        print(
            f"LibEphemeris is on average {1 / summary['mean_ratio']:.1f}x faster than pyswisseph."
        )


if __name__ == "__main__":
    main()
