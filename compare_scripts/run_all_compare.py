#!/usr/bin/env python3
"""
Master Runner for All Comparison Scripts.

Runs all comparison tests between pyswisseph and libephemeris
and produces a comprehensive report.

Usage:
    python run_all_compare.py [OPTIONS]

Options:
    -v, --verbose    Show detailed output for each test
    --quick          Run only quick tests (skip slow ones)
    --module NAME    Run only specific module
    --list           List available modules
    -h, --help       Show this help message
"""

import sys
import os
import subprocess
import time
from dataclasses import dataclass
from typing import List, Tuple, Optional


@dataclass
class ModuleResult:
    """Result from running a comparison module."""

    name: str
    passed: int
    total: int
    duration: float
    error: Optional[str] = None

    @property
    def pass_rate(self) -> float:
        return (self.passed / self.total * 100) if self.total > 0 else 0.0

    @property
    def status(self) -> str:
        if self.error:
            return "ERROR"
        elif self.passed == self.total:
            return "PASS"
        else:
            return "FAIL"


# ============================================================================
# MODULE DEFINITIONS
# ============================================================================

COMPARISON_MODULES = [
    # (module_name, script_name, description, is_quick)
    ("time", "compare_time.py", "Time functions (julday, deltat, sidtime)", True),
    ("planets", "compare_planets.py", "Planetary calculations", True),
    ("houses", "compare_houses.py", "House systems", True),
    ("sidereal", "compare_sidereal.py", "Sidereal/ayanamsha modes", True),
    ("lunar", "compare_lunar.py", "Lunar nodes and Lilith", True),
    ("crossings", "compare_crossings.py", "Sun/Moon crossings", True),
    ("crossings_ext", "compare_crossings_ext.py", "Extended crossings", True),
    ("observations", "compare_observations.py", "Observation modes", True),
    ("minor_bodies", "compare_minor_bodies.py", "Asteroids", True),
    ("eclipses", "compare_eclipses.py", "Solar and lunar eclipses", False),
    ("occultations", "compare_occultations.py", "Lunar occultations", False),
    ("rise_transit", "compare_rise_transit.py", "Rise/set/transit", False),
    ("fixedstars", "compare_fixedstars.py", "Fixed stars", True),
    ("heliacal", "compare_heliacal.py", "Heliacal events", False),
    ("coordinates", "compare_coordinates.py", "Coordinate transforms", True),
    ("orbital", "compare_orbital.py", "Orbital elements", True),
    ("phenomena", "compare_phenomena.py", "Planetary phenomena", True),
    ("houses_ext", "compare_houses_ext.py", "Extended house functions", True),
    ("utilities", "compare_utilities.py", "Utility functions", True),
]


# ============================================================================
# RUNNER FUNCTIONS
# ============================================================================


def run_module(script_name: str, verbose: bool = False) -> Tuple[int, int, str]:
    """
    Run a comparison module and extract results.

    Returns:
        (passed, total, output)
    """
    cmd = [sys.executable, "-u", f"compare_scripts/{script_name}"]
    if verbose:
        cmd.append("--verbose")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )
        output = result.stdout + result.stderr

        # Parse output for pass/total counts
        # Look for "Passed:" and "Total tests:" in summary section
        passed = 0
        total = 0

        for line in output.split("\n"):
            line = line.strip()
            # Match "Total tests:   39" format
            if line.startswith("Total tests:"):
                try:
                    total = int(line.split(":")[1].strip())
                except (ValueError, IndexError):
                    pass
            # Match "Passed:        24 ✓" format - avoid "Pass rate:"
            elif line.startswith("Passed:") and "✓" in line:
                try:
                    # Handle "Passed:        24 ✓" format
                    parts = line.split(":")[1].strip()
                    # Get the number before any symbols
                    num_str = parts.split()[0]
                    passed = int(num_str)
                except (ValueError, IndexError):
                    pass

        return passed, total, output

    except subprocess.TimeoutExpired:
        return 0, 0, "TIMEOUT: Module exceeded 5 minute limit"
    except Exception as e:
        return 0, 0, f"ERROR: {e}"


def run_all_modules(modules: List[Tuple], verbose: bool = False) -> List[ModuleResult]:
    """Run all specified modules and collect results."""
    results = []

    for name, script, desc, is_quick in modules:
        print(f"\n{'=' * 60}")
        print(f"Running: {name} ({desc})")
        print(f"{'=' * 60}")

        start_time = time.time()
        passed, total, output = run_module(script, verbose)
        duration = time.time() - start_time

        # Check for errors
        error = None
        if "ERROR" in output or "Traceback" in output:
            if passed == 0 and total == 0:
                error = "Module failed to run"

        result = ModuleResult(
            name=name,
            passed=passed,
            total=total,
            duration=duration,
            error=error,
        )
        results.append(result)

        # Print brief status
        status_char = (
            "✓"
            if result.status == "PASS"
            else ("!" if result.status == "ERROR" else "✗")
        )
        print(
            f"  Result: {result.passed}/{result.total} ({result.pass_rate:.1f}%) {status_char}"
        )
        print(f"  Time: {result.duration:.2f}s")

    return results


def print_summary_report(results: List[ModuleResult]):
    """Print comprehensive summary report."""
    # Calculate totals
    total_passed = sum(r.passed for r in results)
    total_tests = sum(r.total for r in results)
    total_time = sum(r.duration for r in results)

    modules_passed = sum(1 for r in results if r.status == "PASS")
    modules_failed = sum(1 for r in results if r.status == "FAIL")
    modules_error = sum(1 for r in results if r.status == "ERROR")

    print("\n")
    print("=" * 80)
    print("       LIBEPHEMERIS VS PYSWISSEPH - COMPREHENSIVE COMPARISON REPORT")
    print("=" * 80)

    # Module summary table
    print("\n" + "-" * 80)
    print(
        f"{'Module':<20} {'Tests':>8} {'Passed':>8} {'Rate':>8} {'Time':>8} {'Status':>8}"
    )
    print("-" * 80)

    for r in results:
        status_symbol = {"PASS": "✓", "FAIL": "✗", "ERROR": "!"}.get(r.status, "?")

        print(
            f"{r.name:<20} {r.total:>8} {r.passed:>8} "
            f"{r.pass_rate:>7.1f}% {r.duration:>7.1f}s {status_symbol:>8}"
        )

    print("-" * 80)
    print(
        f"{'TOTAL':<20} {total_tests:>8} {total_passed:>8} "
        f"{(total_passed / total_tests * 100 if total_tests else 0):>7.1f}% {total_time:>7.1f}s"
    )
    print("=" * 80)

    # Overall summary
    print("\n" + "-" * 40)
    print("OVERALL SUMMARY")
    print("-" * 40)
    print(f"  Modules tested:    {len(results)}")
    print(f"  Modules passed:    {modules_passed} ✓")
    print(f"  Modules failed:    {modules_failed} ✗")
    print(f"  Modules errored:   {modules_error} !")
    print(f"  Total tests:       {total_tests}")
    print(f"  Tests passed:      {total_passed}")
    print(
        f"  Overall pass rate: {(total_passed / total_tests * 100 if total_tests else 0):.1f}%"
    )
    print(f"  Total time:        {total_time:.1f}s")

    # Final verdict
    print("\n" + "=" * 80)
    if modules_error > 0:
        print("  STATUS: SOME MODULES HAD ERRORS - Review output for details")
    elif modules_failed > 0:
        print("  STATUS: SOME TESTS FAILED - Libraries have differences")
    else:
        print("  STATUS: ALL TESTS PASSED! - Excellent compatibility ✓")
    print("=" * 80)

    return total_passed == total_tests and modules_error == 0


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print(__doc__)
    print("\nAvailable modules:")
    for name, script, desc, is_quick in COMPARISON_MODULES:
        quick_mark = "[quick]" if is_quick else "[slow]"
        print(f"  {name:<15} {quick_mark:<8} {desc}")


def list_modules():
    """List available modules."""
    print("\nAvailable comparison modules:")
    print("-" * 60)
    for name, script, desc, is_quick in COMPARISON_MODULES:
        quick_mark = "[quick]" if is_quick else "[slow]"
        print(f"  {name:<15} {quick_mark:<8} {desc}")
    print("-" * 60)
    print(f"\nTotal: {len(COMPARISON_MODULES)} modules")


def main():
    """Main entry point."""
    args = sys.argv[1:]

    # Parse arguments
    verbose = "--verbose" in args or "-v" in args
    quick_only = "--quick" in args

    if "--help" in args or "-h" in args:
        print_help()
        sys.exit(0)

    if "--list" in args:
        list_modules()
        sys.exit(0)

    # Check for specific module
    module_filter = None
    if "--module" in args:
        try:
            idx = args.index("--module")
            module_filter = args[idx + 1]
        except (IndexError, ValueError):
            print("Error: --module requires a module name")
            sys.exit(1)

    # Select modules to run
    modules = COMPARISON_MODULES

    if module_filter:
        modules = [
            (n, s, d, q) for n, s, d, q in COMPARISON_MODULES if n == module_filter
        ]
        if not modules:
            print(f"Error: Unknown module '{module_filter}'")
            print("Use --list to see available modules")
            sys.exit(1)
    elif quick_only:
        modules = [(n, s, d, q) for n, s, d, q in COMPARISON_MODULES if q]

    # Print header
    print("=" * 80)
    print("       LIBEPHEMERIS VS PYSWISSEPH COMPARISON SUITE")
    print("=" * 80)
    print(f"Running {len(modules)} comparison module(s)")
    if quick_only:
        print("(Quick mode - skipping slow tests)")
    if verbose:
        print("(Verbose mode enabled)")

    # Run modules
    results = run_all_modules(modules, verbose)

    # Print summary
    all_passed = print_summary_report(results)

    # Exit with appropriate code
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
