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
import threading
from dataclasses import dataclass
from typing import List, Tuple, Optional, Any, Callable

# Try to import tqdm for progress bars
try:
    from tqdm import tqdm as tqdm_func

    TQDM_AVAILABLE = True
except ImportError:
    tqdm_func: Any = None
    TQDM_AVAILABLE = False


@dataclass
class ErrorInfo:
    """Detailed information about a module error."""

    category: str
    message: str
    suggestion: str
    last_lines: List[str]


# Error patterns for categorization: (pattern, category, suggestion)
ERROR_PATTERNS: List[Tuple[str, str, str]] = [
    (
        "ModuleNotFoundError",
        "import_error",
        "Install missing module with: pip install <module_name>",
    ),
    (
        "ImportError",
        "import_error",
        "Check module installation and PYTHONPATH configuration",
    ),
    (
        "SyntaxError",
        "syntax_error",
        "Fix syntax error in the indicated file and line",
    ),
    (
        "IndentationError",
        "syntax_error",
        "Fix indentation in the indicated file and line",
    ),
    (
        "AssertionError",
        "assertion_error",
        "Check test assertions - expected values may differ from actual",
    ),
    (
        "FileNotFoundError",
        "file_error",
        "Ensure required data files exist (check ephemeris files)",
    ),
    (
        "PermissionError",
        "file_error",
        "Check file permissions for required data files",
    ),
    (
        "TimeoutError",
        "timeout_error",
        "Module execution exceeded time limit - consider optimizing or increasing timeout",
    ),
    (
        "TIMEOUT:",
        "timeout_error",
        "Module execution exceeded time limit - consider running with fewer tests",
    ),
    (
        "MemoryError",
        "resource_error",
        "Not enough memory - try running with fewer concurrent tests",
    ),
    (
        "RecursionError",
        "resource_error",
        "Infinite recursion detected - check for circular dependencies",
    ),
    (
        "AttributeError",
        "attribute_error",
        "Object missing expected attribute - check API compatibility",
    ),
    (
        "TypeError",
        "type_error",
        "Type mismatch in function call - check argument types",
    ),
    (
        "ValueError",
        "value_error",
        "Invalid value passed to function - check input parameters",
    ),
    (
        "KeyError",
        "key_error",
        "Dictionary key not found - check data structure",
    ),
    (
        "ZeroDivisionError",
        "math_error",
        "Division by zero occurred - check calculations",
    ),
    (
        "ConnectionError",
        "network_error",
        "Network connection failed - check internet connectivity",
    ),
]


def analyze_error(output: str) -> Optional[ErrorInfo]:
    """
    Analyze module output to extract detailed error information.

    Args:
        output: The complete output from running a module

    Returns:
        ErrorInfo with category, message, suggestion, and last lines,
        or None if no error detected
    """
    if not output:
        return None

    # Check if there's an error in the output
    has_error = "ERROR" in output or "Traceback" in output or "TIMEOUT:" in output

    if not has_error:
        return None

    # Get last 10 non-empty lines for context
    lines = output.strip().split("\n")
    non_empty_lines = [line for line in lines if line.strip()]
    last_lines = non_empty_lines[-10:] if len(non_empty_lines) > 10 else non_empty_lines

    # Try to categorize the error
    category = "unknown_error"
    suggestion = "Review the error output and module code for issues"
    error_message = "Module failed to run"

    # Check each error pattern
    for pattern, cat, sugg in ERROR_PATTERNS:
        if pattern in output:
            category = cat
            suggestion = sugg

            # Try to extract the specific error message
            for line in reversed(lines):
                if pattern in line:
                    error_message = line.strip()
                    break
            break

    # Special handling for tracebacks - try to get the actual error line
    if "Traceback" in output and error_message == "Module failed to run":
        # Last non-empty line often contains the error
        for line in reversed(lines):
            stripped = line.strip()
            if (
                stripped
                and not stripped.startswith("File ")
                and "Traceback" not in stripped
            ):
                error_message = stripped
                break

    return ErrorInfo(
        category=category,
        message=error_message,
        suggestion=suggestion,
        last_lines=last_lines,
    )


def format_error_report(error_info: ErrorInfo) -> str:
    """
    Format error information for display.

    Args:
        error_info: The ErrorInfo object containing error details

    Returns:
        Formatted string with error details
    """
    lines = [
        f"Error Type: {error_info.category.replace('_', ' ').title()}",
        f"Message: {error_info.message}",
        f"Suggestion: {error_info.suggestion}",
        "",
        "Last 10 lines of output:",
        "-" * 40,
    ]
    lines.extend(f"  {line}" for line in error_info.last_lines)
    lines.append("-" * 40)

    return "\n".join(lines)


@dataclass
class ModuleResult:
    """Result from running a comparison module."""

    name: str
    passed: int
    total: int
    duration: float
    error: Optional[str] = None
    error_info: Optional[ErrorInfo] = None

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
    ("minor_bodies", "compare_minor_bodies.py", "Asteroids, Centaurs, TNOs", True),
    ("hypothetical", "compare_hypothetical.py", "Uranian/hypothetical planets", True),
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


def run_module(
    script_name: str, verbose: bool = False, progress_interval: int = 30
) -> Tuple[int, int, str]:
    """
    Run a comparison module and extract results with progress reporting.

    Args:
        script_name: Name of the script to run
        verbose: Whether to pass --verbose flag
        progress_interval: Seconds between progress reports (default 30)

    Returns:
        (passed, total, output)
    """
    cmd = [sys.executable, "-u", f"compare_scripts/{script_name}"]
    if verbose:
        cmd.append("--verbose")

    timeout = 300  # 5 minute timeout

    try:
        # Use Popen for progress monitoring
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )

        start_time = time.time()
        last_progress_time = start_time
        output_lines: List[str] = []
        stderr_lines: List[str] = []

        # Progress indicator state
        progress_chars = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
        progress_idx = 0

        # Thread to read stdout without blocking
        def read_stdout():
            if process.stdout:
                for line in iter(process.stdout.readline, ""):
                    output_lines.append(line)

        def read_stderr():
            if process.stderr:
                for line in iter(process.stderr.readline, ""):
                    stderr_lines.append(line)

        stdout_thread = threading.Thread(target=read_stdout, daemon=True)
        stderr_thread = threading.Thread(target=read_stderr, daemon=True)
        stdout_thread.start()
        stderr_thread.start()

        # Poll with progress reporting
        while process.poll() is None:
            elapsed = time.time() - start_time

            # Check timeout
            if elapsed >= timeout:
                process.kill()
                stdout_thread.join(timeout=1)
                stderr_thread.join(timeout=1)
                return 0, 0, "TIMEOUT: Module exceeded 5 minute limit"

            # Print progress every progress_interval seconds
            current_time = time.time()
            if current_time - last_progress_time >= progress_interval:
                last_progress_time = current_time
                elapsed_mins = int(elapsed // 60)
                elapsed_secs = int(elapsed % 60)
                remaining = timeout - elapsed
                remaining_mins = int(remaining // 60)
                remaining_secs = int(remaining % 60)

                # Show spinner and progress info
                spinner = progress_chars[progress_idx % len(progress_chars)]
                progress_idx += 1

                print(
                    f"    {spinner} Running {script_name}... "
                    f"[{elapsed_mins:02d}:{elapsed_secs:02d} elapsed, "
                    f"{remaining_mins:02d}:{remaining_secs:02d} remaining]",
                    flush=True,
                )

            # Small sleep to avoid busy-waiting
            time.sleep(0.1)

        # Wait for threads to finish
        stdout_thread.join(timeout=2)
        stderr_thread.join(timeout=2)

        output = "".join(output_lines) + "".join(stderr_lines)

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

    except Exception as e:
        return 0, 0, f"ERROR: {e}"


def run_all_modules(modules: List[Tuple], verbose: bool = False) -> List[ModuleResult]:
    """Run all specified modules and collect results."""
    results = []
    total_modules = len(modules)

    # Create iterator with tqdm if available
    if TQDM_AVAILABLE:
        module_iter = tqdm_func(
            enumerate(modules),
            total=total_modules,
            desc="Running modules",
            unit="module",
        )
    else:
        module_iter = enumerate(modules)

    for idx, (name, script, desc, is_quick) in module_iter:
        module_num = idx + 1
        print(f"\n{'=' * 60}")
        print(f"Running [{module_num}/{total_modules}]: {name} ({desc})")
        print(f"{'=' * 60}")

        start_time = time.time()
        passed, total, output = run_module(script, verbose)
        duration = time.time() - start_time

        # Check for errors with detailed analysis
        error = None
        error_info = None
        if "ERROR" in output or "Traceback" in output:
            if passed == 0 and total == 0:
                error_info = analyze_error(output)
                if error_info:
                    error = error_info.message
                else:
                    error = "Module failed to run"

        result = ModuleResult(
            name=name,
            passed=passed,
            total=total,
            duration=duration,
            error=error,
            error_info=error_info,
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

        # Print detailed error information if available
        if error_info:
            print(f"\n  {format_error_report(error_info)}")

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
