#!/usr/bin/env python3
"""
Shared diagnostic logic for precision tier scripts.

This module provides the core functionality used by the per-tier diagnostic
scripts (tier_diagnostic_base.py, tier_diagnostic_medium.py,
tier_diagnostic_extended.py).

It calculates positions for all supported celestial bodies and prints a
detailed table showing ecliptic coordinates, equatorial coordinates,
velocities, and whether each body uses DE4xx, SPK, Keplerian fallback,
or analytical computation.
"""

from __future__ import annotations

import argparse
import os
import sys

# Ensure libephemeris can be imported from the project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import libephemeris as eph  # noqa: E402
from libephemeris import state  # noqa: E402
from libephemeris.constants import (  # noqa: E402
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_CHIRON,
    SE_PHOLUS,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
    SE_ERIS,
    SE_SEDNA,
    SE_HAUMEA,
    SE_MAKEMAKE,
    SE_IXION,
    SE_ORCUS,
    SE_QUAOAR,
    SE_VARUNA,
    SE_NESSUS,
    SE_ASBOLUS,
    SE_CHARIKLO,
    SE_GONGGONG,
    SE_APOPHIS,
    SE_HYGIEA,
    SE_EROS,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
)
from libephemeris.spk_auto import discover_local_spks  # noqa: E402
from libephemeris.state import TIERS, set_precision_tier  # noqa: E402


# =============================================================================
# BODY DEFINITIONS
# =============================================================================

# (ipl, display_name, category)
# category: "major", "lunar", "minor"
BODIES: list[tuple[int, str, str]] = [
    # Major planets
    (SE_SUN, "Sun", "major"),
    (SE_MOON, "Moon", "major"),
    (SE_MERCURY, "Mercury", "major"),
    (SE_VENUS, "Venus", "major"),
    (SE_MARS, "Mars", "major"),
    (SE_JUPITER, "Jupiter", "major"),
    (SE_SATURN, "Saturn", "major"),
    (SE_URANUS, "Uranus", "major"),
    (SE_NEPTUNE, "Neptune", "major"),
    (SE_PLUTO, "Pluto", "major"),
    # Lunar points
    (SE_MEAN_NODE, "Mean Node", "lunar"),
    (SE_TRUE_NODE, "True Node", "lunar"),
    # Minor bodies (from SPK_BODY_NAME_MAP)
    (SE_CHIRON, "Chiron", "minor"),
    (SE_PHOLUS, "Pholus", "minor"),
    (SE_CERES, "Ceres", "minor"),
    (SE_PALLAS, "Pallas", "minor"),
    (SE_JUNO, "Juno", "minor"),
    (SE_VESTA, "Vesta", "minor"),
    (SE_ERIS, "Eris", "minor"),
    (SE_SEDNA, "Sedna", "minor"),
    (SE_HAUMEA, "Haumea", "minor"),
    (SE_MAKEMAKE, "Makemake", "minor"),
    (SE_IXION, "Ixion", "minor"),
    (SE_ORCUS, "Orcus", "minor"),
    (SE_QUAOAR, "Quaoar", "minor"),
    (SE_VARUNA, "Varuna", "minor"),
    (SE_NESSUS, "Nessus", "minor"),
    (SE_ASBOLUS, "Asbolus", "minor"),
    (SE_CHARIKLO, "Chariklo", "minor"),
    (SE_GONGGONG, "Gonggong", "minor"),
    (SE_APOPHIS, "Apophis", "minor"),
    (SE_HYGIEA, "Hygiea", "minor"),
    (SE_EROS, "Eros", "minor"),
]

MAJOR_IPLS = {
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
}
LUNAR_IPLS = {SE_MEAN_NODE, SE_TRUE_NODE}

# Default dates per tier (label, year, month, day, hour)
DEFAULT_DATES: dict[str, list[tuple[str, int, int, int, float]]] = {
    "base": [
        ("J2000", 2000, 1, 1, 12.0),
        ("1950-06-15", 1950, 6, 15, 0.0),
        ("2024-11-05", 2024, 11, 5, 18.0),
    ],
    "medium": [
        ("J2000", 2000, 1, 1, 12.0),
        ("1900-01-01", 1900, 1, 1, 0.0),
        ("2024-11-05", 2024, 11, 5, 18.0),
        ("1600-01-01", 1600, 1, 1, 0.0),
    ],
    "extended": [
        ("J2000", 2000, 1, 1, 12.0),
        ("1900-01-01", 1900, 1, 1, 0.0),
        ("2024-11-05", 2024, 11, 5, 18.0),
        ("1550-01-01", 1550, 1, 1, 0.0),
        ("0800-01-01", 800, 1, 1, 0.0),
    ],
}


# =============================================================================
# FORMATTING HELPERS
# =============================================================================


def _format_deg(val: float) -> str:
    """Format degrees with 4 decimal places."""
    return f"{val:11.4f}\u00b0"


def _format_speed(val: float) -> str:
    """Format velocity in degrees/day."""
    return f"{val:+10.6f}"


def _format_dist(val: float) -> str:
    """Format distance in AU."""
    return f"{val:12.6f}"


def _format_dist_speed(val: float) -> str:
    """Format distance speed in AU/day."""
    return f"{val:+12.8f}"


def _format_hms(deg: float) -> str:
    """Convert degrees to hours:minutes:seconds string."""
    h_total = deg / 15.0
    h = int(h_total)
    m_total = (h_total - h) * 60.0
    m = int(m_total)
    s = (m_total - m) * 60.0
    return f"{h:2d}h{m:02d}m{s:05.2f}s"


def _format_dms(deg: float) -> str:
    """Convert degrees to signed degrees:arcmin:arcsec string."""
    sign = "+" if deg >= 0 else "-"
    deg = abs(deg)
    d = int(deg)
    m_total = (deg - d) * 60.0
    m = int(m_total)
    s = (m_total - m) * 60.0
    return f"{sign}{d:02d}\u00b0{m:02d}'{s:05.2f}\""


# =============================================================================
# SOURCE DETECTION
# =============================================================================


def _get_source(ipl: int, tier_name: str) -> str:
    """Determine the data source used for a body calculation.

    Returns one of: "DE440s", "DE440", "DE441", "Analytical", "SPK", "Keplerian".
    """
    if ipl in MAJOR_IPLS:
        eph_file = TIERS[tier_name].ephemeris_file
        return eph_file.replace(".bsp", "").upper()
    if ipl in LUNAR_IPLS:
        return "Analytical"
    if ipl in state._SPK_BODY_MAP:
        return "SPK"
    return "Keplerian"


# =============================================================================
# CALCULATION
# =============================================================================


def _calc_body(jd: float, ipl: int) -> dict:
    """Calculate all available data for a body at a given Julian Day.

    Returns dict with keys:
        ecl: (lon, lat, dist, lon_spd, lat_spd, dist_spd) or None
        equ: (ra, dec, dist, ra_spd, dec_spd, dist_spd) or None
        error: str or None
    """
    result: dict = {"ecl": None, "equ": None, "error": None}

    # Ecliptic coordinates + speed
    try:
        pos_ecl, _ = eph.swe_calc_ut(jd, ipl, SEFLG_SPEED)
        result["ecl"] = pos_ecl
    except Exception as e:
        result["error"] = str(e)[:60]
        return result

    # Equatorial coordinates + speed
    try:
        pos_equ, _ = eph.swe_calc_ut(jd, ipl, SEFLG_SPEED | SEFLG_EQUATORIAL)
        result["equ"] = pos_equ
    except Exception:
        # Equatorial may fail for some bodies; ecliptic is still valid
        pass

    return result


# =============================================================================
# TABLE OUTPUT
# =============================================================================


def _print_separator(char: str = "\u2550", width: int = 160) -> None:
    print(char * width)


def _print_header(tier_name: str) -> None:
    """Print the tier header box."""
    tier = TIERS[tier_name]
    start, end = tier.spk_date_range
    print()
    _print_separator("\u2550")
    print(
        f"  TIER: {tier_name} \u2014 {tier.ephemeris_file}"
        f" \u2014 SPK range: {start} to {end}"
    )
    print(f"  {tier.description}")
    _print_separator("\u2550")


def _print_column_headers() -> None:
    """Print the column header row."""
    print(
        f"  {'Body':<14s}"
        f" {'Lon':>11s}"
        f" {'Lat':>11s}"
        f" {'Dist(AU)':>12s}"
        f" {'LonSpd':>10s}"
        f" {'LatSpd':>10s}"
        f" {'DistSpd':>13s}"
        f"  {'RA':>12s}"
        f" {'Dec':>12s}"
        f" {'RASpd':>10s}"
        f" {'DecSpd':>10s}"
        f"  {'Source':<10s}"
    )
    _print_separator("\u2500")


def _print_group_label(label: str) -> None:
    """Print a group separator label."""
    print(f"  \u2500\u2500 {label} " + "\u2500" * (160 - len(label) - 6))


def _print_body_row(name: str, data: dict, source: str) -> None:
    """Print one body's data row."""
    if data["error"]:
        print(f"  {name:<14s}  ** {data['error']:<140s}  {source}")
        return

    ecl = data["ecl"]
    equ = data["equ"]

    lon = _format_deg(ecl[0])
    lat = _format_deg(ecl[1])
    dist = _format_dist(ecl[2])
    lon_spd = _format_speed(ecl[3])
    lat_spd = _format_speed(ecl[4])
    dist_spd = _format_dist_speed(ecl[5])

    if equ:
        ra = _format_hms(equ[0])
        dec = _format_dms(equ[1])
        ra_spd = _format_speed(equ[3])
        dec_spd = _format_speed(equ[4])
    else:
        ra = "     N/A    "
        dec = "     N/A    "
        ra_spd = "       N/A"
        dec_spd = "       N/A"

    print(
        f"  {name:<14s}"
        f" {lon}"
        f" {lat}"
        f" {dist}"
        f" {lon_spd}"
        f" {lat_spd}"
        f" {dist_spd}"
        f"  {ra}"
        f" {dec}"
        f" {ra_spd}"
        f" {dec_spd}"
        f"  {source}"
    )


def _print_summary(source_counts: dict[str, int], error_count: int) -> None:
    """Print the per-date summary line."""
    parts = []
    for src, count in sorted(source_counts.items()):
        parts.append(f"{count} {src}")
    if error_count > 0:
        parts.append(f"{error_count} Errors")
    print(f"\n  Summary: {' | '.join(parts)}")


# =============================================================================
# DATE PARSING
# =============================================================================


def _parse_date_str(date_str: str) -> tuple[str, int, int, int, float]:
    """Parse a YYYY-MM-DD string into (label, year, month, day, hour)."""
    parts = date_str.split("-")
    if len(parts) != 3:
        raise ValueError(f"Invalid date format: {date_str!r} (expected YYYY-MM-DD)")
    year = int(parts[0])
    month = int(parts[1])
    day = int(parts[2])
    return (date_str, year, month, day, 0.0)


def _parse_cli_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Tier diagnostic: calculate all bodies and show data sources.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python scripts/tier_diagnostic_base.py\n"
            "  python scripts/tier_diagnostic_base.py 1985-06-15\n"
            "  python scripts/tier_diagnostic_base.py 2000-01-01 1950-03-15\n"
            "  python scripts/tier_diagnostic_base.py --jd 2451545.0\n"
        ),
    )
    parser.add_argument(
        "dates",
        nargs="*",
        help="Dates in YYYY-MM-DD format (default: tier-specific set)",
    )
    parser.add_argument(
        "--jd",
        type=float,
        nargs="+",
        metavar="JD",
        help="Dates as Julian Day numbers",
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Omit equatorial coordinates (narrower table)",
    )
    return parser.parse_args()


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def run_diagnostic(tier_name: str) -> None:
    """Run the full tier diagnostic.

    Args:
        tier_name: One of "base", "medium", "extended".
    """
    args = _parse_cli_args()

    # Build date list
    dates: list[tuple[str, float]] = []  # (label, jd)

    if args.dates or args.jd:
        # User-provided dates
        if args.dates:
            for d in args.dates:
                info = _parse_date_str(d)
                jd = eph.swe_julday(info[1], info[2], info[3], info[4])
                dates.append((info[0], jd))
        if args.jd:
            for jd_val in args.jd:
                dates.append((f"JD {jd_val}", jd_val))
    else:
        # Default dates for this tier
        for label, y, m, d, h in DEFAULT_DATES[tier_name]:
            jd = eph.swe_julday(y, m, d, h)
            dates.append((label, jd))

    # Setup
    eph.close()
    set_precision_tier(tier_name)
    eph.set_strict_precision(False)
    eph.set_auto_spk_download(False)

    # Discover local SPK files
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    discovered = discover_local_spks(repo_root)

    # Print header
    _print_header(tier_name)

    if discovered:
        spk_bodies = [name for name, st in discovered.items() if st == "registered"]
        if spk_bodies:
            print(f"\n  Local SPK files discovered: {', '.join(spk_bodies)}")

    # Process each date
    for label, jd in dates:
        print(f"\n  Date: {label} (JD {jd:.1f})")
        _print_separator("\u2500")
        if not args.compact:
            _print_column_headers()

        source_counts: dict[str, int] = {}
        error_count = 0
        current_category = ""

        for ipl, name, category in BODIES:
            # Print group separator
            if category != current_category:
                current_category = category
                group_labels = {
                    "major": "Major planets",
                    "lunar": "Lunar points",
                    "minor": "Minor bodies",
                }
                _print_group_label(group_labels[category])

            # Calculate
            data = _calc_body(jd, ipl)
            source = _get_source(ipl, tier_name)

            if data["error"]:
                error_count += 1
            else:
                source_counts[source] = source_counts.get(source, 0) + 1

            # Print row
            if args.compact:
                _print_body_row_compact(name, data, source)
            else:
                _print_body_row(name, data, source)

        _print_summary(source_counts, error_count)

    # Cleanup
    print()
    _print_separator("\u2550")
    eph.close()


def _print_body_row_compact(name: str, data: dict, source: str) -> None:
    """Print a compact body row (ecliptic only, no equatorial)."""
    if data["error"]:
        print(f"  {name:<14s}  ** {data['error']:<80s}  {source}")
        return

    ecl = data["ecl"]
    lon = _format_deg(ecl[0])
    lat = _format_deg(ecl[1])
    dist = _format_dist(ecl[2])
    lon_spd = _format_speed(ecl[3])
    lat_spd = _format_speed(ecl[4])
    dist_spd = _format_dist_speed(ecl[5])

    print(f"  {name:<14s} {lon} {lat} {dist} {lon_spd} {lat_spd} {dist_spd}  {source}")
