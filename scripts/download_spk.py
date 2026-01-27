#!/usr/bin/env python3
"""
Download SPK files for common celestial bodies.

This script pre-downloads SPK (SPICE kernel) files from JPL Horizons for common
minor bodies, allowing users to populate the cache ahead of time rather than
downloading on demand.

Usage:
    python scripts/download_spk.py                    # Download all common bodies
    python scripts/download_spk.py --bodies chiron ceres  # Download specific bodies
    python scripts/download_spk.py --list             # List available bodies
    python scripts/download_spk.py --cache-dir /path  # Use custom cache directory

Requirements:
    pip install astroquery

"""

import argparse
import os
import sys
from typing import Optional

# Ensure libephemeris can be imported from the project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# Body definitions: name -> (horizons_id, ipl, naif_id)
AVAILABLE_BODIES: dict[str, tuple[str, int, int]] = {}


def _init_bodies() -> None:
    """Initialize body definitions from libephemeris constants."""
    global AVAILABLE_BODIES
    try:
        from libephemeris.constants import (
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
            NAIF_CHIRON,
            NAIF_PHOLUS,
            NAIF_CERES,
            NAIF_PALLAS,
            NAIF_JUNO,
            NAIF_VESTA,
            NAIF_ERIS,
            NAIF_SEDNA,
            NAIF_HAUMEA,
            NAIF_MAKEMAKE,
            NAIF_IXION,
            NAIF_ORCUS,
            NAIF_QUAOAR,
        )

        AVAILABLE_BODIES = {
            "chiron": ("2060", SE_CHIRON, NAIF_CHIRON),
            "pholus": ("5145", SE_PHOLUS, NAIF_PHOLUS),
            "ceres": ("1", SE_CERES, NAIF_CERES),
            "pallas": ("2", SE_PALLAS, NAIF_PALLAS),
            "juno": ("3", SE_JUNO, NAIF_JUNO),
            "vesta": ("4", SE_VESTA, NAIF_VESTA),
            "eris": ("136199", SE_ERIS, NAIF_ERIS),
            "sedna": ("90377", SE_SEDNA, NAIF_SEDNA),
            "haumea": ("136108", SE_HAUMEA, NAIF_HAUMEA),
            "makemake": ("136472", SE_MAKEMAKE, NAIF_MAKEMAKE),
            "ixion": ("28978", SE_IXION, NAIF_IXION),
            "orcus": ("90482", SE_ORCUS, NAIF_ORCUS),
            "quaoar": ("50000", SE_QUAOAR, NAIF_QUAOAR),
            "varuna": ("20000", SE_VARUNA, 2020000),
        }
    except ImportError as e:
        print(f"Error importing libephemeris constants: {e}", file=sys.stderr)
        sys.exit(1)


# Common body groups for convenience
COMMON_BODIES = [
    "chiron",
    "pholus",
    "ceres",
    "pallas",
    "juno",
    "vesta",
    "eris",
    "sedna",
]
ALL_BODIES = list(AVAILABLE_BODIES.keys()) if AVAILABLE_BODIES else []


def check_astroquery() -> bool:
    """Check if astroquery is available."""
    try:
        from astroquery.jplhorizons import Horizons

        return True
    except ImportError:
        return False


def download_spk_for_body(
    body_name: str,
    start_date: str,
    end_date: str,
    cache_dir: Optional[str] = None,
    force: bool = False,
    verbose: bool = True,
) -> Optional[str]:
    """
    Download SPK file for a specific body.

    Args:
        body_name: Name of the body (e.g., 'chiron', 'ceres')
        start_date: Start date in YYYY-MM-DD format
        end_date: End date in YYYY-MM-DD format
        cache_dir: Optional custom cache directory
        force: If True, re-download even if file exists
        verbose: If True, print progress messages

    Returns:
        Path to downloaded SPK file, or None on failure
    """
    from libephemeris import spk_auto

    body_name_lower = body_name.lower()
    if body_name_lower not in AVAILABLE_BODIES:
        if verbose:
            print(f"Unknown body: {body_name}", file=sys.stderr)
        return None

    horizons_id, ipl, naif_id = AVAILABLE_BODIES[body_name_lower]

    # Determine cache directory
    if cache_dir is None:
        cache_dir = spk_auto.DEFAULT_AUTO_SPK_DIR

    # Ensure cache directory exists
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir, exist_ok=True)

    # Convert dates to Julian Days for the filename
    from libephemeris.spk_auto import _jd_to_iso_date

    def _iso_to_jd(date_str: str) -> float:
        """Convert ISO date string to Julian Day."""
        parts = date_str.split("-")
        year, month, day = int(parts[0]), int(parts[1]), int(parts[2])

        # Algorithm from Meeus, Astronomical Algorithms
        if month <= 2:
            year -= 1
            month += 12

        a = int(year / 100)
        b = 2 - a + int(a / 4)

        jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + b - 1524.5
        return jd

    jd_start = _iso_to_jd(start_date)
    jd_end = _iso_to_jd(end_date)

    # Generate filename
    filename = spk_auto._generate_spk_cache_filename(horizons_id, jd_start, jd_end)
    output_path = os.path.join(cache_dir, filename)

    # Check if already cached
    if os.path.exists(output_path) and not force:
        if verbose:
            print(f"  {body_name}: Already cached at {output_path}")
        return output_path

    # Download
    if verbose:
        print(f"  {body_name}: Downloading from JPL Horizons...")

    try:
        spk_path = spk_auto.download_spk_from_horizons(
            body_id=horizons_id,
            jd_start=jd_start,
            jd_end=jd_end,
            output_path=output_path,
        )
        if verbose:
            size_mb = os.path.getsize(spk_path) / (1024 * 1024)
            print(f"  {body_name}: Downloaded {size_mb:.2f} MB to {spk_path}")
        return spk_path
    except Exception as e:
        if verbose:
            print(f"  {body_name}: Failed - {e}", file=sys.stderr)
        return None


def list_available_bodies() -> None:
    """Print list of available bodies."""
    print("\nAvailable bodies for SPK download:")
    print("=" * 60)

    # Group by category
    centaurs = ["chiron", "pholus"]
    main_belt = ["ceres", "pallas", "juno", "vesta"]
    tnos = ["eris", "sedna", "haumea", "makemake", "ixion", "orcus", "quaoar", "varuna"]

    print("\nCentaurs:")
    for name in centaurs:
        if name in AVAILABLE_BODIES:
            horizons_id, ipl, naif_id = AVAILABLE_BODIES[name]
            print(f"  {name:12s} - Horizons ID: {horizons_id:8s} (NAIF: {naif_id})")

    print("\nMain Belt Asteroids:")
    for name in main_belt:
        if name in AVAILABLE_BODIES:
            horizons_id, ipl, naif_id = AVAILABLE_BODIES[name]
            print(f"  {name:12s} - Horizons ID: {horizons_id:8s} (NAIF: {naif_id})")

    print("\nTrans-Neptunian Objects (TNOs):")
    for name in tnos:
        if name in AVAILABLE_BODIES:
            horizons_id, ipl, naif_id = AVAILABLE_BODIES[name]
            print(f"  {name:12s} - Horizons ID: {horizons_id:8s} (NAIF: {naif_id})")

    print("\nPreset groups:")
    print("  --common     Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris, Sedna")
    print("  --all        All available bodies")


def list_cache_contents(cache_dir: Optional[str] = None) -> None:
    """List contents of the SPK cache."""
    from libephemeris import spk_auto

    cached = spk_auto.list_cached_spk(cache_dir)

    if not cached:
        print("\nSPK cache is empty.")
        return

    print(f"\nCached SPK files ({len(cached)} files):")
    print("=" * 80)

    total_size = 0.0
    for entry in cached:
        filename = entry["filename"]
        size_mb = float(entry["size_mb"])
        date_start = entry.get("date_start", "Unknown")
        date_end = entry.get("date_end", "Unknown")
        total_size += size_mb

        print(f"  {filename}")
        print(f"    Size: {size_mb:.2f} MB | Coverage: {date_start} to {date_end}")

    print(f"\nTotal cache size: {total_size:.2f} MB")


def main() -> int:
    """Main entry point."""
    # Initialize body definitions
    _init_bodies()

    # Update ALL_BODIES after initialization
    global ALL_BODIES
    ALL_BODIES = list(AVAILABLE_BODIES.keys())

    parser = argparse.ArgumentParser(
        description="Download SPK files for common celestial bodies.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/download_spk.py                        # Download common bodies
  python scripts/download_spk.py --all                  # Download all bodies
  python scripts/download_spk.py --bodies chiron ceres  # Download specific bodies
  python scripts/download_spk.py --list                 # List available bodies
  python scripts/download_spk.py --list-cache           # Show cache contents
  python scripts/download_spk.py --start 2020-01-01 --end 2050-01-01
        """,
    )

    parser.add_argument(
        "--bodies",
        nargs="+",
        help="Specific bodies to download (e.g., chiron ceres eris)",
    )
    parser.add_argument(
        "--common",
        action="store_true",
        help="Download common bodies (Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris, Sedna)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        dest="download_all",
        help="Download all available bodies",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_bodies",
        help="List available bodies and exit",
    )
    parser.add_argument(
        "--list-cache",
        action="store_true",
        help="List cached SPK files and exit",
    )
    parser.add_argument(
        "--start",
        default="2000-01-01",
        help="Start date for SPK coverage (default: 2000-01-01)",
    )
    parser.add_argument(
        "--end",
        default="2100-01-01",
        help="End date for SPK coverage (default: 2100-01-01)",
    )
    parser.add_argument(
        "--cache-dir",
        help="Custom cache directory (default: ~/.libephemeris/spk/)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if file already exists",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress output",
    )

    args = parser.parse_args()

    # Handle list options
    if args.list_bodies:
        list_available_bodies()
        return 0

    if args.list_cache:
        list_cache_contents(args.cache_dir)
        return 0

    # Check astroquery
    if not check_astroquery():
        print(
            "Error: astroquery is required for downloading SPK files.",
            file=sys.stderr,
        )
        print("Install it with: pip install astroquery", file=sys.stderr)
        return 1

    # Determine which bodies to download
    bodies_to_download: list[str] = []

    if args.bodies:
        bodies_to_download = [b.lower() for b in args.bodies]
        # Validate bodies
        for body in bodies_to_download:
            if body not in AVAILABLE_BODIES:
                print(f"Error: Unknown body '{body}'", file=sys.stderr)
                print("Use --list to see available bodies", file=sys.stderr)
                return 1
    elif args.download_all:
        bodies_to_download = ALL_BODIES
    else:
        # Default to common bodies
        bodies_to_download = COMMON_BODIES

    if not bodies_to_download:
        print("No bodies to download.", file=sys.stderr)
        return 1

    verbose = not args.quiet

    if verbose:
        print(f"\nDownloading SPK files for {len(bodies_to_download)} bodies...")
        print(f"Date range: {args.start} to {args.end}")
        if args.cache_dir:
            print(f"Cache directory: {args.cache_dir}")
        print()

    # Download each body
    success_count = 0
    fail_count = 0

    for body_name in bodies_to_download:
        result = download_spk_for_body(
            body_name=body_name,
            start_date=args.start,
            end_date=args.end,
            cache_dir=args.cache_dir,
            force=args.force,
            verbose=verbose,
        )
        if result:
            success_count += 1
        else:
            fail_count += 1

    if verbose:
        print()
        print(f"Completed: {success_count} successful, {fail_count} failed")

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
