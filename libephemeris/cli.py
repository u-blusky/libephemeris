"""
Command-line interface for libephemeris.

This module provides CLI commands for managing libephemeris data files.

Usage:
    libephemeris download:base         Download data for 'base' tier (1850-2150)
    libephemeris download:medium       Download data for 'medium' tier (1550-2650)
    libephemeris download:extended     Download data for 'extended' tier (-13200 to +17191)
    libephemeris status                Show data file status
    libephemeris --version             Show version
    libephemeris --help                Show help
"""

from __future__ import annotations

import argparse
import sys

from . import __version__


# ---------------------------------------------------------------------------
# Tier descriptions for --help output
# ---------------------------------------------------------------------------

_TIER_INFO = {
    "base": {
        "label": "base",
        "ephemeris": "de440s.bsp (~31 MB)",
        "range": "1850-2150 CE",
        "spk_range": "1850-2150",
        "description": "Lightweight tier for modern-era calculations.",
    },
    "medium": {
        "label": "medium",
        "ephemeris": "de440.bsp (~114 MB)",
        "range": "1550-2650 CE",
        "spk_range": "1900-2100",
        "description": "General purpose tier (default). Covers most historical and future dates.",
    },
    "extended": {
        "label": "extended",
        "ephemeris": "de441.bsp (~3.1 GB)",
        "range": "-13200 to +17191 CE",
        "spk_range": "1600-2500 (JPL Horizons limit)",
        "description": "Full extended range for deep historical and far-future research.",
    },
}


# ---------------------------------------------------------------------------
# Command handlers
# ---------------------------------------------------------------------------


def _cmd_download(tier_name: str, args: argparse.Namespace) -> int:
    """Handle a download:<tier> command."""
    from .download import download_for_tier

    try:
        download_for_tier(
            tier_name=tier_name,
            force=args.force,
            show_progress=not args.no_progress,
            quiet=args.quiet,
        )
        return 0
    except KeyboardInterrupt:
        print("\nDownload cancelled.")
        return 130
    except Exception as e:
        if not args.quiet:
            print(f"Error: {e}", file=sys.stderr)
        return 1


def cmd_download_base(args: argparse.Namespace) -> int:
    """Download data for the 'base' tier."""
    return _cmd_download("base", args)


def cmd_download_medium(args: argparse.Namespace) -> int:
    """Download data for the 'medium' tier."""
    return _cmd_download("medium", args)


def cmd_download_extended(args: argparse.Namespace) -> int:
    """Download data for the 'extended' tier."""
    return _cmd_download("extended", args)


def cmd_status(args: argparse.Namespace) -> int:
    """Handle the status command."""
    from .download import print_data_status

    print_data_status()
    return 0


def cmd_version(args: argparse.Namespace) -> int:
    """Handle the version command."""
    print(f"libephemeris {__version__}")
    return 0


# ---------------------------------------------------------------------------
# Parser construction
# ---------------------------------------------------------------------------


def _add_download_flags(parser: argparse.ArgumentParser) -> None:
    """Add common --force / --no-progress / --quiet flags to a download subparser."""
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Force download even if files already exist",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress output",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress all output except errors",
    )


def _make_download_description(tier: str) -> str:
    """Build the long description for a download:<tier> subparser."""
    info = _TIER_INFO[tier]
    return f"""\
Download all data files for the '{tier}' precision tier.

  Ephemeris:  {info["ephemeris"]}
  Date range: {info["range"]}
  SPK range:  {info["spk_range"]}

{info["description"]}

Downloads:
  1. The ephemeris file ({info["ephemeris"].split(" ")[0]})
  2. planet_centers.bsp precision offsets (~25 MB)
  3. SPK kernels for 21 minor bodies (asteroids, centaurs, TNOs)
"""


def create_parser() -> argparse.ArgumentParser:
    """Create the argument parser for the CLI."""
    parser = argparse.ArgumentParser(
        prog="libephemeris",
        description="High-precision astronomical ephemeris library",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  libephemeris download:medium       Download data for the default tier
  libephemeris download:base         Lightweight, modern-era data
  libephemeris download:extended     Full range (-13200 to +17191 CE)
  libephemeris status                Show installed data files
  libephemeris --version             Show version information

For more information, visit: https://github.com/g-battaglia/libephemeris
""",
    )

    parser.add_argument(
        "--version",
        action="store_true",
        help="Show version and exit",
    )

    # Subcommands
    subparsers = parser.add_subparsers(
        title="commands",
        dest="command",
        metavar="<command>",
    )

    # download:base
    dl_base = subparsers.add_parser(
        "download:base",
        help=f"Download data for 'base' tier ({_TIER_INFO['base']['range']})",
        description=_make_download_description("base"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_download_flags(dl_base)
    dl_base.set_defaults(func=cmd_download_base)

    # download:medium
    dl_medium = subparsers.add_parser(
        "download:medium",
        help=f"Download data for 'medium' tier ({_TIER_INFO['medium']['range']})",
        description=_make_download_description("medium"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_download_flags(dl_medium)
    dl_medium.set_defaults(func=cmd_download_medium)

    # download:extended
    dl_extended = subparsers.add_parser(
        "download:extended",
        help=f"Download data for 'extended' tier ({_TIER_INFO['extended']['range']})",
        description=_make_download_description("extended"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_download_flags(dl_extended)
    dl_extended.set_defaults(func=cmd_download_extended)

    # status
    status_parser = subparsers.add_parser(
        "status",
        help="Show data file status",
        description="Show the status of installed data files and current tier.",
    )
    status_parser.set_defaults(func=cmd_status)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Main entry point for the CLI.

    Args:
        argv: Command line arguments (defaults to sys.argv[1:])

    Returns:
        Exit code (0 for success, non-zero for errors)
    """
    parser = create_parser()
    args = parser.parse_args(argv)

    # Handle --version at top level
    if args.version:
        return cmd_version(args)

    # Handle no command
    if not args.command:
        parser.print_help()
        return 0

    # Run the command
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
