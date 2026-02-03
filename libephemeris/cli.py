"""
Command-line interface for libephemeris.

This module provides CLI commands for managing libephemeris, including
downloading optional data files.

Usage:
    libephemeris download-data       Download precision data files
    libephemeris status              Show data file status
    libephemeris --version           Show version
    libephemeris --help              Show help
"""

from __future__ import annotations

import argparse
import sys

from . import __version__


def cmd_download_data(args: argparse.Namespace) -> int:
    """Handle the download-data command."""
    from .download import download_planet_centers

    try:
        download_planet_centers(
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


def cmd_status(args: argparse.Namespace) -> int:
    """Handle the status command."""
    from .download import print_data_status

    print_data_status()
    return 0


def cmd_version(args: argparse.Namespace) -> int:
    """Handle the version command."""
    print(f"libephemeris {__version__}")
    return 0


def create_parser() -> argparse.ArgumentParser:
    """Create the argument parser for the CLI."""
    parser = argparse.ArgumentParser(
        prog="libephemeris",
        description="High-precision astronomical ephemeris library",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  libephemeris download-data     Download optional precision data files
  libephemeris status            Show installed data files
  libephemeris --version         Show version information

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

    # download-data command
    download_parser = subparsers.add_parser(
        "download-data",
        help="Download optional precision data files",
        description="""
Download optional data files that enhance calculation precision.

The main data file is planet_centers.bsp which provides precise planet
center positions for Jupiter, Saturn, Uranus, Neptune, and Pluto.
This enables sub-arcsecond precision for outer planet calculations.

Without these files, libephemeris uses analytical approximations which
are still accurate to ~0.1 arcseconds.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    download_parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Force download even if file already exists",
    )
    download_parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bar",
    )
    download_parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress all output except errors",
    )
    download_parser.set_defaults(func=cmd_download_data)

    # status command
    status_parser = subparsers.add_parser(
        "status",
        help="Show data file status",
        description="Show the status of optional data files.",
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
