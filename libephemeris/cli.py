"""Command-line interface for libephemeris.

This module provides CLI commands for managing libephemeris data files.
Migrated from argparse to click for tab completion and consistency with
the dev CLI (``leph``).

Usage:
    libephemeris download base          Download data for 'base' tier (1850-2150)
    libephemeris download medium        Download data for 'medium' tier (1550-2650)
    libephemeris download extended      Download data for 'extended' tier (-13200 to +17191)
    libephemeris download leb-base      Download LEB binary ephemeris for 'base' (~53 MB)
    libephemeris download leb-medium    Download LEB binary ephemeris for 'medium' (~175 MB)
    libephemeris download leb-extended  Download LEB binary ephemeris for 'extended'
    libephemeris download assist        Download ASSIST n-body data files (~714 MB)
    libephemeris status                 Show data file status
    libephemeris info                   Show version, calc mode, LEB file, tier
    libephemeris --version              Show version
    libephemeris --help                 Show help

Shell completion (add to your shell profile):

    # zsh
    eval "$(_LIBEPHEMERIS_COMPLETE=zsh_source libephemeris)"

    # bash
    eval "$(_LIBEPHEMERIS_COMPLETE=bash_source libephemeris)"
"""

from __future__ import annotations

import sys

import click

from . import __version__
from .cli_shared import TIER_INFO, leb_download_help, tier_download_help


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _handle_download(func, quiet: bool, **kwargs) -> None:  # type: ignore[no-untyped-def]
    """Run a download function with standard error handling."""
    try:
        func(**kwargs)
    except KeyboardInterrupt:
        click.echo("\nDownload cancelled.")
        sys.exit(130)
    except ImportError as e:
        if not quiet:
            click.echo(
                f"Error: {e}\nInstall the nbody extra: pip install libephemeris[nbody]",
                err=True,
            )
        sys.exit(1)
    except (OSError, ValueError, RuntimeError) as e:
        if not quiet:
            click.echo(f"Error: {e}", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Root CLI group
# ---------------------------------------------------------------------------


@click.group(
    name="libephemeris",
    help="High-precision astronomical ephemeris library.\n\n"
    "This CLI manages data files, shows library status, and configures\n"
    "the calculation backend. It is intended for end-users and CI pipelines.\n\n"
    "First-time setup:\n\n"
    "  libephemeris download medium       Download data for the default tier\n"
    "  libephemeris status                Verify everything is installed\n"
    "  libephemeris info                  Show active backend and configuration",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 120},
    epilog="""\
Examples:
  libephemeris download medium        Download data for the default tier
  libephemeris download base          Lightweight, modern-era data
  libephemeris download extended      Full range (-13200 to +17191 CE)
  libephemeris download leb-base      LEB binary ephemeris (~53 MB, ~14x speedup)
  libephemeris download leb-medium    LEB binary ephemeris (~175 MB, ~14x speedup)
  libephemeris download assist        ASSIST n-body data (~714 MB)
  libephemeris status                 Show installed data files
  libephemeris info                   Show version, calc mode, tier info

For more information, visit: https://github.com/g-battaglia/libephemeris
""",
)
@click.version_option(__version__, prog_name="libephemeris")
def cli() -> None:
    """Root CLI group."""


# ---------------------------------------------------------------------------
# download subgroup
# ---------------------------------------------------------------------------


@click.group(
    "download",
    help="Download data files required by libephemeris.\n\n"
    "Three types of data:\n\n"
    "  Tier data (base/medium/extended)  DE440/441 kernels + asteroid SPKs\n"
    "  LEB files (leb-base/medium/ext)   Precomputed Chebyshev (~14x speedup)\n"
    "  ASSIST data                       N-body integration files (~714 MB)\n\n"
    "Most users need only: libephemeris download medium",
)
def download_group() -> None:
    """Download subcommands."""


# --- Common download options ---


def _download_options(f):  # type: ignore[no-untyped-def]
    """Add --force / --no-progress / --quiet to a download command."""
    f = click.option(
        "--force", "-f", is_flag=True, help="Force download even if files already exist"
    )(f)
    f = click.option("--no-progress", is_flag=True, help="Disable progress output")(f)
    f = click.option(
        "--quiet", "-q", is_flag=True, help="Suppress all output except errors"
    )(f)
    return f


# --- Tier data downloads ---


def _make_tier_download(tier: str) -> click.Command:
    """Create a download command for a tier."""

    @click.command(
        tier,
        help=f"Download data for '{tier}' tier ({TIER_INFO[tier]['range']}).\n\n"
        + tier_download_help(tier),
    )
    @_download_options
    def cmd(force: bool, no_progress: bool, quiet: bool) -> None:
        from .download import download_for_tier

        _handle_download(
            download_for_tier,
            quiet=quiet,
            tier_name=tier,
            force=force,
            show_progress=not no_progress,
        )

    return cmd


for _tier in TIER_INFO:
    download_group.add_command(_make_tier_download(_tier))


# --- LEB downloads ---


def _make_leb_download(tier: str) -> click.Command:
    """Create a LEB download command for a tier."""

    @click.command(
        f"leb-{tier}",
        help=f"Download LEB binary ephemeris for '{tier}' tier.\n\n"
        + leb_download_help(tier),
    )
    @_download_options
    def cmd(force: bool, no_progress: bool, quiet: bool) -> None:
        from .download import download_leb_for_tier

        _handle_download(
            download_leb_for_tier,
            quiet=quiet,
            tier_name=tier,
            force=force,
            show_progress=not no_progress,
            activate=False,
        )

    return cmd


for _tier in TIER_INFO:
    download_group.add_command(_make_leb_download(_tier))


# --- ASSIST download ---


@download_group.command()
@_download_options
@click.option("--no-planets", is_flag=True, help="Skip planet ephemeris download")
@click.option("--no-asteroids", is_flag=True, help="Skip asteroid perturbers download")
@click.option(
    "--target-dir",
    type=str,
    default=None,
    help="Directory to save files (default: ~/.libephemeris/assist/)",
)
def assist(
    force: bool,
    no_progress: bool,
    quiet: bool,
    no_planets: bool,
    no_asteroids: bool,
    target_dir: str | None,
) -> None:
    """Download ASSIST n-body data files (~714 MB).

    ASSIST provides sub-arcsecond precision for asteroid orbit propagation by
    including gravitational perturbations from the Sun, Moon, 8 planets, and
    16 massive asteroids.

    \b
    Downloads:
      1. Planet ephemeris (linux_p1550p2650.440, ~98 MB)
      2. Asteroid perturbers (sb441-n16.bsp, ~616 MB)

    Files are saved to ~/.libephemeris/assist/ by default.
    Requires: pip install libephemeris[nbody]
    """
    from .rebound_integration import download_assist_data

    _handle_download(
        download_assist_data,
        quiet=quiet,
        target_dir=target_dir,
        planets=not no_planets,
        asteroids=not no_asteroids,
        force=force,
        show_progress=not no_progress,
    )


cli.add_command(download_group)


# ---------------------------------------------------------------------------
# status command
# ---------------------------------------------------------------------------


@cli.command(short_help="Show installed data files, current tier, and download status.")
def status() -> None:
    """Show installed data files, current precision tier, and download status.

    Lists all data files (BSP kernels, SPKs, LEB files), their sizes,
    and whether they are present on disk.
    """
    from .download import print_data_status

    print_data_status()


# ---------------------------------------------------------------------------
# info command (NEW)
# ---------------------------------------------------------------------------


@cli.command(short_help="Show version, active calc mode, LEB file, and precision tier.")
def info() -> None:
    """Show library version, active calculation mode, LEB file path, and precision tier.

    Displays a quick summary of the current libephemeris configuration:
    which backend is active (skyfield/leb/horizons/auto), the precision tier,
    the configured LEB file (if any), and the data directory path.
    """
    from . import __version__ as ver

    click.echo(f"libephemeris {ver}")
    click.echo()

    # Calculation mode
    try:
        import os

        from .state import get_calc_mode, get_precision_tier

        mode = get_calc_mode()
        click.echo(f"  Calc mode:      {mode}")

        tier = get_precision_tier()
        tier_info = TIER_INFO.get(tier, {})
        click.echo(f"  Precision tier:  {tier} ({tier_info.get('range', 'unknown')})")

        # LEB file: check env var and internal state
        leb_path = os.environ.get("LIBEPHEMERIS_LEB")
        if not leb_path:
            from . import state as _state

            leb_path = getattr(_state, "_LEB_FILE", None)
        if leb_path:
            click.echo(f"  LEB file:        {leb_path}")
        else:
            click.echo("  LEB file:        (none configured)")
    except Exception as e:
        click.echo(f"  (could not read state: {e})")

    # Data directory
    try:
        from .download import get_data_dir

        data_dir = get_data_dir()
        click.echo(f"  Data directory:  {data_dir}")
    except Exception:
        pass

    click.echo()
    click.echo("For detailed file status: libephemeris status")


# ---------------------------------------------------------------------------
# Backward compatibility aliases (download:base -> download base)
# ---------------------------------------------------------------------------
# These allow the old colon-separated syntax to still work for users
# who have scripts or muscle memory using the old argparse-based CLI.

_LEGACY_ALIASES = {
    "download:base": ["download", "base"],
    "download:medium": ["download", "medium"],
    "download:extended": ["download", "extended"],
    "download:leb:base": ["download", "leb-base"],
    "download:leb:medium": ["download", "leb-medium"],
    "download:leb:extended": ["download", "leb-extended"],
    "download:assist": ["download", "assist"],
}


def main(argv: list[str] | None = None) -> None:
    """Main entry point for the CLI.

    Supports both old (colon-separated) and new (space-separated) command syntax.
    """
    args = argv if argv is not None else sys.argv[1:]

    # Rewrite legacy colon-separated commands to space-separated equivalents
    if args:
        first = args[0]
        if first in _LEGACY_ALIASES:
            args = [*_LEGACY_ALIASES[first], *args[1:]]

    cli(args=args, standalone_mode=True)


if __name__ == "__main__":
    main()
