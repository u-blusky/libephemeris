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
    libephemeris download leb2-base     Download LEB2 compressed ephemeris for 'base'
    libephemeris download leb2-medium   Download LEB2 compressed ephemeris for 'medium'
    libephemeris download leb2-extended Download LEB2 compressed ephemeris for 'extended'
    libephemeris download assist        Download ASSIST n-body data files (~714 MB)
    libephemeris status                 Show comprehensive library and data status
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
    "  libephemeris status                Verify everything is installed",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 120},
    epilog="""\
Examples:
  libephemeris download medium        Download data for the default tier
  libephemeris download base          Lightweight, modern-era data
  libephemeris download extended      Full range (-13200 to +17191 CE)
  libephemeris download leb-base      LEB binary ephemeris (~53 MB, ~14x speedup)
  libephemeris download leb-medium    LEB binary ephemeris (~175 MB, ~14x speedup)
  libephemeris download leb2-base     LEB2 compressed ephemeris (smaller, modular)
  libephemeris download assist        ASSIST n-body data (~714 MB)
  libephemeris status                 Show comprehensive library and data status

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
    short_help="Download data files: tier SPKs, LEB, LEB2, ASSIST.",
    help="Download data files required by libephemeris.\n\n"
    "Four types of data:\n\n"
    "  Tier data (base/medium/extended)   DE440/441 kernels + asteroid SPKs\n"
    "  LEB files (leb-base/medium/ext)    Precomputed Chebyshev (~14x speedup)\n"
    "  LEB2 files (leb2-base/medium/ext)  Compressed modular (4-10x smaller)\n"
    "  ASSIST data                        N-body integration files (~714 MB)\n\n"
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
        short_help=f"Download DE kernel + SPKs for '{tier}' tier ({TIER_INFO[tier]['range']}).",
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
        short_help=f"Download LEB1 binary ephemeris for '{tier}' tier (~14x speedup).",
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


# --- LEB2 downloads ---

_LEB2_SIZES = {"base": "~28 MB", "medium": "~99 MB", "extended": "~734 MB"}


def _make_leb2_download(tier: str) -> click.Command:
    """Create a LEB2 download command for a tier."""

    @click.command(
        f"leb2-{tier}",
        short_help=f"Download LEB2 compressed ephemeris for '{tier}' tier ({_LEB2_SIZES.get(tier, '')}).",
        help=f"Download LEB2 compressed modular ephemeris for '{tier}' tier.\n\n"
        f"LEB2 uses error-bounded lossy compression (mantissa truncation + zstd)\n"
        f'to achieve 4-10x smaller files while maintaining <0.001" precision vs LEB1.\n\n'
        f"Downloads 4 group files: core, asteroids, apogee, uranians.\n"
        f"Total size: {_LEB2_SIZES.get(tier, 'varies')}.\n\n"
        f"Files are saved to ~/.libephemeris/leb/ by default.",
    )
    @_download_options
    def cmd(force: bool, no_progress: bool, quiet: bool) -> None:
        from .download import download_leb2_for_tier

        _handle_download(
            download_leb2_for_tier,
            quiet=quiet,
            tier_name=tier,
            force=force,
            show_progress=not no_progress,
            activate=True,
        )

    return cmd


for _tier in TIER_INFO:
    download_group.add_command(_make_leb2_download(_tier))


# --- ASSIST download ---


@download_group.command(
    short_help="Download ASSIST n-body data files (~714 MB, requires libephemeris[nbody]).",
)
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
# status command — comprehensive library and data overview
# ---------------------------------------------------------------------------


@cli.command(
    short_help="Show comprehensive library status: version, config, all data files.",
)
@click.option(
    "--json",
    "as_json",
    is_flag=True,
    help="Output status as machine-readable JSON.",
)
def status(as_json: bool) -> None:
    """Show comprehensive library and data file status.

    Displays version, calculation mode, precision tier, LEB file, data directory,
    and the status of all data files: DE kernels, planet center corrections,
    LEB1 binary ephemeris, LEB2 compressed files, SPK asteroid cache,
    ASSIST n-body data, and IERS Earth orientation data.

    \b
    Use --json for machine-readable output (e.g. for CI pipelines).
    """
    from .download import print_data_status

    print_data_status(as_json=as_json)


# ---------------------------------------------------------------------------
# info command (deprecated — alias for status --brief)
# ---------------------------------------------------------------------------


@cli.command(
    short_help="[Deprecated] Use 'status' instead.",
    deprecated=True,
)
def info() -> None:
    """Show library version, calc mode, and precision tier.

    DEPRECATED: This command is superseded by 'libephemeris status' which
    shows everything 'info' showed plus comprehensive data file status.
    """
    click.echo(
        click.style("Note: ", fg="yellow")
        + "'info' is deprecated. Use 'libephemeris status' instead."
    )
    click.echo()

    from .download import print_data_status

    print_data_status(as_json=False)


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
