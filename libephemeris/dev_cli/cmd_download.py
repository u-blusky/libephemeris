"""Data download commands for local libephemeris development.

Replaces 7 poe tasks: spk:download:*, download:leb:*, download:assist.
"""

from __future__ import annotations

import subprocess
import sys

import click


def _run_python(args: list[str]) -> int:
    """Run a Python module or script and return its exit code."""
    return subprocess.call([sys.executable, *args])


def _python(args: list[str]) -> None:
    """Run a Python module or script and exit with the same status."""
    sys.exit(_run_python(args))


def _libephemeris_download_args(
    *parts: str,
    force: bool,
    no_progress: bool,
    quiet: bool,
) -> list[str]:
    """Build a wrapped `libephemeris download ...` invocation."""
    args = ["-m", "libephemeris.cli", "download", *parts]
    if force:
        args.append("--force")
    if no_progress:
        args.append("--no-progress")
    if quiet:
        args.append("--quiet")
    return args


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group(
    "download",
    short_help="Download runtime and developer data files.",
    help="Download data files needed for libephemeris development.\n\n"
    "\b\n"
    "Categories:\n"
    "  all          Full developer dataset: runtime files + pre-built LEB1 + ASSIST\n"
    "  spk-*        DE kernels + asteroid SPKs for a specific tier\n"
    "  leb-*        Pre-built LEB1 files for verification and tooling\n"
    "  assist       REBOUND/ASSIST n-body files\n\n"
    "Use `leph download all` to prepare a machine for local development, data\n"
    "generation, verification, and packaging workflows.\n\n"
    "\b\n"
    "Examples:\n"
    "  leph download all\n"
    "  leph download spk-medium\n"
    "  leph download leb-medium",
)
def download_group() -> None:
    """Data download commands."""


@download_group.command(
    "all",
    short_help="Download the full local-development dataset for libephemeris.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Re-download files even if they are already present.",
)
@click.option(
    "--no-progress",
    is_flag=True,
    help="Disable progress output for wrapped libephemeris downloads.",
)
@click.option(
    "--quiet",
    is_flag=True,
    help="Suppress wrapper output and pass quiet mode through to subcommands.",
)
def download_all(force: bool, no_progress: bool, quiet: bool) -> None:
    """Download the full developer dataset used for local libephemeris work.

    This is the developer superset of `libephemeris download all`.

    \b
    Includes:
      1. All runtime data for every tier: DE kernels, planet centers,
         asteroid SPKs, IERS files, and LEB2 modules.
      2. Pre-built LEB1 files for base, medium, and extended tiers.
      3. ASSIST n-body data used by REBOUND/ASSIST integrations.

    This command is meant for contributors working on generation, verification,
    packaging, and release tasks. Some generation workflows still require extra
    tools or upstream network calls at generation time (for example `spiceypy`).
    """
    phases = [
        ("runtime dataset", ["all"]),
        ("LEB1 base", ["leb-base"]),
        ("LEB1 medium", ["leb-medium"]),
        ("LEB1 extended", ["leb-extended"]),
        ("ASSIST", ["assist"]),
    ]

    for label, command in phases:
        if not quiet:
            click.echo(f"\n== {label} ==")

        exit_code = _run_python(
            _libephemeris_download_args(
                *command,
                force=force,
                no_progress=no_progress,
                quiet=quiet,
            )
        )
        if exit_code != 0:
            sys.exit(exit_code)


# ===========================================================================
# SPK downloads
# ===========================================================================


@download_group.command(
    "spk-base",
    short_help="Download asteroid SPK kernels for base tier (1850-2150).",
)
def spk_base() -> None:
    """Download asteroid SPK21 kernel files for base tier (1850-2150).

    Downloads SPK kernels for Chiron, Ceres, Pallas, Juno, Vesta and other
    minor bodies. Also downloads DE440s if not already present.
    Files are saved to ~/.libephemeris/ by default.
    """
    _python(
        [
            "-c",
            "import libephemeris as eph; eph.set_precision_tier('base'); eph.ensure_all_ephemerides(show_progress=True)",
        ]
    )


@download_group.command(
    "spk-medium",
    short_help="Download asteroid SPK kernels for medium tier (1550-2650).",
)
def spk_medium() -> None:
    """Download asteroid SPK21 kernel files for medium tier (1550-2650).

    Downloads SPK kernels for all supported minor bodies plus DE440.
    """
    _python(
        [
            "-c",
            "import libephemeris as eph; eph.set_precision_tier('medium'); eph.ensure_all_ephemerides(show_progress=True)",
        ]
    )


@download_group.command(
    "spk-extended",
    short_help="Download max-range asteroid SPK files for extended tier.",
)
def spk_extended() -> None:
    """Download max-range asteroid SPK files for extended tier from JPL Horizons.

    These are custom-generated SPK files with the widest date range available.
    Requires internet access to query JPL Horizons for each body.
    """
    _python(["scripts/download_max_range_spk.py"])


# ===========================================================================
# LEB downloads
# ===========================================================================


@download_group.command(
    "leb-base",
    short_help="Download LEB binary ephemeris for base tier (~5 MB).",
)
def leb_base() -> None:
    """Download pre-built LEB binary ephemeris for base tier (~5 MB).

    Covers 1850-2150 (de440s). Provides ~14x speedup over Skyfield.
    Downloaded from GitHub Releases to data/leb/.
    """
    _python(["-m", "libephemeris.cli", "download:leb:base"])


@download_group.command(
    "leb-medium",
    short_help="Download LEB binary ephemeris for medium tier (~20 MB).",
)
def leb_medium() -> None:
    """Download pre-built LEB binary ephemeris for medium tier (~20 MB).

    Covers 1550-2650 (de440). This is the recommended tier for most users.
    Downloaded from GitHub Releases to data/leb/.
    """
    _python(["-m", "libephemeris.cli", "download:leb:medium"])


@download_group.command(
    "leb-extended",
    short_help="Download LEB binary ephemeris for extended tier (~180 MB).",
)
def leb_extended() -> None:
    """Download pre-built LEB binary ephemeris for extended tier (~180 MB).

    Covers -5000 to 5000 CE (de441). Largest file, widest date range.
    Downloaded from GitHub Releases to data/leb/.
    """
    _python(["-m", "libephemeris.cli", "download:leb:extended"])


# ===========================================================================
# ASSIST download
# ===========================================================================


@download_group.command(
    short_help="Download REBOUND/ASSIST n-body data (~120 MB).",
)
def assist() -> None:
    """Download REBOUND/ASSIST n-body integration data (~120 MB).

    Downloads planet ephemeris and asteroid perturber SPK files needed
    for high-precision TNO/asteroid orbit propagation via REBOUND/ASSIST.
    Files are saved to ~/.libephemeris/assist/.
    Requires: pip install libephemeris[nbody]
    """
    _python(["-m", "libephemeris.cli", "download:assist"])
