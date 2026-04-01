"""Data download commands: SPK kernels, LEB files, ASSIST data.

Replaces 7 poe tasks: spk:download:*, download:leb:*, download:assist.
"""

from __future__ import annotations

import subprocess
import sys

import click


def _python(args: list[str]) -> None:
    """Run a python expression or script."""
    sys.exit(
        subprocess.call(
            [sys.executable, *args],
        )
    )


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group(
    "download",
    help="Download data files needed by libephemeris.\n\n"
    "Three categories of data:\n\n"
    "  SPK kernels   NASA JPL binary ephemeris files (DE440/DE441) + asteroid SPKs\n"
    "  LEB files     Precomputed Chebyshev polynomial approximations (~14x speedup)\n"
    "  ASSIST data   REBOUND/ASSIST n-body integration data for TNOs/asteroids\n\n"
    "SPK files are required for the Skyfield backend. LEB files are optional\n"
    "but provide a major speedup. ASSIST data is only needed for n-body work.\n\n"
    "  leph download spk-medium    # DE440 + asteroid SPKs for medium tier\n"
    "  leph download leb-medium    # Pre-built LEB file (~20 MB)",
)
def download_group() -> None:
    """Data download commands."""


# ===========================================================================
# SPK downloads
# ===========================================================================


@download_group.command("spk-base")
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


@download_group.command("spk-medium")
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


@download_group.command("spk-extended")
def spk_extended() -> None:
    """Download max-range asteroid SPK files for extended tier from JPL Horizons.

    These are custom-generated SPK files with the widest date range available.
    Requires internet access to query JPL Horizons for each body.
    """
    _python(["scripts/download_max_range_spk.py"])


# ===========================================================================
# LEB downloads
# ===========================================================================


@download_group.command("leb-base")
def leb_base() -> None:
    """Download pre-built LEB binary ephemeris for base tier (~5 MB).

    Covers 1850-2150 (de440s). Provides ~14x speedup over Skyfield.
    Downloaded from GitHub Releases to data/leb/.
    """
    _python(["-m", "libephemeris.cli", "download:leb:base"])


@download_group.command("leb-medium")
def leb_medium() -> None:
    """Download pre-built LEB binary ephemeris for medium tier (~20 MB).

    Covers 1550-2650 (de440). This is the recommended tier for most users.
    Downloaded from GitHub Releases to data/leb/.
    """
    _python(["-m", "libephemeris.cli", "download:leb:medium"])


@download_group.command("leb-extended")
def leb_extended() -> None:
    """Download pre-built LEB binary ephemeris for extended tier (~180 MB).

    Covers -5000 to 5000 CE (de441). Largest file, widest date range.
    Downloaded from GitHub Releases to data/leb/.
    """
    _python(["-m", "libephemeris.cli", "download:leb:extended"])


# ===========================================================================
# ASSIST download
# ===========================================================================


@download_group.command()
def assist() -> None:
    """Download REBOUND/ASSIST n-body integration data (~120 MB).

    Downloads planet ephemeris and asteroid perturber SPK files needed
    for high-precision TNO/asteroid orbit propagation via REBOUND/ASSIST.
    Files are saved to ~/.libephemeris/assist/.
    Requires: pip install libephemeris[nbody]
    """
    _python(["-m", "libephemeris.cli", "download:assist"])
