"""Dev CLI for libephemeris — development tools and task runner.

Provides the ``leph`` command with hierarchical subcommands for testing,
LEB generation, diagnostics, data downloads, release management, and more.

Shell completion (add to your shell profile):

    # zsh
    eval "$(_LEPH_COMPLETE=zsh_source leph)"

    # bash
    eval "$(_LEPH_COMPLETE=bash_source leph)"

    # fish
    _LEPH_COMPLETE=fish_source leph | source
"""

from __future__ import annotations

import click

from .cmd_calibrate import calibrate_group
from .cmd_code import code_group
from .cmd_diag import diag_group
from .cmd_download import download_group
from .cmd_generate import generate_group
from .cmd_leb import leb_group
from .cmd_leb2 import leb2_group
from .cmd_manual import manual_group
from .cmd_release import release_group
from .cmd_test import test_group


@click.group(
    name="leph",
    help="libephemeris developer CLI — run tests, generate data, lint, release.\n\n"
    "Every subgroup is self-contained: use TAB completion or --help on any\n"
    "subcommand to discover what's available.\n\n"
    "Quick start:\n\n"
    "  leph test skyfield essential     # Fast sanity check (~490 tests, ~20s)\n"
    "  leph test leb-backend unit-fast  # Recommended daily driver (~1 min)\n"
    "  leph code lint                   # Ruff linter with auto-fix\n",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(package_name="libephemeris", prog_name="leph")
def cli() -> None:
    """Root CLI group."""


# Register all subgroups
cli.add_command(code_group)
cli.add_command(test_group)
cli.add_command(leb_group)
cli.add_command(leb2_group)
cli.add_command(download_group)
cli.add_command(diag_group)
cli.add_command(generate_group)
cli.add_command(calibrate_group)
cli.add_command(release_group)
cli.add_command(manual_group)


def main() -> None:
    """Entry point for the ``leph`` console script."""
    cli()
