"""Status command for the developer CLI."""

from __future__ import annotations

import click


@click.command(
    "status",
    short_help="Show libephemeris data and configuration status.",
)
@click.option(
    "--json",
    "as_json",
    is_flag=True,
    help="Output status as machine-readable JSON.",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase detail: -v shows paths, -vv shows full file lists.",
)
def status(as_json: bool, verbose: int) -> None:
    """Show the same status report exposed by the end-user CLI."""
    from ..download import print_data_status

    print_data_status(as_json=as_json, verbose=verbose)
