"""Code quality commands: lint, format, typecheck.

Replaces poe tasks: lint, format, format-black, typecheck.
"""

from __future__ import annotations

import subprocess
import sys

import click


@click.group(
    "code",
    help="Code quality tools: lint with ruff, format with ruff, type-check with mypy.\n\nRun before committing to catch issues early.",
)
def code_group() -> None:
    """Code quality tools."""


@code_group.command()
def lint() -> None:
    """Run ruff linter on the entire project with auto-fix enabled.

    Checks for style violations, unused imports, undefined names, and more.
    Safe fixes are applied automatically; unsafe fixes require manual review.
    """
    sys.exit(subprocess.call(["ruff", "check", ".", "--fix"]))


@code_group.command()
def format() -> None:
    """Format all Python files with the ruff formatter (line-length 88).

    Equivalent to black but faster. Formats libephemeris/, tests/, scripts/, etc.
    """
    sys.exit(subprocess.call(["ruff", "format", "."]))


@code_group.command("format-black")
def format_black() -> None:
    """Format with black (legacy, kept for compatibility -- prefer 'leph code format')."""
    sys.exit(subprocess.call(["black", "libephemeris", "tests"]))


@code_group.command()
def typecheck() -> None:
    """Run mypy static type checker on the libephemeris package.

    Checks type annotations, function signatures, and catches type errors
    at development time. Uses the mypy config from pyproject.toml.
    """
    sys.exit(subprocess.call(["mypy", "libephemeris"]))
