"""Manual/documentation build commands: EPUB, PDF, pandoc and ebooklib workflows.

Replaces 8 poe tasks: manual:build*, docs:manual:generate*.
"""

from __future__ import annotations

import subprocess
import sys

import click


def _python(args: list[str]) -> None:
    """Run a python script."""
    sys.exit(subprocess.call([sys.executable, *args]))


@click.group(
    "manual",
    help="Build the libephemeris user manual in EPUB and/or PDF format.\n\n"
    "Two workflows are available:\n\n"
    "  Pandoc workflow (build*)      Requires: pandoc + tectonic\n"
    "  Ebooklib workflow (generate*) Requires: ebooklib + markdown (no pandoc)\n\n"
    "Both produce manuals in Italian and English. The ebooklib workflow\n"
    "generates Kobo-compatible EPUBs without external tools.\n\n"
    "  leph manual build             # EPUB + PDF, both languages\n"
    "  leph manual generate-epub     # EPUB only, no pandoc needed",
)
def manual_group() -> None:
    """Manual build commands."""


# ===========================================================================
# Pandoc + tectonic workflow (requires: brew install pandoc tectonic)
# ===========================================================================


@manual_group.command()
def build() -> None:
    """Build all manuals (EPUB + PDF, Italian + English).

    Requires: pandoc and tectonic (brew install pandoc tectonic).
    """
    _python(["scripts/build_manual.py"])


@manual_group.command("build-epub")
def build_epub() -> None:
    """Build all manuals in EPUB format only (Italian + English).

    Requires: pandoc (brew install pandoc).
    """
    _python(["scripts/build_manual.py", "--format", "epub"])


@manual_group.command("build-pdf")
def build_pdf() -> None:
    """Build all manuals in PDF format only (Italian + English).

    Requires: pandoc and tectonic (brew install pandoc tectonic).
    """
    _python(["scripts/build_manual.py", "--format", "pdf"])


@manual_group.command("build-it")
def build_it() -> None:
    """Build Italian manual (EPUB + PDF).

    Requires: pandoc and tectonic.
    """
    _python(["scripts/build_manual.py", "--lang", "it"])


@manual_group.command("build-en")
def build_en() -> None:
    """Build English manual (EPUB + PDF).

    Requires: pandoc and tectonic.
    """
    _python(["scripts/build_manual.py", "--lang", "en"])


# ===========================================================================
# ebooklib workflow (no pandoc required, Kobo-compatible)
# ===========================================================================


@manual_group.command("generate-epub")
def generate_epub() -> None:
    """Generate all manual EPUBs (Italian + English, no pandoc required).

    Uses ebooklib + markdown. Kobo-compatible output.
    Requires: pip install ebooklib markdown pyyaml
    """
    _python(["scripts/generate_manual_epub.py"])


@manual_group.command("generate-epub-it")
def generate_epub_it() -> None:
    """Generate Italian manual EPUB (no pandoc required)."""
    _python(["scripts/generate_manual_epub.py", "--lang", "it"])


@manual_group.command("generate-epub-en")
def generate_epub_en() -> None:
    """Generate English manual EPUB (no pandoc required)."""
    _python(["scripts/generate_manual_epub.py", "--lang", "en"])
