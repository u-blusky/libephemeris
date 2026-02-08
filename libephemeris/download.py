"""
Data file download utilities for libephemeris.

This module provides functionality to download optional data files that enhance
precision for certain calculations. The main file is planet_centers.bsp which
provides precise planet center positions for outer planets (Jupiter-Pluto).

Usage:
    # From command line
    libephemeris download-data

    # From Python
    from libephemeris.download import download_planet_centers
    download_planet_centers()

Without the data files, libephemeris will fall back to analytical approximations
which are still accurate to ~0.1 arcseconds.
"""

from __future__ import annotations

import hashlib
import os
import sys
import tempfile
import urllib.request
from pathlib import Path
from typing import Callable, Optional

from .logging_config import get_logger

# GitHub Releases URL for data files
# The release tag "data-v1" can be updated when new data files are generated
GITHUB_RELEASES_BASE = (
    "https://github.com/g-battaglia/libephemeris/releases/download/data-v1"
)

# Data file definitions: (filename, sha256 hash, description)
DATA_FILES = {
    "planet_centers.bsp": {
        "url": f"{GITHUB_RELEASES_BASE}/planet_centers.bsp",
        "sha256": None,  # Will be set after first release
        "size_mb": 25.4,
        "description": "Precise planet center positions for Jupiter, Saturn, Uranus, Neptune, Pluto (1989-2049)",
    },
}


def get_data_dir() -> Path:
    """Get the data directory path.

    Returns the path to libephemeris/data/ within the package installation.
    Creates the directory if it doesn't exist.

    Returns:
        Path to the data directory
    """
    data_dir = Path(__file__).parent / "data"
    data_dir.mkdir(exist_ok=True)
    return data_dir


def _format_size(size_bytes: int) -> str:
    """Format byte size to human readable string."""
    for unit in ["B", "KB", "MB", "GB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} TB"


class SimpleProgressBar:
    """Simple progress bar that works without external dependencies."""

    def __init__(
        self,
        total: int,
        description: str = "",
        width: int = 40,
        file=None,
    ):
        self.total = total
        self.description = description
        self.width = width
        self.current = 0
        self.file = file or sys.stdout
        self._last_percent = -1

    def update(self, amount: int = 1):
        """Update progress by amount."""
        self.current += amount
        self._render()

    def _render(self):
        """Render the progress bar."""
        if self.total == 0:
            return

        percent = int(100 * self.current / self.total)

        # Only update if percent changed (avoid too many writes)
        if percent == self._last_percent:
            return
        self._last_percent = percent

        filled = int(self.width * self.current / self.total)
        bar = "=" * filled + "-" * (self.width - filled)

        # Format: Description [====----] 50% (12.5 MB / 25.0 MB)
        current_str = _format_size(self.current)
        total_str = _format_size(self.total)

        line = (
            f"\r{self.description} [{bar}] {percent:3d}% ({current_str} / {total_str})"
        )
        self.file.write(line)
        self.file.flush()

    def close(self):
        """Finish the progress bar."""
        self.file.write("\n")
        self.file.flush()


def _get_progress_bar(total: int, description: str) -> SimpleProgressBar:
    """Get the best available progress bar.

    Tries to use rich if available, falls back to simple progress bar.

    Args:
        total: Total size in bytes
        description: Description to show

    Returns:
        Progress bar object with update() and close() methods
    """
    # Try rich first (best experience)
    try:
        from rich.progress import (
            Progress,
            BarColumn,
            DownloadColumn,
            TransferSpeedColumn,
            TimeRemainingColumn,
        )

        class RichProgressWrapper:
            def __init__(self, total: int, description: str):
                self.progress = Progress(
                    "[progress.description]{task.description}",
                    BarColumn(),
                    "[progress.percentage]{task.percentage:>3.0f}%",
                    DownloadColumn(),
                    TransferSpeedColumn(),
                    TimeRemainingColumn(),
                )
                self.progress.start()
                self.task = self.progress.add_task(description, total=total)

            def update(self, amount: int):
                self.progress.update(self.task, advance=amount)

            def close(self):
                self.progress.stop()

        return RichProgressWrapper(total, description)
    except ImportError:
        pass

    # Fall back to simple progress bar
    return SimpleProgressBar(total, description)


def download_file(
    url: str,
    dest_path: Path,
    description: str = "Downloading",
    expected_sha256: Optional[str] = None,
    show_progress: bool = True,
) -> bool:
    """Download a file with progress bar.

    Args:
        url: URL to download from
        dest_path: Destination path for the file
        description: Description to show in progress bar
        expected_sha256: Expected SHA256 hash (optional, for verification)
        show_progress: Whether to show progress bar

    Returns:
        True if download successful, False otherwise

    Raises:
        urllib.error.URLError: If download fails
        ValueError: If hash verification fails
    """
    logger = get_logger()

    # Use a temporary file for atomic download
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_fd, temp_path = tempfile.mkstemp(dir=dest_path.parent, suffix=".download")

    try:
        # Open URL and get content length
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "libephemeris-download/1.0"},
        )

        logger.info("Downloading %s...", description)

        with urllib.request.urlopen(req, timeout=30) as response:
            total_size = int(response.headers.get("Content-Length", 0))

            if total_size > 0:
                size_mb = total_size / (1024 * 1024)
                logger.info("File size: %.1f MB", size_mb)

            # Setup progress bar
            if show_progress and total_size > 0:
                progress = _get_progress_bar(total_size, description)
            else:
                progress = None

            # Download with progress updates
            sha256 = hashlib.sha256()
            chunk_size = 64 * 1024  # 64KB chunks

            with os.fdopen(temp_fd, "wb") as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    sha256.update(chunk)
                    if progress:
                        progress.update(len(chunk))

            if progress:
                progress.close()

        # Verify hash if provided
        if expected_sha256:
            actual_hash = sha256.hexdigest()
            if actual_hash != expected_sha256:
                os.unlink(temp_path)
                raise ValueError(
                    f"Hash mismatch: expected {expected_sha256}, got {actual_hash}"
                )

        # Atomic move to final destination
        os.replace(temp_path, dest_path)
        logger.info("Download complete: %s", dest_path.name)
        return True

    except Exception:
        # Clean up temp file on error
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        raise


def download_planet_centers(
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
) -> Path:
    """Download the planet_centers.bsp data file.

    This file provides precise planet center positions for Jupiter, Saturn,
    Uranus, Neptune, and Pluto. It enables sub-arcsecond precision for
    outer planet calculations.

    Args:
        force: If True, download even if file already exists
        show_progress: If True, show progress bar during download
        quiet: If True, suppress all output except errors

    Returns:
        Path to the downloaded file

    Raises:
        urllib.error.URLError: If download fails
        ValueError: If hash verification fails
    """
    file_info = DATA_FILES["planet_centers.bsp"]
    dest_path = get_data_dir() / "planet_centers.bsp"

    # Check if already exists
    if dest_path.exists() and not force:
        if not quiet:
            print(f"planet_centers.bsp already exists at {dest_path}")
            print("Use --force to re-download.")
        return dest_path

    if not quiet:
        print(f"Downloading planet_centers.bsp (~{file_info['size_mb']:.1f} MB)...")
        print(f"  {file_info['description']}")
        print()

    try:
        download_file(
            url=file_info["url"],
            dest_path=dest_path,
            description="planet_centers.bsp",
            expected_sha256=file_info.get("sha256"),
            show_progress=show_progress,
        )

        if not quiet:
            print()
            print(f"Downloaded to: {dest_path}")
            print(
                "Planet center data is now available for high-precision calculations."
            )

        return dest_path

    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(
                "\nError: Data file not found on server.",
                file=sys.stderr,
            )
            print(
                "The data release may not be available yet.",
                file=sys.stderr,
            )
            print(
                "libephemeris will use analytical approximations (~0.1 arcsec precision).",
                file=sys.stderr,
            )
        raise
    except Exception as e:
        print(f"\nError downloading file: {e}", file=sys.stderr)
        raise


def download_all(
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
) -> list[Path]:
    """Download all optional data files.

    Args:
        force: If True, download even if files already exist
        show_progress: If True, show progress bar during download
        quiet: If True, suppress all output except errors

    Returns:
        List of paths to downloaded files
    """
    paths = []
    paths.append(
        download_planet_centers(force=force, show_progress=show_progress, quiet=quiet)
    )
    return paths


def check_data_status() -> dict[str, dict]:
    """Check the status of all data files.

    Returns:
        Dictionary with file names as keys and status info as values.
        Each status dict contains:
        - exists: bool
        - path: Path or None
        - size: int or None (file size in bytes)
        - description: str
    """
    data_dir = get_data_dir()
    status = {}

    for filename, info in DATA_FILES.items():
        path = data_dir / filename
        exists = path.exists()
        status[filename] = {
            "exists": exists,
            "path": path if exists else None,
            "size": path.stat().st_size if exists else None,
            "description": info["description"],
            "expected_size_mb": info["size_mb"],
        }

    return status


def print_data_status():
    """Print the status of all data files to stdout."""
    status = check_data_status()

    print("libephemeris data files:")
    print()

    for filename, info in status.items():
        if info["exists"]:
            size_str = _format_size(info["size"])
            print(f"  [OK] {filename} ({size_str})")
            print(f"        {info['description']}")
        else:
            print(f"  [--] {filename} (not installed)")
            print(f"        {info['description']}")
            print(f"        Expected size: ~{info['expected_size_mb']:.1f} MB")

    print()
    print("Run 'libephemeris download-data' to download missing files.")
