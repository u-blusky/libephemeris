"""
Data file download utilities for libephemeris.

This module provides functionality to download optional data files that enhance
precision for certain calculations. The main file is planet_centers.bsp which
provides precise planet center positions for outer planets (Jupiter-Pluto).

Usage:
    # From command line
    libephemeris download:medium       Download data for 'medium' tier
    libephemeris download:base         Download data for 'base' tier
    libephemeris download:extended     Download data for 'extended' tier

    # From Python
    from libephemeris.download import download_for_tier
    download_for_tier("medium")

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
from typing import Any, Callable, Optional, Union

from .logging_config import get_logger


def _is_valid_bsp(filepath: str) -> bool:
    """Check if a BSP file can be opened and its segments enumerated by jplephem.

    Uses jplephem.spk.SPK.open() which validates:
    - DAF header (magic bytes, endianness)
    - FTP corruption test string
    - All segment summary records are readable

    This catches truncated files where the header is intact but segment
    records are missing deeper in the file.

    Args:
        filepath: Path to the BSP file to validate

    Returns:
        True if file can be opened and all segments enumerated, False otherwise
    """
    try:
        from jplephem.spk import SPK

        with SPK.open(filepath) as kernel:
            # Force iteration through all segments to detect truncation
            for segment in kernel.segments:
                _ = segment.center, segment.target
        return True
    except Exception:
        return False


# GitHub Releases URL for data files
# The release tag "data-v1" can be updated when new data files are generated
GITHUB_RELEASES_BASE = (
    "https://github.com/g-battaglia/libephemeris/releases/download/data-v1"
)

# Data file definitions: (filename, sha256 hash, description)
DATA_FILES = {
    # Legacy single file (kept for backward compatibility)
    "planet_centers.bsp": {
        "url": f"{GITHUB_RELEASES_BASE}/planet_centers.bsp",
        "sha256": None,  # Will be set after first release
        "size_mb": 25.4,
        "description": "Precise planet center positions for Jupiter, Saturn, Uranus, Neptune, Pluto (1989-2049)",
    },
    # Tier-specific files
    "planet_centers_base.bsp": {
        "url": f"{GITHUB_RELEASES_BASE}/planet_centers_base.bsp",
        "sha256": "a9ec744ff412b095129166587ea0814f81c850faebf92586a738cb5dc103c92a",
        "size_mb": 25.4,
        "description": "Planet centers for 'base' tier (1850-2150)",
    },
    "planet_centers_medium.bsp": {
        "url": f"{GITHUB_RELEASES_BASE}/planet_centers_medium.bsp",
        "sha256": "b4fd366f2d00958ee3dd4a8884164d339674d5f7f2ea25c5f6959705c2b66852",
        "size_mb": 72.6,
        "description": "Planet centers for 'medium' tier (1550-2650)",
    },
    "planet_centers_extended.bsp": {
        "url": f"{GITHUB_RELEASES_BASE}/planet_centers_extended.bsp",
        "sha256": "a07b046b89a9992fc7fda445b00e656341a3bab66a035adb8108de7d4bd69edc",
        "size_mb": 222.6,
        "description": "Planet centers for 'extended' tier (partial -12000 to +17000)",
    },
    # LEB (LibEphemeris Binary) precomputed ephemeris files
    "ephemeris_base.leb": {
        "url": f"{GITHUB_RELEASES_BASE}/ephemeris_base.leb",
        "sha256": "0313211cabd17b398a1182673f164d333a24500314c2bd2e979b5a1a40f5b0fd",
        "size_mb": 53.1,
        "description": "LEB binary ephemeris for 'base' tier (1850-2150, ~14x speedup)",
    },
    "ephemeris_medium.leb": {
        "url": f"{GITHUB_RELEASES_BASE}/ephemeris_medium.leb",
        "sha256": "3cf0665f06e175ed1f8b07febbdd47bf273e1d00527f14cb6e20cc58fd810791",
        "size_mb": 174.6,
        "description": "LEB binary ephemeris for 'medium' tier (1550-2650, ~14x speedup)",
    },
    "ephemeris_extended.leb": {
        "url": f"{GITHUB_RELEASES_BASE}/ephemeris_extended.leb",
        "sha256": None,  # Not yet generated
        "size_mb": None,
        "description": "LEB binary ephemeris for 'extended' tier (-5000 to +5000, ~14x speedup)",
    },
}


def get_data_dir() -> Path:
    """Get the data directory path.

    Returns the path to ~/.libephemeris (or LIBEPHEMERIS_DATA_DIR env var).
    Creates the directory if it doesn't exist.

    Returns:
        Path to the data directory
    """
    from .state import _get_data_dir

    return Path(_get_data_dir())


def _format_size(size_bytes: float) -> str:
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


def _get_progress_bar(total: int, description: str) -> Any:
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
    logger = get_logger()

    # Check if already exists and valid
    if dest_path.exists() and not force:
        if _is_valid_bsp(str(dest_path)):
            if not quiet:
                print(f"planet_centers.bsp already exists at {dest_path}")
                print("Use --force to re-download.")
            return dest_path
        logger.warning("Cached file %s is corrupted, re-downloading", dest_path)
        try:
            os.remove(dest_path)
        except OSError:
            pass

    if not quiet:
        print(f"Downloading planet_centers.bsp (~{file_info['size_mb']:.1f} MB)...")
        print(f"  {file_info['description']}")
        print()

    try:
        download_file(
            url=str(file_info["url"]),
            dest_path=dest_path,
            description="planet_centers.bsp",
            expected_sha256=str(file_info["sha256"])
            if file_info.get("sha256")
            else None,
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
    workspace_root = Path(__file__).parent.parent
    status = {}

    for filename, info in DATA_FILES.items():
        # Tier-specific files are in workspace root, legacy in data dir
        if filename.startswith("planet_centers_"):
            path = workspace_root / filename
        else:
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
    from .state import get_precision_tier

    status = check_data_status()
    current_tier = get_precision_tier()

    print("libephemeris data files:")
    print()

    # Show tier-specific planet_centers first
    tier_files = [
        "planet_centers_base.bsp",
        "planet_centers_medium.bsp",
        "planet_centers_extended.bsp",
    ]
    for filename in tier_files:
        if filename not in status:
            continue
        info = status[filename]
        tier_name = filename.replace("planet_centers_", "").replace(".bsp", "")
        marker = " *" if tier_name == current_tier else ""
        if info["exists"]:
            size_str = _format_size(info["size"])
            print(f"  [OK] {filename} ({size_str}){marker}")
            print(f"        {info['description']}")
        else:
            print(f"  [--] {filename} (not installed){marker}")
            print(f"        {info['description']}")

    # Show legacy planet_centers.bsp
    legacy_file = "planet_centers.bsp"
    if legacy_file in status:
        info = status[legacy_file]
        print()
        if info["exists"]:
            size_str = _format_size(info["size"])
            print(f"  [OK] {legacy_file} ({size_str}) - legacy")
            print(f"        {info['description']}")
        else:
            print(f"  [--] {legacy_file} (not installed) - legacy")

    print()
    print(f"Current tier: {current_tier}")
    print()
    print("Run 'libephemeris download:<tier>' to download all data files.")
    print("Available tiers: base, medium, extended")


def init_all(
    cache_dir: Optional[str] = None,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
    start_year: int = 1550,
    end_year: int = 2650,
) -> dict:
    """Initialize libephemeris with all required data files.

    Downloads:
    1. DE440.bsp planetary ephemeris (~128 MB)
    2. planet_centers.bsp precision data (~25 MB)
    3. SPK kernels for all minor bodies defined in SPK_BODY_NAME_MAP
       (default 1550-2650, 20-year chunks)

    Args:
        cache_dir: Custom SPK cache directory (default: ~/.libephemeris/spk/).
            Also configurable via LIBEPHEMERIS_SPK_DIR env var or
            set_spk_cache_dir().
        force: Re-download existing files
        show_progress: Show progress output
        quiet: Suppress non-error output
        start_year: First year of SPK coverage (default: 1550)
        end_year: Last year of SPK coverage (default: 2650)

    Returns:
        dict with summary:
            - de440: bool (success)
            - planet_centers: bool (success)
            - spk_success: int (number of SPK chunks downloaded)
            - spk_failed: int (number of SPK chunks that failed)
            - spk_skipped: int (number of SPK chunks already cached)
            - spk_out_of_range: int (chunks outside body's Horizons SPK limits)
            - failed_details: list of (body_name, chunk_start, chunk_end, error)
    """
    import time

    logger = get_logger()

    result: dict[str, Any] = {
        "de440": False,
        "planet_centers": False,
        "spk_success": 0,
        "spk_failed": 0,
        "spk_skipped": 0,
        "spk_out_of_range": 0,
        "failed_details": [],
    }

    # -------------------------------------------------------------------------
    # Step 1: DE440.bsp
    # -------------------------------------------------------------------------
    if not quiet:
        print("Step 1/3: DE440.bsp planetary ephemeris")

    try:
        from .state import get_planets

        get_planets()  # Downloads if not present
        result["de440"] = True
        if not quiet:
            print("  [OK] DE440.bsp ready")
    except Exception as e:
        logger.error("Failed to initialize DE440: %s", e)
        if not quiet:
            print(f"  [FAIL] DE440.bsp: {e}", file=sys.stderr)

    # -------------------------------------------------------------------------
    # Step 2: planet_centers.bsp
    # -------------------------------------------------------------------------
    if not quiet:
        print("Step 2/3: planet_centers.bsp precision data")

    try:
        download_planet_centers(
            force=force,
            show_progress=show_progress,
            quiet=quiet,
        )
        result["planet_centers"] = True
        if not quiet:
            print("  [OK] planet_centers.bsp ready")
    except Exception as e:
        logger.error("Failed to download planet_centers.bsp: %s", e)
        if not quiet:
            print(f"  [FAIL] planet_centers.bsp: {e}", file=sys.stderr)

    # -------------------------------------------------------------------------
    # Step 3: SPK for all minor bodies
    # -------------------------------------------------------------------------
    if not quiet:
        print("Step 3/3: SPK kernels for minor bodies")

    from .constants import SPK_BODY_NAME_MAP
    from .spk_auto import (
        download_spk_from_horizons,
        _generate_spk_cache_filename,
        _iso_to_jd,
        ensure_cache_dir,
    )
    from .spk import _get_body_name

    effective_cache_dir = ensure_cache_dir(cache_dir)

    if not quiet:
        print(f"  Cache directory: {effective_cache_dir}")

    CHUNK_SIZE_YEARS = 20
    START_YEAR = start_year
    END_YEAR = end_year
    chunks_per_body = (END_YEAR - START_YEAR) // CHUNK_SIZE_YEARS
    total_bodies = len(SPK_BODY_NAME_MAP)

    if not quiet:
        print(
            f"  Bodies: {total_bodies} | "
            f"Chunks per body: {chunks_per_body} | "
            f"Total: {chunks_per_body * total_bodies}"
        )
        print()

    for body_idx, (ipl, (horizons_id, naif_id)) in enumerate(SPK_BODY_NAME_MAP.items()):
        body_name = _get_body_name(ipl) or f"body_{ipl}"

        if not quiet:
            print(
                f"  [{body_idx + 1}/{total_bodies}] {body_name} "
                f"(Horizons: {horizons_id})"
            )

        body_success = 0
        body_skipped = 0
        body_failed = 0
        body_out_of_range = 0
        chunk_idx = 0

        for chunk_start in range(START_YEAR, END_YEAR, CHUNK_SIZE_YEARS):
            chunk_end = min(chunk_start + CHUNK_SIZE_YEARS, END_YEAR)
            chunk_idx += 1

            jd_start = _iso_to_jd(f"{chunk_start}-01-01")
            jd_end = _iso_to_jd(f"{chunk_end}-01-01")

            filename = _generate_spk_cache_filename(horizons_id, jd_start, jd_end)
            output_path = os.path.join(effective_cache_dir, filename)

            # Show per-chunk progress
            if not quiet:
                sys.stdout.write(
                    f"\r    chunk {chunk_idx}/{chunks_per_body} "
                    f"({chunk_start}-{chunk_end}) ... "
                )
                sys.stdout.flush()

            if os.path.exists(output_path) and not force:
                body_skipped += 1
                result["spk_skipped"] += 1
                continue

            try:
                download_spk_from_horizons(
                    body_id=horizons_id,
                    jd_start=jd_start,
                    jd_end=jd_end,
                    output_path=output_path,
                    ipl=ipl,
                    naif_id=naif_id,
                )
                body_success += 1
                result["spk_success"] += 1

                # Rate limiting: small delay between Horizons requests
                time.sleep(1.5)

            except KeyboardInterrupt:
                if not quiet:
                    sys.stdout.write("\n")
                raise
            except Exception as e:
                error_msg = str(e)

                # Detect Horizons SPK limit errors to skip gracefully
                if "START time outside set SPK limits" in error_msg:
                    # Body data doesn't go back this far; try later chunks
                    body_out_of_range += 1
                    result["spk_out_of_range"] += 1
                    logger.debug(
                        "%s: chunk %d-%d before SPK start, skipping",
                        body_name,
                        chunk_start,
                        chunk_end,
                    )
                    continue

                if "STOP time outside set SPK limits" in error_msg:
                    # Body data ends before this chunk; skip remaining chunks
                    remaining = (END_YEAR - chunk_start) // CHUNK_SIZE_YEARS
                    body_out_of_range += remaining
                    result["spk_out_of_range"] += remaining
                    logger.debug(
                        "%s: chunk %d-%d past SPK end, skipping %d remaining",
                        body_name,
                        chunk_start,
                        chunk_end,
                        remaining,
                    )
                    break

                # Genuine failure
                body_failed += 1
                result["spk_failed"] += 1
                result["failed_details"].append(
                    (body_name, chunk_start, chunk_end, error_msg)
                )
                logger.warning(
                    "Failed to download %s chunk %d-%d: %s",
                    body_name,
                    chunk_start,
                    chunk_end,
                    e,
                )

        if not quiet:
            # Clear the per-chunk line and print summary
            sys.stdout.write("\r" + " " * 60 + "\r")
            parts = []
            if body_success:
                parts.append(f"{body_success} downloaded")
            if body_skipped:
                parts.append(f"{body_skipped} cached")
            if body_out_of_range:
                parts.append(f"{body_out_of_range} out of range")
            if body_failed:
                parts.append(f"{body_failed} failed")
            print(f"    {', '.join(parts)}")

    return result


def download_for_tier(
    tier_name: str,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
) -> dict:
    """Download all data files required for a specific precision tier.

    Sets the precision tier, then downloads:
    1. The tier's ephemeris file (de440s / de440 / de441)
    2. planet_centers.bsp precision offsets
    3. SPK kernels for all 21 minor bodies (full tier date range)

    For the 'extended' tier, SPK files are downloaded as single max-range
    files (1600-2500) rather than chunked, since that's the full Horizons
    range for minor bodies.

    Args:
        tier_name: One of "base", "medium", "extended"
        force: If True, re-download even if files already exist
        show_progress: If True, show progress output
        quiet: If True, suppress non-error output

    Returns:
        Dict with summary (same format as ensure_all_ephemerides)

    Raises:
        ValueError: If tier_name is not valid
    """
    from .state import TIERS, set_precision_tier

    if tier_name not in TIERS:
        valid = ", ".join(sorted(TIERS.keys()))
        raise ValueError(f"Unknown tier '{tier_name}'. Valid tiers: {valid}")

    tier = TIERS[tier_name]

    if not quiet:
        print(f"libephemeris: downloading data for tier '{tier_name}'")
        print(f"  Ephemeris file: {tier.ephemeris_file}")
        print(f"  Description:    {tier.description}")
        print(f"  SPK date range: {tier.spk_date_range[0]} to {tier.spk_date_range[1]}")
        print()

    # Activate the tier so ensure_all_ephemerides uses the right config
    set_precision_tier(tier_name)

    # Step 1: planet_centers file for this tier (workspace root)
    if not quiet:
        print("Step 1/2: planet_centers precision data")

    try:
        _download_planet_centers_for_tier(
            tier_name=tier_name,
            force=force,
            show_progress=show_progress,
            quiet=quiet,
        )
        if not quiet:
            print("  [OK] planet_centers ready")
            print()
    except Exception as e:
        if not quiet:
            print(f"  [WARN] planet_centers: {e}", file=sys.stderr)
            print("  (non-critical, continuing...)")
            print()

    # Step 2: Ephemeris file + SPK kernels (via ensure_all_ephemerides)
    if not quiet:
        print(f"Step 2/2: {tier.ephemeris_file} + SPK kernels for minor bodies")

    from .spk_auto import ensure_all_ephemerides

    results = ensure_all_ephemerides(
        force_download=force,
        show_progress=show_progress and not quiet,
    )

    # Print final summary
    if not quiet:
        summary = results.get("summary", {})
        print()
        print("Download complete:")
        print(f"  Tier:         {tier_name}")
        print(f"  Ephemeris:    {tier.ephemeris_file}")
        if isinstance(summary, dict):
            print(f"  SPK cached:   {summary.get('cached', 0)}")
            print(f"  SPK downloaded: {summary.get('downloaded', 0)}")
            if summary.get("fallback", 0):
                print(f"  SPK fallback: {summary['fallback']}")
            if summary.get("errors", 0):
                print(f"  SPK errors:   {summary['errors']}")

    return results


def _is_valid_leb(filepath: str) -> bool:
    """Check if a LEB file can be opened and parsed by LEBReader.

    Validates magic bytes, version, section directory, and body index.

    Args:
        filepath: Path to the LEB file to validate

    Returns:
        True if file can be opened and parsed, False otherwise
    """
    try:
        from .leb_reader import LEBReader

        reader = LEBReader(filepath)
        reader.close()
        return True
    except Exception:
        return False


def get_leb_dir() -> Path:
    """Get the LEB data directory path.

    Returns the path to ~/.libephemeris/leb/ (or under LIBEPHEMERIS_DATA_DIR).
    Creates the directory if it doesn't exist.

    Returns:
        Path to the LEB data directory
    """
    leb_dir = get_data_dir() / "leb"
    leb_dir.mkdir(parents=True, exist_ok=True)
    return leb_dir


def get_leb_path_for_tier(tier_name: str) -> Path:
    """Get the expected LEB file path for a tier.

    Args:
        tier_name: One of "base", "medium", "extended"

    Returns:
        Path where the LEB file should be stored
    """
    return get_leb_dir() / f"ephemeris_{tier_name}.leb"


def download_leb_for_tier(
    tier_name: str,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
    activate: bool = True,
) -> Path:
    """Download the precomputed LEB binary ephemeris file for a specific tier.

    Downloads the .leb file from GitHub Releases to ~/.libephemeris/leb/.
    After a successful download, optionally activates it via set_leb_file().

    LEB files contain precomputed Chebyshev polynomial approximations for
    all celestial bodies, providing ~14x speedup over the Skyfield pipeline.

    Args:
        tier_name: One of "base", "medium", "extended"
        force: If True, re-download even if file already exists
        show_progress: If True, show progress bar during download
        quiet: If True, suppress all output except errors
        activate: If True, call set_leb_file() after successful download

    Returns:
        Path to the downloaded LEB file

    Raises:
        ValueError: If tier_name is invalid or hash verification fails
        urllib.error.URLError: If download fails
        RuntimeError: If the extended tier LEB is not yet available
    """
    logger = get_logger()
    filename = f"ephemeris_{tier_name}.leb"
    file_info = DATA_FILES.get(filename)

    if file_info is None:
        valid_tiers = ["base", "medium", "extended"]
        raise ValueError(
            f"Unknown tier '{tier_name}'. Valid tiers: {', '.join(valid_tiers)}"
        )

    # Extended tier is not yet generated
    if file_info.get("sha256") is None:
        raise RuntimeError(
            f"LEB file for '{tier_name}' tier is not yet available for download. "
            f"You can generate it locally with: poe leb:generate:{tier_name}:groups"
        )

    dest_path = get_leb_path_for_tier(tier_name)

    # Check if already exists and valid
    if dest_path.exists() and not force:
        if _is_valid_leb(str(dest_path)):
            if not quiet:
                print(f"  {filename} already exists at {dest_path}")
                print("  Use --force to re-download.")
            # Activate if requested
            if activate:
                from .state import set_leb_file

                set_leb_file(str(dest_path))
            return dest_path
        else:
            logger.warning("Cached LEB file %s is corrupted, re-downloading", dest_path)
            try:
                os.remove(dest_path)
            except OSError:
                pass

    size_mb = file_info.get("size_mb", 0)
    if not quiet:
        print(f"  Downloading {filename} (~{size_mb:.0f} MB)...")
        print(f"  {file_info['description']}")
        print()

    try:
        download_file(
            url=str(file_info["url"]),
            dest_path=dest_path,
            description=filename,
            expected_sha256=str(file_info["sha256"]),
            show_progress=show_progress,
        )
    except Exception:
        # Clean up partial file
        if dest_path.exists():
            try:
                os.remove(dest_path)
            except OSError:
                pass
        raise

    # Validate the downloaded file
    if not _is_valid_leb(str(dest_path)):
        try:
            os.remove(dest_path)
        except OSError:
            pass
        raise ValueError(
            f"Downloaded {filename} failed LEB validation — corrupt or incompatible"
        )

    if not quiet:
        print()
        print(f"  Downloaded to: {dest_path}")
        print(f"  LEB binary ephemeris for '{tier_name}' tier is now available.")
        if activate:
            print("  Activated automatically for this session.")

    # Activate if requested
    if activate:
        from .state import set_leb_file

        set_leb_file(str(dest_path))
        logger.info("Activated LEB file: %s", dest_path)

    return dest_path


def _download_planet_centers_for_tier(
    tier_name: str,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
) -> Path:
    """Download planet_centers SPK file for a specific tier.

    Saves to data directory (~/.libephemeris).

    Args:
        tier_name: One of "base", "medium", "extended"
        force: If True, re-download even if file exists
        show_progress: If True, show progress bar
        quiet: If True, suppress output

    Returns:
        Path to the downloaded file

    Raises:
        Exception: If download fails
    """
    import urllib.request

    logger = get_logger()
    filename = f"planet_centers_{tier_name}.bsp"
    file_info = DATA_FILES.get(filename)

    if file_info is None:
        raise ValueError(f"No planet_centers file for tier '{tier_name}'")

    dest_path = get_data_dir() / filename

    if dest_path.exists() and not force:
        if _is_valid_bsp(str(dest_path)):
            if not quiet:
                print(f"  {filename} already exists at {dest_path}")
            return dest_path
        else:
            logger.warning("Cached file %s is corrupted, re-downloading", dest_path)
            try:
                os.remove(dest_path)
            except OSError:
                pass

    url = file_info["url"]
    expected_size_mb = file_info.get("size_mb", 50)

    if not quiet:
        print(f"  Downloading {filename} (~{expected_size_mb:.0f} MB)...")

    dest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_fd, temp_path = tempfile.mkstemp(dir=dest_path.parent, suffix=".download")

    try:
        total_size = 0
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=30) as response:
                total_size = int(response.headers.get("Content-Length", 0))
        except Exception:
            pass

        downloaded = 0
        chunk_size = 1024 * 1024

        with urllib.request.urlopen(url, timeout=300) as response:
            with os.fdopen(temp_fd, "wb") as f:
                if show_progress and total_size > 0 and not quiet:
                    progress = SimpleProgressBar(total_size, f"  {filename}")
                    while True:
                        chunk = response.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        progress.update(len(chunk))
                    print()
                else:
                    while True:
                        chunk = response.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        if not quiet:
                            print(
                                f"\r  Downloaded {downloaded / 1024 / 1024:.1f} MB",
                                end="",
                            )

        if not _is_valid_bsp(temp_path):
            os.unlink(temp_path)
            raise ValueError(
                f"Downloaded file {filename} failed validation - corrupt or incomplete"
            )

        os.replace(temp_path, dest_path)

        if not quiet:
            print(f"  Downloaded to {dest_path}")

        return dest_path

    except Exception:
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        raise
