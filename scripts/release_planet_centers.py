#!/usr/bin/env python3
"""
Upload planet_centers SPK files to GitHub Releases.

This script calculates SHA256 hashes and uploads planet_centers files
to a GitHub release using the gh CLI.

Requirements:
    - gh CLI installed and authenticated (gh auth login)

Usage:
    # Upload all planet_centers files
    python scripts/release_planet_centers.py --version 0.15.0

    # Upload specific tier
    python scripts/release_planet_centers.py --version 0.15.0 --tier extended

    # Dry run (show what would be uploaded)
    python scripts/release_planet_centers.py --version 0.15.0 --dry-run

Files are uploaded to:
    https://github.com/g-battaglia/libephemeris/releases/download/data-v1/
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
from pathlib import Path


TIER_FILES = {
    "base": "planet_centers_base.bsp",
    "medium": "planet_centers_medium.bsp",
    "extended": "planet_centers_extended.bsp",
    "legacy": "planet_centers.bsp",
}

DEFAULT_TIERS = ["base", "medium", "extended"]

RELEASE_TAG = "data-v1"
REPO = "g-battaglia/libephemeris"


def get_data_dir() -> Path:
    """Get the libephemeris data directory.

    Checks LIBEPHEMERIS_DATA_DIR env var, falls back to ~/.libephemeris.
    """
    env_dir = os.environ.get("LIBEPHEMERIS_DATA_DIR")
    if env_dir:
        return Path(env_dir)
    return Path.home() / ".libephemeris"


def find_bsp_file(filename: str) -> Path | None:
    """Find a BSP file in data dir or repo root.

    Search order:
        1. Data directory (~/.libephemeris or LIBEPHEMERIS_DATA_DIR)
        2. Repository root (fallback)

    Args:
        filename: BSP filename to find.

    Returns:
        Path to the file, or None if not found.
    """
    data_path = get_data_dir() / filename
    if data_path.exists():
        return data_path

    repo_root = Path(__file__).parent.parent
    repo_path = repo_root / filename
    if repo_path.exists():
        return repo_path

    return None


def calculate_sha256(filepath: Path) -> str:
    """Calculate SHA256 hash of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def get_file_size_mb(filepath: Path) -> float:
    """Get file size in MB."""
    return filepath.stat().st_size / (1024 * 1024)


def check_gh_auth() -> bool:
    """Check if gh CLI is authenticated."""
    try:
        result = subprocess.run(
            ["gh", "auth", "status"],
            capture_output=True,
            text=True,
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def check_release_exists(tag: str) -> bool:
    """Check if a release with the given tag exists."""
    try:
        result = subprocess.run(
            ["gh", "release", "view", tag, "--repo", REPO],
            capture_output=True,
            text=True,
        )
        return result.returncode == 0
    except Exception:
        return False


def upload_file(filepath: Path, tag: str, dry_run: bool = False) -> bool:
    """Upload a file to a GitHub release."""
    if dry_run:
        print(f"  [DRY RUN] Would upload: {filepath.name}")
        return True

    try:
        result = subprocess.run(
            [
                "gh",
                "release",
                "upload",
                tag,
                str(filepath),
                "--repo",
                REPO,
                "--clobber",
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            print(f"  [OK] Uploaded: {filepath.name}")
            return True
        else:
            print(f"  [FAIL] Upload failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"  [FAIL] Exception: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Upload planet_centers SPK files to GitHub Releases",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python scripts/release_planet_centers.py --version 0.15.0
    python scripts/release_planet_centers.py --version 0.15.0 --tier extended
    python scripts/release_planet_centers.py --version 0.15.0 --dry-run
        """,
    )

    parser.add_argument(
        "--version",
        required=True,
        help="Version string for logging purposes (e.g., 0.15.0)",
    )
    parser.add_argument(
        "--tier",
        choices=["base", "medium", "extended", "legacy", "all"],
        default="all",
        help="Which tier(s) to upload (default: all)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be uploaded without actually uploading",
    )
    parser.add_argument(
        "--tag",
        default=RELEASE_TAG,
        help=f"Release tag (default: {RELEASE_TAG})",
    )

    args = parser.parse_args()

    # Check authentication
    if not args.dry_run and not check_gh_auth():
        print("ERROR: gh CLI is not authenticated.")
        print("Run: gh auth login")
        return 1

    # Determine files to upload
    if args.tier == "all":
        tiers = DEFAULT_TIERS
    else:
        tiers = [args.tier]

    files_to_upload = []
    for tier in tiers:
        filename = TIER_FILES[tier]
        filepath = find_bsp_file(filename)
        if filepath is not None:
            files_to_upload.append((tier, filepath))
        else:
            print(f"[SKIP] {filename} not found in {get_data_dir()} or repo root")

    if not files_to_upload:
        print("ERROR: No files found to upload.")
        return 1

    # Print summary
    print("=" * 70)
    print("PLANET CENTERS RELEASE UPLOAD")
    print("=" * 70)
    print(f"Version: {args.version}")
    print(f"Release tag: {args.tag}")
    print(f"Repository: {REPO}")
    print(f"Dry run: {args.dry_run}")
    print()
    print("Files to upload:")
    print()

    for tier, filepath in files_to_upload:
        size_mb = get_file_size_mb(filepath)
        sha256 = calculate_sha256(filepath)
        print(f"  {filepath.name}")
        print(f"    Size:   {size_mb:.1f} MB")
        print(f"    SHA256: {sha256}")
        print()

    # Check release exists
    if not args.dry_run:
        if not check_release_exists(args.tag):
            print(f"WARNING: Release '{args.tag}' does not exist.")
            print(f"Create it first with: gh release create {args.tag} --repo {REPO}")
            return 1

    # Upload files
    print("Uploading...")
    print()

    success_count = 0
    for tier, filepath in files_to_upload:
        if upload_file(filepath, args.tag, args.dry_run):
            success_count += 1

    # Summary
    print()
    print("=" * 70)
    if args.dry_run:
        print(
            f"DRY RUN COMPLETE: {success_count}/{len(files_to_upload)} files would be uploaded"
        )
    else:
        print(f"UPLOAD COMPLETE: {success_count}/{len(files_to_upload)} files uploaded")
    print("=" * 70)

    # Print URLs
    if not args.dry_run:
        print()
        print("Download URLs:")
        base_url = f"https://github.com/{REPO}/releases/download/{args.tag}"
        for tier, filepath in files_to_upload:
            print(f"  {base_url}/{filepath.name}")

    # Update download.py with SHA256 hashes
    print()
    print("To update download.py with SHA256 hashes, add:")
    print()
    for tier, filepath in files_to_upload:
        sha256 = calculate_sha256(filepath)
        print(f'    "{filepath.name}": {{"sha256": "{sha256}", ...}},')

    return 0


if __name__ == "__main__":
    sys.exit(main())
