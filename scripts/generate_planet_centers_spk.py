#!/usr/bin/env python3
"""
Generate a compact SPK file containing only planet center segments.

This script downloads the full satellite SPK files from JPL NAIF, extracts
only the planet center segments (599, 699, 799, 899, 999), and creates a
single compact SPK file suitable for distribution with libephemeris.

The resulting file is much smaller than the original satellite SPKs
(~5-10 MB vs ~500 MB) because it excludes all satellite data.

BACKGROUND
----------
JPL's DE ephemerides (DE421, DE440, etc.) provide only system barycenters
for the outer planets (Jupiter, Saturn, Uranus, Neptune, Pluto). To get
the actual planet center positions, we need the satellite SPK files which
contain the offset from barycenter to planet center.

NAIF ID Codes:
    5 = Jupiter Barycenter    599 = Jupiter Center
    6 = Saturn Barycenter     699 = Saturn Center
    7 = Uranus Barycenter     799 = Uranus Center
    8 = Neptune Barycenter    899 = Neptune Center
    9 = Pluto Barycenter      999 = Pluto Center

USAGE
-----
    # Generate the compact SPK (requires ~500MB temporary download)
    python scripts/generate_planet_centers_spk.py

    # Or using poe task
    poe generate-planet-centers-spk

OUTPUT
------
    libephemeris/data/planet_centers.bsp

REQUIREMENTS
------------
    - spiceypy >= 6.0.0 (pip install spiceypy)
    - Internet connection (to download source SPK files)
    - ~500 MB temporary disk space

REFERENCES
----------
    - NAIF Generic Kernels: https://naif.jpl.nasa.gov/naif/data_generic.html
    - SPK File Format: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
    - Skyfield SPK Support: https://rhodesmill.org/skyfield/planets.html

Author: libephemeris team
Date: 2024
"""

import os
import ssl
import struct
import sys
import tempfile
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

# =============================================================================
# CONFIGURATION
# =============================================================================

# Source SPK files from JPL NAIF (older, smaller versions)
# These contain planet centers relative to system barycenters
# Format: (url, barycenter_id, center_id)
SPK_SOURCES: Dict[str, Tuple[str, int, int]] = {
    "jupiter": (
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/jup204.bsp",
        5,
        599,
    ),
    "saturn": (
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/sat319.bsp",
        6,
        699,
    ),
    "uranus": (
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/ura083.bsp",
        7,
        799,
    ),
    "neptune": (
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/nep050.bsp",
        8,
        899,
    ),
    "pluto": (
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/plu017.bsp",
        9,
        999,
    ),
}

# Output file path (relative to repository root)
OUTPUT_FILE = "libephemeris/data/planet_centers.bsp"

# Leap seconds kernel URL (required by SPICE for time conversions)
LEAPSECONDS_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"


# =============================================================================
# SSL CONTEXT (for downloads)
# =============================================================================


def _get_ssl_context(verify: bool = True) -> ssl.SSLContext:
    """
    Create an SSL context for HTTPS downloads.

    Some systems have SSL certificate issues with JPL servers.
    This creates either a verified or permissive context.

    Returns:
        SSL context for urllib
    """
    ctx = ssl.create_default_context()
    if not verify:
        # Permissive mode for JPL servers (trusted source)
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
    return ctx


# =============================================================================
# DOWNLOAD FUNCTIONS
# =============================================================================


def download_file(
    url: str,
    dest_path: str,
    ssl_ctx: ssl.SSLContext,
    fallback_ctx: Optional[ssl.SSLContext] = None,
) -> str:
    """
    Download a file from URL to destination path with progress reporting.

    Args:
        url: Source URL
        dest_path: Destination file path
        ssl_ctx: SSL context for HTTPS

    Returns:
        Path to downloaded file

    Raises:
        urllib.error.URLError: If download fails
    """
    filename = os.path.basename(url)
    print(f"  Downloading {filename}...")

    # Get file size
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=30, context=ssl_ctx) as response:
            total_size = int(response.headers.get("Content-Length", 0))
    except Exception:
        if fallback_ctx is None:
            total_size = 0
        else:
            try:
                req = urllib.request.Request(url, method="HEAD")
                with urllib.request.urlopen(
                    req, timeout=30, context=fallback_ctx
                ) as response:
                    total_size = int(response.headers.get("Content-Length", 0))
                print("    Warning: SSL verification failed; using permissive mode")
            except Exception:
                total_size = 0

    # Download with progress
    downloaded = 0
    chunk_size = 1024 * 1024  # 1 MB chunks

    try:
        with urllib.request.urlopen(url, timeout=300, context=ssl_ctx) as response:
            with open(dest_path, "wb") as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        pct = downloaded / total_size * 100
                        print(
                            f"\r    {pct:5.1f}% ({downloaded / 1024 / 1024:.1f} MB)",
                            end="",
                        )
                    else:
                        print(f"\r    {downloaded / 1024 / 1024:.1f} MB", end="")
    except Exception:
        if fallback_ctx is None:
            raise
        print("    Warning: SSL verification failed; using permissive mode")
        with urllib.request.urlopen(url, timeout=300, context=fallback_ctx) as response:
            with open(dest_path, "wb") as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        pct = downloaded / total_size * 100
                        print(
                            f"\r    {pct:5.1f}% ({downloaded / 1024 / 1024:.1f} MB)",
                            end="",
                        )
                    else:
                        print(f"\r    {downloaded / 1024 / 1024:.1f} MB", end="")

    print(f"\r    100.0% ({downloaded / 1024 / 1024:.1f} MB) - Done")
    return dest_path


def download_source_spks(
    cache_dir: str,
    ssl_ctx: ssl.SSLContext,
    fallback_ctx: Optional[ssl.SSLContext] = None,
) -> Dict[str, str]:
    """
    Download all source SPK files to cache directory.

    Args:
        cache_dir: Directory to store downloaded files
        ssl_ctx: SSL context for HTTPS

    Returns:
        Dictionary mapping planet name to local file path
    """
    print("\n=== Step 1: Downloading Source SPK Files ===\n")

    local_files = {}
    for planet, (url, _, _) in SPK_SOURCES.items():
        filename = os.path.basename(url)
        local_path = os.path.join(cache_dir, filename)

        if os.path.exists(local_path):
            print(f"  {filename}: Already cached")
            local_files[planet] = local_path
        else:
            download_file(url, local_path, ssl_ctx, fallback_ctx)
            local_files[planet] = local_path

    return local_files


def download_leapseconds(
    cache_dir: str,
    ssl_ctx: ssl.SSLContext,
    fallback_ctx: Optional[ssl.SSLContext] = None,
) -> str:
    """
    Download leap seconds kernel (required by SPICE).

    Args:
        cache_dir: Directory to store downloaded file
        ssl_ctx: SSL context for HTTPS

    Returns:
        Path to leap seconds kernel
    """
    filename = os.path.basename(LEAPSECONDS_URL)
    local_path = os.path.join(cache_dir, filename)

    if not os.path.exists(local_path):
        print(f"  Downloading leap seconds kernel ({filename})...")
        download_file(LEAPSECONDS_URL, local_path, ssl_ctx, fallback_ctx)

    return local_path


# =============================================================================
# SPK EXTRACTION FUNCTIONS
# =============================================================================


def pack_spk_descriptor(dc: np.ndarray, ic: np.ndarray) -> np.ndarray:
    """
    Pack double and integer components into a 5-element SPK descriptor.

    The spiceypy dafgs() function returns an oversized array, but spksub()
    requires exactly 5 elements. This function creates the correct format.

    SPK descriptors have ND=2 (doubles) and NI=6 (integers), which pack into
    2 + ceil(6/2) = 5 doubles.

    Args:
        dc: Array of 2 doubles (start_et, end_et)
        ic: Array of 6 integers (target, center, frame, type, start_addr, end_addr)

    Returns:
        5-element numpy array containing packed descriptor
    """
    descr = np.zeros(5, dtype=np.float64)
    descr[0] = dc[0]
    descr[1] = dc[1]
    # Pack 6 integers into 3 doubles (each double holds 2 32-bit ints)
    for i in range(3):
        i1 = int(ic[2 * i])
        i2 = int(ic[2 * i + 1])
        packed = struct.pack("<ii", i1, i2)
        descr[2 + i] = struct.unpack("<d", packed)[0]
    return descr


def extract_segment_with_spksub(
    source_file: str,
    output_handle: int,
    target_id: int,
    planet_name: str,
) -> Optional[Tuple[float, float]]:
    """
    Extract a segment from source SPK and write to output using spksub.

    spksub extracts a subset of an SPK segment, preserving the original
    data type and precision.

    Args:
        source_file: Path to source SPK file
        output_handle: Handle to open output SPK file
        target_id: NAIF ID of body to extract (e.g., 599)
        planet_name: Name for logging

    Returns:
        Tuple of (start_et, end_et) or None if failed
    """
    import spiceypy as spice

    print(f"  Extracting {planet_name} center ({target_id})...")

    # Open source file for reading
    source_handle = spice.spklef(source_file)

    try:
        # Get coverage for the target body
        cover = spice.spkcov(source_file, target_id)
        n_windows = spice.wncard(cover)

        if n_windows == 0:
            print(f"    ERROR: No coverage for {target_id} in {source_file}")
            return None

        # Get the coverage window
        start_et, end_et = spice.wnfetd(cover, 0)
        print(f"    Coverage: {start_et:.0f} to {end_et:.0f} ET")

        # Find the segment containing this body
        # We need to search the DAF for the segment descriptor
        spice.dafbfs(source_handle)

        found_segment = False
        while True:
            found = spice.daffna()
            if not found:
                break

            # Get segment summary (descriptor) - spiceypy returns oversized array
            raw_descr = spice.dafgs()

            # Unpack: dc = [start_et, end_et], ic = [target, center, frame, type, start, end]
            dc, ic = spice.dafus(raw_descr, 2, 6)

            if int(ic[0]) == target_id:
                found_segment = True
                seg_start, seg_end = dc[0], dc[1]
                seg_center = int(ic[1])
                seg_type = int(ic[3])

                # Get segment identifier
                ident = spice.dafgn()

                print(f"    Segment: center={seg_center}, type={seg_type}")
                print(f"    Identifier: {ident}")

                # Pack descriptor into 5-element format required by spksub
                descr5 = pack_spk_descriptor(dc, ic)

                # Use spksub to copy the segment to output
                spice.spksub(
                    source_handle,  # Source SPK handle
                    descr5,  # 5-element packed descriptor
                    ident,  # Segment identifier
                    seg_start,  # Begin time
                    seg_end,  # End time
                    output_handle,  # Destination handle
                )

                print(f"    Written: {ident}")
                break

        if not found_segment:
            print(f"    ERROR: Segment for {target_id} not found")
            return None

        return (start_et, end_et)

    finally:
        # Unload source file
        spice.spkuef(source_handle)


def generate_compact_spk(
    source_files: Dict[str, str],
    output_file: str,
    leapseconds_file: str,
) -> str:
    """
    Generate compact SPK containing only planet center segments.

    Args:
        source_files: Dictionary mapping planet name to source SPK path
        output_file: Path for output SPK file
        leapseconds_file: Path to leap seconds kernel

    Returns:
        Path to generated SPK file
    """
    import spiceypy as spice

    print("\n=== Step 2: Extracting Planet Center Segments ===\n")

    # Load leap seconds kernel (required for time conversions)
    spice.furnsh(leapseconds_file)

    # Create output SPK file
    # Internal filename (max 60 chars)
    ifname = "Planet Centers SPK for libephemeris"

    # Number of characters for comments (we won't add any via SPICE)
    ncomch = 0

    output_handle = spice.spkopn(output_file, ifname, ncomch)

    # Track results
    results = {}

    # Extract each planet center
    for planet, (url, center_id, target_id) in SPK_SOURCES.items():
        source_file = source_files[planet]

        try:
            result = extract_segment_with_spksub(
                source_file,
                output_handle,
                target_id,
                planet,
            )
            results[planet] = result
        except Exception as e:
            print(f"    ERROR: {e}")
            results[planet] = None

    # Close the output file
    spice.spkcls(output_handle)

    # Clear SPICE state
    spice.kclear()

    return output_file


# =============================================================================
# MAIN FUNCTION
# =============================================================================


def main():
    """
    Main entry point for SPK generation.

    This function orchestrates the full workflow:
    1. Download source SPK files from JPL NAIF
    2. Download leap seconds kernel
    3. Extract planet center segments using spksub
    4. Generate compact output SPK
    5. Verify the result
    """
    print("=" * 70)
    print("PLANET CENTERS SPK GENERATOR")
    print("=" * 70)
    print()
    print("This script generates a compact SPK file containing only planet")
    print("center positions for Jupiter, Saturn, Uranus, Neptune, and Pluto.")
    print()

    # Check for spiceypy
    try:
        import spiceypy

        print(f"Using spiceypy version {spiceypy.__version__}")
    except ImportError:
        print("ERROR: spiceypy is required but not installed.")
        print("Install with: pip install spiceypy>=6.0.0")
        return 1

    # Determine paths
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent
    output_path = repo_root / OUTPUT_FILE

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create SSL contexts for downloads (verified first, permissive fallback)
    ssl_ctx = _get_ssl_context(verify=True)
    fallback_ctx = _get_ssl_context(verify=False)

    # Use a temporary directory for source files
    # These are large and we don't need to keep them
    with tempfile.TemporaryDirectory(prefix="planet_centers_spk_") as cache_dir:
        print(f"Cache directory: {cache_dir}")

        # Download source SPK files
        source_files = download_source_spks(cache_dir, ssl_ctx, fallback_ctx)

        # Download leap seconds kernel
        leapseconds = download_leapseconds(cache_dir, ssl_ctx, fallback_ctx)

        # Generate compact SPK
        result = generate_compact_spk(source_files, str(output_path), leapseconds)

        print("\n=== Step 3: Verification ===\n")

        # Verify the output
        file_size = os.path.getsize(result)
        print(f"  Output file: {result}")
        print(f"  File size: {file_size / 1024 / 1024:.2f} MB")

        # Verify with spiceypy
        import spiceypy as spice

        # Get list of bodies in the output file
        bodies = spice.spkobj(str(output_path))
        print("\n  Bodies in output file:")
        for body_id in bodies:
            print(f"    {body_id}")

        # Quick verification with Skyfield
        try:
            from skyfield.api import Loader

            loader = Loader(cache_dir)
            spk = loader(str(output_path))

            print("\n  Segments (via Skyfield):")
            for seg in spk.segments:
                print(
                    f"    {seg.center} -> {seg.target}: {seg.center_name} -> {seg.target_name}"
                )

        except Exception as e:
            print(f"\n  Warning: Could not verify with Skyfield: {e}")

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"\nGenerated: {output_path}")
    print(f"Size: {file_size / 1024 / 1024:.2f} MB")
    print("\nTo use with libephemeris, this file will be loaded automatically")
    print("alongside DE421/DE440 to provide planet center positions.")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
