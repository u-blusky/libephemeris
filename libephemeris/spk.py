"""
SPK kernel support for high-precision minor body calculations.

This module provides functionality to:
- Download SPK (SPICE kernel) files from JPL Horizons API
- Register mappings between libephemeris body IDs and SPK targets
- Calculate positions using SPK data instead of Keplerian approximations

Using SPK kernels provides significantly higher precision for asteroids and TNOs:
- Keplerian model: ~10-30 arcseconds (asteroids), ~1-3 arcminutes (TNOs)
- SPK kernel: ~arcseconds to sub-arcsecond (within kernel coverage)

Usage:
    >>> import libephemeris as eph
    >>> # Download and register in one step
    >>> eph.download_and_register_spk(
    ...     body="Chiron",
    ...     ipl=eph.SE_CHIRON,
    ...     start="2000-01-01",
    ...     end="2100-01-01",
    ... )
    >>> # Now calc_ut uses SPK automatically
    >>> pos, _ = eph.calc_ut(2451545.0, eph.SE_CHIRON, eph.SEFLG_SPEED)

References:
    - JPL Horizons API: https://ssd-api.jpl.nasa.gov/doc/horizons.html
    - NAIF SPICE: https://naif.jpl.nasa.gov/naif/
"""

import json
import math
import os
import re
import time
import urllib.error
import urllib.request
from typing import Optional, Union

from skyfield.framelib import ecliptic_frame

from .exceptions import SPKNotFoundError
from .state import get_library_path, get_loader, get_timescale


# =============================================================================
# MODULE-LEVEL STATE (managed by state.py, but defined here for type hints)
# =============================================================================
# Actual state is in state.py to maintain single source of truth

# NAIF ID convention for numbered asteroids: naif_id = asteroid_number + 2000000
NAIF_ASTEROID_OFFSET = 2000000


# =============================================================================
# NAIF ID UTILITIES
# =============================================================================


def _extract_asteroid_number(body: str) -> Optional[int]:
    """
    Extract asteroid catalog number from body string.

    Args:
        body: Body identifier (e.g., "Chiron", "2060", "2060 Chiron", "(2060)")

    Returns:
        Asteroid number if found, None otherwise.

    Examples:
        >>> _extract_asteroid_number("2060")
        2060
        >>> _extract_asteroid_number("2060 Chiron")
        2060
        >>> _extract_asteroid_number("(136199) Eris")
        136199
        >>> _extract_asteroid_number("Chiron")
        None
    """
    # Try to find a number at the start or in parentheses
    # Pattern: optional parentheses, digits, optional name
    match = re.match(r"^\(?(\d+)\)?", body.strip())
    if match:
        return int(match.group(1))
    return None


def _deduce_naif_id(body: str, asteroid_number: Optional[int] = None) -> Optional[int]:
    """
    Deduce NAIF ID from body identifier.

    For numbered asteroids, NAIF ID = asteroid_number + 2000000.
    For unnumbered objects or named objects without numbers, returns None.

    Args:
        body: Body identifier string
        asteroid_number: Optional pre-extracted asteroid number

    Returns:
        NAIF ID if deducible, None otherwise.
    """
    if asteroid_number is None:
        asteroid_number = _extract_asteroid_number(body)

    if asteroid_number is not None:
        return asteroid_number + NAIF_ASTEROID_OFFSET

    return None


# =============================================================================
# BODY NAME UTILITIES
# =============================================================================

# Mapping from body IDs to human-readable names for helpful error messages
# This complements the SPK_BODY_NAME_MAP in constants.py which maps to Horizons IDs
_BODY_NAMES: dict[int, str] = {
    15: "Chiron",  # SE_CHIRON
    16: "Pholus",  # SE_PHOLUS
    17: "Ceres",  # SE_CERES
    18: "Pallas",  # SE_PALLAS
    19: "Juno",  # SE_JUNO
    20: "Vesta",  # SE_VESTA
    10000 + 136199: "Eris",  # SE_ERIS
    10000 + 90377: "Sedna",  # SE_SEDNA
    10000 + 136108: "Haumea",  # SE_HAUMEA
    10000 + 136472: "Makemake",  # SE_MAKEMAKE
    10000 + 28978: "Ixion",  # SE_IXION
    10000 + 90482: "Orcus",  # SE_ORCUS
    10000 + 50000: "Quaoar",  # SE_QUAOAR
    10000 + 20000: "Varuna",  # SE_VARUNA
    10000 + 7066: "Nessus",  # SE_NESSUS
    10000 + 8405: "Asbolus",  # SE_ASBOLUS
    10000 + 10199: "Chariklo",  # SE_CHARIKLO
    10000 + 225088: "Gonggong",  # SE_GONGGONG
    10000 + 99942: "Apophis",  # SE_APOPHIS
    10000 + 10: "Hygiea",  # SE_HYGIEA
    10000 + 433: "Eros",  # SE_EROS
}


def _get_body_name(ipl: int) -> Optional[str]:
    """
    Get the human-readable name for a body ID.

    Args:
        ipl: libephemeris body ID

    Returns:
        Human-readable name if known, None otherwise.
    """
    return _BODY_NAMES.get(ipl)


def _get_horizons_id_for_body(ipl: int) -> Optional[str]:
    """
    Get the JPL Horizons target identifier for a body ID.

    Args:
        ipl: libephemeris body ID

    Returns:
        Horizons ID if available, None otherwise.
    """
    from .constants import SPK_BODY_NAME_MAP

    if ipl in SPK_BODY_NAME_MAP:
        return SPK_BODY_NAME_MAP[ipl][0]
    return None


# =============================================================================
# HORIZONS API DOWNLOAD
# =============================================================================

# JPL Horizons API endpoint
HORIZONS_API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"


def _build_horizons_url(
    body: str,
    start: str,
    end: str,
    center: str = "500@0",
) -> str:
    """
    Build URL for JPL Horizons SPK file request.

    Args:
        body: Target body (name or number)
        start: Start date (YYYY-MM-DD)
        end: End date (YYYY-MM-DD)
        center: Reference center (default: 500@0 = Solar System Barycenter)

    Returns:
        URL string for Horizons API request.
    """
    # URL-encode the body name (handle spaces, special chars)
    # For numbered asteroids, use the number directly
    # For named objects, Horizons accepts the name
    params = {
        "format": "json",
        "COMMAND": f"'{body}'",
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "SPK",
        "CENTER": f"'{center}'",
        "START_TIME": f"'{start}'",
        "STOP_TIME": f"'{end}'",
    }

    query = "&".join(f"{k}={v}" for k, v in params.items())
    return f"{HORIZONS_API_URL}?{query}"


def _sanitize_filename(body: str) -> str:
    """
    Create a safe filename from body name.

    Args:
        body: Body name/identifier

    Returns:
        Sanitized string safe for use in filenames.
    """
    # Remove/replace unsafe characters
    safe = re.sub(r"[^\w\-]", "_", body.lower())
    # Collapse multiple underscores
    safe = re.sub(r"_+", "_", safe)
    # Remove leading/trailing underscores
    safe = safe.strip("_")
    return safe or "body"


def download_spk(
    body: str,
    start: str,
    end: str,
    path: Optional[str] = None,
    center: str = "500@0",
    overwrite: bool = False,
    timeout: int = 120,
) -> str:
    """
    Download SPK ephemeris file from JPL Horizons.

    Requests an SPK (SPICE kernel) file for the specified body from JPL Horizons
    API. The file contains high-precision ephemeris data that can be used for
    accurate position calculations.

    Args:
        body: Target body identifier. Can be:
            - Asteroid number: "2060", "136199"
            - Name: "Chiron", "Eris"
            - Combined: "2060 Chiron", "(136199) Eris"
        start: Start date in YYYY-MM-DD format
        end: End date in YYYY-MM-DD format
        path: Directory to save the file. If None, uses get_library_path()
        center: Reference center for ephemeris (default: "500@0" = SSB).
            Use "500@0" for compatibility with Skyfield/DE kernels.
        overwrite: If True, overwrite existing file. If False, skip if exists.
        timeout: Request timeout in seconds (default: 120)

    Returns:
        str: Full path to the downloaded SPK file

    Raises:
        ValueError: If body not found on Horizons or invalid parameters
        ConnectionError: If network request fails after retries
        IOError: If file cannot be written

    Example:
        >>> path = download_spk("Chiron", "2000-01-01", "2100-01-01")
        >>> print(path)
        /path/to/kernels/chiron_2000_2100.bsp
    """
    # Determine output directory
    if path is None:
        path = get_library_path()

    # Ensure directory exists
    os.makedirs(path, exist_ok=True)

    # Generate filename
    body_safe = _sanitize_filename(body)
    start_short = start.replace("-", "")[:6]  # YYYYMM
    end_short = end.replace("-", "")[:6]
    filename = f"{body_safe}_{start_short}_{end_short}.bsp"
    filepath = os.path.join(path, filename)

    # Check if already exists
    if os.path.exists(filepath) and not overwrite:
        return filepath

    # Build request URL
    url = _build_horizons_url(body, start, end, center)

    # Make request with retry
    max_retries = 2
    last_error = None

    for attempt in range(max_retries + 1):
        try:
            # First request: get SPK file URL from Horizons
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "libephemeris/0.1")

            with urllib.request.urlopen(req, timeout=timeout) as response:
                data = json.loads(response.read().decode("utf-8"))

            # Check for errors in response
            if "error" in data:
                raise ValueError(f"Horizons API error: {data['error']}")

            # Get SPK file URL from response
            if "spk_file_id" not in data or "spk" not in data:
                # Check if it's a "result" response (body not found, etc.)
                if "result" in data:
                    raise ValueError(f"Horizons error: {data['result'][:500]}")
                raise ValueError(
                    f"Unexpected Horizons response format. Keys: {list(data.keys())}"
                )

            spk_url = data["spk"]

            # Download the actual SPK file
            spk_req = urllib.request.Request(spk_url)
            spk_req.add_header("User-Agent", "libephemeris/0.1")

            with urllib.request.urlopen(spk_req, timeout=timeout) as spk_response:
                spk_data = spk_response.read()

            # Write to file
            with open(filepath, "wb") as f:
                f.write(spk_data)

            return filepath

        except urllib.error.HTTPError as e:
            last_error = e
            if e.code == 400:
                # Bad request - likely invalid body name
                try:
                    error_data = json.loads(e.read().decode("utf-8"))
                    msg = error_data.get("message", str(e))
                except Exception:
                    msg = str(e)
                raise ValueError(f"Invalid body '{body}': {msg}") from e
            elif attempt < max_retries:
                time.sleep(1)  # Brief pause before retry
                continue
            else:
                raise ConnectionError(
                    f"Failed to download SPK after {max_retries + 1} attempts: {e}"
                ) from e

        except urllib.error.URLError as e:
            last_error = e
            if attempt < max_retries:
                time.sleep(1)
                continue
            else:
                raise ConnectionError(
                    f"Network error downloading SPK: {e.reason}"
                ) from e

    # Should not reach here, but just in case
    raise ConnectionError(f"Download failed: {last_error}")


# =============================================================================
# SPK REGISTRATION
# =============================================================================


def register_spk_body(
    ipl: int,
    spk_file: str,
    naif_id: Union[int, str],
) -> None:
    """
    Register a mapping between a libephemeris body ID and an SPK kernel target.

    After registration, swe_calc_ut() and swe_calc() will automatically use the
    SPK kernel for this body instead of the Keplerian approximation.

    Args:
        ipl: libephemeris body ID (e.g., SE_CHIRON, SE_ERIS)
        spk_file: Path to the SPK file, or filename if in library path
        naif_id: NAIF ID of the body in the SPK kernel.
            For numbered asteroids: asteroid_number + 2000000
            (e.g., Chiron = 2060 -> naif_id = 2002060)

    Raises:
        SPKNotFoundError: If SPK file not found (with helpful instructions)
        ValueError: If naif_id not found in SPK kernel

    Example:
        >>> register_spk_body(SE_CHIRON, "/path/to/chiron.bsp", 2002060)
        >>> # Now calc_ut(jd, SE_CHIRON, ...) uses SPK data
    """
    from . import state

    # Convert naif_id to int if string
    if isinstance(naif_id, str):
        naif_id = int(naif_id)

    # Resolve file path
    if not os.path.isabs(spk_file):
        # Try library path
        lib_path = get_library_path()
        full_path = os.path.join(lib_path, spk_file)
        if os.path.exists(full_path):
            spk_file = full_path
        elif not os.path.exists(spk_file):
            # Get helpful info for error message
            body_name = _get_body_name(ipl)
            body_id = _get_horizons_id_for_body(ipl)
            raise SPKNotFoundError.from_filepath(
                filepath=spk_file,
                body_name=body_name,
                body_id=body_id,
            )

    if not os.path.exists(spk_file):
        # Get helpful info for error message
        body_name = _get_body_name(ipl)
        body_id = _get_horizons_id_for_body(ipl)
        raise SPKNotFoundError.from_filepath(
            filepath=spk_file,
            body_name=body_name,
            body_id=body_id,
        )

    # Load kernel if not already cached
    state._load_spk_kernel(spk_file)

    # Validate that naif_id exists in kernel
    kernel = state._SPK_KERNELS.get(spk_file)
    if kernel is not None:
        # Check if target exists in kernel
        try:
            _ = kernel[naif_id]
        except KeyError:
            # Try with string name
            try:
                _ = kernel[str(naif_id)]
            except KeyError:
                available = []
                if hasattr(kernel, "names"):
                    available = list(kernel.names())[:10]
                raise ValueError(
                    f"NAIF ID {naif_id} not found in SPK kernel {spk_file}. "
                    f"Available targets (first 10): {available}"
                )

    # Register mapping
    state._SPK_BODY_MAP[ipl] = (spk_file, naif_id)


def unregister_spk_body(ipl: int) -> None:
    """
    Remove SPK registration for a body.

    After unregistration, the body will use Keplerian approximation again.

    Args:
        ipl: libephemeris body ID (e.g., SE_CHIRON)
    """
    from . import state

    if ipl in state._SPK_BODY_MAP:
        del state._SPK_BODY_MAP[ipl]


def get_spk_body_info(ipl: int) -> Optional[tuple[str, int]]:
    """
    Get SPK registration info for a body.

    Args:
        ipl: libephemeris body ID

    Returns:
        Tuple of (spk_file, naif_id) if registered, None otherwise.
    """
    from . import state

    return state._SPK_BODY_MAP.get(ipl)


def list_spk_bodies() -> dict[int, tuple[str, int]]:
    """
    List all registered SPK body mappings.

    Returns:
        Dict mapping ipl -> (spk_file, naif_id) for all registered bodies.
    """
    from . import state

    return dict(state._SPK_BODY_MAP)


def get_spk_coverage(spk_file: str) -> Optional[tuple[float, float]]:
    """
    Get the time coverage of an SPK kernel.

    Args:
        spk_file: Path to SPK file

    Returns:
        Tuple of (start_jd, end_jd) Julian Day coverage, or None if error.
    """
    from . import state

    # Load if not cached
    if spk_file not in state._SPK_KERNELS:
        if not os.path.isabs(spk_file):
            full_path = os.path.join(get_library_path(), spk_file)
            if os.path.exists(full_path):
                spk_file = full_path

        if not os.path.exists(spk_file):
            return None

        state._load_spk_kernel(spk_file)

    kernel = state._SPK_KERNELS.get(spk_file)
    if kernel is None:
        return None

    try:
        if hasattr(kernel, "spk") and hasattr(kernel.spk, "segments"):
            segments = list(kernel.spk.segments)
            if segments:
                start_jd = min(
                    float(s.start_jd) for s in segments if hasattr(s, "start_jd")
                )
                end_jd = max(float(s.end_jd) for s in segments if hasattr(s, "end_jd"))
                return (start_jd, end_jd)
    except Exception:
        pass

    return None


# =============================================================================
# CONVENIENCE FUNCTION
# =============================================================================


def download_and_register_spk(
    body: str,
    ipl: int,
    start: str,
    end: str,
    naif_id: Optional[int] = None,
    path: Optional[str] = None,
    center: str = "500@0",
    overwrite: bool = False,
    timeout: int = 120,
) -> str:
    """
    Download SPK from Horizons and register it for use in calculations.

    Convenience function that combines download_spk() and register_spk_body().

    Args:
        body: Target body identifier (see download_spk)
        ipl: libephemeris body ID (e.g., SE_CHIRON)
        start: Start date (YYYY-MM-DD)
        end: End date (YYYY-MM-DD)
        naif_id: NAIF ID in the kernel. If None, deduced from body number
            using convention: naif_id = asteroid_number + 2000000
        path: Directory to save file (default: get_library_path())
        center: Reference center (default: "500@0" = SSB)
        overwrite: Overwrite existing file if True
        timeout: Request timeout in seconds

    Returns:
        str: Path to downloaded SPK file

    Raises:
        ValueError: If body not found, or naif_id cannot be deduced and not provided

    Example:
        >>> download_and_register_spk(
        ...     body="Chiron",
        ...     ipl=SE_CHIRON,
        ...     start="2000-01-01",
        ...     end="2100-01-01",
        ... )
        >>> # Chiron now uses SPK for calculations
    """
    # Download
    spk_path = download_spk(
        body=body,
        start=start,
        end=end,
        path=path,
        center=center,
        overwrite=overwrite,
        timeout=timeout,
    )

    # Deduce NAIF ID if not provided
    if naif_id is None:
        naif_id = _deduce_naif_id(body)
        if naif_id is None:
            raise ValueError(
                f"Cannot deduce NAIF ID for '{body}'. "
                f"Please provide naif_id explicitly. "
                f"For numbered asteroids: naif_id = asteroid_number + 2000000"
            )

    # Register
    register_spk_body(ipl, spk_path, naif_id)

    return spk_path


# =============================================================================
# SPK POSITION CALCULATION
# =============================================================================


def calc_spk_body_position(
    t,
    ipl: int,
    iflag: int,
) -> Optional[tuple[float, float, float, float, float, float]]:
    """
    Calculate body position using SPK kernel.

    Internal function called by planets._calc_body() when SPK is registered.

    Args:
        t: Skyfield Time object
        ipl: libephemeris body ID
        iflag: Calculation flags (SEFLG_SPEED, SEFLG_HELCTR, etc.)

    Returns:
        Position tuple (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        or None if body not registered or outside coverage.

    Raises:
        ValueError: If JD is outside SPK coverage
    """
    from . import state
    from .constants import SEFLG_HELCTR, SEFLG_SPEED, SEFLG_SIDEREAL
    from .planets import swe_get_ayanamsa_ut

    # Check if body is registered
    if ipl not in state._SPK_BODY_MAP:
        return None

    spk_file, naif_id = state._SPK_BODY_MAP[ipl]

    # Get kernel
    kernel = state._SPK_KERNELS.get(spk_file)
    if kernel is None:
        return None

    # Check coverage
    coverage = get_spk_coverage(spk_file)
    if coverage is not None:
        start_jd, end_jd = coverage
        if t.tt < start_jd or t.tt > end_jd:
            raise ValueError(
                f"JD {t.tt:.1f} outside SPK coverage [{start_jd:.1f}, {end_jd:.1f}] "
                f"for body {ipl}"
            )

    # Get target from kernel
    try:
        target = kernel[naif_id]
    except KeyError:
        return None

    # Get main ephemeris for Earth/Sun
    planets = state.get_planets()

    # Determine observer
    is_heliocentric = bool(iflag & SEFLG_HELCTR)

    if is_heliocentric:
        observer = planets["sun"]
    else:
        observer = planets["earth"]

    # Calculate position
    # Target position relative to SSB
    target_pos = target.at(t)

    # Observer position relative to SSB
    observer_pos = observer.at(t)

    # Relative position (target as seen from observer)
    # We need to compute this manually since SPK targets may not support observe()
    target_xyz = target_pos.position.au
    observer_xyz = observer_pos.position.au

    rel_xyz = [target_xyz[i] - observer_xyz[i] for i in range(3)]

    # Convert to ecliptic coordinates
    # Use Skyfield's ecliptic frame
    from skyfield.positionlib import ICRF

    rel_pos = ICRF(rel_xyz, t=t, center=399)

    # Get ecliptic lat/lon
    ecl_pos = rel_pos.frame_latlon(ecliptic_frame)
    lat = ecl_pos[0].degrees
    lon = ecl_pos[1].degrees
    dist = ecl_pos[2].au

    # Normalize longitude to 0-360
    lon = lon % 360.0

    # Calculate speeds if requested
    speed_lon, speed_lat, speed_dist = 0.0, 0.0, 0.0

    if iflag & SEFLG_SPEED:
        # Numerical differentiation (1 second timestep)
        ts = get_timescale()
        dt = 1.0 / 86400.0  # 1 second in days

        t_next = ts.tt_jd(t.tt + dt)

        # Position at t + dt
        target_pos_next = target.at(t_next)
        observer_pos_next = observer.at(t_next)

        target_xyz_next = target_pos_next.position.au
        observer_xyz_next = observer_pos_next.position.au

        rel_xyz_next = [target_xyz_next[i] - observer_xyz_next[i] for i in range(3)]

        rel_pos_next = ICRF(rel_xyz_next, t=t_next, center=399)
        ecl_pos_next = rel_pos_next.frame_latlon(ecliptic_frame)

        lat_next = ecl_pos_next[0].degrees
        lon_next = ecl_pos_next[1].degrees
        dist_next = ecl_pos_next[2].au

        lon_next = lon_next % 360.0

        # Compute rates (per day)
        speed_lon = (lon_next - lon) / dt
        speed_lat = (lat_next - lat) / dt
        speed_dist = (dist_next - dist) / dt

        # Handle 360 wrap
        if speed_lon > 180.0 / dt:
            speed_lon -= 360.0 / dt
        if speed_lon < -180.0 / dt:
            speed_lon += 360.0 / dt

    # Apply sidereal correction if requested
    if iflag & SEFLG_SIDEREAL:
        ayanamsa = swe_get_ayanamsa_ut(t.ut1)
        lon = (lon - ayanamsa) % 360.0

    return (lon, lat, dist, speed_lon, speed_lat, speed_dist)
