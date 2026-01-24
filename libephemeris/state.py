"""
Global state management for libephemeris.

This module maintains the library's singleton state including:
- Ephemeris data loader (Skyfield Loader)
- Planetary ephemeris (DE421 or other JPL files)
- Timescale (for UTC/TT conversions)
- Observer topocentric location
- Sidereal mode configuration
- Cached angles for Arabic parts calculations

All state is stored in module-level globals to provide SwissEphemeris-compatible
stateful API behavior. This is thread-unsafe by design, matching SwissEph behavior.
"""

import os
from typing import Optional, Union
from skyfield.api import Loader, Topos
from skyfield.timelib import Timescale
from skyfield.jpllib import SpiceKernel

# =============================================================================
# GLOBAL STATE VARIABLES
# =============================================================================

_EPHEMERIS_PATH: Optional[str] = None  # Custom ephemeris directory
_EPHEMERIS_FILE: str = "de421.bsp"  # Ephemeris file to use (default: DE421)
_LOADER: Optional[Loader] = None  # Skyfield data loader
_PLANETS: Optional[SpiceKernel] = None  # Loaded planetary ephemeris
_TS: Optional[Timescale] = None  # Timescale object
_TOPO: Optional[Topos] = None  # Observer location
_SIDEREAL_MODE: Optional[int] = None  # Active sidereal mode ID
_SIDEREAL_AYAN_T0: float = 0.0  # Ayanamsha value at reference epoch
_SIDEREAL_T0: float = 0.0  # Reference epoch (JD) for ayanamsha
_ANGLES_CACHE: dict[str, float] = {}  # Pre-calculated angles {name: longitude}
_TIDAL_ACCELERATION: Optional[float] = None  # Tidal acceleration for Delta T
_DELTA_T_USERDEF: Optional[float] = None  # User-defined Delta T value


def get_loader() -> Loader:
    """
    Get or create the Skyfield data loader.

    Returns:
        Loader: Skyfield Loader instance for downloading/caching ephemeris files

    Note:
        Data files are cached in the parent directory of this module by default.
    """
    global _LOADER
    if _LOADER is None:
        data_dir = os.path.join(os.path.dirname(__file__), "..")
        _LOADER = Loader(data_dir)
    return _LOADER


def get_timescale() -> Timescale:
    """
    Get or create the Skyfield timescale object.

    Returns:
        Timescale: Skyfield timescale for time conversions (UTC, TT, etc.)

    Note:
        Automatically downloads IERS data for Delta T calculations if needed.
    """
    global _TS
    if _TS is None:
        load = get_loader()
        _TS = load.timescale()
    return _TS


def get_planets() -> SpiceKernel:
    """
    Get or load the planetary ephemeris (DE421 by default).

    Returns:
        SpiceKernel: Loaded JPL ephemeris kernel containing planetary positions

    Raises:
        FileNotFoundError: If ephemeris file cannot be found or downloaded

    Note:
        Uses the ephemeris file set via set_ephemeris_file() (default: de421.bsp).
        Searches in _EPHEMERIS_PATH if set, then workspace root, then downloads.
    """
    global _PLANETS
    if _PLANETS is None:
        load = get_loader()

        # Try custom ephemeris path first if set
        if _EPHEMERIS_PATH:
            bsp_path = os.path.join(_EPHEMERIS_PATH, _EPHEMERIS_FILE)
            if os.path.exists(bsp_path):
                _PLANETS = load(bsp_path)
                return _PLANETS

        # Try workspace root
        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        bsp_path = os.path.join(base_dir, _EPHEMERIS_FILE)
        if os.path.exists(bsp_path):
            _PLANETS = load(bsp_path)
        else:
            # Download from internet
            _PLANETS = load(_EPHEMERIS_FILE)
    return _PLANETS


def set_topo(lon: float, lat: float, alt: float) -> None:
    """
    Set observer's topocentric location for planet calculations.

    Args:
        lon: Geographic longitude in degrees (East positive)
        lat: Geographic latitude in degrees (North positive)
        alt: Elevation above sea level in meters

    Note:
        Required for topocentric calculations (SEFLG_TOPOCTR),
        angles (Ascendant, MC), and Arabic parts.
    """
    global _TOPO
    _TOPO = Topos(latitude_degrees=lat, longitude_degrees=lon, elevation_m=alt)


def get_topo() -> Optional[Topos]:
    """
    Get the observer's topocentric location.

    Returns:
        Optional[Topos]: Current observer location or None if not set
    """
    return _TOPO


def set_sid_mode(mode: int, t0: float = 0.0, ayan_t0: float = 0.0) -> None:
    """
    Set the sidereal mode (ayanamsha system) for calculations.

    Args:
        mode: Sidereal mode ID (SE_SIDM_*) or 255 for custom
        t0: Reference epoch (Julian Day) for custom ayanamsha (default: J2000.0)
        ayan_t0: Ayanamsha value at t0 in degrees (for custom mode)

    Note:
        Affects all position calculations when SEFLG_SIDEREAL is set.
        Default is Lahiri (SE_SIDM_LAHIRI) if never set.
    """
    global _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0
    _SIDEREAL_MODE = mode
    _SIDEREAL_T0 = t0 if t0 != 0.0 else 2451545.0
    _SIDEREAL_AYAN_T0 = ayan_t0


def get_sid_mode(full: bool = False) -> Union[int, tuple[int, float, float]]:
    """
    Get the current sidereal mode configuration.

    Args:
        full: If True, return (mode, t0, ayan_t0); if False, return only mode ID

    Returns:
        int or tuple: Sidereal mode ID, or full configuration tuple

    Note:
        Returns SE_SIDM_LAHIRI (1) by default if never set.
    """
    if full:
        return _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0
    return _SIDEREAL_MODE if _SIDEREAL_MODE is not None else 1


def get_angles_cache() -> dict[str, float]:
    """
    Get cached astrological angles for the current calculation context.

    Returns:
        dict: Cached angles {name: longitude_degrees}

    Note:
        Used by Arabic parts calculations which require pre-calculated
        planetary positions and angles.
    """
    return _ANGLES_CACHE


def set_angles_cache(angles: dict[str, float]) -> None:
    """
    Cache pre-calculated angles for use in Arabic parts and other virtual points.

    Args:
        angles: Dictionary of angles {name: longitude_degrees}
                e.g., {"Sun": 120.5, "Moon": 240.3, "Asc": 15.7}

    Note:
        Creates a copy to prevent external mutation of cache.
    """
    global _ANGLES_CACHE
    _ANGLES_CACHE = angles.copy()


def clear_angles_cache() -> None:
    """
    Clear the angles cache.

    Use this between unrelated calculation contexts to prevent stale data.
    """
    global _ANGLES_CACHE
    _ANGLES_CACHE = {}

    _ANGLES_CACHE = {}


def set_ephe_path(path: Optional[str]) -> None:
    """
    Set the path for ephemeris files.

    Args:
        path: Path to directory containing ephemeris files.

    Note:
        This sets the directory where get_planets() will look for the ephemeris
        file specified by set_ephemeris_file(). If the file is not found there,
        it will fall back to the workspace root and then download if needed.
    """
    global _EPHEMERIS_PATH, _PLANETS
    _EPHEMERIS_PATH = path
    # Clear cached planets to force reload from new path
    _PLANETS = None


def set_ephemeris_file(filename: str) -> None:
    """
    Set the ephemeris file to use for planetary calculations.

    Args:
        filename: Name of the JPL ephemeris file (e.g., "de421.bsp", "de422.bsp", "de431.bsp")

    Note:
        This allows using different ephemeris files with varying date ranges:
        - de421.bsp: 1900-2050 (default, 16 MB)
        - de422.bsp: -3000-3000 (623 MB)
        - de430.bsp: 1550-2650 (128 MB)
        - de431.bsp: -13200-17191 (3.4 GB)

        The file will be searched in:
        1. The path set by set_ephe_path() if configured
        2. The workspace root directory
        3. Downloaded from JPL if not found locally

        Changing the ephemeris file clears the cached planets and forces a reload.
    """
    global _EPHEMERIS_FILE, _PLANETS
    _EPHEMERIS_FILE = filename
    # Clear cached planets to force reload with new file
    _PLANETS = None


def set_jpl_file(filename: str) -> None:
    """
    Set the JPL ephemeris file to use for planetary calculations.

    This is an alias for set_ephemeris_file() with a more descriptive name
    emphasizing JPL (Jet Propulsion Laboratory) ephemeris files.

    Args:
        filename: Name of the JPL ephemeris file (e.g., "de421.bsp", "de441.bsp")

    Note:
        Skyfield/libephemeris uses JPL Development Ephemeris (DE) files:
        - de421.bsp: 1900-2050 (default, 16 MB)
        - de422.bsp: -3000-3000 (623 MB)
        - de430.bsp: 1550-2650 (128 MB)
        - de431.bsp: -13200-17191 (3.4 GB)
        - de440.bsp: 1550-2650 (128 MB)
        - de441.bsp: -13200-17191 (3.4 GB)

        If the file is not found locally, Skyfield will attempt to download it.
        Use set_ephe_path() to specify a local directory containing the file.

    Example:
        >>> from libephemeris import set_jpl_file, set_ephe_path
        >>> set_ephe_path("/path/to/ephemeris/files")
        >>> set_jpl_file("de441.bsp")  # Use local de441.bsp file
    """
    set_ephemeris_file(filename)


def set_tid_acc(value: float) -> None:
    """
    Set the tidal acceleration used in Delta T calculations.

    The tidal acceleration of the Moon affects the long-term extrapolation
    of Delta T (TT - UT1), which is important for historical astronomical
    calculations. Different JPL ephemeris files assume different values.

    Args:
        value: Tidal acceleration in arcsec/century^2.
               Use SE_TIDAL_* constants for standard ephemeris values,
               or SE_TIDAL_AUTOMATIC (0.0) to use the default.

    Note:
        - The default value is based on DE431 (-25.80 arcsec/cy^2)
        - This setting primarily affects Delta T calculations for dates
          far from the present (historical studies, ancient astronomy)
        - For modern dates (1900-2100), the effect is negligible
        - Common values:
          - DE421: -25.85 arcsec/cy^2
          - DE431: -25.80 arcsec/cy^2 (default)
          - DE441: -25.936 arcsec/cy^2

    Example:
        >>> from libephemeris import set_tid_acc, get_tid_acc
        >>> from libephemeris import SE_TIDAL_DE421, SE_TIDAL_DE431
        >>> set_tid_acc(SE_TIDAL_DE421)  # Use DE421 tidal acceleration
        >>> get_tid_acc()
        -25.85
        >>> set_tid_acc(SE_TIDAL_DE431)  # Use DE431 tidal acceleration
        >>> get_tid_acc()
        -25.8
    """
    global _TIDAL_ACCELERATION
    _TIDAL_ACCELERATION = value if value != 0.0 else None


def get_tid_acc() -> float:
    """
    Get the current tidal acceleration used in Delta T calculations.

    Returns:
        float: The tidal acceleration in arcsec/century^2.
               Returns SE_TIDAL_DEFAULT (-25.80, based on DE431) if not explicitly set.

    Note:
        The tidal acceleration affects how Delta T is extrapolated for dates
        outside the range of direct observations. This is particularly
        important for historical astronomical calculations.

        If set_tid_acc() was called with SE_TIDAL_AUTOMATIC (0.0),
        this returns the default value (DE431-based).

    Example:
        >>> from libephemeris import get_tid_acc, set_tid_acc, SE_TIDAL_DE441
        >>> get_tid_acc()  # Default value
        -25.8
        >>> set_tid_acc(SE_TIDAL_DE441)
        >>> get_tid_acc()
        -25.936
    """
    from .constants import SE_TIDAL_DEFAULT

    if _TIDAL_ACCELERATION is None:
        return SE_TIDAL_DEFAULT
    return _TIDAL_ACCELERATION


def set_delta_t_userdef(dt: Optional[float]) -> None:
    """
    Set a user-defined Delta T value to use instead of computed values.

    When set, swe_deltat() and swe_deltat_ex() will return this fixed value
    instead of computing Delta T from IERS data. This is useful for:
    - Testing and reproducibility
    - Very ancient dates where Delta T is highly uncertain
    - Very future dates where Delta T cannot be predicted
    - Experimentation with different Delta T assumptions

    Args:
        dt: Delta T value in days (TT - UT1), or None to clear and resume
            using computed values. The value should be in the same units
            as returned by swe_deltat() (days, not seconds).

    Note:
        - To convert from seconds to days, divide by 86400
        - Modern Delta T (year 2000) is approximately 64 seconds (~0.00074 days)
        - For ancient dates, Delta T can be hours or even days
        - Use get_delta_t_userdef() to check if a user-defined value is active

    Example:
        >>> from libephemeris import set_delta_t_userdef, get_delta_t_userdef, swe_deltat
        >>> # Set a fixed Delta T of 65 seconds (in days)
        >>> set_delta_t_userdef(65.0 / 86400.0)
        >>> get_delta_t_userdef()
        7.523148148148148e-04
        >>> swe_deltat(2451545.0)  # Now returns fixed value
        7.523148148148148e-04
        >>> # Clear to resume computed values
        >>> set_delta_t_userdef(None)
        >>> get_delta_t_userdef() is None
        True
    """
    global _DELTA_T_USERDEF
    _DELTA_T_USERDEF = dt


def get_delta_t_userdef() -> Optional[float]:
    """
    Get the current user-defined Delta T value.

    Returns:
        Optional[float]: The user-defined Delta T in days if set,
                        or None if using computed values.

    Note:
        When this returns None, swe_deltat() computes Delta T from IERS data.
        When this returns a float, that value is used directly by swe_deltat().

    Example:
        >>> from libephemeris import set_delta_t_userdef, get_delta_t_userdef
        >>> get_delta_t_userdef()  # No user-defined value
        None
        >>> set_delta_t_userdef(0.001)  # Set ~86 seconds
        >>> get_delta_t_userdef()
        0.001
        >>> set_delta_t_userdef(None)  # Clear
        >>> get_delta_t_userdef()
        None
    """
    return _DELTA_T_USERDEF
