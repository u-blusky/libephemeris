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
import threading
from typing import Optional, Union
from skyfield.api import Loader, Topos
from skyfield.timelib import Timescale
from skyfield.jpllib import SpiceKernel

# =============================================================================
# GLOBAL STATE VARIABLES
# =============================================================================

# Lock to protect global state during context-swap operations.
# This ensures thread-safety when EphemerisContext temporarily swaps
# global state for calculations. Without this lock, concurrent threads
# could interfere with each other's state during the save-set-restore cycle.
# We use RLock (reentrant lock) to allow nested locking in case any
# calculation function internally calls other state-swapping functions.
_CONTEXT_SWAP_LOCK = threading.RLock()

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
_LAPSE_RATE: Optional[float] = None  # Atmospheric lapse rate for refraction


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
    mode = _SIDEREAL_MODE if _SIDEREAL_MODE is not None else 1
    if full:
        return mode, _SIDEREAL_T0, _SIDEREAL_AYAN_T0
    return mode


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


# Default atmospheric lapse rate in K/m (standard atmosphere)
SE_LAPSE_RATE_DEFAULT: float = 0.0065


def set_lapse_rate(lapse_rate: Optional[float]) -> None:
    """
    Set the atmospheric temperature lapse rate for refraction calculations.

    The lapse rate is the rate at which temperature decreases with altitude
    in the atmosphere. It affects the calculation of atmospheric refraction,
    particularly the dip of the horizon for elevated observers.

    Args:
        lapse_rate: Temperature lapse rate dT/dh in degrees Kelvin per meter.
                    Use None to reset to the default value (0.0065 K/m).
                    Typical values range from 0.0034 to 0.010 K/m.

    Note:
        - Standard atmospheric lapse rate is 0.0065 K/m (6.5°C per 1000m)
        - This setting affects refrac_extended() calculations
        - Higher lapse rates result in less refraction of the horizon
        - The lapse rate is used in calculating the dip of the horizon
          for elevated observers

    Example:
        >>> from libephemeris import set_lapse_rate, get_lapse_rate
        >>> get_lapse_rate()  # Default value
        0.0065
        >>> set_lapse_rate(0.005)  # Set custom lapse rate
        >>> get_lapse_rate()
        0.005
        >>> set_lapse_rate(None)  # Reset to default
        >>> get_lapse_rate()
        0.0065
    """
    global _LAPSE_RATE
    _LAPSE_RATE = lapse_rate


def get_lapse_rate() -> float:
    """
    Get the current atmospheric temperature lapse rate.

    Returns:
        float: The lapse rate in degrees Kelvin per meter.
               Returns 0.0065 K/m (standard atmosphere) if not explicitly set.

    Note:
        The lapse rate affects refrac_extended() calculations, particularly
        the dip of the horizon for elevated observers.

    Example:
        >>> from libephemeris import get_lapse_rate, set_lapse_rate
        >>> get_lapse_rate()  # Default value
        0.0065
        >>> set_lapse_rate(0.008)
        >>> get_lapse_rate()
        0.008
    """
    if _LAPSE_RATE is None:
        return SE_LAPSE_RATE_DEFAULT
    return _LAPSE_RATE


def get_library_path() -> str:
    """
    Get the path where ephemeris files are stored.

    For libephemeris, this returns the directory containing the JPL ephemeris
    files (.bsp files used by Skyfield). This is analogous to pyswisseph's
    get_library_path(), which returns the path of the Swiss Ephemeris library.

    Returns:
        str: The absolute path to the directory containing ephemeris files.
             If set_ephe_path() was called, returns that custom path.
             Otherwise returns the default data directory (parent of the
             libephemeris package).

    Note:
        - If a custom ephemeris path was set via set_ephe_path(), that path
          is returned regardless of whether files actually exist there.
        - Otherwise, returns the workspace root directory where de421.bsp
          and other data files are located by default.

    Example:
        >>> from libephemeris import get_library_path, set_ephe_path
        >>> get_library_path()  # Returns default path
        '/path/to/workspace'
        >>> set_ephe_path('/custom/ephemeris/path')
        >>> get_library_path()  # Returns custom path
        '/custom/ephemeris/path'
    """
    if _EPHEMERIS_PATH is not None:
        return os.path.abspath(_EPHEMERIS_PATH)
    # Default: parent directory of this module (workspace root)
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def close() -> None:
    """
    Close all opened ephemeris files and release resources.

    This function closes the SPK kernel file handles and resets all global
    state to its initial values. Call this when you're done using the library
    or want to reload ephemeris files with different settings.

    Note:
        - This closes the underlying file handles in the JPL ephemeris kernel
        - After calling close(), the next ephemeris calculation will
          automatically reload the ephemeris files as needed
        - This is useful for:
          - Freeing memory and file handles in long-running applications
          - Switching to a different ephemeris file
          - Ensuring clean state in test suites

    Example:
        >>> from libephemeris import calc_ut, close, SE_SUN
        >>> pos, _ = calc_ut(2451545.0, SE_SUN, 0)  # Loads ephemeris
        >>> close()  # Close files and reset state
        >>> pos, _ = calc_ut(2451545.0, SE_SUN, 0)  # Reloads ephemeris
    """
    global _EPHEMERIS_PATH, _EPHEMERIS_FILE, _LOADER, _PLANETS, _TS
    global _TOPO, _SIDEREAL_MODE, _SIDEREAL_AYAN_T0, _SIDEREAL_T0
    global _ANGLES_CACHE, _TIDAL_ACCELERATION, _DELTA_T_USERDEF, _LAPSE_RATE

    # Close the SPK kernel file handles if loaded
    if _PLANETS is not None:
        try:
            _PLANETS.close()
        except (AttributeError, Exception):
            # SpiceKernel may not have close() in all versions,
            # or may already be closed
            pass

    # Reset all global state to initial values
    _EPHEMERIS_PATH = None
    _EPHEMERIS_FILE = "de421.bsp"
    _LOADER = None
    _PLANETS = None
    _TS = None
    _TOPO = None
    _SIDEREAL_MODE = None
    _SIDEREAL_AYAN_T0 = 0.0
    _SIDEREAL_T0 = 0.0
    _ANGLES_CACHE = {}
    _TIDAL_ACCELERATION = None
    _DELTA_T_USERDEF = None
    _LAPSE_RATE = None


def get_current_file_data(ifno: int = 0) -> tuple[str, float, float, int]:
    """
    Get information about the currently loaded ephemeris file.

    This function returns details about the JPL ephemeris file currently in use,
    including the file path and the date range covered by the ephemeris.

    Args:
        ifno: File number indicator (for compatibility with pyswisseph):
              - 0: planet file (main ephemeris - DE421, DE431, etc.)
              - 1: moon file (same as planets in JPL ephemeris)
              - 2: main asteroid file (not applicable for JPL ephemeris)
              - 3: other asteroid file (not applicable for JPL ephemeris)
              - 4: star file (not applicable for JPL ephemeris)

    Returns:
        A tuple containing:
            - path (str): Full path to the ephemeris file, or empty string if
                         no ephemeris is loaded or ifno is not applicable
            - start (float): Start date of the ephemeris file (Julian Day)
            - end (float): End date of the ephemeris file (Julian Day)
            - denum (int): JPL Development Ephemeris number (e.g., 421, 431, 441)
                          or 0 if not determinable

    Note:
        - For libephemeris using JPL/Skyfield ephemerides, ifno values 0 and 1
          return the same data since planets and Moon are in the same file.
        - ifno values 2, 3, 4 return empty data as those file types are not
          used by libephemeris.
        - The ephemeris file must have been loaded (by calling calc_ut or
          get_planets) before this function returns meaningful data.

    Example:
        >>> from libephemeris import get_current_file_data, calc_ut, SE_SUN
        >>> calc_ut(2451545.0, SE_SUN, 0)  # Triggers ephemeris loading
        >>> path, start, end, denum = get_current_file_data(0)
        >>> print(f"Using DE{denum}, covering JD {start:.1f} to {end:.1f}")
        Using DE421, covering JD 2414864.5 to 2471184.5
    """
    # Only file types 0 (planets) and 1 (moon) are meaningful for JPL ephemeris
    # Types 2, 3, 4 (asteroids, stars) are not applicable
    if ifno not in (0, 1):
        return ("", 0.0, 0.0, 0)

    # If no ephemeris loaded, return empty data
    if _PLANETS is None:
        return ("", 0.0, 0.0, 0)

    try:
        # Get the file path
        path = str(_PLANETS.path) if hasattr(_PLANETS, "path") else ""
        if not path and hasattr(_PLANETS, "filename"):
            path = str(_PLANETS.filename)

        # Get start/end dates from the SPK segments
        # All segments in a JPL DE file typically have the same date range
        start_jd = 0.0
        end_jd = 0.0

        if hasattr(_PLANETS, "spk") and hasattr(_PLANETS.spk, "segments"):
            segments = list(_PLANETS.spk.segments)
            if segments:
                # Use the first segment's date range as representative
                first_seg = segments[0]
                start_jd = (
                    float(first_seg.start_jd) if hasattr(first_seg, "start_jd") else 0.0
                )
                end_jd = (
                    float(first_seg.end_jd) if hasattr(first_seg, "end_jd") else 0.0
                )

        # Extract DE number from filename (e.g., "de421.bsp" -> 421)
        denum = 0
        filename = _EPHEMERIS_FILE.lower()
        if filename.startswith("de") and ".bsp" in filename:
            try:
                # Extract digits between "de" and ".bsp"
                num_str = filename[2:].split(".")[0]
                denum = int(num_str)
            except (ValueError, IndexError):
                denum = 0

        return (path, start_jd, end_jd, denum)

    except Exception:
        # Return empty data on any error
        return ("", 0.0, 0.0, 0)
