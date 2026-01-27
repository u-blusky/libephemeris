"""
Automatic SPK download and caching for minor bodies using astroquery.

This module provides automatic SPK file management for minor body calculations:
- Automatic download of SPK kernels from JPL Horizons on demand
- Local caching to avoid redundant downloads
- Seamless integration with calc_ut() for transparent SPK usage

The module uses astroquery.jplhorizons to download SPK files, which provides
a more robust and feature-rich interface than direct HTTP requests.

Usage:
    >>> from libephemeris import spk_auto
    >>> # Enable auto-downloading for a body
    >>> spk_auto.enable_auto_spk(
    ...     ipl=SE_CHIRON,
    ...     body_id="2060",  # Horizons identifier
    ...     start="2000-01-01",
    ...     end="2100-01-01",
    ... )
    >>> # Now calc_ut automatically downloads and uses SPK if needed
    >>> pos, _ = calc_ut(2451545.0, SE_CHIRON, 0)

Requirements:
    pip install astroquery

References:
    - astroquery.jplhorizons: https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html
    - JPL Horizons: https://ssd.jpl.nasa.gov/horizons/
"""

import hashlib
import os
import threading
from typing import Optional, Union

from .state import get_library_path


# =============================================================================
# MODULE STATE
# =============================================================================

# Lock for thread-safe operations
_AUTO_SPK_LOCK = threading.RLock()

# Registry of auto-SPK configurations: {ipl: AutoSpkConfig}
_AUTO_SPK_REGISTRY: dict[int, "AutoSpkConfig"] = {}

# Default cache directory name within library path
DEFAULT_CACHE_DIR = "spk_cache"

# Default date range for auto-downloaded SPK files
DEFAULT_START_DATE = "2000-01-01"
DEFAULT_END_DATE = "2100-01-01"


# =============================================================================
# CACHE DIRECTORY MANAGEMENT
# =============================================================================


def ensure_cache_dir(cache_dir: Optional[str] = None) -> str:
    """
    Ensure the SPK cache directory exists, creating it if necessary.

    This function creates the cache directory if it doesn't exist and returns
    the absolute path to the cache directory.

    Args:
        cache_dir: Optional custom cache directory path. If None, uses the
            default cache directory ({library_path}/spk_cache/).

    Returns:
        str: Absolute path to the cache directory.

    Example:
        >>> from libephemeris.spk_auto import ensure_cache_dir
        >>> cache_path = ensure_cache_dir()
        >>> print(cache_path)
        /path/to/libephemeris/spk_cache
        >>> # Custom cache directory
        >>> cache_path = ensure_cache_dir("/custom/cache/path")
    """
    if cache_dir is not None:
        cache_path = os.path.abspath(cache_dir)
    else:
        cache_path = os.path.join(get_library_path(), DEFAULT_CACHE_DIR)

    if not os.path.exists(cache_path):
        os.makedirs(cache_path, exist_ok=True)

    return cache_path


def get_cache_path(body_id: Union[int, str], cache_dir: Optional[str] = None) -> str:
    """
    Get the path where an SPK file for a body should be stored.

    This function returns the expected cache file path for a given body ID.
    The path is constructed using the body_id as part of the filename.
    Note: This does not create the file or directory - use ensure_cache_dir()
    to create the cache directory if needed.

    Args:
        body_id: JPL Horizons body identifier. Can be:
            - Asteroid number (int or str): 2060, "2060", "136199"
            - Name (str): "Chiron", "Eris"
        cache_dir: Optional custom cache directory path. If None, uses the
            default cache directory ({library_path}/spk_cache/).

    Returns:
        str: Full path where the SPK file would be stored.
            The filename format is: {sanitized_body_id}.bsp

    Example:
        >>> from libephemeris.spk_auto import get_cache_path
        >>> path = get_cache_path("2060")
        >>> print(path)
        /path/to/libephemeris/spk_cache/2060.bsp
        >>> path = get_cache_path("Chiron")
        >>> print(path)
        /path/to/libephemeris/spk_cache/chiron.bsp
    """
    if cache_dir is not None:
        dir_path = os.path.abspath(cache_dir)
    else:
        dir_path = os.path.join(get_library_path(), DEFAULT_CACHE_DIR)

    # Sanitize body_id for filename
    body_str = str(body_id).lower()
    safe_body_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in body_str)

    filename = f"{safe_body_id}.bsp"
    return os.path.join(dir_path, filename)


def cache_info(
    cache_dir: Optional[str] = None,
) -> dict[str, Union[int, float, str, list[str]]]:
    """
    Get information about the SPK cache directory.

    This function returns statistics about cached SPK files including the
    number of files, total size, and list of cached filenames.

    Args:
        cache_dir: Optional custom cache directory path. If None, uses the
            default cache directory ({library_path}/spk_cache/).

    Returns:
        dict: Cache information with the following keys:
            - cache_dir (str): Absolute path to the cache directory
            - num_files (int): Number of cached SPK files (.bsp files)
            - total_size_mb (float): Total size of all cached files in megabytes
            - files (list[str]): List of cached SPK filenames

    Example:
        >>> from libephemeris.spk_auto import cache_info
        >>> info = cache_info()
        >>> print(f"Cache at: {info['cache_dir']}")
        >>> print(f"Files cached: {info['num_files']}")
        >>> print(f"Total size: {info['total_size_mb']:.2f} MB")
    """
    if cache_dir is not None:
        dir_path = os.path.abspath(cache_dir)
    else:
        dir_path = os.path.join(get_library_path(), DEFAULT_CACHE_DIR)

    if not os.path.exists(dir_path):
        return {
            "cache_dir": dir_path,
            "num_files": 0,
            "total_size_mb": 0.0,
            "files": [],
        }

    files = [f for f in os.listdir(dir_path) if f.endswith(".bsp")]
    total_size = sum(os.path.getsize(os.path.join(dir_path, f)) for f in files)

    return {
        "cache_dir": dir_path,
        "num_files": len(files),
        "total_size_mb": round(total_size / (1024 * 1024), 2),
        "files": files,
    }


# =============================================================================
# CONFIGURATION CLASS
# =============================================================================


class AutoSpkConfig:
    """Configuration for automatic SPK download for a minor body.

    Attributes:
        ipl: libephemeris body ID (e.g., SE_CHIRON)
        body_id: JPL Horizons identifier (e.g., "2060", "Chiron")
        naif_id: NAIF ID for SPK lookup (e.g., 2002060)
        start: Start date for SPK coverage (YYYY-MM-DD)
        end: End date for SPK coverage (YYYY-MM-DD)
        cache_dir: Directory for storing downloaded SPK files
        enabled: Whether auto-download is enabled
        spk_path: Path to downloaded SPK file (None if not yet downloaded)
    """

    def __init__(
        self,
        ipl: int,
        body_id: str,
        naif_id: Optional[int] = None,
        start: str = DEFAULT_START_DATE,
        end: str = DEFAULT_END_DATE,
        cache_dir: Optional[str] = None,
    ):
        self.ipl = ipl
        self.body_id = body_id
        self.naif_id = naif_id
        self.start = start
        self.end = end
        self.cache_dir = cache_dir
        self.enabled = True
        self.spk_path: Optional[str] = None

    def get_cache_dir(self) -> str:
        """Get the cache directory, creating if necessary."""
        from . import state

        cache_path: str
        if self.cache_dir:
            cache_path = self.cache_dir
        elif state.get_spk_cache_dir() is not None:
            cache_path = state.get_spk_cache_dir()  # type: ignore
        else:
            cache_path = os.path.join(get_library_path(), DEFAULT_CACHE_DIR)

        if not os.path.exists(cache_path):
            os.makedirs(cache_path, exist_ok=True)

        return cache_path

    def get_cache_filename(self) -> str:
        """Generate a unique cache filename based on configuration."""
        # Create a hash of the configuration for unique identification
        config_str = f"{self.body_id}_{self.start}_{self.end}"
        config_hash = hashlib.md5(config_str.encode()).hexdigest()[:8]

        # Sanitize body_id for filename
        safe_body_id = "".join(
            c if c.isalnum() or c in "-_" else "_" for c in self.body_id.lower()
        )

        return f"{safe_body_id}_{config_hash}.bsp"

    def get_cache_path(self) -> str:
        """Get the full path where the SPK file should be cached."""
        return os.path.join(self.get_cache_dir(), self.get_cache_filename())


# =============================================================================
# ASTROQUERY-BASED DOWNLOAD
# =============================================================================


def _download_spk_astroquery(
    body_id: str,
    start: str,
    end: str,
    output_path: str,
    location: str = "@0",
) -> str:
    """
    Download SPK file using astroquery.jplhorizons.

    Args:
        body_id: JPL Horizons target identifier (e.g., "2060", "Chiron")
        start: Start date (YYYY-MM-DD)
        end: End date (YYYY-MM-DD)
        output_path: Path to save the SPK file
        location: Observer location code (default: "@0" = solar system barycenter)

    Returns:
        str: Path to the downloaded SPK file

    Raises:
        ImportError: If astroquery is not installed
        ValueError: If body not found or download fails
    """
    try:
        from astroquery.jplhorizons import Horizons
    except ImportError as e:
        raise ImportError(
            "astroquery is required for automatic SPK downloads. "
            "Install it with: pip install astroquery"
        ) from e

    # Create Horizons query object
    obj = Horizons(id=body_id, location=location, epochs={"start": start, "stop": end})

    try:
        # Download SPK file
        # Note: This uses Horizons' SPK file generation capability
        obj.download_spk(output_path)
        return output_path
    except Exception as e:
        raise ValueError(
            f"Failed to download SPK for '{body_id}' from JPL Horizons: {e}"
        ) from e


def _check_astroquery_available() -> bool:
    """Check if astroquery is available."""
    try:
        from astroquery.jplhorizons import Horizons

        return True
    except ImportError:
        return False


# =============================================================================
# PUBLIC API
# =============================================================================


def enable_auto_spk(
    ipl: int,
    body_id: str,
    naif_id: Optional[int] = None,
    start: str = DEFAULT_START_DATE,
    end: str = DEFAULT_END_DATE,
    cache_dir: Optional[str] = None,
) -> None:
    """
    Enable automatic SPK download for a minor body.

    When enabled, the library will automatically download and cache an SPK
    kernel from JPL Horizons when the body's position is calculated.

    Args:
        ipl: libephemeris body ID (e.g., SE_CHIRON, SE_ERIS)
        body_id: JPL Horizons target identifier. Can be:
            - Asteroid number: "2060", "136199"
            - Name: "Chiron", "Eris"
            - Designation: "2003 UB313"
        naif_id: NAIF ID for SPK lookup. If None, will be deduced from body_id
            using convention: asteroid_number + 2000000
        start: Start date for SPK coverage (YYYY-MM-DD). Default: "2000-01-01"
        end: End date for SPK coverage (YYYY-MM-DD). Default: "2100-01-01"
        cache_dir: Directory for storing downloaded SPK files.
            Default: {library_path}/spk_cache/

    Example:
        >>> from libephemeris import spk_auto, SE_CHIRON
        >>> spk_auto.enable_auto_spk(
        ...     ipl=SE_CHIRON,
        ...     body_id="2060",
        ...     start="2000-01-01",
        ...     end="2100-01-01",
        ... )
        >>> # Future calculations for Chiron will use SPK automatically

    Note:
        This function only registers the configuration. The actual download
        occurs lazily when the body's position is first calculated, or
        immediately if download_now() is called.
    """
    if not _check_astroquery_available():
        raise ImportError(
            "astroquery is required for automatic SPK downloads. "
            "Install it with: pip install astroquery"
        )

    # Import state to get configuration settings
    from . import state

    # Use parameter if provided, otherwise fall back to global config
    if cache_dir is None:
        cache_dir = state.get_spk_cache_dir()

    config = AutoSpkConfig(
        ipl=ipl,
        body_id=body_id,
        naif_id=naif_id,
        start=start,
        end=end,
        cache_dir=cache_dir,
    )

    with _AUTO_SPK_LOCK:
        _AUTO_SPK_REGISTRY[ipl] = config


def disable_auto_spk(ipl: int) -> None:
    """
    Disable automatic SPK download for a body.

    This removes the auto-SPK configuration but does not delete any
    already-downloaded SPK files or unregister them from the SPK registry.

    Args:
        ipl: libephemeris body ID (e.g., SE_CHIRON)
    """
    with _AUTO_SPK_LOCK:
        if ipl in _AUTO_SPK_REGISTRY:
            del _AUTO_SPK_REGISTRY[ipl]


def is_auto_spk_enabled(ipl: int) -> bool:
    """
    Check if automatic SPK download is enabled for a body.

    Args:
        ipl: libephemeris body ID

    Returns:
        bool: True if auto-SPK is enabled for this body
    """
    with _AUTO_SPK_LOCK:
        return ipl in _AUTO_SPK_REGISTRY and _AUTO_SPK_REGISTRY[ipl].enabled


def get_auto_spk_config(ipl: int) -> Optional[AutoSpkConfig]:
    """
    Get the auto-SPK configuration for a body.

    Args:
        ipl: libephemeris body ID

    Returns:
        AutoSpkConfig or None if not configured
    """
    with _AUTO_SPK_LOCK:
        return _AUTO_SPK_REGISTRY.get(ipl)


def list_auto_spk_bodies() -> dict[int, AutoSpkConfig]:
    """
    List all bodies with auto-SPK enabled.

    Returns:
        Dict mapping ipl -> AutoSpkConfig
    """
    with _AUTO_SPK_LOCK:
        return dict(_AUTO_SPK_REGISTRY)


def download_now(ipl: int, force: bool = False) -> str:
    """
    Download SPK for a body immediately.

    This triggers immediate download of the SPK file for a body that has
    auto-SPK enabled, rather than waiting for the first calculation.

    Args:
        ipl: libephemeris body ID (must have auto-SPK enabled)
        force: If True, re-download even if cached file exists

    Returns:
        str: Path to the downloaded SPK file

    Raises:
        ValueError: If auto-SPK is not enabled for this body
        ImportError: If astroquery is not installed
    """
    config = get_auto_spk_config(ipl)
    if config is None:
        raise ValueError(
            f"Auto-SPK not enabled for body {ipl}. Call enable_auto_spk() first."
        )

    return _ensure_spk_downloaded(config, force=force)


def _ensure_spk_downloaded(config: AutoSpkConfig, force: bool = False) -> str:
    """
    Ensure SPK is downloaded and registered for a body.

    Internal function that handles the download and registration logic.

    Args:
        config: Auto-SPK configuration
        force: If True, re-download even if cached

    Returns:
        str: Path to the SPK file
    """
    from . import spk

    cache_path = config.get_cache_path()

    # Check if already cached
    if os.path.exists(cache_path) and not force:
        config.spk_path = cache_path
    else:
        # Download using astroquery
        with _AUTO_SPK_LOCK:
            # Double-check after acquiring lock
            if not os.path.exists(cache_path) or force:
                _download_spk_astroquery(
                    body_id=config.body_id,
                    start=config.start,
                    end=config.end,
                    output_path=cache_path,
                )
        config.spk_path = cache_path

    # Deduce NAIF ID if not provided
    naif_id = config.naif_id
    if naif_id is None:
        naif_id = spk._deduce_naif_id(config.body_id)
        if naif_id is None:
            raise ValueError(
                f"Cannot deduce NAIF ID for '{config.body_id}'. "
                "Please provide naif_id explicitly in enable_auto_spk()."
            )
        config.naif_id = naif_id

    # Register the SPK body if not already registered
    from . import state

    if config.ipl not in state._SPK_BODY_MAP:
        spk.register_spk_body(config.ipl, cache_path, naif_id)

    return cache_path


def try_auto_download(ipl: int) -> Optional[str]:
    """
    Try to auto-download SPK for a body if configured.

    This function is called internally during position calculations
    to transparently download SPK files when needed.

    Args:
        ipl: libephemeris body ID

    Returns:
        Path to SPK file if downloaded, None if not configured or failed
    """
    config = get_auto_spk_config(ipl)
    if config is None or not config.enabled:
        return None

    try:
        return _ensure_spk_downloaded(config)
    except Exception:
        # Log error but don't raise - fall back to Keplerian
        return None


def clear_cache(ipl: Optional[int] = None) -> int:
    """
    Clear cached SPK files.

    Args:
        ipl: If provided, only clear cache for this body.
             If None, clear all cached SPK files.

    Returns:
        int: Number of files deleted
    """
    deleted = 0

    with _AUTO_SPK_LOCK:
        if ipl is not None:
            config = _AUTO_SPK_REGISTRY.get(ipl)
            if config and config.spk_path and os.path.exists(config.spk_path):
                os.remove(config.spk_path)
                config.spk_path = None
                deleted = 1
        else:
            # Clear all cached files
            for config in _AUTO_SPK_REGISTRY.values():
                if config.spk_path and os.path.exists(config.spk_path):
                    try:
                        os.remove(config.spk_path)
                        config.spk_path = None
                        deleted += 1
                    except OSError:
                        pass

    return deleted


def get_cache_info() -> dict[str, Union[int, float, str, list[str]]]:
    """
    Get information about the SPK cache.

    Returns:
        Dict with cache information:
            - cache_dir: Path to cache directory
            - num_files: Number of cached SPK files
            - total_size_mb: Total size in MB
            - files: List of cached file names
    """
    cache_dir = os.path.join(get_library_path(), DEFAULT_CACHE_DIR)

    if not os.path.exists(cache_dir):
        return {
            "cache_dir": cache_dir,
            "num_files": 0,
            "total_size_mb": 0.0,
            "files": [],
        }

    files = [f for f in os.listdir(cache_dir) if f.endswith(".bsp")]
    total_size = sum(os.path.getsize(os.path.join(cache_dir, f)) for f in files)

    return {
        "cache_dir": cache_dir,
        "num_files": len(files),
        "total_size_mb": round(total_size / (1024 * 1024), 2),
        "files": files,
    }


# =============================================================================
# CACHE MANAGEMENT FUNCTIONS
# =============================================================================


def list_cached_spk(
    cache_dir: Optional[str] = None,
) -> list[dict[str, Union[str, float, int, None]]]:
    """
    List all cached SPK files with their date ranges and sizes.

    This function scans the cache directory for SPK files and returns detailed
    information about each file including its size, date range coverage, and
    last access time.

    Args:
        cache_dir: Optional custom cache directory path. If None, checks both:
            - The default cache directory ({library_path}/spk_cache/)
            - The auto SPK directory (~/.libephemeris/spk/)

    Returns:
        list: A list of dictionaries, each containing:
            - filename (str): Name of the SPK file
            - path (str): Full path to the file
            - size_bytes (int): File size in bytes
            - size_mb (float): File size in megabytes
            - jd_start (float or None): Start Julian Day from SPK coverage (if readable)
            - jd_end (float or None): End Julian Day from SPK coverage (if readable)
            - date_start (str or None): Start date in YYYY-MM-DD format (if readable)
            - date_end (str or None): End date in YYYY-MM-DD format (if readable)
            - last_accessed (float): Last access time as Unix timestamp

    Example:
        >>> from libephemeris.spk_auto import list_cached_spk
        >>> cached_files = list_cached_spk()
        >>> for f in cached_files:
        ...     print(f"{f['filename']}: {f['size_mb']:.2f} MB")
        ...     if f['date_start']:
        ...         print(f"  Coverage: {f['date_start']} to {f['date_end']}")
    """
    from .spk import get_spk_coverage

    result = []
    dirs_to_check = []

    if cache_dir is not None:
        dirs_to_check = [os.path.abspath(cache_dir)]
    else:
        # Check both default cache directories
        default_cache = os.path.join(get_library_path(), DEFAULT_CACHE_DIR)
        dirs_to_check = [default_cache, DEFAULT_AUTO_SPK_DIR]

    for dir_path in dirs_to_check:
        if not os.path.exists(dir_path):
            continue

        for filename in os.listdir(dir_path):
            if not filename.endswith(".bsp"):
                continue

            filepath = os.path.join(dir_path, filename)
            stat_info = os.stat(filepath)

            file_info: dict[str, Union[str, float, int, None]] = {
                "filename": filename,
                "path": filepath,
                "size_bytes": stat_info.st_size,
                "size_mb": round(stat_info.st_size / (1024 * 1024), 2),
                "jd_start": None,
                "jd_end": None,
                "date_start": None,
                "date_end": None,
                "last_accessed": stat_info.st_atime,
            }

            # Try to get actual coverage from SPK file
            try:
                coverage = get_spk_coverage(filepath)
                if coverage is not None:
                    jd_start, jd_end = coverage
                    file_info["jd_start"] = jd_start
                    file_info["jd_end"] = jd_end
                    file_info["date_start"] = _jd_to_iso_date(jd_start)
                    file_info["date_end"] = _jd_to_iso_date(jd_end)
            except Exception:
                # If we can't read coverage, try to extract from filename
                parts = filename[:-4].split("_")  # Remove .bsp
                if len(parts) >= 3:
                    try:
                        file_jd_start = float(parts[-2])
                        file_jd_end = float(parts[-1])
                        file_info["jd_start"] = file_jd_start
                        file_info["jd_end"] = file_jd_end
                        file_info["date_start"] = _jd_to_iso_date(file_jd_start)
                        file_info["date_end"] = _jd_to_iso_date(file_jd_end)
                    except ValueError:
                        pass

            result.append(file_info)

    return result


def clear_spk_cache(cache_dir: Optional[str] = None) -> int:
    """
    Delete all cached SPK files.

    This function removes all .bsp files from the cache directory. It does not
    delete the directory itself, only the SPK files within it.

    Warning: This operation cannot be undone. All cached SPK files will need
    to be re-downloaded if needed.

    Args:
        cache_dir: Optional custom cache directory path. If None, clears both:
            - The default cache directory ({library_path}/spk_cache/)
            - The auto SPK directory (~/.libephemeris/spk/)

    Returns:
        int: Number of files deleted

    Example:
        >>> from libephemeris.spk_auto import clear_spk_cache, get_cache_size
        >>> print(f"Cache size before: {get_cache_size():.2f} MB")
        Cache size before: 15.32 MB
        >>> deleted = clear_spk_cache()
        >>> print(f"Deleted {deleted} files")
        Deleted 5 files
        >>> print(f"Cache size after: {get_cache_size():.2f} MB")
        Cache size after: 0.00 MB
    """
    deleted_count = 0
    dirs_to_clear = []

    if cache_dir is not None:
        dirs_to_clear = [os.path.abspath(cache_dir)]
    else:
        # Clear both default cache directories
        default_cache = os.path.join(get_library_path(), DEFAULT_CACHE_DIR)
        dirs_to_clear = [default_cache, DEFAULT_AUTO_SPK_DIR]

    for dir_path in dirs_to_clear:
        if not os.path.exists(dir_path):
            continue

        for filename in os.listdir(dir_path):
            if not filename.endswith(".bsp"):
                continue

            filepath = os.path.join(dir_path, filename)
            try:
                os.remove(filepath)
                deleted_count += 1
            except OSError:
                # Skip files that can't be deleted (e.g., permission issues)
                pass

    # Also clear the spk_path from any registered configurations
    with _AUTO_SPK_LOCK:
        for config in _AUTO_SPK_REGISTRY.values():
            if config.spk_path and not os.path.exists(config.spk_path):
                config.spk_path = None

    return deleted_count


def get_cache_size(cache_dir: Optional[str] = None) -> float:
    """
    Get the total disk usage of the SPK cache in megabytes.

    This function calculates the combined size of all .bsp files in the
    cache directory.

    Args:
        cache_dir: Optional custom cache directory path. If None, calculates
            the combined size of both:
            - The default cache directory ({library_path}/spk_cache/)
            - The auto SPK directory (~/.libephemeris/spk/)

    Returns:
        float: Total size of cached SPK files in megabytes

    Example:
        >>> from libephemeris.spk_auto import get_cache_size
        >>> size = get_cache_size()
        >>> print(f"Total cache size: {size:.2f} MB")
        Total cache size: 15.32 MB
        >>> if size > 100:
        ...     print("Consider pruning old cache files")
    """
    total_size = 0
    dirs_to_check = []

    if cache_dir is not None:
        dirs_to_check = [os.path.abspath(cache_dir)]
    else:
        # Check both default cache directories
        default_cache = os.path.join(get_library_path(), DEFAULT_CACHE_DIR)
        dirs_to_check = [default_cache, DEFAULT_AUTO_SPK_DIR]

    for dir_path in dirs_to_check:
        if not os.path.exists(dir_path):
            continue

        for filename in os.listdir(dir_path):
            if not filename.endswith(".bsp"):
                continue

            filepath = os.path.join(dir_path, filename)
            try:
                total_size += os.path.getsize(filepath)
            except OSError:
                pass

    return round(total_size / (1024 * 1024), 2)


def prune_old_cache(
    max_age_days: int,
    cache_dir: Optional[str] = None,
) -> int:
    """
    Remove SPK files that have not been accessed recently.

    This function removes cached SPK files whose last access time exceeds
    the specified age threshold. This helps manage disk space by removing
    files that are no longer in active use.

    Args:
        max_age_days: Maximum age in days since last access. Files not accessed
            within this period will be deleted. Must be a positive integer.
        cache_dir: Optional custom cache directory path. If None, prunes both:
            - The default cache directory ({library_path}/spk_cache/)
            - The auto SPK directory (~/.libephemeris/spk/)

    Returns:
        int: Number of files deleted

    Raises:
        ValueError: If max_age_days is not a positive integer

    Example:
        >>> from libephemeris.spk_auto import prune_old_cache, list_cached_spk
        >>> # Remove files not accessed in the last 30 days
        >>> pruned = prune_old_cache(max_age_days=30)
        >>> print(f"Removed {pruned} old cache files")
        Removed 3 old cache files
        >>> # Check remaining files
        >>> remaining = list_cached_spk()
        >>> print(f"Remaining cached files: {len(remaining)}")

    Note:
        The access time is based on the file's last access timestamp (atime),
        which may not be updated on all file systems. Some file systems mount
        with 'noatime' for performance, which would prevent this from working
        correctly.
    """
    import time

    if not isinstance(max_age_days, int) or max_age_days <= 0:
        raise ValueError("max_age_days must be a positive integer")

    deleted_count = 0
    current_time = time.time()
    max_age_seconds = max_age_days * 24 * 60 * 60
    dirs_to_prune = []

    if cache_dir is not None:
        dirs_to_prune = [os.path.abspath(cache_dir)]
    else:
        # Prune both default cache directories
        default_cache = os.path.join(get_library_path(), DEFAULT_CACHE_DIR)
        dirs_to_prune = [default_cache, DEFAULT_AUTO_SPK_DIR]

    for dir_path in dirs_to_prune:
        if not os.path.exists(dir_path):
            continue

        for filename in os.listdir(dir_path):
            if not filename.endswith(".bsp"):
                continue

            filepath = os.path.join(dir_path, filename)
            try:
                stat_info = os.stat(filepath)
                last_access = stat_info.st_atime
                age_seconds = current_time - last_access

                if age_seconds > max_age_seconds:
                    os.remove(filepath)
                    deleted_count += 1
            except OSError:
                # Skip files that can't be accessed or deleted
                pass

    # Also clear the spk_path from any registered configurations for deleted files
    with _AUTO_SPK_LOCK:
        for config in _AUTO_SPK_REGISTRY.values():
            if config.spk_path and not os.path.exists(config.spk_path):
                config.spk_path = None

    return deleted_count


# =============================================================================
# PRESET CONFIGURATIONS
# =============================================================================


def enable_common_bodies(
    start: str = DEFAULT_START_DATE,
    end: str = DEFAULT_END_DATE,
    cache_dir: Optional[str] = None,
) -> None:
    """
    Enable auto-SPK for commonly used minor bodies.

    This is a convenience function that enables automatic SPK downloads
    for popular asteroids and TNOs used in astrological and astronomical
    calculations.

    Bodies enabled:
        - Chiron (2060)
        - Pholus (5145)
        - Ceres (1)
        - Pallas (2)
        - Juno (3)
        - Vesta (4)
        - Eris (136199)
        - Sedna (90377)

    Args:
        start: Start date for SPK coverage (YYYY-MM-DD)
        end: End date for SPK coverage (YYYY-MM-DD)
        cache_dir: Directory for storing downloaded SPK files

    Example:
        >>> from libephemeris import spk_auto
        >>> spk_auto.enable_common_bodies()
        >>> # Now all common minor bodies use SPK automatically
    """
    from .constants import (
        SE_CHIRON,
        SE_PHOLUS,
        SE_CERES,
        SE_PALLAS,
        SE_JUNO,
        SE_VESTA,
        SE_ERIS,
        SE_SEDNA,
        NAIF_CHIRON,
        NAIF_PHOLUS,
        NAIF_CERES,
        NAIF_PALLAS,
        NAIF_JUNO,
        NAIF_VESTA,
        NAIF_ERIS,
        NAIF_SEDNA,
    )

    common_bodies = [
        (SE_CHIRON, "2060", NAIF_CHIRON),
        (SE_PHOLUS, "5145", NAIF_PHOLUS),
        (SE_CERES, "1", NAIF_CERES),
        (SE_PALLAS, "2", NAIF_PALLAS),
        (SE_JUNO, "3", NAIF_JUNO),
        (SE_VESTA, "4", NAIF_VESTA),
        (SE_ERIS, "136199", NAIF_ERIS),
        (SE_SEDNA, "90377", NAIF_SEDNA),
    ]

    for ipl, body_id, naif_id in common_bodies:
        enable_auto_spk(
            ipl=ipl,
            body_id=body_id,
            naif_id=naif_id,
            start=start,
            end=end,
            cache_dir=cache_dir,
        )


def disable_all() -> None:
    """
    Disable auto-SPK for all bodies.

    This clears the auto-SPK registry but does not delete cached files
    or unregister already-registered SPK bodies.
    """
    with _AUTO_SPK_LOCK:
        _AUTO_SPK_REGISTRY.clear()


# =============================================================================
# DEFAULT CACHE DIRECTORY
# =============================================================================

# Default cache directory for auto_get_spk
DEFAULT_AUTO_SPK_DIR = os.path.join(os.path.expanduser("~"), ".libephemeris", "spk")


def _jd_to_iso_date(jd: float) -> str:
    """
    Convert Julian Day to ISO date string (YYYY-MM-DD).

    Args:
        jd: Julian Day number

    Returns:
        str: ISO date string (YYYY-MM-DD)
    """
    # Julian Day to calendar date conversion
    # Algorithm from Meeus, Astronomical Algorithms
    jd = jd + 0.5
    z = int(jd)
    f = jd - z

    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - int(alpha / 4)

    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    day = b - d - int(30.6001 * e)

    if e < 14:
        month = e - 1
    else:
        month = e - 13

    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    return f"{year:04d}-{month:02d}-{day:02d}"


def _generate_spk_cache_filename(
    body_id: Union[int, str], jd_start: float, jd_end: float
) -> str:
    """
    Generate a cache filename for an SPK file.

    The filename includes body identifier and truncated Julian Day range
    to enable cache lookup.

    Args:
        body_id: Body identifier
        jd_start: Start Julian Day
        jd_end: End Julian Day

    Returns:
        str: Filename like "2060_2458849_2462502.bsp"
    """
    # Sanitize body_id for filename
    body_str = str(body_id).lower()
    safe_body_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in body_str)

    # Truncate JD to integer for filename (sufficient precision for cache lookup)
    jd_start_int = int(jd_start)
    jd_end_int = int(jd_end)

    return f"{safe_body_id}_{jd_start_int}_{jd_end_int}.bsp"


def _find_covering_spk(
    body_id: Union[int, str],
    jd_start: float,
    jd_end: float,
    cache_dir: str,
) -> Optional[str]:
    """
    Find an existing SPK file in the cache that covers the requested date range.

    Searches the cache directory for any SPK file for the given body that
    covers the requested Julian Day range.

    Args:
        body_id: Body identifier
        jd_start: Start Julian Day
        jd_end: End Julian Day
        cache_dir: Cache directory path

    Returns:
        str: Path to covering SPK file if found, None otherwise
    """
    if not os.path.exists(cache_dir):
        return None

    # Sanitize body_id to match filename pattern
    body_str = str(body_id).lower()
    safe_body_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in body_str)

    # Look for files matching this body
    for filename in os.listdir(cache_dir):
        if not filename.endswith(".bsp"):
            continue

        # Check if filename matches our body
        if not filename.startswith(safe_body_id + "_"):
            continue

        # Extract JD range from filename (format: body_jdstart_jdend.bsp)
        parts = filename[:-4].split("_")  # Remove .bsp
        if len(parts) >= 3:
            try:
                # Get the last two parts as JD values
                file_jd_start = float(parts[-2])
                file_jd_end = float(parts[-1])

                # Check if this file covers our requested range
                # We compare truncated integer values since filenames use truncated JDs
                # Add 1 to file_jd_end to account for the truncation (covers the full day)
                if file_jd_start <= int(jd_start) and file_jd_end >= int(jd_end):
                    return os.path.join(cache_dir, filename)
            except ValueError:
                continue

    return None


def is_spk_cached(
    body_id: Union[int, str],
    jd_start: float,
    jd_end: float,
    cache_dir: Optional[str] = None,
) -> bool:
    """
    Check if a cached SPK file exists and covers the requested date range.

    This function searches the cache directory for SPK files matching the given
    body_id and verifies that the file covers the requested Julian Day range by
    parsing the SPK file header (not just relying on filename patterns).

    Args:
        body_id: JPL Horizons body identifier. Can be:
            - Asteroid number (int): 2060, 136199
            - Name (str): "Chiron", "Eris"
            - Combined: "2060 Chiron"
        jd_start: Start of the date range (Julian Day)
        jd_end: End of the date range (Julian Day)
        cache_dir: Directory to search for SPK files.
            Default: ~/.libephemeris/spk/

    Returns:
        bool: True if a cached SPK file exists that covers the requested range,
              False otherwise.

    Example:
        >>> from libephemeris.spk_auto import is_spk_cached
        >>> # Check if Chiron SPK is cached for 2020-2030
        >>> jd_start = 2458849.5  # 2020-01-01
        >>> jd_end = 2462502.5    # 2030-01-01
        >>> if is_spk_cached("2060", jd_start, jd_end):
        ...     print("Chiron SPK is cached")
        ... else:
        ...     print("Need to download Chiron SPK")
    """
    # Determine cache directory
    if cache_dir is None:
        cache_dir = DEFAULT_AUTO_SPK_DIR

    if not os.path.exists(cache_dir):
        return False

    # Sanitize body_id to match filename pattern
    body_str = str(body_id).lower()
    safe_body_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in body_str)

    # Import get_spk_coverage here to avoid circular imports
    from .spk import get_spk_coverage

    # Look for files matching this body
    for filename in os.listdir(cache_dir):
        if not filename.endswith(".bsp"):
            continue

        # Check if filename matches our body
        if not filename.startswith(safe_body_id + "_"):
            continue

        # Found a potential match - parse the actual SPK file to get coverage
        spk_path = os.path.join(cache_dir, filename)

        try:
            coverage = get_spk_coverage(spk_path)
            if coverage is not None:
                file_jd_start, file_jd_end = coverage
                # Check if this file covers our requested range
                if file_jd_start <= jd_start and file_jd_end >= jd_end:
                    return True
        except Exception:
            # If we can't read the file, skip it and try others
            continue

    return False


def auto_get_spk(
    body_id: Union[int, str],
    jd_start: float,
    jd_end: float,
    cache_dir: Optional[str] = None,
    ipl: Optional[int] = None,
    naif_id: Optional[int] = None,
) -> str:
    """
    Automatically get SPK file for a body and date range, downloading if needed.

    This function checks if an SPK file for the requested body and date range
    exists in the local cache directory. If not, it automatically downloads
    the file from JPL Horizons using the Horizons class.

    If ``ipl`` is provided, the SPK body is automatically registered with
    libephemeris after download, making it immediately usable by ``calc_ut()``.

    Args:
        body_id: JPL Horizons body identifier. Can be:
            - Asteroid number (int): 2060, 136199
            - Name (str): "Chiron", "Eris"
            - Combined: "2060 Chiron"
        jd_start: Start of the date range (Julian Day)
        jd_end: End of the date range (Julian Day)
        cache_dir: Directory for storing/finding SPK files.
            Default: Uses set_spk_cache_dir() value if set, otherwise
            ~/.libephemeris/spk/
        ipl: Optional libephemeris body ID (e.g., SE_CHIRON). If provided,
            the SPK body is automatically registered after download.
        naif_id: Optional NAIF ID for SPK lookup. If not provided but ipl is,
            the NAIF ID is deduced from body_id using the convention:
            naif_id = asteroid_number + 2000000

    Returns:
        str: Path to the SPK file (either existing cached or newly downloaded)

    Raises:
        ImportError: If astroquery is not installed
        ValueError: If body not found, download fails, or naif_id cannot be
            deduced when ipl is provided

    Note:
        The date range may be extended by the value set via set_spk_date_padding()
        to provide a buffer around the requested range.

    Example:
        >>> from libephemeris.spk_auto import auto_get_spk
        >>> from libephemeris.constants import SE_CHIRON
        >>> # Get SPK for Chiron covering 2020-2030 and register it
        >>> jd_start = 2458849.5  # 2020-01-01
        >>> jd_end = 2462502.5    # 2030-01-01
        >>> spk_path = auto_get_spk("2060", jd_start, jd_end, ipl=SE_CHIRON)
        >>> print(spk_path)
        /home/user/.libephemeris/spk/2060_2458849_2462502.bsp
        >>> # Now calc_ut(jd, SE_CHIRON, ...) uses SPK data automatically
    """
    # Check if astroquery is available
    if not _check_astroquery_available():
        raise ImportError(
            "astroquery is required for automatic SPK downloads. "
            "Install it with: pip install astroquery"
        )

    # Validate inputs
    if jd_end <= jd_start:
        raise ValueError(
            f"jd_end ({jd_end}) must be greater than jd_start ({jd_start})"
        )

    # Import state to get configuration settings
    from . import state

    # Determine cache directory (parameter > global config > default)
    if cache_dir is None:
        cache_dir = state.get_spk_cache_dir()
    if cache_dir is None:
        cache_dir = DEFAULT_AUTO_SPK_DIR

    # Apply date padding from global configuration
    date_padding = state.get_spk_date_padding()
    if date_padding > 0:
        jd_start = jd_start - date_padding
        jd_end = jd_end + date_padding

    # Create cache directory if it doesn't exist
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir, exist_ok=True)

    # Check if we already have a cached SPK that covers this range
    existing_spk = _find_covering_spk(body_id, jd_start, jd_end, cache_dir)
    if existing_spk is not None:
        # Even if cached, register if ipl is provided
        if ipl is not None:
            _register_spk_after_download(existing_spk, body_id, ipl, naif_id)
        return existing_spk

    # Need to download - convert Julian Days to date strings
    start_date = _jd_to_iso_date(jd_start)
    end_date = _jd_to_iso_date(jd_end)

    # Generate cache filename
    filename = _generate_spk_cache_filename(body_id, jd_start, jd_end)
    output_path = os.path.join(cache_dir, filename)

    # Download using astroquery
    with _AUTO_SPK_LOCK:
        # Double-check after acquiring lock (another thread may have downloaded)
        if os.path.exists(output_path):
            # Register if ipl is provided
            if ipl is not None:
                _register_spk_after_download(output_path, body_id, ipl, naif_id)
            return output_path

        # Check again for covering SPK (may have been downloaded with different range)
        existing_spk = _find_covering_spk(body_id, jd_start, jd_end, cache_dir)
        if existing_spk is not None:
            # Register if ipl is provided
            if ipl is not None:
                _register_spk_after_download(existing_spk, body_id, ipl, naif_id)
            return existing_spk

        # Download SPK file
        _download_spk_astroquery(
            body_id=str(body_id),
            start=start_date,
            end=end_date,
            output_path=output_path,
        )

    # Register the SPK body if ipl is provided
    if ipl is not None:
        _register_spk_after_download(output_path, body_id, ipl, naif_id)

    return output_path


def _register_spk_after_download(
    spk_path: str,
    body_id: Union[int, str],
    ipl: int,
    naif_id: Optional[int] = None,
) -> None:
    """
    Register an SPK body after download.

    Internal helper function that registers the SPK body with libephemeris
    so it can be used by calc_ut().

    Args:
        spk_path: Path to the downloaded SPK file
        body_id: JPL Horizons body identifier (used for NAIF ID deduction)
        ipl: libephemeris body ID (e.g., SE_CHIRON)
        naif_id: Optional NAIF ID. If not provided, deduced from body_id.

    Raises:
        ValueError: If naif_id cannot be deduced and is not provided
    """
    from . import spk

    # Deduce NAIF ID if not provided
    if naif_id is None:
        naif_id = spk._deduce_naif_id(str(body_id))
        if naif_id is None:
            raise ValueError(
                f"Cannot deduce NAIF ID for '{body_id}'. "
                f"Please provide naif_id explicitly. "
                f"For numbered asteroids: naif_id = asteroid_number + 2000000"
            )

    # Register the SPK body
    spk.register_spk_body(ipl, spk_path, naif_id)


# =============================================================================
# JPL HORIZONS SPK DOWNLOAD (PUBLIC API)
# =============================================================================


def download_spk_from_horizons(
    body_id: Union[int, str],
    jd_start: float,
    jd_end: float,
    output_path: str,
    location: str = "@0",
    ipl: Optional[int] = None,
    naif_id: Optional[int] = None,
) -> str:
    """
    Download an SPK file from JPL Horizons for a specified body and date range.

    This function uses astroquery.jplhorizons to download an SPK (SPICE kernel)
    file for the specified body and Julian Day range. The file can then be used
    for high-precision ephemeris calculations.

    If ``ipl`` is provided, the SPK body is automatically registered with
    libephemeris after download, making it immediately usable by ``calc_ut()``.

    Args:
        body_id: JPL Horizons body identifier. Can be:
            - Asteroid number (int or str): 2060, "2060", "136199"
            - Name (str): "Chiron", "Eris"
            - Designation (str): "2003 UB313"
        jd_start: Start of the date range (Julian Day)
        jd_end: End of the date range (Julian Day)
        output_path: Full path where the SPK file should be saved (including filename)
        location: Observer location code for Horizons (default: "@0" = solar system barycenter)
        ipl: Optional libephemeris body ID (e.g., SE_CHIRON). If provided,
            the SPK body is automatically registered after download.
        naif_id: Optional NAIF ID for SPK lookup. If not provided but ipl is,
            the NAIF ID is deduced from body_id using the convention:
            naif_id = asteroid_number + 2000000

    Returns:
        str: Path to the downloaded SPK file (same as output_path)

    Raises:
        ImportError: If astroquery is not installed
        ValueError: If body is not found on Horizons, date range is invalid,
            the date range is too large for Horizons to process, or naif_id
            cannot be deduced when ipl is provided
        ConnectionError: If network request fails

    Example:
        >>> from libephemeris.spk_auto import download_spk_from_horizons
        >>> from libephemeris.constants import SE_CHIRON
        >>> # Download Chiron SPK for 2020-2030 and register it
        >>> jd_start = 2458849.5  # 2020-01-01
        >>> jd_end = 2462502.5    # 2030-01-01
        >>> path = download_spk_from_horizons(
        ...     body_id="2060",
        ...     jd_start=jd_start,
        ...     jd_end=jd_end,
        ...     output_path="/path/to/chiron.bsp",
        ...     ipl=SE_CHIRON,  # Automatically registers the SPK
        ... )
        >>> print(path)
        /path/to/chiron.bsp
        >>> # Now calc_ut(jd, SE_CHIRON, ...) uses SPK data automatically

    Notes:
        - Horizons has limits on the date range for SPK generation. For most
          bodies, a 100-year range is typically acceptable.
        - For numbered asteroids, use the asteroid number (e.g., "2060" for Chiron).
        - For major planets, use the body name (e.g., "Mars", "Jupiter").
    """
    # Check if astroquery is available
    if not _check_astroquery_available():
        raise ImportError(
            "astroquery is required for downloading SPK files from Horizons. "
            "Install it with: pip install astroquery"
        )

    # Validate inputs
    if jd_end <= jd_start:
        raise ValueError(
            f"jd_end ({jd_end}) must be greater than jd_start ({jd_start})"
        )

    # Convert Julian Days to ISO date strings
    start_date = _jd_to_iso_date(jd_start)
    end_date = _jd_to_iso_date(jd_end)

    # Ensure output directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Import astroquery here to handle import errors gracefully
    try:
        from astroquery.jplhorizons import Horizons
    except ImportError as e:
        raise ImportError(
            "astroquery is required for downloading SPK files from Horizons. "
            "Install it with: pip install astroquery"
        ) from e

    # Create Horizons query object
    body_id_str = str(body_id)
    try:
        obj = Horizons(
            id=body_id_str,
            location=location,
            epochs={"start": start_date, "stop": end_date},
        )
    except Exception as e:
        raise ValueError(f"Failed to create Horizons query for '{body_id}': {e}") from e

    # Download SPK file with error handling
    try:
        obj.download_spk(output_path)
    except Exception as e:
        error_msg = str(e).lower()

        # Handle body not found
        if "unknown target" in error_msg or "no matches found" in error_msg:
            raise ValueError(
                f"Body '{body_id}' not found on JPL Horizons. "
                "Check the body identifier and try again."
            ) from e

        # Handle date range too large
        if (
            "time span" in error_msg
            or "too large" in error_msg
            or "exceeds" in error_msg
        ):
            raise ValueError(
                f"Date range too large for JPL Horizons. "
                f"Requested range: JD {jd_start} to {jd_end} "
                f"({start_date} to {end_date}). "
                "Try a smaller date range."
            ) from e

        # Handle network failures
        if (
            "connection" in error_msg
            or "timeout" in error_msg
            or "network" in error_msg
            or "urlopen" in error_msg
            or "http" in error_msg
        ):
            raise ConnectionError(
                f"Network error downloading SPK from Horizons: {e}"
            ) from e

        # Re-raise with more context for other errors
        raise ValueError(
            f"Failed to download SPK for '{body_id}' from JPL Horizons: {e}"
        ) from e

    # Register the SPK body if ipl is provided
    if ipl is not None:
        _register_spk_after_download(output_path, body_id, ipl, naif_id)

    return output_path
