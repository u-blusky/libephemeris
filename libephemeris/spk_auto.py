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
        if self.cache_dir:
            cache_path = self.cache_dir
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
