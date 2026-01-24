"""
Utility functions for libephemeris.

Provides helper functions compatible with pyswisseph API including
angular calculations and other mathematical utilities.
"""

import math
from typing import Tuple


def cotrans(
    coord: Tuple[float, float, float], obliquity: float
) -> Tuple[float, float, float]:
    """
    Transform coordinates between ecliptic and equatorial systems.

    Compatible with pyswisseph's swe.cotrans() function.

    The direction of transformation depends on the sign of obliquity:
    - Positive obliquity: ecliptic (lon, lat) → equatorial (RA, Dec)
    - Negative obliquity: equatorial (RA, Dec) → ecliptic (lon, lat)

    Args:
        coord: Tuple of (longitude/RA, latitude/Dec, distance) in degrees
        obliquity: Obliquity of the ecliptic in degrees.
                   Positive for ecliptic→equatorial, negative for equatorial→ecliptic.

    Returns:
        Tuple of (transformed_lon/RA, transformed_lat/Dec, distance)
        Distance is unchanged by the transformation.

    Examples:
        >>> # Ecliptic to equatorial (positive obliquity)
        >>> cotrans((0.0, 0.0, 1.0), 23.4)
        (0.0, 0.0, 1.0)
        >>> # Equatorial to ecliptic (negative obliquity)
        >>> cotrans((0.0, 0.0, 1.0), -23.4)
        (0.0, 0.0, 1.0)
    """
    lon = coord[0]
    lat = coord[1]
    dist = coord[2]

    # Convert to radians
    lon_rad = math.radians(lon)
    lat_rad = math.radians(lat)
    eps_rad = math.radians(obliquity)

    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # Calculate the new latitude (Dec for ecl→eq, β for eq→ecl)
    # sin(new_lat) = sin(lat) * cos(eps) - cos(lat) * sin(eps) * sin(lon)
    sin_new_lat = sin_lat * cos_eps - cos_lat * sin_eps * sin_lon
    # Clamp to avoid numerical issues with asin
    sin_new_lat = max(-1.0, min(1.0, sin_new_lat))
    new_lat_rad = math.asin(sin_new_lat)

    # Calculate the new longitude (RA for ecl→eq, λ for eq→ecl)
    # tan(new_lon) = (sin(lon) * cos(eps) + tan(lat) * sin(eps)) / cos(lon)
    tan_lat = math.tan(lat_rad)
    y = sin_lon * cos_eps + tan_lat * sin_eps
    x = cos_lon

    new_lon_rad = math.atan2(y, x)

    # Convert back to degrees
    new_lon = math.degrees(new_lon_rad)
    new_lat = math.degrees(new_lat_rad)

    # Normalize longitude to [0, 360)
    new_lon = new_lon % 360.0

    return (new_lon, new_lat, dist)


def difdeg2n(p1: float, p2: float) -> float:
    """
    Calculate distance in degrees p1 - p2 normalized to [-180;180].

    Compatible with pyswisseph's swe.difdeg2n() function.
    Computes the signed angular difference, handling 360° wrapping.

    Args:
        p1: First angle in degrees
        p2: Second angle in degrees

    Returns:
        Normalized difference in range [-180, 180]

    Examples:
        >>> difdeg2n(10, 20)
        -10.0
        >>> difdeg2n(350, 10)
        -20.0
        >>> difdeg2n(10, 350)
        20.0
        >>> difdeg2n(180, 0)
        180.0
    """
    diff = (p1 - p2) % 360.0
    if diff > 180.0:
        diff -= 360.0
    return diff


def swe_calc_angles(jd_ut: float, lat: float, lon: float):
    """
    Pre-calculate and cache astrological angles and planet positions
    for use with Arabic parts.

    Args:
        jd_ut: Julian Day (UT)
        lat: Latitude (degrees)
        lon: Longitude (degrees)

    Returns:
        Dictionary with calculated positions
    """
    from .state import set_angles_cache, set_topo
    from .angles import calc_angles
    from .planets import swe_calc_ut
    from .constants import SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS

    # Set observer location
    set_topo(lon, lat, 0)

    # Calculate angles
    angles_dict = calc_angles(jd_ut, lat, lon)

    # Calculate and add planet positions for Arabic parts
    sun_pos, _ = swe_calc_ut(jd_ut, SE_SUN, 0)
    moon_pos, _ = swe_calc_ut(jd_ut, SE_MOON, 0)
    mercury_pos, _ = swe_calc_ut(jd_ut, SE_MERCURY, 0)
    venus_pos, _ = swe_calc_ut(jd_ut, SE_VENUS, 0)

    angles_dict["Sun"] = sun_pos[0]
    angles_dict["Moon"] = moon_pos[0]
    angles_dict["Mercury"] = mercury_pos[0]
    angles_dict["Venus"] = venus_pos[0]

    # Cache for Arabic parts
    set_angles_cache(angles_dict)

    return angles_dict
