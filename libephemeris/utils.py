"""
Utility functions for libephemeris.

Provides helper functions compatible with pyswisseph API including
angular calculations and other mathematical utilities.
"""

import math
from typing import Tuple


def cotrans_sp(
    coord: Tuple[float, float, float],
    speed: Tuple[float, float, float],
    obliquity: float,
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """
    Transform coordinates and velocities between ecliptic and equatorial systems.

    This function extends cotrans() to also transform velocity (speed) components,
    using analytical derivatives of the coordinate transformation equations.

    The direction of transformation depends on the sign of obliquity:
    - Positive obliquity: ecliptic (lon, lat) -> equatorial (RA, Dec)
    - Negative obliquity: equatorial (RA, Dec) -> ecliptic (lon, lat)

    Args:
        coord: Tuple of (longitude/RA, latitude/Dec, distance) in degrees
        speed: Tuple of (lon_speed, lat_speed, dist_speed) in degrees/day
        obliquity: Obliquity of the ecliptic in degrees.
                   Positive for ecliptic->equatorial, negative for equatorial->ecliptic.

    Returns:
        Tuple of (coord_transformed, speed_transformed) where:
        - coord_transformed: (new_lon/RA, new_lat/Dec, distance)
        - speed_transformed: (new_lon_speed, new_lat_speed, dist_speed)
        Distance and distance speed are unchanged by the transformation.

    Examples:
        >>> # Ecliptic to equatorial (positive obliquity)
        >>> coord, speed = cotrans_sp((90.0, 0.0, 1.0), (1.0, 0.0, 0.0), 23.4)
        >>> # Returns transformed position and velocity
    """
    lon = coord[0]
    lat = coord[1]
    dist = coord[2]
    lon_speed = speed[0]
    lat_speed = speed[1]
    dist_speed = speed[2]

    # Convert to radians
    lon_rad = math.radians(lon)
    lat_rad = math.radians(lat)
    eps_rad = math.radians(obliquity)

    # Precompute trig values
    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # Calculate the new latitude (Dec for ecl->eq, beta for eq->ecl)
    # sin(new_lat) = sin(lat) * cos(eps) - cos(lat) * sin(eps) * sin(lon)
    sin_new_lat = sin_lat * cos_eps - cos_lat * sin_eps * sin_lon
    # Clamp to avoid numerical issues with asin
    sin_new_lat = max(-1.0, min(1.0, sin_new_lat))
    new_lat_rad = math.asin(sin_new_lat)
    cos_new_lat = math.cos(new_lat_rad)

    # Calculate the new longitude (RA for ecl->eq, lambda for eq->ecl)
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

    # --- Velocity transformation ---
    # Convert speed to radians/day for the calculation
    lon_speed_rad = math.radians(lon_speed)
    lat_speed_rad = math.radians(lat_speed)

    # Derivative of new latitude:
    # d/dt[sin(new_lat)] = cos(new_lat) * d(new_lat)/dt
    # d/dt[sin(lat)*cos(eps) - cos(lat)*sin(eps)*sin(lon)]
    #   = cos(lat)*cos(eps)*d(lat)/dt + sin(lat)*sin(eps)*sin(lon)*d(lat)/dt
    #     - cos(lat)*sin(eps)*cos(lon)*d(lon)/dt
    #
    # new_lat_speed = (1/cos(new_lat)) * [
    #   (cos(lat)*cos(eps) + sin(lat)*sin(eps)*sin(lon)) * lat_speed
    #   - cos(lat)*sin(eps)*cos(lon) * lon_speed
    # ]
    if abs(cos_new_lat) > 1e-10:
        new_lat_speed_rad = (
            (cos_lat * cos_eps + sin_lat * sin_eps * sin_lon) * lat_speed_rad
            - cos_lat * sin_eps * cos_lon * lon_speed_rad
        ) / cos_new_lat
    else:
        # At poles, latitude speed is undefined; use 0
        new_lat_speed_rad = 0.0

    # Derivative of new longitude:
    # new_lon = atan2(y, x) where y = sin(lon)*cos(eps) + tan(lat)*sin(eps), x = cos(lon)
    # d(new_lon)/dt = (x * dy/dt - y * dx/dt) / (x^2 + y^2)
    #
    # dx/dt = -sin(lon) * d(lon)/dt
    # dy/dt = cos(lon)*cos(eps)*d(lon)/dt + sec^2(lat)*sin(eps)*d(lat)/dt
    #       = cos(lon)*cos(eps)*d(lon)/dt + sin(eps)/(cos^2(lat))*d(lat)/dt
    dx_dt = -sin_lon * lon_speed_rad
    cos_lat_sq = cos_lat * cos_lat
    if abs(cos_lat_sq) > 1e-10:
        dy_dt = (
            cos_lon * cos_eps * lon_speed_rad + (sin_eps / cos_lat_sq) * lat_speed_rad
        )
    else:
        # At poles of input coordinates
        dy_dt = cos_lon * cos_eps * lon_speed_rad

    denom = x * x + y * y
    if abs(denom) > 1e-10:
        new_lon_speed_rad = (x * dy_dt - y * dx_dt) / denom
    else:
        new_lon_speed_rad = 0.0

    # Convert speeds back to degrees/day
    new_lon_speed = math.degrees(new_lon_speed_rad)
    new_lat_speed = math.degrees(new_lat_speed_rad)

    return (
        (new_lon, new_lat, dist),
        (new_lon_speed, new_lat_speed, dist_speed),
    )


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
