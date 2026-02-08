"""
Utility functions for Moshier ephemeris calculations.

Provides coordinate conversions, aberration corrections, and other mathematical
utilities needed for semi-analytical ephemeris computations.

This module has NO dependencies on Skyfield or SPK files.
Only numpy is used for numerical operations.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np
from numpy.typing import NDArray

# =============================================================================
# CONSTANTS
# =============================================================================

# Speed of light in AU/day
C_LIGHT_AU_DAY: float = 173.1446326846693

# Astronomical constants
AU_KM: float = 149597870.7  # 1 AU in km
EARTH_RADIUS_KM: float = 6378.137  # Earth equatorial radius in km

# Julian day constants
J2000: float = 2451545.0  # JD of J2000.0 epoch
JD_PER_CENTURY: float = 36525.0  # Julian days per Julian century

# Pi constants
TWO_PI: float = 2.0 * math.pi
DEG_TO_RAD: float = math.pi / 180.0
RAD_TO_DEG: float = 180.0 / math.pi
ARCSEC_TO_RAD: float = math.pi / (180.0 * 3600.0)
RAD_TO_ARCSEC: float = (180.0 * 3600.0) / math.pi


# =============================================================================
# ANGLE NORMALIZATION
# =============================================================================


def normalize_angle(angle: float) -> float:
    """Normalize angle to range [0, 360) degrees.

    Args:
        angle: Angle in degrees.

    Returns:
        Normalized angle in [0, 360) degrees.
    """
    angle = angle % 360.0
    if angle < 0:
        angle += 360.0
    return angle


def normalize_radians(angle: float) -> float:
    """Normalize angle to range [0, 2*pi) radians.

    Args:
        angle: Angle in radians.

    Returns:
        Normalized angle in [0, 2*pi) radians.
    """
    angle = angle % TWO_PI
    if angle < 0:
        angle += TWO_PI
    return angle


# =============================================================================
# COORDINATE CONVERSIONS
# =============================================================================


def spherical_to_cartesian(
    lon: float, lat: float, dist: float
) -> Tuple[float, float, float]:
    """Convert spherical coordinates to Cartesian.

    Args:
        lon: Longitude in degrees.
        lat: Latitude in degrees.
        dist: Distance (radius).

    Returns:
        Tuple of (x, y, z) Cartesian coordinates.
    """
    lon_rad = lon * DEG_TO_RAD
    lat_rad = lat * DEG_TO_RAD

    cos_lat = math.cos(lat_rad)
    x = dist * cos_lat * math.cos(lon_rad)
    y = dist * cos_lat * math.sin(lon_rad)
    z = dist * math.sin(lat_rad)

    return x, y, z


def cartesian_to_spherical(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Convert Cartesian coordinates to spherical.

    Args:
        x: X coordinate.
        y: Y coordinate.
        z: Z coordinate.

    Returns:
        Tuple of (longitude, latitude, distance) where angles are in degrees.
    """
    dist = math.sqrt(x * x + y * y + z * z)

    if dist < 1e-30:
        return 0.0, 0.0, 0.0

    lon = math.atan2(y, x) * RAD_TO_DEG
    lat = math.asin(z / dist) * RAD_TO_DEG

    # Normalize longitude to [0, 360)
    lon = normalize_angle(lon)

    return lon, lat, dist


def ecliptic_to_equatorial(
    lon: float,
    lat: float,
    obliquity: float,
) -> Tuple[float, float]:
    """Convert ecliptic coordinates to equatorial.

    Args:
        lon: Ecliptic longitude in degrees.
        lat: Ecliptic latitude in degrees.
        obliquity: Obliquity of the ecliptic in degrees.

    Returns:
        Tuple of (right_ascension, declination) in degrees.
    """
    lon_rad = lon * DEG_TO_RAD
    lat_rad = lat * DEG_TO_RAD
    eps_rad = obliquity * DEG_TO_RAD

    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # sin(dec) = sin(lat)*cos(eps) + cos(lat)*sin(eps)*sin(lon)
    sin_dec = sin_lat * cos_eps + cos_lat * sin_eps * sin_lon
    sin_dec = max(-1.0, min(1.0, sin_dec))  # Clamp for safety
    dec = math.asin(sin_dec)

    # cos(ra) = cos(lat)*cos(lon) / cos(dec)
    # sin(ra) = [cos(lat)*sin(lon)*cos(eps) - sin(lat)*sin(eps)] / cos(dec)
    cos_dec = math.cos(dec)
    if abs(cos_dec) > 1e-10:
        x = cos_lat * cos_lon
        y = cos_lat * sin_lon * cos_eps - sin_lat * sin_eps
        ra = math.atan2(y, x)
    else:
        ra = 0.0

    ra_deg = normalize_angle(ra * RAD_TO_DEG)
    dec_deg = dec * RAD_TO_DEG

    return ra_deg, dec_deg


def equatorial_to_ecliptic(
    ra: float,
    dec: float,
    obliquity: float,
) -> Tuple[float, float]:
    """Convert equatorial coordinates to ecliptic.

    Args:
        ra: Right ascension in degrees.
        dec: Declination in degrees.
        obliquity: Obliquity of the ecliptic in degrees.

    Returns:
        Tuple of (ecliptic_longitude, ecliptic_latitude) in degrees.
    """
    ra_rad = ra * DEG_TO_RAD
    dec_rad = dec * DEG_TO_RAD
    eps_rad = obliquity * DEG_TO_RAD

    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_dec = math.cos(dec_rad)
    sin_dec = math.sin(dec_rad)
    cos_ra = math.cos(ra_rad)
    sin_ra = math.sin(ra_rad)

    # sin(lat) = sin(dec)*cos(eps) - cos(dec)*sin(eps)*sin(ra)
    sin_lat = sin_dec * cos_eps - cos_dec * sin_eps * sin_ra
    sin_lat = max(-1.0, min(1.0, sin_lat))
    lat = math.asin(sin_lat)

    # cos(lon) = cos(dec)*cos(ra) / cos(lat)
    # sin(lon) = [cos(dec)*sin(ra)*cos(eps) + sin(dec)*sin(eps)] / cos(lat)
    cos_lat = math.cos(lat)
    if abs(cos_lat) > 1e-10:
        x = cos_dec * cos_ra
        y = cos_dec * sin_ra * cos_eps + sin_dec * sin_eps
        lon = math.atan2(y, x)
    else:
        lon = 0.0

    lon_deg = normalize_angle(lon * RAD_TO_DEG)
    lat_deg = lat * RAD_TO_DEG

    return lon_deg, lat_deg


def rectangular_to_spherical_velocity(
    x: float,
    y: float,
    z: float,
    vx: float,
    vy: float,
    vz: float,
) -> Tuple[float, float, float, float, float, float]:
    """Convert rectangular position/velocity to spherical.

    Args:
        x, y, z: Position in rectangular coordinates (AU).
        vx, vy, vz: Velocity in rectangular coordinates (AU/day).

    Returns:
        Tuple of (lon, lat, dist, dlon, dlat, ddist) where:
        - lon, lat in degrees
        - dist in AU
        - dlon, dlat in degrees/day
        - ddist in AU/day
    """
    # Position
    r = math.sqrt(x * x + y * y + z * z)
    rho = math.sqrt(x * x + y * y)

    if r < 1e-30:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    lon = math.atan2(y, x)
    lat = math.asin(z / r)

    # Velocity derivatives
    # dr/dt = (x*vx + y*vy + z*vz) / r
    dr_dt = (x * vx + y * vy + z * vz) / r

    # drho/dt = (x*vx + y*vy) / rho
    if rho > 1e-10:
        # dlon/dt = (x*vy - y*vx) / rho^2
        dlon_dt = (x * vy - y * vx) / (rho * rho)

        # dlat/dt = (vz*r - z*dr_dt) / (r^2 * sqrt(1 - z^2/r^2))
        # Simplified: dlat/dt = (vz - z*dr_dt/r) / rho
        # Actually: d/dt(asin(z/r)) = (vz*r - z*dr_dt) / (r^2 * rho / r)
        #                          = (vz*r - z*dr_dt) / (r * rho)
        dlat_dt = (vz * r - z * dr_dt) / (r * rho)
    else:
        dlon_dt = 0.0
        dlat_dt = 0.0

    lon_deg = normalize_angle(lon * RAD_TO_DEG)
    lat_deg = lat * RAD_TO_DEG
    dlon_deg = dlon_dt * RAD_TO_DEG
    dlat_deg = dlat_dt * RAD_TO_DEG

    return lon_deg, lat_deg, r, dlon_deg, dlat_deg, dr_dt


# =============================================================================
# ABERRATION
# =============================================================================


def annual_aberration(
    lon: float,
    lat: float,
    earth_x: float,
    earth_y: float,
    earth_z: float,
    earth_vx: float,
    earth_vy: float,
    earth_vz: float,
) -> Tuple[float, float]:
    """Calculate annual aberration correction.

    The annual aberration is the apparent displacement of a celestial body
    due to the finite speed of light and the motion of the Earth.

    Args:
        lon: Ecliptic longitude in degrees.
        lat: Ecliptic latitude in degrees.
        earth_x, earth_y, earth_z: Earth heliocentric position (AU).
        earth_vx, earth_vy, earth_vz: Earth heliocentric velocity (AU/day).

    Returns:
        Tuple of (delta_lon, delta_lat) corrections in degrees.
    """
    # Convert to radians
    lon_rad = lon * DEG_TO_RAD
    lat_rad = lat * DEG_TO_RAD

    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)

    # Aberration formula (first-order approximation)
    # delta_lon = (-vx*sin(lon) + vy*cos(lon)) / (c * cos(lat))
    # delta_lat = (-vx*cos(lon)*sin(lat) - vy*sin(lon)*sin(lat) + vz*cos(lat)) / c

    if abs(cos_lat) > 1e-10:
        delta_lon = (-earth_vx * sin_lon + earth_vy * cos_lon) / (
            C_LIGHT_AU_DAY * cos_lat
        )
    else:
        delta_lon = 0.0

    delta_lat = (
        -earth_vx * cos_lon * sin_lat
        - earth_vy * sin_lon * sin_lat
        + earth_vz * cos_lat
    ) / C_LIGHT_AU_DAY

    return delta_lon * RAD_TO_DEG, delta_lat * RAD_TO_DEG


# =============================================================================
# LIGHT-TIME CORRECTION
# =============================================================================


def light_time_correction(distance_au: float) -> float:
    """Calculate light-time delay.

    Args:
        distance_au: Distance in AU.

    Returns:
        Light-time delay in days.
    """
    return distance_au / C_LIGHT_AU_DAY


def iterative_light_time_correction(
    jd_tt: float,
    earth_pos: Tuple[float, float, float],
    target_func,
    iterations: int = 3,
) -> Tuple[Tuple[float, float, float], float]:
    """Calculate light-time corrected position using iterative method.

    Light travels at a finite speed (~8 minutes/AU). When we observe a planet
    at time t, we see it where it was at time t - Δt, where Δt is the light
    travel time. This function iterates to find the correct antedated position.

    The Swiss Ephemeris Moshier algorithm uses 3 iterations, which achieves
    convergence to ~10^-10 AU accuracy for all solar system objects.

    Args:
        jd_tt: Julian Day in Terrestrial Time (observation time).
        earth_pos: Observer position as (x, y, z) in AU (heliocentric ecliptic).
        target_func: Function(jd_tt) -> (x, y, z) returning target heliocentric
                    position in AU at given time.
        iterations: Number of light-time iterations (default 3, like Swiss Eph).

    Returns:
        Tuple containing:
        - (x, y, z): Light-time corrected geocentric position in AU.
        - light_time: Light travel time in days.

    Example:
        >>> def mars_helio(jd):
        ...     return calc_vsop87_heliocentric_cartesian(jd, MARS)
        >>> corrected_pos, lt = iterative_light_time_correction(
        ...     jd_tt, earth_pos, mars_helio
        ... )
    """
    # Initial estimate: target position at observation time
    target_pos = target_func(jd_tt)

    # Geocentric vector (initial)
    dx = target_pos[0] - earth_pos[0]
    dy = target_pos[1] - earth_pos[1]
    dz = target_pos[2] - earth_pos[2]
    distance = math.sqrt(dx * dx + dy * dy + dz * dz)

    # Light-time in days
    light_time = distance / C_LIGHT_AU_DAY

    # Iterate to refine
    for _ in range(iterations):
        # Antedated time
        t_corrected = jd_tt - light_time

        # Target position at corrected time
        target_pos = target_func(t_corrected)

        # Recompute geocentric vector
        dx = target_pos[0] - earth_pos[0]
        dy = target_pos[1] - earth_pos[1]
        dz = target_pos[2] - earth_pos[2]
        distance = math.sqrt(dx * dx + dy * dy + dz * dz)

        # Update light-time
        light_time = distance / C_LIGHT_AU_DAY

    return (dx, dy, dz), light_time


def annual_aberration_cartesian(
    target_direction: Tuple[float, float, float],
    earth_velocity: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """Calculate stellar aberration correction using Cartesian vectors.

    Annual aberration is caused by Earth's orbital velocity (~30 km/s) combined
    with the finite speed of light. The apparent position of a star or planet
    is shifted in the direction of Earth's motion by up to κ ≈ 20.49552 arcsec.

    This function uses the exact vector formula (relativistic to first order):
        Δr = (v/c) - (r·v/c) * r

    where r is the unit direction vector to the target and v is Earth's velocity.

    Args:
        target_direction: Unit direction vector (x, y, z) toward target (geocentric).
        earth_velocity: Earth heliocentric velocity (vx, vy, vz) in AU/day.

    Returns:
        Aberration correction (dx, dy, dz) to be ADDED to the apparent position.
        Units are dimensionless (fraction of unit vector).

    Note:
        The correction is typically ~20 arcsec maximum (~10^-4 radians).
        To apply: apparent_direction = geometric_direction + correction
        Then renormalize to unit length.

    References:
        - Explanatory Supplement to the Astronomical Almanac (3rd ed), Section 7.2
        - Urban & Seidelmann (2013), "Aberration"
    """
    # Normalize target direction (should already be normalized, but ensure)
    rx, ry, rz = target_direction
    r_mag = math.sqrt(rx * rx + ry * ry + rz * rz)
    if r_mag < 1e-30:
        return (0.0, 0.0, 0.0)
    rx /= r_mag
    ry /= r_mag
    rz /= r_mag

    # Earth velocity in units of c (dimensionless)
    vx = earth_velocity[0] / C_LIGHT_AU_DAY
    vy = earth_velocity[1] / C_LIGHT_AU_DAY
    vz = earth_velocity[2] / C_LIGHT_AU_DAY

    # Dot product: r · v/c
    r_dot_v = rx * vx + ry * vy + rz * vz

    # Aberration formula: Δr = v/c - (r · v/c) * r
    # This is the first-order (classical) Bradley aberration formula
    dx = vx - r_dot_v * rx
    dy = vy - r_dot_v * ry
    dz = vz - r_dot_v * rz

    return (dx, dy, dz)


def apply_aberration_to_position(
    position: Tuple[float, float, float],
    earth_velocity: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """Apply stellar aberration correction to a geocentric position.

    Convenience function that applies aberration to a position vector,
    returning the apparent (aberrated) position.

    Args:
        position: Geocentric position (x, y, z) in AU.
        earth_velocity: Earth heliocentric velocity (vx, vy, vz) in AU/day.

    Returns:
        Apparent (aberrated) position (x, y, z) in AU.
        The distance is preserved; only direction is affected.
    """
    # Get distance
    x, y, z = position
    distance = math.sqrt(x * x + y * y + z * z)
    if distance < 1e-30:
        return position

    # Unit direction
    ux, uy, uz = x / distance, y / distance, z / distance

    # Get aberration correction
    dx, dy, dz = annual_aberration_cartesian((ux, uy, uz), earth_velocity)

    # Apply correction to unit direction
    ax = ux + dx
    ay = uy + dy
    az = uz + dz

    # Renormalize
    a_mag = math.sqrt(ax * ax + ay * ay + az * az)
    if a_mag < 1e-30:
        return position
    ax /= a_mag
    ay /= a_mag
    az /= a_mag

    # Scale back to original distance
    return (ax * distance, ay * distance, az * distance)


def aberration_in_longitude_latitude(
    lon: float,
    lat: float,
    sun_lon: float,
    obliquity: float,
) -> Tuple[float, float]:
    """Calculate aberration correction in ecliptic longitude and latitude.

    This is the classical Bradley aberration formula for ecliptic coordinates,
    using the approximation that Earth's velocity is directed along the ecliptic
    at longitude (sun_lon + 90°).

    The constant of aberration κ = 20.49552 arcseconds.

    Args:
        lon: Ecliptic longitude in degrees.
        lat: Ecliptic latitude in degrees.
        sun_lon: Apparent Sun longitude in degrees (geocentric).
        obliquity: Obliquity of the ecliptic in degrees.

    Returns:
        Tuple of (delta_lon, delta_lat) corrections in degrees.
        Add these to the geometric position to get apparent position.

    Note:
        This is a simplified formula assuming circular Earth orbit.
        For higher precision, use annual_aberration_cartesian with
        actual Earth velocity vectors.
    """
    # Constant of aberration in degrees
    KAPPA = 20.49552 / 3600.0  # arcsec -> degrees

    # Convert to radians
    lon_rad = lon * DEG_TO_RAD
    lat_rad = lat * DEG_TO_RAD
    sun_rad = sun_lon * DEG_TO_RAD
    eps_rad = obliquity * DEG_TO_RAD

    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_sun = math.cos(sun_rad)
    sin_sun = math.sin(sun_rad)
    cos_eps = math.cos(eps_rad)

    # Bradley aberration formulas (simplified, assuming circular orbit)
    # Δλ = -κ · (cos(λ☉) · cos(λ) + sin(λ☉) · sin(λ) · cos(ε)) / cos(β)
    # Δβ = -κ · sin(β) · (sin(λ☉) · cos(λ) · cos(ε) - cos(λ☉) · sin(λ))
    #      - κ · cos(β) · sin(λ☉) · sin(ε)
    # However, a more standard form:
    # Δλ = κ · (-cos(λ☉ - λ) + sin(λ☉) · sin(λ) · (1 - cos(ε))) / cos(β)
    # For simplicity, use the vector approach which is more accurate.

    # Simplified: velocity direction is perpendicular to Sun direction
    # v_direction = sun_lon + 90°
    v_lon_rad = sun_rad + math.pi / 2

    # Δλ = (κ / cos(β)) · cos(v_lon - λ)
    if abs(cos_lat) > 1e-10:
        delta_lon = KAPPA * math.cos(v_lon_rad - lon_rad) / cos_lat
    else:
        delta_lon = 0.0

    # Δβ = κ · sin(β) · sin(v_lon - λ)
    delta_lat = KAPPA * sin_lat * math.sin(v_lon_rad - lon_rad)

    return delta_lon, delta_lat


# =============================================================================
# CARTESIAN COORDINATE UTILITIES
# =============================================================================


def spherical_to_cartesian_velocity(
    lon: float,
    lat: float,
    dist: float,
    dlon: float,
    dlat: float,
    ddist: float,
) -> Tuple[float, float, float, float, float, float]:
    """Convert spherical coordinates and velocities to Cartesian.

    Args:
        lon: Longitude in degrees.
        lat: Latitude in degrees.
        dist: Distance in AU.
        dlon: Longitude velocity in degrees/day.
        dlat: Latitude velocity in degrees/day.
        ddist: Distance velocity in AU/day.

    Returns:
        Tuple of (x, y, z, vx, vy, vz) in AU and AU/day.
    """
    lon_rad = lon * DEG_TO_RAD
    lat_rad = lat * DEG_TO_RAD
    dlon_rad = dlon * DEG_TO_RAD  # rad/day
    dlat_rad = dlat * DEG_TO_RAD  # rad/day

    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)

    # Position
    x = dist * cos_lat * cos_lon
    y = dist * cos_lat * sin_lon
    z = dist * sin_lat

    # Velocity (chain rule derivatives)
    # dx/dt = cos(lat)·cos(lon)·dr/dt - r·cos(lat)·sin(lon)·dlon/dt
    #         - r·sin(lat)·cos(lon)·dlat/dt
    vx = (
        cos_lat * cos_lon * ddist
        - dist * cos_lat * sin_lon * dlon_rad
        - dist * sin_lat * cos_lon * dlat_rad
    )
    vy = (
        cos_lat * sin_lon * ddist
        + dist * cos_lat * cos_lon * dlon_rad
        - dist * sin_lat * sin_lon * dlat_rad
    )
    vz = sin_lat * ddist + dist * cos_lat * dlat_rad

    return x, y, z, vx, vy, vz


def compute_earth_velocity(
    jd_tt: float,
    h: float = 0.001,
) -> Tuple[float, float, float]:
    """Compute Earth's heliocentric velocity by numerical differentiation.

    Uses central differences on the VSOP87 Earth position.

    Args:
        jd_tt: Julian Day in Terrestrial Time.
        h: Step size in days (default 0.001 = ~1.44 minutes).

    Returns:
        Tuple of (vx, vy, vz) velocity in AU/day (heliocentric ecliptic).

    Note:
        This function requires the vsop87 module. Import it locally to
        avoid circular imports.
    """
    # Use lazy import to avoid circular dependency
    from .vsop87 import calc_earth_heliocentric

    # Get positions at t±h
    lon_plus, lat_plus, r_plus = calc_earth_heliocentric(jd_tt + h)
    lon_minus, lat_minus, r_minus = calc_earth_heliocentric(jd_tt - h)

    # Convert to Cartesian
    x_plus, y_plus, z_plus = spherical_to_cartesian(lon_plus, lat_plus, r_plus)
    x_minus, y_minus, z_minus = spherical_to_cartesian(lon_minus, lat_minus, r_minus)

    # Central differences
    vx = (x_plus - x_minus) / (2 * h)
    vy = (y_plus - y_minus) / (2 * h)
    vz = (z_plus - z_minus) / (2 * h)

    return vx, vy, vz


def apply_light_time_and_aberration(
    jd_tt: float,
    geometric_pos: Tuple[float, float, float],
    earth_pos: Tuple[float, float, float],
    earth_velocity: Tuple[float, float, float],
    target_func=None,
    iterations: int = 3,
) -> Tuple[Tuple[float, float, float], float]:
    """Apply both light-time and aberration corrections.

    This is the main entry point for computing apparent positions in the
    Moshier ephemeris. It combines:
    1. Light-time correction (iterative)
    2. Annual stellar aberration

    Args:
        jd_tt: Julian Day in Terrestrial Time (observation time).
        geometric_pos: Uncorrected geocentric position (x, y, z) in AU.
        earth_pos: Earth heliocentric position (x, y, z) in AU.
        earth_velocity: Earth heliocentric velocity (vx, vy, vz) in AU/day.
        target_func: Optional function(jd) -> (x, y, z) for light-time iteration.
                    If None, only simple light-time delay is applied.
        iterations: Number of light-time iterations (default 3).

    Returns:
        Tuple containing:
        - (x, y, z): Apparent position with all corrections applied.
        - light_time: Light travel time in days.
    """
    if target_func is not None:
        # Full iterative light-time correction
        corrected_pos, light_time = iterative_light_time_correction(
            jd_tt, earth_pos, target_func, iterations
        )
    else:
        # Simple light-time (no iteration, geometric position assumed correct)
        x, y, z = geometric_pos
        distance = math.sqrt(x * x + y * y + z * z)
        light_time = distance / C_LIGHT_AU_DAY
        corrected_pos = geometric_pos

    # Apply aberration
    apparent_pos = apply_aberration_to_position(corrected_pos, earth_velocity)

    return apparent_pos, light_time


# =============================================================================
# TIME CONVERSIONS
# =============================================================================


def jd_to_julian_centuries(jd_tt: float) -> float:
    """Convert Julian Day (TT) to Julian centuries from J2000.0.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Julian centuries from J2000.0.
    """
    return (jd_tt - J2000) / JD_PER_CENTURY


def jd_to_julian_millennia(jd_tt: float) -> float:
    """Convert Julian Day (TT) to Julian millennia from J2000.0.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Julian millennia from J2000.0.
    """
    return (jd_tt - J2000) / (JD_PER_CENTURY * 10.0)


def julian_centuries_to_jd(t: float) -> float:
    """Convert Julian centuries from J2000.0 to Julian Day.

    Args:
        t: Julian centuries from J2000.0.

    Returns:
        Julian Day (TT).
    """
    return t * JD_PER_CENTURY + J2000


# =============================================================================
# NUMERICAL DIFFERENTIATION
# =============================================================================


def numerical_derivative(
    func,
    jd_tt: float,
    body_id: int,
    h: float = 0.001,
) -> Tuple[float, float, float]:
    """Compute numerical derivatives for velocity.

    Uses central differences to compute velocity from position function.

    Args:
        func: Position function that takes (jd_tt, body_id) and returns
              (lon, lat, dist).
        jd_tt: Julian Day (TT) at which to compute derivative.
        body_id: Body identifier.
        h: Step size in days (default 0.001 = ~1.44 minutes).

    Returns:
        Tuple of (dlon_dt, dlat_dt, ddist_dt) in degrees/day and AU/day.
    """
    # Forward and backward positions
    pos_plus = func(jd_tt + h, body_id)
    pos_minus = func(jd_tt - h, body_id)

    # Central differences
    dlon = pos_plus[0] - pos_minus[0]
    dlat = pos_plus[1] - pos_minus[1]
    ddist = pos_plus[2] - pos_minus[2]

    # Handle longitude wrap-around
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0

    return dlon / (2 * h), dlat / (2 * h), ddist / (2 * h)


# =============================================================================
# POLYNOMIAL EVALUATION
# =============================================================================


def evaluate_polynomial(coeffs: Tuple[float, ...] | list[float], x: float) -> float:
    """Evaluate polynomial using Horner's method.

    Args:
        coeffs: Polynomial coefficients [a0, a1, a2, ...] for a0 + a1*x + a2*x^2 + ...
        x: Independent variable.

    Returns:
        Polynomial value.
    """
    result = 0.0
    for c in reversed(coeffs):
        result = result * x + c
    return result


def evaluate_polynomial_derivative(
    coeffs: Tuple[float, ...] | list[float], x: float
) -> float:
    """Evaluate polynomial derivative using Horner's method.

    Args:
        coeffs: Polynomial coefficients [a0, a1, a2, ...].
        x: Independent variable.

    Returns:
        Derivative value.
    """
    if len(coeffs) <= 1:
        return 0.0

    # Derivative coefficients: [a1, 2*a2, 3*a3, ...]
    deriv_coeffs = [i * c for i, c in enumerate(coeffs) if i > 0]
    return evaluate_polynomial(deriv_coeffs, x)


# =============================================================================
# TRIGONOMETRIC SERIES EVALUATION
# =============================================================================


def evaluate_trig_series(
    terms: NDArray[np.float64],
    t: float,
) -> float:
    """Evaluate a trigonometric series.

    Each term is [amplitude, phase, frequency] for:
    amplitude * cos(phase + frequency * t)

    Args:
        terms: Array of shape (n, 3) with [amplitude, phase, frequency] rows.
        t: Time variable (e.g., Julian centuries).

    Returns:
        Sum of all terms.
    """
    if len(terms) == 0:
        return 0.0

    amplitudes = terms[:, 0]
    phases = terms[:, 1]
    frequencies = terms[:, 2]

    arguments = phases + frequencies * t
    return float(np.sum(amplitudes * np.cos(arguments)))


def evaluate_trig_series_with_derivative(
    terms: NDArray[np.float64],
    t: float,
) -> Tuple[float, float]:
    """Evaluate a trigonometric series and its derivative.

    Args:
        terms: Array of shape (n, 3) with [amplitude, phase, frequency] rows.
        t: Time variable.

    Returns:
        Tuple of (value, derivative).
    """
    if len(terms) == 0:
        return 0.0, 0.0

    amplitudes = terms[:, 0]
    phases = terms[:, 1]
    frequencies = terms[:, 2]

    arguments = phases + frequencies * t

    value = float(np.sum(amplitudes * np.cos(arguments)))
    derivative = float(np.sum(-amplitudes * frequencies * np.sin(arguments)))

    return value, derivative
