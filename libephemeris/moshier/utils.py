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
