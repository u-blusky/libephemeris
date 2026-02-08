"""
ELP 2000-82B lunar theory implementation for Moshier ephemeris.

Implements a truncated version of the ELP 2000-82B lunar theory for
calculating the geocentric position of the Moon.

The ELP (Éphéméride Lunaire Parisienne) theory was developed by
M. Chapront-Touzé and J. Chapront at the Bureau des Longitudes, Paris.

This implementation uses truncated series based on Jean Meeus's
"Astronomical Algorithms" (Chapter 47), which provides arcsecond-level
accuracy sufficient for most applications.

This module has NO dependencies on Skyfield or SPK files.
Only numpy is used for numerical operations.

References:
- Chapront-Touzé, M. and Chapront, J. (1983), "The lunar ephemeris ELP 2000",
  A&A 124, 50-62
- Meeus, J. (1991), "Astronomical Algorithms", Chapter 47
"""

from __future__ import annotations

import math
from typing import Tuple

from .utils import (
    ARCSEC_TO_RAD,
    DEG_TO_RAD,
    J2000,
    JD_PER_CENTURY,
    RAD_TO_DEG,
    jd_to_julian_centuries,
    normalize_angle,
)

# =============================================================================
# CONSTANTS
# =============================================================================

# Swiss Ephemeris Moon body ID
MOSHIER_MOON = 1

# Mean distance to Moon in km
MOON_MEAN_DISTANCE_KM = 385000.56

# Earth radius in km (for parallax)
EARTH_RADIUS_KM = 6378.137

# 1 AU in km
AU_KM = 149597870.7


# =============================================================================
# FUNDAMENTAL ARGUMENTS
# =============================================================================


def _fundamental_arguments(T: float) -> Tuple[float, float, float, float, float]:
    """Calculate fundamental arguments for lunar theory.

    Args:
        T: Julian centuries from J2000.0.

    Returns:
        Tuple of (L', D, M, M', F) in degrees:
        - L': Moon's mean longitude
        - D: Mean elongation of Moon from Sun
        - M: Sun's mean anomaly
        - M': Moon's mean anomaly
        - F: Moon's argument of latitude
    """
    # Moon's mean longitude (L')
    Lp = (
        218.3164477
        + 481267.88123421 * T
        - 0.0015786 * T**2
        + T**3 / 538841.0
        - T**4 / 65194000.0
    )

    # Mean elongation of Moon from Sun (D)
    D = (
        297.8501921
        + 445267.1114034 * T
        - 0.0018819 * T**2
        + T**3 / 545868.0
        - T**4 / 113065000.0
    )

    # Sun's mean anomaly (M)
    M = 357.5291092 + 35999.0502909 * T - 0.0001536 * T**2 + T**3 / 24490000.0

    # Moon's mean anomaly (M')
    Mp = (
        134.9633964
        + 477198.8675055 * T
        + 0.0087414 * T**2
        + T**3 / 69699.0
        - T**4 / 14712000.0
    )

    # Moon's argument of latitude (F)
    F = (
        93.2720950
        + 483202.0175233 * T
        - 0.0036539 * T**2
        - T**3 / 3526000.0
        + T**4 / 863310000.0
    )

    return Lp, D, M, Mp, F


def _mean_longitude_ascending_node(T: float) -> float:
    """Calculate mean longitude of Moon's ascending node.

    Args:
        T: Julian centuries from J2000.0.

    Returns:
        Omega in degrees.
    """
    Omega = (
        125.0445479
        - 1934.1362891 * T
        + 0.0020754 * T**2
        + T**3 / 467441.0
        - T**4 / 60616000.0
    )
    return Omega


# =============================================================================
# LONGITUDE SERIES (sigma l)
# =============================================================================

# Terms for Moon's longitude
# Each term: (D, M, M', F, sin_coeff, cos_coeff)
# D, M, M', F are integer multipliers for fundamental arguments
# sin_coeff is in units of 0.000001 degrees

_LONGITUDE_TERMS = [
    # D    M   M'   F     sine coefficient
    (0, 0, 1, 0, 6288774),
    (2, 0, -1, 0, 1274027),
    (2, 0, 0, 0, 658314),
    (0, 0, 2, 0, 213618),
    (0, 1, 0, 0, -185116),
    (0, 0, 0, 2, -114332),
    (2, 0, -2, 0, 58793),
    (2, -1, -1, 0, 57066),
    (2, 0, 1, 0, 53322),
    (2, -1, 0, 0, 45758),
    (0, 1, -1, 0, -40923),
    (1, 0, 0, 0, -34720),
    (0, 1, 1, 0, -30383),
    (2, 0, 0, -2, 15327),
    (0, 0, 1, 2, -12528),
    (0, 0, 1, -2, 10980),
    (4, 0, -1, 0, 10675),
    (0, 0, 3, 0, 10034),
    (4, 0, -2, 0, 8548),
    (2, 1, -1, 0, -7888),
    (2, 1, 0, 0, -6766),
    (1, 0, -1, 0, -5163),
    (1, 1, 0, 0, 4987),
    (2, -1, 1, 0, 4036),
    (2, 0, 2, 0, 3994),
    (4, 0, 0, 0, 3861),
    (2, 0, -3, 0, 3665),
    (0, 1, -2, 0, -2689),
    (2, 0, -1, 2, -2602),
    (2, -1, -2, 0, 2390),
    (1, 0, 1, 0, -2348),
    (2, -2, 0, 0, 2236),
    (0, 1, 2, 0, -2120),
    (0, 2, 0, 0, -2069),
    (2, -2, -1, 0, 2048),
    (2, 0, 1, -2, -1773),
    (2, 0, 0, 2, -1595),
    (4, -1, -1, 0, 1215),
    (0, 0, 2, 2, -1110),
    (3, 0, -1, 0, -892),
    (2, 1, 1, 0, -810),
    (4, -1, -2, 0, 759),
    (0, 2, -1, 0, -713),
    (2, 2, -1, 0, -700),
    (2, 1, -2, 0, 691),
    (2, -1, 0, -2, 596),
    (4, 0, 1, 0, 549),
    (0, 0, 4, 0, 537),
    (4, -1, 0, 0, 520),
    (1, 0, -2, 0, -487),
    (2, 1, 0, -2, -399),
    (0, 0, 2, -2, -381),
    (1, 1, 1, 0, 351),
    (3, 0, -2, 0, -340),
    (4, 0, -3, 0, 330),
    (2, -1, 2, 0, 327),
    (0, 2, 1, 0, -323),
    (1, 1, -1, 0, 299),
    (2, 0, 3, 0, 294),
]


# =============================================================================
# LATITUDE SERIES (sigma b)
# =============================================================================

# Terms for Moon's latitude
# Each term: (D, M, M', F, sin_coeff)

_LATITUDE_TERMS = [
    # D    M   M'   F     sine coefficient
    (0, 0, 0, 1, 5128122),
    (0, 0, 1, 1, 280602),
    (0, 0, 1, -1, 277693),
    (2, 0, 0, -1, 173237),
    (2, 0, -1, 1, 55413),
    (2, 0, -1, -1, 46271),
    (2, 0, 0, 1, 32573),
    (0, 0, 2, 1, 17198),
    (2, 0, 1, -1, 9266),
    (0, 0, 2, -1, 8822),
    (2, -1, 0, -1, 8216),
    (2, 0, -2, -1, 4324),
    (2, 0, 1, 1, 4200),
    (2, 1, 0, -1, -3359),
    (2, -1, -1, 1, 2463),
    (2, -1, 0, 1, 2211),
    (2, -1, -1, -1, 2065),
    (0, 1, -1, -1, -1870),
    (4, 0, -1, -1, 1828),
    (0, 1, 0, 1, -1794),
    (0, 0, 0, 3, -1749),
    (0, 1, -1, 1, -1565),
    (1, 0, 0, 1, -1491),
    (0, 1, 1, 1, -1475),
    (0, 1, 1, -1, -1410),
    (0, 1, 0, -1, -1344),
    (1, 0, 0, -1, -1335),
    (0, 0, 3, 1, 1107),
    (4, 0, 0, -1, 1021),
    (4, 0, -1, 1, 833),
    (0, 0, 1, -3, 777),
    (4, 0, -2, 1, 671),
    (2, 0, 0, -3, 607),
    (2, 0, 2, -1, 596),
    (2, -1, 1, -1, 491),
    (2, 0, -2, 1, -451),
    (0, 0, 3, -1, 439),
    (2, 0, 2, 1, 422),
    (2, 0, -3, -1, 421),
    (2, 1, -1, 1, -366),
    (2, 1, 0, 1, -351),
    (4, 0, 0, 1, 331),
    (2, -1, 1, 1, 315),
    (2, -2, 0, -1, 302),
    (0, 0, 1, 3, -283),
    (2, 1, 1, -1, -229),
    (1, 1, 0, -1, 223),
    (1, 1, 0, 1, 223),
    (0, 1, -2, -1, -220),
    (2, 1, -1, -1, -220),
    (1, 0, 1, 1, -185),
    (2, -1, -2, -1, 181),
    (0, 1, 2, 1, -177),
    (4, 0, -2, -1, 176),
    (4, -1, -1, -1, 166),
    (1, 0, 1, -1, -164),
    (4, 0, 1, -1, 132),
    (1, 0, -1, -1, -119),
    (4, -1, 0, -1, 115),
    (2, -2, 0, 1, 107),
]


# =============================================================================
# DISTANCE SERIES (sigma r)
# =============================================================================

# Terms for Moon's distance
# Each term: (D, M, M', F, cos_coeff)

_DISTANCE_TERMS = [
    # D    M   M'   F     cosine coefficient
    (0, 0, 1, 0, -20905355),
    (2, 0, -1, 0, -3699111),
    (2, 0, 0, 0, -2955968),
    (0, 0, 2, 0, -569925),
    (0, 1, 0, 0, 48888),
    (0, 0, 0, 2, -3149),
    (2, 0, -2, 0, 246158),
    (2, -1, -1, 0, -152138),
    (2, 0, 1, 0, -170733),
    (2, -1, 0, 0, -204586),
    (0, 1, -1, 0, -129620),
    (1, 0, 0, 0, 108743),
    (0, 1, 1, 0, 104755),
    (2, 0, 0, -2, 10321),
    (0, 0, 1, 2, 0),
    (0, 0, 1, -2, 79661),
    (4, 0, -1, 0, -34782),
    (0, 0, 3, 0, -23210),
    (4, 0, -2, 0, -21636),
    (2, 1, -1, 0, 24208),
    (2, 1, 0, 0, 30824),
    (1, 0, -1, 0, -8379),
    (1, 1, 0, 0, -16675),
    (2, -1, 1, 0, -12831),
    (2, 0, 2, 0, -10445),
    (4, 0, 0, 0, -11650),
    (2, 0, -3, 0, 14403),
    (0, 1, -2, 0, -7003),
    (2, 0, -1, 2, 0),
    (2, -1, -2, 0, 10056),
    (1, 0, 1, 0, 6322),
    (2, -2, 0, 0, -9884),
    (0, 1, 2, 0, 5751),
    (0, 2, 0, 0, 0),
    (2, -2, -1, 0, -4950),
    (2, 0, 1, -2, 4130),
    (2, 0, 0, 2, 0),
    (4, -1, -1, 0, -3958),
    (0, 0, 2, 2, 0),
    (3, 0, -1, 0, 3258),
    (2, 1, 1, 0, 2616),
    (4, -1, -2, 0, -1897),
    (0, 2, -1, 0, -2117),
    (2, 2, -1, 0, 2354),
    (2, 1, -2, 0, 0),
    (2, -1, 0, -2, 0),
    (4, 0, 1, 0, -1423),
    (0, 0, 4, 0, -1117),
    (4, -1, 0, 0, -1571),
    (1, 0, -2, 0, -1739),
]


# =============================================================================
# ECCENTRICITY CORRECTION
# =============================================================================


def _eccentricity_correction(T: float, M_power: int) -> float:
    """Calculate eccentricity correction factor.

    The amplitude of terms involving M (Sun's mean anomaly) must be
    multiplied by E or E^2 depending on the power of M.

    Args:
        T: Julian centuries from J2000.0.
        M_power: Absolute value of M coefficient (0, 1, or 2).

    Returns:
        Correction factor.
    """
    if M_power == 0:
        return 1.0

    E = 1.0 - 0.002516 * T - 0.0000074 * T**2

    if M_power == 1:
        return E
    elif M_power >= 2:
        return E * E
    return 1.0


# =============================================================================
# MAIN CALCULATION
# =============================================================================


def calc_moon_position(jd_tt: float) -> Tuple[float, float, float]:
    """Calculate geocentric ecliptic position of the Moon.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (longitude, latitude, distance) where:
        - longitude, latitude in degrees (J2000.0 ecliptic)
        - distance in km
    """
    # Julian centuries from J2000.0
    T = jd_to_julian_centuries(jd_tt)

    # Fundamental arguments
    Lp, D, M, Mp, F = _fundamental_arguments(T)

    # Convert to radians
    D_rad = D * DEG_TO_RAD
    M_rad = M * DEG_TO_RAD
    Mp_rad = Mp * DEG_TO_RAD
    F_rad = F * DEG_TO_RAD

    # Sum longitude terms
    sigma_l = 0.0
    for term in _LONGITUDE_TERMS:
        d, m, mp, f, coeff = term
        arg = d * D_rad + m * M_rad + mp * Mp_rad + f * F_rad
        E = _eccentricity_correction(T, abs(m))
        sigma_l += coeff * E * math.sin(arg)

    # Sum latitude terms
    sigma_b = 0.0
    for term in _LATITUDE_TERMS:
        d, m, mp, f, coeff = term
        arg = d * D_rad + m * M_rad + mp * Mp_rad + f * F_rad
        E = _eccentricity_correction(T, abs(m))
        sigma_b += coeff * E * math.sin(arg)

    # Sum distance terms
    sigma_r = 0.0
    for term in _DISTANCE_TERMS:
        d, m, mp, f, coeff = term
        arg = d * D_rad + m * M_rad + mp * Mp_rad + f * F_rad
        E = _eccentricity_correction(T, abs(m))
        sigma_r += coeff * E * math.cos(arg)

    # Additional corrections for Venus and Jupiter
    A1 = (119.75 + 131.849 * T) * DEG_TO_RAD
    A2 = (53.09 + 479264.290 * T) * DEG_TO_RAD
    A3 = (313.45 + 481266.484 * T) * DEG_TO_RAD

    sigma_l += 3958 * math.sin(A1)
    sigma_l += 1962 * math.sin(Lp * DEG_TO_RAD - F_rad)
    sigma_l += 318 * math.sin(A2)

    sigma_b += -2235 * math.sin(Lp * DEG_TO_RAD)
    sigma_b += 382 * math.sin(A3)
    sigma_b += 175 * math.sin(A1 - F_rad)
    sigma_b += 175 * math.sin(A1 + F_rad)
    sigma_b += 127 * math.sin(Lp * DEG_TO_RAD - Mp_rad)
    sigma_b += -115 * math.sin(Lp * DEG_TO_RAD + Mp_rad)

    # Calculate final coordinates
    # sigma_l and sigma_b are in units of 0.000001 degrees
    longitude = Lp + sigma_l / 1000000.0
    latitude = sigma_b / 1000000.0
    # sigma_r is in km, add to mean distance
    distance_km = 385000.56 + sigma_r / 1000.0

    # Normalize longitude
    longitude = normalize_angle(longitude)

    return longitude, latitude, distance_km


def calc_position(
    jd_tt: float,
    body_id: int,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate geocentric ecliptic position of the Moon.

    This is the main entry point for lunar calculations, matching the
    signature of vsop87.calc_position.

    Args:
        jd_tt: Julian Day in Terrestrial Time.
        body_id: Body ID (must be MOSHIER_MOON = 1).

    Returns:
        Tuple of (lon, lat, dist, dlon, dlat, ddist) where:
        - lon, lat in degrees (ecliptic J2000.0)
        - dist in AU
        - dlon, dlat in degrees/day
        - ddist in AU/day

    Raises:
        ValueError: If body_id is not the Moon.
    """
    if body_id != MOSHIER_MOON:
        raise ValueError(f"Body {body_id} is not the Moon (expected {MOSHIER_MOON})")

    lon, lat, dist_km = calc_moon_position(jd_tt)

    # Convert distance from km to AU
    dist_au = dist_km / AU_KM

    # Calculate velocities using numerical differentiation
    h = 0.001  # 1.44 minutes step

    pos_plus = calc_moon_position(jd_tt + h)
    pos_minus = calc_moon_position(jd_tt - h)

    dlon = pos_plus[0] - pos_minus[0]
    dlat = pos_plus[1] - pos_minus[1]
    ddist_km = pos_plus[2] - pos_minus[2]

    # Handle longitude wrap-around
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0

    dlon /= 2 * h
    dlat /= 2 * h
    ddist_au = (ddist_km / AU_KM) / (2 * h)

    return lon, lat, dist_au, dlon, dlat, ddist_au


def is_moon_body(body_id: int) -> bool:
    """Check if a body ID is the Moon.

    Args:
        body_id: Swiss Ephemeris body ID.

    Returns:
        True if body is the Moon.
    """
    return body_id == MOSHIER_MOON
