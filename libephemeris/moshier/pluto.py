"""
Pluto ephemeris for Moshier calculations.

Implements a polynomial ephemeris for Pluto based on the theory by
Chapront and Francou, as adapted by Meeus in "Astronomical Algorithms".

This provides heliocentric rectangular coordinates of Pluto valid from
approximately 1885 to 2099 with accuracy of about 0.07 degrees in longitude.

For dates outside this range, a simplified Keplerian solution is used.

This module has NO dependencies on Skyfield or SPK files.
Only numpy is used for numerical operations.

References:
- Meeus, J. (1991), "Astronomical Algorithms", Chapter 37
- Chapront, J. and Francou, G. (1995), Pluto theory tables
"""

from __future__ import annotations

import math
from typing import Tuple

from .utils import (
    DEG_TO_RAD,
    J2000,
    JD_PER_CENTURY,
    RAD_TO_DEG,
    cartesian_to_spherical,
    jd_to_julian_centuries,
    normalize_angle,
)
from .vsop87 import calc_earth_heliocentric

# =============================================================================
# CONSTANTS
# =============================================================================

# Swiss Ephemeris Pluto body ID
MOSHIER_PLUTO = 9

# Pluto orbital elements at J2000.0
PLUTO_EPOCH_JD = 2451545.0  # J2000.0
PLUTO_SEMI_MAJOR_AXIS = 39.48211675  # AU
PLUTO_ECCENTRICITY = 0.2488273
PLUTO_INCLINATION = 17.14001206  # degrees
PLUTO_LONG_ASC_NODE = 110.30393684  # degrees
PLUTO_ARG_PERIHELION = 224.06891629  # degrees
PLUTO_MEAN_LONGITUDE = 238.92881030  # degrees

# Daily rates
PLUTO_MEAN_MOTION = 0.00397570  # degrees/day


# =============================================================================
# PLUTO PERIODIC TERMS (Meeus, Chapter 37)
# =============================================================================

# Each term: (J, S, P, lon_sin, lon_cos, lat_sin, lat_cos, rad_sin, rad_cos)
# J, S, P are integer arguments
# lon/lat in 0.000001 degrees, rad in 0.0000001 AU

_PLUTO_TERMS = [
    # J   S   P     lon_sin   lon_cos    lat_sin   lat_cos    rad_sin   rad_cos
    (0, 0, 1, -19799805, 19850055, -5452852, -14974862, 66865439, 68951812),
    (0, 0, 2, 897144, -4954829, 3527812, 1672790, -11827535, -332538),
    (0, 0, 3, 611149, 1211027, -1050748, 327647, 1593179, -1438890),
    (0, 0, 4, -341243, -189585, 178690, -292153, -18444, 483220),
    (0, 0, 5, 129287, -34992, 18763, 100340, -65977, -85431),
    (0, 0, 6, -38164, 30893, -30697, -25823, 31174, -6032),
    (0, 1, -1, 20442, -9987, 4878, 11248, -5794, 22161),
    (0, 1, 0, -4063, -5071, 226, -64, 4601, 4032),
    (0, 1, 1, -6016, -3336, 2030, -836, -1729, 234),
    (0, 1, 2, -3956, 3039, 69, -604, -415, 702),
    (0, 1, 3, -667, 3572, -247, -567, 239, 723),
    (0, 2, -2, 1276, 501, -57, 1, 67, -67),
    (0, 2, -1, 1152, -917, -122, 175, 1034, -451),
    (0, 2, 0, 630, -1277, -49, -164, -129, 504),
    (1, -1, 0, 2571, -459, -197, 199, 480, -231),
    (1, -1, 1, 899, -1449, -25, 217, 2, -441),
    (1, 0, -3, -1016, 1043, 589, -248, -3359, 265),
    (1, 0, -2, -2343, -1012, -269, 711, 7856, -7832),
    (1, 0, -1, 7042, 788, 185, 193, 36, 45763),
    (1, 0, 0, 1199, -338, 315, 807, 8663, 8547),
    (1, 0, 1, 418, -67, -130, -43, -809, -769),
    (1, 0, 2, 120, -274, 5, 3, 263, -144),
    (1, 0, 3, -60, -159, 2, 17, -126, 32),
    (1, 0, 4, -82, -29, 2, 5, -35, -16),
    (1, 1, -3, -36, -29, 2, 3, -19, -4),
    (1, 1, -2, -40, 7, 3, 1, -15, 8),
    (1, 1, -1, -14, 22, 2, -1, -4, 12),
    (1, 1, 0, 4, 13, 1, -1, 5, 6),
    (1, 1, 1, 5, 2, 0, -1, 3, 1),
    (1, 1, 3, -1, 0, 0, 0, 1, -1),
    (2, 0, -6, 2, 0, 0, -2, 2, 2),
    (2, 0, -5, -4, 5, 2, 2, -2, -2),
    (2, 0, -4, 4, -7, -7, 0, 14, 13),
    (2, 0, -3, 14, 24, 10, -8, -63, 66),
    (2, 0, -2, -49, -34, -3, 20, 136, -236),
    (2, 0, -1, 163, -48, 6, 5, 273, 1065),
    (2, 0, 0, 9, -24, 14, 17, 251, 149),
    (2, 0, 1, -4, 1, -2, 0, -25, -9),
    (2, 0, 2, -3, 1, 0, 0, 9, -2),
    (2, 0, 3, 1, 3, 0, 0, -8, 7),
    (3, 0, -2, -3, -1, 0, 1, 2, -10),
    (3, 0, -1, 5, -3, 0, 0, 19, 35),
    (3, 0, 0, 0, 0, 1, 0, 10, 3),
]


# =============================================================================
# PLUTO CALCULATION
# =============================================================================


def _calc_pluto_periodic(jd_tt: float) -> Tuple[float, float, float]:
    """Calculate Pluto heliocentric coordinates using periodic terms.

    Valid approximately 1885-2099.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (longitude, latitude, radius) in degrees and AU.
    """
    # Julian centuries from J2000.0
    T = jd_to_julian_centuries(jd_tt)

    # Calculate the arguments J, S, P
    # These are mean longitudes of Jupiter, Saturn, and Pluto
    J = (34.35 + 3034.9057 * T) * DEG_TO_RAD
    S = (50.08 + 1222.1138 * T) * DEG_TO_RAD
    P = (238.96 + 144.9600 * T) * DEG_TO_RAD

    # Sum the periodic terms
    sum_lon = 0.0
    sum_lat = 0.0
    sum_rad = 0.0

    for term in _PLUTO_TERMS:
        j_mult, s_mult, p_mult = term[0], term[1], term[2]
        lon_sin, lon_cos = term[3], term[4]
        lat_sin, lat_cos = term[5], term[6]
        rad_sin, rad_cos = term[7], term[8]

        alpha = j_mult * J + s_mult * S + p_mult * P
        sin_alpha = math.sin(alpha)
        cos_alpha = math.cos(alpha)

        sum_lon += lon_sin * sin_alpha + lon_cos * cos_alpha
        sum_lat += lat_sin * sin_alpha + lat_cos * cos_alpha
        sum_rad += rad_sin * sin_alpha + rad_cos * cos_alpha

    # Base values
    lon_base = 238.958116 + 144.96 * T
    lat_base = -3.908239
    rad_base = 40.7241346

    # Final values
    longitude = lon_base + sum_lon * 1e-6
    latitude = lat_base + sum_lat * 1e-6
    radius = rad_base + sum_rad * 1e-7

    return normalize_angle(longitude), latitude, radius


def _calc_pluto_keplerian(jd_tt: float) -> Tuple[float, float, float]:
    """Calculate Pluto position using Keplerian elements.

    This is a fallback for dates outside the periodic term range.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (longitude, latitude, radius) in degrees and AU.
    """
    # Days from J2000.0
    d = jd_tt - J2000

    # Mean anomaly
    n = PLUTO_MEAN_MOTION  # degrees/day
    M = PLUTO_MEAN_LONGITUDE - PLUTO_ARG_PERIHELION - PLUTO_LONG_ASC_NODE + n * d
    M = normalize_angle(M) * DEG_TO_RAD

    # Solve Kepler's equation (Newton-Raphson)
    e = PLUTO_ECCENTRICITY
    E = M
    for _ in range(10):
        delta = E - e * math.sin(E) - M
        E -= delta / (1 - e * math.cos(E))
        if abs(delta) < 1e-12:
            break

    # True anomaly
    cos_E = math.cos(E)
    sin_E = math.sin(E)
    sqrt_1_e2 = math.sqrt(1 - e * e)
    true_anom = math.atan2(sqrt_1_e2 * sin_E, cos_E - e)

    # Radius
    a = PLUTO_SEMI_MAJOR_AXIS
    r = a * (1 - e * cos_E)

    # Argument of latitude
    omega = PLUTO_ARG_PERIHELION * DEG_TO_RAD
    u = omega + true_anom

    # Orbital plane to ecliptic
    Omega = PLUTO_LONG_ASC_NODE * DEG_TO_RAD
    i = PLUTO_INCLINATION * DEG_TO_RAD

    cos_u = math.cos(u)
    sin_u = math.sin(u)
    cos_i = math.cos(i)
    sin_i = math.sin(i)
    cos_Omega = math.cos(Omega)
    sin_Omega = math.sin(Omega)

    # Heliocentric rectangular coordinates
    x = r * (cos_Omega * cos_u - sin_Omega * sin_u * cos_i)
    y = r * (sin_Omega * cos_u + cos_Omega * sin_u * cos_i)
    z = r * sin_u * sin_i

    # Convert to spherical
    longitude, latitude, radius = cartesian_to_spherical(x, y, z)

    return longitude, latitude, radius


def calc_pluto_heliocentric(jd_tt: float) -> Tuple[float, float, float]:
    """Calculate heliocentric ecliptic coordinates of Pluto.

    Uses periodic terms for 1885-2099, Keplerian elements otherwise.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (longitude, latitude, radius) in degrees and AU.
    """
    # Check if within periodic term range
    year = 2000.0 + (jd_tt - J2000) / 365.25

    if 1885.0 <= year <= 2099.0:
        return _calc_pluto_periodic(jd_tt)
    else:
        return _calc_pluto_keplerian(jd_tt)


def calc_pluto_geocentric(jd_tt: float) -> Tuple[float, float, float]:
    """Calculate geocentric ecliptic coordinates of Pluto.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (longitude, latitude, distance) in degrees and AU.
    """
    # Heliocentric coordinates of Pluto
    pL, pB, pR = calc_pluto_heliocentric(jd_tt)

    # Heliocentric coordinates of Earth
    eL, eB, eR = calc_earth_heliocentric(jd_tt)

    # Convert to rectangular
    pL_rad = pL * DEG_TO_RAD
    pB_rad = pB * DEG_TO_RAD
    eL_rad = eL * DEG_TO_RAD
    eB_rad = eB * DEG_TO_RAD

    # Pluto heliocentric rectangular
    px = pR * math.cos(pB_rad) * math.cos(pL_rad)
    py = pR * math.cos(pB_rad) * math.sin(pL_rad)
    pz = pR * math.sin(pB_rad)

    # Earth heliocentric rectangular
    ex = eR * math.cos(eB_rad) * math.cos(eL_rad)
    ey = eR * math.cos(eB_rad) * math.sin(eL_rad)
    ez = eR * math.sin(eB_rad)

    # Geocentric rectangular
    gx = px - ex
    gy = py - ey
    gz = pz - ez

    # Convert back to spherical
    lon, lat, dist = cartesian_to_spherical(gx, gy, gz)

    return lon, lat, dist


# =============================================================================
# PUBLIC API
# =============================================================================


def calc_position(
    jd_tt: float,
    body_id: int,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate geocentric ecliptic position of Pluto.

    This is the main entry point for Pluto calculations, matching the
    signature of vsop87.calc_position.

    Args:
        jd_tt: Julian Day in Terrestrial Time.
        body_id: Body ID (must be MOSHIER_PLUTO = 9).

    Returns:
        Tuple of (lon, lat, dist, dlon, dlat, ddist) where:
        - lon, lat in degrees (ecliptic J2000.0)
        - dist in AU
        - dlon, dlat in degrees/day
        - ddist in AU/day

    Raises:
        ValueError: If body_id is not Pluto.
    """
    if body_id != MOSHIER_PLUTO:
        raise ValueError(f"Body {body_id} is not Pluto (expected {MOSHIER_PLUTO})")

    lon, lat, dist = calc_pluto_geocentric(jd_tt)

    # Calculate velocities using numerical differentiation
    h = 0.01  # Larger step for slow-moving Pluto

    pos_plus = calc_pluto_geocentric(jd_tt + h)
    pos_minus = calc_pluto_geocentric(jd_tt - h)

    dlon = pos_plus[0] - pos_minus[0]
    dlat = pos_plus[1] - pos_minus[1]
    ddist = pos_plus[2] - pos_minus[2]

    # Handle longitude wrap-around
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0

    dlon /= 2 * h
    dlat /= 2 * h
    ddist /= 2 * h

    return lon, lat, dist, dlon, dlat, ddist


def is_pluto_body(body_id: int) -> bool:
    """Check if a body ID is Pluto.

    Args:
        body_id: Swiss Ephemeris body ID.

    Returns:
        True if body is Pluto.
    """
    return body_id == MOSHIER_PLUTO
