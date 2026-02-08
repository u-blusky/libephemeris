"""
ELP 2000-82B lunar theory implementation for Moshier ephemeris.

Implements a truncated version of the ELP 2000-82B lunar theory for
calculating the geocentric position of the Moon.

The ELP (Éphéméride Lunaire Parisienne) theory was developed by
M. Chapront-Touzé and J. Chapront at the Bureau des Longitudes, Paris.

This implementation uses truncated series based on the Moshier implementation
used in Swiss Ephemeris, which provides arcsecond-level accuracy (0.5-5 arcsec)
sufficient for most applications over the range -3000 to +3000 CE.

This module has NO dependencies on Skyfield or SPK files.
Only standard library math module is used for numerical operations.

References:
- Chapront-Touzé, M. and Chapront, J. (1983), "The lunar ephemeris ELP 2000",
  A&A 124, 50-62
- Chapront-Touzé, M. and Chapront, J. (1988), "ELP 2000-85: a semi-analytical
  lunar ephemeris adequate for historical times", A&A 190, 342
- Meeus, J. (1991), "Astronomical Algorithms", Chapter 47
- Moshier, S.L., "Astronomical Ephemeris and Reduction Library"
"""

from __future__ import annotations

import math
from typing import Tuple

from .elp82b_data import (
    # Polynomial coefficients for fundamental arguments
    MOON_MEAN_LON_COEFFS,
    MEAN_ELONGATION_COEFFS,
    SUN_MEAN_ANOMALY_COEFFS,
    MOON_MEAN_ANOMALY_COEFFS,
    MOON_ARG_LAT_COEFFS,
    ASCENDING_NODE_COEFFS,
    # Main series terms
    LONGITUDE_MAIN_TERMS,
    LATITUDE_MAIN_TERMS,
    DISTANCE_MAIN_TERMS,
    # Planetary perturbation terms
    LONGITUDE_PLANETARY_TERMS,
    LATITUDE_PLANETARY_TERMS,
    DISTANCE_PLANETARY_TERMS,
    # Constants
    MOON_MEAN_DISTANCE_KM,
    AU_KM,
    EARTH_ECCENTRICITY_COEFFS,
)
from .utils import (
    DEG_TO_RAD,
    J2000,
    JD_PER_CENTURY,
    normalize_angle,
)

# =============================================================================
# CONSTANTS
# =============================================================================

# Swiss Ephemeris Moon body ID
MOSHIER_MOON = 1

# Mean longitudes of planets for planetary perturbations (in degrees at J2000)
# These are used for the planetary perturbation terms
_PLANET_MEAN_LON_COEFFS = {
    # Mercury
    "Me": (252.2509, 149472.6747),
    # Venus
    "Ve": (181.9798, 58517.8156),
    # Earth
    "Ea": (100.4664, 35999.3729),
    # Mars
    "Ma": (355.4330, 19140.2993),
    # Jupiter
    "Ju": (34.3515, 3034.9057),
    # Saturn
    "Sa": (50.0775, 1222.1138),
}


# =============================================================================
# FUNDAMENTAL ARGUMENTS
# =============================================================================


def _evaluate_polynomial(coeffs: Tuple[float, ...], T: float) -> float:
    """Evaluate polynomial using Horner's method.

    Args:
        coeffs: Polynomial coefficients [a0, a1, a2, ...].
        T: Independent variable.

    Returns:
        Polynomial value.
    """
    result = 0.0
    for c in reversed(coeffs):
        result = result * T + c
    return result


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
    Lp = _evaluate_polynomial(MOON_MEAN_LON_COEFFS, T)
    D = _evaluate_polynomial(MEAN_ELONGATION_COEFFS, T)
    M = _evaluate_polynomial(SUN_MEAN_ANOMALY_COEFFS, T)
    Mp = _evaluate_polynomial(MOON_MEAN_ANOMALY_COEFFS, T)
    F = _evaluate_polynomial(MOON_ARG_LAT_COEFFS, T)

    return Lp, D, M, Mp, F


def _mean_longitude_ascending_node(T: float) -> float:
    """Calculate mean longitude of Moon's ascending node.

    Args:
        T: Julian centuries from J2000.0.

    Returns:
        Omega in degrees.
    """
    return _evaluate_polynomial(ASCENDING_NODE_COEFFS, T)


def _planet_mean_longitudes(
    T: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate mean longitudes of planets for perturbation terms.

    Args:
        T: Julian centuries from J2000.0.

    Returns:
        Tuple of (Me, Ve, Ea, Ma, Ju, Sa) in degrees.
    """
    Me = _PLANET_MEAN_LON_COEFFS["Me"][0] + _PLANET_MEAN_LON_COEFFS["Me"][1] * T
    Ve = _PLANET_MEAN_LON_COEFFS["Ve"][0] + _PLANET_MEAN_LON_COEFFS["Ve"][1] * T
    Ea = _PLANET_MEAN_LON_COEFFS["Ea"][0] + _PLANET_MEAN_LON_COEFFS["Ea"][1] * T
    Ma = _PLANET_MEAN_LON_COEFFS["Ma"][0] + _PLANET_MEAN_LON_COEFFS["Ma"][1] * T
    Ju = _PLANET_MEAN_LON_COEFFS["Ju"][0] + _PLANET_MEAN_LON_COEFFS["Ju"][1] * T
    Sa = _PLANET_MEAN_LON_COEFFS["Sa"][0] + _PLANET_MEAN_LON_COEFFS["Sa"][1] * T

    return Me, Ve, Ea, Ma, Ju, Sa


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

    E = (
        EARTH_ECCENTRICITY_COEFFS[0]
        + EARTH_ECCENTRICITY_COEFFS[1] * T
        + EARTH_ECCENTRICITY_COEFFS[2] * T * T
    )

    if M_power == 1:
        return E
    elif M_power >= 2:
        return E * E
    return 1.0


# =============================================================================
# SERIES SUMMATION
# =============================================================================


def _sum_longitude_series(
    T: float,
    D_rad: float,
    M_rad: float,
    Mp_rad: float,
    F_rad: float,
) -> float:
    """Sum longitude series terms.

    Args:
        T: Julian centuries from J2000.0.
        D_rad: Mean elongation in radians.
        M_rad: Sun's mean anomaly in radians.
        Mp_rad: Moon's mean anomaly in radians.
        F_rad: Moon's argument of latitude in radians.

    Returns:
        Sum of longitude terms in units of 0.000001 degrees.
    """
    sigma_l = 0.0

    for term in LONGITUDE_MAIN_TERMS:
        d, m, mp, f, coeff = term
        if coeff == 0:
            continue
        arg = d * D_rad + m * M_rad + mp * Mp_rad + f * F_rad
        E = _eccentricity_correction(T, abs(m))
        sigma_l += coeff * E * math.sin(arg)

    return sigma_l


def _sum_latitude_series(
    T: float,
    D_rad: float,
    M_rad: float,
    Mp_rad: float,
    F_rad: float,
) -> float:
    """Sum latitude series terms.

    Args:
        T: Julian centuries from J2000.0.
        D_rad: Mean elongation in radians.
        M_rad: Sun's mean anomaly in radians.
        Mp_rad: Moon's mean anomaly in radians.
        F_rad: Moon's argument of latitude in radians.

    Returns:
        Sum of latitude terms in units of 0.000001 degrees.
    """
    sigma_b = 0.0

    for term in LATITUDE_MAIN_TERMS:
        d, m, mp, f, coeff = term
        if coeff == 0:
            continue
        arg = d * D_rad + m * M_rad + mp * Mp_rad + f * F_rad
        E = _eccentricity_correction(T, abs(m))
        sigma_b += coeff * E * math.sin(arg)

    return sigma_b


def _sum_distance_series(
    T: float,
    D_rad: float,
    M_rad: float,
    Mp_rad: float,
    F_rad: float,
) -> float:
    """Sum distance series terms.

    Args:
        T: Julian centuries from J2000.0.
        D_rad: Mean elongation in radians.
        M_rad: Sun's mean anomaly in radians.
        Mp_rad: Moon's mean anomaly in radians.
        F_rad: Moon's argument of latitude in radians.

    Returns:
        Sum of distance terms in km.
    """
    sigma_r = 0.0

    for term in DISTANCE_MAIN_TERMS:
        d, m, mp, f, coeff = term
        if coeff == 0:
            continue
        arg = d * D_rad + m * M_rad + mp * Mp_rad + f * F_rad
        E = _eccentricity_correction(T, abs(m))
        sigma_r += coeff * E * math.cos(arg)

    return sigma_r


def _sum_planetary_longitude(
    T: float,
    D_rad: float,
    M_rad: float,
    Mp_rad: float,
    F_rad: float,
    planet_lons: Tuple[float, float, float, float, float, float],
) -> float:
    """Sum planetary perturbation terms for longitude.

    Args:
        T: Julian centuries from J2000.0.
        D_rad: Mean elongation in radians.
        M_rad: Sun's mean anomaly in radians.
        Mp_rad: Moon's mean anomaly in radians.
        F_rad: Moon's argument of latitude in radians.
        planet_lons: Tuple of (Me, Ve, Ea, Ma, Ju, Sa) in radians.

    Returns:
        Sum of planetary longitude perturbations in units of 0.000001 degrees.
    """
    Me_rad, Ve_rad, Ea_rad, Ma_rad, Ju_rad, Sa_rad = planet_lons

    sigma_pl = 0.0

    for term in LONGITUDE_PLANETARY_TERMS:
        d, m, mp, f, me, ve, ea, ma, ju, sa, coeff = term
        if coeff == 0:
            continue
        arg = (
            d * D_rad
            + m * M_rad
            + mp * Mp_rad
            + f * F_rad
            + me * Me_rad
            + ve * Ve_rad
            + ea * Ea_rad
            + ma * Ma_rad
            + ju * Ju_rad
            + sa * Sa_rad
        )
        sigma_pl += coeff * math.sin(arg)

    return sigma_pl


def _sum_planetary_latitude(
    T: float,
    D_rad: float,
    M_rad: float,
    Mp_rad: float,
    F_rad: float,
    planet_lons: Tuple[float, float, float, float, float, float],
) -> float:
    """Sum planetary perturbation terms for latitude.

    Args:
        T: Julian centuries from J2000.0.
        D_rad: Mean elongation in radians.
        M_rad: Sun's mean anomaly in radians.
        Mp_rad: Moon's mean anomaly in radians.
        F_rad: Moon's argument of latitude in radians.
        planet_lons: Tuple of (Me, Ve, Ea, Ma, Ju, Sa) in radians.

    Returns:
        Sum of planetary latitude perturbations in units of 0.000001 degrees.
    """
    Me_rad, Ve_rad, Ea_rad, Ma_rad, Ju_rad, Sa_rad = planet_lons

    sigma_pb = 0.0

    for term in LATITUDE_PLANETARY_TERMS:
        d, m, mp, f, me, ve, ea, ma, ju, sa, coeff = term
        if coeff == 0:
            continue
        arg = (
            d * D_rad
            + m * M_rad
            + mp * Mp_rad
            + f * F_rad
            + me * Me_rad
            + ve * Ve_rad
            + ea * Ea_rad
            + ma * Ma_rad
            + ju * Ju_rad
            + sa * Sa_rad
        )
        sigma_pb += coeff * math.sin(arg)

    return sigma_pb


def _sum_planetary_distance(
    T: float,
    D_rad: float,
    M_rad: float,
    Mp_rad: float,
    F_rad: float,
    planet_lons: Tuple[float, float, float, float, float, float],
) -> float:
    """Sum planetary perturbation terms for distance.

    Args:
        T: Julian centuries from J2000.0.
        D_rad: Mean elongation in radians.
        M_rad: Sun's mean anomaly in radians.
        Mp_rad: Moon's mean anomaly in radians.
        F_rad: Moon's argument of latitude in radians.
        planet_lons: Tuple of (Me, Ve, Ea, Ma, Ju, Sa) in radians.

    Returns:
        Sum of planetary distance perturbations in km.
    """
    Me_rad, Ve_rad, Ea_rad, Ma_rad, Ju_rad, Sa_rad = planet_lons

    sigma_pr = 0.0

    for term in DISTANCE_PLANETARY_TERMS:
        d, m, mp, f, me, ve, ea, ma, ju, sa, coeff = term
        if coeff == 0:
            continue
        arg = (
            d * D_rad
            + m * M_rad
            + mp * Mp_rad
            + f * F_rad
            + me * Me_rad
            + ve * Ve_rad
            + ea * Ea_rad
            + ma * Ma_rad
            + ju * Ju_rad
            + sa * Sa_rad
        )
        sigma_pr += coeff * math.cos(arg)

    return sigma_pr


# =============================================================================
# MAIN CALCULATION
# =============================================================================


def calc_moon_position(jd_tt: float) -> Tuple[float, float, float]:
    """Calculate geocentric ecliptic position of the Moon.

    Uses the ELP 2000-82B theory with truncated series for high precision.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (longitude, latitude, distance) where:
        - longitude, latitude in degrees (J2000.0 ecliptic)
        - distance in km
    """
    # Julian centuries from J2000.0
    T = (jd_tt - J2000) / JD_PER_CENTURY

    # Fundamental arguments in degrees
    Lp, D, M, Mp, F = _fundamental_arguments(T)

    # Convert to radians
    D_rad = D * DEG_TO_RAD
    M_rad = M * DEG_TO_RAD
    Mp_rad = Mp * DEG_TO_RAD
    F_rad = F * DEG_TO_RAD
    Lp_rad = Lp * DEG_TO_RAD

    # Get planet mean longitudes for perturbations
    Me, Ve, Ea, Ma, Ju, Sa = _planet_mean_longitudes(T)
    planet_lons_rad = (
        Me * DEG_TO_RAD,
        Ve * DEG_TO_RAD,
        Ea * DEG_TO_RAD,
        Ma * DEG_TO_RAD,
        Ju * DEG_TO_RAD,
        Sa * DEG_TO_RAD,
    )

    # Sum main series
    sigma_l = _sum_longitude_series(T, D_rad, M_rad, Mp_rad, F_rad)
    sigma_b = _sum_latitude_series(T, D_rad, M_rad, Mp_rad, F_rad)
    sigma_r = _sum_distance_series(T, D_rad, M_rad, Mp_rad, F_rad)

    # Add planetary perturbations
    sigma_l += _sum_planetary_longitude(T, D_rad, M_rad, Mp_rad, F_rad, planet_lons_rad)
    sigma_b += _sum_planetary_latitude(T, D_rad, M_rad, Mp_rad, F_rad, planet_lons_rad)
    sigma_r += _sum_planetary_distance(T, D_rad, M_rad, Mp_rad, F_rad, planet_lons_rad)

    # Additional corrections for Venus and Jupiter (from Meeus)
    A1 = (119.75 + 131.849 * T) * DEG_TO_RAD
    A2 = (53.09 + 479264.290 * T) * DEG_TO_RAD
    A3 = (313.45 + 481266.484 * T) * DEG_TO_RAD

    # Longitude corrections
    sigma_l += 3958 * math.sin(A1)
    sigma_l += 1962 * math.sin(Lp_rad - F_rad)
    sigma_l += 318 * math.sin(A2)

    # Latitude corrections
    sigma_b += -2235 * math.sin(Lp_rad)
    sigma_b += 382 * math.sin(A3)
    sigma_b += 175 * math.sin(A1 - F_rad)
    sigma_b += 175 * math.sin(A1 + F_rad)
    sigma_b += 127 * math.sin(Lp_rad - Mp_rad)
    sigma_b += -115 * math.sin(Lp_rad + Mp_rad)

    # Calculate final coordinates
    # sigma_l and sigma_b are in units of 0.000001 degrees
    longitude = Lp + sigma_l / 1000000.0
    latitude = sigma_b / 1000000.0
    # sigma_r is in meters originally but our data is in km
    # Distance terms are in units of 0.001 km
    distance_km = MOON_MEAN_DISTANCE_KM + sigma_r / 1000.0

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
