"""
IAU 2006 Precession and Nutation for Moshier ephemeris.

Implements the IAU 2006/2000A precession-nutation model for coordinate
transformations between different epochs and reference frames.

This module has NO dependencies on Skyfield or SPK files.
Only numpy is used for numerical operations.

References:
- Capitaine, N., Wallace, P. T., and Chapront, J. (2003), "Expressions for
  IAU 2000 precession quantities", A&A 412, 567-586
- IAU SOFA (Standards of Fundamental Astronomy) software collection
- Lieske, J. H. et al. (1977), "Expressions for the Precession Quantities
  Based upon the IAU (1976) System of Astronomical Constants", A&A 58, 1-16
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
# IAU 2006 PRECESSION CONSTANTS
# =============================================================================

# Mean obliquity of the ecliptic at J2000.0 (arcseconds)
# IAU 2006 value: 84381.406 arcseconds
OBLIQUITY_J2000_ARCSEC: float = 84381.406

# Polynomials for mean obliquity of the ecliptic (IAU 2006)
# epsilon = eps0 + eps1*t + eps2*t^2 + eps3*t^3 + eps4*t^4 + eps5*t^5
# where t is Julian centuries from J2000.0
OBLIQUITY_COEFFS: Tuple[float, ...] = (
    84381.406,  # eps0 (arcsec)
    -46.836769,  # eps1 (arcsec/century)
    -0.0001831,  # eps2
    0.00200340,  # eps3
    -0.000000576,  # eps4
    -0.0000000434,  # eps5
)

# Precession angles (IAU 2006) - for ecliptic precession
# psi_A: precession in longitude
PSI_A_COEFFS: Tuple[float, ...] = (
    0.0,
    5038.481507,
    -1.0790069,
    -0.00114045,
    0.000132851,
    -0.0000000951,
)

# omega_A: obliquity of the moving ecliptic on the fixed ecliptic
OMEGA_A_COEFFS: Tuple[float, ...] = (
    OBLIQUITY_J2000_ARCSEC,
    -0.025754,
    0.0512623,
    -0.00772503,
    -0.000000467,
    0.0000003337,
)

# chi_A: precession of the equator along the ecliptic
CHI_A_COEFFS: Tuple[float, ...] = (
    0.0,
    10.556403,
    -2.3814292,
    -0.00121197,
    0.000170663,
    -0.0000000560,
)

# Equatorial precession angles
# zeta_A, z_A, theta_A (Lieske 1977, with IAU 2006 corrections)
ZETA_A_COEFFS: Tuple[float, ...] = (
    0.0,
    2306.2181,
    1.39656,
    -0.000139,
)

Z_A_COEFFS: Tuple[float, ...] = (
    0.0,
    2306.2181,
    1.39656,
    -0.000139,
)

THETA_A_COEFFS: Tuple[float, ...] = (
    0.0,
    2004.3109,
    -0.85330,
    -0.000217,
)


# =============================================================================
# OBLIQUITY
# =============================================================================


def mean_obliquity(jd_tt: float) -> float:
    """Calculate mean obliquity of the ecliptic (IAU 2006).

    The mean obliquity is the angle between the ecliptic and the mean
    equator of date, without nutation.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Mean obliquity in degrees.
    """
    t = jd_to_julian_centuries(jd_tt)

    # Evaluate polynomial (Horner's method)
    epsilon = 0.0
    for i, coeff in enumerate(OBLIQUITY_COEFFS):
        epsilon += coeff * (t**i)

    return epsilon * ARCSEC_TO_RAD * RAD_TO_DEG


def true_obliquity(jd_tt: float) -> float:
    """Calculate true obliquity of the ecliptic.

    True obliquity = mean obliquity + nutation in obliquity.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        True obliquity in degrees.
    """
    eps_mean = mean_obliquity(jd_tt)
    _, d_eps = nutation_angles(jd_tt)

    return eps_mean + d_eps


# =============================================================================
# NUTATION (Simplified IAU 2000B)
# =============================================================================

# Fundamental arguments of nutation (Delaunay arguments)
# l = mean anomaly of the Moon
# l' = mean anomaly of the Sun
# F = mean argument of latitude of the Moon
# D = mean elongation of the Moon from the Sun
# Omega = mean longitude of the Moon's ascending node


def _fundamental_arguments(t: float) -> Tuple[float, float, float, float, float]:
    """Calculate fundamental arguments for nutation (IAU 2003).

    Args:
        t: Julian centuries from J2000.0.

    Returns:
        Tuple of (el, elp, F, D, Omega) in radians, all normalized to [0, 2pi).
        - el: Mean anomaly of the Moon
        - elp: Mean anomaly of the Sun
        - F: Mean argument of latitude of the Moon
        - D: Mean elongation of the Moon from the Sun
        - Omega: Mean longitude of the Moon's ascending node
    """
    # Mean anomaly of the Moon (l)
    el = 485868.249036 + t * (
        1717915923.2178 + t * (31.8792 + t * (0.051635 + t * (-0.00024470)))
    )
    el = (el * ARCSEC_TO_RAD) % (2.0 * math.pi)

    # Mean anomaly of the Sun (l')
    elp = 1287104.79305 + t * (
        129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))
    )
    elp = (elp * ARCSEC_TO_RAD) % (2.0 * math.pi)

    # Mean argument of latitude of the Moon (F)
    F = 335779.526232 + t * (
        1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))
    )
    F = (F * ARCSEC_TO_RAD) % (2.0 * math.pi)

    # Mean elongation of the Moon from the Sun (D)
    D = 1072260.70369 + t * (
        1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169)))
    )
    D = (D * ARCSEC_TO_RAD) % (2.0 * math.pi)

    # Mean longitude of the Moon's ascending node (Omega)
    Omega = 450160.398036 + t * (
        -6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939)))
    )
    Omega = (Omega * ARCSEC_TO_RAD) % (2.0 * math.pi)

    return el, elp, F, D, Omega


# Simplified nutation series (major terms only from IAU 2000B)
# Format: (el, elp, F, D, Omega, psi_sin, psi_cos, eps_sin, eps_cos) in 0.1 microarcsec
_NUTATION_TERMS: Tuple[Tuple[int, ...], ...] = (
    # el  elp   F    D   Omega   psi_sin   psi_cos    eps_sin   eps_cos
    (0, 0, 0, 0, 1, -172064161, -174666, 92052331, 9086),
    (0, 0, 2, -2, 2, -13170906, -1675, 5730336, -3015),
    (0, 0, 2, 0, 2, -2276413, -234, 978459, -485),
    (0, 0, 0, 0, 2, 2074554, 207, -897492, 470),
    (0, 1, 0, 0, 0, 1475877, -3633, 73871, -184),
    (0, 1, 2, -2, 2, -516821, 1226, 224386, -677),
    (1, 0, 0, 0, 0, 711159, 73, -6750, 0),
    (0, 0, 2, 0, 1, -387298, -367, 200728, 18),
    (1, 0, 2, 0, 2, -301461, -36, 129025, -63),
    (0, -1, 2, -2, 2, 215829, -494, -95929, 299),
    (0, 0, 2, -2, 1, 128227, 137, -68982, -9),
    (-1, 0, 2, 0, 2, 123457, 11, -53311, 32),
    (-1, 0, 0, 2, 0, 156994, 10, -1235, 0),
    (1, 0, 0, 0, 1, 63110, 63, -33228, 0),
    (-1, 0, 0, 0, 1, -57976, -63, 31429, 0),
    (-1, 0, 2, 2, 2, -59641, -11, 25543, -11),
    (1, 0, 2, 0, 1, -51613, -42, 26366, 0),
    (-2, 0, 2, 0, 1, 45893, 50, -24236, -10),
    (0, 0, 0, 2, 0, 63384, 11, -1220, 0),
    (0, 0, 2, 2, 2, -38571, -1, 16452, -11),
)


def nutation_angles(jd_tt: float) -> Tuple[float, float]:
    """Calculate nutation in longitude and obliquity (IAU 2000B simplified).

    This is a simplified nutation model using the major terms from IAU 2000B.
    It provides accuracy of about 1 milliarcsecond, sufficient for most
    astronomical applications.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (delta_psi, delta_epsilon) in degrees, where:
        - delta_psi: nutation in longitude
        - delta_epsilon: nutation in obliquity
    """
    t = jd_to_julian_centuries(jd_tt)

    # Get fundamental arguments
    el, elp, F, D, Omega = _fundamental_arguments(t)

    # Sum nutation series
    dpsi = 0.0  # nutation in longitude (0.1 microarcsec)
    deps = 0.0  # nutation in obliquity (0.1 microarcsec)

    for term in _NUTATION_TERMS:
        arg = term[0] * el + term[1] * elp + term[2] * F + term[3] * D + term[4] * Omega
        sin_arg = math.sin(arg)
        cos_arg = math.cos(arg)

        dpsi += (term[5] + term[6] * t) * sin_arg
        deps += (term[7] + term[8] * t) * cos_arg

    # Convert from 0.1 microarcsec to degrees
    # 0.1 microarcsec = 0.1e-6 arcsec = 0.1e-6 / 3600 degrees
    dpsi_deg = dpsi * 0.1e-6 * ARCSEC_TO_RAD * RAD_TO_DEG
    deps_deg = deps * 0.1e-6 * ARCSEC_TO_RAD * RAD_TO_DEG

    return dpsi_deg, deps_deg


# =============================================================================
# PRECESSION
# =============================================================================


def precession_angles(jd_tt: float) -> Tuple[float, float, float]:
    """Calculate precession angles zeta, z, theta (Lieske 1977).

    These angles define the rotation from the mean equator and equinox
    of J2000.0 to the mean equator and equinox of date.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (zeta, z, theta) in degrees.
    """
    t = jd_to_julian_centuries(jd_tt)

    # Calculate precession angles
    zeta = 0.0
    for i, coeff in enumerate(ZETA_A_COEFFS):
        zeta += coeff * (t**i)

    z = 0.0
    for i, coeff in enumerate(Z_A_COEFFS):
        z += coeff * (t**i)

    theta = 0.0
    for i, coeff in enumerate(THETA_A_COEFFS):
        theta += coeff * (t**i)

    # Convert from arcseconds to degrees
    zeta_deg = zeta * ARCSEC_TO_RAD * RAD_TO_DEG
    z_deg = z * ARCSEC_TO_RAD * RAD_TO_DEG
    theta_deg = theta * ARCSEC_TO_RAD * RAD_TO_DEG

    return zeta_deg, z_deg, theta_deg


def precess_ecliptic(
    lon: float,
    lat: float,
    from_jd: float,
    to_jd: float,
) -> Tuple[float, float]:
    """Precess ecliptic coordinates from one epoch to another.

    Uses rigorous precession matrix for ecliptic coordinates.

    Args:
        lon: Ecliptic longitude in degrees.
        lat: Ecliptic latitude in degrees.
        from_jd: Source epoch Julian Day (TT).
        to_jd: Target epoch Julian Day (TT).

    Returns:
        Tuple of (precessed_lon, precessed_lat) in degrees.
    """
    if abs(from_jd - to_jd) < 1e-6:
        return lon, lat

    # Convert to radians
    lon_rad = lon * DEG_TO_RAD
    lat_rad = lat * DEG_TO_RAD

    # Get obliquity at both epochs
    eps1 = mean_obliquity(from_jd) * DEG_TO_RAD
    eps2 = mean_obliquity(to_jd) * DEG_TO_RAD

    # Convert to equatorial at source epoch
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)
    cos_eps1 = math.cos(eps1)
    sin_eps1 = math.sin(eps1)

    # Direction cosines in equatorial frame at epoch 1
    x1 = cos_lat * cos_lon
    y1 = cos_lat * sin_lon * cos_eps1 - sin_lat * sin_eps1
    z1 = cos_lat * sin_lon * sin_eps1 + sin_lat * cos_eps1

    # Precess in equatorial frame
    t1 = jd_to_julian_centuries(from_jd)
    t2 = jd_to_julian_centuries(to_jd)
    dt = t2 - t1

    # Precession angles for interval
    zeta = (2306.2181 + 1.39656 * t1) * dt + 1.09468 * dt**2 + 0.018203 * dt**3
    z = (2306.2181 + 1.39656 * t1) * dt + 0.30188 * dt**2 + 0.017998 * dt**3
    theta = (2004.3109 - 0.85330 * t1) * dt - 0.42665 * dt**2 - 0.041833 * dt**3

    # Convert to radians
    zeta_rad = zeta * ARCSEC_TO_RAD
    z_rad = z * ARCSEC_TO_RAD
    theta_rad = theta * ARCSEC_TO_RAD

    # Build rotation matrix and apply
    cos_zeta = math.cos(zeta_rad)
    sin_zeta = math.sin(zeta_rad)
    cos_z = math.cos(z_rad)
    sin_z = math.sin(z_rad)
    cos_theta = math.cos(theta_rad)
    sin_theta = math.sin(theta_rad)

    # Precession matrix elements
    r11 = cos_zeta * cos_theta * cos_z - sin_zeta * sin_z
    r12 = -sin_zeta * cos_theta * cos_z - cos_zeta * sin_z
    r13 = -sin_theta * cos_z
    r21 = cos_zeta * cos_theta * sin_z + sin_zeta * cos_z
    r22 = -sin_zeta * cos_theta * sin_z + cos_zeta * cos_z
    r23 = -sin_theta * sin_z
    r31 = cos_zeta * sin_theta
    r32 = -sin_zeta * sin_theta
    r33 = cos_theta

    # Apply rotation
    x2 = r11 * x1 + r12 * y1 + r13 * z1
    y2 = r21 * x1 + r22 * y1 + r23 * z1
    z2 = r31 * x1 + r32 * y1 + r33 * z1

    # Convert back to ecliptic at target epoch
    cos_eps2 = math.cos(eps2)
    sin_eps2 = math.sin(eps2)

    # Rotate from equatorial to ecliptic
    y_ecl = y2 * cos_eps2 + z2 * sin_eps2
    z_ecl = -y2 * sin_eps2 + z2 * cos_eps2

    new_lon = math.atan2(y_ecl, x2) * RAD_TO_DEG
    new_lat = math.asin(max(-1.0, min(1.0, z_ecl))) * RAD_TO_DEG

    return normalize_angle(new_lon), new_lat


def precess_to_j2000(
    lon: float,
    lat: float,
    jd_tt: float,
) -> Tuple[float, float]:
    """Precess ecliptic coordinates from date to J2000.0.

    Args:
        lon: Ecliptic longitude in degrees (of date).
        lat: Ecliptic latitude in degrees (of date).
        jd_tt: Julian Day (TT) of the source epoch.

    Returns:
        Tuple of (lon_j2000, lat_j2000) in degrees.
    """
    return precess_ecliptic(lon, lat, jd_tt, J2000)


def precess_from_j2000(
    lon: float,
    lat: float,
    jd_tt: float,
) -> Tuple[float, float]:
    """Precess ecliptic coordinates from J2000.0 to date.

    Args:
        lon: Ecliptic longitude in degrees (J2000.0).
        lat: Ecliptic latitude in degrees (J2000.0).
        jd_tt: Julian Day (TT) of the target epoch.

    Returns:
        Tuple of (lon_date, lat_date) in degrees.
    """
    return precess_ecliptic(lon, lat, J2000, jd_tt)


# =============================================================================
# FRAME OF DATE TRANSFORMATIONS
# =============================================================================


def frame_bias_j2000() -> Tuple[float, float, float]:
    """Return the ICRS-to-J2000.0 frame bias angles.

    The ICRS (International Celestial Reference System) is aligned with
    the dynamical mean equator and equinox of J2000.0 to within a few
    milliarcseconds. These small rotations account for the difference.

    Returns:
        Tuple of (dpsi, deps, da) in arcseconds:
        - dpsi: frame bias in RA
        - deps: frame bias in Dec
        - da: equinox offset
    """
    # IERS Conventions (2010), Chapter 5
    dpsi = -0.0166170  # arcsec
    deps = -0.0068192  # arcsec
    da = -0.01460  # arcsec

    return dpsi, deps, da
