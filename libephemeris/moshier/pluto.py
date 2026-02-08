"""
Pluto ephemeris for Moshier calculations.

Implements the Moshier analytical theory for Pluto based on a trigonometric
fit to JPL DE404 ephemeris. This provides heliocentric ecliptic coordinates
of Pluto valid from approximately -3000 CE to +3000 CE with accuracy of
about 1-3 arcseconds.

The algorithm evaluates Poisson series (trigonometric terms multiplied by
powers of time) for longitude, latitude, and radius. The fundamental
arguments are the mean longitudes of the major planets.

This module has NO dependencies on Skyfield or SPK files.
Only numpy is used for numerical operations.

References:
- Swiss Ephemeris swemplan.c (Astrodienst AG)
- Moshier, S.L., DE404 ephemeris fit
"""

from __future__ import annotations

import math
from typing import List, Tuple

from .pluto_data import (
    FREQS,
    PHASES,
    PLU404_DISTANCE,
    PLU404_MAX_HARMONIC,
    PLUTO_ARG_TABLE,
    PLUTO_LAT_TABLE,
    PLUTO_LON_TABLE,
    PLUTO_RAD_TABLE,
)
from .utils import (
    DEG_TO_RAD,
    J2000,
    RAD_TO_DEG,
    cartesian_to_spherical,
    normalize_angle,
)
from .vsop87 import calc_earth_heliocentric

# =============================================================================
# CONSTANTS
# =============================================================================

# Swiss Ephemeris Pluto body ID
MOSHIER_PLUTO = 9

# Time scale: 10000 Julian years in days
TIMESCALE = 3652500.0

# Arcseconds to radians
STR = 4.8481368110953599359e-6  # arcseconds to radians

# Modulo for 360 degrees in arcseconds
ARCSEC_360 = 1296000.0


def _mods3600(x: float) -> float:
    """Reduce angle in arcseconds to range [0, 360 degrees)."""
    return x - ARCSEC_360 * math.floor(x / ARCSEC_360)


# =============================================================================
# SINE/COSINE TABLE FOR MULTIPLE ANGLES
# =============================================================================


def _sscc(arg: float, n: int) -> Tuple[List[float], List[float]]:
    """Compute sin and cos for multiple angles.

    Builds tables of sin(k*arg) and cos(k*arg) for k = 1 to n.

    Args:
        arg: Fundamental angle in radians.
        n: Maximum harmonic number.

    Returns:
        Tuple of (sin_table, cos_table) where each is a list of length n.
    """
    if n <= 0:
        return [], []

    ss = [0.0] * n
    cc = [0.0] * n

    su = math.sin(arg)
    cu = math.cos(arg)
    ss[0] = su  # sin(L)
    cc[0] = cu  # cos(L)

    if n >= 2:
        # sin(2L) = 2*sin(L)*cos(L)
        # cos(2L) = cos^2(L) - sin^2(L)
        sv = 2.0 * su * cu
        cv = cu * cu - su * su
        ss[1] = sv
        cc[1] = cv

        # Recurrence for higher harmonics
        for i in range(2, n):
            s = su * cv + cu * sv
            cv = cu * cv - su * sv
            sv = s
            ss[i] = sv
            cc[i] = cv

    return ss, cc


# =============================================================================
# MAIN PLUTO CALCULATION
# =============================================================================


def _calc_pluto_moshier(jd_tt: float) -> Tuple[float, float, float]:
    """Calculate heliocentric ecliptic coordinates of Pluto using Moshier theory.

    This implements the algorithm from Swiss Ephemeris swemplan.c, evaluating
    the Poisson series for Pluto's longitude, latitude, and radius.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (longitude, latitude, radius) where:
        - longitude is in radians (ecliptic J2000.0)
        - latitude is in radians (ecliptic J2000.0)
        - radius is in AU
    """
    # Time in units of 10000 Julian years from J2000.0
    T = (jd_tt - J2000) / TIMESCALE

    # Build sin/cos tables for each planet's mean longitude
    ss_all: List[List[float]] = []
    cc_all: List[List[float]] = []

    for i in range(9):
        max_harm = PLU404_MAX_HARMONIC[i]
        if max_harm > 0:
            # Mean longitude at time T
            sr = (_mods3600(FREQS[i] * T) + PHASES[i]) * STR
            ss_tab, cc_tab = _sscc(sr, max_harm)
            ss_all.append(ss_tab)
            cc_all.append(cc_tab)
        else:
            ss_all.append([])
            cc_all.append([])

    # Pointers into the data tables
    p_idx = 0  # Index into argument table
    pl_idx = 0  # Index into longitude table
    pb_idx = 0  # Index into latitude table
    pr_idx = 0  # Index into radius table

    # Accumulated sums (in arcseconds for lon/lat, relative units for radius)
    sl = 0.0
    sb = 0.0
    sr = 0.0

    # Process the argument table
    args = PLUTO_ARG_TABLE
    lon_tbl = PLUTO_LON_TABLE
    lat_tbl = PLUTO_LAT_TABLE
    rad_tbl = PLUTO_RAD_TABLE

    while True:
        # Number of periodic arguments
        np = args[p_idx]
        p_idx += 1

        if np < 0:
            # End of table
            break

        if np == 0:
            # Polynomial term
            nt = args[p_idx]
            p_idx += 1

            # Longitude polynomial
            cu = lon_tbl[pl_idx]
            pl_idx += 1
            for _ in range(nt):
                cu = cu * T + lon_tbl[pl_idx]
                pl_idx += 1
            sl += _mods3600(cu)

            # Latitude polynomial
            cu = lat_tbl[pb_idx]
            pb_idx += 1
            for _ in range(nt):
                cu = cu * T + lat_tbl[pb_idx]
                pb_idx += 1
            sb += cu

            # Radius polynomial
            cu = rad_tbl[pr_idx]
            pr_idx += 1
            for _ in range(nt):
                cu = cu * T + rad_tbl[pr_idx]
                pr_idx += 1
            sr += cu

            continue

        # Trigonometric term: build the argument as sum of planet harmonics
        k1 = 0
        cv = 0.0
        sv = 0.0

        for _ in range(np):
            # Harmonic number
            j = args[p_idx]
            p_idx += 1
            # Planet index (1-based in table)
            m = args[p_idx] - 1
            p_idx += 1

            if j != 0:
                k = abs(j)
                k -= 1  # Convert to 0-based index

                # Get sin(k*angle) and cos(k*angle) for planet m
                if k < len(ss_all[m]):
                    su = ss_all[m][k]
                    if j < 0:
                        su = -su
                    cu = cc_all[m][k]
                else:
                    # Handle case where harmonic exceeds table
                    continue

                if k1 == 0:
                    # First angle
                    sv = su
                    cv = cu
                    k1 = 1
                else:
                    # Combine angles: sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
                    t = su * cv + cu * sv
                    cv = cu * cv - su * sv
                    sv = t

        # Highest power of T
        nt = args[p_idx]
        p_idx += 1

        # Longitude contribution
        cu = lon_tbl[pl_idx]
        pl_idx += 1
        su = lon_tbl[pl_idx]
        pl_idx += 1
        for _ in range(nt):
            cu = cu * T + lon_tbl[pl_idx]
            pl_idx += 1
            su = su * T + lon_tbl[pl_idx]
            pl_idx += 1
        sl += cu * cv + su * sv

        # Latitude contribution
        cu = lat_tbl[pb_idx]
        pb_idx += 1
        su = lat_tbl[pb_idx]
        pb_idx += 1
        for _ in range(nt):
            cu = cu * T + lat_tbl[pb_idx]
            pb_idx += 1
            su = su * T + lat_tbl[pb_idx]
            pb_idx += 1
        sb += cu * cv + su * sv

        # Radius contribution
        cu = rad_tbl[pr_idx]
        pr_idx += 1
        su = rad_tbl[pr_idx]
        pr_idx += 1
        for _ in range(nt):
            cu = cu * T + rad_tbl[pr_idx]
            pr_idx += 1
            su = su * T + rad_tbl[pr_idx]
            pr_idx += 1
        sr += cu * cv + su * sv

    # Convert results
    # Longitude and latitude: arcseconds to radians
    lon_rad = STR * sl
    lat_rad = STR * sb
    # Radius: scale factor applied
    radius = STR * PLU404_DISTANCE * sr + PLU404_DISTANCE

    return lon_rad, lat_rad, radius


def calc_pluto_heliocentric(jd_tt: float) -> Tuple[float, float, float]:
    """Calculate heliocentric ecliptic coordinates of Pluto.

    Uses the Moshier DE404 fit for the full range -3000 to +3000 CE.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (longitude, latitude, radius) in degrees and AU.
    """
    lon_rad, lat_rad, radius = _calc_pluto_moshier(jd_tt)

    # Convert to degrees
    lon_deg = normalize_angle(lon_rad * RAD_TO_DEG)
    lat_deg = lat_rad * RAD_TO_DEG

    return lon_deg, lat_deg, radius


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
