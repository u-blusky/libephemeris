"""
Fast calculation pipeline using precomputed .leb binary ephemeris data.

This module reimplements the Skyfield pipeline using LEBReader as the data
source, handling coordinate transforms, light-time correction, aberration,
and flag dispatch.

Three pipelines:
    A: ICRS barycentric bodies (Sun-Pluto, Earth, Chiron, Ceres-Vesta)
    B: Ecliptic direct bodies (lunar nodes, Lilith variants)
    C: Heliocentric bodies (Uranians, Transpluto)
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Dict, Optional, Tuple

from .constants import (
    SE_EARTH,
    SE_MOON,
    SE_SUN,
    SEFLG_BARYCTR,
    SEFLG_EQUATORIAL,
    SEFLG_HELCTR,
    SEFLG_ICRS,
    SEFLG_J2000,
    SEFLG_MOSEPH,
    SEFLG_NOABERR,
    SEFLG_NOGDEFL,
    SEFLG_NONUT,
    SEFLG_RADIANS,
    SEFLG_SIDEREAL,
    SEFLG_SPEED,
    SEFLG_TOPOCTR,
    SEFLG_TRUEPOS,
    SEFLG_XYZ,
)
from .leb_format import (
    COORD_ECLIPTIC,
    COORD_HELIO_ECL,
    COORD_ICRS_BARY,
    COORD_ICRS_BARY_SYSTEM,
)

if TYPE_CHECKING:
    from .leb_reader import LEBReader

# =============================================================================
# CONSTANTS
# =============================================================================

C_LIGHT_AU_DAY = 173.1446326846693  # Speed of light in AU/day
J2000 = 2451545.0  # J2000.0 epoch in JD
OBLIQUITY_J2000_DEG = 23.4392911  # Mean obliquity at J2000 (degrees)
OBLIQUITY_J2000_RAD = math.radians(OBLIQUITY_J2000_DEG)

# IAU 2006 mean obliquity polynomial coefficients (arcseconds)
# eps = 84381.406 - 46.836769*T - 0.0001831*T^2 + 0.00200340*T^3 ...
_OBLIQUITY_COEFFS = (
    84381.406,
    -46.836769,
    -0.0001831,
    0.00200340,
    -0.000000576,
    -0.0000000434,
)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================


def _mean_obliquity_iau2006(jd_tt: float) -> float:
    """Compute IAU 2006 mean obliquity in degrees.

    Args:
        jd_tt: Julian Day in TT.

    Returns:
        Mean obliquity in degrees.
    """
    T = (jd_tt - J2000) / 36525.0
    eps_arcsec = _OBLIQUITY_COEFFS[0]
    T_power = T
    for i in range(1, len(_OBLIQUITY_COEFFS)):
        eps_arcsec += _OBLIQUITY_COEFFS[i] * T_power
        T_power *= T
    return eps_arcsec / 3600.0


def _vec3_sub(a: Tuple[float, ...], b: Tuple[float, ...]) -> Tuple[float, float, float]:
    """Subtract two 3-vectors."""
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec3_dist(v: Tuple[float, float, float]) -> float:
    """Euclidean distance of a 3-vector."""
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def _cartesian_to_spherical(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Convert Cartesian to spherical (lon, lat, dist) in degrees and AU.

    Returns:
        (longitude_deg, latitude_deg, distance_au)
    """
    dist = math.sqrt(x * x + y * y + z * z)
    if dist == 0.0:
        return (0.0, 0.0, 0.0)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, z / dist))))
    return (lon, lat, dist)


def _cartesian_velocity_to_spherical(
    x: float,
    y: float,
    z: float,
    vx: float,
    vy: float,
    vz: float,
) -> Tuple[float, float, float]:
    """Convert Cartesian velocity to spherical velocity.

    Given position (x, y, z) and velocity (vx, vy, vz), compute the time
    derivatives of (longitude, latitude, distance) using the standard
    Cartesian-to-spherical Jacobian.

    Args:
        x, y, z: Position in AU (Cartesian).
        vx, vy, vz: Velocity in AU/day (Cartesian).

    Returns:
        (dlon_deg_day, dlat_deg_day, ddist_au_day)
    """
    r_xy_sq = x * x + y * y
    r_sq = r_xy_sq + z * z
    r = math.sqrt(r_sq)
    r_xy = math.sqrt(r_xy_sq)

    if r == 0.0 or r_xy == 0.0:
        # At the pole or origin — angular rates undefined
        ddist = math.sqrt(vx * vx + vy * vy + vz * vz) if r == 0.0 else 0.0
        return (0.0, 0.0, ddist)

    dlon_rad = (x * vy - y * vx) / r_xy_sq  # rad/day
    dlat_rad = (vz * r_xy_sq - z * (x * vx + y * vy)) / (r_sq * r_xy)  # rad/day
    ddist = (x * vx + y * vy + z * vz) / r  # AU/day

    dlon_deg = math.degrees(dlon_rad)  # deg/day
    dlat_deg = math.degrees(dlat_rad)  # deg/day

    return (dlon_deg, dlat_deg, ddist)


def _rotate_equatorial_to_ecliptic(
    x: float, y: float, z: float, eps_rad: float
) -> Tuple[float, float, float]:
    """Rotate from equatorial to ecliptic frame.

    Args:
        x, y, z: Equatorial Cartesian coordinates.
        eps_rad: Obliquity in radians.

    Returns:
        (x_ecl, y_ecl, z_ecl) in ecliptic frame.
    """
    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    return (
        x,
        y * cos_eps + z * sin_eps,
        -y * sin_eps + z * cos_eps,
    )


def _rotate_icrs_to_ecliptic_j2000(
    x: float, y: float, z: float
) -> Tuple[float, float, float]:
    """Rotate ICRS (≈equatorial J2000) to ecliptic J2000."""
    return _rotate_equatorial_to_ecliptic(x, y, z, OBLIQUITY_J2000_RAD)


def _apply_aberration(
    geo: Tuple[float, float, float],
    earth_vel: Tuple[float, float, float],
    light_time: float = 0.0,
) -> Tuple[float, float, float]:
    """Apply relativistic aberration to a geometric position vector.

    Uses the full special-relativistic formula matching Skyfield's
    ``add_aberration()``.  This includes the Lorentz factor γ for
    exact agreement with the reference pipeline (< 0.001 mas).

    When *light_time* is zero or negative, falls back to a first-order
    approximation (classical Bradley aberration) for backward
    compatibility with callers that do not track light-time.

    Args:
        geo: Light-time-corrected geocentric position vector (ICRS, AU).
        earth_vel: Earth barycentric velocity vector (ICRS, AU/day).
        light_time: One-way light travel time in days (observer→target).

    Returns:
        Aberrated position vector (same frame, AU).
    """
    dist = _vec3_dist(geo)
    if dist == 0.0:
        return geo

    # ── Full special-relativistic aberration (Skyfield formula) ──────────
    if light_time > 0.0:
        p1mag = light_time * C_LIGHT_AU_DAY  # distance in AU
        vemag = math.sqrt(earth_vel[0] ** 2 + earth_vel[1] ** 2 + earth_vel[2] ** 2)
        if vemag == 0.0 or p1mag == 0.0:
            return geo

        beta = vemag / C_LIGHT_AU_DAY
        dot = geo[0] * earth_vel[0] + geo[1] * earth_vel[1] + geo[2] * earth_vel[2]
        cosd = dot / (p1mag * vemag)
        gammai = math.sqrt(1.0 - beta * beta)  # inverse Lorentz factor
        p = beta * cosd
        q = (1.0 + p / (1.0 + gammai)) * light_time
        r = 1.0 + p
        if abs(r) < 1e-30:
            r = 1.0

        return (
            (gammai * geo[0] + q * earth_vel[0]) / r,
            (gammai * geo[1] + q * earth_vel[1]) / r,
            (gammai * geo[2] + q * earth_vel[2]) / r,
        )

    # ── Fallback: first-order Bradley aberration ────────────────────────
    ux = geo[0] / dist
    uy = geo[1] / dist
    uz = geo[2] / dist

    vx = earth_vel[0] / C_LIGHT_AU_DAY
    vy = earth_vel[1] / C_LIGHT_AU_DAY
    vz = earth_vel[2] / C_LIGHT_AU_DAY

    dot = ux * vx + uy * vy + uz * vz

    ax = ux + vx - ux * dot
    ay = uy + vy - uy * dot
    az = uz + vz - uz * dot

    a_dist = math.sqrt(ax * ax + ay * ay + az * az)
    if a_dist == 0.0:
        return geo

    return (ax / a_dist * dist, ay / a_dist * dist, az / a_dist * dist)


# Gravitational deflection constants (matching Skyfield's relativity module)
_GS = 1.32712440017987e20  # heliocentric gravitational constant (m^3/s^2)
_C_MS = 299792458.0  # speed of light (m/s)
_AU_M = 149597870700  # 1 AU in metres

# Deflector reciprocal masses (solar mass / deflector mass)
_DEFLECTORS: Tuple[Tuple[int, float], ...] = (
    (SE_SUN, 1.0),  # Sun
    (5, 1047.3486),  # Jupiter barycenter
    (6, 3497.898),  # Saturn barycenter
)

# Body IDs in the LEB file for deflector barycenters
# (Jupiter=5, Saturn=6 map to body_ids 5, 6 in LEB — their barycentric ICRS)
_DEFLECTOR_LEB_IDS = {SE_SUN: SE_SUN, 5: 5, 6: 6}

# =============================================================================
# CENTER-OF-BODY (COB) CORRECTION FOR SYSTEM BARYCENTERS
# =============================================================================

# Barycenter names for outer planets (used by get_cob_offset / SPK centers)
_SYSTEM_BARY_NAMES: Dict[int, str] = {
    5: "jupiter barycenter",
    6: "saturn barycenter",
    7: "uranus barycenter",
    8: "neptune barycenter",
    9: "pluto barycenter",
}


def _apply_cob_correction(
    pos: Tuple[float, float, float],
    ipl: int,
    jd_tt: float,
) -> Tuple[float, float, float]:
    """Apply center-of-body correction to a system barycenter position.

    Converts system barycenter ICRS position to planet center ICRS position
    by adding the COB offset. Uses SPK planet_centers segments where available
    (high precision), falling back to analytical moon-theory COB corrections.

    This matches the behavior of _SpkCenterTarget in planets.py.

    Args:
        pos: System barycenter ICRS position (x, y, z) in AU.
        ipl: Body ID (5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune, 9=Pluto).
        jd_tt: Julian Day TT at which to evaluate the COB offset.

    Returns:
        Planet center ICRS position (x, y, z) in AU.
    """
    from .planets import _PLANET_CENTER_NAIF_IDS
    from .state import get_planet_center_segment, get_timescale

    bary_name = _SYSTEM_BARY_NAMES.get(ipl)
    if bary_name is None:
        return pos  # Not an outer planet, no COB needed

    ts = get_timescale()
    t = ts.tt_jd(jd_tt)

    # Map body_id to planet name for NAIF lookup
    planet_name = {5: "jupiter", 6: "saturn", 7: "uranus", 8: "neptune", 9: "pluto"}[
        ipl
    ]

    # Try SPK center offset first (high precision)
    if planet_name in _PLANET_CENTER_NAIF_IDS:
        naif_id = _PLANET_CENTER_NAIF_IDS[planet_name]
        seg = get_planet_center_segment(naif_id)
        if seg is not None:
            try:
                offset_pos = seg.at(t).position.au
                return (
                    pos[0] + float(offset_pos[0]),
                    pos[1] + float(offset_pos[1]),
                    pos[2] + float(offset_pos[2]),
                )
            except Exception:
                pass  # Fall through to analytical COB

    # Fallback: analytical COB from moon theories
    from .moon_theories import get_cob_offset

    offset = get_cob_offset(bary_name, t)
    return (
        pos[0] + offset[0],
        pos[1] + offset[1],
        pos[2] + offset[2],
    )


def _apply_gravitational_deflection(
    geo: Tuple[float, float, float],
    earth_bary: Tuple[float, float, float],
    jd_tt: float,
    light_time: float,
    reader: "LEBReader",
) -> Tuple[float, float, float]:
    """Apply PPN gravitational light deflection by Sun, Jupiter, Saturn.

    Matches Skyfield's ``apparent(deflectors=(10, 599, 699))`` formula.
    Uses LEB data for deflector positions (barycentric ICRS).

    For the Sun, deflection is the dominant correction (~max 1.75" at
    the limb, typically 0.01–4" for planets).  Jupiter and Saturn add
    ~0.01" near their limbs.

    Args:
        geo: Light-time-corrected geocentric ICRS position (AU).
        earth_bary: Earth ICRS barycentric position (AU).
        jd_tt: Observation Julian Day in TT.
        light_time: Light travel time to target (days).
        reader: Open LEBReader for deflector positions.

    Returns:
        Deflection-corrected geocentric ICRS position (AU).
    """
    result = list(geo)
    pmag = math.sqrt(result[0] ** 2 + result[1] ** 2 + result[2] ** 2)
    if pmag == 0.0:
        return geo

    for defl_body_id, rmass in _DEFLECTORS:
        # 1. Deflector barycentric position at observation time
        try:
            defl_pos, _ = reader.eval_body(defl_body_id, jd_tt)
        except (KeyError, ValueError):
            continue

        # 2. Deflector relative to observer
        gpv = (
            defl_pos[0] - earth_bary[0],
            defl_pos[1] - earth_bary[1],
            defl_pos[2] - earth_bary[2],
        )

        # 3. Unit vector to target
        phat = (result[0] / pmag, result[1] / pmag, result[2] / pmag)

        # 4. Light-time difference: when photon passed closest to deflector
        dlt = (phat[0] * gpv[0] + phat[1] * gpv[1] + phat[2] * gpv[2]) / C_LIGHT_AU_DAY

        # 5. Clamp and compute time at closest approach
        tclose_offset = max(0.0, min(dlt, light_time))
        tclose_jd = jd_tt - tclose_offset

        # 6. Deflector position at closest approach
        try:
            defl_close, _ = reader.eval_body(defl_body_id, tclose_jd)
        except (KeyError, ValueError):
            continue

        # 7. pe = observer - deflector (observer relative to deflector)
        pe = (
            earth_bary[0] - defl_close[0],
            earth_bary[1] - defl_close[1],
            earth_bary[2] - defl_close[2],
        )

        # 8. pq = target relative to deflector (from observer frame)
        pq = (result[0] + pe[0], result[1] + pe[1], result[2] + pe[2])

        qmag = math.sqrt(pq[0] ** 2 + pq[1] ** 2 + pq[2] ** 2)
        emag = math.sqrt(pe[0] ** 2 + pe[1] ** 2 + pe[2] ** 2)

        if qmag == 0.0 or emag == 0.0:
            continue

        qhat = (pq[0] / qmag, pq[1] / qmag, pq[2] / qmag)
        ehat = (pe[0] / emag, pe[1] / emag, pe[2] / emag)

        pdotq = phat[0] * qhat[0] + phat[1] * qhat[1] + phat[2] * qhat[2]
        qdote = qhat[0] * ehat[0] + qhat[1] * ehat[1] + qhat[2] * ehat[2]
        edotp = ehat[0] * phat[0] + ehat[1] * phat[1] + ehat[2] * phat[2]

        # Skip if object is on the line-of-sight to the deflector
        if abs(edotp) > 0.99999999999:
            continue

        fac1 = 2.0 * _GS / (_C_MS * _C_MS * emag * _AU_M * rmass)
        fac2 = 1.0 + qdote
        if abs(fac2) < 1e-30:
            continue

        coeff = fac1 / fac2 * pmag
        result[0] += coeff * (pdotq * ehat[0] - edotp * qhat[0])
        result[1] += coeff * (pdotq * ehat[1] - edotp * qhat[1])
        result[2] += coeff * (pdotq * ehat[2] - edotp * qhat[2])

    return (result[0], result[1], result[2])


def _get_skyfield_frame_data(
    jd_tt: float,
) -> Tuple[
    Tuple[Tuple[float, float, float], ...],
    float,
    float,
    float,
]:
    """Get precession-nutation matrix and nutation angles from Skyfield.

    Uses Skyfield's IAU 2000A nutation model, matching the reference pipeline
    in planets.py.  This ensures LEB-vs-Skyfield comparison tests measure only
    Chebyshev fitting error, not nutation model mismatch.

    Args:
        jd_tt: Julian Day in TT.

    Returns:
        (pn_mat, dpsi, deps, eps_true_rad) where:
        - pn_mat: 3x3 ICRS→true-equatorial-of-date rotation as nested tuples
        - dpsi: nutation in longitude (radians, IAU 2000A)
        - deps: nutation in obliquity (radians, IAU 2000A)
        - eps_true_rad: true obliquity (radians)
    """
    from .state import get_timescale

    ts = get_timescale()
    t = ts.tt_jd(jd_tt)

    # PNM matrix: ICRS -> true equatorial of date (N × P × B)
    M = t.M
    pn_mat = (
        (float(M[0][0]), float(M[0][1]), float(M[0][2])),
        (float(M[1][0]), float(M[1][1]), float(M[1][2])),
        (float(M[2][0]), float(M[2][1]), float(M[2][2])),
    )

    # Nutation angles (IAU 2000A, matching Skyfield's own computation)
    dpsi, deps = t._nutation_angles_radians
    eps_true_rad = float(t._mean_obliquity_radians + deps)

    return pn_mat, float(dpsi), float(deps), eps_true_rad


def _mat3_vec3(
    mat: Tuple[Tuple[float, float, float], ...],
    vec: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """Multiply a 3x3 matrix by a 3-vector."""
    return (
        mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2],
        mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2],
        mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2],
    )


def _cotrans(lon: float, lat: float, eps: float) -> Tuple[float, float]:
    """Coordinate transform between ecliptic and equatorial.

    Negative eps: ecliptic -> equatorial.
    Positive eps: equatorial -> ecliptic.

    Args:
        lon: Longitude/RA in degrees.
        lat: Latitude/Dec in degrees.
        eps: Obliquity in degrees (sign determines direction).

    Returns:
        (transformed_lon, transformed_lat) in degrees.
    """
    eps_rad = math.radians(-eps)  # Negate to match pyswisseph convention
    lon_rad = math.radians(lon)
    lat_rad = math.radians(lat)

    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # Spherical rotation
    x = cos_lat * cos_lon
    y = cos_lat * sin_lon * cos_eps - sin_lat * sin_eps
    z = cos_lat * sin_lon * sin_eps + sin_lat * cos_eps

    new_lon = math.degrees(math.atan2(y, x)) % 360.0
    r = math.sqrt(x * x + y * y + z * z)
    new_lat = math.degrees(math.asin(max(-1.0, min(1.0, z / r)))) if r > 0 else 0.0

    return new_lon, new_lat


def _precess_ecliptic(
    lon: float, lat: float, from_jd: float, to_jd: float
) -> Tuple[float, float]:
    """Precess ecliptic coordinates between two epochs.

    Uses the astrometry module's precession.

    Args:
        lon, lat: Ecliptic coordinates in degrees.
        from_jd, to_jd: Julian Days (TT).

    Returns:
        (precessed_lon, precessed_lat) in degrees.
    """
    from .astrometry import _precess_ecliptic as _pe

    return _pe(lon, lat, from_jd, to_jd)


# =============================================================================
# AYANAMSHA
# =============================================================================

# IAU 2006 general precession polynomial (arcsec/century)
_PREC_COEFFS = (5028.796195, 1.1054348, 0.00007964, -0.000023857, -0.0000000383)

# Ayanamsha J2000 offsets for formula-based sidereal modes (degrees).
# Star-based / galactic modes have (0.0) placeholders and require Skyfield.
# These values mirror the ayanamsha_data dict in planets._calc_ayanamsa().
_AYANAMSHA_J2000: Dict[int, float] = {
    0: 24.740300,  # SE_SIDM_FAGAN_BRADLEY
    1: 23.857092,  # SE_SIDM_LAHIRI
    2: 27.815753,  # SE_SIDM_DELUCE
    3: 22.410791,  # SE_SIDM_RAMAN
    4: 20.057541,  # SE_SIDM_USHASHASHI
    5: 23.760240,  # SE_SIDM_KRISHNAMURTI
    6: 28.359679,  # SE_SIDM_DJWHAL_KHUL
    7: 22.478803,  # SE_SIDM_YUKTESHWAR
    8: 22.762137,  # SE_SIDM_JN_BHASIN
    9: 23.533640,  # SE_SIDM_BABYL_KUGLER1
    10: 24.933640,  # SE_SIDM_BABYL_KUGLER2
    11: 25.783640,  # SE_SIDM_BABYL_KUGLER3
    12: 24.733640,  # SE_SIDM_BABYL_HUBER
    13: 24.522528,  # SE_SIDM_BABYL_ETPSC
    14: 24.758924,  # SE_SIDM_ALDEBARAN_15TAU
    15: 20.247788,  # SE_SIDM_HIPPARCHOS
    16: 19.992959,  # SE_SIDM_SASSANIAN
    18: 0.0,  # SE_SIDM_J2000
    19: 1.396581,  # SE_SIDM_J1900
    20: 0.698370,  # SE_SIDM_B1950
    21: 20.895059,  # SE_SIDM_SURYASIDDHANTA
    22: 20.680425,  # SE_SIDM_SURYASIDDHANTA_MSUN
    23: 20.895060,  # SE_SIDM_ARYABHATA
    24: 20.657427,  # SE_SIDM_ARYABHATA_MSUN
    25: 20.103388,  # SE_SIDM_SS_REVATI
    26: 23.005763,  # SE_SIDM_SS_CITRA
    37: 20.575847,  # SE_SIDM_ARYABHATA_522
    38: 24.615753,  # SE_SIDM_BABYL_BRITTON
    41: 25.000019,  # SE_SIDM_GALEQU_FIORENZA
}

# Star-based modes that cannot be computed without Skyfield
_STAR_BASED_MODES = frozenset({17, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 39, 40, 42})


def _calc_ayanamsa_from_leb(
    reader: "LEBReader",
    jd_tt: float,
    sid_mode: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
) -> float:
    """Compute ayanamsa using only .leb data (no Skyfield).

    For formula-based ayanamsha modes only. Star-based and galactic modes
    require Skyfield and will raise KeyError (triggering fallback).

    Args:
        reader: Open LEBReader instance.
        jd_tt: Julian Day in TT.
        sid_mode: Sidereal mode ID (if None, reads from global state).
        sid_t0: Reference epoch JD for custom ayanamsha.
        sid_ayan_t0: Ayanamsha value at reference epoch (degrees).

    Returns:
        Ayanamsa in degrees.

    Raises:
        KeyError: If the active sidereal mode requires Skyfield.
    """
    if sid_mode is None:
        from .state import _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0

        sid_mode = _SIDEREAL_MODE
        sid_t0 = _SIDEREAL_T0
        sid_ayan_t0 = _SIDEREAL_AYAN_T0

    T = (jd_tt - J2000) / 36525.0

    # IAU 2006 general precession
    precession_arcsec = 0.0
    T_power = T
    for coeff in _PREC_COEFFS:
        precession_arcsec += coeff * T_power
        T_power *= T

    # Get reference offset for active sidereal mode
    mode = sid_mode if sid_mode is not None else 1  # Default Lahiri

    if mode in _STAR_BASED_MODES:
        raise KeyError(f"Star-based sidereal mode {mode} requires Skyfield")

    if mode == 255:
        # SE_SIDM_USER: custom user-defined ayanamsha
        ayan_t0 = sid_ayan_t0 if sid_ayan_t0 is not None else 0.0
        t0 = sid_t0 if sid_t0 is not None else J2000
        T0 = (t0 - J2000) / 36525.0

        # Delta precession from user epoch to current epoch
        def _prec(Tc: float) -> float:
            return sum(c * Tc ** (i + 1) for i, c in enumerate(_PREC_COEFFS))

        delta_prec = _prec(T) - _prec(T0)
        mean_aya = ayan_t0 + delta_prec / 3600.0
    elif mode in _AYANAMSHA_J2000:
        aya_j2000 = _AYANAMSHA_J2000[mode]
        mean_aya = aya_j2000 + precession_arcsec / 3600.0
    else:
        raise KeyError(f"Unknown sidereal mode {mode}")

    # Mean ayanamsa (without nutation) - swe_get_ayanamsa_ut() returns mean value
    return mean_aya % 360.0


# =============================================================================
# PIPELINE A: ICRS BARYCENTRIC BODIES
# =============================================================================


def _pipeline_icrs(
    reader: "LEBReader",
    jd_tt: float,
    ipl: int,
    iflag: int,
    want_velocity: bool = False,
    is_system_bary: bool = False,
) -> Tuple[float, ...]:
    """Pipeline A: compute ecliptic coordinates for ICRS barycentric bodies.

    Handles both planet-center bodies (COORD_ICRS_BARY) and system-barycenter
    bodies (COORD_ICRS_BARY_SYSTEM). For system barycenters, COB correction
    is applied at runtime to match Skyfield's _SpkCenterTarget behavior.

    When want_velocity is False (default), returns (lon, lat, dist).
    When want_velocity is True, returns (lon, lat, dist, dlon, dlat, ddist)
    where the velocity components are analytically derived from the Chebyshev
    polynomial derivatives, transformed through the same rotation matrices
    as the position.

    Args:
        reader: Open LEBReader.
        jd_tt: Julian Day TT.
        ipl: Body ID.
        iflag: Flags.
        want_velocity: Whether to compute velocity.
        is_system_bary: If True, stored data is system barycenter; apply COB.

    Returns:
        (lon_deg, lat_deg, dist_au) or
        (lon_deg, lat_deg, dist_au, dlon_deg_day, dlat_deg_day, ddist_au_day)
    """
    # 1. Get body position (and velocity if needed)
    target_pos, target_vel = reader.eval_body(ipl, jd_tt)

    # 1b. For system barycenters, apply COB only for TRUEPOS (no light-time).
    #     For normal path, COB is deferred until after light-time iteration
    #     to match Skyfield's _SpkCenterTarget._observe_from_bcrs() behavior:
    #     iterate light-time on barycenter, apply COB once at retarded time.
    if is_system_bary and (iflag & SEFLG_TRUEPOS):
        target_pos = _apply_cob_correction(target_pos, ipl, jd_tt)

    # Pre-initialize velocity variables to satisfy type checker.
    # These are always set before use when want_velocity=True.
    geo_vel: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    dlon = dlat = ddist = 0.0

    # 2. Observer selection (defer Earth fetch for helio/bary)
    if iflag & SEFLG_HELCTR:
        sun_pos, sun_vel = reader.eval_body(SE_SUN, jd_tt)
        observer = sun_pos
        observer_vel = sun_vel
        earth_vel = (0.0, 0.0, 0.0)  # Not needed (aberration skipped)
    elif iflag & SEFLG_BARYCTR:
        observer = (0.0, 0.0, 0.0)
        observer_vel = (0.0, 0.0, 0.0)
        earth_vel = (0.0, 0.0, 0.0)  # Not needed (aberration skipped)
    else:
        # Geocentric (default) — need Earth position and velocity
        earth_pos, earth_vel = reader.eval_body(SE_EARTH, jd_tt)
        observer = earth_pos
        observer_vel = earth_vel

    # 3. Geometric vector
    geo = _vec3_sub(target_pos, observer)

    # 4. Light-time correction (unless SEFLG_TRUEPOS)
    #    For system barycenters: iterate on raw barycenter positions (smooth),
    #    then apply COB once after convergence at OBSERVER time. This matches
    #    Skyfield's _SpkCenterTarget._observe_from_bcrs() which:
    #      1. Calls barycenter._observe_from_bcrs(observer) to iterate light-time
    #         on the barycenter, returning (pos, vel, t, light_time) where
    #         t = observer.t (the OBSERVATION time, not retarded time).
    #      2. Evaluates center_segment.at(t) at OBSERVER time, not retarded time.
    #    We must match this: COB offset evaluated at jd_tt (observer time).
    retarded_vel = target_vel
    lt = 0.0
    if not (iflag & SEFLG_TRUEPOS):
        for _ in range(3):  # Fixed-point iterations
            dist = _vec3_dist(geo)
            if dist == 0.0:
                break
            lt = dist / C_LIGHT_AU_DAY
            retarded_pos, retarded_vel = reader.eval_body(ipl, jd_tt - lt)
            geo = _vec3_sub(retarded_pos, observer)

        # Apply COB correction at OBSERVER time (jd_tt), matching Skyfield's
        # _SpkCenterTarget._observe_from_bcrs() which evaluates the center
        # segment at observer.t, not at the retarded time.
        if is_system_bary and lt > 0.0:
            retarded_pos_cob = _apply_cob_correction(
                (geo[0] + observer[0], geo[1] + observer[1], geo[2] + observer[2]),
                ipl,
                jd_tt,
            )
            geo = _vec3_sub(retarded_pos_cob, observer)

    if want_velocity:
        geo_vel = _vec3_sub(retarded_vel, observer_vel)

    # 5. Gravitational deflection by Sun, Jupiter, Saturn (PPN formula).
    #    Dominant correction: up to ~4" for Saturn near the Sun's limb.
    #    Skipped for helio/bary/truepos/nogdefl and for the Moon (negligible at
    #    ~0.0026 AU, deflection < 0.000001").
    if not (
        iflag
        & (SEFLG_NOABERR | SEFLG_NOGDEFL | SEFLG_HELCTR | SEFLG_BARYCTR | SEFLG_TRUEPOS)
    ):
        if ipl != SE_MOON and lt > 0.0:
            geo = _apply_gravitational_deflection(geo, observer, jd_tt, lt, reader)

    # 6. Aberration (full special-relativistic, matching Skyfield).
    #    For velocity: the aberration correction depends on Earth's velocity
    #    which changes slowly (~0.017 deg/day²). The velocity component of
    #    aberration is ~1e-8 deg/day — negligible. We skip it.
    if not (iflag & (SEFLG_NOABERR | SEFLG_HELCTR | SEFLG_BARYCTR | SEFLG_TRUEPOS)):
        geo = _apply_aberration(geo, earth_vel, lt)

    # 6. Coordinate transform — apply the same transform to velocity
    if (iflag & SEFLG_EQUATORIAL) and (iflag & SEFLG_J2000):
        # ICRS J2000 equatorial -- geo is already in this frame
        lon_deg, lat_deg, dist = _cartesian_to_spherical(geo[0], geo[1], geo[2])
        if want_velocity:
            dlon, dlat, ddist = _cartesian_velocity_to_spherical(
                geo[0],
                geo[1],
                geo[2],
                geo_vel[0],
                geo_vel[1],
                geo_vel[2],
            )

    elif iflag & SEFLG_EQUATORIAL:
        # True equator of date
        pn_mat, _, _, _ = _get_skyfield_frame_data(jd_tt)
        geo_eq = _mat3_vec3(pn_mat, geo)
        lon_deg, lat_deg, dist = _cartesian_to_spherical(
            geo_eq[0], geo_eq[1], geo_eq[2]
        )
        if want_velocity:
            vel_eq = _mat3_vec3(pn_mat, geo_vel)
            dlon, dlat, ddist = _cartesian_velocity_to_spherical(
                geo_eq[0],
                geo_eq[1],
                geo_eq[2],
                vel_eq[0],
                vel_eq[1],
                vel_eq[2],
            )

    elif iflag & SEFLG_J2000:
        # J2000 ecliptic
        ecl = _rotate_icrs_to_ecliptic_j2000(geo[0], geo[1], geo[2])
        lon_deg, lat_deg, dist = _cartesian_to_spherical(ecl[0], ecl[1], ecl[2])
        if want_velocity:
            vel_ecl = _rotate_icrs_to_ecliptic_j2000(geo_vel[0], geo_vel[1], geo_vel[2])
            dlon, dlat, ddist = _cartesian_velocity_to_spherical(
                ecl[0],
                ecl[1],
                ecl[2],
                vel_ecl[0],
                vel_ecl[1],
                vel_ecl[2],
            )

    else:
        # TRUE ECLIPTIC OF DATE (default) -- most common path
        # Get PNM matrix and true obliquity from Skyfield (IAU 2000A)
        pn_mat, _, _, eps_true_rad = _get_skyfield_frame_data(jd_tt)

        # Step 1: Precess ICRS -> equatorial of date
        geo_eq = _mat3_vec3(pn_mat, geo)

        # Step 2: Rotate equatorial -> ecliptic using true obliquity
        ecl = _rotate_equatorial_to_ecliptic(
            geo_eq[0], geo_eq[1], geo_eq[2], eps_true_rad
        )
        lon_deg, lat_deg, dist = _cartesian_to_spherical(ecl[0], ecl[1], ecl[2])

        if want_velocity:
            vel_eq = _mat3_vec3(pn_mat, geo_vel)
            vel_ecl = _rotate_equatorial_to_ecliptic(
                vel_eq[0], vel_eq[1], vel_eq[2], eps_true_rad
            )
            dlon, dlat, ddist = _cartesian_velocity_to_spherical(
                ecl[0],
                ecl[1],
                ecl[2],
                vel_ecl[0],
                vel_ecl[1],
                vel_ecl[2],
            )

    if want_velocity:
        return lon_deg, lat_deg, dist, dlon, dlat, ddist
    return lon_deg, lat_deg, dist


# =============================================================================
# PIPELINE B: ECLIPTIC DIRECT BODIES
# =============================================================================


def _pipeline_ecliptic(
    reader: "LEBReader",
    jd_tt: float,
    ipl: int,
    iflag: int,
) -> Tuple[float, float, float, float, float, float]:
    """Pipeline B: evaluate ecliptic-direct bodies.

    Returns:
        (lon, lat, dist, dlon, dlat, ddist)
    """
    (lon, lat, dist), (dlon, dlat, ddist) = reader.eval_body(ipl, jd_tt)

    # Coordinate transforms for ecliptic-direct bodies.
    # Input coords are always ecliptic of date.
    if (iflag & SEFLG_EQUATORIAL) and (iflag & SEFLG_J2000):
        # J2000 equatorial: precess ecliptic-of-date -> J2000 ecliptic,
        # then rotate J2000 ecliptic -> J2000 equatorial.
        eps = OBLIQUITY_J2000_DEG

        def _ecl_date_to_eq_j2000(lo: float, la: float) -> tuple[float, float]:
            lo_j, la_j = _precess_ecliptic(lo, la, jd_tt, J2000)
            return _cotrans(lo_j, la_j, -eps)

        # Velocity via finite difference on original ecliptic coords
        dt_step = 0.001  # days
        eq_now_lon, eq_now_lat = _ecl_date_to_eq_j2000(lon, lat)
        eq_fwd_lon, eq_fwd_lat = _ecl_date_to_eq_j2000(
            lon + dlon * dt_step, lat + dlat * dt_step
        )
        d_eq_lon = eq_fwd_lon - eq_now_lon
        if d_eq_lon > 180.0:
            d_eq_lon -= 360.0
        elif d_eq_lon < -180.0:
            d_eq_lon += 360.0
        dlon = d_eq_lon / dt_step
        dlat = (eq_fwd_lat - eq_now_lat) / dt_step
        lon = eq_now_lon
        lat = eq_now_lat

    elif iflag & SEFLG_EQUATORIAL:
        # True equatorial of date: rotate ecliptic-of-date -> equatorial-of-date
        _, _, deps, _ = _get_skyfield_frame_data(jd_tt)
        eps_mean = _mean_obliquity_iau2006(jd_tt)
        eps = eps_mean + math.degrees(deps)

        # Velocity via finite difference on original ecliptic coords
        dt_step = 0.001  # days
        eq_now_lon, eq_now_lat = _cotrans(lon, lat, -eps)
        eq_fwd_lon, eq_fwd_lat = _cotrans(
            lon + dlon * dt_step, lat + dlat * dt_step, -eps
        )
        d_eq_lon = eq_fwd_lon - eq_now_lon
        if d_eq_lon > 180.0:
            d_eq_lon -= 360.0
        elif d_eq_lon < -180.0:
            d_eq_lon += 360.0
        dlon = d_eq_lon / dt_step
        dlat = (eq_fwd_lat - eq_now_lat) / dt_step
        lon = eq_now_lon
        lat = eq_now_lat

    elif iflag & SEFLG_J2000:
        # J2000 ecliptic: precess from ecliptic of date to J2000
        # Velocity must also be transformed via finite difference
        dt_step = 0.001  # days
        j_now_lon, j_now_lat = _precess_ecliptic(lon, lat, jd_tt, J2000)
        j_fwd_lon, j_fwd_lat = _precess_ecliptic(
            lon + dlon * dt_step, lat + dlat * dt_step, jd_tt, J2000
        )
        d_j_lon = j_fwd_lon - j_now_lon
        if d_j_lon > 180.0:
            d_j_lon -= 360.0
        elif d_j_lon < -180.0:
            d_j_lon += 360.0
        dlon = d_j_lon / dt_step
        dlat = (j_fwd_lat - j_now_lat) / dt_step
        lon = j_now_lon
        lat = j_now_lat

    return lon, lat, dist, dlon, dlat, ddist


# =============================================================================
# PIPELINE C: HELIOCENTRIC BODIES
# =============================================================================


def _pipeline_helio(
    reader: "LEBReader",
    jd_tt: float,
    ipl: int,
    iflag: int,
) -> Tuple[float, float, float, float, float, float]:
    """Pipeline C: evaluate heliocentric ecliptic bodies (Uranians, Transpluto).

    Same as Pipeline B -- these bodies are already in heliocentric ecliptic.

    Returns:
        (lon, lat, dist, dlon, dlat, ddist)
    """
    return _pipeline_ecliptic(reader, jd_tt, ipl, iflag)


# =============================================================================
# ENTRY POINTS
# =============================================================================


def fast_calc_ut(
    reader: "LEBReader",
    tjd_ut: float,
    ipl: int,
    iflag: int,
    *,
    sid_mode: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Fast equivalent of swe_calc_ut() using precomputed .leb data.

    Args:
        reader: An open LEBReader instance.
        tjd_ut: Julian Day in Universal Time (UT1).
        ipl: Planet/body ID (SE_SUN, SE_MOON, etc.).
        iflag: Calculation flags (SEFLG_SPEED, SEFLG_HELCTR, etc.).
        sid_mode: Sidereal mode override (for thread-safe context calls).
        sid_t0: Sidereal reference epoch override.
        sid_ayan_t0: Sidereal ayanamsha-at-epoch override.

    Returns:
        Same as swe_calc_ut(): ((lon, lat, dist, dlon, dlat, ddist), iflag)

    Raises:
        KeyError: If body is not in the .leb file (caller should fall back).
        ValueError: If JD is outside the .leb file's range.
    """
    # Flags not supported in LEB mode — fall back to Skyfield
    if iflag & SEFLG_TOPOCTR:
        raise KeyError("SEFLG_TOPOCTR not supported in LEB mode")
    if iflag & SEFLG_XYZ:
        raise KeyError("SEFLG_XYZ not supported in LEB mode")
    if iflag & SEFLG_RADIANS:
        raise KeyError("SEFLG_RADIANS not supported in LEB mode")
    if iflag & SEFLG_NONUT:
        raise KeyError("SEFLG_NONUT not supported in LEB mode")
    if iflag & SEFLG_ICRS:
        raise KeyError("SEFLG_ICRS not supported in LEB mode")

    # Strip SEFLG_MOSEPH (always ignored)
    iflag = iflag & ~SEFLG_MOSEPH

    # Snapshot sidereal state once at entry (thread-safe)
    if sid_mode is None and (iflag & SEFLG_SIDEREAL):
        from .state import _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0

        sid_mode = _SIDEREAL_MODE
        sid_t0 = _SIDEREAL_T0
        sid_ayan_t0 = _SIDEREAL_AYAN_T0

    # Delta-T conversion: UT -> TT
    # Use swe_deltat() for exact match with the Skyfield reference path.
    # The LEB reader's linearly-interpolated sparse table introduces up to
    # ~0.004s error near 1985, which at the Moon's ~0.5"/s speed exceeds
    # the 0.001" target.
    from .time_utils import swe_deltat

    delta_t = swe_deltat(tjd_ut)
    jd_tt = tjd_ut + delta_t

    return _fast_calc_core(
        reader,
        jd_tt,
        tjd_ut,
        ipl,
        iflag,
        sid_mode=sid_mode,
        sid_t0=sid_t0,
        sid_ayan_t0=sid_ayan_t0,
    )


def fast_calc_tt(
    reader: "LEBReader",
    tjd_tt: float,
    ipl: int,
    iflag: int,
    *,
    sid_mode: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Fast equivalent of swe_calc() using precomputed .leb data.

    Args:
        reader: An open LEBReader instance.
        tjd_tt: Julian Day in Terrestrial Time (TT).
        ipl: Planet/body ID.
        iflag: Calculation flags.
        sid_mode: Sidereal mode override (for thread-safe context calls).
        sid_t0: Sidereal reference epoch override.
        sid_ayan_t0: Sidereal ayanamsha-at-epoch override.

    Returns:
        Same as swe_calc(): ((lon, lat, dist, dlon, dlat, ddist), iflag)

    Raises:
        KeyError: If body is not in the .leb file.
        ValueError: If JD is outside the .leb file's range.
    """
    # Flags not supported in LEB mode — fall back to Skyfield
    if iflag & SEFLG_TOPOCTR:
        raise KeyError("SEFLG_TOPOCTR not supported in LEB mode")
    if iflag & SEFLG_XYZ:
        raise KeyError("SEFLG_XYZ not supported in LEB mode")
    if iflag & SEFLG_RADIANS:
        raise KeyError("SEFLG_RADIANS not supported in LEB mode")
    if iflag & SEFLG_NONUT:
        raise KeyError("SEFLG_NONUT not supported in LEB mode")
    if iflag & SEFLG_ICRS:
        raise KeyError("SEFLG_ICRS not supported in LEB mode")

    iflag = iflag & ~SEFLG_MOSEPH

    # Snapshot sidereal state once at entry (thread-safe)
    if sid_mode is None and (iflag & SEFLG_SIDEREAL):
        from .state import _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0

        sid_mode = _SIDEREAL_MODE
        sid_t0 = _SIDEREAL_T0
        sid_ayan_t0 = _SIDEREAL_AYAN_T0

    # For swe_calc, input is already TT
    # We need UT for sidereal ayanamsa, approximate: ut ≈ tt - delta_t(tt)
    tjd_ut = tjd_tt - reader.delta_t(tjd_tt)

    return _fast_calc_core(
        reader,
        tjd_tt,
        tjd_ut,
        ipl,
        iflag,
        sid_mode=sid_mode,
        sid_t0=sid_t0,
        sid_ayan_t0=sid_ayan_t0,
    )


def _fast_calc_core(
    reader: "LEBReader",
    jd_tt: float,
    jd_ut: float,
    ipl: int,
    iflag: int,
    *,
    sid_mode: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Core fast calculation logic shared by fast_calc_ut and fast_calc_tt.

    Args:
        reader: Open LEBReader.
        jd_tt: Julian Day TT (for position computation).
        jd_ut: Julian Day UT (for sidereal ayanamsa).
        ipl: Body ID.
        iflag: Flags.
        sid_mode: Sidereal mode (already snapshot at entry point).
        sid_t0: Sidereal reference epoch.
        sid_ayan_t0: Sidereal ayanamsha at reference epoch.

    Returns:
        ((lon, lat, dist, dlon, dlat, ddist), iflag)
    """
    # Check if body is in the .leb file
    if not reader.has_body(ipl):
        raise KeyError(f"Body {ipl} not in LEB file")

    body = reader._bodies[ipl]

    # Dispatch to appropriate pipeline based on coordinate type
    if body.coord_type == COORD_ICRS_BARY:
        # Pipeline A: ICRS barycentric with analytical velocity
        if iflag & SEFLG_SPEED:
            lon, lat, dist, dlon, dlat, ddist = _pipeline_icrs(
                reader, jd_tt, ipl, iflag, want_velocity=True
            )
        else:
            lon, lat, dist = _pipeline_icrs(reader, jd_tt, ipl, iflag)
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    elif body.coord_type == COORD_ICRS_BARY_SYSTEM:
        # Pipeline A': ICRS system barycenter — COB correction applied at runtime
        if iflag & SEFLG_SPEED:
            lon, lat, dist, dlon, dlat, ddist = _pipeline_icrs(
                reader,
                jd_tt,
                ipl,
                iflag,
                want_velocity=True,
                is_system_bary=True,
            )
        else:
            lon, lat, dist = _pipeline_icrs(
                reader,
                jd_tt,
                ipl,
                iflag,
                is_system_bary=True,
            )
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    elif body.coord_type == COORD_ECLIPTIC:
        # Pipeline B: ecliptic direct
        lon, lat, dist, dlon, dlat, ddist = _pipeline_ecliptic(
            reader, jd_tt, ipl, iflag
        )
        if not (iflag & SEFLG_SPEED):
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    elif body.coord_type == COORD_HELIO_ECL:
        # Pipeline C: heliocentric
        lon, lat, dist, dlon, dlat, ddist = _pipeline_helio(reader, jd_tt, ipl, iflag)
        if not (iflag & SEFLG_SPEED):
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    else:
        raise ValueError(f"Unknown coord_type {body.coord_type}")

    # Sidereal correction
    if (iflag & SEFLG_SIDEREAL) and not (iflag & SEFLG_EQUATORIAL):
        try:
            mean_aya = _calc_ayanamsa_from_leb(
                reader,
                jd_tt,
                sid_mode=sid_mode,
                sid_t0=sid_t0,
                sid_ayan_t0=sid_ayan_t0,
            )
            # True ayanamsa: add nutation in longitude (Δψ) to mean ayanamsa.
            # Ecliptic longitude from pipelines includes nutation, so we must
            # subtract the true ayanamsa (mean + Δψ) for correct cancellation.
            try:
                _, dpsi_sid, _, _ = _get_skyfield_frame_data(jd_tt)
                nutation_deg = math.degrees(dpsi_sid)
                aya = mean_aya + nutation_deg
            except (ValueError, Exception):
                aya = mean_aya
            lon = (lon - aya) % 360.0
            # Sidereal speed correction: subtract precession rate from dlon
            # _PREC_COEFFS are arcsec/century: dP/dT = c0 + 2*c1*T + ...
            # Convert: deg/day = (arcsec/century) / 3600 / 36525
            T = (jd_tt - J2000) / 36525.0
            prec_rate_arcsec_cy = _PREC_COEFFS[0] + 2 * _PREC_COEFFS[1] * T
            prec_rate_deg_day = prec_rate_arcsec_cy / (3600.0 * 36525.0)
            dlon -= prec_rate_deg_day
        except KeyError:
            # Star-based sidereal mode, fall back
            raise

    return (lon, lat, dist, dlon, dlat, ddist), iflag
