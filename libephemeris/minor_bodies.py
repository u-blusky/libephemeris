"""
Minor body calculations for asteroids and Trans-Neptunian Objects (TNOs).

This module computes positions for:
- Main belt asteroids: Ceres, Pallas, Juno, Vesta
- Centaurs: Chiron, Pholus
- Trans-Neptunian Objects (TNOs): Eris, Sedna, Haumea, Makemake, Orcus, Quaoar, Ixion

Method: Keplerian orbital elements with 2-body dynamics (Sun-body only).

IMPORTANT PRECISION LIMITATIONS:
- Uses simplified Keplerian orbits (no perturbations from planets)
- Accuracies: 1-5 arcminutes (1-5') typical for asteroids
- TNOs: Lower precision due to longer periods and perturbations
- For research-grade precision, use full numerical integration (Swiss Ephemeris, JPL Horizons)

This is intentionally simpler than Swiss Ephemeris's full dynamical model for:
- Faster computation
- No external data files required
- Acceptable precision for most astrological applications

Orbital elements source: JPL Small-Body Database (epoch 2025.0, JD 2461000.5)
Algorithm: Standard Keplerian orbital mechanics (Curtis, Vallado)
"""

import math
from dataclasses import dataclass
from typing import Tuple
from .constants import (
    SE_CHIRON,
    SE_PHOLUS,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
    SE_ERIS,
    SE_SEDNA,
    SE_HAUMEA,
    SE_MAKEMAKE,
    SE_IXION,
    SE_ORCUS,
    SE_QUAOAR,
)


@dataclass
class OrbitalElements:
    """
    Classical Keplerian orbital elements for a minor body.

    Attributes:
        name: Body name
        epoch: Reference epoch (Julian Day in TT)
        a: Semi-major axis in AU
        e: Eccentricity (0-1, dimensionless)
        i: Inclination to ecliptic in degrees
        omega: Argument of perihelion (ω) in degrees
        Omega: Longitude of ascending node (Ω) in degrees
        M0: Mean anomaly at epoch in degrees
        n: Mean motion in degrees/day

    Note:
        These are osculating elements at the given epoch.
        They drift over time due to planetary perturbations.
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float


# =============================================================================
# ORBITAL ELEMENTS DATABASE (Epoch 2025.0 = JD 2461000.5)
# =============================================================================
# Source: NASA JPL Small-Body Database (sbdb.api), retrieved 2025
# Epoch: JD 2461000.5 TDB (2025-Sep-19)
#
# Note: Elements are osculating at the given epoch.
# Accuracy degrades ~10-50 arcsec/year due to secular perturbations.
# For dates >10 years from epoch, consider updating elements.

MINOR_BODY_ELEMENTS = {
    SE_CHIRON: OrbitalElements(
        name="Chiron",
        epoch=2461000.5,
        a=13.6922,  # AU - between Saturn and Uranus
        e=0.378979,
        i=6.926,
        omega=339.254,  # argument of perihelion
        Omega=209.298,  # longitude of ascending node
        M0=212.840,
        n=0.01945,  # ~51 year period
    ),
    SE_PHOLUS: OrbitalElements(
        name="Pholus",
        epoch=2461000.5,
        a=20.2834,
        e=0.574745,  # Highly eccentric
        i=24.757,
        omega=354.730,
        Omega=119.290,
        M0=134.471,
        n=0.01079,  # ~91 year period
    ),
    SE_CERES: OrbitalElements(
        name="Ceres",
        epoch=2461000.5,
        a=2.7656,
        e=0.079576,
        i=10.588,
        omega=73.300,
        Omega=80.250,
        M0=231.540,
        n=0.21430,  # ~4.6 year period
    ),
    SE_PALLAS: OrbitalElements(
        name="Pallas",
        epoch=2461000.5,
        a=2.7699,
        e=0.230643,
        i=34.928,  # High inclination
        omega=310.933,
        Omega=172.889,
        M0=211.530,
        n=0.21380,  # ~4.6 year period
    ),
    SE_JUNO: OrbitalElements(
        name="Juno",
        epoch=2461000.5,
        a=2.6709,
        e=0.255826,
        i=12.986,
        omega=247.884,
        Omega=169.820,
        M0=217.591,
        n=0.22580,  # ~4.4 year period
    ),
    SE_VESTA: OrbitalElements(
        name="Vesta",
        epoch=2461000.5,
        a=2.3615,
        e=0.090168,
        i=7.144,
        omega=151.537,
        Omega=103.702,
        M0=26.810,
        n=0.27159,  # ~3.6 year period
    ),
    SE_ERIS: OrbitalElements(
        name="Eris",
        epoch=2461000.5,
        a=67.9964,  # Highly distant
        e=0.436965,
        i=43.869,  # Extreme inclination
        omega=150.732,
        Omega=36.027,
        M0=211.449,
        n=0.001758,  # ~561 year period
    ),
    SE_SEDNA: OrbitalElements(
        name="Sedna",
        epoch=2461000.5,
        a=549.541,  # Extreme distance (detached object)
        e=0.861297,  # Very eccentric
        i=11.926,
        omega=311.010,
        Omega=144.479,
        M0=358.607,
        n=0.0000765,  # ~12,880 year period
    ),
    SE_HAUMEA: OrbitalElements(
        name="Haumea",
        epoch=2461000.5,
        a=43.0055,
        e=0.195775,
        i=28.208,
        omega=240.888,
        Omega=121.797,
        M0=222.328,
        n=0.003495,  # ~282 year period
    ),
    SE_MAKEMAKE: OrbitalElements(
        name="Makemake",
        epoch=2461000.5,
        a=45.5107,
        e=0.160425,
        i=29.032,
        omega=297.075,
        Omega=79.269,
        M0=169.320,
        n=0.003210,  # ~307 year period
    ),
    SE_IXION: OrbitalElements(
        name="Ixion",
        epoch=2461000.5,
        a=39.3505,  # Plutino (2:3 resonance with Neptune)
        e=0.244233,
        i=19.670,
        omega=300.659,
        Omega=71.093,
        M0=294.200,
        n=0.003993,  # ~247 year period
    ),
    SE_ORCUS: OrbitalElements(
        name="Orcus",
        epoch=2461000.5,
        a=39.3358,  # Plutino (anti-Pluto phase)
        e=0.221730,
        i=20.556,
        omega=73.722,
        Omega=268.386,
        M0=188.111,
        n=0.003995,  # ~247 year period
    ),
    SE_QUAOAR: OrbitalElements(
        name="Quaoar",
        epoch=2461000.5,
        a=43.1477,
        e=0.035839,  # Nearly circular
        i=7.991,
        omega=163.923,
        Omega=188.963,
        M0=291.482,
        n=0.003478,  # ~284 year period
    ),
}


def solve_kepler_equation(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Solve Kepler's equation M = E - e·sin(E) for eccentric anomaly E.

    Uses Newton-Raphson iteration for robust convergence.

    Args:
        M: Mean anomaly in radians
        e: Eccentricity (0 ≤ e < 1)
        tol: Convergence tolerance (default 1e-8 ~ 0.002 arcsec)

    Returns:
        float: Eccentric anomaly E in radians

    Algorithm:
        Newton-Raphson: E_{n+1} = E_n - f(E_n)/f'(E_n)
        where f(E) = E - e·sin(E) - M
        and f'(E) = 1 - e·cos(E)

    Note:
        Converges in ~3-6 iterations for typical eccentricities (e < 0.8).
        Initial guess: M for e < 0.8, π for highly eccentric orbits.

    References:
        Curtis "Orbital Mechanics for Engineering Students" §3.1
        Vallado "Fundamentals of Astrodynamics" Algorithm 2
    """
    # FIXME: Precision - For parabolic/hyperbolic orbits (e ≥ 1), use different equation
    # This implementation assumes elliptical orbits (0 ≤ e < 1)
    E = M if e < 0.8 else math.pi

    for _ in range(30):
        f = E - e * math.sin(E) - M
        fp = 1 - e * math.cos(E)
        E_new = E - f / fp

        if abs(E_new - E) < tol:
            return E_new
        E = E_new

    return E


def calc_minor_body_position(
    elements: OrbitalElements, jd_tt: float
) -> Tuple[float, float, float]:
    """
    Calculate heliocentric position using Keplerian orbital elements.

    Propagates orbit from epoch to target time using mean motion.

    Args:
        elements: Orbital elements at epoch
        jd_tt: Target Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple[float, float, float]: (x, y, z) heliocentric position in AU
            Coordinates in ecliptic J2000.0 frame

    Algorithm:
        1. Propagate mean anomaly: M(t) = M0 + n·Δt
        2. Solve Kepler's equation for eccentric anomaly E
        3. Calculate true anomaly ν from E
        4. Compute position in orbital plane
        5. Rotate to ecliptic frame using Ω, i, ω

    FIXME: Precision - 2-body Keplerian propagation ignores:
        - Planetary perturbations (Jupiter, Saturn especially)
        - Non-gravitational forces (radiation pressure, Yarkovsky)
        - Relativistic effects (minor for asteroids)
        Typical errors: 1-5 arcminutes for asteroids, worse for TNOs

    References:
        Curtis §3 (orbital elements)
        Vallado §2.3 (coordinate transformations)
    """
    dt = jd_tt - elements.epoch

    # Propagate mean anomaly
    M = math.radians((elements.M0 + elements.n * dt) % 360.0)

    # Solve Kepler's equation
    E = solve_kepler_equation(M, elements.e)

    # True anomaly from eccentric anomaly
    nu = 2.0 * math.atan2(
        math.sqrt(1 + elements.e) * math.sin(E / 2),
        math.sqrt(1 - elements.e) * math.cos(E / 2),
    )

    # Heliocentric distance
    r = elements.a * (1 - elements.e * math.cos(E))

    # Position in orbital plane (perifocal frame)
    x_orb = r * math.cos(nu)
    y_orb = r * math.sin(nu)

    # Convert Euler angles to radians
    omega_rad = math.radians(elements.omega)
    Omega_rad = math.radians(elements.Omega)
    i_rad = math.radians(elements.i)

    # Precompute trig functions
    cos_omega = math.cos(omega_rad)
    sin_omega = math.sin(omega_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)
    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)

    # Rotation matrix from perifocal to ecliptic (standard transformation)
    P11 = cos_omega * cos_Omega - sin_omega * sin_Omega * cos_i
    P12 = -sin_omega * cos_Omega - cos_omega * sin_Omega * cos_i
    P21 = cos_omega * sin_Omega + sin_omega * cos_Omega * cos_i
    P22 = -sin_omega * sin_Omega + cos_omega * cos_Omega * cos_i
    P31 = sin_omega * sin_i
    P32 = cos_omega * sin_i

    # Transform to heliocentric ecliptic J2000 coordinates
    x = P11 * x_orb + P12 * y_orb
    y = P21 * x_orb + P22 * y_orb
    z = P31 * x_orb + P32 * y_orb

    return x, y, z


def calc_minor_body_heliocentric(
    body_id: int, jd_tt: float
) -> Tuple[float, float, float]:
    """
    Calculate heliocentric ecliptic coordinates for a minor body.

    Args:
        body_id: Minor body identifier (SE_CHIRON, SE_ERIS, etc.)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (-90 to +90)
            - distance: Heliocentric distance in AU

    Raises:
        ValueError: If body_id is not in the database

    Note:
        Returns HELIOCENTRIC coordinates. For geocentric, caller must
        subtract Earth's heliocentric position (see planets.py).

    Precision:
        Asteroids (Ceres, etc.): ~1-3 arcminutes typical
        TNOs (Eris, etc.): ~3-10 arcminutes typical
        Errors increase with time from epoch (2023.0)
    """
    if body_id not in MINOR_BODY_ELEMENTS:
        raise ValueError(f"Unknown minor body ID: {body_id}")

    elements = MINOR_BODY_ELEMENTS[body_id]
    x, y, z = calc_minor_body_position(elements, jd_tt)

    # Convert Cartesian to spherical coordinates
    r = math.sqrt(x**2 + y**2 + z**2)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(z / r))

    return lon, lat, r
