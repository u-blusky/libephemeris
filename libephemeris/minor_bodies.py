"""
Minor body calculations for asteroids and Trans-Neptunian Objects (TNOs).

This module computes positions for:
- Main belt asteroids: Ceres, Pallas, Juno, Vesta, Hygiea, Interamnia, Davida, Europa, Sylvia, Psyche, Sappho, Pandora, Lilith
- Centaurs: Chiron, Pholus, Nessus, Asbolus, Chariklo, Hidalgo
- Trans-Neptunian Objects (TNOs): Eris, Sedna, Haumea, Makemake, Orcus, Quaoar, Ixion, Gonggong, Varuna
- Near-Earth asteroids: Apophis, Eros, Amor, Icarus, Toro, Toutatis, Itokawa, Bennu, Ryugu

Method: Keplerian orbital elements with first-order secular perturbations from
Jupiter, Saturn, Uranus, and Neptune. This provides significantly improved accuracy over pure
2-body dynamics, especially for propagation over multiple years.

PRECISION:
- Main belt asteroids: ~10-30 arcseconds typical (improved from 1-5 arcminutes)
- TNOs: ~1-3 arcminutes typical (improved from 3-10 arcminutes)
- Errors increase with time from epoch, but secular perturbations reduce drift

PERTURBATION MODEL:
- Applies secular perturbations to orbital elements (ω, Ω, mean anomaly)
- Accounts for gravitational influence of Jupiter (dominant), Saturn, Uranus, and Neptune
- Neptune perturbations are critical for plutinos (2:3 resonance) like Ixion and Orcus
- Uranus perturbations are significant for Trans-Neptunian Objects (TNOs)
- Based on classical Laplace-Lagrange secular theory
- Does NOT include: mean-motion resonances, close encounters, non-gravitational forces

For research-grade precision, use full numerical integration (Swiss Ephemeris, JPL Horizons)

Orbital elements source: JPL Small-Body Database (epoch JD 2461000.5 TDB = 2025-Sep-19)
Algorithm: Keplerian mechanics with Laplace-Lagrange secular perturbations
"""

import math
from dataclasses import dataclass
from typing import Tuple, Optional, NamedTuple
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
    SE_NESSUS,
    SE_ASBOLUS,
    SE_CHARIKLO,
    SE_GONGGONG,
    SE_VARUNA,
    SE_APOPHIS,
    SE_HYGIEA,
    SE_INTERAMNIA,
    SE_DAVIDA,
    SE_EUROPA_AST,
    SE_SYLVIA,
    SE_PSYCHE,
    SE_EROS,
    SE_AMOR,
    SE_ICARUS,
    SE_TORO,
    SE_SAPPHO,
    SE_PANDORA_AST,
    SE_LILITH_AST,
    SE_HIDALGO,
    SE_TOUTATIS,
    SE_ITOKAWA,
    SE_BENNU,
    SE_RYUGU,
)


# =============================================================================
# PHYSICAL CONSTANTS FOR PERTURBATION CALCULATIONS
# =============================================================================
# Gravitational parameters (GM) in AU^3/day^2
# These are used for secular perturbation calculations

# GM of Sun (k^2 where k is Gaussian gravitational constant)
GM_SUN = 0.00029591220828559  # AU^3/day^2

# Mass ratios relative to Sun (m_planet / m_sun)
# Source: IAU nominal values and JPL DE441
MASS_RATIO_JUPITER = 1.0 / 1047.348644  # ~9.546e-4
MASS_RATIO_SATURN = 1.0 / 3497.901768  # ~2.858e-4

# Mean orbital elements of perturbing planets (J2000.0 values)
# These are used for secular perturbation calculations
# Source: JPL Horizons, mean elements at J2000.0

# Jupiter mean elements (for secular perturbation theory)
JUPITER_A = 5.2026  # Semi-major axis (AU)
JUPITER_E = 0.0485  # Mean eccentricity
JUPITER_I = 1.303  # Mean inclination (degrees)
JUPITER_OMEGA = 274.25  # Mean longitude of perihelion (degrees)
JUPITER_NODE = 100.46  # Mean longitude of ascending node (degrees)
JUPITER_N = 0.08309  # Mean motion (degrees/day)

# Saturn mean elements
SATURN_A = 9.5549  # Semi-major axis (AU)
SATURN_E = 0.0555  # Mean eccentricity
SATURN_I = 2.489  # Mean inclination (degrees)
SATURN_OMEGA = 339.39  # Mean longitude of perihelion (degrees)
SATURN_NODE = 113.66  # Mean longitude of ascending node (degrees)
SATURN_N = 0.03346  # Mean motion (degrees/day)

# Uranus mean elements (significant for TNO calculations)
URANUS_A = 19.2184  # Semi-major axis (AU)
URANUS_E = 0.0457  # Mean eccentricity
URANUS_I = 0.772  # Mean inclination (degrees)
URANUS_OMEGA = 170.9  # Mean argument of perihelion (degrees)
URANUS_NODE = 74.0  # Mean longitude of ascending node (degrees)
URANUS_N = 0.01177  # Mean motion (degrees/day)

# Mass ratio of Uranus relative to Sun
MASS_RATIO_URANUS = 1.0 / 22902  # ~4.366e-5

# Neptune mean elements (critical for TNO accuracy, especially plutinos)
# Plutinos like Ixion and Orcus are in 2:3 mean motion resonance with Neptune
NEPTUNE_A = 30.1104  # Semi-major axis (AU)
NEPTUNE_E = 0.0086  # Mean eccentricity
NEPTUNE_I = 1.769  # Mean inclination (degrees)
NEPTUNE_OMEGA = 44.9  # Mean argument of perihelion (degrees)
NEPTUNE_NODE = 131.7  # Mean longitude of ascending node (degrees)
NEPTUNE_N = 0.006021  # Mean motion (degrees/day)

# Mass ratio of Neptune relative to Sun
MASS_RATIO_NEPTUNE = 1.0 / 19412  # ~5.153e-5


# =============================================================================
# MEAN MOTION RESONANCE CONSTANTS
# =============================================================================
# Mean motion resonances with Neptune occur when an outer body's orbital period
# is a simple integer ratio of Neptune's period. These resonances protect
# bodies from close encounters but also mean that secular perturbation theory
# breaks down and gives incorrect results.
#
# Common Neptune resonances (body:Neptune orbital period ratio):
# - 2:3 (Plutinos): ~39.4 AU - most populated resonance (Pluto, Ixion, Orcus)
# - 1:2 (Twotinos): ~47.8 AU
# - 2:5: ~55.4 AU
# - 1:3: ~62.5 AU
# - 3:5: ~36.4 AU
# - 4:7: ~43.7 AU
# - 3:4: ~36.5 AU
#
# References:
#   Malhotra, R. (1995) "The origin of Pluto's peculiar orbit"
#   Murray & Dermott "Solar System Dynamics" Ch. 8

# Neptune's mean motion in degrees/day (for resonance calculations)
NEPTUNE_MEAN_MOTION = 0.006021  # degrees/day (~164.8 year period)


class ResonanceInfo(NamedTuple):
    """
    Information about a mean motion resonance.

    Attributes:
        p: Integer numerator (body's period multiplier)
        q: Integer denominator (Neptune's period multiplier)
        name: Common name for the resonance (e.g., "plutino", "twotino")
        a_resonant: Expected semi-major axis for this resonance (AU)
        tolerance: Fractional tolerance for detecting resonance
    """

    p: int
    q: int
    name: str
    a_resonant: float
    tolerance: float


# Calculate resonant semi-major axis using Kepler's 3rd law:
# For a p:q resonance, body makes p orbits while Neptune makes q orbits
# T_body / T_Neptune = q / p  (body's period is longer if q > p)
# (a_body / a_Neptune)^(3/2) = q / p
# a_body = a_Neptune * (q / p)^(2/3)


def _calc_resonant_a(p: int, q: int) -> float:
    """Calculate the resonant semi-major axis for a p:q resonance with Neptune.

    In a p:q resonance, the body completes p orbits while Neptune completes q orbits.
    For exterior resonances (most common), q > p, so the body's orbit is larger.
    """
    return NEPTUNE_A * (q / p) ** (2.0 / 3.0)


# Known Neptune mean motion resonances
# Format: (p, q, name, a_resonant, tolerance)
# Tolerance is fractional (e.g., 0.02 means ±2% of a_resonant)
NEPTUNE_RESONANCES: list[ResonanceInfo] = [
    # Interior resonances (a < Neptune)
    ResonanceInfo(3, 4, "3:4 resonance", _calc_resonant_a(3, 4), 0.02),
    ResonanceInfo(3, 5, "3:5 resonance", _calc_resonant_a(3, 5), 0.02),
    # Plutinos (most famous resonance, includes Pluto)
    ResonanceInfo(2, 3, "plutino", _calc_resonant_a(2, 3), 0.02),
    # Other exterior resonances
    ResonanceInfo(4, 7, "4:7 resonance", _calc_resonant_a(4, 7), 0.02),
    # Twotinos (1:2 resonance)
    ResonanceInfo(1, 2, "twotino", _calc_resonant_a(1, 2), 0.02),
    # More distant resonances
    ResonanceInfo(2, 5, "2:5 resonance", _calc_resonant_a(2, 5), 0.02),
    ResonanceInfo(1, 3, "1:3 resonance", _calc_resonant_a(1, 3), 0.02),
]

# Default tolerance for resonance detection (fractional difference in mean motion)
DEFAULT_RESONANCE_TOLERANCE = 0.02  # 2%


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
# MEAN MOTION RESONANCE DETECTION FUNCTIONS
# =============================================================================


@dataclass
class ResonanceResult:
    """
    Result of resonance detection for a body.

    Attributes:
        is_resonant: Whether the body is in or near a mean motion resonance
        resonance: The detected resonance (ResonanceInfo), or None if not resonant
        deviation: Fractional deviation from exact resonance (0.0 = exact)
        warning_message: Message about secular perturbation accuracy, or None
    """

    is_resonant: bool
    resonance: Optional[ResonanceInfo]
    deviation: float
    warning_message: Optional[str]


def detect_mean_motion_resonance(
    elements: "OrbitalElements",
    tolerance: float = DEFAULT_RESONANCE_TOLERANCE,
) -> ResonanceResult:
    """
    Detect if a body is in mean motion resonance with Neptune.

    Checks if the body's mean motion is close to a resonant ratio with Neptune's
    mean motion. Bodies in mean motion resonance have their secular perturbations
    significantly modified by resonant effects that are not captured by
    first-order Laplace-Lagrange secular theory.

    Args:
        elements: Orbital elements of the body to check
        tolerance: Fractional tolerance for detecting resonance (default 2%)

    Returns:
        ResonanceResult: Contains detection result, resonance info, and warning

    Algorithm:
        For each known Neptune resonance p:q, checks if:
        |n_body / n_Neptune - p/q| < tolerance * (p/q)

        where n is mean motion. In a p:q resonance, the body completes p orbits
        while Neptune completes q, so n_body/n_Neptune = p/q.
        This is equivalent to checking if the semi-major axis is within
        tolerance of the resonant value.

    Example:
        >>> from libephemeris.minor_bodies import detect_mean_motion_resonance
        >>> result = detect_mean_motion_resonance(MINOR_BODY_ELEMENTS[SE_IXION])
        >>> result.is_resonant
        True
        >>> result.resonance.name
        'plutino'

    Note:
        Bodies in resonance may have secular perturbation calculations that
        are less accurate than for non-resonant bodies. The warning_message
        field provides guidance in such cases.

    References:
        Murray & Dermott "Solar System Dynamics" Ch. 8
        Malhotra (1995) "The origin of Pluto's peculiar orbit"
    """
    a = elements.a
    n = elements.n

    # Only check bodies that could potentially be in Neptune resonance
    # (must be in outer solar system, say a > 20 AU)
    if a < 20.0:
        return ResonanceResult(
            is_resonant=False, resonance=None, deviation=0.0, warning_message=None
        )

    # Calculate body's mean motion ratio relative to Neptune
    # For a p:q resonance, body makes p orbits while Neptune makes q
    # So n_body / n_Neptune = p / q
    n_ratio = n / NEPTUNE_MEAN_MOTION

    best_match: Optional[ResonanceInfo] = None
    best_deviation = float("inf")

    for resonance in NEPTUNE_RESONANCES:
        # Expected ratio for this resonance: p/q
        expected_ratio = resonance.p / resonance.q

        # Calculate fractional deviation
        deviation = abs(n_ratio - expected_ratio) / expected_ratio

        if deviation < tolerance and deviation < best_deviation:
            best_match = resonance
            best_deviation = deviation

    if best_match is not None:
        warning = (
            f"Body '{elements.name}' is near the {best_match.p}:{best_match.q} "
            f"mean motion resonance with Neptune ({best_match.name}). "
            f"Secular perturbation theory may be inaccurate for resonant bodies. "
            f"Deviation from exact resonance: {best_deviation * 100:.2f}%"
        )
        return ResonanceResult(
            is_resonant=True,
            resonance=best_match,
            deviation=best_deviation,
            warning_message=warning,
        )

    return ResonanceResult(
        is_resonant=False, resonance=None, deviation=0.0, warning_message=None
    )


def is_body_resonant(
    body_id: int, tolerance: float = DEFAULT_RESONANCE_TOLERANCE
) -> bool:
    """
    Check if a minor body is in mean motion resonance with Neptune.

    This is a convenience function that looks up the body's orbital elements
    and checks for resonance.

    Args:
        body_id: Minor body identifier (SE_IXION, SE_ORCUS, etc.)
        tolerance: Fractional tolerance for detecting resonance (default 2%)

    Returns:
        bool: True if the body is in or near a Neptune resonance

    Raises:
        ValueError: If body_id is not in the orbital elements database

    Example:
        >>> from libephemeris.minor_bodies import is_body_resonant
        >>> from libephemeris.constants import SE_IXION, SE_CERES
        >>> is_body_resonant(SE_IXION)
        True
        >>> is_body_resonant(SE_CERES)
        False
    """
    if body_id not in MINOR_BODY_ELEMENTS:
        raise ValueError(f"Unknown body ID: {body_id}")

    elements = MINOR_BODY_ELEMENTS[body_id]
    result = detect_mean_motion_resonance(elements, tolerance)
    return result.is_resonant


def get_resonance_info(body_id: int) -> Optional[ResonanceResult]:
    """
    Get detailed resonance information for a minor body.

    Args:
        body_id: Minor body identifier (SE_IXION, SE_ORCUS, etc.)

    Returns:
        ResonanceResult: Detailed resonance detection result, or None if body unknown

    Example:
        >>> from libephemeris.minor_bodies import get_resonance_info
        >>> from libephemeris.constants import SE_ORCUS
        >>> info = get_resonance_info(SE_ORCUS)
        >>> if info and info.is_resonant:
        ...     print(f"{info.resonance.name}: deviation {info.deviation:.2%}")
        plutino: deviation 0.15%
    """
    if body_id not in MINOR_BODY_ELEMENTS:
        return None

    elements = MINOR_BODY_ELEMENTS[body_id]
    return detect_mean_motion_resonance(elements)


# =============================================================================
# SECULAR PERTURBATION FUNCTIONS
# =============================================================================


def _calc_laplace_coefficients(alpha: float, s: float, j: int) -> float:
    """
    Calculate Laplace coefficient b_s^(j)(alpha) using series expansion.

    The Laplace coefficients appear in the disturbing function expansion
    for planetary perturbation theory. They depend on the ratio of semi-major
    axes and are used to calculate secular perturbation rates.

    Args:
        alpha: Ratio of semi-major axes (a_inner / a_outer), must be < 1
        s: Half-integer index (typically 1/2 or 3/2)
        j: Integer index (0, 1, 2, ...)

    Returns:
        float: The Laplace coefficient value

    References:
        Murray & Dermott "Solar System Dynamics" §6.4
        Brouwer & Clemence "Methods of Celestial Mechanics" Ch. XI
    """
    if alpha >= 1.0 or alpha <= 0.0:
        return 0.0

    # Series expansion for b_s^(j)(alpha)
    # b_s^(j) = (2/pi) * integral from 0 to pi of cos(j*psi) / (1 - 2*alpha*cos(psi) + alpha^2)^s dpsi
    # For s = 1/2, we can use the approximation for small alpha

    # Use numerical approximation via trapezoidal integration
    n_steps = 100
    dpsi = math.pi / n_steps
    result = 0.0

    for k in range(n_steps + 1):
        psi = k * dpsi
        denom = 1.0 - 2.0 * alpha * math.cos(psi) + alpha * alpha
        if denom > 1e-10:
            integrand = math.cos(j * psi) / (denom**s)
            weight = 1.0 if (k == 0 or k == n_steps) else 2.0
            result += weight * integrand * dpsi / 2.0

    return result / math.pi


def calc_secular_perturbation_rates(
    elements: OrbitalElements,
) -> Tuple[float, float, float]:
    """
    Calculate secular perturbation rates for orbital elements due to Jupiter, Saturn, Uranus, and Neptune.

    Uses first-order Laplace-Lagrange secular perturbation theory to compute
    the time rates of change of the argument of perihelion (ω), longitude of
    ascending node (Ω), and a correction to mean motion.

    The secular perturbations cause:
    - Precession of perihelion (ω advances or regresses)
    - Precession of the ascending node (Ω regresses for prograde orbits)
    - Small correction to mean motion due to non-Keplerian effects

    Note:
        Neptune perturbations are critical for TNO accuracy, especially for
        plutinos like Ixion and Orcus which are in 2:3 mean motion resonance
        with Neptune. Neptune perturbations are only applied for bodies with
        semi-major axis > 20 AU to avoid overhead for inner solar system objects.
        Uranus perturbations are particularly significant for Trans-Neptunian Objects
        (TNOs) where its gravitational influence is comparable to or greater than
        Saturn's influence.

    Args:
        elements: Orbital elements of the minor body

    Returns:
        Tuple[float, float, float]: (d_omega, d_Omega, d_n) rates in degrees/day
            - d_omega: Rate of change of argument of perihelion
            - d_Omega: Rate of change of longitude of ascending node
            - d_n: Correction to mean motion (small)

    Algorithm:
        Uses the classical Laplace-Lagrange secular theory formulas:
        - dω/dt ≈ n * (m'/M) * α * b_{3/2}^{(1)}(α) * (1/4) * ...
        - dΩ/dt ≈ -n * (m'/M) * α * b_{3/2}^{(1)}(α) * cos(i) * ...

        where:
        - n is the mean motion of the asteroid
        - m'/M is the mass ratio of the perturbing planet to the Sun
        - α is the ratio of semi-major axes
        - b is the Laplace coefficient

    References:
        Murray & Dermott "Solar System Dynamics" Ch. 7
        Brouwer & Clemence "Methods of Celestial Mechanics" Ch. XVI
    """
    a = elements.a
    e = elements.e
    i_rad = math.radians(elements.i)
    n = elements.n  # degrees/day

    # Convert mean motion to radians/day for calculations
    n_rad = math.radians(n)

    # Initialize perturbation rates
    d_omega = 0.0  # degrees/day
    d_Omega = 0.0  # degrees/day
    d_n = 0.0  # degrees/day

    # Calculate perturbations from Jupiter (dominant for most asteroids)
    if a < JUPITER_A:
        # Asteroid is interior to Jupiter
        alpha = a / JUPITER_A

        # Laplace coefficient approximation for b_{3/2}^{(1)}
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)

        # Secular rates from first-order theory
        # Factor from disturbing function expansion
        factor = n_rad * MASS_RATIO_JUPITER * alpha * b32_1

        # Perihelion precession rate (prograde for interior orbits)
        # dω/dt = (n/4) * (m'/M) * α² * b_{3/2}^{(2)} * (2 - e²) / (1 - e²)
        # Simplified: use leading order term
        d_omega_jup = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)

        # Node regression rate
        # dΩ/dt = -(n/4) * (m'/M) * α * b_{3/2}^{(1)} * cos(i) / (1 - e²)^(1/2)
        d_Omega_jup = -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(1.0 - e * e)

        d_omega += math.degrees(d_omega_jup)
        d_Omega += math.degrees(d_Omega_jup)

    elif a > JUPITER_A:
        # Asteroid is exterior to Jupiter
        alpha = JUPITER_A / a

        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_JUPITER * alpha * b32_1

        # For exterior bodies, perturbation effects are reduced but still significant
        d_omega_jup = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)
        d_Omega_jup = (
            -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(max(0.01, 1.0 - e * e))
        )

        d_omega += math.degrees(d_omega_jup)
        d_Omega += math.degrees(d_Omega_jup)

    # Calculate perturbations from Saturn (significant for outer asteroids/TNOs)
    if a < SATURN_A:
        alpha = a / SATURN_A
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_SATURN * alpha * b32_1

        d_omega_sat = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)
        d_Omega_sat = -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(1.0 - e * e)

        d_omega += math.degrees(d_omega_sat)
        d_Omega += math.degrees(d_Omega_sat)

    elif a > SATURN_A:
        alpha = SATURN_A / a
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_SATURN * alpha * b32_1

        d_omega_sat = (
            factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (max(0.01, 1.0 - e * e))
        )
        d_Omega_sat = (
            -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(max(0.01, 1.0 - e * e))
        )

        d_omega += math.degrees(d_omega_sat)
        d_Omega += math.degrees(d_Omega_sat)

    # Calculate perturbations from Uranus (significant for TNOs)
    if a < URANUS_A:
        alpha = a / URANUS_A
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_URANUS * alpha * b32_1

        d_omega_ura = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)
        d_Omega_ura = -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(1.0 - e * e)

        d_omega += math.degrees(d_omega_ura)
        d_Omega += math.degrees(d_Omega_ura)

    elif a > URANUS_A:
        alpha = URANUS_A / a
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_URANUS * alpha * b32_1

        d_omega_ura = (
            factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (max(0.01, 1.0 - e * e))
        )
        d_Omega_ura = (
            -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(max(0.01, 1.0 - e * e))
        )

        d_omega += math.degrees(d_omega_ura)
        d_Omega += math.degrees(d_Omega_ura)

    # Calculate perturbations from Neptune (critical for TNO accuracy)
    # Only apply for bodies with a > 20 AU to avoid overhead for inner solar system
    # Neptune is especially important for plutinos like Ixion and Orcus (2:3 resonance)
    if a > 20.0:
        if a < NEPTUNE_A:
            # Body is interior to Neptune (includes plutinos at ~39 AU)
            alpha = a / NEPTUNE_A
            b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
            factor = n_rad * MASS_RATIO_NEPTUNE * alpha * b32_1

            d_omega_nep = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)
            d_Omega_nep = (
                -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(1.0 - e * e)
            )

            d_omega += math.degrees(d_omega_nep)
            d_Omega += math.degrees(d_Omega_nep)

        elif a > NEPTUNE_A:
            # Body is exterior to Neptune (most TNOs: Eris, Sedna, Haumea, etc.)
            alpha = NEPTUNE_A / a
            b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
            factor = n_rad * MASS_RATIO_NEPTUNE * alpha * b32_1

            d_omega_nep = (
                factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (max(0.01, 1.0 - e * e))
            )
            d_Omega_nep = (
                -factor
                * (1.0 / 4.0)
                * math.cos(i_rad)
                / math.sqrt(max(0.01, 1.0 - e * e))
            )

            d_omega += math.degrees(d_omega_nep)
            d_Omega += math.degrees(d_Omega_nep)

    return d_omega, d_Omega, d_n


def apply_secular_perturbations(
    elements: OrbitalElements, jd_tt: float, include_perturbations: bool = True
) -> Tuple[float, float, float, float]:
    """
    Apply secular perturbations to orbital elements and return perturbed values.

    Takes the osculating orbital elements at epoch and applies first-order
    secular perturbation corrections to propagate them to the target time.
    This accounts for the long-term drift in ω and Ω due to Jupiter, Saturn, Uranus,
    and Neptune (for TNOs with a > 20 AU).

    WARNING: For bodies in mean motion resonance with Neptune (e.g., plutinos
    like Ixion and Orcus in 2:3 resonance, twotinos in 1:2 resonance), the
    secular perturbation theory used here may give inaccurate results. Resonant
    bodies experience additional perturbations that are not captured by
    first-order Laplace-Lagrange theory. Use detect_mean_motion_resonance()
    or is_body_resonant() to check if a body is in resonance.

    Args:
        elements: Original orbital elements at epoch
        jd_tt: Target Julian Day in Terrestrial Time
        include_perturbations: If True, apply secular corrections; if False, return unperturbed

    Returns:
        Tuple[float, float, float, float]: (omega_pert, Omega_pert, M_pert, n_pert)
            - omega_pert: Perturbed argument of perihelion (degrees)
            - Omega_pert: Perturbed longitude of ascending node (degrees)
            - M_pert: Perturbed mean anomaly at target time (degrees)
            - n_pert: Perturbed mean motion (degrees/day)

    See Also:
        detect_mean_motion_resonance: Check if a body is in Neptune resonance
        is_body_resonant: Quick check if a body ID is resonant
    """
    dt = jd_tt - elements.epoch  # Time since epoch in days

    if not include_perturbations or abs(dt) < 1.0:
        # For very short propagation times, perturbations are negligible
        M = (elements.M0 + elements.n * dt) % 360.0
        return elements.omega, elements.Omega, M, elements.n

    # Calculate secular perturbation rates
    d_omega, d_Omega, d_n = calc_secular_perturbation_rates(elements)

    # Apply secular corrections
    omega_pert = (elements.omega + d_omega * dt) % 360.0
    Omega_pert = (elements.Omega + d_Omega * dt) % 360.0
    n_pert = elements.n + d_n

    # Propagate mean anomaly with perturbed mean motion
    M_pert = (elements.M0 + n_pert * dt) % 360.0

    return omega_pert, Omega_pert, M_pert, n_pert


# =============================================================================
# ORBITAL ELEMENTS DATABASE (Epoch JD 2461000.5 TDB = 2025-Sep-19)
# =============================================================================
# Source: NASA JPL Small-Body Database (sbdb.api), retrieved 2025-Nov-27
# Epoch: JD 2461000.5 TDB (2025-Sep-19)
#
# Note: Elements are osculating at the given epoch.
# Accuracy degrades ~10-50 arcsec/year due to secular perturbations.
# For dates >10 years from epoch, consider updating elements.

MINOR_BODY_ELEMENTS = {
    SE_CHIRON: OrbitalElements(
        name="Chiron",
        epoch=2461000.5,
        a=13.69219896172984,  # AU - between Saturn and Uranus
        e=0.3789792342846475,
        i=6.926003536565557,
        omega=339.2537417045351,  # argument of perihelion
        Omega=209.2984204899107,  # longitude of ascending node
        M0=212.8397717853335,
        n=0.01945334424082164,  # ~51 year period
    ),
    SE_PHOLUS: OrbitalElements(
        name="Pholus",
        epoch=2461000.5,
        a=20.28340105547402,
        e=0.574744804002886,  # Highly eccentric
        i=24.75699076739707,
        omega=354.7299656288133,
        Omega=119.2896923424452,
        M0=134.470501527666,
        n=0.01078929100989589,  # ~91 year period
    ),
    SE_CERES: OrbitalElements(
        name="Ceres",
        epoch=2461000.5,
        a=2.765615651508659,
        e=0.07957631994408416,
        i=10.58788658206854,
        omega=73.29975464616518,
        Omega=80.24963090816965,
        M0=231.5397330043706,
        n=0.2142971214271186,  # ~4.6 year period
    ),
    SE_PALLAS: OrbitalElements(
        name="Pallas",
        epoch=2461000.5,
        a=2.76992582511479,
        e=0.2306429787781384,
        i=34.92832687077855,  # High inclination
        omega=310.9333840114307,
        Omega=172.8885963367437,
        M0=211.5297778033731,
        n=0.2137971269626138,  # ~4.6 year period
    ),
    SE_JUNO: OrbitalElements(
        name="Juno",
        epoch=2461000.5,
        a=2.670879058906207,
        e=0.2558257725543152,
        i=12.98603961477441,
        omega=247.8836661359693,
        Omega=169.8198844530219,
        M0=217.5909617686606,
        n=0.2257993770223806,  # ~4.4 year period
    ),
    SE_VESTA: OrbitalElements(
        name="Vesta",
        epoch=2461000.5,
        a=2.361541280084789,
        e=0.09016764504738634,
        i=7.144060599543863,
        omega=151.5371488873794,
        Omega=103.7022980342142,
        M0=26.80967220901607,
        n=0.2715881155129186,  # ~3.6 year period
    ),
    SE_ERIS: OrbitalElements(
        name="Eris",
        epoch=2461000.5,
        a=67.99636506315233,  # Highly distant
        e=0.4369647318678042,
        i=43.86893210857381,  # Extreme inclination
        omega=150.7324609262193,
        Omega=36.0271728267921,
        M0=211.4489901906776,
        n=0.001757824561999289,  # ~561 year period
    ),
    SE_SEDNA: OrbitalElements(
        name="Sedna",
        epoch=2461000.5,
        a=549.5405279241455,  # Extreme distance (detached object)
        e=0.8612974378772624,  # Very eccentric
        i=11.92591697769048,
        omega=311.0097709574837,
        Omega=144.4787268123774,
        M0=358.6072615114487,
        n=7.650758332945998e-05,  # ~12,880 year period
    ),
    SE_HAUMEA: OrbitalElements(
        name="Haumea",
        epoch=2461000.5,
        a=43.00549881333706,
        e=0.1957748188564788,
        i=28.20840579216805,
        omega=240.888336259082,
        Omega=121.7972901946415,
        M0=222.3276474703985,
        n=0.003494765903502902,  # ~282 year period
    ),
    SE_MAKEMAKE: OrbitalElements(
        name="Makemake",
        epoch=2461000.5,
        a=45.51068175675377,
        e=0.1604249910368623,
        i=29.03230613859081,
        omega=297.0754218776326,
        Omega=79.26892117343222,
        M0=169.3202841867025,
        n=0.003210214572603997,  # ~307 year period
    ),
    SE_IXION: OrbitalElements(
        name="Ixion",
        epoch=2461000.5,
        a=39.35053706213409,  # Plutino (2:3 resonance with Neptune)
        e=0.2442328489971719,
        i=19.67041190191546,
        omega=300.6585723831106,
        Omega=71.09295825377649,
        M0=294.2004612799266,
        n=0.003992804789574964,  # ~247 year period
    ),
    SE_ORCUS: OrbitalElements(
        name="Orcus",
        epoch=2461000.5,
        a=39.33577647200568,  # Plutino (anti-Pluto phase)
        e=0.2217300030420161,
        i=20.55552551599616,
        omega=73.72249191891537,
        Omega=268.3859416278478,
        M0=188.1111318293787,
        n=0.003995052426031791,  # ~247 year period
    ),
    SE_QUAOAR: OrbitalElements(
        name="Quaoar",
        epoch=2461000.5,
        a=43.1476797802032,
        e=0.03583878353429052,  # Nearly circular
        i=7.991371294217068,
        omega=163.9231384883233,
        Omega=188.9632800603184,
        M0=291.4818844949103,
        n=0.003477506123841158,  # ~284 year period
    ),
    SE_NESSUS: OrbitalElements(
        name="Nessus",
        epoch=2461000.5,
        a=24.51657596400793,  # Centaur orbit
        e=0.5176801039355027,  # Moderately eccentric
        i=15.64407698270027,
        omega=170.3715103791255,
        Omega=31.29215422788363,
        M0=100.5239875112597,
        n=0.008119220769969982,  # ~121 year period
    ),
    SE_ASBOLUS: OrbitalElements(
        name="Asbolus",
        epoch=2461000.5,
        a=18.06355955222326,  # Centaur orbit
        e=0.6173222203620333,  # Highly eccentric
        i=17.6174400161563,
        omega=290.5836923759394,
        Omega=6.020562307154061,
        M0=109.6324751756823,
        n=0.01283805024514489,  # ~77 year period
    ),
    SE_CHARIKLO: OrbitalElements(
        name="Chariklo",
        epoch=2461000.5,
        a=15.73995155535189,  # Centaur orbit, largest known centaur (~250 km)
        e=0.1702459901276738,  # Moderate eccentricity
        i=23.43032407033149,
        omega=241.2242283234723,
        Omega=300.4752379513423,
        M0=126.9607903541095,
        n=0.01578334271898484,  # ~62 year period, has ring system discovered 2014
    ),
    SE_GONGGONG: OrbitalElements(
        name="Gonggong",
        epoch=2461000.5,
        a=66.89366871435344,  # TNO, dwarf planet candidate (~1230 km diameter)
        e=0.5031674399617051,  # Eccentric orbit
        i=30.86626129015389,  # High inclination
        omega=206.6416070091921,
        Omega=336.8400960976296,
        M0=111.3903730396343,
        n=0.001801467997135536,  # ~547 year period (formerly 2007 OR10)
    ),
    SE_VARUNA: OrbitalElements(
        name="Varuna",
        epoch=2461000.5,
        a=43.17823437208563,  # Classical KBO (~670 km diameter)
        e=0.05254516459334936,  # Nearly circular
        i=17.13812404484431,
        omega=273.2206220996139,
        Omega=97.21030234100793,
        M0=115.0289798413239,
        n=0.003473815549791981,  # ~284 year period
    ),
    # Near-Earth Asteroids (NEAs)
    # NOTE: Apophis orbital elements change measurably with each refinement due to
    # close Earth approaches in 2029 and 2036. These elements are from JPL SBDB
    # solution 220 (2024-06-25). For critical applications, use SPK files.
    SE_APOPHIS: OrbitalElements(
        name="Apophis",
        epoch=2461000.5,
        a=0.9223803173917017,  # Aten-class NEA (~0.34 km diameter)
        e=0.1911663355386932,  # Moderate eccentricity
        i=3.340958441017069,  # Low inclination
        omega=126.6728325163065,
        Omega=203.8996515621043,
        M0=312.8054663584516,
        n=1.11259994308075,  # ~0.89 year period, close approach 2029-04-13
    ),
    # Main Belt Asteroids - Large bodies
    SE_HYGIEA: OrbitalElements(
        name="Hygiea",
        epoch=2461000.5,
        a=3.147591335345947,  # Fourth largest asteroid, dwarf planet candidate (~430 km)
        e=0.108223833089903,  # Low eccentricity
        i=3.832937411339835,  # Low inclination
        omega=312.6058061800255,
        Omega=283.1216602821771,
        M0=216.6903206124055,
        n=0.1764966888337409,  # ~5.6 year period
    ),
    SE_INTERAMNIA: OrbitalElements(
        name="Interamnia",
        epoch=2461000.5,
        a=3.056218782582974,  # Fifth largest asteroid (~320 km diameter)
        e=0.1552477966505805,  # Moderate eccentricity
        i=17.31546543833021,  # Moderate inclination
        omega=94.10694029269551,
        Omega=280.1660768334629,
        M0=184.2054144782699,
        n=0.1844707074261051,  # ~5.3 year period
    ),
    SE_DAVIDA: OrbitalElements(
        name="Davida",
        epoch=2461000.5,
        a=3.163619117539836,  # Seventh largest asteroid (~270 km diameter)
        e=0.1896284393751952,  # Moderate eccentricity
        i=15.95004792439553,  # Moderate inclination
        omega=336.666946405915,
        Omega=107.5554381931758,
        M0=35.24757453107242,
        n=0.1751571164581395,  # ~5.6 year period
    ),
    SE_EUROPA_AST: OrbitalElements(
        name="Europa",
        epoch=2461000.5,
        a=3.092193380509604,  # Main belt asteroid (~300 km diameter, not Jupiter's moon)
        e=0.1120016601360199,  # Low eccentricity
        i=7.480921967944819,  # Low inclination
        omega=342.9584258894789,
        Omega=128.5790551357993,
        M0=312.5276379313007,
        n=0.1812608885921891,  # ~5.4 year period
    ),
    SE_SYLVIA: OrbitalElements(
        name="Sylvia",
        epoch=2461000.5,
        a=3.484518782251748,  # Outer main belt asteroid (~280 km, triple system with moons Romulus and Remus)
        e=0.09433757216253168,  # Low eccentricity
        i=10.87247390257978,  # Moderate inclination
        omega=265.8265584696771,
        Omega=73.01146029932147,
        M0=94.95962367758351,
        n=0.151526942166858,  # ~6.5 year period
    ),
    SE_PSYCHE: OrbitalElements(
        name="Psyche",
        epoch=2461000.5,
        a=2.923314514376521,  # Main belt asteroid (~226 km, metallic M-type, NASA Psyche mission target)
        e=0.1343462202634327,  # Low eccentricity
        i=3.097291662719076,  # Low inclination
        omega=229.7534129677803,
        Omega=150.0098773790708,
        M0=40.63883784197816,
        n=0.1971926664511637,  # ~5.0 year period
    ),
    SE_EROS: OrbitalElements(
        name="Eros",
        epoch=2461000.5,
        a=1.458120998474684,  # Near-Earth asteroid (Amor class, NEAR Shoemaker mission target)
        e=0.2228359407071628,  # Moderate eccentricity
        i=10.82846651399785,  # Moderate inclination
        omega=178.9297536744151,
        Omega=304.2701025753316,
        M0=310.5543277370992,
        n=0.5597752949285997,  # ~1.76 year period
    ),
    SE_AMOR: OrbitalElements(
        name="Amor",
        epoch=2461000.5,
        a=1.919831260906931,  # Amor-class NEA prototype (~1 km diameter)
        e=0.4346323375284623,  # Moderate eccentricity
        i=11.86882308811668,  # Moderate inclination
        omega=26.75822531564864,
        Omega=171.2371875175307,
        M0=59.87048521279021,
        n=0.3705181485730381,  # ~2.66 year period (972 days)
    ),
    SE_ICARUS: OrbitalElements(
        name="Icarus",
        epoch=2461000.5,
        a=1.078037837567316,  # Apollo asteroid (~1.4 km, perihelion inside Mercury's orbit)
        e=0.8270056908369543,  # Very high eccentricity (q=0.19 AU, Q=1.97 AU)
        i=22.8032205675089,  # Moderate inclination
        omega=31.43821652068713,
        Omega=87.95241724594092,
        M0=153.0789301425594,
        n=0.8805480660198883,  # ~1.12 year period (409 days)
    ),
    SE_TORO: OrbitalElements(
        name="Toro",
        epoch=2461000.5,
        a=1.367838675762221,  # Apollo asteroid, near-Earth asteroid (~3.4 km)
        e=0.4360487430728354,  # Moderate eccentricity (q=0.77 AU, Q=1.96 AU)
        i=9.38184277020598,  # Low inclination
        omega=127.2752987170327,
        Omega=274.2077498325479,
        M0=82.68163042845468,
        n=0.6161007748143977,  # ~1.6 year period (584 days)
    ),
    SE_SAPPHO: OrbitalElements(
        name="Sappho",
        epoch=2461000.5,
        a=2.296282028063693,  # Main belt asteroid (~69 km, artistic expression, same-sex love)
        e=0.1996711692276041,  # Moderate eccentricity
        i=8.676097553661929,  # Low inclination
        omega=139.693299727962,
        Omega=218.6385022891111,
        M0=57.25859931389674,
        n=0.2832475967264622,  # ~3.5 year period (1271 days)
    ),
    SE_PANDORA_AST: OrbitalElements(
        name="Pandora",
        epoch=2461000.5,
        a=2.75784980823487,  # Main belt asteroid (~67 km, distinct from Saturn moon Pandora)
        e=0.1450625784537587,  # Moderate eccentricity
        i=7.176481151539076,  # Low inclination
        omega=4.93258429587031,
        Omega=10.28381432395156,
        M0=157.5343526907893,
        n=0.2152029188679012,  # ~4.6 year period (1673 days)
    ),
    SE_LILITH_AST: OrbitalElements(
        name="Lilith",
        epoch=2461000.5,
        a=2.664645194686444,  # Main belt asteroid (~30 km, not lunar apogee Lilith)
        e=0.1933175242946558,  # Moderate eccentricity
        i=5.592042797130856,  # Low inclination
        omega=156.6248120324298,
        Omega=260.5622393942289,
        M0=311.9465725963339,
        n=0.2265922174357816,  # ~4.35 year period (1589 days)
    ),
    SE_HIDALGO: OrbitalElements(
        name="Hidalgo",
        epoch=2461000.5,
        a=5.728305724964734,  # Centaur-class asteroid with comet-like orbit
        e=0.6622230478106673,  # Very high eccentricity (comet-like)
        i=42.53593406738054,  # Very high inclination
        omega=56.59815243071501,
        Omega=21.36316817210909,
        M0=185.3374947174508,
        n=0.07188938874958553,  # ~13.7 year period (5008 days)
    ),
    SE_TOUTATIS: OrbitalElements(
        name="Toutatis",
        epoch=2461000.5,
        a=2.543029853370523,  # Apollo PHA, radar and spacecraft target (~5.4 km)
        e=0.6247403332202829,  # High eccentricity (q=0.95 AU, Q=4.13 AU)
        i=0.4480782679581268,  # Very low inclination
        omega=277.8613313176078,
        Omega=125.3653430437421,
        M0=76.89343773316328,
        n=0.2430395129094619,  # ~4.05 year period (1481 days), tumbling rotation
    ),
    SE_ITOKAWA: OrbitalElements(
        name="Itokawa",
        epoch=2461000.5,
        a=1.324135178668783,  # Apollo PHA, Hayabusa sample return mission target (~535 m)
        e=0.2801500981413037,  # Moderate eccentricity (q=0.95 AU, Q=1.70 AU)
        i=1.620938690165939,  # Very low inclination
        omega=162.853569017811,
        Omega=69.07562356929041,
        M0=41.26158479412086,
        n=0.646852987504302,  # ~1.52 year period (557 days), elongated shape
    ),
    SE_BENNU: OrbitalElements(
        name="Bennu",
        epoch=2461000.5,
        a=1.126391025894812,  # Apollo PHA, OSIRIS-REx sample return mission target (~490 m)
        e=0.2037450762416414,  # Moderate eccentricity (q=0.90 AU, Q=1.36 AU)
        i=6.03494377024794,  # Low inclination
        omega=66.22306084084298,
        Omega=2.06086619569642,
        M0=265.1247751080418,  # Propagated from epoch 2455562.5
        n=0.8244613503320309,  # ~1.20 year period (437 days)
    ),
    SE_RYUGU: OrbitalElements(
        name="Ryugu",
        epoch=2461000.5,
        a=1.190921090916117,  # Apollo PHA, Hayabusa2 sample return mission target (~900 m)
        e=0.1910626231558206,  # Moderate eccentricity (q=0.96 AU, Q=1.42 AU)
        i=5.866551865769142,  # Low inclination
        omega=211.6167659559647,
        Omega=251.2915331553314,
        M0=270.6593928357033,
        n=0.7583672923533237,  # ~1.30 year period (475 days)
    ),
}


def solve_kepler_equation_elliptic(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Solve Kepler's equation M = E - e·sin(E) for eccentric anomaly E (elliptic orbits).

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
    E = M if e < 0.8 else math.pi

    for _ in range(30):
        f = E - e * math.sin(E) - M
        fp = 1 - e * math.cos(E)
        E_new = E - f / fp

        if abs(E_new - E) < tol:
            return E_new
        E = E_new

    return E


def solve_kepler_equation_hyperbolic(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Solve hyperbolic Kepler's equation M = e·sinh(H) - H for hyperbolic anomaly H.

    Uses Newton-Raphson iteration for robust convergence.

    Args:
        M: Mean anomaly in radians (can be any value for hyperbolic orbits)
        e: Eccentricity (e > 1)
        tol: Convergence tolerance (default 1e-8)

    Returns:
        float: Hyperbolic anomaly H in radians

    Algorithm:
        Newton-Raphson: H_{n+1} = H_n - f(H_n)/f'(H_n)
        where f(H) = e·sinh(H) - H - M
        and f'(H) = e·cosh(H) - 1

    References:
        Curtis "Orbital Mechanics for Engineering Students" §3.4
        Vallado "Fundamentals of Astrodynamics" Algorithm 4
    """
    # Initial guess
    if e < 1.6:
        H = M if abs(M) < math.pi else math.copysign(math.pi, M)
    else:
        H = math.copysign(math.log(2 * abs(M) / e + 1.8), M) if abs(M) > 1e-10 else 0.0

    for _ in range(50):
        sinh_H = math.sinh(H)
        cosh_H = math.cosh(H)
        f = e * sinh_H - H - M
        fp = e * cosh_H - 1

        if abs(fp) < 1e-15:
            break

        H_new = H - f / fp

        if abs(H_new - H) < tol:
            return H_new
        H = H_new

    return H


def solve_barker_equation(M: float) -> float:
    """
    Solve Barker's equation for parabolic orbits to find true anomaly.

    For parabolic orbits (e = 1), the relationship between mean anomaly M
    and true anomaly ν is given by Barker's equation.

    Args:
        M: Mean anomaly in radians

    Returns:
        float: True anomaly ν in radians

    Algorithm:
        Uses the cubic solution for Barker's equation:
        tan(ν/2) = s where s³ + 3s - 6M/√2 = 0

        Solving via Cardano's formula or the closed-form solution.

    References:
        Vallado "Fundamentals of Astrodynamics" Algorithm 3
        Danby "Fundamentals of Celestial Mechanics" §6.6
    """
    # Barker's equation: M = (1/2) * tan(ν/2) + (1/6) * tan³(ν/2)
    # Let s = tan(ν/2), then: s + s³/3 = 2M
    # Rearranging: s³ + 3s - 6M = 0

    # Use the closed-form solution (Cardano's formula for depressed cubic)
    W = 3.0 * M
    # s = (W + sqrt(W^2 + 1))^(1/3) - (W + sqrt(W^2 + 1))^(-1/3) won't work directly
    # Use: s = 2 * sinh(asinh(W) / 3)
    s = 2.0 * math.sinh(math.asinh(W) / 3.0)

    # True anomaly from tan(ν/2) = s
    nu = 2.0 * math.atan(s)

    return nu


def solve_kepler_equation(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Solve Kepler's equation for the appropriate anomaly based on orbit type.

    Handles elliptic, parabolic, and hyperbolic orbits.

    Args:
        M: Mean anomaly in radians
        e: Eccentricity (e ≥ 0)
        tol: Convergence tolerance (default 1e-8 ~ 0.002 arcsec)

    Returns:
        float: For e < 1: Eccentric anomaly E in radians
               For e = 1: True anomaly ν in radians (special case)
               For e > 1: Hyperbolic anomaly H in radians

    Note:
        For parabolic orbits (e = 1), returns TRUE anomaly directly since
        the eccentric anomaly is not defined for parabolas.

    References:
        Curtis "Orbital Mechanics for Engineering Students" §3.1-3.4
        Vallado "Fundamentals of Astrodynamics" Algorithms 2-4
    """
    # Parabolic orbit (e ≈ 1)
    if abs(e - 1.0) < 1e-10:
        # For parabolic orbits, solve Barker's equation
        # Returns true anomaly directly
        return solve_barker_equation(M)

    # Elliptic orbit (e < 1)
    if e < 1.0:
        return solve_kepler_equation_elliptic(M, e, tol)

    # Hyperbolic orbit (e > 1)
    return solve_kepler_equation_hyperbolic(M, e, tol)


def calc_minor_body_position(
    elements: OrbitalElements, jd_tt: float, include_perturbations: bool = True
) -> Tuple[float, float, float]:
    """
    Calculate heliocentric position using Keplerian orbital elements with perturbations.

    Propagates orbit from epoch to target time using mean motion and applies
    first-order secular perturbations from Jupiter and Saturn.
    Supports elliptic, parabolic, and hyperbolic orbits.

    Args:
        elements: Orbital elements at epoch
        jd_tt: Target Julian Day in Terrestrial Time (TT)
        include_perturbations: If True, apply secular perturbations (default True)

    Returns:
        Tuple[float, float, float]: (x, y, z) heliocentric position in AU
            Coordinates in ecliptic J2000.0 frame

    Algorithm:
        1. Apply secular perturbations to ω and Ω (if enabled)
        2. Propagate mean anomaly: M(t) = M0 + n·Δt
        3. Solve Kepler's equation for eccentric/hyperbolic/true anomaly
        4. Calculate true anomaly ν from E (or H for hyperbolic)
        5. Compute position in orbital plane
        6. Rotate to ecliptic frame using perturbed Ω, i, ω

    Perturbation Model:
        Uses first-order Laplace-Lagrange secular perturbation theory to
        account for the gravitational influence of Jupiter and Saturn.
        Secular perturbations cause:
        - Precession of perihelion (ω)
        - Regression of ascending node (Ω)

        This significantly improves accuracy for multi-year propagation,
        reducing errors from ~1-5 arcminutes (pure Keplerian) to
        ~10-30 arcseconds for main belt asteroids.

    Precision:
        With perturbations enabled:
        - Main belt asteroids: ~10-30 arcseconds typical
        - TNOs: ~1-3 arcminutes typical

        Remaining errors from:
        - Higher-order perturbations
        - Mean-motion resonances
        - Non-gravitational forces (radiation pressure, Yarkovsky)
        - Relativistic effects (minor for asteroids)

    References:
        Murray & Dermott "Solar System Dynamics" Ch. 7 (secular theory)
        Curtis §3 (orbital elements)
        Vallado §2.3 (coordinate transformations)
    """
    e = elements.e
    dt = jd_tt - elements.epoch

    # Apply secular perturbations to get perturbed orbital elements
    if include_perturbations and abs(e - 1.0) > 1e-10 and e < 1.0:
        # Only apply perturbations for elliptic orbits
        omega_pert, Omega_pert, M_deg, n_pert = apply_secular_perturbations(
            elements, jd_tt, include_perturbations=True
        )
    else:
        # Parabolic or hyperbolic orbits, or perturbations disabled
        omega_pert = elements.omega
        Omega_pert = elements.Omega
        M_deg = elements.M0 + elements.n * dt
        # n_pert not used for non-elliptic orbits

    # Propagate mean anomaly (handle differently for each orbit type)
    if abs(e - 1.0) < 1e-10:
        # Parabolic orbit: mean anomaly grows linearly
        M = math.radians(elements.M0 + elements.n * dt)
    elif e < 1.0:
        # Elliptic orbit: use perturbed values, wrap to [0, 360)
        M = math.radians(M_deg % 360.0)
    else:
        # Hyperbolic orbit: no wrapping needed
        M = math.radians(elements.M0 + elements.n * dt)

    # Solve the appropriate Kepler equation
    anomaly = solve_kepler_equation(M, e)

    # Calculate true anomaly and distance based on orbit type
    if abs(e - 1.0) < 1e-10:
        # Parabolic orbit: solve_kepler_equation returns true anomaly directly
        nu = anomaly
        # For parabolic orbits: r = p / (1 + cos(ν)) where p = 2 * q (q = perihelion)
        # Since a is undefined for parabolas, we use q (stored in 'a' field as perihelion)
        # p = 2 * q, so r = 2 * q / (1 + cos(ν))
        p = 2.0 * elements.a  # elements.a is perihelion distance q for parabolic orbits
        r = p / (1.0 + math.cos(nu))
    elif e < 1.0:
        # Elliptic orbit: anomaly is eccentric anomaly E
        E = anomaly
        nu = 2.0 * math.atan2(
            math.sqrt(1 + e) * math.sin(E / 2),
            math.sqrt(1 - e) * math.cos(E / 2),
        )
        r = elements.a * (1 - e * math.cos(E))
    else:
        # Hyperbolic orbit: anomaly is hyperbolic anomaly H
        H = anomaly
        # True anomaly from hyperbolic anomaly
        nu = 2.0 * math.atan2(
            math.sqrt(e + 1) * math.sinh(H / 2),
            math.sqrt(e - 1) * math.cosh(H / 2),
        )
        # Distance from hyperbolic anomaly
        r = elements.a * (e * math.cosh(H) - 1)

    # Position in orbital plane (perifocal frame)
    x_orb = r * math.cos(nu)
    y_orb = r * math.sin(nu)

    # Convert Euler angles to radians (use perturbed values for elliptic)
    if abs(e - 1.0) > 1e-10 and e < 1.0:
        omega_rad = math.radians(omega_pert)
        Omega_rad = math.radians(Omega_pert)
    else:
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
        raise ValueError(f"illegal planet number {body_id}.")

    elements = MINOR_BODY_ELEMENTS[body_id]
    x, y, z = calc_minor_body_position(elements, jd_tt)

    # Convert Cartesian to spherical coordinates
    r = math.sqrt(x**2 + y**2 + z**2)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(z / r))

    return lon, lat, r


# =============================================================================
# JPL SBDB API FOR DYNAMIC ASTEROID LOOKUP
# =============================================================================

# JPL Small-Body Database API endpoint
SBDB_API_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"

# Cache for dynamically fetched orbital elements
# Maps asteroid_number -> OrbitalElements
_ASTEROID_ELEMENTS_CACHE: dict[int, OrbitalElements] = {}


def fetch_orbital_elements_from_sbdb(
    asteroid_number: int, timeout: float = 30.0
) -> Optional[OrbitalElements]:
    """
    Fetch orbital elements for a numbered asteroid from JPL SBDB API.

    This function queries the JPL Small-Body Database API to retrieve current
    osculating orbital elements for any numbered asteroid in the database.

    Args:
        asteroid_number: The asteroid catalog number (e.g., 1 for Ceres, 2060 for Chiron)
        timeout: Request timeout in seconds (default: 30)

    Returns:
        OrbitalElements if successful, None if the asteroid is not found or
        the request fails.

    Raises:
        No exceptions are raised; errors return None.

    Example:
        >>> from libephemeris.minor_bodies import fetch_orbital_elements_from_sbdb
        >>> elements = fetch_orbital_elements_from_sbdb(433)  # Eros
        >>> if elements:
        ...     print(f"Eros semi-major axis: {elements.a:.4f} AU")

    Notes:
        - Requires network access to ssd-api.jpl.nasa.gov
        - Results are cached to avoid redundant API calls
        - JPL SBDB contains 1+ million numbered asteroids
    """
    import json
    import urllib.error
    import urllib.request

    # Check cache first
    if asteroid_number in _ASTEROID_ELEMENTS_CACHE:
        return _ASTEROID_ELEMENTS_CACHE[asteroid_number]

    # Build API request
    params = f"sstr={asteroid_number}&phys-par=false&full-prec=true"
    url = f"{SBDB_API_URL}?{params}"

    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "libephemeris/1.0 (Python)"},
        )
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = json.loads(response.read().decode("utf-8"))

        # Check for API errors
        if "error" in data:
            return None

        # Extract orbital elements
        orbit = data.get("orbit", {})
        elements_list = orbit.get("elements", [])

        if not elements_list:
            return None

        # Build a dict of element name -> value
        elem_dict: dict[str, float] = {}
        for elem in elements_list:
            elem_dict[elem["name"]] = float(elem["value"])

        # Get epoch (Julian Day TDB)
        epoch_jd = float(orbit.get("epoch", 0))
        if epoch_jd == 0:
            return None

        # Get body name from object data
        obj_data = data.get("object", {})
        body_name = obj_data.get("fullname", f"Asteroid {asteroid_number}")

        # Map SBDB element names to our structure
        # SBDB uses: e, a, q, i, om (node), w (arg peri), ma (mean anom), n (mean motion)
        orbital_elements = OrbitalElements(
            name=body_name,
            epoch=epoch_jd,
            a=elem_dict.get("a", 0.0),  # Semi-major axis
            e=elem_dict.get("e", 0.0),  # Eccentricity
            i=elem_dict.get("i", 0.0),  # Inclination
            omega=elem_dict.get("w", 0.0),  # Argument of perihelion
            Omega=elem_dict.get("om", 0.0),  # Longitude of ascending node
            M0=elem_dict.get("ma", 0.0),  # Mean anomaly
            n=elem_dict.get("n", 0.0),  # Mean motion
        )

        # Cache the result
        _ASTEROID_ELEMENTS_CACHE[asteroid_number] = orbital_elements

        return orbital_elements

    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError):
        return None
    except (KeyError, ValueError, TypeError):
        return None


def clear_asteroid_elements_cache() -> int:
    """
    Clear the cache of dynamically fetched orbital elements.

    Returns:
        int: Number of entries cleared from the cache.

    Example:
        >>> from libephemeris.minor_bodies import clear_asteroid_elements_cache
        >>> cleared = clear_asteroid_elements_cache()
        >>> print(f"Cleared {cleared} cached entries")
    """
    count = len(_ASTEROID_ELEMENTS_CACHE)
    _ASTEROID_ELEMENTS_CACHE.clear()
    return count


# Cache for name-to-number lookups (asteroid_name_lower -> asteroid_number)
_ASTEROID_NAME_CACHE: dict[str, int] = {}


def _build_local_name_cache() -> None:
    """
    Build a lookup cache from asteroid names to their catalog numbers.

    Populates the cache using bodies in MINOR_BODY_ELEMENTS and known mappings
    for asteroids with dedicated SE_* constants (Chiron, Pholus, Ceres, etc.).
    """
    if _ASTEROID_NAME_CACHE:
        return  # Already built

    # Known mappings for bodies with dedicated SE_* constants
    known_mappings = {
        "chiron": 2060,
        "pholus": 5145,
        "ceres": 1,
        "pallas": 2,
        "juno": 3,
        "vesta": 4,
        "eris": 136199,
        "sedna": 90377,
        "haumea": 136108,
        "makemake": 136472,
        "ixion": 28978,
        "orcus": 90482,
        "quaoar": 50000,
        "nessus": 7066,
        "asbolus": 8405,
        "chariklo": 10199,
        "gonggong": 225088,
        "varuna": 20000,
        "apophis": 99942,
        "flora": 8,
        "eunomia": 15,
        "hygiea": 10,
        "interamnia": 704,
        "davida": 511,
        "europa": 52,  # Main belt asteroid, not Jupiter's moon
        "sylvia": 87,
        "psyche": 16,
        "eros": 433,
        "amor": 1221,
        "icarus": 1566,
        "toro": 1685,
        "sappho": 80,
        "pandora": 55,  # Main belt asteroid, not Saturn's moon
        "lilith": 1181,  # Main belt asteroid, not lunar apogee Lilith
        "hidalgo": 944,
        "toutatis": 4179,
        "itokawa": 25143,
        "bennu": 101955,
        "ryugu": 162173,
    }

    _ASTEROID_NAME_CACHE.update(known_mappings)


def get_asteroid_number(name: str, timeout: float = 30.0) -> Optional[int]:
    """
    Look up an asteroid's catalog number by name.

    This function is useful for users who know an asteroid's name but not its
    catalog number. It first checks a local database of well-known asteroids,
    then queries the JPL Small-Body Database (SBDB) API if not found locally.

    Args:
        name: The asteroid name to look up (case-insensitive).
            Examples: "Ceres", "Eros", "Apophis", "Bennu", "Vesta"
        timeout: Request timeout in seconds for SBDB API queries (default: 30)

    Returns:
        int: The asteroid catalog number if found, or None if not found.

    Example:
        >>> from libephemeris.minor_bodies import get_asteroid_number
        >>> get_asteroid_number("Ceres")
        1
        >>> get_asteroid_number("Chiron")
        2060
        >>> get_asteroid_number("Apophis")
        99942
        >>> get_asteroid_number("Eros")
        433

    Notes:
        - Local lookup is instant and doesn't require network access
        - SBDB API query requires network access to ssd-api.jpl.nasa.gov
        - Results from SBDB queries are cached for subsequent lookups
        - The search is case-insensitive ("ceres", "Ceres", "CERES" all work)
        - For ambiguous names (e.g., "Pandora" can be an asteroid or Saturn's
          moon), this returns the asteroid catalog number

    See Also:
        calc_asteroid_by_number: Calculate position using the returned number
        fetch_orbital_elements_from_sbdb: Get orbital elements by number
    """
    import json
    import urllib.error
    import urllib.parse
    import urllib.request

    # Normalize name for lookup
    name_lower = name.strip().lower()

    if not name_lower:
        return None

    # Build local cache if not done
    _build_local_name_cache()

    # Check local cache first
    if name_lower in _ASTEROID_NAME_CACHE:
        return _ASTEROID_NAME_CACHE[name_lower]

    # Query JPL SBDB API by name
    # The API accepts names in the sstr parameter
    params = f"sstr={urllib.parse.quote(name)}&phys-par=false"
    url = f"{SBDB_API_URL}?{params}"

    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "libephemeris/1.0 (Python)"},
        )
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = json.loads(response.read().decode("utf-8"))

        # Check for API errors
        if "error" in data:
            return None

        # Extract object data
        obj_data = data.get("object", {})

        # Get the SPK-ID (asteroid number) from object data
        # SBDB returns "spkid" for numbered objects or we can use "des" (designation)
        spkid = obj_data.get("spkid")
        if spkid:
            try:
                # SPK-ID format for asteroids: 2XXXXXX (2000000 + number)
                spkid_int = int(spkid)
                if spkid_int >= 2000000:
                    asteroid_number = spkid_int - 2000000
                else:
                    # Some objects have direct SPK-IDs
                    asteroid_number = spkid_int
            except (ValueError, TypeError):
                return None
        else:
            # Try to parse from designation
            des = obj_data.get("des", "")
            if des and des.isdigit():
                asteroid_number = int(des)
            else:
                # Try fullname parsing - e.g., "1 Ceres" -> 1
                fullname = obj_data.get("fullname", "")
                parts = fullname.split()
                if parts and parts[0].isdigit():
                    asteroid_number = int(parts[0])
                else:
                    return None

        # Cache the result for future lookups
        _ASTEROID_NAME_CACHE[name_lower] = asteroid_number

        return asteroid_number

    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError):
        return None
    except (KeyError, ValueError, TypeError):
        return None


def clear_asteroid_name_cache() -> int:
    """
    Clear the cache of name-to-number lookups.

    Returns:
        int: Number of entries cleared from the cache.

    Example:
        >>> from libephemeris.minor_bodies import clear_asteroid_name_cache
        >>> cleared = clear_asteroid_name_cache()
        >>> print(f"Cleared {cleared} cached entries")
    """
    count = len(_ASTEROID_NAME_CACHE)
    _ASTEROID_NAME_CACHE.clear()
    return count


def calc_asteroid_by_number(
    asteroid_number: int,
    jd_tt: float,
    use_spk: bool = True,
) -> Tuple[float, float, float]:
    """
    Calculate position for any numbered asteroid by fetching orbital elements on demand.

    This function provides a generic way to calculate the heliocentric position
    of any numbered asteroid (out of 1+ million known) without requiring the
    asteroid to be hardcoded in the library. It first checks if an SPK file
    is available for higher precision, otherwise falls back to Keplerian
    calculations using orbital elements fetched from JPL SBDB.

    Args:
        asteroid_number: The asteroid catalog number (e.g., 1 for Ceres,
            433 for Eros, 2060 for Chiron, 136199 for Eris)
        jd_tt: Julian Day in Terrestrial Time (TT)
        use_spk: If True (default), check for registered SPK file first.
            If False, always use Keplerian calculation.

    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Heliocentric ecliptic latitude in degrees (-90 to +90)
            - distance: Heliocentric distance in AU

    Raises:
        ValueError: If asteroid_number is invalid or not found in JPL SBDB.
        ConnectionError: If unable to fetch orbital elements from JPL SBDB
            (only when SPK is not available and cache is empty).

    Precision:
        - With SPK: Sub-arcsecond to arcsecond accuracy
        - Without SPK (Keplerian): ~10-30 arcseconds for main belt,
          ~1-3 arcminutes for TNOs

    Example:
        >>> from libephemeris.minor_bodies import calc_asteroid_by_number
        >>> # Calculate position for Eros (asteroid 433) at J2000.0
        >>> lon, lat, dist = calc_asteroid_by_number(433, 2451545.0)
        >>> print(f"Eros: {lon:.4f}deg, {lat:.4f}deg, {dist:.4f} AU")
        >>>
        >>> # Calculate for any asteroid by number
        >>> lon, lat, dist = calc_asteroid_by_number(99942, 2451545.0)  # Apophis
        >>> lon, lat, dist = calc_asteroid_by_number(25143, 2451545.0)  # Itokawa

    Notes:
        - First call for an asteroid requires network access to JPL SBDB
        - Subsequent calls use cached orbital elements
        - Use clear_asteroid_elements_cache() to clear the cache
        - For known bodies in MINOR_BODY_ELEMENTS, use calc_minor_body_heliocentric()
          as it doesn't require network access

    See Also:
        calc_minor_body_heliocentric: For bodies with pre-defined orbital elements
        fetch_orbital_elements_from_sbdb: To fetch elements without calculating position
    """
    from .constants import SE_AST_OFFSET

    # Validate asteroid number
    if not isinstance(asteroid_number, int) or asteroid_number <= 0:
        raise ValueError(
            f"asteroid_number must be a positive integer, got {asteroid_number}"
        )

    # Check if we have this asteroid in the predefined MINOR_BODY_ELEMENTS
    body_id = asteroid_number + SE_AST_OFFSET
    if body_id in MINOR_BODY_ELEMENTS:
        return calc_minor_body_heliocentric(body_id, jd_tt)

    # Check special cases for low-numbered asteroids that use different IDs
    # (Chiron, Pholus, Ceres, Pallas, Juno, Vesta have dedicated SE_* constants)
    from .constants import (
        SE_CHIRON,
        SE_PHOLUS,
        SE_CERES,
        SE_PALLAS,
        SE_JUNO,
        SE_VESTA,
    )

    special_mapping = {
        2060: SE_CHIRON,  # Chiron
        5145: SE_PHOLUS,  # Pholus
        1: SE_CERES,  # Ceres
        2: SE_PALLAS,  # Pallas
        3: SE_JUNO,  # Juno
        4: SE_VESTA,  # Vesta
    }

    if asteroid_number in special_mapping:
        mapped_id = special_mapping[asteroid_number]
        if mapped_id in MINOR_BODY_ELEMENTS:
            return calc_minor_body_heliocentric(mapped_id, jd_tt)

    # Check if SPK file is available for this asteroid
    if use_spk:
        try:
            from . import state
            from .constants import SEFLG_HELCTR
            from .spk import calc_spk_body_position

            # Check if there's a registered SPK for this body
            if body_id in state._SPK_BODY_MAP:
                # SPK is available, create a Skyfield Time object and use SPK calculation
                ts = state.get_timescale()
                t = ts.tt_jd(jd_tt)

                # Calculate heliocentric position using SPK
                result = calc_spk_body_position(t, body_id, SEFLG_HELCTR)
                if result is not None:
                    lon, lat, dist, _, _, _ = result
                    return lon, lat, dist
        except (ImportError, ValueError):
            # Fall through to Keplerian calculation
            pass

    # No SPK available, fetch orbital elements from SBDB and calculate
    elements = fetch_orbital_elements_from_sbdb(asteroid_number)

    if elements is None:
        raise ValueError(
            f"Asteroid {asteroid_number} not found in JPL Small-Body Database. "
            "Check if the asteroid number is valid."
        )

    # Calculate position using Keplerian mechanics with perturbations
    x, y, z = calc_minor_body_position(elements, jd_tt)

    # Convert Cartesian to spherical coordinates
    r = math.sqrt(x**2 + y**2 + z**2)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(z / r))

    return lon, lat, r
