"""
Hypothetical and fictitious body calculations for libephemeris.

This module provides calculation functions for hypothetical celestial bodies
that are used in various astrological traditions but do not correspond to
known physical objects in the solar system.

Supported Bodies:

Uranian Planets (Hamburg School Astrology):
    The Hamburg School, founded by Alfred Witte in the early 20th century,
    uses eight hypothetical planets beyond Neptune. Their positions are
    calculated using secular polynomial formulas.
    - Cupido: Love, art, marriage, social connections
    - Hades: Decay, hidden matters, the past, investigation
    - Zeus: Controlled energy, machinery, fire, procreation
    - Kronos: Authority, mastery, expertise, government
    - Apollon: Expansion, science, commerce, success
    - Admetos: Endurance, depth, focus, raw materials
    - Vulkanus: Intensity, power, strength, compulsion
    - Poseidon: Spirituality, enlightenment, idealism, truth

Transpluto (Persephone/Isis):
    A hypothetical planet beyond Pluto proposed by various astrologers.
    Its orbit was derived from perturbation analysis of outer planet orbits.
    Position calculated using Keplerian orbital elements.

Other Fictitious Bodies:
    - Vulcan: Hypothetical intra-Mercurial planet (now known not to exist)
    - White Moon (Selena): Second Moon of Earth (symbolic point)
    - Proserpina: Another proposed trans-Plutonian planet
    - Waldemath Black Moon: Second focus of lunar orbit

References:
    - Swiss Ephemeris documentation: "Hypothetical Bodies"
    - Witte, Alfred: "Regelwerk fuer Planetenbilder" (1932)
    - Jacobson: "The Dark Moon Lilith in Astrology" (1961)
    - Landscheidt: "Cosmic Cybernetics" (1973)
"""

import math
from dataclasses import dataclass
from typing import Tuple, Dict, Any
from .constants import SE_FICT_OFFSET


# =============================================================================
# HYPOTHETICAL BODY IDENTIFIERS
# =============================================================================
# Body IDs follow Swiss Ephemeris convention: SE_FICT_OFFSET (40) + index

# Uranian planets (Hamburg School)
SE_CUPIDO: int = SE_FICT_OFFSET + 0  # 40
SE_HADES: int = SE_FICT_OFFSET + 1  # 41
SE_ZEUS: int = SE_FICT_OFFSET + 2  # 42
SE_KRONOS: int = SE_FICT_OFFSET + 3  # 43
SE_APOLLON: int = SE_FICT_OFFSET + 4  # 44
SE_ADMETOS: int = SE_FICT_OFFSET + 5  # 45
SE_VULKANUS: int = SE_FICT_OFFSET + 6  # 46
SE_POSEIDON: int = SE_FICT_OFFSET + 7  # 47

# Other hypothetical bodies
SE_ISIS: int = SE_FICT_OFFSET + 8  # 48 - Transpluto (also called Isis/Persephone)
SE_TRANSPLUTO: int = SE_ISIS  # Alias
SE_NIBIRU: int = SE_FICT_OFFSET + 9  # 49 - Sitchin's hypothetical planet
SE_HARRINGTON: int = SE_FICT_OFFSET + 10  # 50 - Harrington's Planet X
SE_NEPTUNE_LEVERRIER: int = SE_FICT_OFFSET + 11  # 51 - Leverrier's Neptune position
SE_NEPTUNE_ADAMS: int = SE_FICT_OFFSET + 12  # 52 - Adams' Neptune position
SE_PLUTO_LOWELL: int = SE_FICT_OFFSET + 13  # 53 - Lowell's Pluto position
SE_PLUTO_PICKERING: int = SE_FICT_OFFSET + 14  # 54 - Pickering's Pluto position

# Additional hypothetical bodies
SE_VULCAN: int = SE_FICT_OFFSET + 15  # 55 - Intra-Mercurial planet
SE_WHITE_MOON: int = SE_FICT_OFFSET + 16  # 56 - Selena (second Moon of Earth)
SE_PROSERPINA: int = SE_FICT_OFFSET + 17  # 57 - Trans-Plutonian
SE_WALDEMATH: int = SE_FICT_OFFSET + 18  # 58 - Waldemath's second moon

# Pyswisseph-compatible aliases (without SE_ prefix)
CUPIDO: int = SE_CUPIDO
HADES: int = SE_HADES
ZEUS: int = SE_ZEUS
KRONOS: int = SE_KRONOS
APOLLON: int = SE_APOLLON
ADMETOS: int = SE_ADMETOS
VULKANUS: int = SE_VULKANUS
POSEIDON: int = SE_POSEIDON
ISIS: int = SE_ISIS
TRANSPLUTO: int = SE_TRANSPLUTO
NIBIRU: int = SE_NIBIRU
VULCAN: int = SE_VULCAN
WHITE_MOON: int = SE_WHITE_MOON
PROSERPINA: int = SE_PROSERPINA
WALDEMATH: int = SE_WALDEMATH


# =============================================================================
# URANIAN PLANET PARAMETERS (Hamburg School)
# =============================================================================
# Formulas derived from Alfred Witte's work and later refinements.
# Position = L0 + n * T + p * sin(A + B * T)
# where T = Julian centuries from J2000.0 (JD 2451545.0)
# L0 = mean longitude at J2000.0 (degrees)
# n = mean motion (degrees per Julian century)
# p = amplitude of oscillation (degrees)
# A = phase at epoch (degrees)
# B = oscillation frequency (degrees per century)
#
# Note: Different sources give slightly different parameters. These values
# match the Swiss Ephemeris implementation for compatibility.


@dataclass
class UranianElements:
    """
    Orbital elements for Uranian (Hamburg School) hypothetical planets.

    The position is calculated as:
        longitude = L0 + n * T + amplitude * sin(radians(phase + phase_rate * T))

    where T is Julian centuries from J2000.0

    Attributes:
        name: Name of the hypothetical body
        L0: Mean longitude at J2000.0 (degrees)
        n: Mean motion (degrees per Julian century)
        amplitude: Oscillation amplitude (degrees)
        phase: Phase of oscillation at J2000.0 (degrees)
        phase_rate: Rate of phase change (degrees per century)
    """

    name: str
    L0: float
    n: float
    amplitude: float = 0.0
    phase: float = 0.0
    phase_rate: float = 0.0


# Uranian planet parameters from Swiss Ephemeris seorbel.txt
# Based on Witte's original formulas, refined over time
URANIAN_ELEMENTS: Dict[int, UranianElements] = {
    SE_CUPIDO: UranianElements(
        name="Cupido",
        L0=241.2067,  # Mean longitude at J2000.0
        n=1.091437,  # ~0.0109 deg/year = ~330 year period
        amplitude=1.5,
        phase=21.6,
        phase_rate=83.0,
    ),
    SE_HADES: UranianElements(
        name="Hades",
        L0=176.0581,
        n=0.736380,  # ~0.0074 deg/year = ~487 year period
        amplitude=1.5,
        phase=258.0,
        phase_rate=56.0,
    ),
    SE_ZEUS: UranianElements(
        name="Zeus",
        L0=32.1893,
        n=0.532664,  # ~0.0053 deg/year = ~676 year period
        amplitude=1.0,
        phase=141.0,
        phase_rate=41.0,
    ),
    SE_KRONOS: UranianElements(
        name="Kronos",
        L0=213.2096,
        n=0.420481,  # ~0.0042 deg/year = ~857 year period
        amplitude=1.0,
        phase=98.0,
        phase_rate=32.0,
    ),
    SE_APOLLON: UranianElements(
        name="Apollon",
        L0=71.8925,
        n=0.341403,  # ~0.0034 deg/year = ~1055 year period
        amplitude=1.0,
        phase=326.0,
        phase_rate=26.0,
    ),
    SE_ADMETOS: UranianElements(
        name="Admetos",
        L0=142.2269,
        n=0.283756,  # ~0.0028 deg/year = ~1269 year period
        amplitude=1.0,
        phase=250.0,
        phase_rate=22.0,
    ),
    SE_VULKANUS: UranianElements(
        name="Vulkanus",
        L0=195.6753,
        n=0.240116,  # ~0.0024 deg/year = ~1500 year period
        amplitude=1.0,
        phase=178.0,
        phase_rate=18.0,
    ),
    SE_POSEIDON: UranianElements(
        name="Poseidon",
        L0=274.4073,
        n=0.207016,  # ~0.0021 deg/year = ~1739 year period
        amplitude=1.0,
        phase=105.0,
        phase_rate=16.0,
    ),
}


@dataclass
class HypotheticalElements:
    """
    Orbital elements for other hypothetical bodies.

    These use Keplerian or simplified orbital mechanics.

    Attributes:
        name: Name of the hypothetical body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        M0: Mean anomaly at epoch (degrees)
        n: Mean motion (degrees per day)
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


# Other hypothetical body elements
# Transpluto (Isis) elements from Swiss Ephemeris
HYPOTHETICAL_ELEMENTS: Dict[int, HypotheticalElements] = {
    SE_ISIS: HypotheticalElements(
        name="Transpluto/Isis",
        epoch=2451545.0,  # J2000.0
        a=77.775,  # AU - very distant orbit
        e=0.2,  # Moderate eccentricity
        i=0.0,  # Assumed zero inclination
        omega=228.3,  # Argument of perihelion at J2000
        Omega=0.0,  # Assumed zero ascending node
        M0=350.8,  # Mean anomaly at J2000
        n=0.000470,  # Very slow motion (~2100 year period)
    ),
    SE_VULCAN: HypotheticalElements(
        name="Vulcan",
        epoch=2451545.0,  # J2000.0
        a=0.14,  # AU - inside Mercury's orbit
        e=0.0,  # Assumed circular
        i=0.0,  # Assumed zero inclination
        omega=0.0,
        Omega=0.0,
        M0=0.0,  # Unknown, set to Sun position
        n=18.0,  # Very fast motion (~20 day period)
    ),
    SE_PROSERPINA: HypotheticalElements(
        name="Proserpina",
        epoch=2451545.0,
        a=81.0,  # AU
        e=0.15,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=0.0,
        n=0.000445,  # ~2200 year period
    ),
}


# =============================================================================
# HYPOTHETICAL BODY NAME MAPPING
# =============================================================================

HYPOTHETICAL_NAMES: Dict[int, str] = {
    SE_CUPIDO: "Cupido",
    SE_HADES: "Hades",
    SE_ZEUS: "Zeus",
    SE_KRONOS: "Kronos",
    SE_APOLLON: "Apollon",
    SE_ADMETOS: "Admetos",
    SE_VULKANUS: "Vulkanus",
    SE_POSEIDON: "Poseidon",
    SE_ISIS: "Transpluto",
    SE_NIBIRU: "Nibiru",
    SE_HARRINGTON: "Harrington",
    SE_NEPTUNE_LEVERRIER: "Neptune-Leverrier",
    SE_NEPTUNE_ADAMS: "Neptune-Adams",
    SE_PLUTO_LOWELL: "Pluto-Lowell",
    SE_PLUTO_PICKERING: "Pluto-Pickering",
    SE_VULCAN: "Vulcan",
    SE_WHITE_MOON: "White Moon (Selena)",
    SE_PROSERPINA: "Proserpina",
    SE_WALDEMATH: "Waldemath",
}


# =============================================================================
# CALCULATION FUNCTIONS
# =============================================================================


def is_hypothetical_body(ipl: int) -> bool:
    """
    Check if a body ID corresponds to a hypothetical body.

    Args:
        ipl: Planet/body ID

    Returns:
        True if the body is a hypothetical body, False otherwise.

    Example:
        >>> from libephemeris.hypothetical import is_hypothetical_body, SE_CUPIDO
        >>> is_hypothetical_body(SE_CUPIDO)
        True
        >>> is_hypothetical_body(0)  # SE_SUN
        False
    """
    return SE_FICT_OFFSET <= ipl < SE_FICT_OFFSET + 30


def get_hypothetical_name(ipl: int) -> str:
    """
    Get the name of a hypothetical body.

    Args:
        ipl: Planet/body ID

    Returns:
        Name of the hypothetical body, or "Unknown" if not found.

    Example:
        >>> from libephemeris.hypothetical import get_hypothetical_name, SE_CUPIDO
        >>> get_hypothetical_name(SE_CUPIDO)
        'Cupido'
    """
    return HYPOTHETICAL_NAMES.get(ipl, f"Unknown ({ipl})")


def calc_uranian_longitude(ipl: int, jd_tt: float) -> float:
    """
    Calculate the ecliptic longitude of a Uranian (Hamburg School) planet.

    Uses the secular polynomial formula:
        longitude = L0 + n * T + amplitude * sin(radians(phase + phase_rate * T))

    where T is Julian centuries from J2000.0

    Args:
        ipl: Uranian planet ID (SE_CUPIDO through SE_POSEIDON)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Ecliptic longitude in degrees (0-360)

    Raises:
        ValueError: If ipl is not a valid Uranian planet ID.

    Example:
        >>> from libephemeris.hypothetical import calc_uranian_longitude, SE_CUPIDO
        >>> lon = calc_uranian_longitude(SE_CUPIDO, 2451545.0)  # J2000.0
        >>> print(f"Cupido longitude: {lon:.4f}")
    """
    if ipl not in URANIAN_ELEMENTS:
        raise ValueError(f"Body ID {ipl} is not a valid Uranian planet")

    elements = URANIAN_ELEMENTS[ipl]

    # Calculate T = Julian centuries from J2000.0
    T = (jd_tt - 2451545.0) / 36525.0

    # Mean longitude
    longitude = elements.L0 + elements.n * T

    # Add periodic oscillation if amplitude is non-zero
    if elements.amplitude != 0.0:
        arg = elements.phase + elements.phase_rate * T
        longitude += elements.amplitude * math.sin(math.radians(arg))

    # Normalize to 0-360 degrees
    longitude = longitude % 360.0

    return longitude


def calc_uranian_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the full position of a Uranian (Hamburg School) planet.

    Returns ecliptic longitude, latitude, distance, and their rates of change.
    Uranian planets are assumed to be on the ecliptic (latitude = 0) and at
    a fixed distance (based on their orbital period).

    Args:
        ipl: Uranian planet ID (SE_CUPIDO through SE_POSEIDON)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees
            - latitude: Ecliptic latitude in degrees (always 0)
            - distance: Distance in AU (estimated from mean motion)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0)

    Raises:
        ValueError: If ipl is not a valid Uranian planet ID.

    Example:
        >>> from libephemeris.hypothetical import calc_uranian_position, SE_KRONOS
        >>> pos = calc_uranian_position(SE_KRONOS, 2451545.0)
        >>> print(f"Kronos at {pos[0]:.4f} deg, dist {pos[2]:.2f} AU")
    """
    if ipl not in URANIAN_ELEMENTS:
        raise ValueError(f"Body ID {ipl} is not a valid Uranian planet")

    elements = URANIAN_ELEMENTS[ipl]

    # Calculate longitude
    longitude = calc_uranian_longitude(ipl, jd_tt)

    # Latitude is assumed to be 0 (on ecliptic)
    latitude = 0.0

    # Estimate distance from mean motion using Kepler's 3rd law
    # n (deg/century) -> n (deg/day) = n / 36525
    # n (rad/day) = n (deg/day) * pi / 180
    # For circular orbit: n = sqrt(GM / a^3), where GM = k^2 (Gaussian constant)
    # k = 0.01720209895 AU^(3/2) / day
    # n = k / a^(3/2)  => a = (k / n)^(2/3)
    n_deg_per_day = elements.n / 36525.0
    n_rad_per_day = math.radians(n_deg_per_day)
    k = 0.01720209895  # Gaussian gravitational constant

    if n_rad_per_day > 0:
        distance = (k / n_rad_per_day) ** (2.0 / 3.0)
    else:
        distance = 100.0  # Default large distance

    # Calculate velocity (numerical differentiation)
    dt = 1.0 / 86400.0  # 1 second in days
    lon_next = calc_uranian_longitude(ipl, jd_tt + dt)

    # Handle wrap-around
    dlon = (lon_next - longitude) / dt
    if dlon > 180.0 / dt:
        dlon -= 360.0 / dt
    elif dlon < -180.0 / dt:
        dlon += 360.0 / dt

    # Latitude and distance are constant
    dlat = 0.0
    ddist = 0.0

    return (longitude, latitude, distance, dlon, dlat, ddist)


def calc_transpluto_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Transpluto (Isis/Persephone).

    Uses Keplerian orbital mechanics with the elements stored in
    HYPOTHETICAL_ELEMENTS[SE_ISIS].

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_transpluto_position
        >>> pos = calc_transpluto_position(2451545.0)  # J2000.0
        >>> print(f"Transpluto at {pos[0]:.4f} deg")
    """
    return _calc_keplerian_position(SE_ISIS, jd_tt)


def _solve_kepler_equation(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Solve Kepler's equation M = E - e * sin(E) for eccentric anomaly E.

    Uses Newton-Raphson iteration.

    Args:
        M: Mean anomaly in radians
        e: Eccentricity (0 <= e < 1)
        tol: Convergence tolerance

    Returns:
        Eccentric anomaly E in radians
    """
    # Initial guess
    E = M if e < 0.8 else math.pi

    for _ in range(30):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)
        E_new = E - f / fp

        if abs(E_new - E) < tol:
            return E_new
        E = E_new

    return E


def _calc_keplerian_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate position using Keplerian orbital elements.

    Internal function for hypothetical bodies with Keplerian orbits.

    Args:
        ipl: Body ID
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
    """
    if ipl not in HYPOTHETICAL_ELEMENTS:
        raise ValueError(f"Body ID {ipl} does not have Keplerian elements")

    elements = HYPOTHETICAL_ELEMENTS[ipl]

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = elements.M0 + elements.n * dt
    M_rad = math.radians(M % 360.0)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    u = nu + math.radians(elements.omega)

    # Convert to ecliptic coordinates
    # For zero inclination and ascending node, this simplifies considerably
    i_rad = math.radians(elements.i)
    Omega_rad = math.radians(elements.Omega)

    # Position in orbital plane
    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    # Rotate to ecliptic frame
    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    # Convert to spherical coordinates
    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    # Calculate velocity via numerical differentiation
    dt_step = 1.0 / 86400.0  # 1 second
    pos_next = _calc_keplerian_position_raw(ipl, jd_tt + dt_step)

    dlon = (pos_next[0] - longitude) / dt_step
    # Handle wrap-around
    if dlon > 180.0 / dt_step:
        dlon -= 360.0 / dt_step
    elif dlon < -180.0 / dt_step:
        dlon += 360.0 / dt_step

    dlat = (pos_next[1] - latitude) / dt_step
    ddist = (pos_next[2] - distance) / dt_step

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_keplerian_position_raw(ipl: int, jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Keplerian position without velocity (helper for differentiation).
    """
    if ipl not in HYPOTHETICAL_ELEMENTS:
        raise ValueError(f"Body ID {ipl} does not have Keplerian elements")

    elements = HYPOTHETICAL_ELEMENTS[ipl]

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = elements.M0 + elements.n * dt
    M_rad = math.radians(M % 360.0)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude
    u = nu + math.radians(elements.omega)

    # Convert to ecliptic coordinates
    i_rad = math.radians(elements.i)
    Omega_rad = math.radians(elements.Omega)

    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    return (longitude, latitude, distance)


def calc_white_moon_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of the White Moon (Selena).

    The White Moon is a symbolic point representing the second Moon of Earth.
    Its position is calculated as the point 180 degrees opposite to the
    Black Moon Lilith (mean lunar apogee).

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)

    Note:
        This is a symbolic calculation. The White Moon is not a physical body.
    """
    # Import lunar module to get mean Lilith position
    from . import lunar

    # Mean Lilith longitude
    lilith_lon = lunar.calc_mean_lilith(jd_tt)

    # White Moon is opposite to Black Moon
    longitude = (lilith_lon + 180.0) % 360.0

    # Latitude, distance, and velocities are simplified
    latitude = 0.0
    distance = 0.0  # Symbolic point, no meaningful distance

    # Velocity is opposite to Lilith velocity
    dt = 1.0 / 86400.0
    lilith_next = lunar.calc_mean_lilith(jd_tt + dt)
    lilith_dlon = (lilith_next - lilith_lon) / dt
    if lilith_dlon > 180.0 / dt:
        lilith_dlon -= 360.0 / dt
    elif lilith_dlon < -180.0 / dt:
        lilith_dlon += 360.0 / dt

    dlon = lilith_dlon  # Same rate, just offset by 180 degrees
    dlat = 0.0
    ddist = 0.0

    return (longitude, latitude, distance, dlon, dlat, ddist)


def calc_waldemath_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of the Waldemath Black Moon.

    The Waldemath Black Moon is the second focus of the lunar orbit ellipse.
    Its position is calculated as twice the Earth-Moon distance in the
    direction opposite to the Moon.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)

    Note:
        This is a geometric point (second focus of ellipse), not a physical body.
    """
    # Import lunar module
    from . import lunar

    # Get mean apogee (Lilith) position - this approximates the apse line direction
    lilith_lon = lunar.calc_mean_lilith(jd_tt)

    # Waldemath is in the same direction as Lilith (apogee direction)
    # but represents the second focus of the ellipse
    longitude = lilith_lon

    # For an ellipse with semi-major axis a and eccentricity e,
    # the distance from center to focus is c = a * e
    # The lunar orbit has a ~ 384,400 km, e ~ 0.0549
    # Distance from Earth to second focus ~ 2 * c ~ 42,200 km
    # In AU: 42,200 km / 149,597,870.7 km/AU ~ 0.00028 AU
    distance = 0.00028  # AU (geometric, not meaningful for astrology)

    latitude = 0.0

    # Calculate velocity
    dt = 1.0 / 86400.0
    lilith_next = lunar.calc_mean_lilith(jd_tt + dt)
    dlon = (lilith_next - lilith_lon) / dt
    if dlon > 180.0 / dt:
        dlon -= 360.0 / dt
    elif dlon < -180.0 / dt:
        dlon += 360.0 / dt

    dlat = 0.0
    ddist = 0.0

    return (longitude, latitude, distance, dlon, dlat, ddist)


def calc_hypothetical_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of any hypothetical body.

    This is the main entry point for hypothetical body calculations.
    Routes the request to the appropriate calculation function based on body ID.

    Args:
        ipl: Hypothetical body ID (SE_CUPIDO through SE_WALDEMATH)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees
            - latitude: Ecliptic latitude in degrees
            - distance: Distance in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Raises:
        ValueError: If ipl is not a valid hypothetical body ID.

    Example:
        >>> from libephemeris.hypothetical import calc_hypothetical_position, SE_CUPIDO
        >>> pos = calc_hypothetical_position(SE_CUPIDO, 2451545.0)
        >>> print(f"Cupido: {pos[0]:.4f} deg")
    """
    if not is_hypothetical_body(ipl):
        raise ValueError(f"Body ID {ipl} is not a hypothetical body")

    # Uranian planets (Hamburg School)
    if ipl in URANIAN_ELEMENTS:
        return calc_uranian_position(ipl, jd_tt)

    # Transpluto / Isis
    if ipl == SE_ISIS:
        return calc_transpluto_position(jd_tt)

    # White Moon (Selena)
    if ipl == SE_WHITE_MOON:
        return calc_white_moon_position(jd_tt)

    # Waldemath Black Moon
    if ipl == SE_WALDEMATH:
        return calc_waldemath_position(jd_tt)

    # Other Keplerian bodies
    if ipl in HYPOTHETICAL_ELEMENTS:
        return _calc_keplerian_position(ipl, jd_tt)

    # Unknown hypothetical body - return zeros
    return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def list_hypothetical_bodies() -> Dict[int, str]:
    """
    Get a dictionary of all supported hypothetical body IDs and names.

    Returns:
        Dictionary mapping body ID to body name.

    Example:
        >>> from libephemeris.hypothetical import list_hypothetical_bodies
        >>> bodies = list_hypothetical_bodies()
        >>> for id, name in bodies.items():
        ...     print(f"{id}: {name}")
    """
    return HYPOTHETICAL_NAMES.copy()
