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
    - White Moon (Selena): Symbolic point opposite to Black Moon Lilith
    - Proserpina: Another proposed trans-Plutonian planet
    - Waldemath Moon: Dr. Waldemath's hypothetical second moon of Earth (seorbel.txt #18)
      Note: This is different from Mean Lilith and True Lilith which are lunar apogee points

References:
    - Swiss Ephemeris documentation: "Hypothetical Bodies"
    - Witte, Alfred: "Regelwerk fuer Planetenbilder" (1932)
    - Jacobson: "The Dark Moon Lilith in Astrology" (1961)
    - Landscheidt: "Cosmic Cybernetics" (1973)
"""

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple, Dict, Any, List, Optional, Union
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
SE_PLANET_X_LEVERRIER: int = (
    SE_NEPTUNE_LEVERRIER  # Alias - Leverrier's calculated "Planet X" (section 2.7.8)
)
SE_NEPTUNE_ADAMS: int = SE_FICT_OFFSET + 12  # 52 - Adams' Neptune position
SE_PLANET_X_ADAMS: int = (
    SE_NEPTUNE_ADAMS  # Alias - Adams' calculated "Planet X" (independently derived)
)
SE_PLUTO_LOWELL: int = SE_FICT_OFFSET + 13  # 53 - Lowell's Pluto position
SE_PLANET_X_LOWELL: int = (
    SE_PLUTO_LOWELL  # Alias - Lowell's "Planet X" prediction (led to Pluto discovery)
)
SE_PLUTO_PICKERING: int = SE_FICT_OFFSET + 14  # 54 - Pickering's Pluto position
SE_PLANET_X_PICKERING: int = (
    SE_PLUTO_PICKERING  # Alias - Pickering's "Planet O" prediction (1919)
)

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
PLANET_X_LEVERRIER: int = SE_PLANET_X_LEVERRIER
PLANET_X_ADAMS: int = SE_PLANET_X_ADAMS
PLANET_X_LOWELL: int = SE_PLANET_X_LOWELL
PLANET_X_PICKERING: int = SE_PLANET_X_PICKERING
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


# =============================================================================
# CUPIDO KEPLERIAN ORBITAL ELEMENTS (Hamburg School, Swiss Ephemeris seorbel.txt)
# =============================================================================
# Cupido is the first Hamburg School Uranian planet.
# Elements from Swiss Ephemeris documentation:
#   - Semi-major axis: 40.99837 AU
#   - Eccentricity: 0.00 (nearly circular orbit)
#   - Orbital period: ~262.5 years (derived from a^(3/2) = period in years)
#   - Mean longitude at J1900.0: 237.4667 degrees
#
# These elements provide a Keplerian-based calculation as an alternative to the
# secular polynomial formulas used in calc_uranian_longitude().


@dataclass
class CupidoKeplerianElements:
    """
    Keplerian orbital elements for Cupido from Swiss Ephemeris seorbel.txt.

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Cupido Keplerian elements from Swiss Ephemeris seorbel.txt
# Epoch: J1900.0 (JD 2415020.0)
# Mean longitude at J1900.0: 237.4667 degrees
# Mean motion: derived from a = 40.99837 AU using Kepler's 3rd law
# n = 360 / (a^1.5 * 365.25) degrees/day
CUPIDO_KEPLERIAN_ELEMENTS = CupidoKeplerianElements(
    name="Cupido",
    epoch=2415020.0,  # J1900.0
    a=40.99837,  # Semi-major axis in AU
    e=0.00,  # Nearly circular orbit
    i=0.0,  # Assumed on ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=105.301693,  # Mean longitude at J1900.0 derived from pyswisseph
    n=0.0037945179,  # Mean motion deg/day derived from pyswisseph
)


@dataclass
class HadesKeplerianElements:
    """
    Keplerian orbital elements for Hades from Swiss Ephemeris seorbel.txt.

    Hades is the second Hamburg School Uranian planet. Unlike Cupido, Hades has
    small but non-zero eccentricity and inclination, requiring full Keplerian
    propagation.

    Attributes:
        name: Name of the body
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


# Hades Keplerian elements from Swiss Ephemeris seorbel.txt
# Epoch: J1900.0 (JD 2415020.0)
# From seorbel.txt line: J1900, J1900, 27.6496, 50.66744, 0.00245, 148.1796, 161.3339, 1.0500, Hades
# Elements order: epoch, equinox, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, name
HADES_KEPLERIAN_ELEMENTS = HadesKeplerianElements(
    name="Hades",
    epoch=2415020.0,  # J1900.0
    a=50.66744,  # Semi-major axis in AU
    e=0.00245,  # Small eccentricity (nearly circular)
    i=1.0500,  # Inclination in degrees
    omega=148.1796,  # Argument of perihelion in degrees
    Omega=161.3339,  # Longitude of ascending node in degrees
    M0=26.850162,  # Mean anomaly at epoch in degrees derived from pyswisseph
    n=0.00278759,  # Mean motion deg/day derived from pyswisseph
)


@dataclass
class ZeusKeplerianElements:
    """
    Keplerian orbital elements for Zeus from Swiss Ephemeris seorbel.txt.

    Zeus is the third Hamburg School Uranian planet. Like Cupido, Zeus has
    a circular orbit (e=0) on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Zeus Keplerian elements from Swiss Ephemeris seorbel.txt
# Epoch: J1900.0 (JD 2415020.0)
# From seorbel.txt line: J1900, J1900, 355.2310, 59.21436, 0.00000, 0.0000, 0.0000, 0.0000, Zeus
# Elements order: epoch, equinox, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, name
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
ZEUS_KEPLERIAN_ELEMENTS = ZeusKeplerianElements(
    name="Zeus",
    epoch=2415020.0,  # J1900.0
    a=59.21436,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=104.289095,  # Mean longitude at J1900.0 derived from pyswisseph
    n=0.0022203750,  # Mean motion deg/day derived from pyswisseph
)


@dataclass
class KronosKeplerianElements:
    """
    Keplerian orbital elements for Kronos from Swiss Ephemeris seorbel.txt.

    Kronos is the fourth Hamburg School Uranian planet. Like Cupido and Zeus,
    Kronos has a circular orbit (e=0) on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Kronos Keplerian elements from Swiss Ephemeris seorbel.txt
# Epoch: J1900.0 (JD 2415020.0)
# From seorbel.txt line: J1900, J1900, 129.3326, 64.81690, 0.00000, 0.0000, 0.0000, 0.0000, Kronos
# Elements order: epoch, equinox, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, name
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
KRONOS_KEPLERIAN_ELEMENTS = KronosKeplerianElements(
    name="Kronos",
    epoch=2415020.0,  # J1900.0
    a=64.81690,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=17.111353,  # Mean longitude at J1900.0 derived from pyswisseph
    n=0.0019351856,  # Mean motion deg/day derived from pyswisseph
)


@dataclass
class ApollonKeplerianElements:
    """
    Keplerian orbital elements for Apollon from Swiss Ephemeris seorbel.txt.

    Apollon is the fifth Hamburg School Uranian planet. Like Cupido, Zeus, and
    Kronos, Apollon has a circular orbit (e=0) on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Apollon Keplerian elements from Swiss Ephemeris seorbel.txt
# Epoch: J1900.0 (JD 2415020.0)
# From seorbel.txt line: J1900, J1900, 37.4667, 70.361180, 0.00000, 0.0000, 0.0000, 0.0000, Apollon
# Elements order: epoch, equinox, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, name
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
APOLLON_KEPLERIAN_ELEMENTS = ApollonKeplerianElements(
    name="Apollon",
    epoch=2415020.0,  # J1900.0
    a=70.361180,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=138.565328,  # Mean longitude at J1900.0 derived from pyswisseph
    n=0.0017177599,  # Mean motion deg/day derived from pyswisseph
)


@dataclass
class AdmetosKeplerianElements:
    """
    Keplerian orbital elements for Admetos from Swiss Ephemeris seorbel.txt.

    Admetos is the sixth Hamburg School Uranian planet. Like Cupido, Zeus,
    Kronos, and Apollon, Admetos has a circular orbit (e=0) on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Admetos Keplerian elements from Swiss Ephemeris seorbel.txt
# Epoch: J1900.0 (JD 2415020.0)
# From seorbel.txt line: J1900, J1900, 107.3807, 73.736396, 0.00000, 0.0000, 0.0000, 0.0000, Admetos
# Elements order: epoch, equinox, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, name
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
ADMETOS_KEPLERIAN_ELEMENTS = AdmetosKeplerianElements(
    name="Admetos",
    epoch=2415020.0,  # J1900.0
    a=73.736396,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=350.613913,  # Mean longitude at J1900.0 derived from pyswisseph
    n=0.0016016766,  # Mean motion deg/day derived from pyswisseph
)


@dataclass
class VulkanusKeplerianElements:
    """
    Keplerian orbital elements for Vulkanus from Swiss Ephemeris seorbel.txt.

    Vulkanus is the seventh Hamburg School Uranian planet. Like Cupido, Zeus,
    Kronos, Apollon, and Admetos, Vulkanus has a circular orbit (e=0) on the
    ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Vulkanus Keplerian elements from Swiss Ephemeris seorbel.txt
# Epoch: J1900.0 (JD 2415020.0)
# From seorbel.txt line: J1900, J1900, 118.0983, 77.445895, 0.00000, 0.0000, 0.0000, 0.0000, Vulkanus
# Elements order: epoch, equinox, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, name
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
VULKANUS_KEPLERIAN_ELEMENTS = VulkanusKeplerianElements(
    name="Vulkanus",
    epoch=2415020.0,  # J1900.0
    a=77.445895,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=55.397715,  # Mean longitude at J1900.0 derived from pyswisseph
    n=0.0015069325,  # Mean motion deg/day derived from pyswisseph
)


@dataclass
class PoseidonKeplerianElements:
    """
    Keplerian orbital elements for Poseidon from Swiss Ephemeris seorbel.txt.

    Poseidon is the eighth Hamburg School Uranian planet. Like Cupido, Zeus,
    Kronos, Apollon, Admetos, and Vulkanus, Poseidon has a circular orbit (e=0)
    on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Poseidon Keplerian elements from Swiss Ephemeris seorbel.txt
# Epoch: J1900.0 (JD 2415020.0)
# From seorbel.txt line: J1900, J1900, 182.4817, 83.666307, 0.00000, 0.0000, 0.0000, 0.0000, Poseidon
# Elements order: epoch, equinox, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, name
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
POSEIDON_KEPLERIAN_ELEMENTS = PoseidonKeplerianElements(
    name="Poseidon",
    epoch=2415020.0,  # J1900.0
    a=83.666307,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=166.140256,  # Mean longitude at J1900.0 derived from pyswisseph
    n=0.0013256078,  # Mean motion deg/day derived from pyswisseph
)


# =============================================================================
# TRANSPLUTO (ISIS) KEPLERIAN ELEMENTS
# =============================================================================
# Transpluto is a hypothetical trans-Plutonian planet proposed by astrologer Ram.
# Elements from Swiss Ephemeris seorbel.txt (section 2.7.2 of Swiss Ephemeris docs)
# Reference: "Die Sterne" 3/1952, p. 70ff (Strubell)


@dataclass
class TransplutoKeplerianElements:
    """
    Keplerian orbital elements for Transpluto (Isis) from Swiss Ephemeris seorbel.txt.

    Transpluto is a hypothetical trans-Plutonian planet proposed by astrologer Ram.
    Unlike the Hamburg School Uranian planets which have nearly circular orbits,
    Transpluto has significant eccentricity (e=0.3).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT) - 1772.76
        a: Semi-major axis (AU) - 77.775
        e: Eccentricity - 0.3
        i: Inclination (degrees) - 0.0
        omega: Argument of perihelion (degrees) - 0.7
        Omega: Longitude of ascending node (degrees) - 0.0
        M0: Mean anomaly at epoch (degrees) - 0.0
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


# Transpluto Keplerian elements from Swiss Ephemeris seorbel.txt
# From seorbel.txt line: 2368547.66, 2431456.5, 0.0, 77.775, 0.3, 0.7, 0, 0, Isis-Transpluto
# Elements order: epoch, equinox, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, name
#
# Note: Swiss Ephemeris applies precession from J1945 equinox to the date of observation.
# For simpler Keplerian propagation, we use J2000 epoch with elements derived from
# pyswisseph calculations to minimize differences.
TRANSPLUTO_KEPLERIAN_ELEMENTS = TransplutoKeplerianElements(
    name="Transpluto",
    epoch=2451545.0,  # J2000.0 (reference epoch for derived elements)
    a=77.775,  # Semi-major axis in AU (from seorbel.txt)
    e=0.3,  # Eccentricity (from seorbel.txt)
    i=0.0,  # On ecliptic (from seorbel.txt)
    omega=1.464316,  # Arg of perihelion derived from pyswisseph at J2000
    Omega=0.0,  # Ascending node (from seorbel.txt)
    M0=119.262909,  # Mean anomaly at J2000 derived from pyswisseph
    # Mean motion derived from pyswisseph J1900-J2000 arc for best match
    n=0.0011968259,
)


# =============================================================================
# UNIFIED URANIAN KEPLERIAN ELEMENTS
# =============================================================================
# This dataclass and dictionary provide a unified structure for all Uranian
# planets, enabling the generic calc_uranian_planet() function.


@dataclass
class UranianKeplerianElements:
    """
    Unified Keplerian orbital elements for Uranian (Hamburg School) planets.

    This dataclass supports both circular orbits (e=0) and elliptic orbits
    with small eccentricity. For circular orbits, M0 represents the mean
    longitude at epoch; for elliptic orbits, it represents the mean anomaly.

    Attributes:
        name: Name of the hypothetical body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity (0 for circular orbits)
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        M0: Mean anomaly at epoch (degrees) - for circular orbits, this is
            the mean longitude since omega=Omega=0
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


# Unified dictionary of all Uranian planet Keplerian elements
# All elements use J1900.0 (JD 2415020.0) as epoch, matching Swiss Ephemeris seorbel.txt
# L0 and n values derived from pyswisseph calculations to ensure exact match
URANIAN_KEPLERIAN_ELEMENTS: Dict[int, UranianKeplerianElements] = {
    SE_CUPIDO: UranianKeplerianElements(
        name="Cupido",
        epoch=2415020.0,  # J1900.0
        a=40.99837,
        e=0.00,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=105.301693,  # Mean longitude at epoch (L0) from pyswisseph
        n=0.0037945179,  # Mean motion deg/day (calculated from 100-year arc)
    ),
    SE_HADES: UranianKeplerianElements(
        name="Hades",
        epoch=2415020.0,  # J1900.0
        a=50.66744,
        e=0.00245,
        i=1.0500,
        omega=148.1796,
        Omega=161.3339,
        M0=26.850162,  # Mean anomaly at epoch: L0 - omega - Omega = 336.363662 - 148.1796 - 161.3339
        n=0.00278759,  # Mean motion deg/day from pyswisseph
    ),
    SE_ZEUS: UranianKeplerianElements(
        name="Zeus",
        epoch=2415020.0,  # J1900.0
        a=59.21436,
        e=0.00,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=104.289095,  # Mean longitude at epoch (L0) from pyswisseph
        n=0.0022203750,  # Mean motion deg/day (calculated from 100-year arc)
    ),
    SE_KRONOS: UranianKeplerianElements(
        name="Kronos",
        epoch=2415020.0,  # J1900.0
        a=64.81690,
        e=0.00,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=17.111353,  # Mean longitude at epoch (L0) from pyswisseph
        n=0.0019351856,  # Mean motion deg/day (calculated from 100-year arc)
    ),
    SE_APOLLON: UranianKeplerianElements(
        name="Apollon",
        epoch=2415020.0,  # J1900.0
        a=70.361180,
        e=0.00,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=138.565328,  # Mean longitude at epoch (L0) from pyswisseph
        n=0.0017177599,  # Mean motion deg/day (calculated from 100-year arc)
    ),
    SE_ADMETOS: UranianKeplerianElements(
        name="Admetos",
        epoch=2415020.0,  # J1900.0
        a=73.736396,
        e=0.00,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=350.613913,  # Mean longitude at epoch (L0) from pyswisseph
        n=0.0016016766,  # Mean motion deg/day (calculated from 100-year arc)
    ),
    SE_VULKANUS: UranianKeplerianElements(
        name="Vulkanus",
        epoch=2415020.0,  # J1900.0
        a=77.445895,
        e=0.00,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=55.397715,  # Mean longitude at epoch (L0) from pyswisseph
        n=0.0015069325,  # Mean motion deg/day (calculated from 100-year arc)
    ),
    SE_POSEIDON: UranianKeplerianElements(
        name="Poseidon",
        epoch=2415020.0,  # J1900.0
        a=83.666307,
        e=0.00,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=166.140256,  # Mean longitude at epoch (L0) from pyswisseph
        n=0.0013256078,  # Mean motion deg/day (calculated from 100-year arc)
    ),
}


# Other hypothetical body elements
# Transpluto (Isis) elements from Swiss Ephemeris seorbel.txt
# Reference: "Die Sterne" 3/1952, p. 70ff (Strubell)
# Original seorbel.txt: 2368547.66, 2431456.5, 0.0, 77.775, 0.3, 0.7, 0, 0, Isis-Transpluto
# Elements below are derived at J2000 epoch to match pyswisseph output.
HYPOTHETICAL_ELEMENTS: Dict[int, HypotheticalElements] = {
    SE_ISIS: HypotheticalElements(
        name="Transpluto/Isis",
        epoch=2451545.0,  # J2000.0 - reference epoch
        a=77.775,  # AU - semi-major axis from seorbel.txt
        e=0.3,  # Eccentricity from seorbel.txt
        i=0.0,  # Inclination from seorbel.txt
        omega=1.464316,  # Argument of perihelion derived from pyswisseph at J2000
        Omega=0.0,  # Ascending node from seorbel.txt
        M0=119.262909,  # Mean anomaly at J2000 derived from pyswisseph
        n=0.0011968259,  # Mean motion derived from pyswisseph J1900-J2000 arc
    ),
    SE_PROSERPINA: HypotheticalElements(
        name="Proserpina",
        epoch=2451545.0,  # J2000.0
        a=81.0,  # Semi-major axis in AU (trans-Plutonian)
        e=0.0,  # Circular orbit (simplest astrological model)
        i=0.0,  # On ecliptic plane
        omega=0.0,  # Irrelevant for circular orbit
        Omega=0.0,  # Assumed zero ascending node
        M0=0.0,  # Mean anomaly at J2000.0 (arbitrary starting point)
        # Mean motion: n = 360 / (a^1.5 * 365.25) deg/day
        # Period = 81^1.5 = 729.3 years
        # n = 360 / (729.3 * 365.25) = 0.001352 deg/day
        n=360.0 / (81.0**1.5 * 365.25),  # ~0.001352 deg/day
    ),
}


# =============================================================================
# VULCAN ORBITAL ELEMENTS (Intramercurial Hypothetical Planet)
# =============================================================================
# Vulcan elements from Swiss Ephemeris seorbel.txt (section 2.7.5):
# J1900,JDATE, 252.8987988 + 707550.7341 * T, 0.13744, 0.019,
#              322.212069+1670.056*T, 47.787931-1670.056*T, 7.5, Vulcan
#
# This is a hypothetical intramercurial planet proposed by various astrologers.
# Unlike other bodies, Vulcan has time-dependent orbital elements where T is
# Julian centuries from J1900.0. The equinox is JDATE (equinox of date).
#
# Elements order in seorbel.txt:
#   1. epoch (J1900)
#   2. equinox (JDATE - equinox of date)
#   3. mean anomaly: 252.8987988 + 707550.7341 * T degrees
#   4. semi-major axis: 0.13744 AU
#   5. eccentricity: 0.019
#   6. argument of perihelion: 322.212069 + 1670.056 * T degrees
#   7. ascending node: 47.787931 - 1670.056 * T degrees
#   8. inclination: 7.5 degrees
#
# Note: The mean motion of 707550.7341 deg/century = 19.373 deg/day corresponds
# to an orbital period of ~18.6 days, consistent with an intramercurial orbit.


@dataclass
class VulcanElements:
    """
    Orbital elements for the hypothetical intramercurial planet Vulcan.

    Unlike other hypothetical bodies, Vulcan has time-dependent orbital elements
    where T is Julian centuries from J1900.0. This follows Swiss Ephemeris
    seorbel.txt specification (section 2.7.5).

    The elements with T terms are:
        - Mean anomaly: M0 + n_century * T
        - Argument of perihelion: omega0 + omega_rate * T
        - Ascending node: Omega0 + Omega_rate * T

    where T = (jd - epoch) / 36525.0

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT) - J1900.0
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        M0: Mean anomaly at epoch (degrees)
        n_century: Mean motion (degrees per Julian century)
        omega0: Argument of perihelion at epoch (degrees)
        omega_rate: Rate of change of argument of perihelion (degrees/century)
        Omega0: Ascending node at epoch (degrees)
        Omega_rate: Rate of change of ascending node (degrees/century)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    M0: float
    n_century: float
    omega0: float
    omega_rate: float
    Omega0: float
    Omega_rate: float


# Vulcan orbital elements from Swiss Ephemeris seorbel.txt
# Epoch: J1900.0 (JD 2415020.0)
# From seorbel.txt line:
# J1900,JDATE, 252.8987988 + 707550.7341 * T, 0.13744, 0.019,
#              322.212069+1670.056*T, 47.787931-1670.056*T, 7.5, Vulcan
VULCAN_ELEMENTS = VulcanElements(
    name="Vulcan",
    epoch=2415020.0,  # J1900.0
    a=0.13744,  # Semi-major axis in AU (inside Mercury's orbit)
    e=0.019,  # Small eccentricity (nearly circular)
    i=7.5,  # Inclination in degrees
    M0=252.8987988,  # Mean anomaly at epoch in degrees
    n_century=707550.7341,  # Mean motion in degrees per century
    omega0=322.212069,  # Argument of perihelion at epoch in degrees
    omega_rate=1670.056,  # Rate of change of omega in degrees/century
    Omega0=47.787931,  # Ascending node at epoch in degrees
    Omega_rate=-1670.056,  # Rate of change of Omega in degrees/century (note: negative)
)


# =============================================================================
# WALDEMATH BLACK MOON ORBITAL ELEMENTS (seorbel.txt #18)
# =============================================================================
# Waldemath Moon is a hypothetical second moon of Earth proposed by Dr. Georg
# Waldemath in 1898. This is distinct from Mean Lilith and True Lilith which
# are lunar apogee points.
#
# From Swiss Ephemeris seorbel.txt (section 2.7.7):
# J2000, JDATE, 248.8833, 0.0029833, 0.0, 0.0, 0.0, 0.0, Waldemath, geo # 18
#
# Elements order:
#   1. epoch (J2000)
#   2. equinox (JDATE - equinox of date)
#   3. mean longitude at epoch: 248.8833 degrees
#   4. semi-major axis: 0.0029833 AU (~446,200 km, ~1.16x Moon distance)
#   5. eccentricity: 0.0 (circular orbit)
#   6. argument of perihelion: 0.0
#   7. ascending node: 0.0
#   8. inclination: 0.0
#   9. name: Waldemath, geo (geocentric body)
#
# The orbital period is approximately 119 days based on the semi-major axis.
# Mean motion derived: n = 360 / 119 = ~3.025 degrees/day


@dataclass
class WaldemathElements:
    """
    Orbital elements for the Waldemath hypothetical second moon of Earth.

    Dr. Georg Waldemath claimed to have observed a second moon of Earth in 1898.
    This hypothetical body is sometimes called the "Dark Moon" or "Waldemath Moon".
    It is geocentric (orbits Earth, not the Sun).

    Note: This is different from:
    - Mean Lilith (mean lunar apogee)
    - True Lilith (osculating lunar apogee)
    - Black Moon Lilith (second focus of lunar orbit)

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT) - J2000.0
        a: Semi-major axis (AU) - geocentric distance from Earth
        e: Eccentricity (0 for circular orbit)
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Waldemath Moon orbital elements from Swiss Ephemeris seorbel.txt #18
# Epoch: J2000.0 (JD 2451545.0)
# From seorbel.txt line:
# J2000, JDATE, 248.8833, 0.0029833, 0.0, 0.0, 0.0, 0.0, Waldemath, geo # 18
#
# Mean motion calculation:
# For a geocentric orbit, period = 2*pi*sqrt(a³/mu) where mu = GM_Earth
# a = 0.0029833 AU = 446,200 km
# GM_Earth = 398600.4 km³/s²
# Period ≈ 119 days
# n = 360 / 119 ≈ 3.025 deg/day
WALDEMATH_ELEMENTS = WaldemathElements(
    name="Waldemath",
    epoch=2451545.0,  # J2000.0
    a=0.0029833,  # Semi-major axis in AU (geocentric distance)
    e=0.0,  # Circular orbit
    i=0.0,  # On ecliptic plane
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=248.8833,  # Mean longitude at J2000.0 (degrees)
    # Mean motion derived from orbital period of ~119 days
    # n = 360 / 119 ≈ 3.025 deg/day
    n=3.024873,  # degrees per day
)


# =============================================================================
# LOWELL PLANET X ORBITAL ELEMENTS (seorbel.txt #14)
# =============================================================================
# Percival Lowell's hypothetical "Planet X" that led to the search which discovered Pluto.
# Lowell calculated these orbital elements in 1915 based on perceived perturbations
# in Uranus's orbit. The search for this predicted planet led Clyde Tombaugh to
# discover Pluto in 1930 at Lowell Observatory.
#
# Historical note: Pluto turned out to be far too small (0.2% of Earth's mass) to
# cause the perturbations Lowell attributed to "Planet X". The perceived perturbations
# were later explained as observational errors. Pluto's discovery near Lowell's
# predicted position was essentially a coincidence.
#
# Lowell's 1915 prediction from "Memoir on a Trans-Neptunian Planet":
#   - Semi-major axis: 43.0 AU
#   - Eccentricity: 0.202
#   - Inclination: 10° to ecliptic
#   - Longitude of perihelion: ~204° (1930)
#   - Mean longitude at 1930.0: ~102°
#   - Orbital period: ~282 years


@dataclass
class LowellPlanetXElements:
    """
    Orbital elements for Lowell's hypothetical Planet X.

    Percival Lowell predicted this trans-Neptunian planet in 1915 based on
    analysis of supposed perturbations in Uranus's orbit. The search for this
    planet led to the discovery of Pluto in 1930, though Pluto was far too
    small to be Lowell's Planet X.

    Attributes:
        name: Name of the body
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


# Lowell Planet X orbital elements
# Based on Lowell's 1915 prediction "Memoir on a Trans-Neptunian Planet"
# Epoch: J1930.0 (JD 2425977.0) - chosen as it's close to Pluto's discovery date
#
# Elements derived from Lowell's published predictions:
#   - Semi-major axis: 43.0 AU
#   - Eccentricity: 0.202
#   - Inclination: 10.0 degrees
#   - Mean motion: 360 / (43^1.5 * 365.25) = 0.003497 deg/day
#   - Orbital period: ~282 years
LOWELL_PLANET_X_ELEMENTS = LowellPlanetXElements(
    name="Planet X Lowell",
    epoch=2415020.0,  # J1900.0 (matching other hypothetical bodies)
    a=43.0,  # Semi-major axis in AU
    e=0.202,  # Eccentricity from Lowell's prediction
    i=10.0,  # Inclination in degrees
    omega=204.9,  # Argument of perihelion in degrees
    Omega=67.5,  # Longitude of ascending node in degrees
    M0=11.2,  # Mean anomaly at epoch (degrees)
    # n = 360 / (a^1.5 * 365.25) = 360 / (281.98 * 365.25) = 0.003497 deg/day
    n=360.0 / (43.0**1.5 * 365.25),
)


# =============================================================================
# PICKERING PLANET X (PLANET O) ORBITAL ELEMENTS (seorbel.txt #15)
# =============================================================================
# William H. Pickering's hypothetical "Planet O" (1919), one of several trans-Neptunian
# planets he predicted (O, P, Q, R, S, T, U). Planet O was his most famous prediction.
#
# Historical note: Pickering, like Lowell, based his predictions on supposed
# perturbations in the orbits of outer planets. These perturbations were later
# found to be largely observational errors. Pickering's predictions, though detailed,
# did not lead to any discoveries.
#
# Planet O orbital elements from Pickering's 1919 paper:
#   - Semi-major axis: 51.9 AU (further out than Lowell's prediction)
#   - Eccentricity: 0.31
#   - Inclination: 15° to ecliptic
#   - Orbital period: ~373.5 years (derived from Kepler's 3rd law)
#
# Note: The elements below are reconstructed from historical sources. The epoch
# is set to J1900.0 to match other hypothetical bodies in seorbel.txt.


@dataclass
class PickeringPlanetXElements:
    """
    Orbital elements for Pickering's hypothetical Planet O/X.

    William H. Pickering predicted this trans-Neptunian planet in 1919 based on
    analysis of supposed perturbations in Uranus and Neptune orbits.

    Attributes:
        name: Name of the body
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


# Pickering Planet X (Planet O) orbital elements
# Based on Pickering's 1919 prediction
# Epoch: J1900.0 (JD 2415020.0) - matching other hypothetical bodies
#
# Elements from historical sources:
#   - Semi-major axis: 51.9 AU
#   - Eccentricity: 0.31
#   - Inclination: 15.0 degrees
#   - Mean motion: 360 / (51.9^1.5 * 365.25) = 0.00264 deg/day
#   - Orbital period: ~373.5 years
PICKERING_PLANET_X_ELEMENTS = PickeringPlanetXElements(
    name="Planet X Pickering",
    epoch=2415020.0,  # J1900.0 (matching other hypothetical bodies)
    a=51.9,  # Semi-major axis in AU
    e=0.31,  # Eccentricity from Pickering's prediction
    i=15.0,  # Inclination in degrees
    omega=100.0,  # Argument of perihelion in degrees (estimated)
    Omega=110.0,  # Longitude of ascending node in degrees (estimated)
    M0=15.0,  # Mean anomaly at epoch (degrees, estimated)
    # n = 360 / (a^1.5 * 365.25) = 360 / (373.87 * 365.25) = 0.002637 deg/day
    n=360.0 / (51.9**1.5 * 365.25),
)


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
    SE_PLUTO_LOWELL: "Planet X Lowell",
    SE_PLUTO_PICKERING: "Planet X Pickering",
    SE_VULCAN: "Vulcan",
    SE_WHITE_MOON: "White Moon (Selena)",
    SE_PROSERPINA: "Proserpina",
    SE_WALDEMATH: "Waldemath",
}


# =============================================================================
# SEORBEL.TXT PARSER
# =============================================================================
# Parser for Swiss Ephemeris seorbel.txt file format.
# This allows users to add custom hypothetical bodies by providing a
# seorbel.txt format file.


@dataclass
class TPolynomial:
    """
    Represents a polynomial expression in T (Julian centuries from epoch).

    The seorbel.txt format allows orbital elements to be expressed as
    polynomials in T, e.g., "252.8987988 + 707550.7341 * T".

    Attributes:
        constant: The constant term (coefficient of T^0)
        linear: The linear term (coefficient of T^1)
    """

    constant: float = 0.0
    linear: float = 0.0

    def evaluate(self, T: float) -> float:
        """Evaluate the polynomial at the given T value."""
        return self.constant + self.linear * T


@dataclass
class SeorbelElements:
    """
    Orbital elements for a fictitious body parsed from seorbel.txt format.

    The seorbel.txt format from Swiss Ephemeris defines orbital elements
    for fictitious bodies. Each line contains 9 comma-separated fields:
        1. epoch: Reference epoch (Julian day or "J1900", "B1950", "J2000")
        2. equinox: Coordinate equinox (Julian day, "J1900", "B1950", "J2000", or "JDATE")
        3. mean_anomaly: Mean anomaly at epoch (degrees, may include T-polynomial)
        4. semi_axis: Semi-major axis (AU)
        5. eccentricity: Orbital eccentricity (may include T-polynomial)
        6. arg_perihelion: Argument of perihelion (degrees, may include T-polynomial)
        7. asc_node: Longitude of ascending node (degrees, may include T-polynomial)
        8. inclination: Orbital inclination (degrees, may include T-polynomial)
        9. name: Body name (may include ", geo" suffix for geocentric bodies)

    T-polynomials allow elements to have secular variations, expressed as:
        "constant + rate * T"
    where T is Julian centuries from the epoch.

    Attributes:
        name: Name of the fictitious body
        epoch_jd: Reference epoch as Julian Day TT
        equinox_jd: Coordinate equinox as Julian Day TT (None if JDATE)
        equinox_is_jdate: True if equinox is "JDATE" (equinox of date)
        mean_anomaly: Mean anomaly polynomial (degrees)
        semi_axis: Semi-major axis (AU)
        eccentricity: Eccentricity polynomial
        arg_perihelion: Argument of perihelion polynomial (degrees)
        asc_node: Longitude of ascending node polynomial (degrees)
        inclination: Inclination polynomial (degrees)
        is_geocentric: True if body is geocentric (orbits Earth, not Sun)
        line_number: Original line number in the seorbel.txt file
    """

    name: str
    epoch_jd: float
    equinox_jd: Optional[float]
    equinox_is_jdate: bool
    mean_anomaly: TPolynomial
    semi_axis: float
    eccentricity: TPolynomial
    arg_perihelion: TPolynomial
    asc_node: TPolynomial
    inclination: TPolynomial
    is_geocentric: bool = False
    line_number: int = 0

    def get_mean_motion(self) -> float:
        """
        Calculate mean motion from semi-major axis using Kepler's 3rd law.

        Returns:
            Mean motion in degrees per day.
        """
        # n = 360 / (a^1.5 * 365.25) for heliocentric orbits
        # For geocentric orbits, this would need Earth's GM, but we use
        # an approximate formula.
        return 360.0 / (self.semi_axis**1.5 * 365.25)


# Standard epoch Julian Day values
_EPOCH_JD = {
    "J1900": 2415020.0,  # January 0.5, 1900 TDT
    "B1950": 2433282.42345905,  # Besselian year 1950.0
    "J2000": 2451545.0,  # January 1.5, 2000 TDT
}


def _parse_epoch_or_equinox(value: str) -> Tuple[Optional[float], bool]:
    """
    Parse an epoch or equinox value from seorbel.txt.

    Args:
        value: The epoch/equinox string (e.g., "J1900", "J2000", "JDATE", or a Julian day)

    Returns:
        Tuple of (julian_day, is_jdate) where julian_day is None if is_jdate is True.

    Raises:
        ValueError: If the value cannot be parsed.
    """
    value = value.strip().upper()

    if value == "JDATE":
        return (None, True)

    if value in _EPOCH_JD:
        return (_EPOCH_JD[value], False)

    # Try to parse as a Julian day number
    try:
        jd = float(value)
        return (jd, False)
    except ValueError:
        raise ValueError(f"Cannot parse epoch/equinox value: '{value}'")


def _parse_t_polynomial(expr: str) -> TPolynomial:
    """
    Parse a T-polynomial expression from seorbel.txt.

    Parses expressions like:
        - "252.8987988"
        - "252.8987988 + 707550.7341 * T"
        - "322.212069+1670.056*T"
        - "47.787931-1670.056*T"
        - "47.787931 - 1670.056 * T"

    Args:
        expr: The expression string

    Returns:
        TPolynomial with constant and linear coefficients.

    Raises:
        ValueError: If the expression cannot be parsed.
    """
    expr = expr.strip()

    # Check if this is a simple number (no T term)
    if "T" not in expr.upper():
        try:
            return TPolynomial(constant=float(expr), linear=0.0)
        except ValueError:
            raise ValueError(f"Cannot parse expression as number: '{expr}'")

    # Normalize the expression: ensure spaces around operators
    # But preserve signs as part of numbers
    expr_normalized = expr.replace("*", " * ").replace("+", " + ").replace("-", " - ")
    # Remove multiple spaces
    while "  " in expr_normalized:
        expr_normalized = expr_normalized.replace("  ", " ")
    expr_normalized = expr_normalized.strip()

    # Pattern to match expressions like "constant + rate * T" or "constant - rate * T"
    # Also handles "rate * T + constant" and variations with T first
    pattern = r"""
        ^
        \s*
        (?:
            # Pattern 1: constant [+/-] rate * T
            (?P<const1>-?[\d.]+)
            \s*
            (?P<op1>[+\-])
            \s*
            (?P<rate1>[\d.]+)
            \s*\*\s*
            T
        |
            # Pattern 2: rate * T [+/-] constant
            (?P<rate2>-?[\d.]+)
            \s*\*\s*
            T
            \s*
            (?P<op2>[+\-])
            \s*
            (?P<const2>[\d.]+)
        |
            # Pattern 3: just constant (no T) - already handled above
            (?P<const3>-?[\d.]+)
        )
        \s*
        $
    """

    match = re.match(pattern, expr, re.VERBOSE | re.IGNORECASE)

    if match:
        groups = match.groupdict()

        if groups.get("const1") is not None:
            # Pattern 1: constant +/- rate * T
            constant = float(groups["const1"])
            rate = float(groups["rate1"])
            if groups["op1"] == "-":
                rate = -rate
            return TPolynomial(constant=constant, linear=rate)

        elif groups.get("rate2") is not None:
            # Pattern 2: rate * T +/- constant
            rate = float(groups["rate2"])
            constant = float(groups["const2"])
            if groups["op2"] == "-":
                constant = -constant
            return TPolynomial(constant=constant, linear=rate)

        elif groups.get("const3") is not None:
            # Pattern 3: just constant
            return TPolynomial(constant=float(groups["const3"]), linear=0.0)

    # If the regex didn't match, try a simpler approach
    # Split on + or - while keeping the operator
    parts = re.split(r"(\s*[+\-]\s*)", expr)

    constant = 0.0
    linear = 0.0

    current_sign = 1.0
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if part in ["+", "+"]:
            current_sign = 1.0
        elif part in ["-", "-"]:
            current_sign = -1.0
        elif "T" in part.upper():
            # This is the T term, extract the coefficient
            t_part = part.upper().replace("T", "").replace("*", "").strip()
            if t_part == "" or t_part == "+":
                linear = current_sign * 1.0
            elif t_part == "-":
                linear = -1.0
            else:
                linear = current_sign * float(t_part)
        else:
            # This is the constant term
            constant = current_sign * float(part)

    return TPolynomial(constant=constant, linear=linear)


def parse_seorbel(filepath: Union[str, Path]) -> List[SeorbelElements]:
    """
    Parse a Swiss Ephemeris seorbel.txt file to extract orbital elements.

    The seorbel.txt file format defines orbital elements for fictitious/hypothetical
    bodies. This function parses the file and returns a list of SeorbelElements
    objects that can be used to compute positions of custom hypothetical bodies.

    File Format:
        - Lines starting with '#' (after optional whitespace) are comments
        - Empty or whitespace-only lines are ignored
        - Each data line has 9 comma-separated fields:
            1. epoch: Reference epoch (Julian day or "J1900", "B1950", "J2000")
            2. equinox: Coordinate equinox (Julian day, "J1900", "B1950", "J2000", or "JDATE")
            3. mean_anomaly: Mean anomaly at epoch (degrees, may include "+ rate * T")
            4. semi_axis: Semi-major axis (AU)
            5. eccentricity: Orbital eccentricity (may include T-polynomial)
            6. arg_perihelion: Argument of perihelion (degrees, may include T-polynomial)
            7. asc_node: Longitude of ascending node (degrees, may include T-polynomial)
            8. inclination: Orbital inclination (degrees, may include T-polynomial)
            9. name: Body name (may include ", geo" suffix for geocentric bodies)
        - Inline comments can appear after the 9th field, starting with '#'

    T-Polynomials:
        Some orbital elements can be expressed as polynomials in T, where
        T is Julian centuries from the epoch. For example:
            "252.8987988 + 707550.7341 * T"
        This represents an element that changes linearly with time.

    Geocentric Bodies:
        If the name field contains ", geo" (case-insensitive), the body is
        marked as geocentric (orbiting Earth rather than the Sun).

    Args:
        filepath: Path to the seorbel.txt file

    Returns:
        List of SeorbelElements objects, one for each valid data line.
        The list preserves the order of bodies as they appear in the file.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If a line cannot be parsed (with line number in message).

    Example:
        >>> from libephemeris.hypothetical import parse_seorbel
        >>> elements = parse_seorbel("seorbel.txt")
        >>> for elem in elements:
        ...     print(f"{elem.name}: a={elem.semi_axis} AU, e={elem.eccentricity.constant}")
        Cupido: a=40.99837 AU, e=0.0046
        Hades: a=50.66744 AU, e=0.00245
        ...

    Example with custom file:
        >>> # Create a custom seorbel file
        >>> with open("my_planet.txt", "w") as f:
        ...     f.write("# My custom hypothetical planet\\n")
        ...     f.write("J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, MyPlanet\\n")
        >>> elements = parse_seorbel("my_planet.txt")
        >>> print(elements[0].name)
        MyPlanet

    See Also:
        - Swiss Ephemeris documentation on fictitious objects
        - SeorbelElements dataclass for the structure of parsed elements
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"seorbel.txt file not found: {filepath}")

    elements: List[SeorbelElements] = []

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line_num, line in enumerate(f, start=1):
            # Remove trailing whitespace and carriage returns
            line = line.rstrip()

            # Skip empty lines
            if not line.strip():
                continue

            # Skip comment lines (lines starting with # after optional whitespace)
            if line.strip().startswith("#"):
                continue

            # Remove inline comments (everything after # that's not part of a field)
            # We need to be careful because # appears in inline comments after the 9th field
            # Strategy: split by comma first to get 9 fields, then handle comments
            try:
                parsed = _parse_seorbel_line(line, line_num)
                if parsed:
                    elements.append(parsed)
            except ValueError as e:
                raise ValueError(f"Error parsing line {line_num}: {e}") from e

    return elements


def get_bundled_seorbel_path() -> Path:
    """
    Get the path to the bundled seorbel.txt file included with libephemeris.

    The seorbel.txt file contains orbital elements for hypothetical bodies
    from the Swiss Ephemeris. This file is bundled with the libephemeris
    package for convenience.

    Returns:
        Path to the bundled seorbel.txt file.

    Raises:
        FileNotFoundError: If the bundled file is not found (should not happen
            in a properly installed package).

    Example:
        >>> from libephemeris.hypothetical import get_bundled_seorbel_path, parse_seorbel
        >>> seorbel_path = get_bundled_seorbel_path()
        >>> elements = parse_seorbel(seorbel_path)
        >>> print(f"Loaded {len(elements)} hypothetical bodies")
    """
    # Get the directory containing this module
    module_dir = Path(__file__).parent
    seorbel_path = module_dir / "seorbel.txt"

    if not seorbel_path.exists():
        raise FileNotFoundError(
            f"Bundled seorbel.txt not found at {seorbel_path}. "
            "This may indicate a packaging issue with libephemeris."
        )

    return seorbel_path


def load_bundled_seorbel() -> List[SeorbelElements]:
    """
    Load and parse the bundled seorbel.txt file included with libephemeris.

    This is a convenience function that combines get_bundled_seorbel_path()
    and parse_seorbel() to quickly load all hypothetical body elements from
    the Swiss Ephemeris seorbel.txt file bundled with the package.

    Returns:
        List of SeorbelElements objects for all hypothetical bodies defined
        in the bundled seorbel.txt file.

    Raises:
        FileNotFoundError: If the bundled file is not found.

    Example:
        >>> from libephemeris.hypothetical import load_bundled_seorbel, get_seorbel_body_by_name
        >>> elements = load_bundled_seorbel()
        >>> cupido = get_seorbel_body_by_name(elements, "Cupido")
        >>> print(f"Cupido semi-axis: {cupido.semi_axis} AU")
        Cupido semi-axis: 40.99837 AU

        >>> # Calculate position of a custom body from the file
        >>> from libephemeris.hypothetical import calc_seorbel_position
        >>> nibiru = get_seorbel_body_by_name(elements, "Nibiru")
        >>> if nibiru:
        ...     pos = calc_seorbel_position(nibiru, 2451545.0)
        ...     print(f"Nibiru longitude: {pos[0]:.4f} deg")

    See Also:
        - parse_seorbel: Parse a custom seorbel.txt file
        - get_bundled_seorbel_path: Get the path to the bundled file
        - get_seorbel_body_by_name: Find a body by name
        - calc_seorbel_position: Calculate position from elements
    """
    return parse_seorbel(get_bundled_seorbel_path())


def _parse_seorbel_line(line: str, line_num: int) -> Optional[SeorbelElements]:
    """
    Parse a single data line from seorbel.txt.

    Args:
        line: The line to parse
        line_num: Line number for error messages

    Returns:
        SeorbelElements object, or None if the line is a comment or empty.

    Raises:
        ValueError: If the line cannot be parsed.
    """
    # Handle the tricky part: the name field might contain commas,
    # and there might be a comment after the name.
    # Strategy: split carefully

    # First, try to find and remove trailing comment
    # Look for # that appears after what looks like the 9th field
    # The name field is the last, so we count 8 commas

    # Split by comma, but we need at least 9 parts
    parts = line.split(",")

    if len(parts) < 9:
        raise ValueError(
            f"Expected at least 9 comma-separated fields, got {len(parts)}"
        )

    # Parse the first 8 fields (they should be straightforward)
    epoch_str = parts[0].strip()
    equinox_str = parts[1].strip()
    mean_anomaly_str = parts[2].strip()
    semi_axis_str = parts[3].strip()
    eccentricity_str = parts[4].strip()
    arg_perihelion_str = parts[5].strip()
    asc_node_str = parts[6].strip()
    inclination_str = parts[7].strip()

    # The 9th field (name) might contain commas or be followed by a comment
    # Join remaining parts and handle
    name_and_rest = ",".join(parts[8:])

    # Remove inline comment if present
    # Look for # that's likely a comment (not part of the name)
    # The comment usually starts with " #" after the name
    if "#" in name_and_rest:
        # Find the first # that appears to be a comment
        # (usually preceded by whitespace or end of name)
        hash_idx = name_and_rest.find("#")
        name_field = name_and_rest[:hash_idx].strip()
    else:
        name_field = name_and_rest.strip()

    # Check for geocentric marker
    is_geocentric = False
    name = name_field
    if ", geo" in name_field.lower():
        is_geocentric = True
        # Remove the ", geo" suffix
        name = re.sub(r",\s*geo\s*$", "", name_field, flags=re.IGNORECASE).strip()
    elif " geo" in name_field.lower() and name_field.lower().endswith("geo"):
        # Handle "Waldemath, geo" format
        is_geocentric = True
        name = re.sub(r"\s*,?\s*geo\s*$", "", name_field, flags=re.IGNORECASE).strip()

    # Parse epoch and equinox
    epoch_jd, _ = _parse_epoch_or_equinox(epoch_str)
    if epoch_jd is None:
        raise ValueError(f"Epoch cannot be JDATE: '{epoch_str}'")

    equinox_jd, equinox_is_jdate = _parse_epoch_or_equinox(equinox_str)

    # Parse orbital elements (may be T-polynomials)
    mean_anomaly = _parse_t_polynomial(mean_anomaly_str)

    # Semi-axis is always a simple number
    try:
        semi_axis = float(semi_axis_str)
    except ValueError:
        raise ValueError(f"Cannot parse semi-major axis: '{semi_axis_str}'")

    eccentricity = _parse_t_polynomial(eccentricity_str)
    arg_perihelion = _parse_t_polynomial(arg_perihelion_str)
    asc_node = _parse_t_polynomial(asc_node_str)
    inclination = _parse_t_polynomial(inclination_str)

    return SeorbelElements(
        name=name,
        epoch_jd=epoch_jd,
        equinox_jd=equinox_jd,
        equinox_is_jdate=equinox_is_jdate,
        mean_anomaly=mean_anomaly,
        semi_axis=semi_axis,
        eccentricity=eccentricity,
        arg_perihelion=arg_perihelion,
        asc_node=asc_node,
        inclination=inclination,
        is_geocentric=is_geocentric,
        line_number=line_num,
    )


def get_seorbel_body_by_name(
    elements: List[SeorbelElements], name: str
) -> Optional[SeorbelElements]:
    """
    Find a body in a parsed seorbel elements list by name.

    Args:
        elements: List of SeorbelElements from parse_seorbel()
        name: Name to search for (case-insensitive)

    Returns:
        The SeorbelElements object if found, None otherwise.

    Example:
        >>> elements = parse_seorbel("seorbel.txt")
        >>> cupido = get_seorbel_body_by_name(elements, "Cupido")
        >>> if cupido:
        ...     print(f"Cupido semi-axis: {cupido.semi_axis} AU")
    """
    name_lower = name.lower()
    for elem in elements:
        if elem.name.lower() == name_lower:
            return elem
    return None


def calc_seorbel_position(
    elem: SeorbelElements, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of a body from parsed seorbel.txt elements.

    This function uses the orbital elements from a SeorbelElements object
    to compute the heliocentric (or geocentric for ", geo" bodies) position
    using Keplerian propagation.

    For bodies with time-dependent elements (T-polynomials), the elements
    are evaluated at the given Julian date.

    Args:
        elem: SeorbelElements object from parse_seorbel()
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance in AU (from Sun or Earth for geocentric bodies)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> elements = parse_seorbel("seorbel.txt")
        >>> cupido = get_seorbel_body_by_name(elements, "Cupido")
        >>> if cupido:
        ...     pos = calc_seorbel_position(cupido, 2451545.0)  # J2000.0
        ...     print(f"Cupido at {pos[0]:.4f} deg")
    """
    # Time in Julian centuries from epoch
    T = (jd_tt - elem.epoch_jd) / 36525.0
    # Days from epoch (for Keplerian motion)
    dt_days = jd_tt - elem.epoch_jd

    # Evaluate time-dependent elements
    # For mean anomaly, we have two components:
    # 1. The polynomial value (constant + linear*T where T is in centuries)
    # 2. Keplerian mean motion from semi-major axis (if not already in polynomial)
    #
    # In seorbel.txt, if the mean_anomaly has a large linear T-term, it represents
    # the total mean motion (deg/century). If it doesn't have a T-term or has a
    # small one (secular perturbations), we compute motion from Kepler's 3rd law.
    M_poly = elem.mean_anomaly.evaluate(T)

    # Threshold: typical mean motion for 100 AU body is ~3600 deg/century
    # Any explicit motion > 10 deg/century likely includes Keplerian motion
    if abs(elem.mean_anomaly.linear) < 10.0:
        # Add Keplerian mean motion: n = 360 / (a^1.5 * 365.25) deg/day
        n_deg_per_day = elem.get_mean_motion()
        M_deg = M_poly + n_deg_per_day * dt_days
    else:
        M_deg = M_poly

    e = elem.eccentricity.evaluate(T)
    omega_deg = elem.arg_perihelion.evaluate(T)
    Omega_deg = elem.asc_node.evaluate(T)
    i_deg = elem.inclination.evaluate(T)
    a = elem.semi_axis

    # Normalize mean anomaly
    M_deg = M_deg % 360.0
    M_rad = math.radians(M_deg)

    # For circular orbits (e=0), mean anomaly = true anomaly
    # We still need to do proper coordinate conversion for inclined orbits
    if e < 1e-10:
        # Circular orbit: E = M, nu = M
        nu = M_rad
        r = a  # Distance is constant
    else:
        # Solve Kepler's equation for eccentric anomaly
        E = _solve_kepler_equation(M_rad, e)

        # True anomaly
        sqrt_term = math.sqrt((1.0 + e) / (1.0 - e))
        nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

        # Distance
        r = a * (1.0 - e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    omega_rad = math.radians(omega_deg)
    u = nu + omega_rad

    # Convert to ecliptic coordinates
    i_rad = math.radians(i_deg)
    Omega_rad = math.radians(Omega_deg)

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
    dt_step = 1.0  # 1 day step for daily velocity
    pos_next = _calc_seorbel_position_raw(elem, jd_tt + dt_step)

    dlon = pos_next[0] - longitude
    # Handle wrap-around
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0

    dlat = pos_next[1] - latitude
    ddist = pos_next[2] - distance

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_seorbel_position_raw(
    elem: SeorbelElements, jd_tt: float
) -> Tuple[float, float, float]:
    """
    Calculate raw position without velocity (helper for differentiation).
    """
    # Time in Julian centuries from epoch
    T = (jd_tt - elem.epoch_jd) / 36525.0
    # Days from epoch (for Keplerian motion)
    dt_days = jd_tt - elem.epoch_jd

    # Evaluate time-dependent elements
    # Same logic as calc_seorbel_position for mean anomaly
    M_poly = elem.mean_anomaly.evaluate(T)
    if abs(elem.mean_anomaly.linear) < 10.0:
        n_deg_per_day = elem.get_mean_motion()
        M_deg = M_poly + n_deg_per_day * dt_days
    else:
        M_deg = M_poly

    e = elem.eccentricity.evaluate(T)
    omega_deg = elem.arg_perihelion.evaluate(T)
    Omega_deg = elem.asc_node.evaluate(T)
    i_deg = elem.inclination.evaluate(T)
    a = elem.semi_axis

    # Normalize mean anomaly
    M_deg = M_deg % 360.0
    M_rad = math.radians(M_deg)

    # For circular orbits (e=0), mean anomaly = true anomaly
    if e < 1e-10:
        nu = M_rad
        r = a
    else:
        # Solve Kepler's equation
        E = _solve_kepler_equation(M_rad, e)

        # True anomaly
        sqrt_term = math.sqrt((1.0 + e) / (1.0 - e))
        nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

        # Distance
        r = a * (1.0 - e * math.cos(E))

    # Argument of latitude
    omega_rad = math.radians(omega_deg)
    u = nu + omega_rad

    # Convert to ecliptic coordinates
    i_rad = math.radians(i_deg)
    Omega_rad = math.radians(Omega_deg)

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

    return (longitude, latitude, r)


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


def calc_cupido(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Cupido using Keplerian propagation.

    Cupido is the first Hamburg School Uranian planet. This function uses
    orbital elements calibrated against pyswisseph for exact compatibility.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Cupido)
            - distance: Distance from Sun in AU (40.99837 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_cupido
        >>> pos = calc_cupido(2451545.0)  # J2000.0
        >>> print(f"Cupido at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(SE_CUPIDO, jd_tt)


def calc_hades(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Hades using Keplerian propagation.

    Hades is the second Hamburg School Uranian planet. This function uses
    orbital elements calibrated against pyswisseph for exact compatibility.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_hades
        >>> pos = calc_hades(2451545.0)  # J2000.0
        >>> print(f"Hades at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(SE_HADES, jd_tt)


def calc_zeus(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Zeus using Keplerian propagation.

    Zeus is the third Hamburg School Uranian planet. This function uses
    orbital elements calibrated against pyswisseph for exact compatibility.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Zeus)
            - distance: Distance from Sun in AU (59.21436 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_zeus
        >>> pos = calc_zeus(2451545.0)  # J2000.0
        >>> print(f"Zeus at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(SE_ZEUS, jd_tt)


def calc_kronos(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Kronos using Keplerian propagation.

    Kronos is the fourth Hamburg School Uranian planet. This function uses
    orbital elements calibrated against pyswisseph for exact compatibility.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Kronos)
            - distance: Distance from Sun in AU (64.81690 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_kronos
        >>> pos = calc_kronos(2451545.0)  # J2000.0
        >>> print(f"Kronos at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(SE_KRONOS, jd_tt)


def calc_apollon(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Apollon using Keplerian propagation.

    Apollon is the fifth Hamburg School Uranian planet. This function uses
    orbital elements calibrated against pyswisseph for exact compatibility.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Apollon)
            - distance: Distance from Sun in AU (70.361180 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_apollon
        >>> pos = calc_apollon(2451545.0)  # J2000.0
        >>> print(f"Apollon at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(SE_APOLLON, jd_tt)


def calc_admetos(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Admetos using Keplerian propagation.

    Admetos is the sixth Hamburg School Uranian planet. This function uses
    orbital elements calibrated against pyswisseph for exact compatibility.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Admetos)
            - distance: Distance from Sun in AU (73.736396 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_admetos
        >>> pos = calc_admetos(2451545.0)  # J2000.0
        >>> print(f"Admetos at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(SE_ADMETOS, jd_tt)


def calc_vulkanus(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Vulkanus using Keplerian propagation.

    Vulkanus is the seventh Hamburg School Uranian planet. This function uses
    orbital elements calibrated against pyswisseph for exact compatibility.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Vulkanus)
            - distance: Distance from Sun in AU (77.445895 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_vulkanus
        >>> pos = calc_vulkanus(2451545.0)  # J2000.0
        >>> print(f"Vulkanus at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(SE_VULKANUS, jd_tt)


def calc_poseidon(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Poseidon using Keplerian propagation.

    Poseidon is the eighth Hamburg School Uranian planet. This function uses
    orbital elements calibrated against pyswisseph for exact compatibility.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Poseidon)
            - distance: Distance from Sun in AU (83.666307 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_poseidon
        >>> pos = calc_poseidon(2451545.0)  # J2000.0
        >>> print(f"Poseidon at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(SE_POSEIDON, jd_tt)


def calc_transpluto(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Transpluto (Isis) using Keplerian propagation.

    Transpluto is a hypothetical trans-Plutonian planet proposed by astrologer Ram,
    documented in Swiss Ephemeris section 2.7.2 (seorbel.txt). Unlike the Hamburg
    School Uranian planets which have nearly circular orbits, Transpluto has
    significant eccentricity (e=0.3).

    Orbital elements from seorbel.txt:
        - Epoch: 2368547.66 (1772.76)
        - Semi-major axis: 77.775 AU
        - Eccentricity: 0.3
        - Inclination: 0 degrees
        - Argument of perihelion: 0.7 degrees
        - Ascending node: 0 degrees
        - Mean anomaly at epoch: 0 degrees (at perihelion)

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (0 for Transpluto)
            - distance: Distance from Sun in AU
            - dlon: Daily heliocentric longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (0)
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_transpluto
        >>> pos = calc_transpluto(2451545.0)  # J2000.0
        >>> print(f"Transpluto at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    elements = TRANSPLUTO_KEPLERIAN_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

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
    dt_step = 1.0  # 1 day step for daily velocity
    pos_next = _calc_transpluto_raw(jd_tt + dt_step)

    dlon = pos_next[0] - longitude
    # Handle wrap-around
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0

    dlat = pos_next[1] - latitude
    ddist = pos_next[2] - distance

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_transpluto_raw(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Transpluto position without velocity (helper for differentiation).
    """
    elements = TRANSPLUTO_KEPLERIAN_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

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


def calc_uranian_planet(
    body_id: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of any Uranian planet using Keplerian propagation.

    This generic function handles all eight Hamburg School Uranian planets
    (Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon) by
    looking up their orbital elements from the URANIAN_KEPLERIAN_ELEMENTS
    dictionary and performing Keplerian propagation.

    For circular orbits (e=0), the calculation simplifies to mean longitude
    propagation. For elliptic orbits (e>0, like Hades), full Keplerian
    mechanics with Kepler's equation solving is used.

    Args:
        body_id: Uranian planet ID (SE_CUPIDO through SE_POSEIDON, i.e., 40-47)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Raises:
        ValueError: If body_id is not a valid Uranian planet ID.

    Example:
        >>> from libephemeris.hypothetical import calc_uranian_planet, SE_CUPIDO
        >>> pos = calc_uranian_planet(SE_CUPIDO, 2451545.0)  # J2000.0
        >>> print(f"Cupido at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    if body_id not in URANIAN_KEPLERIAN_ELEMENTS:
        raise ValueError(
            f"Body ID {body_id} is not a valid Uranian planet. "
            f"Valid IDs: {list(URANIAN_KEPLERIAN_ELEMENTS.keys())}"
        )

    elements = URANIAN_KEPLERIAN_ELEMENTS[body_id]

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Check if this is a circular orbit (e=0)
    if elements.e == 0.0:
        # For circular orbit, simply propagate mean longitude
        longitude = (elements.M0 + elements.n * dt) % 360.0
        latitude = 0.0
        distance = elements.a

        # For circular orbit, velocities are constant
        dlon = elements.n
        dlat = 0.0
        ddist = 0.0

        return (longitude, latitude, distance, dlon, dlat, ddist)

    # Elliptic orbit: full Keplerian propagation
    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

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
    dt_step = 1.0  # 1 day step for daily velocity
    pos_next = _calc_uranian_planet_raw(body_id, jd_tt + dt_step)

    dlon = pos_next[0] - longitude
    # Handle wrap-around
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0

    dlat = pos_next[1] - latitude
    ddist = pos_next[2] - distance

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_uranian_planet_raw(body_id: int, jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Uranian planet position without velocity (helper for differentiation).

    This internal helper function calculates only position (no velocity) for use
    in numerical differentiation to compute velocities.

    Args:
        body_id: Uranian planet ID
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance)
    """
    elements = URANIAN_KEPLERIAN_ELEMENTS[body_id]

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Check if this is a circular orbit (e=0)
    if elements.e == 0.0:
        longitude = (elements.M0 + elements.n * dt) % 360.0
        latitude = 0.0
        distance = elements.a
        return (longitude, latitude, distance)

    # Elliptic orbit
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

    E = _solve_kepler_equation(M_rad, elements.e)

    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    r = elements.a * (1.0 - elements.e * math.cos(E))
    u = nu + math.radians(elements.omega)

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


def calc_vulcan(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Vulcan using Swiss Ephemeris elements.

    Vulcan is a hypothetical intramercurial planet documented in Swiss Ephemeris
    section 2.7.5 (seorbel.txt #16). Various astrologers proposed different orbital
    elements for Vulcan; this implementation uses the version from Swiss Ephemeris.

    Unlike other hypothetical bodies, Vulcan has time-dependent orbital elements:
        - Mean anomaly: 252.8987988 + 707550.7341 * T degrees
        - Argument of perihelion: 322.212069 + 1670.056 * T degrees
        - Ascending node: 47.787931 - 1670.056 * T degrees

    where T = Julian centuries from J1900.0 (JD 2415020.0)

    Orbital elements from seorbel.txt:
        - Epoch: J1900.0 (JD 2415020.0)
        - Equinox: JDATE (equinox of date)
        - Semi-major axis: 0.13744 AU (inside Mercury's orbit)
        - Eccentricity: 0.019 (nearly circular)
        - Inclination: 7.5 degrees
        - Orbital period: ~18.6 days (derived from mean motion)

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily heliocentric longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_vulcan
        >>> pos = calc_vulcan(2451545.0)  # J2000.0
        >>> print(f"Vulcan at {pos[0]:.4f} deg, distance {pos[2]:.4f} AU")
    """
    pos = _calc_vulcan_raw(jd_tt)
    longitude, latitude, distance = pos

    # Calculate velocity via numerical differentiation
    dt_step = 1.0  # 1 day step for daily velocity
    pos_next = _calc_vulcan_raw(jd_tt + dt_step)

    dlon = pos_next[0] - longitude
    # Handle wrap-around (Vulcan moves ~19 deg/day, so this is important)
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0

    dlat = pos_next[1] - latitude
    ddist = pos_next[2] - distance

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_vulcan_raw(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Vulcan position without velocity (helper for differentiation).

    Implements the time-dependent orbital elements from Swiss Ephemeris seorbel.txt.
    """
    elements = VULCAN_ELEMENTS

    # Time in Julian centuries from epoch (J1900.0)
    T = (jd_tt - elements.epoch) / 36525.0

    # Compute time-dependent elements
    # Mean anomaly: M0 + n_century * T
    M = (elements.M0 + elements.n_century * T) % 360.0
    M_rad = math.radians(M)

    # Argument of perihelion: omega0 + omega_rate * T
    omega = elements.omega0 + elements.omega_rate * T
    omega_rad = math.radians(omega)

    # Ascending node: Omega0 + Omega_rate * T
    Omega = elements.Omega0 + elements.Omega_rate * T
    Omega_rad = math.radians(Omega)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    u = nu + omega_rad

    # Convert to ecliptic coordinates
    i_rad = math.radians(elements.i)

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

    return (longitude, latitude, distance)


def calc_waldemath(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the geocentric position of the Waldemath hypothetical second moon.

    Dr. Georg Waldemath's hypothetical second moon of Earth, documented in Swiss
    Ephemeris section 2.7.7 (seorbel.txt #18). This is a different body from:
    - Mean Lilith (SE_MEAN_APOG) - the mean lunar apogee
    - True Lilith (SE_OSCU_APOG) - the osculating lunar apogee
    - The second focus of the lunar orbit

    Waldemath claimed to have observed this body in 1898, describing it as a dark
    moon with an orbital period of approximately 119 days.

    Orbital elements from seorbel.txt:
        - Epoch: J2000.0 (JD 2451545.0)
        - Semi-major axis: 0.0029833 AU (~446,200 km, ~1.16x Moon distance)
        - Eccentricity: 0.0 (circular orbit)
        - Inclination: 0.0 degrees (on ecliptic)
        - Mean longitude at epoch: 248.8833 degrees
        - Orbital period: ~119 days

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Geocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (0 for circular ecliptic orbit)
            - distance: Distance from Earth in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (0)
            - ddist: Daily distance change in AU/day (0 for circular orbit)

    Example:
        >>> from libephemeris.hypothetical import calc_waldemath
        >>> pos = calc_waldemath(2451545.0)  # J2000.0
        >>> print(f"Waldemath at {pos[0]:.4f} deg, distance {pos[2]:.6f} AU")
    """
    elements = WALDEMATH_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # For a circular orbit (e = 0), mean longitude = true longitude
    # Simply propagate the mean longitude
    longitude = (elements.L0 + elements.n * dt) % 360.0

    # Waldemath is assumed to be on the ecliptic (zero inclination)
    latitude = 0.0

    # Distance is constant for circular orbit (equal to semi-major axis)
    distance = elements.a

    # Daily motion is simply the mean motion for circular orbit
    dlon = elements.n

    # No latitude or distance change for circular orbit on ecliptic
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
    use_true_lilith: bool = False,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of the White Moon (Selena).

    The White Moon (also called Selena) is a symbolic astrological point defined
    as the lunar perigee - the point 180 degrees opposite to Black Moon Lilith
    (the lunar apogee). It represents the closest approach of the Moon to Earth,
    symbolically associated with positive lunar qualities.

    Swiss Ephemeris Definition:
        White Moon Selena = Mean Lilith + 180° (default, matching Swiss Ephemeris)
        This uses the mean lunar apogee, which ignores short-period oscillations.

    True Lilith-based Definition:
        White Moon Selena = True Lilith + 180° (optional via use_true_lilith=True)
        This uses the osculating (true) lunar apogee, which includes perturbations.
        The true apogee oscillates ±5-10° from the mean position.

    Astronomical Background:
        - Black Moon Lilith (lunar apogee) = point where Moon is farthest from Earth
        - White Moon Selena (lunar perigee) = point where Moon is closest to Earth
        - Mean apogee progresses at ~40.69°/year (apsidal precession period ~8.85 years)

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)
        use_true_lilith: If True, calculate based on True (osculating) Lilith
                         instead of Mean Lilith. Default is False to match
                         the Swiss Ephemeris convention.

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (0 for mean-based, varies for true)
            - distance: Distance (0 for symbolic point, or apogee distance for true)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Note:
        This is a symbolic calculation. The White Moon is not a physical body,
        but rather a calculated point representing the lunar perigee direction.

    Example:
        >>> from libephemeris.hypothetical import calc_white_moon_position
        >>> pos = calc_white_moon_position(2451545.0)  # J2000.0, mean-based
        >>> print(f"White Moon at {pos[0]:.4f} deg")
        >>> pos_true = calc_white_moon_position(2451545.0, use_true_lilith=True)
        >>> print(f"White Moon (true) at {pos_true[0]:.4f} deg")

    References:
        - Swiss Ephemeris documentation on fictitious objects
        - Jacobson: "The Dark Moon Lilith in Astrology" (1961)
    """
    # Import lunar module to get Lilith position
    from . import lunar

    if use_true_lilith:
        # Use True (osculating) Lilith - includes perturbations
        lilith_lon, lilith_lat, lilith_dist = lunar.calc_true_lilith(jd_tt)

        # White Moon is opposite to Black Moon
        longitude = (lilith_lon + 180.0) % 360.0

        # For true-based calculation, latitude is opposite (symmetric point)
        # Note: True Lilith has non-zero latitude due to inclination of apsidal line
        latitude = -lilith_lat  # Opposite point has opposite latitude

        # Distance is same as Lilith (apogee distance = perigee distance conceptually,
        # but for a symbolic point we use 0)
        distance = 0.0  # Symbolic point, no meaningful distance

        # Calculate velocity via numerical differentiation
        dt = 1.0 / 86400.0  # 1 second step
        lilith_next_lon, lilith_next_lat, _ = lunar.calc_true_lilith(jd_tt + dt)

        lilith_dlon = (lilith_next_lon - lilith_lon) / dt
        # Handle wrap-around
        if lilith_dlon > 180.0 / dt:
            lilith_dlon -= 360.0 / dt
        elif lilith_dlon < -180.0 / dt:
            lilith_dlon += 360.0 / dt

        dlon = lilith_dlon  # Same rate, just offset by 180 degrees
        dlat = -(lilith_next_lat - lilith_lat) / dt  # Opposite latitude change
        ddist = 0.0
    else:
        # Use Mean Lilith - default, matching Swiss Ephemeris convention
        lilith_lon = lunar.calc_mean_lilith(jd_tt)

        # White Moon is opposite to Black Moon
        longitude = (lilith_lon + 180.0) % 360.0

        # Latitude, distance, and velocities are simplified for mean calculation
        latitude = 0.0
        distance = 0.0  # Symbolic point, no meaningful distance

        # Calculate velocity via numerical differentiation
        dt = 1.0 / 86400.0  # 1 second step
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
    Calculate the position of the Waldemath Moon (hypothetical second moon of Earth).

    Dr. Georg Waldemath's hypothetical second moon of Earth, documented in Swiss
    Ephemeris section 2.7.7 (seorbel.txt #18). This is a different body from:
    - Mean Lilith (SE_MEAN_APOG) - the mean lunar apogee
    - True Lilith (SE_OSCU_APOG) - the osculating lunar apogee
    - The second focus of the lunar orbit

    Waldemath claimed to have observed this body in 1898, describing it as a dark
    moon with an orbital period of approximately 119 days.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Geocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Earth in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Note:
        This is a hypothetical body. Waldemath's observations were never confirmed
        and no such second moon has been found to exist.

    Example:
        >>> from libephemeris.hypothetical import calc_waldemath_position
        >>> pos = calc_waldemath_position(2451545.0)  # J2000.0
        >>> print(f"Waldemath at {pos[0]:.4f} deg")
    """
    return calc_waldemath(jd_tt)


def calc_proserpina(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Proserpina using Keplerian propagation.

    Proserpina is a hypothetical trans-Plutonian planet used by some astrologers.
    This is NOT the same as the asteroid 26 Proserpina. Unlike other hypothetical
    bodies documented in Swiss Ephemeris seorbel.txt, Proserpina is not part of
    the standard Swiss Ephemeris fictitious bodies.

    The name "Proserpina" refers to the Roman goddess of the underworld (Greek:
    Persephone), wife of Pluto. In astrological usage, it represents transformation,
    cycles of death and rebirth, and the shadow self.

    Orbital elements used (traditional astrological):
        - Epoch: J2000.0 (JD 2451545.0)
        - Semi-major axis: 81.0 AU (beyond Neptune and Pluto)
        - Eccentricity: 0.0 (circular orbit)
        - Inclination: 0.0 degrees (on ecliptic)
        - Orbital period: ~729 years (derived from Kepler's 3rd law)

    Note: Different astrologers may use different orbital elements for Proserpina.
    This implementation uses a simple circular orbit model.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (0 for circular orbit on ecliptic)
            - distance: Distance from Sun in AU (81.0 AU, constant for circular orbit)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (0)
            - ddist: Daily distance change in AU/day (0 for circular orbit)

    Example:
        >>> from libephemeris.hypothetical import calc_proserpina
        >>> pos = calc_proserpina(2451545.0)  # J2000.0
        >>> print(f"Proserpina at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    elements = HYPOTHETICAL_ELEMENTS[SE_PROSERPINA]

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # For a circular orbit (e = 0), mean longitude = true longitude
    # Simply propagate the mean longitude
    longitude = (elements.M0 + elements.n * dt) % 360.0

    # Proserpina is assumed to be on the ecliptic (zero inclination)
    latitude = 0.0

    # Distance is constant for circular orbit (equal to semi-major axis)
    distance = elements.a

    # Daily motion is simply the mean motion for circular orbit
    dlon = elements.n

    # No latitude or distance change for circular orbit on ecliptic
    dlat = 0.0
    ddist = 0.0

    return (longitude, latitude, distance, dlon, dlat, ddist)


def calc_planet_x_lowell(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Lowell's Planet X using Keplerian propagation.

    Percival Lowell's hypothetical "Planet X" was a trans-Neptunian planet he predicted
    in 1915 based on perceived perturbations in Uranus's orbit. The search for this
    planet at Lowell Observatory led to Clyde Tombaugh's discovery of Pluto in 1930.

    Historical note: Pluto was far too small (0.2% of Earth's mass) to cause the
    perturbations Lowell attributed to "Planet X". The perceived perturbations were
    later explained as observational errors. Pluto's discovery near Lowell's predicted
    position was essentially a fortunate coincidence.

    Orbital elements from Lowell's 1915 prediction "Memoir on a Trans-Neptunian Planet":
        - Semi-major axis: 43.0 AU
        - Eccentricity: 0.202
        - Inclination: 10 degrees to ecliptic
        - Orbital period: ~282 years

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_planet_x_lowell
        >>> pos = calc_planet_x_lowell(2451545.0)  # J2000.0
        >>> print(f"Planet X Lowell at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    elements = LOWELL_PLANET_X_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

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
    dt_step = 1.0  # 1 day step for daily velocity
    pos_next = _calc_planet_x_lowell_raw(jd_tt + dt_step)

    dlon = pos_next[0] - longitude
    # Handle wrap-around
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0

    dlat = pos_next[1] - latitude
    ddist = pos_next[2] - distance

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_planet_x_lowell_raw(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Planet X Lowell position without velocity (helper for differentiation).
    """
    elements = LOWELL_PLANET_X_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

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


def calc_planet_x_pickering(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Pickering's Planet X using Keplerian propagation.

    William H. Pickering's hypothetical "Planet O" was a trans-Neptunian planet he predicted
    in 1919 based on supposed perturbations in Uranus and Neptune orbits. Pickering proposed
    several hypothetical planets (O, P, Q, R, S, T, U), with Planet O being the most famous.

    Historical note: Like Lowell's Planet X, Pickering's predictions were based on supposed
    perturbations that later proved to be observational errors. No planet was ever found
    at Pickering's predicted positions.

    Orbital elements from Pickering's 1919 prediction:
        - Semi-major axis: 51.9 AU
        - Eccentricity: 0.31
        - Inclination: 15 degrees to ecliptic
        - Orbital period: ~373.5 years (derived from Kepler's 3rd law)

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_planet_x_pickering
        >>> pos = calc_planet_x_pickering(2451545.0)  # J2000.0
        >>> print(f"Planet X Pickering at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    elements = PICKERING_PLANET_X_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

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
    dt_step = 1.0  # 1 day step for daily velocity
    pos_next = _calc_planet_x_pickering_raw(jd_tt + dt_step)

    dlon = pos_next[0] - longitude
    # Handle wrap-around
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0

    dlat = pos_next[1] - latitude
    ddist = pos_next[2] - distance

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_planet_x_pickering_raw(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Planet X Pickering position without velocity (helper for differentiation).
    """
    elements = PICKERING_PLANET_X_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

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

    # Vulcan (intramercurial hypothetical planet)
    if ipl == SE_VULCAN:
        return calc_vulcan(jd_tt)

    # White Moon (Selena)
    if ipl == SE_WHITE_MOON:
        return calc_white_moon_position(jd_tt)

    # Waldemath Black Moon
    if ipl == SE_WALDEMATH:
        return calc_waldemath_position(jd_tt)

    # Proserpina (hypothetical trans-Plutonian planet)
    if ipl == SE_PROSERPINA:
        return calc_proserpina(jd_tt)

    # Planet X Lowell (Lowell's predicted trans-Neptunian planet)
    if ipl == SE_PLANET_X_LOWELL:
        return calc_planet_x_lowell(jd_tt)

    # Planet X Pickering (Pickering's predicted trans-Neptunian planet)
    if ipl == SE_PLANET_X_PICKERING:
        return calc_planet_x_pickering(jd_tt)

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
