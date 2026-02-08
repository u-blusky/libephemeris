"""
Moshier semi-analytical ephemeris package.

This package provides semi-analytical algorithms for calculating celestial
body positions WITHOUT any dependency on Skyfield or SPK files. It implements:

- VSOP87: Planetary theories for Mercury through Neptune
- ELP 2000-82B: Lunar theory
- Pluto theory: Chapront-Francou periodic terms and Keplerian fallback
- IAU 2006 precession and nutation

The calculations are activated when the SEFLG_MOSEPH flag (4) is passed
to swe_calc_ut() or swe_calc(). The range is approximately -3000 CE to +3000 CE,
much wider than the JPL DE440 ephemeris range (1550-2650 CE).

Usage:
    from libephemeris.moshier import calc_position

    # Calculate Sun position
    lon, lat, dist, dlon, dlat, ddist = calc_position(jd_tt, SE_SUN)

    # Check if a body is supported
    from libephemeris.moshier import is_moshier_body
    if is_moshier_body(body_id):
        position = calc_position(jd_tt, body_id)

Accuracy:
    - Planets (Mercury-Neptune): ~1 arcsecond (truncated VSOP87)
    - Moon: ~1 arcsecond (ELP 2000-82B truncated)
    - Pluto (1885-2099): ~0.1 degrees
    - Pluto (outside range): ~1 degree

References:
    - Bretagnon, P. and Francou, G. (1988), "Planetary theories in rectangular
      and spherical variables. VSOP87 solutions", A&A 202, 309-315
    - Chapront-Touzé, M. and Chapront, J. (1983), "The lunar ephemeris ELP 2000",
      A&A 124, 50-62
    - Meeus, J. (1991), "Astronomical Algorithms"
    - Capitaine, N. et al. (2003), "Expressions for IAU 2000 precession quantities"
"""

from __future__ import annotations

from typing import Tuple

# Import body ID constants
from .elp82b import MOSHIER_MOON
from .pluto import MOSHIER_PLUTO
from .vsop87 import (
    MOSHIER_EARTH,
    MOSHIER_JUPITER,
    MOSHIER_MARS,
    MOSHIER_MERCURY,
    MOSHIER_NEPTUNE,
    MOSHIER_SATURN,
    MOSHIER_SUN,
    MOSHIER_URANUS,
    MOSHIER_VENUS,
)

# Import calculation functions from submodules
from .vsop87 import calc_position as vsop_calc_position
from .vsop87 import is_vsop_body
from .elp82b import calc_position as elp_calc_position
from .elp82b import is_moon_body
from .pluto import calc_position as pluto_calc_position
from .pluto import is_pluto_body

# Import precession/nutation functions
from .precession import (
    mean_obliquity,
    true_obliquity,
    nutation_angles,
    precession_angles,
    precess_ecliptic,
    precess_to_j2000,
    precess_from_j2000,
)

# Import utility functions
from .utils import (
    # Constants
    J2000,
    JD_PER_CENTURY,
    AU_KM,
    C_LIGHT_AU_DAY,
    DEG_TO_RAD,
    RAD_TO_DEG,
    ARCSEC_TO_RAD,
    RAD_TO_ARCSEC,
    # Angle normalization
    normalize_angle,
    normalize_radians,
    # Coordinate conversions
    spherical_to_cartesian,
    cartesian_to_spherical,
    ecliptic_to_equatorial,
    equatorial_to_ecliptic,
    rectangular_to_spherical_velocity,
    # Aberration
    annual_aberration,
    light_time_correction,
    # Time conversions
    jd_to_julian_centuries,
    jd_to_julian_millennia,
    julian_centuries_to_jd,
    # Numerical utilities
    numerical_derivative,
    evaluate_polynomial,
    evaluate_polynomial_derivative,
)

# =============================================================================
# SUPPORTED BODIES
# =============================================================================

# All body IDs supported by the Moshier ephemeris
MOSHIER_BODIES = frozenset(
    [
        MOSHIER_SUN,
        MOSHIER_MOON,
        MOSHIER_MERCURY,
        MOSHIER_VENUS,
        MOSHIER_MARS,
        MOSHIER_JUPITER,
        MOSHIER_SATURN,
        MOSHIER_URANUS,
        MOSHIER_NEPTUNE,
        MOSHIER_PLUTO,
        MOSHIER_EARTH,
    ]
)


def is_moshier_body(body_id: int) -> bool:
    """Check if a body is supported by the Moshier ephemeris.

    Args:
        body_id: Swiss Ephemeris body ID.

    Returns:
        True if the body can be calculated with Moshier algorithms.
    """
    return body_id in MOSHIER_BODIES


# =============================================================================
# MAIN CALCULATION FUNCTION
# =============================================================================


def calc_position(
    jd_tt: float,
    body_id: int,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate geocentric ecliptic position for a body using Moshier algorithms.

    This is the main entry point for all Moshier calculations. It dispatches
    to the appropriate theory based on the body ID.

    Args:
        jd_tt: Julian Day in Terrestrial Time.
        body_id: Swiss Ephemeris body ID (SE_SUN, SE_MOON, SE_MERCURY, etc.)

    Returns:
        Tuple of (lon, lat, dist, dlon, dlat, ddist) where:
        - lon: Ecliptic longitude in degrees (J2000.0)
        - lat: Ecliptic latitude in degrees (J2000.0)
        - dist: Distance in AU (or km for Moon, converted to AU)
        - dlon: Longitude velocity in degrees/day
        - dlat: Latitude velocity in degrees/day
        - ddist: Distance velocity in AU/day

    Raises:
        ValueError: If the body ID is not supported.

    Example:
        >>> from libephemeris.moshier import calc_position
        >>> from libephemeris.constants import SE_MARS
        >>> lon, lat, dist, dlon, dlat, ddist = calc_position(2451545.0, SE_MARS)
        >>> print(f"Mars longitude: {lon:.4f}°")
    """
    # Dispatch to appropriate theory
    if is_moon_body(body_id):
        return elp_calc_position(jd_tt, body_id)
    elif is_pluto_body(body_id):
        return pluto_calc_position(jd_tt, body_id)
    elif is_vsop_body(body_id):
        return vsop_calc_position(jd_tt, body_id)
    else:
        raise ValueError(
            f"Body {body_id} is not supported by Moshier ephemeris. "
            f"Supported bodies: Sun (0), Moon (1), Mercury (2), Venus (3), "
            f"Mars (4), Jupiter (5), Saturn (6), Uranus (7), Neptune (8), "
            f"Pluto (9), Earth (14)."
        )


# =============================================================================
# MODULE EXPORTS
# =============================================================================

__all__ = [
    # Main function
    "calc_position",
    "is_moshier_body",
    # Body ID constants
    "MOSHIER_SUN",
    "MOSHIER_MOON",
    "MOSHIER_MERCURY",
    "MOSHIER_VENUS",
    "MOSHIER_EARTH",
    "MOSHIER_MARS",
    "MOSHIER_JUPITER",
    "MOSHIER_SATURN",
    "MOSHIER_URANUS",
    "MOSHIER_NEPTUNE",
    "MOSHIER_PLUTO",
    "MOSHIER_BODIES",
    # Theory-specific check functions
    "is_vsop_body",
    "is_moon_body",
    "is_pluto_body",
    # Precession/Nutation
    "mean_obliquity",
    "true_obliquity",
    "nutation_angles",
    "precession_angles",
    "precess_ecliptic",
    "precess_to_j2000",
    "precess_from_j2000",
    # Constants
    "J2000",
    "JD_PER_CENTURY",
    "AU_KM",
    "C_LIGHT_AU_DAY",
    "DEG_TO_RAD",
    "RAD_TO_DEG",
    "ARCSEC_TO_RAD",
    "RAD_TO_ARCSEC",
    # Utility functions
    "normalize_angle",
    "normalize_radians",
    "spherical_to_cartesian",
    "cartesian_to_spherical",
    "ecliptic_to_equatorial",
    "equatorial_to_ecliptic",
    "rectangular_to_spherical_velocity",
    "annual_aberration",
    "light_time_correction",
    "jd_to_julian_centuries",
    "jd_to_julian_millennia",
    "julian_centuries_to_jd",
    "numerical_derivative",
    "evaluate_polynomial",
    "evaluate_polynomial_derivative",
]
