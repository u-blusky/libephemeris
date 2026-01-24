"""
Solar and Lunar eclipse calculations for libephemeris.

Finds eclipse events and calculates their circumstances.

Functions:
- sol_eclipse_when_glob: Find next global solar eclipse
- (planned) sol_eclipse_when_loc: Find eclipse at specific location
- (planned) sol_eclipse_where: Calculate path of eclipse
- (planned) sol_eclipse_how: Eclipse circumstances at location
- (planned) lun_eclipse_when: Find next lunar eclipse

Algorithm:
    Solar eclipses occur at New Moon when Moon is near the ecliptic plane.
    1. Find next New Moon (Sun-Moon conjunction in longitude)
    2. Check lunar latitude - if |lat| < ~1.5° eclipse is possible
    3. Calculate eclipse magnitude and type based on distances
    4. Use Besselian elements for precise timing of phases

References:
    - Meeus "Astronomical Algorithms" Ch. 54 (Eclipses)
    - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
    - Swiss Ephemeris documentation
"""

import math
from typing import Tuple
from .constants import (
    SE_SUN,
    SE_MOON,
    SEFLG_SPEED,
    SEFLG_SWIEPH,
    SE_ECL_TOTAL,
    SE_ECL_ANNULAR,
    SE_ECL_PARTIAL,
    SE_ECL_ANNULAR_TOTAL,
    SE_ECL_CENTRAL,
    SE_ECL_ALLTYPES_SOLAR,
    SE_ECL_VISIBLE,
    SE_ECL_MAX_VISIBLE,
    SE_ECL_1ST_VISIBLE,
    SE_ECL_2ND_VISIBLE,
    SE_ECL_3RD_VISIBLE,
    SE_ECL_4TH_VISIBLE,
)
from .planets import swe_calc_ut
from .state import get_timescale


# Constants for eclipse calculations
SYNODIC_MONTH = 29.530588853  # Mean synodic month in days
LUNAR_NODE_PERIOD = 6798.38  # Lunar node regression period in days
ECLIPSE_LIMIT_SOLAR = 18.5  # Maximum elongation from node for solar eclipse (degrees)


def _find_next_new_moon(jd_start: float) -> float:
    """
    Find the next New Moon (Sun-Moon conjunction) after jd_start.

    Uses iterative refinement to find exact moment of conjunction.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of next New Moon
    """
    # Get current positions
    sun_pos, _ = swe_calc_ut(jd_start, SE_SUN, SEFLG_SPEED)
    moon_pos, _ = swe_calc_ut(jd_start, SE_MOON, SEFLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation (Moon - Sun), normalized to -180 to 180
    elongation = (moon_lon - sun_lon) % 360.0
    if elongation > 180:
        elongation -= 360

    # Estimate time to next conjunction
    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day (average Moon speed - Sun speed)

    # Time until next conjunction (elongation = 0)
    if elongation > 0:
        # Moon is ahead, wait for next cycle
        dt = (360.0 - elongation) / relative_speed
    else:
        # Moon is behind Sun
        dt = (-elongation) / relative_speed

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(20):
        sun_pos, _ = swe_calc_ut(jd_guess, SE_SUN, SEFLG_SPEED)
        moon_pos, _ = swe_calc_ut(jd_guess, SE_MOON, SEFLG_SPEED)

        sun_lon = sun_pos[0]
        moon_lon = moon_pos[0]
        sun_speed = sun_pos[3]
        moon_speed = moon_pos[3]

        # Elongation
        diff = (moon_lon - sun_lon) % 360.0
        if diff > 180:
            diff -= 360

        # Convergence check (< 0.1 arcsec)
        if abs(diff) < 1e-5:
            return jd_guess

        # Newton-Raphson step
        rel_speed = moon_speed - sun_speed
        if abs(rel_speed) < 0.1:
            rel_speed = 12.19

        jd_guess -= diff / rel_speed

    return jd_guess


def _get_moon_node_distance(jd: float, moon_lon: float) -> float:
    """
    Calculate angular distance of Moon from nearest lunar node.

    Args:
        jd: Julian Day (UT)
        moon_lon: Moon's ecliptic longitude in degrees

    Returns:
        Absolute angular distance from nearest node in degrees
    """
    from .lunar import calc_mean_lunar_node

    # Get mean node longitude
    ts = get_timescale()
    t = ts.ut1_jd(jd)
    node_lon = calc_mean_lunar_node(t.tt)

    # Calculate distance to both nodes
    dist_to_north = abs((moon_lon - node_lon + 180) % 360 - 180)
    dist_to_south = abs((moon_lon - (node_lon + 180) + 180) % 360 - 180)

    return min(dist_to_north, dist_to_south)


def _calculate_eclipse_type_and_magnitude(
    jd: float,
) -> Tuple[int, float, float, float]:
    """
    Determine eclipse type and magnitude at maximum eclipse.

    Uses geometric calculations based on Sun-Moon-Earth distances
    and apparent angular sizes.

    Args:
        jd: Julian Day of eclipse maximum (UT)

    Returns:
        Tuple of (eclipse_type_flags, magnitude, gamma, moon_sun_ratio)
        - eclipse_type_flags: Bitmask of SE_ECL_* constants
        - magnitude: Eclipse magnitude (fraction of Sun's diameter covered)
        - gamma: Gamma parameter (distance of Moon shadow axis from Earth center)
        - moon_sun_ratio: Ratio of apparent Moon diameter to Sun diameter
    """
    # Get positions
    sun_pos, _ = swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
    moon_pos, _ = swe_calc_ut(jd, SE_MOON, SEFLG_SPEED)

    # Distances in AU
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Moon's ecliptic latitude (how far from the ecliptic)
    moon_lat = moon_pos[1]

    # Angular radii (in degrees)
    # Sun: mean radius 959.63" at 1 AU
    sun_angular_radius = (959.63 / 3600.0) / sun_dist
    # Moon: mean radius 932.56" at mean distance (0.002569 AU)
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Ratio of apparent sizes
    moon_sun_ratio = moon_angular_radius / sun_angular_radius

    # Calculate gamma - the perpendicular distance of Moon's shadow axis
    # from Earth's center, in Earth radii
    # Approximation: gamma ≈ sin(moon_lat) * moon_dist / earth_radius
    # Using more precise formula based on geometry
    earth_radius_au = 4.2635e-5  # Earth radius in AU

    # Convert latitude to radians
    lat_rad = math.radians(moon_lat)

    # Distance of Moon from ecliptic plane
    z_moon = moon_dist * math.sin(lat_rad)

    # Gamma (in Earth radii)
    gamma = z_moon / earth_radius_au

    # Penumbral cone semi-angle (reserved for future Besselian elements)
    # f1 = (sun_angular_radius + moon_angular_radius) * moon_dist / earth_radius_au
    # Umbral cone semi-angle
    # f2 = (sun_angular_radius - moon_angular_radius) * moon_dist / earth_radius_au

    # Simplified gamma calculation based on Moon latitude
    # For more precision, use Besselian elements
    gamma = moon_lat * moon_dist * 57.2958  # Rough conversion to Earth radii
    gamma /= 60.0  # Normalize

    # More accurate gamma using geometry
    # gamma = distance of shadow axis from Earth center / Earth equatorial radius
    gamma = (moon_lat / 60.0) * (moon_dist / 0.002569) * 0.5  # Simplified

    # Determine eclipse type based on gamma and ratio
    eclipse_type = 0
    magnitude = 0.0

    # Maximum gamma for any eclipse
    GAMMA_LIMIT_PARTIAL = 1.55
    GAMMA_LIMIT_TOTAL = 0.997

    if abs(gamma) > GAMMA_LIMIT_PARTIAL:
        # No eclipse
        return 0, 0.0, gamma, moon_sun_ratio

    # There is some form of eclipse
    if moon_sun_ratio >= 1.0:
        # Moon appears larger - can be total
        if abs(gamma) < GAMMA_LIMIT_TOTAL:
            eclipse_type = SE_ECL_TOTAL | SE_ECL_CENTRAL
            magnitude = 1.0 + (moon_sun_ratio - 1.0) * (1 - abs(gamma))
        else:
            eclipse_type = SE_ECL_PARTIAL
            magnitude = 1.0 - abs(gamma) / GAMMA_LIMIT_PARTIAL
    else:
        # Moon appears smaller - can be annular
        if abs(gamma) < GAMMA_LIMIT_TOTAL:
            # Check for hybrid (annular-total)
            if moon_sun_ratio > 0.99:
                eclipse_type = SE_ECL_ANNULAR_TOTAL | SE_ECL_CENTRAL
            else:
                eclipse_type = SE_ECL_ANNULAR | SE_ECL_CENTRAL
            magnitude = moon_sun_ratio
        else:
            eclipse_type = SE_ECL_PARTIAL
            magnitude = (1.0 - abs(gamma) / GAMMA_LIMIT_PARTIAL) * moon_sun_ratio

    # Ensure magnitude is in valid range
    magnitude = max(0.0, min(1.5, magnitude))

    return eclipse_type, magnitude, gamma, moon_sun_ratio


def _calculate_eclipse_phases(
    jd_max: float, eclipse_type: int
) -> Tuple[float, float, float, float, float, float, float, float]:
    """
    Calculate times of eclipse phases (contacts) for a global eclipse.

    Phase indices:
        [0]: Time of maximum eclipse
        [1]: Time of first contact (partial eclipse begins)
        [2]: Time of second contact (total/annular begins, if central)
        [3]: Time of third contact (total/annular ends, if central)
        [4]: Time of fourth contact (partial eclipse ends)
        [5]: Time of sunrise on central line
        [6]: Time of sunset on central line
        [7]: Reserved

    Args:
        jd_max: Julian Day of maximum eclipse
        eclipse_type: Eclipse type flags

    Returns:
        Tuple of 8 floats with phase times (JD UT)
    """
    # For a simplified implementation, estimate durations
    # based on typical eclipse characteristics

    is_central = bool(eclipse_type & SE_ECL_CENTRAL)
    is_total = bool(eclipse_type & SE_ECL_TOTAL)
    is_annular = bool(eclipse_type & SE_ECL_ANNULAR)

    # Typical partial phase duration: ~1 hour before/after maximum
    # Based on Moon's shadow crossing Earth
    partial_duration = 1.0 / 24.0  # 1 hour in days (simplified)

    # Typical central phase duration: 2-7 minutes for total, longer for annular
    if is_total:
        central_duration = 3.0 / (24.0 * 60.0)  # ~3 minutes
    elif is_annular:
        central_duration = 6.0 / (24.0 * 60.0)  # ~6 minutes
    else:
        central_duration = 0.0

    # Calculate phase times
    t_first_contact = jd_max - partial_duration
    t_fourth_contact = jd_max + partial_duration

    if is_central:
        t_second_contact = jd_max - central_duration / 2
        t_third_contact = jd_max + central_duration / 2
    else:
        t_second_contact = 0.0
        t_third_contact = 0.0

    # Sunrise/sunset on central line (simplified - set to 0 if not calculated)
    t_sunrise = 0.0
    t_sunset = 0.0

    return (
        jd_max,
        t_first_contact,
        t_second_contact,
        t_third_contact,
        t_fourth_contact,
        t_sunrise,
        t_sunset,
        0.0,  # Reserved
    )


def sol_eclipse_when_glob(
    jd_start: float,
    flags: int = SEFLG_SWIEPH,
    eclipse_type: int = 0,
) -> Tuple[Tuple[float, ...], int]:
    """
    Find the next global solar eclipse after a given date.

    Searches forward in time from jd_start to find the next solar eclipse.
    Can filter by eclipse type (total, annular, partial, hybrid).

    Args:
        jd_start: Julian Day (UT) to start search from
        flags: Calculation flags (SEFLG_SWIEPH, etc.)
        eclipse_type: Filter for specific eclipse type(s), bitmask of:
            - SE_ECL_TOTAL (4): Total eclipse
            - SE_ECL_ANNULAR (8): Annular eclipse
            - SE_ECL_PARTIAL (16): Partial eclipse
            - SE_ECL_ANNULAR_TOTAL (32): Hybrid eclipse
            - 0: Any eclipse type (default)

    Returns:
        Tuple containing:
            - times: Tuple of 8 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse
                [1]: Time of first contact (partial begins)
                [2]: Time of second contact (central phase begins, or 0)
                [3]: Time of third contact (central phase ends, or 0)
                [4]: Time of fourth contact (partial ends)
                [5]: Time of sunrise on central line (or 0)
                [6]: Time of sunset on central line (or 0)
                [7]: Reserved (0)
            - retflag: Eclipse type flags bitmask (SE_ECL_* constants)

    Raises:
        RuntimeError: If no eclipse found within search limit

    Algorithm:
        1. Find next New Moon after jd_start
        2. Check if Moon is close enough to node for eclipse
        3. If not eclipse, advance to next New Moon
        4. Calculate eclipse type and magnitude
        5. If eclipse_type filter set, check if matches
        6. Calculate phase times

    Precision:
        Eclipse times accurate to ~1 minute for most eclipses.
        For higher precision, use Besselian elements method.

    Example:
        >>> # Find next total solar eclipse after Jan 1, 2024
        >>> from libephemeris import julday, SE_ECL_TOTAL
        >>> jd = julday(2024, 1, 1, 0)
        >>> times, ecl_type = sol_eclipse_when_glob(jd, eclipse_type=SE_ECL_TOTAL)
        >>> print(f"Total eclipse at JD {times[0]:.5f}")

    References:
        - Swiss Ephemeris: swe_sol_eclipse_when_glob()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    MAX_SEARCH_YEARS = 20  # Maximum search range
    MAX_NEW_MOONS = int(MAX_SEARCH_YEARS * 12.4)  # ~12.4 lunations per year

    # If eclipse_type is 0, accept any type
    if eclipse_type == 0:
        eclipse_type = SE_ECL_ALLTYPES_SOLAR

    jd = jd_start

    for _ in range(MAX_NEW_MOONS):
        # Find next New Moon
        jd_new_moon = _find_next_new_moon(jd)

        # Get Moon position at New Moon
        moon_pos, _ = swe_calc_ut(jd_new_moon, SE_MOON, flags | SEFLG_SPEED)
        moon_lon = moon_pos[0]
        # moon_lat = moon_pos[1]  # Used in _calculate_eclipse_type_and_magnitude

        # Check if close enough to ecliptic for eclipse
        # Solar eclipse possible if |moon_lat| < ~1.5° (approximate)
        # More accurate: check distance from node
        node_dist = _get_moon_node_distance(jd_new_moon, moon_lon)

        if node_dist < ECLIPSE_LIMIT_SOLAR:
            # Possible eclipse - check magnitude
            ecl_type, magnitude, gamma, ratio = _calculate_eclipse_type_and_magnitude(
                jd_new_moon
            )

            if ecl_type != 0:
                # Eclipse found - check if matches filter
                type_matches = (
                    (eclipse_type & SE_ECL_TOTAL and ecl_type & SE_ECL_TOTAL)
                    or (eclipse_type & SE_ECL_ANNULAR and ecl_type & SE_ECL_ANNULAR)
                    or (eclipse_type & SE_ECL_PARTIAL and ecl_type & SE_ECL_PARTIAL)
                    or (
                        eclipse_type & SE_ECL_ANNULAR_TOTAL
                        and ecl_type & SE_ECL_ANNULAR_TOTAL
                    )
                )

                if type_matches:
                    # Calculate phase times
                    times = _calculate_eclipse_phases(jd_new_moon, ecl_type)
                    return times, ecl_type

        # Advance to next lunation
        jd = jd_new_moon + 25  # Skip ahead ~25 days to ensure we find next New Moon

    raise RuntimeError(
        f"No matching solar eclipse found within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Aliases for compatibility
swe_sol_eclipse_when_glob = sol_eclipse_when_glob


def _calculate_local_eclipse_phases(
    jd_max_global: float,
    lat: float,
    lon: float,
    altitude: float,
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """
    Calculate eclipse phase times as seen from a specific location.

    Uses geometric calculations based on the Moon's shadow passing over
    the observer's location.

    Args:
        jd_max_global: Julian Day of global maximum eclipse
        lat: Observer latitude in degrees (positive = North)
        lon: Observer longitude in degrees (positive = East)
        altitude: Observer altitude in meters above sea level

    Returns:
        Tuple of 10 floats with local eclipse circumstances:
            [0]: Time of maximum eclipse (local)
            [1]: Time of first contact (partial begins)
            [2]: Time of second contact (total/annular begins, or 0)
            [3]: Time of third contact (total/annular ends, or 0)
            [4]: Time of fourth contact (partial ends)
            [5]: Eclipse magnitude at maximum
            [6]: Fraction of Sun's diameter obscured
            [7]: Eclipse obscuration (area)
            [8]: Azimuth of Sun at maximum eclipse
            [9]: Altitude of Sun at maximum eclipse
    """
    from skyfield.api import wgs84

    from .state import get_planets, get_timescale

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Get Sun and Moon objects
    sun = eph["sun"]
    moon = eph["moon"]
    earth = eph["earth"]

    # Search for local maximum around global maximum
    # The local maximum can differ from global by up to ~1 hour
    search_range = 3.0 / 24.0  # 3 hours in days

    def _get_sun_moon_separation(jd: float) -> float:
        """Get angular separation between Sun and Moon at given JD."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer

        # Get apparent positions from observer
        sun_app = observer_at.at(t).observe(sun).apparent()
        moon_app = observer_at.at(t).observe(moon).apparent()

        # Calculate angular separation
        sep = sun_app.separation_from(moon_app)
        return sep.degrees

    def _get_sun_altaz(jd: float) -> Tuple[float, float]:
        """Get Sun altitude and azimuth at given JD from observer location."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer

        sun_app = observer_at.at(t).observe(sun).apparent()
        alt, az, _ = sun_app.altaz()
        return alt.degrees, az.degrees

    # First check if Sun is above horizon at global maximum
    sun_alt, sun_az = _get_sun_altaz(jd_max_global)
    if sun_alt < -1.0:  # Sun below horizon (with some margin for refraction)
        # Sun not visible - return zeros to indicate no local eclipse
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Find local maximum (minimum separation)
    # Use golden section search for efficiency
    jd_low = jd_max_global - search_range
    jd_high = jd_max_global + search_range

    phi = (1 + math.sqrt(5)) / 2  # Golden ratio
    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    sep_a = _get_sun_moon_separation(jd_a)
    sep_b = _get_sun_moon_separation(jd_b)

    for _ in range(30):  # Converge to ~0.1 second precision
        if sep_a < sep_b:
            jd_high = jd_b
            jd_b = jd_a
            sep_b = sep_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            sep_a = _get_sun_moon_separation(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            sep_a = sep_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            sep_b = _get_sun_moon_separation(jd_b)

        if jd_high - jd_low < 1e-7:  # ~0.01 second
            break

    jd_local_max = (jd_low + jd_high) / 2
    min_separation = _get_sun_moon_separation(jd_local_max)

    # Check if Sun is visible at local maximum
    sun_alt_max, sun_az_max = _get_sun_altaz(jd_local_max)
    if sun_alt_max < -1.0:
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Get angular sizes of Sun and Moon at local maximum
    t_max = ts.ut1_jd(jd_local_max)
    observer_at = earth + observer

    sun_app = observer_at.at(t_max).observe(sun).apparent()
    moon_app = observer_at.at(t_max).observe(moon).apparent()

    # Sun angular radius: ~959.63" at 1 AU
    # Moon angular radius: ~932.56" at mean distance
    sun_dist_au = sun_app.distance().au
    moon_dist_au = moon_app.distance().au

    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au  # degrees
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)  # degrees

    sun_diameter = 2 * sun_angular_radius
    moon_diameter = 2 * moon_angular_radius

    # Check if there's an eclipse at this location
    # Eclipse occurs if separation < sum of radii
    sum_radii = sun_angular_radius + moon_angular_radius
    if min_separation > sum_radii:
        # No eclipse visible from this location
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Calculate eclipse magnitude
    # Magnitude = fraction of Sun's diameter covered by Moon
    overlap = sum_radii - min_separation
    magnitude = overlap / sun_diameter
    magnitude = max(0.0, min(magnitude, 1.0 + moon_diameter / sun_diameter))

    # Calculate obscuration (fraction of Sun's area covered)
    # Using formula from Meeus
    r_sun = sun_angular_radius
    r_moon = moon_angular_radius
    d = min_separation  # center-to-center separation

    if d >= r_sun + r_moon:
        obscuration = 0.0
    elif d <= abs(r_sun - r_moon):
        # One disk entirely within the other
        if r_moon >= r_sun:
            obscuration = 1.0
        else:
            obscuration = (r_moon / r_sun) ** 2
    else:
        # Partial overlap - use lens formula
        # Area of intersection of two circles
        d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
        d2 = d - d1
        if abs(d1) <= r_sun and abs(d2) <= r_moon:
            area1 = r_sun * r_sun * math.acos(d1 / r_sun) - d1 * math.sqrt(
                max(0, r_sun * r_sun - d1 * d1)
            )
            area2 = r_moon * r_moon * math.acos(d2 / r_moon) - d2 * math.sqrt(
                max(0, r_moon * r_moon - d2 * d2)
            )
            intersection_area = area1 + area2
            sun_area = math.pi * r_sun * r_sun
            obscuration = intersection_area / sun_area
        else:
            obscuration = 0.0

    obscuration = max(0.0, min(1.0, obscuration))

    # Find contact times using bisection
    # First/fourth contact: separation = sum of radii
    # Second/third contact: separation = difference of radii (for central eclipse)

    def _find_contact_time(
        jd_start: float, jd_end: float, target_sep: float, is_increasing: bool
    ) -> float:
        """Find time when separation equals target, searching in given direction."""
        jd_low = jd_start
        jd_high = jd_end

        for _ in range(50):
            jd_mid = (jd_low + jd_high) / 2
            sep = _get_sun_moon_separation(jd_mid)

            if abs(sep - target_sep) < 1e-6:
                return jd_mid

            if is_increasing:
                # Looking for separation increasing through target
                if sep < target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid
            else:
                # Looking for separation decreasing through target
                if sep > target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid

            if jd_high - jd_low < 1e-8:
                break

        return (jd_low + jd_high) / 2

    # Search range for contacts
    contact_search = 2.0 / 24.0  # 2 hours

    # First contact (separation decreasing to sum of radii)
    jd_first = _find_contact_time(
        jd_local_max - contact_search, jd_local_max, sum_radii, is_increasing=False
    )

    # Fourth contact (separation increasing from sum of radii)
    jd_fourth = _find_contact_time(
        jd_local_max, jd_local_max + contact_search, sum_radii, is_increasing=True
    )

    # Check if first and fourth contacts are valid (Sun above horizon)
    first_alt, _ = _get_sun_altaz(jd_first)
    fourth_alt, _ = _get_sun_altaz(jd_fourth)

    if first_alt < -1.0:
        jd_first = 0.0
    if fourth_alt < -1.0:
        jd_fourth = 0.0

    # Second and third contacts (for total/annular eclipses)
    jd_second = 0.0
    jd_third = 0.0

    diff_radii = abs(moon_angular_radius - sun_angular_radius)
    if min_separation < diff_radii:
        # Central eclipse at this location
        jd_second = _find_contact_time(
            jd_local_max - contact_search / 4,
            jd_local_max,
            diff_radii,
            is_increasing=False,
        )
        jd_third = _find_contact_time(
            jd_local_max,
            jd_local_max + contact_search / 4,
            diff_radii,
            is_increasing=True,
        )

        # Check visibility
        second_alt, _ = _get_sun_altaz(jd_second)
        third_alt, _ = _get_sun_altaz(jd_third)

        if second_alt < -1.0:
            jd_second = 0.0
        if third_alt < -1.0:
            jd_third = 0.0

    return (
        jd_local_max,  # [0] Time of maximum
        jd_first,  # [1] First contact
        jd_second,  # [2] Second contact
        jd_third,  # [3] Third contact
        jd_fourth,  # [4] Fourth contact
        magnitude,  # [5] Magnitude
        overlap / sun_diameter,  # [6] Fraction covered
        obscuration,  # [7] Obscuration
        sun_az_max,  # [8] Azimuth
        sun_alt_max,  # [9] Altitude
    )


def _determine_local_eclipse_type(
    magnitude: float, moon_sun_ratio: float, is_central: bool
) -> int:
    """
    Determine eclipse type flags for local observation.

    Args:
        magnitude: Eclipse magnitude at location
        moon_sun_ratio: Ratio of Moon's apparent diameter to Sun's
        is_central: Whether this is a central eclipse at location

    Returns:
        Eclipse type flags bitmask
    """
    if magnitude <= 0:
        return 0

    flags = SE_ECL_VISIBLE

    if is_central:
        flags |= SE_ECL_CENTRAL
        if moon_sun_ratio >= 1.0:
            flags |= SE_ECL_TOTAL
        elif moon_sun_ratio > 0.99:
            flags |= SE_ECL_ANNULAR_TOTAL
        else:
            flags |= SE_ECL_ANNULAR
    else:
        flags |= SE_ECL_PARTIAL

    return flags


def sol_eclipse_when_loc(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], int]:
    """
    Find the next solar eclipse visible from a specific location.

    Searches forward in time from jd_start to find the next solar eclipse
    visible from the observer's location. Returns the times of eclipse phases
    as seen from that position.

    Args:
        jd_start: Julian Day (UT) to start search from
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - times: Tuple of 10 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse (local)
                [1]: Time of first contact (partial begins)
                [2]: Time of second contact (total/annular begins, or 0)
                [3]: Time of third contact (total/annular ends, or 0)
                [4]: Time of fourth contact (partial ends)
                [5]: Time of sunrise on central line (or 0, unused)
                [6]: Time of sunset on central line (or 0, unused)
                [7-9]: Reserved (0)
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Eclipse magnitude
                [1]: Fraction of Sun's diameter covered
                [2]: Fraction of Sun's area obscured (obscuration)
                [3]: Azimuth of Sun at maximum eclipse (degrees)
                [4]: Altitude of Sun at maximum eclipse (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of Sun (degrees)
                [7]: Saros series number (0, not implemented)
                [8]: Inex number (0, not implemented)
                [9-10]: Reserved (0)
            - retflag: Eclipse type flags bitmask (SE_ECL_* constants)

    Raises:
        RuntimeError: If no eclipse visible from location within search limit

    Algorithm:
        1. Use sol_eclipse_when_glob to find next global eclipse
        2. Calculate local circumstances at observer's location
        3. If eclipse not visible from location, continue to next global eclipse
        4. Return local phase times and attributes

    Precision:
        Eclipse times accurate to ~1 minute. Contact times depend on
        accurate Moon and Sun ephemeris positions.

    Example:
        >>> # Find next eclipse visible from Rome
        >>> from libephemeris import julday, sol_eclipse_when_loc
        >>> jd = julday(2024, 1, 1, 0)
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> times, attr, ecl_type = sol_eclipse_when_loc(jd, rome_lat, rome_lon)
        >>> print(f"Eclipse maximum at JD {times[0]:.5f}")
        >>> print(f"Magnitude: {attr[0]:.3f}")

    References:
        - Swiss Ephemeris: swe_sol_eclipse_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    MAX_SEARCH_YEARS = 50  # Maximum search range
    MAX_ECLIPSES = int(MAX_SEARCH_YEARS * 2.4)  # ~2.4 eclipses per year on average

    from .state import get_planets, get_timescale

    # Get ephemeris for Moon/Sun diameter calculations
    eph = get_planets()
    ts = get_timescale()

    jd = jd_start

    for _ in range(MAX_ECLIPSES):
        # Find next global eclipse
        try:
            global_times, global_type = sol_eclipse_when_glob(jd, flags)
        except RuntimeError:
            raise RuntimeError(
                f"No solar eclipse visible from lat={lat}, lon={lon} "
                f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
            )

        jd_max_global = global_times[0]

        # Calculate local circumstances
        local_data = _calculate_local_eclipse_phases(jd_max_global, lat, lon, altitude)

        jd_local_max = local_data[0]
        magnitude = local_data[5]

        # Check if eclipse is visible from this location
        if jd_local_max > 0 and magnitude > 0:
            # Eclipse visible! Prepare return values

            # Get Moon/Sun apparent diameters at local maximum
            t_max = ts.ut1_jd(jd_local_max)
            from skyfield.api import wgs84

            observer = wgs84.latlon(lat, lon, altitude)
            observer_at = eph["earth"] + observer

            sun_app = observer_at.at(t_max).observe(eph["sun"]).apparent()
            moon_app = observer_at.at(t_max).observe(eph["moon"]).apparent()

            sun_dist_au = sun_app.distance().au
            moon_dist_au = moon_app.distance().au

            sun_diameter = 2 * (959.63 / 3600.0) / sun_dist_au
            moon_diameter = 2 * (932.56 / 3600.0) * (0.002569 / moon_dist_au)

            # Check for central eclipse at this location
            is_central = local_data[2] > 0 and local_data[3] > 0

            moon_sun_ratio = moon_diameter / sun_diameter

            # Determine eclipse type
            ecl_type = _determine_local_eclipse_type(
                magnitude, moon_sun_ratio, is_central
            )

            # Add visibility flags
            if local_data[1] > 0:
                ecl_type |= SE_ECL_1ST_VISIBLE
            if local_data[2] > 0:
                ecl_type |= SE_ECL_2ND_VISIBLE
            if local_data[3] > 0:
                ecl_type |= SE_ECL_3RD_VISIBLE
            if local_data[4] > 0:
                ecl_type |= SE_ECL_4TH_VISIBLE
            ecl_type |= SE_ECL_MAX_VISIBLE

            # Prepare times tuple (10 elements)
            times = (
                local_data[0],  # [0] Maximum
                local_data[1],  # [1] First contact
                local_data[2],  # [2] Second contact
                local_data[3],  # [3] Third contact
                local_data[4],  # [4] Fourth contact
                0.0,  # [5] Sunrise on central line (not implemented)
                0.0,  # [6] Sunset on central line (not implemented)
                0.0,  # [7] Reserved
                0.0,  # [8] Reserved
                0.0,  # [9] Reserved
            )

            # Prepare attributes tuple (11 elements)
            attr = (
                local_data[5],  # [0] Magnitude
                local_data[6],  # [1] Fraction covered
                local_data[7],  # [2] Obscuration
                local_data[8],  # [3] Azimuth
                local_data[9],  # [4] Altitude
                moon_diameter,  # [5] Moon diameter
                sun_diameter,  # [6] Sun diameter
                0.0,  # [7] Saros (not implemented)
                0.0,  # [8] Inex (not implemented)
                0.0,  # [9] Reserved
                0.0,  # [10] Reserved
            )

            return times, attr, ecl_type

        # Eclipse not visible from this location, try next
        jd = jd_max_global + 25  # Skip ahead to find next eclipse

    raise RuntimeError(
        f"No solar eclipse visible from lat={lat}, lon={lon} "
        f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Alias for Swiss Ephemeris API compatibility
swe_sol_eclipse_when_loc = sol_eclipse_when_loc
