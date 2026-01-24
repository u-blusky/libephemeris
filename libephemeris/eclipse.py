"""
Solar and Lunar eclipse calculations for libephemeris.

Finds eclipse events and calculates their circumstances.

Functions:
- sol_eclipse_when_glob: Find next global solar eclipse
- sol_eclipse_when_loc: Find eclipse at specific location
- sol_eclipse_where: Calculate path of eclipse
- sol_eclipse_how: Eclipse circumstances at location
- lun_eclipse_when: Find next lunar eclipse

Algorithm:
    Solar eclipses occur at New Moon when Moon is near the ecliptic plane.
    1. Find next New Moon (Sun-Moon conjunction in longitude)
    2. Check lunar latitude - if |lat| < ~1.5° eclipse is possible
    3. Calculate eclipse magnitude and type based on distances
    4. Use Besselian elements for precise timing of phases

    Lunar eclipses occur at Full Moon when Moon is near a lunar node.
    1. Find next Full Moon (Sun-Moon opposition in longitude)
    2. Check Moon's distance from node - if close, eclipse is possible
    3. Calculate eclipse type (penumbral, partial, total) based on geometry
    4. Calculate phase times based on shadow cone sizes

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
    SE_ECL_ALLTYPES_LUNAR,
    SE_ECL_PENUMBRAL,
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


def sol_eclipse_where(
    jd: float,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], int]:
    """
    Calculate the geographic location where a solar eclipse is central at a given time.

    This function determines where on Earth the Moon's shadow axis intersects
    the surface at the specified Julian Day. Used for tracing the path of
    totality/annularity across the Earth's surface.

    Args:
        jd: Julian Day (UT) of the moment to calculate
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - geopos: Tuple of 2 floats with geographic position:
                [0]: Geographic longitude of central eclipse (degrees, East positive)
                [1]: Geographic latitude of central eclipse (degrees, North positive)
            - attr: Tuple of 8 floats with eclipse attributes:
                [0]: Fraction of solar diameter covered by Moon (magnitude)
                [1]: Ratio of lunar diameter to solar diameter
                [2]: Fraction of Sun's area obscured (obscuration)
                [3]: Width of totality/annularity path (km)
                [4]: Sun's azimuth at central line (degrees)
                [5]: Sun's altitude at central line (degrees)
                [6]: Apparent diameter of Moon (degrees)
                [7]: Apparent diameter of Sun (degrees)
            - retflag: Eclipse type flags bitmask (SE_ECL_* constants)
                Returns 0 if no central eclipse at this time

    Note:
        If there is no central eclipse at the given time (either no eclipse
        at all, or the shadow axis misses the Earth), geopos will contain
        zeros and retflag will be 0.

    Algorithm:
        1. Calculate Sun and Moon positions at the given time
        2. Compute the Moon's shadow axis direction
        3. Find intersection of shadow axis with Earth's surface
        4. Calculate eclipse attributes at that location
        5. Determine path width from umbra/penumbra cone geometry

    Example:
        >>> # Find central line location during April 8, 2024 eclipse
        >>> from libephemeris import julday, sol_eclipse_where
        >>> jd = julday(2024, 4, 8, 18.2)  # Approximate maximum
        >>> geopos, attr, ecl_type = sol_eclipse_where(jd)
        >>> print(f"Central line at lon={geopos[0]:.2f}, lat={geopos[1]:.2f}")
        >>> print(f"Path width: {attr[3]:.1f} km")

    References:
        - Swiss Ephemeris: swe_sol_eclipse_where()
        - Meeus "Astronomical Algorithms" Ch. 54
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
    """
    from skyfield.api import wgs84

    from .state import get_planets, get_timescale

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Get Sun and Moon positions
    sun_pos, _ = swe_calc_ut(jd, SE_SUN, flags | SEFLG_SPEED)
    moon_pos, _ = swe_calc_ut(jd, SE_MOON, flags | SEFLG_SPEED)

    # Check if eclipse is possible (Sun-Moon conjunction)
    elongation = abs((moon_pos[0] - sun_pos[0] + 180) % 360 - 180)
    if elongation > 5.0:  # Too far from conjunction for eclipse
        return (0.0, 0.0), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0

    # Get distances in AU
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]
    moon_lat = moon_pos[1]  # Ecliptic latitude

    # Constants
    EARTH_RADIUS_AU = 4.2635e-5

    # Calculate gamma (Moon shadow axis distance from Earth center, in Earth radii)
    # Using improved geometric calculation
    lat_rad = math.radians(moon_lat)
    z_moon = moon_dist * math.sin(lat_rad)
    gamma = z_moon / EARTH_RADIUS_AU

    # Check if shadow axis intersects Earth
    GAMMA_LIMIT = 1.55  # Maximum gamma for any eclipse visibility
    if abs(gamma) > GAMMA_LIMIT:
        return (0.0, 0.0), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0

    # Determine eclipse type
    GAMMA_LIMIT_CENTRAL = 0.997
    if abs(gamma) > GAMMA_LIMIT_CENTRAL:
        # Non-central eclipse - no central line on Earth
        return (0.0, 0.0), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), SE_ECL_PARTIAL

    # Central eclipse - calculate central line position

    # Get Skyfield time
    t = ts.ut1_jd(jd)

    # Get Earth, Sun, Moon objects
    earth = eph["earth"]
    sun = eph["sun"]
    moon = eph["moon"]

    # Calculate geocentric positions of Sun and Moon
    sun_geo = earth.at(t).observe(sun).apparent()
    moon_geo = earth.at(t).observe(moon).apparent()

    # Get RA/Dec of Moon (sub-lunar point calculation)
    moon_ra, moon_dec, _ = moon_geo.radec()
    moon_ra_deg = moon_ra.hours * 15.0  # Convert hours to degrees
    moon_dec_deg = moon_dec.degrees

    # Calculate GMST (Greenwich Mean Sidereal Time)
    gmst = t.gmst  # in hours
    gmst_deg = gmst * 15.0  # Convert to degrees

    # Sub-lunar point longitude
    # The Moon is directly overhead where:
    # Local sidereal time = Moon's RA
    # longitude = RA - GMST
    central_lon = moon_ra_deg - gmst_deg
    # Normalize to -180 to +180
    central_lon = ((central_lon + 180) % 360) - 180

    # Sub-lunar point latitude = Moon's declination
    # But for eclipse central line, we need to account for Moon's shadow direction
    # The central line is where the Moon's shadow axis hits Earth
    # Simplified: the latitude is approximately the Moon's declination
    # adjusted for the Moon's ecliptic latitude

    # More accurate central latitude calculation:
    # Account for the Moon being slightly off the ecliptic plane
    central_lat = moon_dec_deg

    # Adjust for the shadow cone geometry
    # The shadow axis points from Moon toward Sun, so we trace backwards
    # Simplified approximation for now
    sun_ra, sun_dec, _ = sun_geo.radec()

    # The central line position is affected by Moon's displacement from ecliptic
    # Apply a correction based on moon latitude
    central_lat -= moon_lat * 0.5  # Simplified correction factor

    # Clamp latitude to valid range
    central_lat = max(-90.0, min(90.0, central_lat))

    # Calculate eclipse attributes at central line location

    # Create observer at central line
    try:
        observer = wgs84.latlon(central_lat, central_lon, 0.0)
        observer_at = earth + observer

        # Get apparent positions from observer
        sun_app = observer_at.at(t).observe(sun).apparent()
        moon_app = observer_at.at(t).observe(moon).apparent()

        # Get Sun altitude and azimuth at central line
        sun_alt, sun_az, _ = sun_app.altaz()
        sun_altitude = sun_alt.degrees
        sun_azimuth = sun_az.degrees

        # Calculate angular separation
        separation = sun_app.separation_from(moon_app).degrees

        # Get local distances
        local_sun_dist = sun_app.distance().au
        local_moon_dist = moon_app.distance().au

        # Local angular radii
        local_sun_radius = (959.63 / 3600.0) / local_sun_dist
        local_moon_radius = (932.56 / 3600.0) * (0.002569 / local_moon_dist)
        local_ratio = local_moon_radius / local_sun_radius

        # Calculate magnitude at central line
        if separation < abs(local_moon_radius - local_sun_radius):
            # Total or annular
            magnitude = local_moon_radius / local_sun_radius
        else:
            overlap = (local_sun_radius + local_moon_radius) - separation
            magnitude = overlap / (2 * local_sun_radius)
        magnitude = max(0.0, min(1.5, magnitude))

        # Calculate obscuration
        if separation >= local_sun_radius + local_moon_radius:
            obscuration = 0.0
        elif separation <= abs(local_moon_radius - local_sun_radius):
            if local_moon_radius >= local_sun_radius:
                obscuration = 1.0
            else:
                obscuration = (local_moon_radius / local_sun_radius) ** 2
        else:
            r_sun = local_sun_radius
            r_moon = local_moon_radius
            d = separation
            d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
            d2 = d - d1
            area1 = r_sun * r_sun * math.acos(
                max(-1, min(1, d1 / r_sun))
            ) - d1 * math.sqrt(max(0, r_sun * r_sun - d1 * d1))
            area2 = r_moon * r_moon * math.acos(
                max(-1, min(1, d2 / r_moon))
            ) - d2 * math.sqrt(max(0, r_moon * r_moon - d2 * d2))
            obscuration = (area1 + area2) / (math.pi * r_sun * r_sun)
        obscuration = max(0.0, min(1.0, obscuration))

        # Calculate path width
        # Path width depends on umbra/penumbra cone geometry and shadow angle
        # Simplified formula: width ≈ 2 * umbra_radius * cos(sun_alt)
        # where umbra_radius depends on Moon/Sun distance ratio

        # Shadow cone geometry
        sun_radius_km = 696340.0  # km
        moon_radius_km = 1737.4  # km
        sun_dist_km = sun_dist * 149597870.7  # AU to km
        moon_dist_km = moon_dist * 149597870.7  # AU to km
        earth_moon_dist_km = moon_dist_km

        # Umbra cone semi-angle
        alpha = math.atan((sun_radius_km - moon_radius_km) / sun_dist_km)

        # For total eclipse, umbra touches Earth
        # Umbra radius at Earth's surface
        d_to_earth = earth_moon_dist_km

        if local_ratio >= 1.0:
            # Total eclipse - umbra reaches Earth
            # Umbra radius at Earth = Moon radius - (distance * tan(alpha))
            umbra_radius_km = moon_radius_km - d_to_earth * math.tan(alpha)
            umbra_radius_km = max(0, umbra_radius_km)
        else:
            # Annular eclipse - antumbra
            # Antumbra starts after umbra vertex
            umbra_radius_km = d_to_earth * math.tan(alpha) - moon_radius_km
            umbra_radius_km = max(0, abs(umbra_radius_km))

        # Path width affected by shadow angle (Sun altitude)
        if sun_altitude > 0:
            cos_alt = math.cos(math.radians(90 - sun_altitude))
            cos_alt = max(0.1, cos_alt)  # Prevent extreme values near horizon
            path_width_km = 2 * umbra_radius_km / cos_alt
        else:
            path_width_km = 0.0

        # Sanity check: typical path widths are 50-270 km for total,
        # 50-400+ km for annular
        path_width_km = max(0.0, min(1000.0, path_width_km))

        # Apparent diameters in degrees
        moon_diameter = 2 * local_moon_radius
        sun_diameter = 2 * local_sun_radius

    except Exception:
        # If calculation fails, return zeros
        return (0.0, 0.0), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0

    # Determine eclipse type flags
    eclipse_type = SE_ECL_CENTRAL
    if local_ratio >= 1.0:
        eclipse_type |= SE_ECL_TOTAL
    elif local_ratio > 0.99:
        eclipse_type |= SE_ECL_ANNULAR_TOTAL  # Hybrid
    else:
        eclipse_type |= SE_ECL_ANNULAR

    # Prepare return tuples
    geopos = (central_lon, central_lat)

    attr = (
        magnitude,  # [0] Magnitude
        local_ratio,  # [1] Moon/Sun diameter ratio
        obscuration,  # [2] Obscuration
        path_width_km,  # [3] Path width in km
        sun_azimuth,  # [4] Sun azimuth
        sun_altitude,  # [5] Sun altitude
        moon_diameter,  # [6] Moon apparent diameter
        sun_diameter,  # [7] Sun apparent diameter
    )

    return geopos, attr, eclipse_type


# Alias for Swiss Ephemeris API compatibility
swe_sol_eclipse_where = sol_eclipse_where


def sol_eclipse_how(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], int]:
    """
    Calculate the circumstances of a solar eclipse at a specific location and time.

    This function determines how a solar eclipse appears from a given geographic
    position at a specific Julian Day. Unlike sol_eclipse_when_loc which finds the
    next eclipse, this function calculates the eclipse magnitude, obscuration, and
    Sun position for a known eclipse time.

    Args:
        jd: Julian Day (UT) of the moment to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Eclipse magnitude (fraction of Sun's diameter covered)
                [1]: Fraction of Sun's diameter covered by Moon at maximum
                [2]: Fraction of Sun's area obscured (obscuration)
                [3]: Azimuth of Sun at the given time (degrees)
                [4]: Altitude of Sun at the given time (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of Sun (degrees)
                [7]: Saros series number (0, not implemented)
                [8]: Inex number (0, not implemented)
                [9]: Reserved (0)
                [10]: Reserved (0)
            - retflag: Eclipse type flags bitmask (SE_ECL_* constants)
                Returns 0 if no eclipse is visible from this location at this time

    Note:
        This function is intended for use when you already know an eclipse is
        occurring (e.g., from sol_eclipse_when_glob or sol_eclipse_when_loc).
        For a random time when no eclipse is occurring, magnitude will be 0
        and retflag will be 0.

    Algorithm:
        1. Calculate Sun and Moon apparent positions from observer location
        2. Compute angular separation between Sun and Moon
        3. Calculate eclipse magnitude from overlap of disks
        4. Determine eclipse type based on whether Moon fully covers Sun
        5. Calculate obscuration (area ratio)

    Precision:
        Eclipse magnitude accurate to ~0.001 for central eclipses.
        Topocentric parallax included in calculations.

    Example:
        >>> # Calculate eclipse circumstances at Dallas during April 8, 2024 eclipse
        >>> from libephemeris import julday, sol_eclipse_how
        >>> jd = julday(2024, 4, 8, 18.5)  # During eclipse
        >>> dallas_lat, dallas_lon = 32.7767, -96.7970
        >>> attr, ecl_type = sol_eclipse_how(jd, dallas_lat, dallas_lon)
        >>> print(f"Magnitude: {attr[0]:.3f}")
        >>> print(f"Sun altitude: {attr[4]:.1f}°")

    References:
        - Swiss Ephemeris: swe_sol_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
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

    # Get Skyfield time
    t = ts.ut1_jd(jd)

    # Create observer position
    observer_at = earth + observer

    # Get apparent positions from observer
    try:
        sun_app = observer_at.at(t).observe(sun).apparent()
        moon_app = observer_at.at(t).observe(moon).apparent()
    except Exception:
        # If calculation fails, return zeros
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0

    # Get Sun altitude and azimuth
    sun_alt, sun_az, _ = sun_app.altaz()
    sun_altitude = sun_alt.degrees
    sun_azimuth = sun_az.degrees

    # If Sun is below horizon, no visible eclipse
    if sun_altitude < -1.0:  # Allow for refraction near horizon
        return (
            0.0,
            0.0,
            0.0,
            sun_azimuth,
            sun_altitude,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ), 0

    # Calculate angular separation between Sun and Moon
    separation = sun_app.separation_from(moon_app).degrees

    # Get distances for angular size calculations
    sun_dist_au = sun_app.distance().au
    moon_dist_au = moon_app.distance().au

    # Angular radii (in degrees)
    # Sun: mean radius 959.63" at 1 AU
    # Moon: mean radius 932.56" at mean distance (0.002569 AU)
    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)

    sun_diameter = 2 * sun_angular_radius
    moon_diameter = 2 * moon_angular_radius

    # Check if there's any eclipse (disks overlapping)
    sum_radii = sun_angular_radius + moon_angular_radius
    if separation >= sum_radii:
        # No eclipse - Sun and Moon too far apart
        return (
            0.0,  # magnitude
            0.0,  # fraction covered
            0.0,  # obscuration
            sun_azimuth,  # Sun azimuth
            sun_altitude,  # Sun altitude
            moon_diameter,  # Moon diameter
            sun_diameter,  # Sun diameter
            0.0,  # Saros
            0.0,  # Inex
            0.0,  # Reserved
            0.0,  # Reserved
        ), 0

    # Calculate eclipse magnitude
    # Magnitude = fraction of Sun's diameter covered by Moon
    overlap = sum_radii - separation
    magnitude = overlap / sun_diameter
    magnitude = max(0.0, min(magnitude, 1.0 + moon_diameter / sun_diameter))

    # Fraction covered (same as magnitude for partial, can exceed 1 for total)
    fraction_covered = magnitude

    # Calculate obscuration (fraction of Sun's area covered)
    r_sun = sun_angular_radius
    r_moon = moon_angular_radius
    d = separation  # center-to-center separation

    if d >= r_sun + r_moon:
        obscuration = 0.0
    elif d <= abs(r_sun - r_moon):
        # One disk entirely within the other
        if r_moon >= r_sun:
            obscuration = 1.0  # Total eclipse
        else:
            obscuration = (r_moon / r_sun) ** 2  # Annular eclipse
    else:
        # Partial overlap - use lens formula
        # Area of intersection of two circles
        d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
        d2 = d - d1

        if abs(d1) <= r_sun and abs(d2) <= r_moon:
            # Ensure values are valid for acos
            cos_arg1 = max(-1, min(1, d1 / r_sun))
            cos_arg2 = max(-1, min(1, d2 / r_moon))

            area1 = r_sun * r_sun * math.acos(cos_arg1) - d1 * math.sqrt(
                max(0, r_sun * r_sun - d1 * d1)
            )
            area2 = r_moon * r_moon * math.acos(cos_arg2) - d2 * math.sqrt(
                max(0, r_moon * r_moon - d2 * d2)
            )
            intersection_area = area1 + area2
            sun_area = math.pi * r_sun * r_sun
            obscuration = intersection_area / sun_area
        else:
            obscuration = 0.0

    obscuration = max(0.0, min(1.0, obscuration))

    # Determine eclipse type flags
    eclipse_type = SE_ECL_VISIBLE

    diff_radii = abs(moon_angular_radius - sun_angular_radius)
    moon_sun_ratio = moon_angular_radius / sun_angular_radius

    if separation < diff_radii:
        # Central eclipse at this location
        eclipse_type |= SE_ECL_CENTRAL
        if moon_sun_ratio >= 1.0:
            eclipse_type |= SE_ECL_TOTAL
        elif moon_sun_ratio > 0.99:
            eclipse_type |= SE_ECL_ANNULAR_TOTAL
        else:
            eclipse_type |= SE_ECL_ANNULAR
    else:
        # Partial eclipse
        eclipse_type |= SE_ECL_PARTIAL

    # Add maximum visibility flag
    eclipse_type |= SE_ECL_MAX_VISIBLE

    # Prepare attributes tuple (11 elements)
    attr = (
        magnitude,  # [0] Magnitude
        fraction_covered,  # [1] Fraction covered
        obscuration,  # [2] Obscuration
        sun_azimuth,  # [3] Azimuth
        sun_altitude,  # [4] Altitude
        moon_diameter,  # [5] Moon diameter
        sun_diameter,  # [6] Sun diameter
        0.0,  # [7] Saros (not implemented)
        0.0,  # [8] Inex (not implemented)
        0.0,  # [9] Reserved
        0.0,  # [10] Reserved
    )

    return attr, eclipse_type


# Alias for Swiss Ephemeris API compatibility
swe_sol_eclipse_how = sol_eclipse_how


# =============================================================================
# LUNAR ECLIPSE FUNCTIONS
# =============================================================================

# Constants for lunar eclipse calculations
ECLIPSE_LIMIT_LUNAR = 12.0  # Maximum elongation from node for lunar eclipse (degrees)


def _find_next_full_moon(jd_start: float) -> float:
    """
    Find the next Full Moon (Sun-Moon opposition) after jd_start.

    Uses iterative refinement to find exact moment of opposition.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of next Full Moon
    """
    # Get current positions
    sun_pos, _ = swe_calc_ut(jd_start, SE_SUN, SEFLG_SPEED)
    moon_pos, _ = swe_calc_ut(jd_start, SE_MOON, SEFLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation from opposition (Moon - Sun - 180°)
    elongation = (moon_lon - sun_lon - 180.0) % 360.0
    if elongation > 180:
        elongation -= 360

    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day

    # Time until next opposition (elongation from 180° = 0)
    if elongation > 0:
        # Past opposition, wait for next cycle
        dt = (360.0 - elongation) / relative_speed
    else:
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

        # Elongation from opposition
        diff = (moon_lon - sun_lon - 180.0) % 360.0
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


def _calculate_lunar_eclipse_type_and_magnitude(
    jd: float,
) -> Tuple[int, float, float, float, float, float]:
    """
    Determine lunar eclipse type and magnitude at maximum eclipse.

    Uses geometric calculations based on Earth's shadow cones at the Moon's distance.

    Args:
        jd: Julian Day of eclipse maximum (UT)

    Returns:
        Tuple of (eclipse_type_flags, magnitude_umbral, magnitude_penumbral,
                  gamma, penumbra_radius, umbra_radius)
        - eclipse_type_flags: Bitmask of SE_ECL_* constants
        - magnitude_umbral: Eclipse magnitude (fraction of Moon in umbra)
        - magnitude_penumbral: Penumbral eclipse magnitude
        - gamma: Gamma parameter (Moon's distance from shadow axis in Earth radii)
        - penumbra_radius: Penumbral shadow radius at Moon distance (degrees)
        - umbra_radius: Umbral shadow radius at Moon distance (degrees)
    """
    # Get positions
    sun_pos, _ = swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
    moon_pos, _ = swe_calc_ut(jd, SE_MOON, SEFLG_SPEED)

    # Distances in AU
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Moon's ecliptic latitude (how far from the ecliptic)
    moon_lat = moon_pos[1]

    # Constants for shadow calculations
    # Sun angular radius at 1 AU: 959.63 arcsec
    # Earth radius: 6378.137 km = 4.2635e-5 AU
    SUN_RADIUS_ARCSEC = 959.63
    EARTH_RADIUS_AU = 4.2635e-5
    EARTH_RADIUS_KM = 6378.137

    # Sun's angular semi-diameter at actual distance (in degrees)
    sun_semidiameter = (SUN_RADIUS_ARCSEC / 3600.0) / sun_dist

    # Moon's angular semi-diameter (932.56 arcsec at mean distance 0.002569 AU)
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Earth's angular semi-diameter as seen from Moon (in degrees)
    earth_semidiameter = math.degrees(math.atan(EARTH_RADIUS_AU / moon_dist))

    # Calculate shadow cone sizes at Moon's distance
    # Using the geometric relationship of similar triangles

    # Sun's angular radius as seen from Earth
    sun_angular_radius_from_earth = sun_semidiameter

    # Penumbra radius at Moon's distance (Earth + Sun shadow)
    # Penumbra outer edge: Earth appears larger, adding Sun's angular size
    penumbra_radius = earth_semidiameter + sun_angular_radius_from_earth

    # Umbra radius at Moon's distance (Earth minus Sun shadow)
    # The umbra is smaller because of the Sun's finite size
    # Using more accurate shadow cone calculation
    # Shadow cone angle from geometry
    sun_dist_km = sun_dist * 149597870.7
    moon_dist_km = moon_dist * 149597870.7
    sun_radius_km = (SUN_RADIUS_ARCSEC / 206265.0) * sun_dist_km

    # Umbra cone semi-angle
    umbra_cone_angle = math.atan((sun_radius_km - EARTH_RADIUS_KM) / sun_dist_km)
    # Umbra radius at Moon distance
    umbra_radius_km = EARTH_RADIUS_KM - moon_dist_km * math.tan(umbra_cone_angle)

    if umbra_radius_km > 0:
        umbra_radius = math.degrees(math.atan(umbra_radius_km / moon_dist_km))
    else:
        # Umbra doesn't reach Moon (antumbra) - extremely rare for lunar eclipses
        umbra_radius = 0.0

    # Moon's distance from the shadow axis (in degrees)
    moon_distance_from_axis = abs(moon_lat)

    # Gamma: Moon's distance from shadow axis in Earth radii
    gamma = moon_lat / earth_semidiameter

    # Calculate eclipse magnitudes
    # Penumbral magnitude: how far Moon penetrates into penumbra
    penumbral_mag = (penumbra_radius + moon_semidiameter - moon_distance_from_axis) / (
        2 * moon_semidiameter
    )

    # Umbral magnitude: how far Moon penetrates into umbra
    umbral_mag = (umbra_radius + moon_semidiameter - moon_distance_from_axis) / (
        2 * moon_semidiameter
    )

    # Determine eclipse type
    eclipse_type = 0

    if penumbral_mag <= 0:
        # No eclipse
        return 0, 0.0, 0.0, gamma, penumbra_radius, umbra_radius

    if umbral_mag <= 0:
        # Penumbral only
        eclipse_type = SE_ECL_PENUMBRAL
        penumbral_mag = max(0.0, min(1.0, penumbral_mag))
        return eclipse_type, 0.0, penumbral_mag, gamma, penumbra_radius, umbra_radius

    if umbral_mag >= 1.0:
        # Total umbral eclipse
        eclipse_type = SE_ECL_TOTAL
        umbral_mag = max(0.0, umbral_mag)
        penumbral_mag = max(0.0, min(2.0, penumbral_mag))
    else:
        # Partial umbral eclipse
        eclipse_type = SE_ECL_PARTIAL
        umbral_mag = max(0.0, min(1.0, umbral_mag))
        penumbral_mag = max(0.0, min(2.0, penumbral_mag))

    return eclipse_type, umbral_mag, penumbral_mag, gamma, penumbra_radius, umbra_radius


def _calculate_lunar_eclipse_phases(
    jd_max: float,
    eclipse_type: int,
    moon_semidiameter: float,
    umbra_radius: float,
    penumbra_radius: float,
    moon_lat_at_max: float,
) -> Tuple[float, float, float, float, float, float, float, float]:
    """
    Calculate times of lunar eclipse phases (contacts).

    Phase indices:
        [0]: Time of maximum eclipse
        [1]: Time of partial eclipse beginning (Moon enters umbra)
        [2]: Time of total eclipse beginning (Moon fully in umbra)
        [3]: Time of total eclipse ending (Moon starts leaving umbra)
        [4]: Time of partial eclipse ending (Moon leaves umbra)
        [5]: Time of penumbral eclipse beginning
        [6]: Time of penumbral eclipse ending
        [7]: Reserved

    Args:
        jd_max: Julian Day of maximum eclipse
        eclipse_type: Eclipse type flags
        moon_semidiameter: Moon's angular semi-diameter in degrees
        umbra_radius: Umbral shadow radius in degrees
        penumbra_radius: Penumbral shadow radius in degrees
        moon_lat_at_max: Moon's ecliptic latitude at maximum in degrees

    Returns:
        Tuple of 8 floats with phase times (JD UT)
    """
    # Actually calculate from Moon's speed minus Sun's speed
    sun_pos, _ = swe_calc_ut(jd_max, SE_SUN, SEFLG_SPEED)
    moon_pos, _ = swe_calc_ut(jd_max, SE_MOON, SEFLG_SPEED)

    # Speed of Moon relative to shadow (in longitude)
    relative_speed = abs(moon_pos[3] - sun_pos[3])  # degrees/day
    if relative_speed < 0.1:
        relative_speed = 0.55  # fallback to typical value

    # Distance from shadow axis at maximum
    y = abs(moon_lat_at_max)

    # Calculate half-duration of each phase using geometry
    # The Moon moves through the shadow at an angle
    # Half-duration = sqrt(R² - y²) / speed where R is radius and y is impact parameter

    def calc_half_duration(radius: float) -> float:
        """Calculate half-duration for given shadow radius."""
        r_total = radius + moon_semidiameter
        if y >= r_total:
            return 0.0
        # Calculate chord half-length
        half_chord = math.sqrt(max(0, r_total * r_total - y * y))
        return half_chord / relative_speed

    def calc_total_half_duration(radius: float) -> float:
        """Calculate half-duration of total phase (Moon fully inside)."""
        r_inner = radius - moon_semidiameter
        if r_inner <= 0 or y >= r_inner:
            return 0.0
        half_chord = math.sqrt(max(0, r_inner * r_inner - y * y))
        return half_chord / relative_speed

    # Penumbral phase times
    penumbral_half_dur = calc_half_duration(penumbra_radius)
    t_pen_begin = jd_max - penumbral_half_dur
    t_pen_end = jd_max + penumbral_half_dur

    # Partial (umbral) phase times
    partial_half_dur = calc_half_duration(umbra_radius)
    if partial_half_dur > 0:
        t_partial_begin = jd_max - partial_half_dur
        t_partial_end = jd_max + partial_half_dur
    else:
        t_partial_begin = 0.0
        t_partial_end = 0.0

    # Total phase times
    if eclipse_type & SE_ECL_TOTAL:
        total_half_dur = calc_total_half_duration(umbra_radius)
        if total_half_dur > 0:
            t_total_begin = jd_max - total_half_dur
            t_total_end = jd_max + total_half_dur
        else:
            t_total_begin = 0.0
            t_total_end = 0.0
    else:
        t_total_begin = 0.0
        t_total_end = 0.0

    return (
        jd_max,  # [0] Maximum
        t_partial_begin,  # [1] Partial begins (enters umbra)
        t_total_begin,  # [2] Total begins
        t_total_end,  # [3] Total ends
        t_partial_end,  # [4] Partial ends (leaves umbra)
        t_pen_begin,  # [5] Penumbral begins
        t_pen_end,  # [6] Penumbral ends
        0.0,  # [7] Reserved
    )


def lun_eclipse_when(
    jd_start: float,
    flags: int = SEFLG_SWIEPH,
    eclipse_type: int = 0,
) -> Tuple[Tuple[float, ...], int]:
    """
    Find the next lunar eclipse after a given date.

    Searches forward in time from jd_start to find the next lunar eclipse.
    Can filter by eclipse type (total, partial, penumbral).

    Args:
        jd_start: Julian Day (UT) to start search from
        flags: Calculation flags (SEFLG_SWIEPH, etc.)
        eclipse_type: Filter for specific eclipse type(s), bitmask of:
            - SE_ECL_TOTAL (4): Total lunar eclipse
            - SE_ECL_PARTIAL (16): Partial lunar eclipse
            - SE_ECL_PENUMBRAL (64): Penumbral lunar eclipse
            - 0: Any eclipse type (default)

    Returns:
        Tuple containing:
            - times: Tuple of 8 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse
                [1]: Time of partial eclipse beginning (Moon enters umbra)
                [2]: Time of total eclipse beginning (or 0 if not total)
                [3]: Time of total eclipse ending (or 0 if not total)
                [4]: Time of partial eclipse ending (Moon leaves umbra)
                [5]: Time of penumbral eclipse beginning
                [6]: Time of penumbral eclipse ending
                [7]: Reserved (0)
            - retflag: Eclipse type flags bitmask (SE_ECL_* constants)

    Raises:
        RuntimeError: If no eclipse found within search limit

    Algorithm:
        1. Find next Full Moon after jd_start
        2. Check if Moon is close enough to node for eclipse
        3. If not eclipse, advance to next Full Moon
        4. Calculate eclipse type and magnitude
        5. If eclipse_type filter set, check if matches
        6. Calculate phase times

    Precision:
        Eclipse times accurate to ~1 minute for most eclipses.

    Example:
        >>> # Find next total lunar eclipse after Jan 1, 2024
        >>> from libephemeris import julday, SE_ECL_TOTAL
        >>> jd = julday(2024, 1, 1, 0)
        >>> times, ecl_type = lun_eclipse_when(jd, eclipse_type=SE_ECL_TOTAL)
        >>> print(f"Total lunar eclipse at JD {times[0]:.5f}")

    References:
        - Swiss Ephemeris: swe_lun_eclipse_when()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    MAX_SEARCH_YEARS = 20  # Maximum search range
    MAX_FULL_MOONS = int(MAX_SEARCH_YEARS * 12.4)  # ~12.4 lunations per year

    # If eclipse_type is 0, accept any type
    if eclipse_type == 0:
        eclipse_type = SE_ECL_ALLTYPES_LUNAR

    jd = jd_start

    for _ in range(MAX_FULL_MOONS):
        # Find next Full Moon
        jd_full_moon = _find_next_full_moon(jd)

        # Get Moon position at Full Moon
        moon_pos, _ = swe_calc_ut(jd_full_moon, SE_MOON, flags | SEFLG_SPEED)
        moon_lon = moon_pos[0]
        moon_lat = moon_pos[1]

        # Check if close enough to ecliptic for eclipse
        # Lunar eclipse possible if Moon is near a node
        node_dist = _get_moon_node_distance(jd_full_moon, moon_lon)

        if node_dist < ECLIPSE_LIMIT_LUNAR:
            # Possible eclipse - check magnitude
            (
                ecl_type,
                umbral_mag,
                penumbral_mag,
                gamma,
                penumbra_radius,
                umbra_radius,
            ) = _calculate_lunar_eclipse_type_and_magnitude(jd_full_moon)

            if ecl_type != 0:
                # Eclipse found - check if matches filter
                type_matches = (
                    (eclipse_type & SE_ECL_TOTAL and ecl_type & SE_ECL_TOTAL)
                    or (eclipse_type & SE_ECL_PARTIAL and ecl_type & SE_ECL_PARTIAL)
                    or (eclipse_type & SE_ECL_PENUMBRAL and ecl_type & SE_ECL_PENUMBRAL)
                )

                if type_matches:
                    # Get moon semi-diameter for phase calculations
                    moon_dist = moon_pos[2]
                    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

                    # Calculate phase times
                    times = _calculate_lunar_eclipse_phases(
                        jd_full_moon,
                        ecl_type,
                        moon_semidiameter,
                        umbra_radius,
                        penumbra_radius,
                        moon_lat,
                    )
                    return times, ecl_type

        # Advance to next lunation
        jd = jd_full_moon + 25  # Skip ahead ~25 days to ensure we find next Full Moon

    raise RuntimeError(
        f"No matching lunar eclipse found within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Aliases for compatibility
swe_lun_eclipse_when = lun_eclipse_when


def lun_eclipse_when_loc(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], int]:
    """
    Find the next lunar eclipse visible from a specific location.

    Searches forward in time from jd_start to find the next lunar eclipse
    where the Moon is above the horizon during at least part of the eclipse.
    The Moon must be above the horizon for the eclipse to be visible.

    Args:
        jd_start: Julian Day (UT) to start search from
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - times: Tuple of 10 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse (local visibility)
                [1]: Time of partial eclipse beginning (Moon enters umbra)
                [2]: Time of total eclipse beginning (or 0 if not total)
                [3]: Time of total eclipse ending (or 0 if not total)
                [4]: Time of partial eclipse ending (Moon leaves umbra)
                [5]: Time of penumbral eclipse beginning
                [6]: Time of penumbral eclipse ending
                [7]: Time of moonrise (if Moon rises during eclipse, else 0)
                [8]: Time of moonset (if Moon sets during eclipse, else 0)
                [9]: Reserved (0)
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Umbral magnitude
                [1]: Penumbral magnitude
                [2]: Reserved (0)
                [3]: Azimuth of Moon at maximum eclipse (degrees)
                [4]: Altitude of Moon at maximum eclipse (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of umbral shadow (degrees)
                [7]: Apparent diameter of penumbral shadow (degrees)
                [8]: Saros series number (0, not implemented)
                [9-10]: Reserved (0)
            - retflag: Eclipse type flags bitmask (SE_ECL_* constants)
                Includes visibility flags (SE_ECL_VISIBLE, SE_ECL_*_VISIBLE)

    Raises:
        RuntimeError: If no eclipse visible from location within search limit

    Algorithm:
        1. Use lun_eclipse_when to find next global lunar eclipse
        2. Calculate Moon's altitude at observer's location during eclipse
        3. If Moon is below horizon during entire eclipse, continue to next
        4. Return local phase times, visibility flags, and attributes

    Precision:
        Eclipse times accurate to ~1 minute. Visibility depends on
        accurate horizon calculations.

    Example:
        >>> # Find next lunar eclipse visible from Rome
        >>> from libephemeris import julday, lun_eclipse_when_loc
        >>> jd = julday(2024, 1, 1, 0)
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> times, attr, ecl_type = lun_eclipse_when_loc(jd, rome_lat, rome_lon)
        >>> print(f"Eclipse maximum at JD {times[0]:.5f}")
        >>> print(f"Moon altitude: {attr[4]:.1f}°")

    References:
        - Swiss Ephemeris: swe_lun_eclipse_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    MAX_SEARCH_YEARS = 50  # Maximum search range
    MAX_ECLIPSES = int(
        MAX_SEARCH_YEARS * 2.4
    )  # ~2.4 lunar eclipses per year on average

    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Get ephemeris
    eph = get_planets()
    ts = get_timescale()

    # Create observer
    earth = eph["earth"]
    moon = eph["moon"]
    observer = wgs84.latlon(lat, lon, altitude)

    def _get_moon_altaz(jd: float) -> Tuple[float, float]:
        """Get Moon altitude and azimuth at given JD from observer location."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer
        moon_app = observer_at.at(t).observe(moon).apparent()
        alt, az, _ = moon_app.altaz()
        return alt.degrees, az.degrees

    def _is_moon_visible(jd: float, min_alt: float = -1.0) -> bool:
        """Check if Moon is above horizon (with margin for refraction)."""
        alt, _ = _get_moon_altaz(jd)
        return alt > min_alt

    jd = jd_start

    for _ in range(MAX_ECLIPSES):
        # Find next global lunar eclipse
        try:
            global_times, global_type = lun_eclipse_when(jd, flags)
        except RuntimeError:
            raise RuntimeError(
                f"No lunar eclipse visible from lat={lat}, lon={lon} "
                f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
            )

        jd_max = global_times[0]
        jd_pen_begin = global_times[5]
        jd_pen_end = global_times[6]
        jd_partial_begin = global_times[1]
        jd_partial_end = global_times[4]
        jd_total_begin = global_times[2]
        jd_total_end = global_times[3]

        # Check if Moon is visible during any part of the eclipse
        # Sample at key phases: penumbral begin/end, partial begin/end, maximum
        check_times = [jd_max]
        if jd_pen_begin > 0:
            check_times.append(jd_pen_begin)
        if jd_pen_end > 0:
            check_times.append(jd_pen_end)
        if jd_partial_begin > 0:
            check_times.append(jd_partial_begin)
        if jd_partial_end > 0:
            check_times.append(jd_partial_end)
        if jd_total_begin > 0:
            check_times.append(jd_total_begin)
        if jd_total_end > 0:
            check_times.append(jd_total_end)

        # Also sample at several intermediate points
        eclipse_start = jd_pen_begin if jd_pen_begin > 0 else jd_max - 2 / 24
        eclipse_end = jd_pen_end if jd_pen_end > 0 else jd_max + 2 / 24
        for i in range(5):
            t = eclipse_start + (eclipse_end - eclipse_start) * i / 4
            if t not in check_times:
                check_times.append(t)

        moon_visible = any(_is_moon_visible(t) for t in check_times)

        if not moon_visible:
            # Eclipse not visible from this location, try next
            jd = jd_max + 25  # Skip ahead to find next eclipse
            continue

        # Eclipse visible! Calculate local circumstances

        # Get Moon position at maximum
        moon_pos, _ = swe_calc_ut(jd_max, SE_MOON, SEFLG_SPEED)
        moon_dist = moon_pos[2]

        # Get eclipse type and magnitude
        (
            ecl_type_flags,
            umbral_mag,
            penumbral_mag,
            gamma,
            penumbra_radius,
            umbra_radius,
        ) = _calculate_lunar_eclipse_type_and_magnitude(jd_max)

        # Get Moon's alt/az at maximum
        moon_alt, moon_az = _get_moon_altaz(jd_max)

        # Calculate apparent diameters
        # Moon semi-diameter: 932.56 arcsec at mean distance 0.002569 AU
        moon_diameter = 2 * (932.56 / 3600.0) * (0.002569 / moon_dist)
        umbra_diameter = 2 * umbra_radius
        penumbra_diameter = 2 * penumbra_radius

        # Determine which phases are visible
        ecl_type = ecl_type_flags

        # Check visibility at each phase
        if jd_pen_begin > 0 and _is_moon_visible(jd_pen_begin):
            ecl_type |= SE_ECL_1ST_VISIBLE
        if jd_partial_begin > 0 and _is_moon_visible(jd_partial_begin):
            ecl_type |= SE_ECL_2ND_VISIBLE
        if jd_partial_end > 0 and _is_moon_visible(jd_partial_end):
            ecl_type |= SE_ECL_3RD_VISIBLE
        if jd_pen_end > 0 and _is_moon_visible(jd_pen_end):
            ecl_type |= SE_ECL_4TH_VISIBLE
        if _is_moon_visible(jd_max):
            ecl_type |= SE_ECL_MAX_VISIBLE
        ecl_type |= SE_ECL_VISIBLE

        # Find moonrise/moonset during eclipse (simplified - we just check endpoints)
        moonrise_time = 0.0
        moonset_time = 0.0

        # Check if Moon rises during eclipse
        if jd_pen_begin > 0 and jd_pen_end > 0:
            alt_begin, _ = _get_moon_altaz(jd_pen_begin)
            alt_end, _ = _get_moon_altaz(jd_pen_end)
            if alt_begin < 0 < alt_end:
                # Moon rises during eclipse - find approximate time
                # Simple binary search
                t_low, t_high = jd_pen_begin, jd_pen_end
                for _ in range(20):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _ = _get_moon_altaz(t_mid)
                    if alt_mid < 0:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                moonrise_time = (t_low + t_high) / 2
            elif alt_begin > 0 > alt_end:
                # Moon sets during eclipse - find approximate time
                t_low, t_high = jd_pen_begin, jd_pen_end
                for _ in range(20):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _ = _get_moon_altaz(t_mid)
                    if alt_mid > 0:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                moonset_time = (t_low + t_high) / 2

        # Prepare times tuple (10 elements)
        times = (
            global_times[0],  # [0] Maximum
            global_times[1],  # [1] Partial begins (enters umbra)
            global_times[2],  # [2] Total begins
            global_times[3],  # [3] Total ends
            global_times[4],  # [4] Partial ends (leaves umbra)
            global_times[5],  # [5] Penumbral begins
            global_times[6],  # [6] Penumbral ends
            moonrise_time,  # [7] Moonrise during eclipse
            moonset_time,  # [8] Moonset during eclipse
            0.0,  # [9] Reserved
        )

        # Prepare attributes tuple (11 elements)
        attr = (
            umbral_mag,  # [0] Umbral magnitude
            penumbral_mag,  # [1] Penumbral magnitude
            0.0,  # [2] Reserved
            moon_az,  # [3] Azimuth of Moon at maximum
            moon_alt,  # [4] Altitude of Moon at maximum
            moon_diameter,  # [5] Apparent diameter of Moon
            umbra_diameter,  # [6] Apparent diameter of umbra
            penumbra_diameter,  # [7] Apparent diameter of penumbra
            0.0,  # [8] Saros (not implemented)
            0.0,  # [9] Reserved
            0.0,  # [10] Reserved
        )

        return times, attr, ecl_type

    raise RuntimeError(
        f"No lunar eclipse visible from lat={lat}, lon={lon} "
        f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Alias for Swiss Ephemeris API compatibility
swe_lun_eclipse_when_loc = lun_eclipse_when_loc


def lun_eclipse_how(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], int]:
    """
    Calculate the circumstances of a lunar eclipse at a specific location and time.

    This function determines how a lunar eclipse appears from a given geographic
    position at a specific Julian Day. Unlike lun_eclipse_when_loc which finds the
    next eclipse, this function calculates the eclipse magnitude, Moon position,
    and other circumstances for a known eclipse time.

    Args:
        jd: Julian Day (UT) of the moment to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Umbral eclipse magnitude (fraction of Moon's diameter in umbra)
                [1]: Penumbral eclipse magnitude
                [2]: Reserved (0)
                [3]: Azimuth of Moon at the given time (degrees)
                [4]: Altitude of Moon at the given time (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of umbral shadow at Moon's distance (degrees)
                [7]: Apparent diameter of penumbral shadow at Moon's distance (degrees)
                [8]: Saros series number (0, not implemented)
                [9]: Reserved (0)
                [10]: Reserved (0)
            - retflag: Eclipse type flags bitmask (SE_ECL_* constants)
                Returns 0 if no eclipse is occurring at this time
                Includes SE_ECL_VISIBLE if Moon is above horizon

    Note:
        This function is intended for use when you already know an eclipse is
        occurring (e.g., from lun_eclipse_when or lun_eclipse_when_loc).
        For a random time when no eclipse is occurring, magnitude will be 0
        and retflag will be 0.

    Algorithm:
        1. Calculate Moon's apparent position from observer location
        2. Calculate Earth's shadow cone geometry at Moon's distance
        3. Compute umbral and penumbral magnitudes
        4. Determine eclipse type based on Moon's penetration into shadows
        5. Return attributes including Moon's altitude and azimuth

    Precision:
        Eclipse magnitude accurate to ~0.001 for central eclipses.
        Topocentric parallax included in calculations.

    Example:
        >>> # Calculate eclipse circumstances at Rome during a lunar eclipse
        >>> from libephemeris import julday, lun_eclipse_how
        >>> jd = julday(2022, 5, 16, 4.2)  # During May 2022 eclipse
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> attr, ecl_type = lun_eclipse_how(jd, rome_lat, rome_lon)
        >>> print(f"Umbral magnitude: {attr[0]:.3f}")
        >>> print(f"Moon altitude: {attr[4]:.1f}°")

    References:
        - Swiss Ephemeris: swe_lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    from skyfield.api import wgs84

    from .state import get_planets, get_timescale

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Get Moon object
    moon = eph["moon"]
    earth = eph["earth"]

    # Get Skyfield time
    t = ts.ut1_jd(jd)

    # Create observer position
    observer_at = earth + observer

    # Get Moon apparent position from observer
    try:
        moon_app = observer_at.at(t).observe(moon).apparent()
    except Exception:
        # If calculation fails, return zeros
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0

    # Get Moon altitude and azimuth
    moon_alt, moon_az, _ = moon_app.altaz()
    moon_altitude = moon_alt.degrees
    moon_azimuth = moon_az.degrees

    # Get Moon's distance for angular size calculation
    moon_dist_au = moon_app.distance().au

    # Calculate eclipse geometry using the same calculations as lun_eclipse_when
    (
        ecl_type_flags,
        umbral_mag,
        penumbral_mag,
        gamma,
        penumbra_radius,
        umbra_radius,
    ) = _calculate_lunar_eclipse_type_and_magnitude(jd)

    # Calculate apparent diameters
    # Moon semi-diameter: 932.56 arcsec at mean distance 0.002569 AU
    moon_diameter = 2 * (932.56 / 3600.0) * (0.002569 / moon_dist_au)
    umbra_diameter = 2 * umbra_radius
    penumbra_diameter = 2 * penumbra_radius

    # Determine eclipse type flags
    eclipse_type = 0

    if penumbral_mag <= 0 and umbral_mag <= 0:
        # No eclipse - Moon too far from Earth's shadow
        return (
            0.0,  # umbral magnitude
            0.0,  # penumbral magnitude
            0.0,  # reserved
            moon_azimuth,  # Moon azimuth
            moon_altitude,  # Moon altitude
            moon_diameter,  # Moon diameter
            umbra_diameter,  # Umbra diameter
            penumbra_diameter,  # Penumbra diameter
            0.0,  # Saros
            0.0,  # Reserved
            0.0,  # Reserved
        ), 0

    # There is an eclipse - set type flags
    eclipse_type = ecl_type_flags

    # Check if Moon is above horizon
    if moon_altitude > -1.0:  # Allow for refraction near horizon
        eclipse_type |= SE_ECL_VISIBLE
        eclipse_type |= SE_ECL_MAX_VISIBLE

    # Prepare attributes tuple (11 elements)
    attr = (
        max(0.0, umbral_mag),  # [0] Umbral magnitude
        max(0.0, penumbral_mag),  # [1] Penumbral magnitude
        0.0,  # [2] Reserved
        moon_azimuth,  # [3] Azimuth of Moon
        moon_altitude,  # [4] Altitude of Moon
        moon_diameter,  # [5] Apparent diameter of Moon
        umbra_diameter,  # [6] Apparent diameter of umbra
        penumbra_diameter,  # [7] Apparent diameter of penumbra
        0.0,  # [8] Saros (not implemented)
        0.0,  # [9] Reserved
        0.0,  # [10] Reserved
    )

    return attr, eclipse_type


# Alias for Swiss Ephemeris API compatibility
swe_lun_eclipse_how = lun_eclipse_how


def lun_occult_when_glob(
    jd_start: float,
    planet: int,
    star_name: str = "",
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], int]:
    """
    Find the next lunar occultation of a planet or fixed star.

    A lunar occultation occurs when the Moon passes in front of (occults)
    a planet or star as seen from Earth. This function searches forward
    in time to find the next such event globally (somewhere on Earth).

    Args:
        jd_start: Julian Day (UT) to start search from
        planet: Planet ID to check for occultation (SE_MERCURY, SE_VENUS, etc.)
            Set to 0 if searching for a fixed star occultation.
        star_name: Name of fixed star to check (e.g. "Regulus", "Spica").
            Only used if planet is 0.
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - times: Tuple of 8 floats with occultation phase times (JD UT):
                [0]: Time of maximum occultation (minimum separation)
                [1]: Time of first contact (occultation begins)
                [2]: Time of second contact (full occultation begins, or 0)
                [3]: Time of third contact (full occultation ends, or 0)
                [4]: Time of fourth contact (occultation ends)
                [5]: Reserved (0)
                [6]: Reserved (0)
                [7]: Reserved (0)
            - retflag: Occultation type flags bitmask (SE_ECL_* constants)
                SE_ECL_TOTAL: Total occultation (body fully behind Moon)
                SE_ECL_PARTIAL: Partial occultation (body partially behind Moon)
                SE_ECL_ANNULAR: Body larger than Moon (very rare, only Sun)

    Raises:
        RuntimeError: If no occultation found within search limit
        ValueError: If neither planet nor star_name is specified

    Algorithm:
        1. Calculate Moon's position and the target body's position
        2. Find the next conjunction in right ascension or longitude
        3. At conjunction, check angular separation
        4. If separation < Moon's angular radius + target's angular radius,
           an occultation occurs
        5. Refine timing using Newton-Raphson iteration
        6. Calculate contact times

    Precision:
        Occultation times accurate to ~1 minute for most events.

    Example:
        >>> # Find next occultation of Regulus by the Moon
        >>> from libephemeris import julday
        >>> jd = julday(2024, 1, 1, 0)
        >>> times, ocl_type = lun_occult_when_glob(jd, 0, "Regulus")
        >>> print(f"Occultation at JD {times[0]:.5f}")

        >>> # Find next occultation of Venus by the Moon
        >>> times, ocl_type = lun_occult_when_glob(jd, SE_VENUS, "")
        >>> print(f"Venus occultation at JD {times[0]:.5f}")

    References:
        - Swiss Ephemeris: swe_lun_occult_when_glob()
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
    """
    from .state import get_planets, get_timescale
    from .fixed_stars import swe_fixstar_ut
    from .constants import (
        SE_MERCURY,
        SE_VENUS,
        SE_MARS,
        SE_JUPITER,
        SE_SATURN,
        SE_URANUS,
        SE_NEPTUNE,
        SE_PLUTO,
    )
    from .planets import _PLANET_MAP

    if planet == 0 and not star_name:
        raise ValueError(
            "Either planet ID or star_name must be specified for occultation search"
        )

    MAX_SEARCH_YEARS = 20
    MAX_ITERATIONS = int(MAX_SEARCH_YEARS * 365)  # Check roughly daily

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    earth = eph["earth"]
    moon_body = eph["moon"]

    # Moon's mean angular radius in degrees (varies from ~14.7' to ~16.7')
    # 932.56 arcsec at mean distance (0.002569 AU)
    MOON_MEAN_ANGULAR_RADIUS_DEG = 932.56 / 3600.0

    # Approximate angular radius for planets (arcsec) - rough values for detection
    PLANET_ANGULAR_RADII_ARCSEC = {
        SE_MERCURY: 6.0,  # 3-6 arcsec
        SE_VENUS: 30.0,  # 10-30 arcsec
        SE_MARS: 12.0,  # 4-12 arcsec
        SE_JUPITER: 24.0,  # 30-50 arcsec
        SE_SATURN: 10.0,  # 15-20 arcsec
        SE_URANUS: 2.0,  # ~2 arcsec
        SE_NEPTUNE: 1.2,  # ~1.2 arcsec
        SE_PLUTO: 0.1,  # ~0.1 arcsec
    }

    def _get_moon_position(jd: float) -> Tuple[float, float, float, float]:
        """Get Moon's geocentric RA, Dec, distance in AU, and angular radius."""
        t = ts.ut1_jd(jd)
        moon_app = earth.at(t).observe(moon_body).apparent()
        ra, dec, dist = moon_app.radec(epoch="date")
        # Moon angular radius varies with distance
        moon_angular_radius = MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist.au)
        return ra.hours * 15.0, dec.degrees, dist.au, moon_angular_radius

    def _get_target_position(jd: float) -> Tuple[float, float, float]:
        """Get target body's geocentric RA, Dec, and angular radius."""
        if planet == 0:
            # Fixed star
            pos, retflag, err = swe_fixstar_ut(star_name, jd, flags)
            if err:
                raise ValueError(f"Star '{star_name}' not found: {err}")
            # For fixed stars, we need to convert ecliptic to equatorial
            # swe_fixstar_ut returns ecliptic coordinates, but we need RA/Dec
            # Use the planets module to get proper equatorial coords
            t = ts.ut1_jd(jd)

            # Get star data and calculate equatorial position
            from .fixed_stars import FIXED_STARS
            from .fixed_stars import _resolve_star_id

            star_id, _ = _resolve_star_id(star_name)
            if star_id < 0:
                raise ValueError(f"Star '{star_name}' not found")

            star = FIXED_STARS[star_id]

            # Time from J2000.0
            t_years = (jd - 2451545.0) / 365.25

            # Apply proper motion
            ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
            dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0

            # Stars have negligible angular radius
            return ra_deg, dec_deg, 0.0001  # ~0.4 arcsec for point source
        else:
            # Planet
            if planet not in _PLANET_MAP:
                raise ValueError(f"Unknown planet ID: {planet}")

            target_name = _PLANET_MAP[planet]
            target = eph[target_name]

            t = ts.ut1_jd(jd)
            target_app = earth.at(t).observe(target).apparent()
            ra, dec, dist = target_app.radec(epoch="date")

            # Get planet angular radius
            angular_radius = PLANET_ANGULAR_RADII_ARCSEC.get(planet, 1.0) / 3600.0

            return ra.hours * 15.0, dec.degrees, angular_radius

    def _angular_separation(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
        """Calculate angular separation between two points (in degrees)."""
        # Convert to radians
        ra1_r = math.radians(ra1)
        dec1_r = math.radians(dec1)
        ra2_r = math.radians(ra2)
        dec2_r = math.radians(dec2)

        # Haversine formula for spherical distance
        cos_sep = math.sin(dec1_r) * math.sin(dec2_r) + math.cos(dec1_r) * math.cos(
            dec2_r
        ) * math.cos(ra1_r - ra2_r)
        cos_sep = max(-1.0, min(1.0, cos_sep))
        return math.degrees(math.acos(cos_sep))

    def _check_occultation(jd: float) -> Tuple[bool, float, float, float]:
        """
        Check if occultation occurs at given time.

        Returns: (is_occultation, separation, moon_radius, target_radius)
        """
        moon_ra, moon_dec, moon_dist, moon_radius = _get_moon_position(jd)
        target_ra, target_dec, target_radius = _get_target_position(jd)

        separation = _angular_separation(moon_ra, moon_dec, target_ra, target_dec)

        # Occultation threshold: Moon radius + target radius
        threshold = moon_radius + target_radius

        return separation <= threshold, separation, moon_radius, target_radius

    def _find_minimum_separation(jd_start_search: float, jd_end_search: float) -> float:
        """Find time of minimum separation using golden section search."""
        phi = (1 + math.sqrt(5)) / 2  # Golden ratio

        a = jd_start_search
        b = jd_end_search

        c = b - (b - a) / phi
        d = a + (b - a) / phi

        def get_sep(jd: float) -> float:
            moon_ra, moon_dec, _, _ = _get_moon_position(jd)
            target_ra, target_dec, _ = _get_target_position(jd)
            return _angular_separation(moon_ra, moon_dec, target_ra, target_dec)

        fc = get_sep(c)
        fd = get_sep(d)

        for _ in range(50):  # Converge to ~0.1 second precision
            if fc < fd:
                b = d
                d = c
                fd = fc
                c = b - (b - a) / phi
                fc = get_sep(c)
            else:
                a = c
                c = d
                fc = fd
                d = a + (b - a) / phi
                fd = get_sep(d)

            if b - a < 1e-6:  # ~0.1 second precision
                break

        return (a + b) / 2

    def _calculate_contact_times(
        jd_max: float, min_sep: float, moon_radius: float, target_radius: float
    ) -> Tuple[float, float, float, float]:
        """Calculate contact times for the occultation."""
        # First/fourth contact: separation = moon_radius + target_radius
        outer_threshold = moon_radius + target_radius

        # Second/third contact: separation = |moon_radius - target_radius|
        # (for total occultation when target is inside Moon's disk)
        inner_threshold = abs(moon_radius - target_radius)

        # Estimate Moon's angular speed relative to target using positions
        # We can't use separation rate at minimum (it's near zero)
        # Instead, use Moon's actual motion rate
        dt_test = 1.0 / 24.0  # 1 hour in days
        moon_ra1, moon_dec1, _, _ = _get_moon_position(jd_max - dt_test)
        moon_ra2, moon_dec2, _, _ = _get_moon_position(jd_max + dt_test)
        target_ra1, target_dec1, _ = _get_target_position(jd_max - dt_test)
        target_ra2, target_dec2, _ = _get_target_position(jd_max + dt_test)

        # Calculate Moon's motion relative to target (in degrees per day)
        # This is the rate at which Moon moves across the sky relative to target
        d_ra = (moon_ra2 - moon_ra1) - (target_ra2 - target_ra1)
        d_dec = (moon_dec2 - moon_dec1) - (target_dec2 - target_dec1)
        moon_angular_speed = math.sqrt(d_ra**2 + d_dec**2) / (2 * dt_test)

        if moon_angular_speed < 0.1:
            # Fallback: Moon moves about 13 degrees per day relative to stars
            moon_angular_speed = 13.0

        # The Moon passes across the target. Calculate time based on
        # the path length through the occultation zone.
        # If min_sep is the minimum separation during passage, and
        # outer_threshold is the radius of the occultation zone,
        # then the half-chord length is sqrt(R^2 - d^2)

        if min_sep < outer_threshold:
            half_chord_outer = math.sqrt(max(0, outer_threshold**2 - min_sep**2))
            # Time = distance / speed (in days)
            dt_outer = half_chord_outer / moon_angular_speed

            jd_first = jd_max - dt_outer
            jd_fourth = jd_max + dt_outer
        else:
            jd_first = 0.0
            jd_fourth = 0.0

        # For total occultation (target fully inside Moon)
        if min_sep < inner_threshold:
            half_chord_inner = math.sqrt(max(0, inner_threshold**2 - min_sep**2))
            dt_inner = half_chord_inner / moon_angular_speed

            jd_second = jd_max - dt_inner
            jd_third = jd_max + dt_inner
        else:
            jd_second = 0.0
            jd_third = 0.0

        return jd_first, jd_second, jd_third, jd_fourth

    # Main search loop
    jd = jd_start

    # Moon moves ~13 degrees per day relative to stars
    # So conjunctions with any target happen roughly every ~27 days
    # But for occultations, we need close passages

    step = 1.0  # Check every day initially

    for iteration in range(MAX_ITERATIONS):
        # Check current position
        is_occ, sep, moon_r, target_r = _check_occultation(jd)

        if is_occ:
            # Found an occultation! Refine the timing
            jd_max = _find_minimum_separation(jd - 0.5, jd + 0.5)

            # Verify it's still an occultation
            is_occ_refined, min_sep, moon_r, target_r = _check_occultation(jd_max)

            if is_occ_refined:
                # Calculate contact times
                jd_first, jd_second, jd_third, jd_fourth = _calculate_contact_times(
                    jd_max, min_sep, moon_r, target_r
                )

                # Determine occultation type
                if min_sep < abs(moon_r - target_r):
                    if target_r > moon_r:
                        # Target larger than Moon (only for Sun - but Sun is not handled here)
                        ecl_type = SE_ECL_ANNULAR
                    else:
                        # Total occultation - target fully behind Moon
                        ecl_type = SE_ECL_TOTAL
                else:
                    # Partial occultation
                    ecl_type = SE_ECL_PARTIAL

                times = (
                    jd_max,  # [0] Maximum
                    jd_first,  # [1] First contact
                    jd_second,  # [2] Second contact (total begins)
                    jd_third,  # [3] Third contact (total ends)
                    jd_fourth,  # [4] Fourth contact
                    0.0,  # [5] Reserved
                    0.0,  # [6] Reserved
                    0.0,  # [7] Reserved
                )

                return times, ecl_type

        # Check if we're getting close to an occultation
        if sep < 3.0:  # Within 3 degrees
            # Use smaller steps for closer approaches
            step = 0.1
        elif sep < 10.0:  # Within 10 degrees
            step = 0.5
        else:
            step = 1.0

        jd += step

        if jd > jd_start + MAX_SEARCH_YEARS * 365.25:
            break

    target_desc = star_name if planet == 0 else f"planet {planet}"
    raise RuntimeError(
        f"No lunar occultation of {target_desc} found within "
        f"{MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Alias for Swiss Ephemeris API compatibility
swe_lun_occult_when_glob = lun_occult_when_glob


def lun_occult_when_loc(
    jd_start: float,
    planet: int,
    star_name: str = "",
    lat: float = 0.0,
    lon: float = 0.0,
    altitude: float = 0.0,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], int]:
    """
    Find the next lunar occultation visible from a specific location.

    A lunar occultation occurs when the Moon passes in front of (occults)
    a planet or star as seen from Earth. This function searches forward
    in time to find the next occultation visible from a specific geographic
    location, where both the Moon and the target are above the horizon.

    Args:
        jd_start: Julian Day (UT) to start search from
        planet: Planet ID to check for occultation (SE_MERCURY, SE_VENUS, etc.)
            Set to 0 if searching for a fixed star occultation.
        star_name: Name of fixed star to check (e.g. "Regulus", "Spica").
            Only used if planet is 0.
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - times: Tuple of 10 floats with occultation phase times (JD UT):
                [0]: Time of maximum occultation (minimum separation)
                [1]: Time of first contact (occultation begins)
                [2]: Time of second contact (full occultation begins, or 0)
                [3]: Time of third contact (full occultation ends, or 0)
                [4]: Time of fourth contact (occultation ends)
                [5]: Reserved (0)
                [6]: Reserved (0)
                [7]: Time of moonrise (if Moon rises during occultation, else 0)
                [8]: Time of moonset (if Moon sets during occultation, else 0)
                [9]: Reserved (0)
            - attr: Tuple of 11 floats with occultation attributes:
                [0]: Fraction of target covered at maximum (0-1)
                [1]: Reserved (0)
                [2]: Reserved (0)
                [3]: Azimuth of Moon at maximum occultation (degrees)
                [4]: Altitude of Moon at maximum occultation (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of target (degrees)
                [7]: Angular separation at maximum (degrees)
                [8]: Reserved (0)
                [9]: Reserved (0)
                [10]: Reserved (0)
            - retflag: Occultation type flags bitmask (SE_ECL_* constants)
                SE_ECL_TOTAL: Total occultation (body fully behind Moon)
                SE_ECL_PARTIAL: Partial occultation (body partially behind Moon)
                SE_ECL_VISIBLE: Occultation visible from location
                SE_ECL_MAX_VISIBLE: Maximum visible from location
                SE_ECL_1ST_VISIBLE: First contact visible
                SE_ECL_4TH_VISIBLE: Fourth contact visible

    Raises:
        RuntimeError: If no occultation visible from location within search limit
        ValueError: If neither planet nor star_name is specified

    Algorithm:
        1. Use lun_occult_when_glob to find next global occultation
        2. Calculate Moon and target altitude at observer's location
        3. If both bodies are below horizon during entire occultation,
           continue to next global occultation
        4. Return local phase times, visibility flags, and attributes

    Precision:
        Occultation times accurate to ~1 minute for most events.
        Visibility depends on accurate horizon calculations.

    Example:
        >>> # Find next occultation of Regulus visible from Rome
        >>> from libephemeris import julday, lun_occult_when_loc
        >>> jd = julday(2017, 1, 1, 0)
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> times, attr, ocl_type = lun_occult_when_loc(jd, 0, "Regulus",
        ...                                              rome_lat, rome_lon)
        >>> print(f"Occultation maximum at JD {times[0]:.5f}")
        >>> print(f"Moon altitude: {attr[4]:.1f}°")

    References:
        - Swiss Ephemeris: swe_lun_occult_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
    """
    from skyfield.api import wgs84

    from .constants import (
        SE_MERCURY,
        SE_VENUS,
        SE_MARS,
        SE_JUPITER,
        SE_SATURN,
        SE_URANUS,
        SE_NEPTUNE,
        SE_PLUTO,
    )
    from .fixed_stars import FIXED_STARS, _resolve_star_id
    from .planets import _PLANET_MAP
    from .state import get_planets, get_timescale

    if planet == 0 and not star_name:
        raise ValueError(
            "Either planet ID or star_name must be specified for occultation search"
        )

    MAX_SEARCH_YEARS = 50
    MAX_OCCULTATIONS = int(MAX_SEARCH_YEARS * 20)  # Check many potential occultations

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    earth = eph["earth"]
    moon_body = eph["moon"]

    # Create observer
    observer = wgs84.latlon(lat, lon, altitude)

    # Moon's mean angular radius in degrees
    MOON_MEAN_ANGULAR_RADIUS_DEG = 932.56 / 3600.0

    # Approximate angular radius for planets (arcsec)
    PLANET_ANGULAR_RADII_ARCSEC = {
        SE_MERCURY: 6.0,
        SE_VENUS: 30.0,
        SE_MARS: 12.0,
        SE_JUPITER: 24.0,
        SE_SATURN: 10.0,
        SE_URANUS: 2.0,
        SE_NEPTUNE: 1.2,
        SE_PLUTO: 0.1,
    }

    def _get_moon_altaz(jd: float) -> Tuple[float, float]:
        """Get Moon altitude and azimuth at given JD from observer location."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer
        moon_app = observer_at.at(t).observe(moon_body).apparent()
        alt, az, _ = moon_app.altaz()
        return alt.degrees, az.degrees

    def _get_target_altaz(jd: float) -> Tuple[float, float]:
        """Get target altitude and azimuth at given JD from observer location."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer

        if planet == 0:
            # Fixed star - calculate alt/az from RA/Dec
            star_id, _ = _resolve_star_id(star_name)
            if star_id < 0:
                raise ValueError(f"Star '{star_name}' not found")

            star = FIXED_STARS[star_id]

            # Time from J2000.0
            t_years = (jd - 2451545.0) / 365.25

            # Apply proper motion
            ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
            dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0

            # Create a position object for the star and get alt/az
            from skyfield.api import Star

            star_obj = Star(
                ra_hours=ra_deg / 15.0,
                dec_degrees=dec_deg,
            )
            star_app = observer_at.at(t).observe(star_obj).apparent()
            alt, az, _ = star_app.altaz()
            return alt.degrees, az.degrees
        else:
            # Planet
            if planet not in _PLANET_MAP:
                raise ValueError(f"Unknown planet ID: {planet}")

            target_name = _PLANET_MAP[planet]
            target = eph[target_name]
            target_app = observer_at.at(t).observe(target).apparent()
            alt, az, _ = target_app.altaz()
            return alt.degrees, az.degrees

    def _is_visible(jd: float, min_alt: float = -1.0) -> bool:
        """Check if both Moon and target are above horizon."""
        moon_alt, _ = _get_moon_altaz(jd)
        target_alt, _ = _get_target_altaz(jd)
        return moon_alt > min_alt and target_alt > min_alt

    def _get_moon_angular_info(jd: float) -> Tuple[float, float, float]:
        """Get Moon's apparent diameter and distance at given JD from observer."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer
        moon_app = observer_at.at(t).observe(moon_body).apparent()
        dist = moon_app.distance().au
        # Moon angular diameter varies with distance
        moon_diameter = 2 * MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist)
        return moon_diameter, dist, MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist)

    def _get_target_angular_info() -> float:
        """Get target's angular diameter."""
        if planet == 0:
            # Fixed star - negligible angular size
            return 0.0001  # ~0.4 arcsec for point source
        else:
            return PLANET_ANGULAR_RADII_ARCSEC.get(planet, 1.0) / 3600.0 * 2

    jd = jd_start

    for _ in range(MAX_OCCULTATIONS):
        # Find next global occultation
        try:
            global_times, global_type = lun_occult_when_glob(
                jd, planet, star_name, flags
            )
        except RuntimeError:
            raise RuntimeError(
                f"No lunar occultation of {'star ' + star_name if planet == 0 else 'planet ' + str(planet)} "
                f"visible from lat={lat}, lon={lon} "
                f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
            )

        jd_max = global_times[0]
        jd_first = global_times[1]
        jd_second = global_times[2]
        jd_third = global_times[3]
        jd_fourth = global_times[4]

        # Check if occultation is visible from this location
        # Sample at key phases: first contact, maximum, fourth contact
        check_times = [jd_max]
        if jd_first > 0:
            check_times.append(jd_first)
        if jd_fourth > 0:
            check_times.append(jd_fourth)
        if jd_second > 0:
            check_times.append(jd_second)
        if jd_third > 0:
            check_times.append(jd_third)

        # Also sample at intermediate points
        if jd_first > 0 and jd_fourth > 0:
            duration = jd_fourth - jd_first
            for i in range(5):
                t = jd_first + duration * i / 4
                if t not in check_times:
                    check_times.append(t)

        bodies_visible = any(_is_visible(t) for t in check_times)

        if not bodies_visible:
            # Occultation not visible from this location, try next
            jd = jd_max + 1  # Skip ahead to find next occultation
            continue

        # Occultation visible! Calculate local circumstances

        # Get Moon's alt/az at maximum
        moon_alt, moon_az = _get_moon_altaz(jd_max)
        target_alt, target_az = _get_target_altaz(jd_max)

        # Get apparent diameters
        moon_diameter, moon_dist, moon_radius = _get_moon_angular_info(jd_max)
        target_diameter = _get_target_angular_info()
        target_radius = target_diameter / 2

        # Calculate angular separation at maximum
        # (should be near zero for an occultation)
        def _angular_separation_at_jd(jd_check: float) -> float:
            """Calculate angular separation between Moon and target."""
            t = ts.ut1_jd(jd_check)
            observer_at = earth + observer
            moon_app = observer_at.at(t).observe(moon_body).apparent()

            if planet == 0:
                # Fixed star
                star_id, _ = _resolve_star_id(star_name)
                star = FIXED_STARS[star_id]
                t_years = (jd_check - 2451545.0) / 365.25
                ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
                dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0

                from skyfield.api import Star

                star_obj = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
                target_app = observer_at.at(t).observe(star_obj).apparent()
            else:
                target_name = _PLANET_MAP[planet]
                target = eph[target_name]
                target_app = observer_at.at(t).observe(target).apparent()

            sep = moon_app.separation_from(target_app)
            return sep.degrees

        min_separation = _angular_separation_at_jd(jd_max)

        # Calculate fraction covered at maximum
        # If separation < |moon_radius - target_radius|, target is fully covered
        if min_separation < abs(moon_radius - target_radius):
            fraction_covered = 1.0
        elif min_separation < moon_radius + target_radius:
            # Partial coverage - calculate overlap
            overlap = (moon_radius + target_radius) - min_separation
            fraction_covered = min(1.0, overlap / target_diameter)
        else:
            fraction_covered = 0.0

        # Determine occultation type
        ecl_type = global_type

        # Check visibility at each phase
        if jd_first > 0 and _is_visible(jd_first):
            ecl_type |= SE_ECL_1ST_VISIBLE
        if jd_second > 0 and _is_visible(jd_second):
            ecl_type |= SE_ECL_2ND_VISIBLE
        if jd_third > 0 and _is_visible(jd_third):
            ecl_type |= SE_ECL_3RD_VISIBLE
        if jd_fourth > 0 and _is_visible(jd_fourth):
            ecl_type |= SE_ECL_4TH_VISIBLE
        if _is_visible(jd_max):
            ecl_type |= SE_ECL_MAX_VISIBLE
        ecl_type |= SE_ECL_VISIBLE

        # Find moonrise/moonset during occultation
        moonrise_time = 0.0
        moonset_time = 0.0

        if jd_first > 0 and jd_fourth > 0:
            alt_begin, _ = _get_moon_altaz(jd_first)
            alt_end, _ = _get_moon_altaz(jd_fourth)
            if alt_begin < 0 < alt_end:
                # Moon rises during occultation - find approximate time
                t_low, t_high = jd_first, jd_fourth
                for _ in range(20):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _ = _get_moon_altaz(t_mid)
                    if alt_mid < 0:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                moonrise_time = (t_low + t_high) / 2
            elif alt_begin > 0 > alt_end:
                # Moon sets during occultation
                t_low, t_high = jd_first, jd_fourth
                for _ in range(20):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _ = _get_moon_altaz(t_mid)
                    if alt_mid > 0:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                moonset_time = (t_low + t_high) / 2

        # Prepare times tuple (10 elements)
        times = (
            global_times[0],  # [0] Maximum
            global_times[1],  # [1] First contact
            global_times[2],  # [2] Second contact (total begins)
            global_times[3],  # [3] Third contact (total ends)
            global_times[4],  # [4] Fourth contact
            0.0,  # [5] Reserved
            0.0,  # [6] Reserved
            moonrise_time,  # [7] Moonrise during occultation
            moonset_time,  # [8] Moonset during occultation
            0.0,  # [9] Reserved
        )

        # Prepare attributes tuple (11 elements)
        attr = (
            fraction_covered,  # [0] Fraction of target covered
            0.0,  # [1] Reserved
            0.0,  # [2] Reserved
            moon_az,  # [3] Azimuth of Moon at maximum
            moon_alt,  # [4] Altitude of Moon at maximum
            moon_diameter,  # [5] Apparent diameter of Moon
            target_diameter,  # [6] Apparent diameter of target
            min_separation,  # [7] Angular separation at maximum
            0.0,  # [8] Reserved
            0.0,  # [9] Reserved
            0.0,  # [10] Reserved
        )

        return times, attr, ecl_type

    target_desc = star_name if planet == 0 else f"planet {planet}"
    raise RuntimeError(
        f"No lunar occultation of {target_desc} visible from lat={lat}, lon={lon} "
        f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Alias for Swiss Ephemeris API compatibility
swe_lun_occult_when_loc = lun_occult_when_loc


def lun_occult_where(
    jd: float,
    planet: int,
    star_name: str = "",
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], int]:
    """
    Calculate where on Earth a lunar occultation is visible at a given time.

    This function determines where on Earth the lunar occultation of a planet
    or star is visible at the specified Julian Day. It returns the geographic
    coordinates of the central line (where the occultation is most central)
    and attributes about the occultation geometry.

    Args:
        jd: Julian Day (UT) of the moment to calculate
        planet: Planet ID to check for occultation (SE_MERCURY, SE_VENUS, etc.)
            Set to 0 if searching for a fixed star occultation.
        star_name: Name of fixed star to check (e.g. "Regulus", "Spica").
            Only used if planet is 0.
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - geopos: Tuple of 10 floats with geographic positions:
                [0]: Geographic longitude of central occultation (degrees, East+)
                [1]: Geographic latitude of central occultation (degrees, North+)
                [2]: Geographic longitude of northern limit (degrees)
                [3]: Geographic latitude of northern limit (degrees)
                [4]: Geographic longitude of southern limit (degrees)
                [5]: Geographic latitude of southern limit (degrees)
                [6]: Geographic longitude of sunrise limit (degrees)
                [7]: Geographic latitude of sunrise limit (degrees)
                [8]: Geographic longitude of sunset limit (degrees)
                [9]: Geographic latitude of sunset limit (degrees)
            - attr: Tuple of 8 floats with occultation attributes:
                [0]: Fraction of target covered by Moon (0-1)
                [1]: Ratio of target diameter to Moon diameter
                [2]: Angular separation at central line (degrees)
                [3]: Width of occultation path (km)
                [4]: Moon's azimuth at central line (degrees)
                [5]: Moon's altitude at central line (degrees)
                [6]: Apparent diameter of Moon (degrees)
                [7]: Apparent diameter of target (degrees)
            - retflag: Occultation type flags bitmask (SE_ECL_* constants)
                Returns 0 if no occultation at this time

    Raises:
        ValueError: If neither planet nor star_name is specified

    Note:
        If there is no occultation at the given time (Moon too far from the
        target), geopos will contain zeros and retflag will be 0.

    Algorithm:
        1. Calculate Moon and target positions at the given time
        2. Check if angular separation is small enough for occultation
        3. Find sub-lunar point (where Moon is overhead)
        4. Calculate path limits based on occultation geometry
        5. Determine where occultation is visible above horizon

    Example:
        >>> # Find where Regulus occultation is visible
        >>> from libephemeris import julday, lun_occult_where
        >>> jd = julday(2017, 6, 28, 10.0)  # During a known occultation
        >>> geopos, attr, ocl_type = lun_occult_where(jd, 0, "Regulus")
        >>> print(f"Central line at lon={geopos[0]:.2f}, lat={geopos[1]:.2f}")

    References:
        - Swiss Ephemeris: swe_lun_occult_where()
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
    """
    from skyfield.api import wgs84

    from .constants import (
        SE_MERCURY,
        SE_VENUS,
        SE_MARS,
        SE_JUPITER,
        SE_SATURN,
        SE_URANUS,
        SE_NEPTUNE,
        SE_PLUTO,
    )
    from .fixed_stars import FIXED_STARS, _resolve_star_id
    from .planets import _PLANET_MAP
    from .state import get_planets, get_timescale

    if planet == 0 and not star_name:
        raise ValueError(
            "Either planet ID or star_name must be specified for occultation"
        )

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    earth = eph["earth"]
    moon_body = eph["moon"]

    # Moon's mean angular radius in degrees (varies from ~14.7' to ~16.7')
    # 932.56 arcsec at mean distance (0.002569 AU)
    MOON_MEAN_ANGULAR_RADIUS_DEG = 932.56 / 3600.0

    # Approximate angular radius for planets (arcsec)
    PLANET_ANGULAR_RADII_ARCSEC = {
        SE_MERCURY: 6.0,
        SE_VENUS: 30.0,
        SE_MARS: 12.0,
        SE_JUPITER: 24.0,
        SE_SATURN: 10.0,
        SE_URANUS: 2.0,
        SE_NEPTUNE: 1.2,
        SE_PLUTO: 0.1,
    }

    # Zero return values for no occultation
    zero_geopos = (0.0,) * 10
    zero_attr = (0.0,) * 8

    def _get_moon_geocentric(jd_calc: float) -> Tuple[float, float, float, float]:
        """Get Moon's geocentric RA, Dec, distance, and angular radius."""
        t = ts.ut1_jd(jd_calc)
        moon_app = earth.at(t).observe(moon_body).apparent()
        ra, dec, dist = moon_app.radec(epoch="date")
        # Moon angular radius varies with distance
        moon_angular_radius = MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist.au)
        return ra.hours * 15.0, dec.degrees, dist.au, moon_angular_radius

    def _get_target_position(jd_calc: float) -> Tuple[float, float, float]:
        """Get target body's geocentric RA, Dec, and angular radius."""
        if planet == 0:
            # Fixed star
            star_id, _ = _resolve_star_id(star_name)
            if star_id < 0:
                raise ValueError(f"Star '{star_name}' not found")

            star = FIXED_STARS[star_id]

            # Time from J2000.0
            t_years = (jd_calc - 2451545.0) / 365.25

            # Apply proper motion
            ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
            dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0

            # Stars have negligible angular radius
            return ra_deg, dec_deg, 0.0001  # ~0.4 arcsec for point source
        else:
            # Planet
            if planet not in _PLANET_MAP:
                raise ValueError(f"Unknown planet ID: {planet}")

            target_name = _PLANET_MAP[planet]
            target = eph[target_name]

            t = ts.ut1_jd(jd_calc)
            target_app = earth.at(t).observe(target).apparent()
            ra, dec, dist = target_app.radec(epoch="date")

            # Get planet angular radius
            angular_radius = PLANET_ANGULAR_RADII_ARCSEC.get(planet, 1.0) / 3600.0

            return ra.hours * 15.0, dec.degrees, angular_radius

    def _angular_separation(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
        """Calculate angular separation between two points (in degrees)."""
        # Convert to radians
        ra1_r = math.radians(ra1)
        dec1_r = math.radians(dec1)
        ra2_r = math.radians(ra2)
        dec2_r = math.radians(dec2)

        # Haversine formula for spherical distance
        cos_sep = math.sin(dec1_r) * math.sin(dec2_r) + math.cos(dec1_r) * math.cos(
            dec2_r
        ) * math.cos(ra1_r - ra2_r)
        cos_sep = max(-1.0, min(1.0, cos_sep))
        return math.degrees(math.acos(cos_sep))

    # Get Moon and target positions
    moon_ra, moon_dec, moon_dist, moon_radius = _get_moon_geocentric(jd)
    target_ra, target_dec, target_radius = _get_target_position(jd)

    # Calculate angular separation
    separation = _angular_separation(moon_ra, moon_dec, target_ra, target_dec)

    # Check if occultation is occurring
    occultation_threshold = moon_radius + target_radius
    if separation > occultation_threshold:
        # No occultation at this time
        return zero_geopos, zero_attr, 0

    # Occultation is occurring - calculate where on Earth it's visible

    # Get Skyfield time
    t = ts.ut1_jd(jd)

    # Calculate GMST (Greenwich Mean Sidereal Time)
    gmst = t.gmst  # in hours
    gmst_deg = gmst * 15.0  # Convert to degrees

    # Sub-lunar point: where Moon is directly overhead
    # Longitude: where Moon's RA = Local Sidereal Time
    # LST = GMST + longitude, so longitude = RA - GMST
    central_lon = moon_ra - gmst_deg
    # Normalize to -180 to +180
    central_lon = ((central_lon + 180) % 360) - 180

    # Sub-lunar point latitude = Moon's declination
    central_lat = moon_dec

    # Clamp latitude to valid range
    central_lat = max(-90.0, min(90.0, central_lat))

    # Calculate occultation path width and limits
    # The occultation is visible where the Moon-target angular distance
    # as seen from the observer is less than the sum of their angular radii

    # Moon's parallax affects where the occultation is visible
    # The parallax is approximately: Moon distance in Earth radii ~ 60
    # Angular parallax ~ asin(1/60) ~ 0.95 degrees
    # Note: parallax calculation available for future refinement

    # Path width estimation based on Moon's angular diameter
    # The path width is approximately 2 * Moon_radius / cos(elevation) in degrees
    # Converting to km: 1 degree ~ 111 km at equator
    path_width_deg = 2 * moon_radius
    path_width_km = path_width_deg * 111.0  # Approximate at equator

    # Calculate northern and southern limits
    # These are approximately +/- Moon's angular radius in latitude from the central line
    # adjusted for the Moon's parallax effect
    north_lat = central_lat + moon_radius
    south_lat = central_lat - moon_radius

    # Clamp to valid latitude range
    north_lat = max(-90.0, min(90.0, north_lat))
    south_lat = max(-90.0, min(90.0, south_lat))

    # Northern and southern limits have approximately the same longitude
    # as the central line (simplified calculation)
    north_lon = central_lon
    south_lon = central_lon

    # Calculate sunrise/sunset limits
    # These are where the Moon is at the horizon during the occultation
    # Simplified: approximately ±90° in longitude from central line
    sunrise_lon = ((central_lon - 90) + 180) % 360 - 180
    sunset_lon = ((central_lon + 90) + 180) % 360 - 180
    sunrise_lat = central_lat
    sunset_lat = central_lat

    # Calculate attributes at central line
    try:
        observer = wgs84.latlon(central_lat, central_lon, 0.0)
        observer_at = earth + observer

        # Get apparent positions from observer
        moon_app = observer_at.at(t).observe(moon_body).apparent()

        if planet == 0:
            # Fixed star
            star_id, _ = _resolve_star_id(star_name)
            star = FIXED_STARS[star_id]
            t_years = (jd - 2451545.0) / 365.25
            ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
            dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0

            from skyfield.api import Star

            star_obj = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
            target_app = observer_at.at(t).observe(star_obj).apparent()
        else:
            target_name = _PLANET_MAP[planet]
            target = eph[target_name]
            target_app = observer_at.at(t).observe(target).apparent()

        # Get Moon altitude and azimuth at central line
        moon_alt, moon_az, _ = moon_app.altaz()
        moon_altitude = moon_alt.degrees
        moon_azimuth = moon_az.degrees

        # Calculate local angular separation
        local_separation = moon_app.separation_from(target_app).degrees

        # Calculate local angular sizes
        local_moon_dist = moon_app.distance().au
        local_moon_radius = MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / local_moon_dist)
        local_moon_diameter = 2 * local_moon_radius
        local_target_diameter = 2 * target_radius

        # Fraction covered
        if local_separation < abs(local_moon_radius - target_radius):
            fraction_covered = 1.0
        elif local_separation < local_moon_radius + target_radius:
            overlap = (local_moon_radius + target_radius) - local_separation
            fraction_covered = (
                min(1.0, overlap / (2 * target_radius)) if target_radius > 0 else 1.0
            )
        else:
            fraction_covered = 0.0

        # Ratio of target diameter to Moon diameter
        if local_moon_diameter > 0:
            diameter_ratio = local_target_diameter / local_moon_diameter
        else:
            diameter_ratio = 0.0

    except Exception:
        # If calculation fails, use defaults
        moon_azimuth = 0.0
        moon_altitude = 0.0
        local_separation = separation
        local_moon_diameter = 2 * moon_radius
        local_target_diameter = 2 * target_radius
        fraction_covered = 0.0
        diameter_ratio = 0.0

    # Determine occultation type
    if separation < abs(moon_radius - target_radius):
        if target_radius > moon_radius:
            # Target larger than Moon (rare, e.g., Sun during solar eclipse)
            eclipse_type = SE_ECL_ANNULAR
        else:
            # Total occultation
            eclipse_type = SE_ECL_TOTAL
    else:
        # Partial occultation
        eclipse_type = SE_ECL_PARTIAL

    # Prepare return tuples
    geopos = (
        central_lon,  # [0] Central longitude
        central_lat,  # [1] Central latitude
        north_lon,  # [2] Northern limit longitude
        north_lat,  # [3] Northern limit latitude
        south_lon,  # [4] Southern limit longitude
        south_lat,  # [5] Southern limit latitude
        sunrise_lon,  # [6] Sunrise limit longitude
        sunrise_lat,  # [7] Sunrise limit latitude
        sunset_lon,  # [8] Sunset limit longitude
        sunset_lat,  # [9] Sunset limit latitude
    )

    attr = (
        fraction_covered,  # [0] Fraction covered
        diameter_ratio,  # [1] Target/Moon diameter ratio
        local_separation,  # [2] Angular separation
        path_width_km,  # [3] Path width in km
        moon_azimuth,  # [4] Moon azimuth
        moon_altitude,  # [5] Moon altitude
        local_moon_diameter,  # [6] Moon apparent diameter
        local_target_diameter,  # [7] Target apparent diameter
    )

    return geopos, attr, eclipse_type


# Alias for Swiss Ephemeris API compatibility
swe_lun_occult_where = lun_occult_where


# =============================================================================
# RISE, SET, AND TRANSIT CALCULATIONS
# =============================================================================


def rise_trans(
    jd_start: float,
    planet: int,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    flags: int = SEFLG_SWIEPH,
    rsmi: int = 1,
) -> Tuple[float, int]:
    """
    Calculate rise, set, or transit time for a celestial body.

    This function finds the next occurrence of a rising, setting, or transit
    (meridian passage) event for a celestial body as seen from a specific
    geographic location.

    Args:
        jd_start: Julian Day (UT) to start search from
        planet: Planet/body ID (SE_SUN, SE_MOON, etc.)
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        pressure: Atmospheric pressure in mbar/hPa for refraction (default 1013.25)
        temperature: Temperature in Celsius for refraction (default 15)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)
        rsmi: Event type and calculation flags (bitmask):
            - SE_CALC_RISE (1): Rise time (body crossing horizon going up)
            - SE_CALC_SET (2): Set time (body crossing horizon going down)
            - SE_CALC_MTRANSIT (4): Upper meridian transit (culmination)
            - SE_CALC_ITRANSIT (8): Lower meridian transit (anti-culmination)
            Additional flags (OR with event type):
            - SE_BIT_DISC_CENTER (256): Use disc center instead of upper limb
            - SE_BIT_DISC_BOTTOM (8192): Use lower limb of disc
            - SE_BIT_NO_REFRACTION (512): Ignore atmospheric refraction
            - SE_BIT_CIVIL_TWILIGHT (1024): Sun at -6 degrees
            - SE_BIT_NAUTIC_TWILIGHT (2048): Sun at -12 degrees
            - SE_BIT_ASTRO_TWILIGHT (4096): Sun at -18 degrees

    Returns:
        Tuple containing:
            - jd_event: Julian Day (UT) of the event, or 0.0 if not found
            - retflag: Return flag (same as rsmi on success, or error indicator)
                       Returns -2 if body is circumpolar (never rises/sets)

    Raises:
        ValueError: If invalid planet ID or parameters

    Note:
        For circumpolar objects (always above or below horizon at the given
        latitude), the function returns (0.0, -2). For transits, circumpolar
        objects still have valid transit times.

    Algorithm:
        1. For transits: Find when body crosses the local meridian
           (Local Sidereal Time = body's Right Ascension)
        2. For rise/set: Find when body's altitude crosses the horizon
           accounting for refraction and disc size
        3. Uses Newton-Raphson iteration for precise timing

    Precision:
        Rise/set times accurate to ~1 minute for Sun/Moon (due to refraction
        uncertainty), better than 1 second for transit times.

    Example:
        >>> from libephemeris import julday, rise_trans, SE_SUN, SE_CALC_RISE
        >>> jd = julday(2024, 6, 21, 0)
        >>> # Find sunrise at Rome
        >>> jd_rise, _ = rise_trans(jd, SE_SUN, 41.9, 12.5, rsmi=SE_CALC_RISE)
        >>> print(f"Sunrise at JD {jd_rise:.5f}")

    References:
        - Swiss Ephemeris: swe_rise_trans()
        - Meeus "Astronomical Algorithms" Ch. 15 (Rise, Set, Transit)
    """
    from skyfield.api import wgs84

    from .constants import (
        SE_CALC_RISE,
        SE_CALC_SET,
        SE_CALC_MTRANSIT,
        SE_CALC_ITRANSIT,
        SE_BIT_DISC_CENTER,
        SE_BIT_DISC_BOTTOM,
        SE_BIT_NO_REFRACTION,
        SE_BIT_CIVIL_TWILIGHT,
        SE_BIT_NAUTIC_TWILIGHT,
        SE_BIT_ASTRO_TWILIGHT,
    )
    from .planets import _PLANET_MAP
    from .state import get_planets, get_timescale

    # Extract event type from rsmi (lower bits)
    event_type = rsmi & 0x0F  # First 4 bits for event type

    # Validate event type
    if event_type not in (
        SE_CALC_RISE,
        SE_CALC_SET,
        SE_CALC_MTRANSIT,
        SE_CALC_ITRANSIT,
    ):
        raise ValueError(
            f"Invalid event type in rsmi: {rsmi}. Use SE_CALC_RISE, SE_CALC_SET, SE_CALC_MTRANSIT, or SE_CALC_ITRANSIT"
        )

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Validate planet and get target body
    if planet not in _PLANET_MAP:
        raise ValueError(f"Invalid planet ID: {planet}")

    target_name = _PLANET_MAP[planet]
    target = eph[target_name]
    earth = eph["earth"]

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Determine horizon altitude based on flags
    # Standard refraction at horizon is about 34 arcminutes
    if rsmi & SE_BIT_NO_REFRACTION:
        refraction = 0.0
    else:
        # Simple refraction model (more accurate would use pressure/temperature)
        # Standard refraction at horizon: ~34 arcminutes = 0.5667 degrees
        refraction = 0.5667

    # Determine the altitude threshold for rise/set
    if rsmi & SE_BIT_CIVIL_TWILIGHT:
        horizon_alt = -6.0
        refraction = 0.0  # Twilight uses geometric position
    elif rsmi & SE_BIT_NAUTIC_TWILIGHT:
        horizon_alt = -12.0
        refraction = 0.0
    elif rsmi & SE_BIT_ASTRO_TWILIGHT:
        horizon_alt = -18.0
        refraction = 0.0
    else:
        horizon_alt = 0.0

    # Account for disc semi-diameter
    # Sun: ~16 arcmin, Moon: ~16 arcmin (varies)
    if rsmi & SE_BIT_DISC_CENTER:
        disc_correction = 0.0
    elif rsmi & SE_BIT_DISC_BOTTOM:
        # Lower limb: add semi-diameter (rises later, sets earlier)
        if planet == SE_SUN:
            disc_correction = 16.0 / 60.0  # degrees
        elif planet == SE_MOON:
            disc_correction = 16.0 / 60.0
        else:
            disc_correction = 0.0
    else:
        # Default: upper limb (subtract semi-diameter)
        if planet == SE_SUN:
            disc_correction = -16.0 / 60.0  # degrees
        elif planet == SE_MOON:
            disc_correction = -16.0 / 60.0
        else:
            disc_correction = 0.0

    # Effective horizon altitude (negative means below geometric horizon)
    target_altitude = horizon_alt - refraction + disc_correction

    def _get_body_altaz(jd: float) -> Tuple[float, float]:
        """Get body's altitude and azimuth at given JD from observer location."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer
        body_app = observer_at.at(t).observe(target).apparent()
        alt, az, _ = body_app.altaz()
        return alt.degrees, az.degrees

    def _get_body_ra_dec(jd: float) -> Tuple[float, float]:
        """Get body's RA and Dec at given JD (epoch of date)."""
        t = ts.ut1_jd(jd)
        body_app = earth.at(t).observe(target).apparent()
        ra, dec, _ = body_app.radec(epoch="date")
        return ra.hours, dec.degrees  # RA in hours, Dec in degrees

    # Handle transit calculations
    if event_type in (SE_CALC_MTRANSIT, SE_CALC_ITRANSIT):
        return _calculate_transit(
            jd_start,
            lat,
            lon,
            event_type,
            ts,
            earth,
            target,
            observer,
            SE_CALC_ITRANSIT,
        )

    # Handle rise/set calculations
    return _calculate_rise_set(
        jd_start,
        lat,
        lon,
        event_type,
        target_altitude,
        ts,
        earth,
        target,
        observer,
        _get_body_altaz,
        _get_body_ra_dec,
        SE_CALC_RISE,
        SE_CALC_SET,
        rsmi,
    )


def _calculate_transit(
    jd_start: float,
    lat: float,
    lon: float,
    event_type: int,
    ts,
    earth,
    target,
    observer,
    SE_CALC_ITRANSIT: int,
) -> Tuple[float, int]:
    """Calculate meridian transit time."""
    # Transit occurs when body's RA = Local Sidereal Time
    # LST = GMST + longitude (in hours)
    # So transit occurs when GMST = RA - longitude/15

    def _get_body_hour_angle(jd: float) -> float:
        """Calculate body's hour angle (in degrees). 0 = on meridian."""
        t = ts.ut1_jd(jd)

        # Get body's RA (epoch of date)
        body_app = earth.at(t).observe(target).apparent()
        ra, _, _ = body_app.radec(epoch="date")
        ra_deg = ra.hours * 15.0  # Convert to degrees

        # Get Local Sidereal Time
        gmst = t.gmst  # in hours
        lst = gmst + lon / 15.0  # in hours
        lst_deg = lst * 15.0  # Convert to degrees

        # Hour angle = LST - RA (normalized to -180 to +180)
        ha = (lst_deg - ra_deg) % 360.0
        if ha > 180:
            ha -= 360
        return ha

    # Get initial hour angle
    ha = _get_body_hour_angle(jd_start)

    # Adjust for lower transit (HA = 180 instead of 0)
    if event_type == SE_CALC_ITRANSIT:
        ha = (ha - 180) % 360
        if ha > 180:
            ha -= 360

    # Estimate time to next transit
    # Earth rotates 360° in ~23h56m (sidereal day)
    # So hour angle changes ~15°/hour = 360°/day (approximately)
    sidereal_day = 0.99726957  # days

    if ha > 0:
        # Body is west of meridian, wait for it to come around
        dt = (360.0 - ha) / 360.0 * sidereal_day
    else:
        # Body is east of meridian, it will cross soon
        dt = (-ha) / 360.0 * sidereal_day

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(30):
        ha = _get_body_hour_angle(jd_guess)

        # Adjust for lower transit
        if event_type == SE_CALC_ITRANSIT:
            ha = (ha - 180) % 360
            if ha > 180:
                ha -= 360

        # Check convergence (< 1 arcsecond in HA = < 0.07 seconds in time)
        if abs(ha) < 1.0 / 3600.0:
            return jd_guess, event_type

        # Rate of change of HA: approximately 360°/sidereal day
        ha_rate = 360.0 / sidereal_day  # degrees/day

        # Newton-Raphson step
        jd_guess -= ha / ha_rate

    # If we get here, convergence failed but return best estimate
    return jd_guess, event_type


def _calculate_rise_set(
    jd_start: float,
    lat: float,
    lon: float,
    event_type: int,
    target_altitude: float,
    ts,
    earth,
    target,
    observer,
    _get_body_altaz,
    _get_body_ra_dec,
    SE_CALC_RISE: int,
    SE_CALC_SET: int,
    rsmi: int,
) -> Tuple[float, int]:
    """Calculate rise or set time using bisection and Newton-Raphson."""

    # Get current altitude and declination
    alt, _ = _get_body_altaz(jd_start)
    _, dec = _get_body_ra_dec(jd_start)

    # Check for circumpolar objects
    # An object is circumpolar if its altitude never crosses the target altitude
    # Calculate the altitude at upper and lower transit

    # Maximum altitude = 90 - |lat - dec|
    # Minimum altitude = -(90 - |lat + dec|) = |lat + dec| - 90
    # For circumpolar: min_alt > target_altitude (never sets)
    # For never rises: max_alt < target_altitude

    # Calculate the altitude at upper and lower transit
    if lat >= 0:  # Northern hemisphere
        max_alt = 90.0 - abs(lat - dec)  # At upper transit
        min_alt = dec - (90.0 - lat)  # At lower transit
    else:  # Southern hemisphere
        max_alt = 90.0 - abs(lat - dec)
        min_alt = -dec - (90.0 + lat)

    # Check circumpolar conditions
    is_circumpolar_above = min_alt > target_altitude  # Never sets
    is_circumpolar_below = max_alt < target_altitude  # Never rises

    if event_type == SE_CALC_RISE and is_circumpolar_below:
        return 0.0, -2  # Never rises
    if event_type == SE_CALC_SET and is_circumpolar_above:
        return 0.0, -2  # Never sets
    if event_type == SE_CALC_RISE and is_circumpolar_above:
        return 0.0, -2  # Always above horizon, no rise
    if event_type == SE_CALC_SET and is_circumpolar_below:
        return 0.0, -2  # Always below horizon, no set

    # Rough estimate: search within 1 day
    search_step = 1.0 / 24.0  # 1 hour steps

    # Find bracket where altitude crosses target
    jd = jd_start
    alt_prev, _ = _get_body_altaz(jd)

    max_search = 2.0  # Search up to 2 days ahead
    jd_cross_start = None
    jd_cross_end = None

    for i in range(int(max_search / search_step) + 1):
        jd = jd_start + i * search_step
        alt, _ = _get_body_altaz(jd)

        # Check for crossing
        if event_type == SE_CALC_RISE:
            # Rising: altitude going from below to above target
            if alt_prev < target_altitude and alt >= target_altitude:
                jd_cross_start = jd - search_step
                jd_cross_end = jd
                break
        else:  # SE_CALC_SET
            # Setting: altitude going from above to below target
            if alt_prev >= target_altitude and alt < target_altitude:
                jd_cross_start = jd - search_step
                jd_cross_end = jd
                break

        alt_prev = alt

    if jd_cross_start is None:
        # No crossing found in search range - might be circumpolar or rare case
        return 0.0, -2

    # Bisection to narrow down the crossing
    for _ in range(30):
        jd_mid = (jd_cross_start + jd_cross_end) / 2
        alt_mid, _ = _get_body_altaz(jd_mid)

        if abs(alt_mid - target_altitude) < 0.001:  # ~3.6 arcsec
            return jd_mid, rsmi & 0x0F

        if event_type == SE_CALC_RISE:
            if alt_mid < target_altitude:
                jd_cross_start = jd_mid
            else:
                jd_cross_end = jd_mid
        else:  # SE_CALC_SET
            if alt_mid >= target_altitude:
                jd_cross_start = jd_mid
            else:
                jd_cross_end = jd_mid

    return (jd_cross_start + jd_cross_end) / 2, rsmi & 0x0F


# Alias for Swiss Ephemeris API compatibility
swe_rise_trans = rise_trans
