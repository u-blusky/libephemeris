"""
Crossing event calculations for libephemeris.

Finds exact times when the Sun or Moon cross specific ecliptic longitudes.
Uses Newton-Raphson iteration for sub-arcsecond precision.

Functions:
- swe_solcross_ut: Sun crossing events (e.g., ingresses, equinoxes)
- swe_mooncross_ut: Moon crossing events (for lunar mansion calculations)
- swe_cross_ut: Generic planet crossing

Precision: Newton-Raphson convergence tolerance
Tolerance: 0.001 arcsecond for Sun (sub-milliarcsecond), 0.05 arcsecond for Moon, 0.1 arcsecond for planets
Iterations: Adaptive based on planet speed:
    - 50 for Sun/Moon/fast planets (speed >= 0.1°/day)
    - 60 for slow planets like Saturn (0.01 <= speed < 0.1°/day)
    - 80 for very slow planets like Pluto (0.001 <= speed < 0.01°/day)
    - 100 near retrograde stations (speed < 0.001°/day)
Typical convergence: 5-8 iterations for Sun, 7-12 for Moon

Algorithm: Initial linear estimate + Newton-Raphson refinement
References: Meeus "Astronomical Algorithms" Ch. 5 (interpolation)
"""

from .constants import SEFLG_SWIEPH, SEFLG_SPEED, SEFLG_HELCTR, SE_SUN, SE_MOON
from .planets import swe_calc_ut, swe_calc

# Newton-Raphson convergence constants
# 0.1 arcsecond tolerance for pyswisseph compatibility
NR_TOLERANCE = 0.1 / 3600.0  # 0.1 arcsecond in degrees
# Tighter tolerance for Sun: pyswisseph achieves < 0.001 arcsec (sub-milliarcsecond)
NR_TOLERANCE_SUN = 0.001 / 3600.0  # 0.001 arcsecond in degrees
# Moon tolerance: 0.05 arcsecond - sub-arcsecond precision with fast convergence
# due to Moon's high speed (~13°/day). Tighter than generic planets.
NR_TOLERANCE_MOON = 0.05 / 3600.0  # 0.05 arcsecond in degrees
NR_MAX_ITER_SUN = 50  # Max iterations for Sun
NR_MAX_ITER_MOON = 50  # Max iterations for Moon
NR_MAX_ITER_PLANET = 50  # Max iterations for generic planets
NR_MAX_ITER_HELIO = 60  # Max iterations for heliocentric (slow planets)
NR_MAX_ITER_VERY_SLOW = 80  # Max iterations for very slow planets (Pluto, TNOs)
NR_MAX_ITER_STATION = 100  # Max iterations near retrograde stations


def _get_adaptive_max_iterations(speed: float) -> int:
    """
    Calculate adaptive iteration limit based on planet speed.

    For slow planets like Pluto (~0.004°/day) or near retrograde stations
    (speed approaching 0), more iterations are needed for Newton-Raphson
    convergence.

    Args:
        speed: Planet speed in degrees/day

    Returns:
        int: Maximum iterations to use

    Speed thresholds:
        |speed| >= 0.1: 50 iterations (normal planets)
        0.01 <= |speed| < 0.1: 60 iterations (slow planets like outer planets)
        0.001 <= |speed| < 0.01: 80 iterations (very slow: Pluto, TNOs)
        |speed| < 0.001: 100 iterations (near station/retrograde)
    """
    abs_speed = abs(speed)

    if abs_speed >= 0.1:
        return NR_MAX_ITER_PLANET
    elif abs_speed >= 0.01:
        return NR_MAX_ITER_HELIO
    elif abs_speed >= 0.001:
        return NR_MAX_ITER_VERY_SLOW
    else:
        return NR_MAX_ITER_STATION


def swe_solcross_ut(x2cross: float, jd_ut: float, flag: int = SEFLG_SWIEPH) -> float:
    """
    Find when the Sun crosses a specific ecliptic longitude.

    Searches FORWARD in time for the next crossing after jd_ut.

    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        float: Julian Day of crossing (UT)

    Raises:
        RuntimeError: If convergence fails or calculation error occurs

    Algorithm:
        1. Get current Sun position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.001 arcsecond (sub-milliarcsecond precision)

    Precision:
        Typically < 0.001 arcsecond (< 0.03 seconds of time for Sun)

    Example:
        >>> # Find next Aries ingress (0°)
        >>> jd_ingress = swe_solcross_ut(0.0, jd_now)
        >>> # Find summer solstice (90°)
        >>> jd_solstice = swe_solcross_ut(90.0, jd_now)
    """
    x2cross = x2cross % 360.0

    try:
        pos, _ = swe_calc_ut(jd_ut, SE_SUN, flag | SEFLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except Exception as e:
        raise RuntimeError(f"Failed to calculate Sun position: {e}")

    # Calculate angular distance to target (always forward)
    diff = (x2cross - lon_start) % 360.0

    # Handle retrograde motion
    if speed < 0 and diff > 0:
        diff -= 360.0

    # If already very close, look for next complete crossing
    if abs(diff) < 1e-5:
        if speed > 0:
            diff += 360.0
        else:
            diff -= 360.0

    # Initial time estimate (linear approximation)
    if speed == 0:
        speed = 0.9856  # Average Sun motion ~1°/day

    dt_guess = diff / speed
    jd_guess = jd_ut + dt_guess

    # Newton-Raphson iteration
    jd = jd_guess
    for iteration in range(NR_MAX_ITER_SUN):
        try:
            pos, _ = swe_calc_ut(jd, SE_SUN, flag | SEFLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(
                f"Failed to calculate Sun position during iteration: {e}"
            )

        # Angular difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.001 arcsecond = 0.001/3600 degree for Sun)
        if abs(diff) < NR_TOLERANCE_SUN:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.01:
            speed = 0.9856

        jd += diff / speed

        # Safety: prevent divergence
        if abs(jd - jd_guess) > 366:
            raise RuntimeError("Solar crossing search diverged")

    raise RuntimeError("Maximum iterations reached in solar crossing calculation")


def swe_solcross(x2cross: float, jd_tt: float, flag: int = SEFLG_SWIEPH) -> float:
    """
    Find when the Sun crosses a specific ecliptic longitude (TT version).

    This is the Terrestrial Time version of swe_solcross_ut(). Takes Julian Day
    in TT (Terrestrial Time, also known as Ephemeris Time) instead of UT.

    Searches FORWARD in time for the next crossing after jd_tt.

    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        jd_tt: Julian Day in Terrestrial Time (TT/ET) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        float: Julian Day of crossing (TT)

    Raises:
        RuntimeError: If convergence fails or calculation error occurs

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use swe_solcross_ut() instead.

    Algorithm:
        1. Get current Sun position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.001 arcsecond (sub-milliarcsecond precision)

    Precision:
        Typically < 0.001 arcsecond (< 0.03 seconds of time for Sun)

    Example:
        >>> # Find next Aries ingress (0°) using TT
        >>> jd_ingress_tt = swe_solcross(0.0, jd_tt_now)
        >>> # Find summer solstice (90°)
        >>> jd_solstice_tt = swe_solcross(90.0, jd_tt_now)
    """
    x2cross = x2cross % 360.0

    try:
        pos, _ = swe_calc(jd_tt, SE_SUN, flag | SEFLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except Exception as e:
        raise RuntimeError(f"Failed to calculate Sun position: {e}")

    # Calculate angular distance to target (always forward)
    diff = (x2cross - lon_start) % 360.0

    # Handle retrograde motion
    if speed < 0 and diff > 0:
        diff -= 360.0

    # If already very close, look for next complete crossing
    if abs(diff) < 1e-5:
        if speed > 0:
            diff += 360.0
        else:
            diff -= 360.0

    # Initial time estimate (linear approximation)
    if speed == 0:
        speed = 0.9856  # Average Sun motion ~1°/day

    dt_guess = diff / speed
    jd_guess = jd_tt + dt_guess

    jd = jd_guess
    for iteration in range(NR_MAX_ITER_SUN):
        try:
            pos, _ = swe_calc(jd, SE_SUN, flag | SEFLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(
                f"Failed to calculate Sun position during iteration: {e}"
            )

        # Angular difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.001 arcsecond for Sun)
        if abs(diff) < NR_TOLERANCE_SUN:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.01:
            speed = 0.9856

        jd += diff / speed

        # Safety: prevent divergence
        if abs(jd - jd_guess) > 366:
            raise RuntimeError("Solar crossing search diverged")

    raise RuntimeError("Maximum iterations reached in solar crossing calculation")


def swe_mooncross_ut(x2cross: float, jd_ut: float, flag: int = SEFLG_SWIEPH) -> float:
    """
    Find when the Moon crosses a specific ecliptic longitude.

    Searches FORWARD in time for the next crossing after jd_ut.

    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        float: Julian Day of crossing (UT)

    Raises:
        RuntimeError: If convergence fails or calculation error occurs

    Note:
        Moon moves ~13° per day (27.3 day cycle).
        More variable speed than Sun due to orbit eccentricity.

    Precision:
        Typically < 0.05 arcsecond (sub-arcsecond) for fast convergence

    Example:
        >>> # Find next new moon (Sun-Moon conjunction at same longitude)
        >>> sun_pos, _ = swe_calc_ut(jd_now, SE_SUN, SEFLG_SWIEPH)
        >>> jd_new_moon = swe_mooncross_ut(sun_pos[0], jd_now)
    """
    x2cross = x2cross % 360.0

    try:
        pos, _ = swe_calc_ut(jd_ut, SE_MOON, flag | SEFLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except Exception as e:
        raise RuntimeError(f"Failed to calculate Moon position: {e}")

    # Calculate initial guess for NEXT crossing
    diff = (x2cross - lon_start) % 360.0

    if speed < 0 and diff > 0:
        diff -= 360.0

    if abs(diff) < 1e-5:
        if speed > 0:
            diff += 360.0
        else:
            diff -= 360.0

    # Initial time estimate
    if speed == 0:
        speed = 13.176  # Average Moon motion ~13.18°/day

    dt_guess = diff / speed
    jd_guess = jd_ut + dt_guess

    # Newton-Raphson iteration
    jd = jd_guess
    for iteration in range(NR_MAX_ITER_MOON):
        try:
            pos, _ = swe_calc_ut(jd, SE_MOON, flag | SEFLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(
                f"Failed to calculate Moon position during iteration: {e}"
            )

        # Difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.05 arcsecond for Moon)
        if abs(diff) < NR_TOLERANCE_MOON:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.1:
            speed = 13.176

        jd += diff / speed

        # Safety check
        if abs(jd - jd_guess) > 31:  # More than a month
            raise RuntimeError("Moon crossing search diverged")

    raise RuntimeError("Maximum iterations reached in moon crossing calculation")


def swe_mooncross(x2cross: float, jd_tt: float, flag: int = SEFLG_SWIEPH) -> float:
    """
    Find when the Moon crosses a specific ecliptic longitude (TT version).

    This is the Terrestrial Time version of swe_mooncross_ut(). Takes Julian Day
    in TT (Terrestrial Time, also known as Ephemeris Time) instead of UT.

    Searches FORWARD in time for the next crossing after jd_tt.

    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        jd_tt: Julian Day in Terrestrial Time (TT/ET) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        float: Julian Day of crossing (TT)

    Raises:
        RuntimeError: If convergence fails or calculation error occurs

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use swe_mooncross_ut() instead.

        Moon moves ~13° per day (27.3 day cycle).
        More variable speed than Sun due to orbit eccentricity.

    Algorithm:
        1. Get current Moon position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.05 arcsecond (sub-arcsecond precision)

    Precision:
        Typically < 0.05 arcsecond (sub-arcsecond) for fast convergence

    Example:
        >>> # Find next new moon (Sun-Moon conjunction at same longitude) using TT
        >>> sun_pos, _ = swe_calc(jd_tt_now, SE_SUN, SEFLG_SWIEPH)
        >>> jd_new_moon_tt = swe_mooncross(sun_pos[0], jd_tt_now)
    """
    x2cross = x2cross % 360.0

    try:
        pos, _ = swe_calc(jd_tt, SE_MOON, flag | SEFLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except Exception as e:
        raise RuntimeError(f"Failed to calculate Moon position: {e}")

    # Calculate initial guess for NEXT crossing
    diff = (x2cross - lon_start) % 360.0

    if speed < 0 and diff > 0:
        diff -= 360.0

    if abs(diff) < 1e-5:
        if speed > 0:
            diff += 360.0
        else:
            diff -= 360.0

    # Initial time estimate
    if speed == 0:
        speed = 13.176  # Average Moon motion ~13.18°/day

    dt_guess = diff / speed
    jd_guess = jd_tt + dt_guess

    jd = jd_guess
    for iteration in range(NR_MAX_ITER_MOON):
        try:
            pos, _ = swe_calc(jd, SE_MOON, flag | SEFLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(
                f"Failed to calculate Moon position during iteration: {e}"
            )

        # Difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.05 arcsecond for Moon)
        if abs(diff) < NR_TOLERANCE_MOON:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.1:
            speed = 13.176

        jd += diff / speed

        # Safety check
        if abs(jd - jd_guess) > 31:  # More than a month
            raise RuntimeError("Moon crossing search diverged")

    raise RuntimeError("Maximum iterations reached in moon crossing calculation")


def swe_mooncross_node_ut(jd_ut: float, flag: int = SEFLG_SWIEPH) -> float:
    """
    Find when the Moon crosses its own orbital node (ascending or descending).

    The Moon crosses a node when its ecliptic latitude becomes zero - i.e., when
    it crosses the ecliptic plane. This is important for eclipse calculations,
    as eclipses can only occur when the Sun and Moon are near the lunar nodes.

    Searches FORWARD in time for the next node crossing after jd_ut.

    Args:
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        float: Julian Day of node crossing (UT)

    Raises:
        RuntimeError: If convergence fails or calculation error occurs

    Algorithm:
        1. Get current Moon latitude and latitude velocity
        2. Linear estimate: dt = -latitude / latitude_velocity
        3. Refine with Newton-Raphson until latitude is ~0
        4. Converge to < 0.05 arcsecond (sub-arcsecond precision)

    Note:
        The function finds the NEXT node crossing regardless of whether it's
        ascending (latitude going from negative to positive) or descending
        (latitude going from positive to negative).

        Moon crosses each node approximately every 13.6 days (half the nodal
        month of ~27.2 days).

    Example:
        >>> # Find next lunar node crossing
        >>> jd_node = swe_mooncross_node_ut(jd_now)
        >>> # Check if ascending or descending by examining latitude velocity
        >>> pos, _ = swe_calc_ut(jd_node, SE_MOON, SEFLG_SPEED)
        >>> is_ascending = pos[4] > 0  # positive lat velocity = ascending
    """
    # Half nodal month - time between successive node crossings
    HALF_NODAL_MONTH = 13.6

    try:
        pos, _ = swe_calc_ut(jd_ut, SE_MOON, flag | SEFLG_SPEED)
        lat = pos[1]  # ecliptic latitude
        lat_speed = pos[4]  # latitude velocity in degrees/day
    except Exception as e:
        raise RuntimeError(f"Failed to calculate Moon position: {e}")

    # If latitude velocity is zero or very small, use average value
    # Moon's latitude varies between about ±5.15° with period ~27.2 days
    # Average latitude speed at zero crossing: ~1.0°/day
    if abs(lat_speed) < 0.1:
        lat_speed = 1.0 if lat >= 0 else -1.0

    # Initial time estimate to reach latitude = 0
    dt_guess = -lat / lat_speed

    # If dt_guess is negative or very small, the crossing is behind us or
    # we're right at one. We need to find the NEXT crossing in the future.
    # The strategy: check at multiple points to bracket the next crossing
    if dt_guess < 0.1:  # Less than ~2.4 hours into future
        # Search in steps to find where latitude changes sign
        # The next crossing is within 0 to ~13.6 days
        jd_search_start = jd_ut + 0.5  # Start half day ahead to avoid current crossing

        # Get sign of latitude at search start
        pos_start, _ = swe_calc_ut(jd_search_start, SE_MOON, flag | SEFLG_SPEED)
        lat_sign = 1 if pos_start[1] >= 0 else -1

        # Step forward in 2-day increments to find sign change
        for step in range(8):  # Up to 16 days
            jd_check = jd_search_start + step * 2.0
            pos_check, _ = swe_calc_ut(jd_check, SE_MOON, flag | SEFLG_SPEED)
            current_sign = 1 if pos_check[1] >= 0 else -1

            if current_sign != lat_sign:
                # Found a sign change - crossing is between jd_check-2 and jd_check
                # Use the midpoint as starting guess and refine
                jd_guess = jd_check - 1.0
                break
        else:
            # Fallback: use half nodal month
            jd_guess = jd_ut + HALF_NODAL_MONTH
    else:
        jd_guess = jd_ut + dt_guess

    jd = jd_guess
    for iteration in range(NR_MAX_ITER_MOON):
        try:
            pos, _ = swe_calc_ut(jd, SE_MOON, flag | SEFLG_SPEED)
            lat = pos[1]
            lat_speed = pos[4]
        except Exception as e:
            raise RuntimeError(
                f"Failed to calculate Moon position during iteration: {e}"
            )

        # Check convergence (< 0.05 arcsecond for Moon)
        if abs(lat) < NR_TOLERANCE_MOON:
            # Make sure we found a crossing that's actually in the future
            if jd > jd_ut + 0.001:  # At least ~1.4 minutes in future
                return jd
            # If we found a crossing too close to start, look for next one
            jd = jd + HALF_NODAL_MONTH / 2
            continue

        # Newton-Raphson step
        if abs(lat_speed) < 0.1:
            lat_speed = 1.0 if lat >= 0 else -1.0

        jd -= lat / lat_speed

        # Safety check: should find a crossing within reasonable time
        if jd < jd_ut:
            # We went backward, push forward
            jd = jd_ut + HALF_NODAL_MONTH / 2
        elif abs(jd - jd_ut) > 30:
            raise RuntimeError("Moon node crossing search diverged")

    raise RuntimeError("Maximum iterations reached in moon node crossing calculation")


def swe_mooncross_node(jd_tt: float, flag: int = SEFLG_SWIEPH) -> float:
    """
    Find when the Moon crosses its own orbital node (TT version).

    This is the Terrestrial Time version of swe_mooncross_node_ut(). Takes Julian Day
    in TT (Terrestrial Time, also known as Ephemeris Time) instead of UT.

    The Moon crosses a node when its ecliptic latitude becomes zero - i.e., when
    it crosses the ecliptic plane. This is important for eclipse calculations,
    as eclipses can only occur when the Sun and Moon are near the lunar nodes.

    Searches FORWARD in time for the next node crossing after jd_tt.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT/ET) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        float: Julian Day of node crossing (TT)

    Raises:
        RuntimeError: If convergence fails or calculation error occurs

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use swe_mooncross_node_ut() instead.

    Example:
        >>> # Find next lunar node crossing using TT
        >>> jd_node_tt = swe_mooncross_node(jd_tt_now)
    """
    # Half nodal month - time between successive node crossings
    HALF_NODAL_MONTH = 13.6

    try:
        pos, _ = swe_calc(jd_tt, SE_MOON, flag | SEFLG_SPEED)
        lat = pos[1]  # ecliptic latitude
        lat_speed = pos[4]  # latitude velocity in degrees/day
    except Exception as e:
        raise RuntimeError(f"Failed to calculate Moon position: {e}")

    # If latitude velocity is zero or very small, use average value
    if abs(lat_speed) < 0.1:
        lat_speed = 1.0 if lat >= 0 else -1.0

    # Initial time estimate to reach latitude = 0
    dt_guess = -lat / lat_speed

    # If dt_guess is negative or very small, the crossing is behind us
    if dt_guess < 0.1:  # Less than ~2.4 hours into future
        # Search in steps to find where latitude changes sign
        jd_search_start = jd_tt + 0.5

        pos_start, _ = swe_calc(jd_search_start, SE_MOON, flag | SEFLG_SPEED)
        lat_sign = 1 if pos_start[1] >= 0 else -1

        for step in range(8):  # Up to 16 days
            jd_check = jd_search_start + step * 2.0
            pos_check, _ = swe_calc(jd_check, SE_MOON, flag | SEFLG_SPEED)
            current_sign = 1 if pos_check[1] >= 0 else -1

            if current_sign != lat_sign:
                jd_guess = jd_check - 1.0
                break
        else:
            jd_guess = jd_tt + HALF_NODAL_MONTH
    else:
        jd_guess = jd_tt + dt_guess

    jd = jd_guess
    for iteration in range(NR_MAX_ITER_MOON):
        try:
            pos, _ = swe_calc(jd, SE_MOON, flag | SEFLG_SPEED)
            lat = pos[1]
            lat_speed = pos[4]
        except Exception as e:
            raise RuntimeError(
                f"Failed to calculate Moon position during iteration: {e}"
            )

        # Check convergence (< 0.05 arcsecond for Moon)
        if abs(lat) < NR_TOLERANCE_MOON:
            if jd > jd_tt + 0.001:
                return jd
            jd = jd + HALF_NODAL_MONTH / 2
            continue

        # Newton-Raphson step
        if abs(lat_speed) < 0.1:
            lat_speed = 1.0 if lat >= 0 else -1.0

        jd -= lat / lat_speed

        # Safety check
        if jd < jd_tt:
            jd = jd_tt + HALF_NODAL_MONTH / 2
        elif abs(jd - jd_tt) > 30:
            raise RuntimeError("Moon node crossing search diverged")

    raise RuntimeError("Maximum iterations reached in moon node crossing calculation")


def swe_cross_ut(
    planet_id: int, x2cross: float, jd_ut: float, flag: int = SEFLG_SWIEPH
) -> float:
    """
    Find when any planet crosses a specific ecliptic longitude.

    Generic version for all planets (Mercury, Venus, Mars, etc.).

    Args:
        planet_id: Planet ID (SE_MERCURY, SE_VENUS, etc.)
        x2cross: Target ecliptic longitude in degrees (0-360)
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags

    Returns:
        float: Julian Day of crossing (UT)

    Raises:
        RuntimeError: If convergence fails or calculation error occurs

    Note:
        Uses adaptive iteration count based on typical planet speed.
        Slower planets (Jupiter, Saturn) may need more iterations.

    Precision:
        Typically < 0.1 arcsecond

    Example:
        >>> # Mars ingress into Aries
        >>> jd_mars_aries = swe_cross_ut(SE_MARS, 0.0, jd_now)
    """
    x2cross = x2cross % 360.0

    try:
        pos, _ = swe_calc_ut(jd_ut, planet_id, flag | SEFLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]
    except Exception as e:
        raise RuntimeError(f"Failed to calculate planet position: {e}")

    # Estimate typical speed if near zero
    typical_speeds = {
        2: 1.4,  # Mercury
        3: 1.2,  # Venus
        4: 0.5,  # Mars
        5: 0.08,  # Jupiter
        6: 0.03,  # Saturn
        7: 0.01,  # Uranus
        8: 0.006,  # Neptune
    }
    speed_default = typical_speeds.get(planet_id, 0.5)

    # Calculate initial guess
    # For forward-only search, always calculate forward angular distance
    diff = (x2cross - lon_start) % 360.0

    # When retrograde, the planet will eventually turn direct and reach the target
    # So we estimate using typical (prograde) speed, not current retrograde speed
    if speed < 0:
        # Use average speed for time estimate (planet will turn direct)
        effective_speed = speed_default
    else:
        effective_speed = speed if abs(speed) > 0.001 else speed_default

    if abs(diff) < 1e-5:
        diff = 360.0  # Already at target, look for next crossing

    dt_guess = diff / effective_speed
    jd_guess = jd_ut + dt_guess

    # Adaptive iteration count based on planet speed
    max_iter = _get_adaptive_max_iterations(speed)

    # Newton-Raphson iteration
    jd = jd_guess
    for iteration in range(max_iter):
        try:
            pos, _ = swe_calc_ut(jd, planet_id, flag | SEFLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(
                f"Failed to calculate planet position during iteration: {e}"
            )

        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.1 arcsecond)
        if abs(diff) < NR_TOLERANCE:
            return jd

        # Update max_iter based on current speed (may change near stations)
        max_iter = max(max_iter, _get_adaptive_max_iterations(speed))

        if abs(speed) < 0.001:
            speed = speed_default

        jd += diff / speed

        # Safety: longer range for slower planets, also account for retrograde
        # Inner planets (Mercury, Venus) can cross same degree multiple times per year
        if abs(speed_default) < 0.1:
            max_range = 500  # Slow outer planets
        elif planet_id in (2, 3):  # Mercury, Venus
            max_range = 500  # Fast inner planets with multiple crossings/year
        else:
            max_range = 400
        if abs(jd - jd_ut) > max_range:  # Use jd_ut not jd_guess
            raise RuntimeError("Planet crossing search diverged")

    raise RuntimeError("Maximum iterations reached in planet crossing calculation")


def swe_helio_cross_ut(
    planet_id: int, x2cross: float, jd_ut: float, flag: int = SEFLG_SWIEPH
) -> float:
    """
    Find when a planet crosses a specific heliocentric ecliptic longitude.

    Calculates the time when a planet, as seen from the Sun (heliocentric
    coordinates), crosses a specific ecliptic longitude. Useful for heliocentric
    astrology calculations.

    Searches FORWARD in time for the next crossing after jd_ut.

    Args:
        planet_id: Planet ID (SE_MERCURY, SE_VENUS, SE_EARTH, SE_MARS, etc.)
                   Note: SE_EARTH can be used to find when Earth crosses a longitude
                   as seen from the Sun.
        x2cross: Target heliocentric ecliptic longitude in degrees (0-360)
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.). SEFLG_HELCTR is automatically added.

    Returns:
        float: Julian Day of crossing (UT)

    Raises:
        RuntimeError: If convergence fails or calculation error occurs

    Note:
        Heliocentric positions show where planets are relative to the Sun,
        not Earth. This is useful for:
        - Heliocentric astrology systems
        - Planetary synodic cycles
        - Finding planetary positions in solar-centered coordinates

    Algorithm:
        1. Get current heliocentric position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.1 arcsecond

    Example:
        >>> # Find when Mars crosses 0° heliocentric longitude
        >>> jd_cross = swe_helio_cross_ut(SE_MARS, 0.0, jd_now)
        >>> # Find when Earth crosses 90° as seen from the Sun
        >>> jd_earth_cross = swe_helio_cross_ut(SE_EARTH, 90.0, jd_now)
    """
    x2cross = x2cross % 360.0

    # Always add SEFLG_HELCTR for heliocentric calculations
    helio_flag = flag | SEFLG_HELCTR | SEFLG_SPEED

    try:
        pos, _ = swe_calc_ut(jd_ut, planet_id, helio_flag)
        lon_start = pos[0]
        speed = pos[3]
    except Exception as e:
        raise RuntimeError(f"Failed to calculate heliocentric planet position: {e}")

    # Estimate typical heliocentric speed if near zero
    # Heliocentric speeds are different from geocentric due to no retrograde
    typical_speeds = {
        2: 4.09,  # Mercury (heliocentric, no retro)
        3: 1.60,  # Venus
        14: 0.986,  # Earth
        4: 0.524,  # Mars
        5: 0.083,  # Jupiter
        6: 0.034,  # Saturn
        7: 0.012,  # Uranus
        8: 0.006,  # Neptune
        9: 0.004,  # Pluto
    }
    speed_default = typical_speeds.get(planet_id, 0.5)

    # Calculate initial guess
    diff = (x2cross - lon_start) % 360.0

    # Heliocentric planets don't go retrograde (except for very minor perturbations)
    # so we always search forward
    if diff < 1e-5:
        diff += 360.0

    if abs(speed) < 0.0001:
        speed = speed_default

    dt_guess = diff / speed
    jd_guess = jd_ut + dt_guess

    # Adaptive iteration count based on planet speed
    max_iter = _get_adaptive_max_iterations(speed)

    jd = jd_guess
    for iteration in range(max_iter):
        try:
            pos, _ = swe_calc_ut(jd, planet_id, helio_flag)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(
                f"Failed to calculate heliocentric position during iteration: {e}"
            )

        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.1 arcsecond = 0.1/3600 degree)
        if abs(diff) < NR_TOLERANCE:
            return jd

        # Update max_iter based on current speed (may change during iteration)
        max_iter = max(max_iter, _get_adaptive_max_iterations(speed))

        if abs(speed) < 0.0001:
            speed = speed_default

        jd += diff / speed

        # Safety: longer range for slower planets
        max_range = 500 if abs(speed_default) < 0.05 else 400
        if abs(jd - jd_guess) > max_range:
            raise RuntimeError("Heliocentric crossing search diverged")

    raise RuntimeError(
        "Maximum iterations reached in heliocentric crossing calculation"
    )


def swe_helio_cross(
    planet_id: int, x2cross: float, jd_tt: float, flag: int = SEFLG_SWIEPH
) -> float:
    """
    Find when a planet crosses a specific heliocentric ecliptic longitude (TT version).

    This is the Terrestrial Time version of swe_helio_cross_ut(). Takes Julian Day
    in TT (Terrestrial Time, also known as Ephemeris Time) instead of UT.

    Calculates the time when a planet, as seen from the Sun (heliocentric
    coordinates), crosses a specific ecliptic longitude.

    Searches FORWARD in time for the next crossing after jd_tt.

    Args:
        planet_id: Planet ID (SE_MERCURY, SE_VENUS, SE_EARTH, SE_MARS, etc.)
        x2cross: Target heliocentric ecliptic longitude in degrees (0-360)
        jd_tt: Julian Day in Terrestrial Time (TT/ET) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.). SEFLG_HELCTR is automatically added.

    Returns:
        float: Julian Day of crossing (TT)

    Raises:
        RuntimeError: If convergence fails or calculation error occurs

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use swe_helio_cross_ut() instead.

    Algorithm:
        1. Get current heliocentric position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.1 arcsecond

    Example:
        >>> # Find when Mars crosses 0° heliocentric longitude using TT
        >>> jd_cross_tt = swe_helio_cross(SE_MARS, 0.0, jd_tt_now)
    """
    x2cross = x2cross % 360.0

    # Always add SEFLG_HELCTR for heliocentric calculations
    helio_flag = flag | SEFLG_HELCTR | SEFLG_SPEED

    try:
        pos, _ = swe_calc(jd_tt, planet_id, helio_flag)
        lon_start = pos[0]
        speed = pos[3]
    except Exception as e:
        raise RuntimeError(f"Failed to calculate heliocentric planet position: {e}")

    # Estimate typical heliocentric speed if near zero
    typical_speeds = {
        2: 4.09,  # Mercury
        3: 1.60,  # Venus
        14: 0.986,  # Earth
        4: 0.524,  # Mars
        5: 0.083,  # Jupiter
        6: 0.034,  # Saturn
        7: 0.012,  # Uranus
        8: 0.006,  # Neptune
        9: 0.004,  # Pluto
    }
    speed_default = typical_speeds.get(planet_id, 0.5)

    # Calculate initial guess
    diff = (x2cross - lon_start) % 360.0

    if diff < 1e-5:
        diff += 360.0

    if abs(speed) < 0.0001:
        speed = speed_default

    dt_guess = diff / speed
    jd_guess = jd_tt + dt_guess

    # Adaptive iteration count based on planet speed
    max_iter = _get_adaptive_max_iterations(speed)

    jd = jd_guess
    for iteration in range(max_iter):
        try:
            pos, _ = swe_calc(jd, planet_id, helio_flag)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(
                f"Failed to calculate heliocentric position during iteration: {e}"
            )

        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.1 arcsecond = 0.1/3600 degree)
        if abs(diff) < NR_TOLERANCE:
            return jd

        # Update max_iter based on current speed (may change during iteration)
        max_iter = max(max_iter, _get_adaptive_max_iterations(speed))

        if abs(speed) < 0.0001:
            speed = speed_default

        jd += diff / speed

        # Safety: longer range for slower planets
        max_range = 500 if abs(speed_default) < 0.05 else 400
        if abs(jd - jd_guess) > max_range:
            raise RuntimeError("Heliocentric crossing search diverged")

    raise RuntimeError(
        "Maximum iterations reached in heliocentric crossing calculation"
    )
