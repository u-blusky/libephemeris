"""
Crossing event calculations for libephemeris.

Finds exact times when the Sun or Moon cross specific ecliptic longitudes.
Uses Newton-Raphson iteration for sub-arcsecond precision, with Brent's method
as a fallback near retrograde stations.

Functions:
- swe_solcross_ut: Sun crossing events (e.g., ingresses, equinoxes)
- swe_mooncross_ut: Moon crossing events (for lunar mansion calculations)
- swe_cross_ut: Generic planet crossing

Precision: Newton-Raphson convergence tolerance
Tolerance: 0.001 arcsecond for Sun (sub-milliarcsecond), 0.05 arcsecond for Moon, 0.1 arcsecond for planets
Iterations: Adaptive based on planet speed:
    - 50 for Sun/fast planets (speed >= 0.1°/day)
    - 30 for Moon (rapid convergence due to high speed ~13°/day)
    - 60 for slow planets like Saturn (0.01 <= speed < 0.1°/day)
    - 80 for very slow planets like Pluto (0.001 <= speed < 0.01°/day)
    - 100 near retrograde stations (speed < 0.001°/day)
Typical convergence: 5-8 iterations for Sun, 7-12 for Moon (up to 20 near nodes)

Station Handling:
    When a planet is near a retrograde station (velocity ~0°/day), Newton-Raphson
    can have convergence problems due to dividing by near-zero derivatives. In these
    cases, the algorithm automatically switches to Brent's method, which is more
    robust as it only requires bracketing the root, not computing derivatives.

Algorithm: Initial linear estimate + Newton-Raphson refinement + Brent's fallback
References:
    - Meeus "Astronomical Algorithms" Ch. 5 (interpolation)
    - Brent, R.P. (1973) "Algorithms for Minimization without Derivatives"
"""

from typing import Callable, Tuple

from .constants import SEFLG_SWIEPH, SEFLG_SPEED, SEFLG_HELCTR, SE_SUN, SE_MOON
from .planets import swe_calc_ut, swe_calc

# Station detection threshold: speed below this indicates proximity to retrograde station
# At stations, Newton-Raphson can fail due to near-zero derivative (speed)
# Typical station speeds: Mercury ~0.05°/day slowing to 0, outer planets much slower
STATION_SPEED_THRESHOLD = 0.001  # degrees/day

# Newton-Raphson convergence constants
# 0.1 arcsecond tolerance for pyswisseph compatibility
NR_TOLERANCE = 0.1 / 3600.0  # 0.1 arcsecond in degrees
# Tighter tolerance for Sun: pyswisseph achieves < 0.001 arcsec (sub-milliarcsecond)
NR_TOLERANCE_SUN = 0.001 / 3600.0  # 0.001 arcsecond in degrees
# Moon tolerance: 0.05 arcsecond - sub-arcsecond precision with fast convergence
# due to Moon's high speed (~13°/day). Tighter than generic planets.
NR_TOLERANCE_MOON = 0.05 / 3600.0  # 0.05 arcsecond in degrees
NR_MAX_ITER_SUN = 50  # Max iterations for Sun
# Moon iterations: 30 is sufficient given rapid convergence.
# - Longitude crossings: Moon moves ~13°/day, typical 7-12 iterations
# - Node crossings: latitude speed ~1°/day at nodes, up to 15-20 iterations
# - 30 provides 1.5-2.5x safety margin for edge cases near nodes
# See CALCS.md for Moon motion details (±5.15° latitude, ~27 day cycle)
NR_MAX_ITER_MOON = 30  # Max iterations for Moon
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


def _is_near_station(speed: float) -> bool:
    """
    Check if a planet is near a retrograde station.

    Near stations, the planet's speed approaches zero, causing Newton-Raphson
    to have convergence issues (division by near-zero). This function detects
    when we should switch to a more robust method like Brent's.

    Args:
        speed: Planet speed in degrees/day

    Returns:
        bool: True if near station (speed below threshold)
    """
    return abs(speed) < STATION_SPEED_THRESHOLD


def _brent_find_crossing(
    get_position_func: Callable[[float], Tuple[float, float]],
    x2cross: float,
    jd_a: float,
    jd_b: float,
    tolerance: float,
    max_iter: int = 100,
) -> float:
    """
    Find exact crossing time using Brent's method.

    Brent's method combines bisection, secant, and inverse quadratic interpolation.
    It's more robust than Newton-Raphson when derivatives are unreliable (near stations).

    Args:
        get_position_func: Function(jd) -> (longitude, speed) in degrees
        x2cross: Target ecliptic longitude in degrees
        jd_a: Start of bracket (Julian Day)
        jd_b: End of bracket (Julian Day)
        tolerance: Convergence tolerance in degrees
        max_iter: Maximum iterations

    Returns:
        float: Julian Day of crossing

    Raises:
        RuntimeError: If no root is bracketed or convergence fails

    Algorithm:
        Brent's method is guaranteed to converge if the function changes sign
        in [a, b]. It uses:
        - Bisection for safety
        - Secant method for faster convergence when applicable
        - Inverse quadratic interpolation when even faster

    References:
        Brent, R. P. (1973). Algorithms for Minimization without Derivatives.
    """

    def f(jd: float) -> float:
        """Return angular difference from target (signed, -180 to 180)."""
        lon, _ = get_position_func(jd)
        diff = (lon - x2cross) % 360.0
        if diff > 180:
            diff -= 360
        return diff

    fa = f(jd_a)
    fb = f(jd_b)

    # Check if root is bracketed
    if fa * fb > 0:
        # Root not bracketed - try to expand the bracket
        # This can happen if we guessed the bracket wrong
        raise RuntimeError(
            f"Brent's method: root not bracketed. f(a)={fa:.6f}, f(b)={fb:.6f}"
        )

    # Ensure |f(a)| >= |f(b)| (b is the better guess)
    if abs(fa) < abs(fb):
        jd_a, jd_b = jd_b, jd_a
        fa, fb = fb, fa

    c = jd_a
    fc = fa
    mflag = True
    d = 0.0  # Will be set before use

    for _ in range(max_iter):
        if abs(fb) < tolerance:
            return jd_b

        if fa != fc and fb != fc:
            # Inverse quadratic interpolation
            s = (
                jd_a * fb * fc / ((fa - fb) * (fa - fc))
                + jd_b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
            )
        else:
            # Secant method
            s = jd_b - fb * (jd_b - jd_a) / (fb - fa)

        # Conditions for using bisection instead
        cond1 = not (
            (3 * jd_a + jd_b) / 4 < s < jd_b or jd_b < s < (3 * jd_a + jd_b) / 4
        )
        cond2 = mflag and abs(s - jd_b) >= abs(jd_b - c) / 2
        cond3 = not mflag and abs(s - jd_b) >= abs(c - d) / 2
        cond4 = mflag and abs(jd_b - c) < tolerance
        cond5 = not mflag and abs(c - d) < tolerance

        if cond1 or cond2 or cond3 or cond4 or cond5:
            # Bisection
            s = (jd_a + jd_b) / 2
            mflag = True
        else:
            mflag = False

        fs = f(s)
        d = c
        c = jd_b
        fc = fb

        if fa * fs < 0:
            jd_b = s
            fb = fs
        else:
            jd_a = s
            fa = fs

        # Ensure |f(a)| >= |f(b)|
        if abs(fa) < abs(fb):
            jd_a, jd_b = jd_b, jd_a
            fa, fb = fb, fa

    # Return best estimate if max iterations reached
    return jd_b


def _find_bracket_for_crossing(
    get_position_func: Callable[[float], Tuple[float, float]],
    x2cross: float,
    jd_start: float,
    jd_end: float,
    num_samples: int = 20,
) -> Tuple[float, float]:
    """
    Find a bracket [jd_a, jd_b] where the crossing occurs.

    Samples the position at regular intervals to find where the difference
    from target changes sign (indicating a crossing).

    Args:
        get_position_func: Function(jd) -> (longitude, speed) in degrees
        x2cross: Target ecliptic longitude in degrees
        jd_start: Start of search interval
        jd_end: End of search interval
        num_samples: Number of samples to take

    Returns:
        Tuple[float, float]: (jd_a, jd_b) bracket containing the crossing

    Raises:
        RuntimeError: If no crossing is found in the interval
    """
    step = (jd_end - jd_start) / num_samples

    def get_diff(jd: float) -> float:
        lon, _ = get_position_func(jd)
        diff = (lon - x2cross) % 360.0
        if diff > 180:
            diff -= 360
        return diff

    prev_jd = jd_start
    prev_diff = get_diff(jd_start)

    for i in range(1, num_samples + 1):
        curr_jd = jd_start + i * step
        curr_diff = get_diff(curr_jd)

        # Check for sign change (crossing)
        if prev_diff * curr_diff <= 0:
            return (prev_jd, curr_jd)

        prev_jd = curr_jd
        prev_diff = curr_diff

    raise RuntimeError(
        f"No crossing found in interval [{jd_start}, {jd_end}] for target {x2cross}°"
    )


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
    # Geocentric average speeds (°/day) - slower planets need more iterations
    # Note: geocentric speeds are affected by retrograde motion, values are
    # approximate averages during direct motion
    typical_speeds = {
        2: 1.4,  # Mercury
        3: 1.2,  # Venus
        4: 0.5,  # Mars
        5: 0.08,  # Jupiter
        6: 0.03,  # Saturn
        7: 0.01,  # Uranus
        8: 0.006,  # Neptune
        9: 0.004,  # Pluto (very slow, ~0.004°/day)
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

    # Helper function for Brent's method fallback
    def get_position(jd_time: float) -> Tuple[float, float]:
        pos_result, _ = swe_calc_ut(jd_time, planet_id, flag | SEFLG_SPEED)
        return pos_result[0], pos_result[3]

    # Check if we're near a retrograde station - use Brent's method for robustness
    if _is_near_station(speed):
        # Near station: Newton-Raphson may fail due to division by near-zero speed
        # Use Brent's method which only requires bracketing, not derivatives
        try:
            # Estimate search window based on typical speeds
            # At stations, planet barely moves, so we need a wider bracket
            search_window = (
                max(30.0, abs(diff / speed_default) * 1.5)
                if speed_default > 0
                else 60.0
            )
            jd_bracket_start = jd_ut
            jd_bracket_end = jd_ut + search_window

            # Find bracket containing the crossing
            jd_a, jd_b = _find_bracket_for_crossing(
                get_position, x2cross, jd_bracket_start, jd_bracket_end
            )

            # Use Brent's method to find the exact crossing
            return _brent_find_crossing(
                get_position, x2cross, jd_a, jd_b, NR_TOLERANCE, max_iter
            )
        except RuntimeError:
            # If Brent's method fails, fall through to Newton-Raphson as last resort
            pass

    # Newton-Raphson iteration
    jd = jd_guess
    station_fallback_triggered = False

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

        # Detect if we've encountered a station during iteration
        if _is_near_station(speed) and not station_fallback_triggered:
            station_fallback_triggered = True
            # Switch to Brent's method for robustness
            try:
                # Create a bracket around current position
                bracket_size = max(
                    5.0, abs(diff / speed_default) if speed_default > 0 else 10.0
                )
                jd_a, jd_b = _find_bracket_for_crossing(
                    get_position, x2cross, jd - bracket_size / 2, jd + bracket_size
                )
                return _brent_find_crossing(
                    get_position,
                    x2cross,
                    jd_a,
                    jd_b,
                    NR_TOLERANCE,
                    max_iter - iteration,
                )
            except RuntimeError:
                # If Brent fails, continue with Newton-Raphson
                pass

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

    # Helper function for Brent's method fallback (heliocentric)
    def get_helio_position(jd_time: float) -> Tuple[float, float]:
        pos_result, _ = swe_calc_ut(jd_time, planet_id, helio_flag)
        return pos_result[0], pos_result[3]

    # Check if we're dealing with a very slow planet - use Brent's method for robustness
    # Heliocentric planets don't have true stations but very slow planets can still
    # benefit from the more robust method
    if _is_near_station(speed):
        try:
            search_window = (
                max(60.0, abs(diff / speed_default) * 1.5)
                if speed_default > 0
                else 120.0
            )
            jd_bracket_start = jd_ut
            jd_bracket_end = jd_ut + search_window

            jd_a, jd_b = _find_bracket_for_crossing(
                get_helio_position, x2cross, jd_bracket_start, jd_bracket_end
            )

            return _brent_find_crossing(
                get_helio_position, x2cross, jd_a, jd_b, NR_TOLERANCE, max_iter
            )
        except RuntimeError:
            pass  # Fall through to Newton-Raphson

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

    # Helper function for Brent's method fallback (heliocentric TT)
    def get_helio_position_tt(jd_time: float) -> Tuple[float, float]:
        pos_result, _ = swe_calc(jd_time, planet_id, helio_flag)
        return pos_result[0], pos_result[3]

    # Check if we're dealing with a very slow planet - use Brent's method for robustness
    if _is_near_station(speed):
        try:
            search_window = (
                max(60.0, abs(diff / speed_default) * 1.5)
                if speed_default > 0
                else 120.0
            )
            jd_bracket_start = jd_tt
            jd_bracket_end = jd_tt + search_window

            jd_a, jd_b = _find_bracket_for_crossing(
                get_helio_position_tt, x2cross, jd_bracket_start, jd_bracket_end
            )

            return _brent_find_crossing(
                get_helio_position_tt, x2cross, jd_a, jd_b, NR_TOLERANCE, max_iter
            )
        except RuntimeError:
            pass  # Fall through to Newton-Raphson

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


# =============================================================================
# RETROGRADE STATION HANDLING
# =============================================================================
#
# Retrograde stations occur when a planet appears to stop and reverse direction
# as seen from Earth. At station points, the planet's velocity approaches zero,
# which can cause numerical instabilities in calculations.
#
# Two types of stations:
#   - Retrograde Station (SR): Planet slows, stops, and begins retrograde motion
#   - Direct Station (SD): Planet stops retrograde motion and resumes direct motion
#
# The functions below provide robust handling for calculations near these points.
# =============================================================================

# Tolerance for station finding (velocity threshold)
# A planet is considered "stationary" when |velocity| < this value
STATION_VELOCITY_TOLERANCE = 1e-6  # degrees/day


def is_retrograde(planet_id: int, jd_ut: float, flag: int = SEFLG_SWIEPH) -> bool:
    """
    Check if a planet is currently in retrograde motion.

    A planet is retrograde when its geocentric longitude is decreasing
    (negative velocity in longitude).

    Args:
        planet_id: Planet ID (SE_MERCURY, SE_VENUS, SE_MARS, etc.)
        jd_ut: Julian Day (UT) to check
        flag: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        bool: True if planet is retrograde (velocity < 0), False otherwise

    Note:
        Sun and Moon never go retrograde from geocentric perspective.
        Only planets beyond Earth (Mars-Pluto) and Mercury/Venus can appear retrograde.

    Example:
        >>> is_retrograde(SE_MERCURY, jd_now)
        True  # Mercury is currently retrograde
    """
    if planet_id in (SE_SUN, SE_MOON):
        # Sun and Moon are never retrograde from geocentric view
        return False

    try:
        pos, _ = swe_calc_ut(jd_ut, planet_id, flag | SEFLG_SPEED)
        return pos[3] < 0
    except Exception:
        return False


def get_station_type(planet_id: int, jd_ut: float, flag: int = SEFLG_SWIEPH) -> str:
    """
    Determine if a planet is near a station and what type.

    Args:
        planet_id: Planet ID (SE_MERCURY, SE_VENUS, SE_MARS, etc.)
        jd_ut: Julian Day (UT) to check
        flag: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        str: One of:
            - "direct": Planet is in direct (prograde) motion
            - "retrograde": Planet is in retrograde motion
            - "stationary_retrograde": Near station, about to go retrograde
            - "stationary_direct": Near station, about to resume direct motion

    Example:
        >>> get_station_type(SE_MARS, jd_station)
        'stationary_retrograde'
    """
    if planet_id in (SE_SUN, SE_MOON):
        return "direct"

    try:
        pos, _ = swe_calc_ut(jd_ut, planet_id, flag | SEFLG_SPEED)
        speed = pos[3]

        if _is_near_station(speed):
            # Near station - determine type by checking velocity trend
            # Look at velocity slightly before and after
            dt = 1.0  # 1 day
            pos_before, _ = swe_calc_ut(jd_ut - dt, planet_id, flag | SEFLG_SPEED)
            pos_after, _ = swe_calc_ut(jd_ut + dt, planet_id, flag | SEFLG_SPEED)

            speed_before = pos_before[3]
            speed_after = pos_after[3]

            # If speed is decreasing (becoming more negative), heading to retrograde
            if speed_before > speed_after:
                return "stationary_retrograde"
            else:
                return "stationary_direct"
        elif speed < 0:
            return "retrograde"
        else:
            return "direct"

    except Exception:
        return "direct"


def _brent_find_station(
    get_speed_func,
    jd_a: float,
    jd_b: float,
    tolerance: float = STATION_VELOCITY_TOLERANCE,
    max_iter: int = 100,
) -> float:
    """
    Find exact station time using Brent's method for root finding.

    Finds when velocity = 0 (station point) within the bracket [jd_a, jd_b].

    Args:
        get_speed_func: Function(jd) -> velocity in degrees/day
        jd_a: Start of bracket (Julian Day)
        jd_b: End of bracket (Julian Day)
        tolerance: Convergence tolerance for velocity
        max_iter: Maximum iterations

    Returns:
        float: Julian Day when velocity = 0

    Raises:
        RuntimeError: If root not bracketed or convergence fails
    """
    fa = get_speed_func(jd_a)
    fb = get_speed_func(jd_b)

    # Check if root is bracketed (velocity changes sign)
    if fa * fb > 0:
        raise RuntimeError(
            f"Station not bracketed: speed(a)={fa:.6f}, speed(b)={fb:.6f}"
        )

    # Ensure |f(a)| >= |f(b)|
    if abs(fa) < abs(fb):
        jd_a, jd_b = jd_b, jd_a
        fa, fb = fb, fa

    c = jd_a
    fc = fa
    mflag = True
    d = 0.0

    for _ in range(max_iter):
        if abs(fb) < tolerance:
            return jd_b

        if fa != fc and fb != fc:
            # Inverse quadratic interpolation
            s = (
                jd_a * fb * fc / ((fa - fb) * (fa - fc))
                + jd_b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
            )
        else:
            # Secant method
            s = jd_b - fb * (jd_b - jd_a) / (fb - fa)

        # Conditions for bisection fallback
        cond1 = not (
            (3 * jd_a + jd_b) / 4 < s < jd_b or jd_b < s < (3 * jd_a + jd_b) / 4
        )
        cond2 = mflag and abs(s - jd_b) >= abs(jd_b - c) / 2
        cond3 = not mflag and abs(s - jd_b) >= abs(c - d) / 2
        cond4 = mflag and abs(jd_b - c) < tolerance * 86400
        cond5 = not mflag and abs(c - d) < tolerance * 86400

        if cond1 or cond2 or cond3 or cond4 or cond5:
            s = (jd_a + jd_b) / 2
            mflag = True
        else:
            mflag = False

        fs = get_speed_func(s)
        d = c
        c = jd_b
        fc = fb

        if fa * fs < 0:
            jd_b = s
            fb = fs
        else:
            jd_a = s
            fa = fs

        if abs(fa) < abs(fb):
            jd_a, jd_b = jd_b, jd_a
            fa, fb = fb, fa

    return jd_b


def _find_station_bracket(
    planet_id: int,
    jd_start: float,
    jd_end: float,
    flag: int,
    num_samples: int = 50,
) -> Tuple[float, float]:
    """
    Find a bracket [jd_a, jd_b] containing a station (velocity sign change).

    Args:
        planet_id: Planet ID
        jd_start: Start of search interval
        jd_end: End of search interval
        flag: Calculation flags
        num_samples: Number of samples to take

    Returns:
        Tuple[float, float]: (jd_a, jd_b) bracket containing a station

    Raises:
        RuntimeError: If no station found in interval
    """
    step = (jd_end - jd_start) / num_samples

    pos_prev, _ = swe_calc_ut(jd_start, planet_id, flag | SEFLG_SPEED)
    speed_prev = pos_prev[3]
    jd_prev = jd_start

    for i in range(1, num_samples + 1):
        jd_curr = jd_start + i * step
        pos_curr, _ = swe_calc_ut(jd_curr, planet_id, flag | SEFLG_SPEED)
        speed_curr = pos_curr[3]

        # Check for sign change (station)
        if speed_prev * speed_curr <= 0:
            return (jd_prev, jd_curr)

        jd_prev = jd_curr
        speed_prev = speed_curr

    raise RuntimeError(
        f"No station found for planet {planet_id} in [{jd_start}, {jd_end}]"
    )


def swe_find_station_ut(
    planet_id: int,
    jd_ut: float,
    station_type: str = "any",
    flag: int = SEFLG_SWIEPH,
) -> Tuple[float, str]:
    """
    Find the next planetary station (stationary point) after a given date.

    A station occurs when a planet's apparent velocity reaches zero, marking
    the transition between direct and retrograde motion (or vice versa).

    Args:
        planet_id: Planet ID (SE_MERCURY through SE_PLUTO)
        jd_ut: Julian Day (UT) to start search from
        station_type: Type of station to find:
            - "any": Find next station regardless of type (default)
            - "retrograde" or "SR": Find next station where planet goes retrograde
            - "direct" or "SD": Find next station where planet resumes direct motion
        flag: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple[float, str]: (Julian Day of station, station type "SR" or "SD")

    Raises:
        ValueError: If planet_id is Sun or Moon (never station)
        RuntimeError: If convergence fails or no station found

    Note:
        Sun and Moon never go retrograde from Earth's perspective.

    Station Periods (approximate time between stations):
        - Mercury: ~116 days (synodic period)
        - Venus: ~584 days
        - Mars: ~780 days
        - Jupiter: ~399 days
        - Saturn: ~378 days
        - Uranus: ~370 days
        - Neptune: ~367 days
        - Pluto: ~367 days

    Precision:
        Typically finds station to within ~1 minute of arc in velocity

    Example:
        >>> jd_station, stype = swe_find_station_ut(SE_MERCURY, jd_now)
        >>> print(f"Mercury stations {stype} on JD {jd_station}")
    """
    if planet_id in (SE_SUN, SE_MOON):
        raise ValueError("Sun and Moon do not have retrograde stations")

    # Normalize station_type
    if station_type.upper() in ("SR", "RETROGRADE"):
        station_type = "retrograde"
    elif station_type.upper() in ("SD", "DIRECT"):
        station_type = "direct"
    else:
        station_type = "any"

    # Typical synodic periods (days) - approximate time between same-type stations
    synodic_periods = {
        2: 116,  # Mercury
        3: 584,  # Venus
        4: 780,  # Mars
        5: 399,  # Jupiter
        6: 378,  # Saturn
        7: 370,  # Uranus
        8: 367,  # Neptune
        9: 367,  # Pluto
    }

    # Search window is half synodic period (time between SR and SD)
    search_window = synodic_periods.get(planet_id, 400) / 2

    def get_speed(jd: float) -> float:
        pos, _ = swe_calc_ut(jd, planet_id, flag | SEFLG_SPEED)
        return pos[3]

    # Get current motion direction
    current_speed = get_speed(jd_ut)

    jd_search_start = jd_ut
    max_attempts = 4  # Allow searching up to 2 full synodic periods

    for attempt in range(max_attempts):
        try:
            # Find the next station bracket
            jd_search_end = jd_search_start + search_window
            jd_a, jd_b = _find_station_bracket(
                planet_id, jd_search_start, jd_search_end, flag
            )

            # Find exact station time
            jd_station = _brent_find_station(
                get_speed, jd_a, jd_b, STATION_VELOCITY_TOLERANCE
            )

            # Determine station type by looking at motion before and after
            speed_before = get_speed(jd_station - 1.0)
            speed_after = get_speed(jd_station + 1.0)

            if speed_before > 0 and speed_after < 0:
                found_type = "SR"  # Stationary Retrograde
            else:
                found_type = "SD"  # Stationary Direct

            # Check if this matches requested type
            if station_type == "any":
                return (jd_station, found_type)
            elif station_type == "retrograde" and found_type == "SR":
                return (jd_station, found_type)
            elif station_type == "direct" and found_type == "SD":
                return (jd_station, found_type)

            # Not the right type, search for next station
            jd_search_start = jd_station + 1.0

        except RuntimeError:
            # No station found in this window, extend search
            jd_search_start = jd_search_start + search_window

    raise RuntimeError(f"Could not find {station_type} station for planet {planet_id}")


def swe_next_retrograde_ut(
    planet_id: int, jd_ut: float, flag: int = SEFLG_SWIEPH
) -> Tuple[float, float]:
    """
    Find the next retrograde period for a planet.

    Returns the start (SR) and end (SD) Julian Days for the next retrograde period.

    Args:
        planet_id: Planet ID (SE_MERCURY through SE_PLUTO)
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags

    Returns:
        Tuple[float, float]: (JD of retrograde start, JD of retrograde end)

    Raises:
        ValueError: If planet is Sun or Moon
        RuntimeError: If stations cannot be found

    Example:
        >>> jd_sr, jd_sd = swe_next_retrograde_ut(SE_MERCURY, jd_now)
        >>> print(f"Mercury Rx: {jd_sr} to {jd_sd}")
    """
    # Check if currently retrograde
    if is_retrograde(planet_id, jd_ut, flag):
        # Find the end of current retrograde (SD), then find next SR
        jd_sd, _ = swe_find_station_ut(planet_id, jd_ut, "direct", flag)
        jd_sr, _ = swe_find_station_ut(planet_id, jd_sd + 1.0, "retrograde", flag)
        jd_sd_next, _ = swe_find_station_ut(planet_id, jd_sr + 1.0, "direct", flag)
        return (jd_sr, jd_sd_next)
    else:
        # Find next SR and then the following SD
        jd_sr, _ = swe_find_station_ut(planet_id, jd_ut, "retrograde", flag)
        jd_sd, _ = swe_find_station_ut(planet_id, jd_sr + 1.0, "direct", flag)
        return (jd_sr, jd_sd)


def calc_velocity_at_station(
    planet_id: int,
    jd_station: float,
    flag: int = SEFLG_SWIEPH,
    dt: float = 1.0,
) -> Tuple[float, float, float]:
    """
    Calculate stable velocity components near a station point.

    Near stations, the standard 1-second numerical differentiation can have
    numerical noise issues. This function uses a wider timestep for more
    stable velocity estimation.

    Args:
        planet_id: Planet ID
        jd_station: Julian Day near the station
        flag: Calculation flags
        dt: Timestep in days for differentiation (default 1.0 day)

    Returns:
        Tuple[float, float, float]: (velocity_lon, velocity_lat, velocity_dist)
            in degrees/day, degrees/day, AU/day respectively

    Note:
        The velocity_lon will be very close to 0 at the actual station point.
        This function is useful for getting stable readings when the standard
        swe_calc_ut velocity might show numerical noise.

    Example:
        >>> vlon, vlat, vdist = calc_velocity_at_station(SE_MARS, jd_station)
        >>> print(f"Station velocity: {vlon:.6f} deg/day")  # Should be ~0
    """
    # Get positions at t-dt and t+dt
    pos_before, _ = swe_calc_ut(jd_station - dt, planet_id, flag)
    pos_after, _ = swe_calc_ut(jd_station + dt, planet_id, flag)

    # Calculate velocity using central difference
    lon_before, lat_before, dist_before = pos_before[0], pos_before[1], pos_before[2]
    lon_after, lat_after, dist_after = pos_after[0], pos_after[1], pos_after[2]

    # Handle longitude wraparound
    dlon = lon_after - lon_before
    if dlon > 180:
        dlon -= 360
    elif dlon < -180:
        dlon += 360

    v_lon = dlon / (2 * dt)
    v_lat = (lat_after - lat_before) / (2 * dt)
    v_dist = (dist_after - dist_before) / (2 * dt)

    return (v_lon, v_lat, v_dist)


def get_station_info(planet_id: int, jd_ut: float, flag: int = SEFLG_SWIEPH) -> dict:
    """
    Get comprehensive information about the nearest station.

    Finds the nearest station to the given date and returns detailed
    information about it.

    Args:
        planet_id: Planet ID (SE_MERCURY through SE_PLUTO)
        jd_ut: Julian Day (UT)
        flag: Calculation flags

    Returns:
        dict: Station information containing:
            - "jd_station": Julian Day of nearest station
            - "station_type": "SR" (stationary retrograde) or "SD" (stationary direct)
            - "days_to_station": Days until station (negative if past)
            - "longitude_at_station": Ecliptic longitude at station
            - "is_currently_retrograde": Current retrograde status
            - "velocity": Current velocity in deg/day

    Example:
        >>> info = get_station_info(SE_MERCURY, jd_now)
        >>> print(f"Next station: {info['station_type']} in {info['days_to_station']:.1f} days")
    """
    if planet_id in (SE_SUN, SE_MOON):
        raise ValueError("Sun and Moon do not have stations")

    # Get current state
    pos, _ = swe_calc_ut(jd_ut, planet_id, flag | SEFLG_SPEED)
    current_velocity = pos[3]
    current_retrograde = current_velocity < 0

    # Find next station
    try:
        jd_next, stype_next = swe_find_station_ut(planet_id, jd_ut, "any", flag)

        # Get position at station
        pos_station, _ = swe_calc_ut(jd_next, planet_id, flag)

        return {
            "jd_station": jd_next,
            "station_type": stype_next,
            "days_to_station": jd_next - jd_ut,
            "longitude_at_station": pos_station[0],
            "is_currently_retrograde": current_retrograde,
            "velocity": current_velocity,
        }
    except RuntimeError as e:
        # Return partial info if station search fails
        return {
            "jd_station": None,
            "station_type": None,
            "days_to_station": None,
            "longitude_at_station": None,
            "is_currently_retrograde": current_retrograde,
            "velocity": current_velocity,
            "error": str(e),
        }
