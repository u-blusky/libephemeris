"""
Time conversion utilities for libephemeris.

Implements standard astronomical time functions for conversions between:
- Calendar dates and Julian Day numbers
- Gregorian and Julian calendar systems
- UT1 (Universal Time) and TT (Terrestrial Time)

Functions match the Swiss Ephemeris API for compatibility.
All algorithms follow Meeus "Astronomical Algorithms" (1998).
"""

from .constants import SE_GREG_CAL, SE_JUL_CAL, SEFLG_JPLEPH, SEFLG_SWIEPH, SEFLG_MOSEPH
from .state import get_timescale

# Julian Day of Gregorian calendar reform: Oct 15, 1582
JD_GREGORIAN_REFORM = 2299161


def swe_julday(
    year: int, month: int, day: int, hour: float, gregflag: int = SE_GREG_CAL
) -> float:
    """
    Convert calendar date to Julian Day number.

    Args:
        year: Calendar year (negative for BCE)
        month: Month (1-12)
        day: Day of month (1-31)
        hour: Decimal hour (0.0-23.999...)
        gregflag: SE_GREG_CAL (1) for Gregorian, SE_JUL_CAL (0) for Julian

    Returns:
        float: Julian Day number (days since JD 0.0 = noon Jan 1, 4713 BCE)

    Note:
        Transition date: Oct 15, 1582 (Gregorian) = Oct 5, 1582 (Julian)
        JD 2451545.0 = Jan 1, 2000 12:00 TT (J2000.0 epoch)
    """
    if month <= 2:
        year -= 1
        month += 12

    a = int(year / 100)

    if gregflag == SE_GREG_CAL:
        b = 2 - a + int(a / 4)
    else:
        b = 0

    jd = (
        int(365.25 * (year + 4716))
        + int(30.6001 * (month + 1))
        + day
        + hour / 24.0
        + b
        - 1524.5
    )
    return jd


def swe_revjul(jd: float, gregflag: int = SE_GREG_CAL) -> tuple[int, int, int, float]:
    """
    Convert Julian Day number to calendar date.

    Args:
        jd: Julian Day number
        gregflag: SE_GREG_CAL (1) for Gregorian, SE_JUL_CAL (0) for Julian

    Returns:
        tuple: (year, month, day, hour) where:
            - year: Calendar year
            - month: Month (1-12)
            - day: Integer day of month
            - hour: Decimal hour (0.0-23.999...)

    Note:
        Automatic Gregorian calendar used for JD >= 2299161 (Oct 15, 1582)
        unless Julian calendar explicitly requested.
    """
    jd = jd + 0.5
    z = int(jd)
    f = jd - z

    if z < 2299161:
        a = z
    else:
        if gregflag == SE_GREG_CAL:
            alpha = int((z - 1867216.25) / 36524.25)
            a = z + 1 + alpha - int(alpha / 4)
        else:
            a = z

    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    day = b - d - int(30.6001 * e) + f
    if e < 14:
        month = e - 1
    else:
        month = e - 13

    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    d_int = int(day)
    d_frac = day - d_int
    hour = d_frac * 24.0
    day = d_int

    return year, month, day, hour


def swe_deltat(tjd: float) -> float:
    """
    Calculate Delta T (TT - UT1) for a given Julian Day.

    Args:
        tjd: Julian Day number in UT1

    Returns:
        float: Delta T in days (TT - UT1)

    Note:
        Delta T accounts for Earth's irregular rotation and is required
        for high-precision planetary calculations. Values are obtained
        from IERS (International Earth Rotation Service) data.

        For modern dates: ~0.0008 days (~69 seconds as of 2024)
        For historical dates: Calculated from polynomial models
    """
    ts = get_timescale()
    t = ts.ut1_jd(tjd)
    # Skyfield returns delta_t in seconds, convert to days for Swiss Ephemeris API compatibility
    delta_t_seconds = float(t.delta_t)
    return delta_t_seconds / 86400.0


def swe_deltat_ex(tjd: float, ephe_flag: int = SEFLG_SWIEPH) -> tuple[float, str]:
    """
    Calculate Delta T (TT - UT1) with explicit ephemeris source specification.

    Extended version of swe_deltat() that allows specifying the ephemeris
    source for Delta T calculation.

    Args:
        tjd: Julian Day number in UT1
        ephe_flag: Ephemeris selection flag:
            - SEFLG_SWIEPH (2): Use Swiss Ephemeris/Skyfield (default)
            - SEFLG_JPLEPH (1): Use JPL ephemeris
            - SEFLG_MOSEPH (4): Use Moshier ephemeris (not supported)

    Returns:
        tuple: (delta_t, serr) where:
            - delta_t: Delta T in days (TT - UT1)
            - serr: Error message string (empty if no error/warning)

    Note:
        Since libephemeris uses Skyfield which internally uses JPL data,
        SEFLG_SWIEPH and SEFLG_JPLEPH produce identical results.

        SEFLG_MOSEPH is not supported and will return the default Delta T
        with a warning message.

    Example:
        >>> from libephemeris import swe_deltat_ex, SEFLG_SWIEPH, SEFLG_JPLEPH
        >>> dt, err = swe_deltat_ex(2451545.0, SEFLG_SWIEPH)
        >>> print(f"Delta T: {dt * 86400:.2f} seconds")
        Delta T: 63.83 seconds
    """
    ts = get_timescale()
    t = ts.ut1_jd(tjd)
    # Skyfield returns delta_t in seconds, convert to days for Swiss Ephemeris API compatibility
    delta_t_seconds = float(t.delta_t)
    delta_t = delta_t_seconds / 86400.0

    # Determine if there's a warning based on the flag
    serr = ""

    # Check for valid ephemeris flags
    ephe_selection = ephe_flag & (SEFLG_JPLEPH | SEFLG_SWIEPH | SEFLG_MOSEPH)

    if ephe_selection == SEFLG_MOSEPH:
        # Moshier ephemeris is not supported, return default with warning
        serr = "Warning: SEFLG_MOSEPH not supported, using default Delta T"
    elif ephe_selection == 0:
        # No ephemeris flag specified, use default (no warning needed)
        pass
    # SEFLG_SWIEPH and SEFLG_JPLEPH both use Skyfield/JPL internally
    # so no warning is needed for those

    return delta_t, serr


def date_conversion(
    year: int, month: int, day: int, hour: float, calendar: str
) -> tuple[int, int, int, float]:
    """
    Convert a date between Julian and Gregorian calendars.

    The function automatically detects the input calendar based on the date:
    - Dates before Oct 15, 1582 are assumed to be Julian calendar
    - Dates from Oct 15, 1582 onwards are assumed to be Gregorian calendar

    Args:
        year: Calendar year
        month: Month (1-12)
        day: Day of month (1-31)
        hour: Decimal hour (0.0-23.999...)
        calendar: Target calendar - 'j' for Julian or 'g' for Gregorian

    Returns:
        tuple: (year, month, day, hour) in the requested calendar

    Raises:
        ValueError: If calendar is not 'j' or 'g'

    Note:
        The Gregorian calendar reform occurred on Oct 15, 1582. On this date,
        the Julian calendar was 10 days behind. The function uses Julian Day
        numbers as an intermediate representation to convert between calendars.

    Example:
        >>> # Convert first Gregorian date to Julian
        >>> date_conversion(1582, 10, 15, 12.0, 'j')
        (1582, 10, 5, 12.0)
        >>> # Convert Julian date to Gregorian
        >>> date_conversion(1582, 10, 5, 12.0, 'g')
        (1582, 10, 15, 12.0)
    """
    calendar = calendar.lower()
    if calendar not in ("j", "g"):
        raise ValueError(f"calendar must be 'j' or 'g', got: {calendar!r}")

    # Determine the input calendar based on the date
    # First convert to JD using Gregorian to check the date
    jd_as_greg = swe_julday(year, month, day, hour, SE_GREG_CAL)

    # If the date is before Oct 15, 1582, assume Julian input
    if jd_as_greg < JD_GREGORIAN_REFORM:
        input_cal = SE_JUL_CAL
        jd = swe_julday(year, month, day, hour, SE_JUL_CAL)
    else:
        input_cal = SE_GREG_CAL
        jd = jd_as_greg

    # Determine target calendar flag
    target_cal = SE_JUL_CAL if calendar == "j" else SE_GREG_CAL

    # If input and target are the same, just return the original values
    if input_cal == target_cal:
        return year, month, day, hour

    # Convert via Julian Day to target calendar
    return swe_revjul(jd, target_cal)


def day_of_week(jd: float) -> int:
    """
    Calculate the day of the week for a given Julian Day number.

    Uses the formula: floor(jd + 0.5) % 7 to get 0=Monday convention.
    This matches pyswisseph's day_of_week function.

    Args:
        jd: Julian Day number

    Returns:
        int: Day of the week where:
            - 0 = Monday
            - 1 = Tuesday
            - 2 = Wednesday
            - 3 = Thursday
            - 4 = Friday
            - 5 = Saturday
            - 6 = Sunday

    Example:
        >>> from libephemeris import day_of_week, swe_julday
        >>> # J2000.0 epoch: Jan 1, 2000 was a Saturday
        >>> jd = swe_julday(2000, 1, 1, 12.0)
        >>> day_of_week(jd)
        5
    """
    import math

    return int(math.floor(jd + 0.5)) % 7
