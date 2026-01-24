"""
Time conversion utilities for libephemeris.

Implements standard astronomical time functions for conversions between:
- Calendar dates and Julian Day numbers
- Gregorian and Julian calendar systems
- UT1 (Universal Time) and TT (Terrestrial Time)

Functions match the Swiss Ephemeris API for compatibility.
All algorithms follow Meeus "Astronomical Algorithms" (1998).
"""

from typing import Any

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


def utc_to_jd(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    calendar: int = SE_GREG_CAL,
) -> tuple[float, float]:
    """
    Convert UTC date/time to Julian Day numbers, properly handling leap seconds.

    Unlike julday() which assumes UT1 input, this function takes UTC input and
    correctly accounts for the difference between UTC and UT1 (DUT1) and between
    UTC and TT (including leap seconds).

    Args:
        year: Calendar year (negative for BCE)
        month: Month (1-12)
        day: Day of month (1-31)
        hour: Hour (0-23)
        minute: Minute (0-59)
        second: Second (0-60, allowing for leap seconds)
        calendar: SE_GREG_CAL (1) for Gregorian, SE_JUL_CAL (0) for Julian

    Returns:
        tuple: (jd_et, jd_ut) where:
            - jd_et: Julian Day in TT (Terrestrial Time / Ephemeris Time)
            - jd_ut: Julian Day in UT1 (Universal Time)

    Note:
        - UTC includes leap seconds while UT1 follows Earth's rotation
        - |UTC - UT1| is always < 0.9 seconds by definition
        - TT = TAI + 32.184 seconds, where TAI is atomic time
        - For dates before 1972 (when UTC was standardized), the function
          treats the input as UT1 approximation and still provides proper
          TT/UT1 conversion using historical Delta T values

    Example:
        >>> from libephemeris import utc_to_jd, SE_GREG_CAL
        >>> # J2000.0 epoch: Jan 1, 2000 12:00:00 UTC
        >>> jd_tt, jd_ut = utc_to_jd(2000, 1, 1, 12, 0, 0.0, SE_GREG_CAL)
        >>> print(f"JD(TT): {jd_tt:.6f}, JD(UT1): {jd_ut:.6f}")
        JD(TT): 2451545.000743, JD(UT1): 2451545.000004
    """
    ts = get_timescale()

    if calendar == SE_JUL_CAL:
        # Convert Julian calendar date to Gregorian for Skyfield
        # Skyfield's utc() expects proleptic Gregorian calendar
        decimal_hour = hour + minute / 60.0 + second / 3600.0
        jd = swe_julday(year, month, day, decimal_hour, SE_JUL_CAL)
        g_year, g_month, g_day, g_hour = swe_revjul(jd, SE_GREG_CAL)
        # Extract hour, minute, second from decimal hour
        g_minute_frac = (g_hour % 1) * 60
        g_second = (g_minute_frac % 1) * 60
        g_hour_int = int(g_hour)
        g_minute_int = int(g_minute_frac)
        t = ts.utc(g_year, g_month, g_day, g_hour_int, g_minute_int, g_second)
    else:
        # Gregorian calendar - use directly with Skyfield
        t = ts.utc(year, month, day, hour, minute, second)

    # Explicit float cast to satisfy type checker (Skyfield uses lazy reify decorator)
    return float(t.tt), float(t.ut1)


def jdet_to_utc(
    jd_et: float, calendar: int = SE_GREG_CAL
) -> tuple[int, int, int, int, int, float]:
    """
    Convert Julian Day in Ephemeris Time (TT/ET) to UTC date/time.

    This is the inverse of utc_to_jd(): it converts a Julian Day number
    in Terrestrial Time (TT, formerly known as Ephemeris Time/ET) back to
    a UTC calendar date and time. The conversion properly accounts for
    Delta-T and leap seconds.

    Args:
        jd_et: Julian Day number in TT (Terrestrial Time / Ephemeris Time)
        calendar: SE_GREG_CAL (1) for Gregorian, SE_JUL_CAL (0) for Julian

    Returns:
        tuple: (year, month, day, hour, minute, second) where:
            - year: Calendar year (negative for BCE)
            - month: Month (1-12)
            - day: Day of month (1-31)
            - hour: Hour (0-23)
            - minute: Minute (0-59)
            - second: Second (0.0-59.999..., or 60.x during leap second)

    Note:
        - TT (Terrestrial Time) is the modern successor to Ephemeris Time (ET)
        - TT = TAI + 32.184 seconds, where TAI is International Atomic Time
        - UTC may include leap seconds (second = 60) on certain dates
        - Delta-T (TT - UT1) is automatically applied using IERS data

    Example:
        >>> from libephemeris import jdet_to_utc, utc_to_jd, SE_GREG_CAL
        >>> # Convert J2000.0 epoch (JD 2451545.0 TT) to UTC
        >>> year, month, day, hour, minute, second = jdet_to_utc(2451545.0)
        >>> print(f"{year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second:05.2f}")
        2000-01-01 11:58:55.82
    """
    ts = get_timescale()

    # Create a Skyfield Time object from TT Julian Day
    t = ts.tt_jd(jd_et)

    # Get UTC components from Skyfield (handles leap seconds automatically)
    # The .utc attribute returns a tuple: (year, month, day, hour, minute, second)
    # We cast to Any to work around Skyfield's reify decorator type annotation issues
    utc_data: Any = t.utc
    g_year = int(utc_data[0])
    g_month = int(utc_data[1])
    g_day = int(utc_data[2])
    g_hour = int(utc_data[3])
    g_minute = int(utc_data[4])
    g_second = float(utc_data[5])

    if calendar == SE_JUL_CAL:
        # Convert Gregorian date to Julian calendar
        # First compute the JD for this Gregorian date
        decimal_hour = g_hour + g_minute / 60.0 + g_second / 3600.0
        jd_greg = swe_julday(g_year, g_month, g_day, decimal_hour, SE_GREG_CAL)
        # Convert to Julian calendar
        j_year, j_month, j_day, j_decimal_hour = swe_revjul(jd_greg, SE_JUL_CAL)
        # Extract time components from decimal hour
        j_hour = int(j_decimal_hour)
        j_minute_frac = (j_decimal_hour - j_hour) * 60.0
        j_minute = int(j_minute_frac)
        j_second = (j_minute_frac - j_minute) * 60.0
        return j_year, j_month, j_day, j_hour, j_minute, j_second

    return g_year, g_month, g_day, g_hour, g_minute, g_second


def jdut1_to_utc(
    jd_ut1: float, calendar: int = SE_GREG_CAL
) -> tuple[int, int, int, int, int, float]:
    """
    Convert Julian Day in UT1 (Universal Time) to UTC date/time.

    Converts a Julian Day number in UT1 back to a UTC calendar date and time.
    The difference between UT1 and UTC is always less than 0.9 seconds by
    definition (maintained by adding leap seconds to UTC).

    Args:
        jd_ut1: Julian Day number in UT1 (Universal Time)
        calendar: SE_GREG_CAL (1) for Gregorian, SE_JUL_CAL (0) for Julian

    Returns:
        tuple: (year, month, day, hour, minute, second) where:
            - year: Calendar year (negative for BCE)
            - month: Month (1-12)
            - day: Day of month (1-31)
            - hour: Hour (0-23)
            - minute: Minute (0-59)
            - second: Second (0.0-59.999..., or 60.x during leap second)

    Note:
        - UT1 is based on Earth's rotation and is not perfectly uniform
        - UTC is atomic time adjusted to stay within 0.9s of UT1
        - The difference DUT1 = UT1 - UTC is published by IERS
        - For high-precision astronomical work, this difference matters

    Example:
        >>> from libephemeris import jdut1_to_utc, utc_to_jd, SE_GREG_CAL
        >>> # Get JD(UT1) for a date, then convert back
        >>> jd_tt, jd_ut1 = utc_to_jd(2020, 6, 15, 14, 30, 0.0)
        >>> year, month, day, hour, minute, second = jdut1_to_utc(jd_ut1)
        >>> print(f"{year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second:05.2f}")
        2020-06-15 14:30:00.00
    """
    ts = get_timescale()

    # Create a Skyfield Time object from UT1 Julian Day
    t = ts.ut1_jd(jd_ut1)

    # Get UTC components from Skyfield (handles leap seconds automatically)
    # The .utc attribute returns a tuple: (year, month, day, hour, minute, second)
    # We cast to Any to work around Skyfield's reify decorator type annotation issues
    utc_data: Any = t.utc
    g_year = int(utc_data[0])
    g_month = int(utc_data[1])
    g_day = int(utc_data[2])
    g_hour = int(utc_data[3])
    g_minute = int(utc_data[4])
    g_second = float(utc_data[5])

    if calendar == SE_JUL_CAL:
        # Convert Gregorian date to Julian calendar
        # First compute the JD for this Gregorian date
        decimal_hour = g_hour + g_minute / 60.0 + g_second / 3600.0
        jd_greg = swe_julday(g_year, g_month, g_day, decimal_hour, SE_GREG_CAL)
        # Convert to Julian calendar
        j_year, j_month, j_day, j_decimal_hour = swe_revjul(jd_greg, SE_JUL_CAL)
        # Extract time components from decimal hour
        j_hour = int(j_decimal_hour)
        j_minute_frac = (j_decimal_hour - j_hour) * 60.0
        j_minute = int(j_minute_frac)
        j_second = (j_minute_frac - j_minute) * 60.0
        return j_year, j_month, j_day, j_hour, j_minute, j_second

    return g_year, g_month, g_day, g_hour, g_minute, g_second


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


def time_equ(jd: float) -> float:
    """
    Calculate the Equation of Time.

    The Equation of Time is the difference between apparent solar time
    (sundial time) and mean solar time (clock time). It arises from:
    1. Earth's orbital eccentricity (elliptical orbit)
    2. The obliquity of the ecliptic (Earth's axial tilt)

    Args:
        jd: Julian Day number in UT

    Returns:
        float: Equation of Time in days (multiply by 1440 for minutes)
               Positive values mean the sundial is ahead of the clock.

    Note:
        The equation of time varies from approximately -14 to +16 minutes
        throughout the year. The extremes occur around:
        - Early November: ~+16 minutes (sundial ahead)
        - Early February: ~-14 minutes (sundial behind)
        - Mid-April, mid-June, early September, late December: ~0 minutes

    Example:
        >>> from libephemeris import time_equ, swe_julday
        >>> # Calculate equation of time for Jan 1, 2000
        >>> jd = swe_julday(2000, 1, 1, 12.0)
        >>> eot = time_equ(jd)
        >>> eot_minutes = eot * 1440
        >>> print(f"Equation of Time: {eot_minutes:.2f} minutes")
        Equation of Time: -3.05 minutes
    """
    from typing import Any

    from .state import get_planets, get_timescale

    ts = get_timescale()
    planets = get_planets()

    # Create time object
    t = ts.ut1_jd(jd)

    # Get Earth and Sun
    earth = planets["earth"]
    sun = planets["sun"]

    # Calculate apparent position of Sun from Earth
    # Type ignore needed for Skyfield's lazy reify decorator
    astrometric: Any = earth.at(t).observe(sun)
    apparent = astrometric.apparent()

    # Get apparent right ascension in degrees
    ra, _, _ = apparent.radec()
    apparent_ra_deg = float(ra.hours) * 15.0  # Convert hours to degrees

    # Calculate Sun's mean longitude using algorithm from
    # Meeus "Astronomical Algorithms" Chapter 28
    # T is Julian centuries from J2000.0
    T = (jd - 2451545.0) / 36525.0

    # Mean longitude of the Sun (in degrees), from Meeus Ch. 25
    L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T**2

    # Normalize to 0-360 degrees
    L0 = L0 % 360.0
    if L0 < 0:
        L0 += 360.0

    # Normalize apparent RA to 0-360 degrees
    apparent_ra_deg = apparent_ra_deg % 360.0
    if apparent_ra_deg < 0:
        apparent_ra_deg += 360.0

    # Equation of time: E = L0 - 0.0057183° - α (all in degrees)
    # The constant 0.0057183° corrects for aberration and nutation
    # already included in apparent position, so we skip it.
    # E = Mean solar longitude - Apparent right ascension (in degrees)
    eot_deg = L0 - apparent_ra_deg

    # Normalize to ±180 degrees
    while eot_deg > 180.0:
        eot_deg -= 360.0
    while eot_deg < -180.0:
        eot_deg += 360.0

    # Convert degrees to time (360° = 24h = 1 day)
    eot_days = eot_deg / 360.0

    return eot_days


def utc_time_zone(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    timezone_offset: float,
) -> tuple[int, int, int, int, int, float]:
    """
    Apply a timezone offset to a UTC date/time and return the local date/time.

    Converts a UTC date/time to local time by adding the specified timezone
    offset. Handles all date boundary crossings (day, month, year) correctly.

    Args:
        year: Calendar year (negative for BCE)
        month: Month (1-12)
        day: Day of month (1-31)
        hour: Hour (0-23)
        minute: Minute (0-59)
        second: Second (0.0-59.999...)
        timezone_offset: Timezone offset in hours from UTC.
            Positive values for timezones east of UTC (e.g., +1 for CET, +9 for JST)
            Negative values for timezones west of UTC (e.g., -5 for EST, -8 for PST)

    Returns:
        tuple: (year, month, day, hour, minute, second) in local time where:
            - year: Calendar year
            - month: Month (1-12)
            - day: Day of month (1-31)
            - hour: Hour (0-23)
            - minute: Minute (0-59)
            - second: Second (0.0-59.999...)

    Example:
        >>> from libephemeris import utc_time_zone
        >>> # Convert 2024-01-15 10:30:00 UTC to CET (UTC+1)
        >>> utc_time_zone(2024, 1, 15, 10, 30, 0.0, 1)
        (2024, 1, 15, 11, 30, 0.0)
        >>> # Convert 2024-01-15 02:00:00 UTC to EST (UTC-5)
        >>> utc_time_zone(2024, 1, 15, 2, 0, 0.0, -5)
        (2024, 1, 14, 21, 0, 0.0)
    """
    # Convert UTC time to decimal hours
    decimal_hour = hour + minute / 60.0 + second / 3600.0

    # Convert to Julian Day
    jd = swe_julday(year, month, day, decimal_hour, SE_GREG_CAL)

    # Add timezone offset (convert hours to days)
    jd_local = jd + timezone_offset / 24.0

    # Convert back to calendar date
    local_year, local_month, local_day, local_decimal_hour = swe_revjul(
        jd_local, SE_GREG_CAL
    )

    # Extract time components from decimal hour with rounding to avoid
    # floating-point precision issues (e.g., 11.4999999 -> 11.5)
    # Round to millisecond precision (3 decimal places in seconds)
    # This is sufficient for timezone conversions and avoids floating-point errors
    total_seconds = local_decimal_hour * 3600.0
    total_seconds = round(total_seconds, 3)

    local_hour = int(total_seconds // 3600)
    remaining = total_seconds - local_hour * 3600
    local_minute = int(remaining // 60)
    local_second = remaining - local_minute * 60

    # Handle edge case where rounding pushes us to 60 seconds
    if local_second >= 60.0:
        local_second -= 60.0
        local_minute += 1
    if local_minute >= 60:
        local_minute -= 60
        local_hour += 1
    if local_hour >= 24:
        local_hour -= 24
        # Day already handled by swe_revjul, this shouldn't happen
        # but included for safety

    return local_year, local_month, local_day, local_hour, local_minute, local_second
