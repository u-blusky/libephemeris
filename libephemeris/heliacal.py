"""
Heliacal event calculations for libephemeris.

Calculates heliacal risings and settings, visual limiting magnitude,
and related heliacal phenomena for celestial bodies.

Functions:
- heliacal_ut: Find heliacal rising/setting time for a body
- swe_heliacal_ut: Alias for heliacal_ut
- heliacal_pheno_ut: Calculate detailed heliacal phenomena
- swe_heliacal_pheno_ut: Alias for heliacal_pheno_ut
- vis_limit_mag: Calculate visual limiting magnitude
- swe_vis_limit_mag: Alias for vis_limit_mag

Historical Note:
    Heliacal risings were crucial for ancient calendars. The heliacal
    rising of Sirius marked the Egyptian new year and predicted the
    Nile flood. Babylonians used heliacal events to track planetary
    positions without modern instruments.

References:
    - Swiss Ephemeris documentation section 5.1
    - Schoch "Planets in Mesopotamian Astral Science"
    - Ptolemy's criteria for heliacal visibility
    - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
"""

import math
from typing import Tuple

from .constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_FIXSTAR_OFFSET,
    SEFLG_SPEED,
    SEFLG_SWIEPH,
)

# Inner planets (orbit inside Earth's orbit)
# These have both inferior conjunction (between Earth and Sun)
# and superior conjunction (behind the Sun)
INNER_PLANETS = {SE_MERCURY, SE_VENUS}

# Outer planets (orbit outside Earth's orbit)
# These only have one type of conjunction (behind the Sun)
# SE_EVENING_FIRST and SE_MORNING_LAST are not applicable to these


def is_fixed_star(body: int) -> bool:
    """
    Check if a body ID corresponds to a fixed star.

    Args:
        body: Body ID constant

    Returns:
        True if the body is a fixed star (ID >= SE_FIXSTAR_OFFSET), False otherwise
    """
    return body >= SE_FIXSTAR_OFFSET


def _get_star_magnitude(star_id: int) -> float:
    """
    Get the visual magnitude of a fixed star from the catalog.

    Args:
        star_id: Star ID (SE_* constant, e.g., SE_SIRIUS)

    Returns:
        Visual magnitude of the star. Brighter stars have lower/negative values.
        Returns 6.0 (faint) if star not found.
    """
    from .fixed_stars import STAR_CATALOG

    for entry in STAR_CATALOG:
        if entry.id == star_id:
            return entry.magnitude
    return 6.0  # Default to faint star if not found


def is_inner_planet(body: int) -> bool:
    """
    Check if a body is an inner planet (Mercury or Venus).

    Inner planets orbit inside Earth's orbit and have both inferior
    and superior conjunctions with the Sun. They can appear as both
    morning and evening stars.

    Args:
        body: Planet ID constant (SE_MERCURY, SE_VENUS, etc.)

    Returns:
        True if the body is Mercury or Venus, False otherwise
    """
    return body in INNER_PLANETS


def heliacal_ut(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 0.5,
    body: int = SE_SUN,
    event_type: int = 1,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[float, int]:
    """
    Calculate heliacal rising or setting time for a celestial body.

    Heliacal events are the first/last visibility of a celestial body
    at dawn or dusk. These were fundamental for ancient calendars:
    - Heliacal rising: First morning visibility after a period of invisibility
    - Heliacal setting: Last evening visibility before becoming invisible

    Args:
        jd_start: Julian Day (UT) to start search from
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        pressure: Atmospheric pressure in mbar/hPa for refraction (default 1013.25)
        temperature: Temperature in Celsius for refraction (default 15)
        humidity: Relative humidity 0.0-1.0 for atmospheric extinction (default 0.5)
        body: Planet/body ID (SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN,
              or fixed stars). Note: SE_SUN and SE_MOON are not valid for heliacal events.
        event_type: Type of heliacal event:
            - SE_HELIACAL_RISING (1): Morning first visibility (heliacal rising)
            - SE_HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
            - SE_EVENING_FIRST (3): First evening visibility (after superior conjunction)
              Note: Only valid for inner planets (Mercury, Venus)
            - SE_MORNING_LAST (4): Last morning visibility (before superior conjunction)
              Note: Only valid for inner planets (Mercury, Venus)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - jd_event: Julian Day (UT) of the heliacal event, or 0.0 if not found
            - retflag: Return flag (event_type on success, negative on error)

    Raises:
        ValueError: If invalid body ID, event_type, or if SE_EVENING_FIRST/SE_MORNING_LAST
            is used with an outer planet (Mars, Jupiter, Saturn, etc.)

    Algorithm:
        The algorithm searches for the moment when:
        1. The body is at a specific altitude above the horizon (arcus visionis)
        2. The Sun is at twilight position (typically -6 to -12 below horizon)
        3. The body's apparent magnitude is brighter than the sky's limiting magnitude

        For heliacal rising (morning first):
        - Search forward for when the body first becomes visible at dawn
        - Body must be above horizon while Sun is still below
        - Sky must be dark enough for the body to be seen

        For heliacal setting (evening last):
        - Search forward for when the body is last visible at dusk
        - Body must be above horizon while Sun is setting
        - Sky brightness must not overwhelm the body's light

    Historical Note:
        Heliacal risings were crucial for ancient calendars. The heliacal
        rising of Sirius marked the Egyptian new year and predicted the
        Nile flood. Babylonians used heliacal events to track planetary
        positions without modern instruments.

    Example:
        >>> from libephemeris import julday, heliacal_ut, SE_VENUS, SE_HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Find next heliacal rising of Venus from Rome
        >>> jd_event, flag = heliacal_ut(jd, 41.9, 12.5, body=SE_VENUS,
        ...                              event_type=SE_HELIACAL_RISING)
        >>> print(f"Heliacal rising at JD {jd_event:.5f}")

    References:
        - Swiss Ephemeris: swe_heliacal_ut()
        - Schoch "Planets in Mesopotamian Astral Science"
        - Ptolemy's criteria for heliacal visibility
    """
    from .constants import (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    )
    from .planets import _PLANET_MAP, swe_pheno_ut
    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Validate event type
    if event_type not in (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {event_type}. Use SE_HELIACAL_RISING, "
            "SE_HELIACAL_SETTING, SE_EVENING_FIRST, or SE_MORNING_LAST."
        )

    # SE_EVENING_FIRST and SE_MORNING_LAST are only valid for inner planets
    # (Mercury and Venus) because they relate to superior conjunction visibility.
    # Outer planets only have one type of conjunction and these events don't apply.
    # Fixed stars only have heliacal rising and setting.
    if event_type in (SE_EVENING_FIRST, SE_MORNING_LAST):
        if is_fixed_star(body):
            raise ValueError(
                "SE_EVENING_FIRST and SE_MORNING_LAST are not valid for fixed stars. "
                "For fixed stars, use SE_HELIACAL_RISING or SE_HELIACAL_SETTING."
            )
        if not is_inner_planet(body):
            raise ValueError(
                "SE_EVENING_FIRST and SE_MORNING_LAST are only valid for inner planets "
                "(Mercury, Venus). For outer planets, use SE_HELIACAL_RISING or "
                "SE_HELIACAL_SETTING."
            )

    # Sun and Moon are not valid for heliacal events
    if body == SE_SUN:
        raise ValueError("SE_SUN is not valid for heliacal calculations")
    if body == SE_MOON:
        raise ValueError("SE_MOON is not valid for heliacal calculations")

    # Check if this is a fixed star
    is_star = is_fixed_star(body)

    # Validate body - must be either a known planet or a fixed star
    if not is_star and body not in _PLANET_MAP:
        raise ValueError(f"illegal planet number {body}.")

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    sun = eph["sun"]
    earth = eph["earth"]

    # Get target - either planet or star
    target = None
    star_object = None
    star_magnitude = 0.0  # Default, will be set for stars

    if is_star:
        # For fixed stars, create a Skyfield Star object
        from skyfield.api import Star
        from .fixed_stars import STAR_CATALOG

        star_magnitude = _get_star_magnitude(body)

        # Find star data from catalog
        star_data = None
        for entry in STAR_CATALOG:
            if entry.id == body:
                star_data = entry.data
                break

        if star_data is None:
            raise ValueError(f"Star ID {body} not found in catalog")

        # Create Skyfield Star object for position calculations
        star_object = Star(
            ra_hours=star_data.ra_j2000 / 15.0,
            dec_degrees=star_data.dec_j2000,
            ra_mas_per_year=star_data.pm_ra * 1000.0,
            dec_mas_per_year=star_data.pm_dec * 1000.0,
        )
    else:
        # For planets, get the target from ephemeris
        target_name = _PLANET_MAP[body]
        from .planets import _PLANET_FALLBACK

        try:
            target = eph[target_name]
        except KeyError:
            if target_name in _PLANET_FALLBACK:
                target = eph[_PLANET_FALLBACK[target_name]]
            else:
                raise

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Note: Arcus visionis (minimum altitude for visibility) and
    # sun altitude thresholds are used in the visibility check functions.
    # Typical values: arcus visionis ~7 for bright objects, sun altitude ~-8.

    def _get_altitudes(jd: float) -> Tuple[float, float, float]:
        """Get Sun altitude, body altitude, and body azimuth at given JD."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer

        # Sun position
        sun_app = observer_at.at(t).observe(sun).apparent()
        sun_alt, _, _ = sun_app.altaz()

        # Body position (handle both planets and stars)
        if is_star and star_object is not None:
            body_app = observer_at.at(t).observe(star_object).apparent()
        else:
            body_app = observer_at.at(t).observe(target).apparent()
        body_alt, body_az, _ = body_app.altaz()

        return sun_alt.degrees, body_alt.degrees, body_az.degrees

    def _get_elongation(jd: float) -> float:
        """Get the elongation of body from Sun in degrees."""
        if not is_star:
            try:
                pheno, _ = swe_pheno_ut(jd, body, flags)
                return pheno[2]  # Elongation
            except Exception:
                pass

        # For stars or fallback: calculate elongation manually
        t = ts.ut1_jd(jd)
        sun_app = earth.at(t).observe(sun).apparent()
        if is_star and star_object is not None:
            body_app = earth.at(t).observe(star_object).apparent()
        else:
            body_app = earth.at(t).observe(target).apparent()
        return body_app.separation_from(sun_app).degrees

    def _get_body_magnitude(jd: float) -> float:
        """Get the visual magnitude of the body."""
        if is_star:
            # For fixed stars, return the catalog magnitude
            return star_magnitude
        try:
            pheno, _ = swe_pheno_ut(jd, body, flags)
            return pheno[4]  # Visual magnitude
        except Exception:
            return 0.0  # Default to bright magnitude

    def _calculate_limiting_magnitude(sun_alt: float, body_alt: float) -> float:
        """
        Calculate the limiting magnitude based on sky brightness.

        The limiting magnitude depends on:
        - Sun's altitude below horizon (darker = higher limit)
        - Body's altitude (extinction increases near horizon)
        - Atmospheric conditions (humidity, pressure)

        This is a simplified model based on Schaefer (1990).
        """
        # Sky brightness model (simplified)
        # When Sun is below -18, sky is fully dark (limiting mag ~6.5)
        # When Sun is at -6 (civil twilight), limiting mag is ~3
        # When Sun is at 0 (horizon), limiting mag is ~-2

        if sun_alt >= 0:
            return -2.0  # Daylight - essentially nothing visible
        elif sun_alt >= -6:
            # Civil twilight
            return -2.0 + (sun_alt + 6) * (-3.0 - (-2.0)) / (-6)
        elif sun_alt >= -12:
            # Nautical twilight
            return 3.0 + (sun_alt + 12) * (5.0 - 3.0) / 6
        elif sun_alt >= -18:
            # Astronomical twilight
            return 5.0 + (sun_alt + 18) * (6.5 - 5.0) / 6
        else:
            # Full darkness
            return 6.5

        # Atmospheric extinction correction (simplified)
        # Add extinction based on altitude (airmass)
        # At low altitudes, extinction can be 0.5-1.0 magnitudes

    def _is_body_visible(jd: float) -> Tuple[bool, float, float, float]:
        """
        Check if body is visible at given time.

        Returns: (is_visible, sun_alt, body_alt, elongation)
        """
        sun_alt, body_alt, body_az = _get_altitudes(jd)
        elongation = _get_elongation(jd)

        # Body must be above minimum altitude
        if body_alt < 0:
            return False, sun_alt, body_alt, elongation

        # Sun must be below horizon
        if sun_alt > 0:
            return False, sun_alt, body_alt, elongation

        # Check limiting magnitude vs body magnitude
        limiting_mag = _calculate_limiting_magnitude(sun_alt, body_alt)
        body_mag = _get_body_magnitude(jd)

        # Body is visible if its magnitude is brighter (lower) than limiting magnitude
        # Also require sufficient elongation from Sun
        min_elongation = 10.0  # Minimum elongation for visibility (degrees)

        is_visible = (body_mag <= limiting_mag) and (elongation >= min_elongation)

        return is_visible, sun_alt, body_alt, elongation

    def _find_twilight_time(jd: float, sun_target_alt: float, rising: bool) -> float:
        """
        Find when Sun crosses target altitude (morning or evening).

        Args:
            jd: Starting JD
            sun_target_alt: Target Sun altitude (negative for below horizon)
            rising: True for morning (Sun rising), False for evening (Sun setting)
        """
        # Search within one day
        for _ in range(50):
            sun_alt, _, _ = _get_altitudes(jd)

            # Check if we're at the right phase of day
            if rising:
                # Morning: looking for Sun rising through target altitude
                if abs(sun_alt - sun_target_alt) < 0.01:
                    return jd
            else:
                # Evening: looking for Sun setting through target altitude
                if abs(sun_alt - sun_target_alt) < 0.01:
                    return jd

            # Adjust time based on Sun's position
            # Sun moves ~15/hour = 0.25/minute = 360/day
            sun_rate = 360.0 / 1.0  # degrees per day
            diff = sun_target_alt - sun_alt

            # Crude estimate of time adjustment
            dt = diff / sun_rate

            # Limit step size
            dt = max(-0.1, min(0.1, dt))

            if abs(dt) < 1e-6:
                return jd

            jd += dt

        return jd

    def _search_heliacal_rising(jd_start: float) -> float:
        """
        Search for heliacal rising (morning first visibility).

        The body becomes visible in the morning before sunrise after
        a period of being hidden in the Sun's glare.
        """
        max_days = 400  # Search up to more than a year

        # Use coarse search first (check once per day at approximate twilight)
        for day in range(max_days):
            jd_day = jd_start + day

            # Check at approximately 5 AM local time (rough morning twilight estimate)
            # This is a coarse check - we'll refine if we find visibility
            for hour in [4, 5, 6]:  # Check a few morning hours
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                # Look for conditions: Sun below horizon but not too deep, body above
                if -15 < sun_alt < -3 and body_alt > 1:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)

                    if visible:
                        # Found a potential heliacal rising
                        # Refine the time using binary search
                        return _refine_heliacal_time(jd_check, is_morning=True)

        return 0.0  # Not found

    def _search_heliacal_setting(jd_start: float) -> float:
        """
        Search for heliacal setting (evening last visibility).

        The body is last visible in the evening after sunset before
        becoming hidden in the Sun's glare.
        """
        max_days = 400

        # First, find when the body is visible in the evening
        first_visible_jd = 0.0

        for day in range(max_days):
            jd_day = jd_start + day

            # Check at a few evening hours (coarse search)
            for hour in [18, 19, 20]:
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                # Look for conditions: Sun below horizon, body above
                if -15 < sun_alt < -3 and body_alt > 1:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)

                    if visible:
                        first_visible_jd = jd_check
                        break

            if first_visible_jd > 0:
                break

        if first_visible_jd == 0:
            return 0.0

        # Now find the LAST day of visibility
        last_visible_jd = first_visible_jd
        start_day = int(first_visible_jd - jd_start)

        for day in range(start_day, start_day + 120):  # Check up to 120 days ahead
            jd_day = jd_start + day
            found_visible = False

            for hour in [18, 19, 20]:
                jd_check = jd_day + hour / 24.0
                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                if -15 < sun_alt < -3 and body_alt > 1:
                    visible, _, _, _ = _is_body_visible(jd_check)
                    if visible:
                        last_visible_jd = jd_check
                        found_visible = True
                        break

            if not found_visible and day > start_day:
                # No visibility found today - return the last visible day
                return _refine_heliacal_time(last_visible_jd, is_morning=False)

        # If still visible after 120 days, return the last checked
        return _refine_heliacal_time(last_visible_jd, is_morning=False)

    def _search_evening_first(jd_start: float) -> float:
        """
        Search for evening first visibility (after superior conjunction).

        The body appears in the evening sky for the first time after
        passing behind the Sun (superior conjunction for inferior planets,
        or conjunction for superior planets).
        """
        max_days = 400
        was_invisible = False

        for day in range(max_days):
            jd_day = jd_start + day

            # Check at a few evening hours (coarse search)
            for hour in [18, 19, 20]:
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                if -15 < sun_alt < -3:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)
                    if not visible and body_alt < 5:
                        was_invisible = True
                    elif visible and was_invisible:
                        # First evening visibility after being invisible
                        return _refine_heliacal_time(jd_check, is_morning=False)

        return 0.0

    def _search_morning_last(jd_start: float) -> float:
        """
        Search for morning last visibility (before superior conjunction).

        The body is last visible in the morning sky before passing
        behind the Sun.
        """
        max_days = 400
        last_visible_jd = 0.0
        found_visible = False

        for day in range(max_days):
            jd_day = jd_start + day

            # Check at a few morning hours (coarse search)
            for hour in [4, 5, 6]:
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                if -15 < sun_alt < -3:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)

                    if visible:
                        last_visible_jd = jd_check
                        found_visible = True
                    elif found_visible and not visible:
                        # Found the transition to invisibility
                        return _refine_heliacal_time(last_visible_jd, is_morning=True)

        return last_visible_jd if last_visible_jd > 0 else 0.0

    def _refine_heliacal_time(jd_approx: float, is_morning: bool) -> float:
        """
        Refine the heliacal event time using binary search.

        Find the exact moment when the body becomes just visible/invisible.
        """
        # Use binary search to find the transition point
        jd_low = jd_approx - 0.1  # ~2.4 hours before
        jd_high = jd_approx + 0.1  # ~2.4 hours after

        for _ in range(30):  # ~30 iterations gives very high precision
            jd_mid = (jd_low + jd_high) / 2

            visible, sun_alt, body_alt, elong = _is_body_visible(jd_mid)

            # For heliacal rising: looking for first visibility
            # For heliacal setting: looking for last visibility
            if is_morning:
                # Morning: visibility increases as time progresses (Sun rises)
                # But we want first visibility, so search backwards
                if visible:
                    jd_high = jd_mid  # Look earlier
                else:
                    jd_low = jd_mid  # Look later
            else:
                # Evening: visibility decreases as time progresses (sky darkens)
                # For last visibility, we want the last moment still visible
                if visible:
                    jd_low = jd_mid  # Look later for last visibility
                else:
                    jd_high = jd_mid  # Look earlier

            if jd_high - jd_low < 1e-6:  # ~0.1 second precision
                break

        return (jd_low + jd_high) / 2

    # Main search logic based on event type
    if event_type == SE_HELIACAL_RISING:
        jd_event = _search_heliacal_rising(jd_start)
    elif event_type == SE_HELIACAL_SETTING:
        jd_event = _search_heliacal_setting(jd_start)
    elif event_type == SE_EVENING_FIRST:
        jd_event = _search_evening_first(jd_start)
    elif event_type == SE_MORNING_LAST:
        jd_event = _search_morning_last(jd_start)
    else:
        jd_event = 0.0

    if jd_event > 0:
        return jd_event, event_type
    else:
        return 0.0, -1  # Not found


def swe_heliacal_ut(
    jd_start: float,
    geopos: tuple,
    datm: tuple,
    dobs: tuple,
    object_name: str,
    event_type: int,
    hel_flag: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], int]:
    """
    Calculate heliacal rising or setting time for a celestial body.

    This function implements the Swiss Ephemeris swe_heliacal_ut() API,
    finding the Julian day of the next heliacal phenomenon after a given
    start date. It works between geographic latitudes 60S - 60N.

    Args:
        jd_start: Julian Day (UT) to start search from for the heliacal event.
        geopos: Geographic position as a sequence of at least 3 values:
            - [0]: Geographic longitude in degrees (east positive)
            - [1]: Geographic latitude in degrees (north positive)
            - [2]: Altitude above sea level in meters (eye height)
        datm: Atmospheric conditions as a sequence of at least 4 values:
            - [0]: Atmospheric pressure in mbar/hPa (default: 1013.25)
            - [1]: Atmospheric temperature in degrees Celsius (default: 15)
            - [2]: Relative humidity in percent (0-100) (default: 40)
            - [3]: If >= 1: Meteorological Range in km
                   If 0 < value < 1: Total atmospheric coefficient (ktot)
                   If 0: Calculate ktot from other parameters
        dobs: Observer description as a sequence of at least 6 values:
            - [0]: Age of observer in years (default: 36)
            - [1]: Snellen ratio of observer's eyes (default: 1.0 = normal)
            - [2]: 0 = monocular, 1 = binocular (used if SE_HELFLAG_OPTICAL_PARAMS)
            - [3]: Telescope magnification (0 = naked eye)
            - [4]: Optical aperture (telescope diameter) in mm
            - [5]: Optical transmission coefficient
        object_name: Name of the celestial body. Can be:
            - Planet name: "Sun", "Moon", "Mercury", "Venus", "Mars",
              "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"
            - Fixed star name: "Sirius", "Aldebaran", "Regulus", etc.
            Note: Sun and Moon are not valid for heliacal calculations
              and will raise a ValueError.
        event_type: Type of heliacal event:
            - SE_HELIACAL_RISING (1): Morning first visibility (heliacal rising)
              Exists for all visible planets and stars.
            - SE_HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
              Exists for all visible planets and stars.
            - SE_EVENING_FIRST (3): First evening visibility after superior
              conjunction. Only valid for inner planets (Mercury, Venus).
            - SE_MORNING_LAST (4): Last morning visibility before superior
              conjunction. Only valid for inner planets (Mercury, Venus).
        hel_flag: Calculation flags (bitmap). Contains ephemeris flags like
            SEFLG_SWIEPH, plus heliacal-specific flags:
            - SE_HELFLAG_OPTICAL_PARAMS (512): Use optical instrument parameters
              from dobs[2-5]. Without this flag, those values are ignored.
            - SE_HELFLAG_NO_DETAILS (1024): Provide date only, skip visibility
              start/optimum/end details. Makes calculation faster.
            - SE_HELFLAG_VISLIM_DARK (4096): Function behaves as if Sun at nadir.
            - SE_HELFLAG_VISLIM_NOMOON (8192): Exclude Moon's brightness
              contribution (useful for calculating epoch heliacal dates).

    Returns:
        Tuple containing:
            - dret: Tuple of 50 floats with heliacal event data:
                - [0]: Start of visibility (Julian day number)
                - [1]: Optimum visibility (Julian day number),
                       zero if SE_HELFLAG_NO_DETAILS is set
                - [2]: End of visibility (Julian day number),
                       zero if SE_HELFLAG_NO_DETAILS is set
                - [3-49]: Reserved for future use (currently 0.0)
            - retflag: Return flag (event_type on success, negative on error)

    Raises:
        ValueError: If invalid object_name, body ID, or event_type.
            Also raised if object_name is "Sun" or "Moon" which are not
            valid for heliacal calculations.

    Algorithm:
        The function searches for the moment when:
        1. The body is at a specific altitude above the horizon (arcus visionis)
        2. The Sun is at twilight position (typically -6 to -12 degrees below)
        3. The body's apparent magnitude is brighter than sky's limiting magnitude

        For heliacal rising (morning first):
        - Search forward for when body first becomes visible at dawn
        - Body must be above horizon while Sun is still below
        - Sky must be dark enough for the body to be seen

        For heliacal setting (evening last):
        - Search forward for when body is last visible at dusk
        - Body must be above horizon while Sun is setting
        - Sky brightness must not overwhelm the body's light

    Historical Note:
        Heliacal risings were crucial for ancient calendars. The heliacal
        rising of Sirius marked the Egyptian new year and predicted the
        Nile flood. Babylonians used heliacal events to track planetary
        positions without modern instruments.

    Example:
        >>> from libephemeris import julday, swe_heliacal_ut
        >>> from libephemeris import SE_HELIACAL_RISING, SEFLG_SWIEPH
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Geographic position: Rome (lon, lat, altitude)
        >>> geopos = (12.5, 41.9, 0)
        >>> # Atmospheric conditions: standard
        >>> datm = (1013.25, 15.0, 40.0, 0.0)
        >>> # Observer: age 36, normal vision
        >>> dobs = (36.0, 1.0, 0, 0, 0, 0)
        >>> # Find next heliacal rising of Venus
        >>> dret, flag = swe_heliacal_ut(jd, geopos, datm, dobs,
        ...                               "Venus", SE_HELIACAL_RISING)
        >>> if flag > 0:
        ...     print(f"Venus heliacal rising at JD {dret[0]:.5f}")

    See Also:
        - heliacal_ut: Internal function using planet ID instead of name
        - heliacal_pheno_ut: Detailed heliacal phenomena calculation
        - vis_limit_mag: Calculate limiting visual magnitude

    References:
        - Swiss Ephemeris documentation section 8.17
        - Schoch "Planets in Mesopotamian Astral Science"
        - Ptolemy's criteria for heliacal visibility
    """
    from .constants import (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
        SE_HELFLAG_NO_DETAILS,
    )

    # Validate event type
    if event_type not in (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {event_type}. Use SE_HELIACAL_RISING, "
            "SE_HELIACAL_SETTING, SE_EVENING_FIRST, or SE_MORNING_LAST."
        )

    # Parse geopos
    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    altitude = geopos[2] if len(geopos) > 2 else 0.0

    # Parse datm with defaults per Swiss Ephemeris documentation
    pressure = datm[0] if len(datm) > 0 and datm[0] > 0 else 1013.25
    temperature = datm[1] if len(datm) > 1 else 15.0
    humidity_pct = datm[2] if len(datm) > 2 else 40.0
    # datm[3] is meteorological range / ktot, handled internally

    # Convert humidity from percent to 0-1 range for internal use
    humidity = humidity_pct / 100.0 if humidity_pct > 1.0 else humidity_pct

    # Parse object_name to get body ID
    body_id = _parse_object_name(object_name)

    # Call the internal heliacal_ut function
    jd_event, retflag = heliacal_ut(
        jd_start=jd_start,
        lat=lat,
        lon=lon,
        altitude=altitude,
        pressure=pressure,
        temperature=temperature,
        humidity=humidity,
        body=body_id,
        event_type=event_type,
        flags=hel_flag,
    )

    # Build the result array (50 elements as per Swiss Ephemeris)
    dret = [0.0] * 50

    if jd_event > 0:
        dret[0] = jd_event  # Start visibility

        # Calculate optimum and end visibility if details requested
        if not (hel_flag & SE_HELFLAG_NO_DETAILS):
            # For detailed calculation, we estimate optimum and end times
            # based on typical visibility window durations
            # This is an approximation; full implementation would require
            # more complex calculations
            dret[1] = jd_event  # Optimum (same as start for now)
            dret[2] = jd_event  # End (same as start for now)

    return tuple(dret), retflag


def _parse_object_name(object_name: str) -> int:
    """
    Parse an object name string to get the corresponding body ID.

    Args:
        object_name: Name of planet or fixed star (e.g., "Venus", "Sirius")

    Returns:
        Body ID (SE_* constant) for planets

    Raises:
        ValueError: If object name is not recognized or not valid for heliacal
    """
    # Import planet constants
    from .constants import (
        SE_SUN,
        SE_MOON,
        SE_MERCURY,
        SE_VENUS,
        SE_MARS,
        SE_JUPITER,
        SE_SATURN,
        SE_URANUS,
        SE_NEPTUNE,
        SE_PLUTO,
    )

    # Normalize name
    name_upper = object_name.upper().strip()

    # Planet name to ID mapping
    planet_map = {
        "SUN": SE_SUN,
        "MOON": SE_MOON,
        "MERCURY": SE_MERCURY,
        "VENUS": SE_VENUS,
        "MARS": SE_MARS,
        "JUPITER": SE_JUPITER,
        "SATURN": SE_SATURN,
        "URANUS": SE_URANUS,
        "NEPTUNE": SE_NEPTUNE,
        "PLUTO": SE_PLUTO,
    }

    # Check if it's a planet
    if name_upper in planet_map:
        body_id = planet_map[name_upper]

        # Sun and Moon are not valid for heliacal calculations
        if body_id == SE_SUN:
            raise ValueError("Sun is not valid for heliacal calculations")
        if body_id == SE_MOON:
            raise ValueError("Moon is not valid for heliacal calculations")

        return body_id

    # Try to parse as planet number
    try:
        body_id = int(object_name)
        if body_id == SE_SUN:
            raise ValueError("Sun is not valid for heliacal calculations")
        if body_id == SE_MOON:
            raise ValueError("Moon is not valid for heliacal calculations")
        return body_id
    except ValueError:
        pass

    # Try to resolve as a fixed star name
    from .fixed_stars import resolve_star_name

    star_id = resolve_star_name(object_name)
    if star_id is not None:
        return star_id

    # Object not recognized
    raise ValueError(
        f"Object '{object_name}' not recognized. "
        "Use planet names (Mercury, Venus, Mars, Jupiter, Saturn, etc.), "
        "planet IDs (2-9 for Mercury-Pluto), or fixed star names "
        "(Sirius, Regulus, Aldebaran, etc.)."
    )


def heliacal_pheno_ut(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 0.5,
    body: int = SE_SUN,
    event_type: int = 1,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[Tuple[float, ...], int]:
    """
    Provides data relevant for the calculation of heliacal risings and settings.

    This function calculates detailed phenomena associated with heliacal events,
    including altitudes, azimuths, arcus visionis, magnitude, visibility times,
    and other parameters used in heliacal visibility calculations.

    Args:
        jd: Julian Day (UT) for the calculation
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        pressure: Atmospheric pressure in mbar/hPa for refraction (default 1013.25)
        temperature: Temperature in Celsius for refraction (default 15)
        humidity: Relative humidity 0.0-1.0 for atmospheric extinction (default 0.5)
        body: Planet/body ID (SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN, etc.)
        event_type: Type of heliacal event:
            - SE_HELIACAL_RISING (1): Morning first visibility (heliacal rising)
            - SE_HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
            - SE_EVENING_FIRST (3): First evening visibility (after superior conjunction)
            - SE_MORNING_LAST (4): Last morning visibility (before superior conjunction)
        flags: Calculation flags (SEFLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - dret: Tuple of 50 floats with heliacal phenomena data:
                - 0: AltO [deg] topocentric altitude of object (unrefracted)
                - 1: AppAltO [deg] apparent altitude of object (refracted)
                - 2: GeoAltO [deg] geocentric altitude of object
                - 3: AziO [deg] azimuth of object
                - 4: AltS [deg] topocentric altitude of Sun
                - 5: AziS [deg] azimuth of Sun
                - 6: TAVact [deg] actual topocentric arcus visionis
                - 7: ARCVact [deg] actual (geocentric) arcus visionis
                - 8: DAZact [deg] actual difference between object's and sun's azimuth
                - 9: ARCLact [deg] actual longitude difference between object and sun
                - 10: kact [-] extinction coefficient
                - 11: minTAV [deg] smallest topocentric arcus visionis
                - 12: TfirstVR [JDN] first time object is visible, according to VR
                - 13: TbVR [JDN] optimum time the object is visible, according to VR
                - 14: TlastVR [JDN] last time object is visible, according to VR
                - 15: TbYallop [JDN] best time the object is visible, according to Yallop
                - 16: WMoon [deg] crescent width of Moon
                - 17: qYal [-] q-test value of Yallop
                - 18: qCrit [-] q-test criterion of Yallop
                - 19: ParO [deg] parallax of object
                - 20: Magn [-] magnitude of object
                - 21: RiseO [JDN] rise/set time of object
                - 22: RiseS [JDN] rise/set time of Sun
                - 23: Lag [JDN] rise/set time of object minus rise/set time of Sun
                - 24: TvisVR [JDN] visibility duration
                - 25: LMoon [deg] crescent length of Moon
                - 26-49: Reserved for future use
            - retflag: Return flag (flags on success, negative on error)

    Raises:
        ValueError: If invalid body ID or event_type

    Example:
        >>> from libephemeris import julday, heliacal_pheno_ut, SE_VENUS, SE_HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Get heliacal phenomena for Venus at Rome
        >>> dret, flag = heliacal_pheno_ut(jd, 41.9, 12.5, body=SE_VENUS,
        ...                                 event_type=SE_HELIACAL_RISING)
        >>> print(f"Object altitude: {dret[0]:.2f}, Sun altitude: {dret[4]:.2f}")

    References:
        - Swiss Ephemeris: swe_heliacal_pheno_ut()
        - Schoch "Planets in Mesopotamian Astral Science"
    """
    from .constants import (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    )
    from .planets import _PLANET_MAP, swe_pheno_ut
    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Validate event type
    if event_type not in (
        SE_HELIACAL_RISING,
        SE_HELIACAL_SETTING,
        SE_EVENING_FIRST,
        SE_MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {event_type}. Use SE_HELIACAL_RISING, "
            "SE_HELIACAL_SETTING, SE_EVENING_FIRST, or SE_MORNING_LAST."
        )

    # Check if this is a fixed star
    is_star = is_fixed_star(body)

    # Validate body - must be either a known planet or a fixed star
    if not is_star and body not in _PLANET_MAP:
        raise ValueError(f"illegal planet number {body}.")

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    sun = eph["sun"]
    earth = eph["earth"]

    # Get target - either planet or star
    target = None
    star_object = None
    star_magnitude = 0.0

    if is_star:
        # For fixed stars, create a Skyfield Star object
        from skyfield.api import Star
        from .fixed_stars import STAR_CATALOG

        star_magnitude = _get_star_magnitude(body)

        # Find star data from catalog
        star_data = None
        for entry in STAR_CATALOG:
            if entry.id == body:
                star_data = entry.data
                break

        if star_data is None:
            raise ValueError(f"Star ID {body} not found in catalog")

        # Create Skyfield Star object for position calculations
        star_object = Star(
            ra_hours=star_data.ra_j2000 / 15.0,
            dec_degrees=star_data.dec_j2000,
            ra_mas_per_year=star_data.pm_ra * 1000.0,
            dec_mas_per_year=star_data.pm_dec * 1000.0,
        )
    else:
        # For planets, get the target from ephemeris
        target_name = _PLANET_MAP[body]
        from .planets import _PLANET_FALLBACK

        try:
            target = eph[target_name]
        except KeyError:
            if target_name in _PLANET_FALLBACK:
                target = eph[_PLANET_FALLBACK[target_name]]
            else:
                raise

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Initialize result array with 50 zeros
    dret = [0.0] * 50

    # Calculate positions at the given time
    t = ts.ut1_jd(jd)
    observer_at = earth + observer

    # Calculate Sun position
    sun_app = observer_at.at(t).observe(sun).apparent()
    sun_alt_topo, sun_az, _ = sun_app.altaz()
    sun_alt_deg = sun_alt_topo.degrees
    sun_az_deg = sun_az.degrees

    # Calculate body position (handle both planets and stars)
    if is_star and star_object is not None:
        body_app = observer_at.at(t).observe(star_object).apparent()
    else:
        body_app = observer_at.at(t).observe(target).apparent()
    body_alt_topo, body_az, body_dist = body_app.altaz()
    body_alt_deg = body_alt_topo.degrees
    body_az_deg = body_az.degrees

    # Get geocentric altitude (without refraction or topocentric correction)
    if is_star and star_object is not None:
        body_geo = earth.at(t).observe(star_object).apparent()
    else:
        body_geo = earth.at(t).observe(target).apparent()
    body_geo_ra, body_geo_dec, body_geo_dist = body_geo.radec()

    # Calculate geocentric altitude using hour angle
    # First get the local sidereal time
    gast = t.gast
    lst = gast + lon / 15.0  # Local sidereal time in hours
    ra_hours = body_geo_ra.hours
    ha = (lst - ra_hours) * 15.0  # Hour angle in degrees

    # Geocentric altitude calculation
    dec_rad = math.radians(body_geo_dec.degrees)
    lat_rad = math.radians(lat)
    ha_rad = math.radians(ha)

    sin_alt = math.sin(lat_rad) * math.sin(dec_rad) + math.cos(lat_rad) * math.cos(
        dec_rad
    ) * math.cos(ha_rad)
    sin_alt = max(-1.0, min(1.0, sin_alt))
    geo_alt_deg = math.degrees(math.asin(sin_alt))

    # Calculate atmospheric refraction
    # Use simplified formula: R = 1.02 / tan(h + 10.3/(h + 5.11)) in arcminutes
    if body_alt_deg > -1:
        refraction = 1.02 / math.tan(
            math.radians(body_alt_deg + 10.3 / (body_alt_deg + 5.11))
        )
        refraction /= 60.0  # Convert arcminutes to degrees
    else:
        refraction = 0.5  # Near horizon approximation

    # Apparent altitude (with refraction)
    app_alt_deg = body_alt_deg + refraction

    # Calculate arcus visionis (altitude difference between body and Sun)
    # Topocentric arcus visionis
    tav_act = body_alt_deg - sun_alt_deg

    # Geocentric arcus visionis
    arcv_act = geo_alt_deg - sun_alt_deg

    # Azimuth difference
    daz_act = body_az_deg - sun_az_deg
    # Normalize to -180 to +180
    while daz_act > 180:
        daz_act -= 360
    while daz_act < -180:
        daz_act += 360

    # Get elongation (longitude difference) from Sun
    # For fixed stars, always calculate manually since swe_pheno_ut doesn't support them
    if is_star:
        # Calculate elongation manually for stars
        sun_geo = earth.at(t).observe(sun).apparent()
        elongation = body_geo.separation_from(sun_geo).degrees
        magnitude = star_magnitude
        phase_angle = 0.0  # Stars don't have phase angle
    else:
        try:
            pheno, _ = swe_pheno_ut(jd, body, flags)
            elongation = pheno[2]  # Elongation
            magnitude = pheno[4]  # Visual magnitude
            phase_angle = pheno[1]  # Phase angle
        except Exception:
            # Calculate elongation manually
            sun_geo = earth.at(t).observe(sun).apparent()
            elongation = body_geo.separation_from(sun_geo).degrees
            magnitude = 0.0
            phase_angle = 0.0

    # Calculate extinction coefficient (simplified Schaefer model)
    # k = k_rayleigh + k_aerosol + k_ozone
    # Simplified: k increases with humidity and decreases with altitude
    # Typical values: 0.2 (excellent) to 0.5 (poor)
    pressure_factor = pressure / 1013.25
    humidity_factor = 1.0 + 0.5 * humidity
    altitude_factor = math.exp(-altitude / 8500.0)  # Scale height ~8.5 km
    k_act = 0.25 * pressure_factor * humidity_factor * altitude_factor

    # Minimum topocentric arcus visionis (simplified)
    # Depends on magnitude and atmospheric conditions
    # Brighter objects need less arcus visionis
    if magnitude < 0:
        min_tav = 5.0 + 0.5 * magnitude  # Bright objects
    elif magnitude < 2:
        min_tav = 6.0 + magnitude * 0.5
    else:
        min_tav = 7.0 + magnitude * 0.7

    # Parallax of object (in degrees)
    # Fixed stars have essentially zero parallax
    if is_star:
        parallax = 0.0
    else:
        # Parallax = arcsin(Earth_radius / distance)
        # Earth radius ~ 6371 km, distance in AU (1 AU ~ 149,597,870.7 km)
        earth_radius_au = 6371.0 / 149597870.7
        if body_geo_dist.au > 0:
            parallax = math.degrees(math.asin(earth_radius_au / body_geo_dist.au))
        else:
            parallax = 0.0

    # Calculate rise/set times for object and Sun
    # For morning events, we look for rise times; for evening, set times
    is_morning = event_type in (SE_HELIACAL_RISING, SE_MORNING_LAST)

    # Use a simple estimate for rise/set time based on altitude
    # Rise/set occurs when altitude crosses 0 (corrected for refraction)
    # Time rate: approximately 1 degree of altitude per 4 minutes (at mid-latitudes)
    rise_set_correction = -0.833  # Standard refraction + semidiameter for Sun

    # Estimate object rise/set time
    if body_alt_deg != 0:
        # Rough estimate: time to rise/set based on altitude rate
        # Altitude rate ~ 15 cos(lat) per hour at the horizon
        alt_rate = 15.0 * math.cos(lat_rad)  # degrees per hour
        if alt_rate > 0:
            if is_morning:
                # Object rising: how long until it crosses horizon
                time_to_horizon = (
                    body_alt_deg - rise_set_correction
                ) / alt_rate  # hours
                rise_o = jd - time_to_horizon / 24.0
            else:
                # Object setting
                time_to_horizon = (
                    body_alt_deg - rise_set_correction
                ) / alt_rate  # hours
                rise_o = jd + time_to_horizon / 24.0
        else:
            rise_o = jd
    else:
        rise_o = jd

    # Estimate Sun rise/set time
    if sun_alt_deg != 0:
        alt_rate = 15.0 * math.cos(lat_rad)
        if alt_rate > 0:
            if is_morning:
                time_to_horizon = (sun_alt_deg - rise_set_correction) / alt_rate
                rise_s = jd - time_to_horizon / 24.0
            else:
                time_to_horizon = (sun_alt_deg - rise_set_correction) / alt_rate
                rise_s = jd + time_to_horizon / 24.0
        else:
            rise_s = jd
    else:
        rise_s = jd

    # Lag time (object rise - Sun rise)
    lag = rise_o - rise_s

    # Visibility duration estimate (simplified)
    # Based on how long the object is above horizon while Sun is below
    if sun_alt_deg < -6 and body_alt_deg > 0:
        # Object visible during civil twilight or darker
        # Estimate visibility window
        tvis_vr = abs(sun_alt_deg + 6) / 15.0 / 24.0  # days
    else:
        tvis_vr = 0.0

    # For Moon-specific calculations
    w_moon = 0.0  # Crescent width
    l_moon = 0.0  # Crescent length
    illumination = 0.0

    if body == SE_MOON:
        # Calculate Moon phase and crescent geometry
        try:
            moon_pheno, _ = swe_pheno_ut(jd, SE_MOON, flags)
            phase = moon_pheno[0]  # Phase 0-1
            illumination = phase * 100.0  # Percentage

            # Crescent width approximation (Danjon's formula)
            # W = 15 * (1 - cos(phase_angle/2)) arcminutes (approximate)
            if len(moon_pheno) > 1:
                pa_rad = math.radians(moon_pheno[1])  # Phase angle in radians
                w_moon = 15.0 * (1 - math.cos(pa_rad / 2)) / 60.0  # In degrees

            # Crescent length (semicircle approximation)
            # L = pi * D / 2 where D is diameter
            moon_diameter = 0.5  # Approximately 0.5 degrees
            l_moon = math.pi * moon_diameter / 2
        except Exception:
            pass

    # Yallop q-test (for lunar crescent visibility)
    # q = (ARCV - (11.8371 - 6.3226*W + 0.7319*W^2 - 0.1018*W^3)) / 10
    # Simplified version
    if body == SE_MOON:
        w = w_moon * 60.0  # Convert to arcminutes for formula
        q_criterion = 11.8371 - 6.3226 * w + 0.7319 * w**2 - 0.1018 * w**3
        q_yallop = (arcv_act - q_criterion) / 10.0
        q_crit = q_criterion
    else:
        q_yallop = 0.0
        q_crit = 0.0

    # First, best, and last visibility times (VR = visibility rule)
    # Simplified: based on Sun altitude thresholds
    # First visibility: Sun at about -10
    # Best visibility: Sun at about -8
    # Last visibility: Sun at about -6
    sun_rate = 15.0 * math.cos(lat_rad)  # degrees per hour
    if sun_rate > 0:
        if is_morning:
            t_first_vr = jd + (sun_alt_deg + 10) / sun_rate / 24.0
            t_best_vr = jd + (sun_alt_deg + 8) / sun_rate / 24.0
            t_last_vr = jd + (sun_alt_deg + 6) / sun_rate / 24.0
        else:
            t_first_vr = jd + (-6 - sun_alt_deg) / sun_rate / 24.0
            t_best_vr = jd + (-8 - sun_alt_deg) / sun_rate / 24.0
            t_last_vr = jd + (-10 - sun_alt_deg) / sun_rate / 24.0
    else:
        t_first_vr = jd
        t_best_vr = jd
        t_last_vr = jd

    # Best time according to Yallop (for Moon)
    t_b_yallop = t_best_vr  # Use same as best VR for simplicity

    # Fill in the result array
    dret[0] = body_alt_deg  # AltO - topocentric altitude (unrefracted)
    dret[1] = app_alt_deg  # AppAltO - apparent altitude (refracted)
    dret[2] = geo_alt_deg  # GeoAltO - geocentric altitude
    dret[3] = body_az_deg  # AziO - azimuth of object
    dret[4] = sun_alt_deg  # AltS - topocentric altitude of Sun
    dret[5] = sun_az_deg  # AziS - azimuth of Sun
    dret[6] = tav_act  # TAVact - topocentric arcus visionis
    dret[7] = arcv_act  # ARCVact - geocentric arcus visionis
    dret[8] = daz_act  # DAZact - azimuth difference
    dret[9] = elongation  # ARCLact - elongation from Sun
    dret[10] = k_act  # kact - extinction coefficient
    dret[11] = min_tav  # minTAV - minimum topocentric arcus visionis
    dret[12] = t_first_vr  # TfirstVR - first visibility time
    dret[13] = t_best_vr  # TbVR - best visibility time
    dret[14] = t_last_vr  # TlastVR - last visibility time
    dret[15] = t_b_yallop  # TbYallop - best time according to Yallop
    dret[16] = w_moon  # WMoon - crescent width
    dret[17] = q_yallop  # qYal - Yallop q-test value
    dret[18] = q_crit  # qCrit - Yallop criterion
    dret[19] = parallax  # ParO - parallax of object
    dret[20] = magnitude  # Magn - magnitude
    dret[21] = rise_o  # RiseO - rise/set time of object
    dret[22] = rise_s  # RiseS - rise/set time of Sun
    dret[23] = lag  # Lag - time difference
    dret[24] = tvis_vr  # TvisVR - visibility duration
    dret[25] = l_moon  # LMoon - crescent length
    dret[26] = phase_angle  # CVAact (using phase angle)
    dret[27] = illumination  # Illum - illumination percentage
    # dret[28] onwards are reserved, already 0.0

    return tuple(dret), flags


# Alias for Swiss Ephemeris API compatibility
swe_heliacal_pheno_ut = heliacal_pheno_ut


def vis_limit_mag(
    jd: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    flags: int = SEFLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate the limiting visual magnitude for observing a celestial body.

    This function determines whether a celestial body (planet, star, etc.)
    is visible given the current sky brightness conditions, atmospheric
    parameters, and observer characteristics. It returns both the visibility
    status and detailed information about the observation conditions.

    Args:
        jd: Julian Day (UT) for the observation time
        geopos: Geographic position as a sequence:
            - [0]: Geographic longitude in degrees (east positive)
            - [1]: Geographic latitude in degrees (north positive)
            - [2]: Altitude above sea level in meters
        atmo: Atmospheric conditions as a sequence:
            - [0]: Atmospheric pressure in mbar/hPa
            - [1]: Atmospheric temperature in degrees Celsius
            - [2]: Relative humidity in percent (0-100)
            - [3]: If >= 1: Meteorological Range in km
                   If 0-1: Total atmospheric coefficient (ktot)
                   If 0: Compute ktot from other parameters
        observer: Observer characteristics as a sequence:
            - [0]: Age of observer in years (default 36)
            - [1]: Snellen ratio of observer's eyes (default 1.0 = normal)
            For optical instruments (when HELFLAG_OPTICAL_PARAMS is set):
            - [2]: 0 = monocular, 1 = binocular
            - [3]: Telescope magnification (0 = naked eye)
            - [4]: Optical aperture (telescope diameter) in mm
            - [5]: Optical transmission coefficient
        objname: Name of the object to observe. Can be:
            - Planet name (e.g., "Venus", "Mars", "Jupiter")
            - Fixed star name (e.g., "Sirius", "Aldebaran")
            - Planet number as string (e.g., "2" for Venus)
        flags: Calculation flags combining ephemeris and heliacal flags:
            - SEFLG_SWIEPH, SEFLG_JPLEPH, etc. for ephemeris
            - HELFLAG_OPTICAL_PARAMS: Use optical instrument parameters
            - HELFLAG_NO_DETAILS: Skip detailed calculations
            - HELFLAG_VISLIM_DARK: Assume Sun at nadir (dark sky)
            - HELFLAG_VISLIM_NOMOON: Exclude Moon's brightness contribution

    Returns:
        Tuple containing:
            - result: Visibility status code:
                - (-2): Object is below horizon
                - (0): OK, photopic vision (bright conditions)
                - (1): OK, scotopic vision (dark conditions)
                - (2): OK, near limit between photopic/scotopic
            - dret: Tuple of 8 floats with observation details:
                - [0]: Limiting visual magnitude (object visible if mag < this)
                - [1]: Altitude of object in degrees
                - [2]: Azimuth of object in degrees
                - [3]: Altitude of Sun in degrees
                - [4]: Azimuth of Sun in degrees
                - [5]: Altitude of Moon in degrees
                - [6]: Azimuth of Moon in degrees
                - [7]: Magnitude of object

    Raises:
        ValueError: If objname is empty or invalid

    Algorithm:
        The limiting magnitude calculation is based on Schaefer's model (1990)
        which considers:
        1. Sky brightness from Sun, Moon, zodiacal light, airglow
        2. Atmospheric extinction based on airmass and conditions
        3. Observer's eye adaptation (scotopic vs photopic)
        4. Optional optical instrument characteristics

        The sky background brightness varies with:
        - Sun altitude (twilight contribution)
        - Moon altitude and phase (moonlight)
        - Atmospheric scattering

    Example:
        >>> from libephemeris import julday, vis_limit_mag
        >>> jd = julday(2024, 8, 15, 22.0)
        >>> # Rome location
        >>> geopos = (12.5, 41.9, 0)
        >>> # Standard atmosphere
        >>> atmo = (1013.25, 15.0, 50.0, 0.0)
        >>> # Normal observer
        >>> observer = (36, 1.0)
        >>> result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")
        >>> if dret[0] > dret[7]:
        ...     print("Venus is visible")
        ... else:
        ...     print("Venus is not visible")

    References:
        - Swiss Ephemeris: swe_vis_limit_mag()
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
    """
    from .planets import _PLANET_MAP, swe_pheno_ut
    from .fixed_stars import swe_fixstar2_ut, swe_fixstar2_mag
    from .state import get_planets, get_timescale
    from .constants import (
        SE_SUN,
        SE_MOON,
        SE_HELFLAG_VISLIM_DARK,
        SE_HELFLAG_VISLIM_NOMOON,
        SE_HELFLAG_BELOW_HORIZON,
        SE_HELFLAG_PHOTOPIC,
        SE_HELFLAG_SCOTOPIC,
        SE_HELFLAG_MIXED,
    )
    from skyfield.api import wgs84

    if not objname:
        raise ValueError("objname cannot be empty")

    # Parse geographic position
    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    alt_m = geopos[2] if len(geopos) > 2 else 0.0

    # Parse atmospheric conditions
    pressure = atmo[0] if len(atmo) > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 50.0
    met_range = atmo[3] if len(atmo) > 3 else 0.0

    # Parse observer data
    observer_age = observer[0] if len(observer) > 0 else 36.0
    snellen_ratio = observer[1] if len(observer) > 1 else 1.0

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Create observer location
    obs_location = wgs84.latlon(lat, lon, alt_m)
    earth = eph["earth"]
    sun = eph["sun"]
    moon = eph["moon"]

    t = ts.ut1_jd(jd)
    observer_at = earth + obs_location

    # Calculate Sun position
    sun_app = observer_at.at(t).observe(sun).apparent()
    sun_alt_deg, sun_az_deg, _ = sun_app.altaz()
    sun_alt = sun_alt_deg.degrees
    sun_az = sun_az_deg.degrees

    # Calculate Moon position
    moon_app = observer_at.at(t).observe(moon).apparent()
    moon_alt_deg, moon_az_deg, _ = moon_app.altaz()
    moon_alt = moon_alt_deg.degrees
    moon_az = moon_az_deg.degrees

    # Determine if objname is a planet ID or name
    body_id = None
    is_fixed_star = False

    # Try parsing as integer (planet ID)
    try:
        body_id = int(objname)
    except ValueError:
        # Try to find planet by name
        name_upper = objname.upper().strip()
        planet_names = {
            "SUN": SE_SUN,
            "MOON": SE_MOON,
            "MERCURY": 2,
            "VENUS": 3,
            "MARS": 4,
            "JUPITER": 5,
            "SATURN": 6,
            "URANUS": 7,
            "NEPTUNE": 8,
            "PLUTO": 9,
        }
        if name_upper in planet_names:
            body_id = planet_names[name_upper]
        else:
            # Assume it's a fixed star
            is_fixed_star = True

    # Calculate object position and magnitude
    obj_alt = 0.0
    obj_az = 0.0
    obj_mag = 0.0

    if is_fixed_star:
        # Fixed star calculation
        try:
            star_name_out, star_result, retflag, error = swe_fixstar2_ut(
                objname, jd, flags & 0xFF
            )
            if error:
                raise ValueError(f"could not find star name {objname.lower()}: {error}")

            # star_result is (lon, lat, dist, lon_speed, lat_speed, dist_speed)
            # We need to convert ecliptic to horizontal
            star_lon = star_result[0]
            star_lat = star_result[1]

            # Get star magnitude
            star_name_mag, star_mag_val, mag_error = swe_fixstar2_mag(objname)
            if not mag_error:
                obj_mag = star_mag_val
            else:
                obj_mag = 2.0  # Default magnitude if not found

            # Convert ecliptic to equatorial then to horizontal
            # Simplified: use azalt function if available
            from .utils import azalt, SE_ECL2HOR

            hor_result = azalt(
                jd,
                SE_ECL2HOR,
                lat,
                lon,
                alt_m,
                pressure,
                temperature,
                (star_lon, star_lat, 1.0),
            )
            obj_az = hor_result[0]
            obj_alt = hor_result[1]

        except ValueError:
            raise
        except Exception as e:
            # Star not found or other error
            raise ValueError(f"could not find star name {objname.lower()}: {e}")
    else:
        # Planet calculation
        if body_id is None:
            raise ValueError(f"Unknown object: {objname}")

        # Get planet name from _PLANET_MAP
        if body_id in _PLANET_MAP:
            target_name = _PLANET_MAP[body_id]
            # Try planet center first, fall back to barycenter if not available
            from .planets import _PLANET_FALLBACK

            try:
                target = eph[target_name]
            except KeyError:
                if target_name in _PLANET_FALLBACK:
                    target = eph[_PLANET_FALLBACK[target_name]]
                else:
                    raise

            # Calculate position
            body_app = observer_at.at(t).observe(target).apparent()
            body_alt_deg, body_az_deg, _ = body_app.altaz()
            obj_alt = body_alt_deg.degrees
            obj_az = body_az_deg.degrees

            # Get magnitude from pheno
            try:
                pheno_result, _ = swe_pheno_ut(jd, body_id, flags)
                obj_mag = pheno_result[4]  # Visual magnitude
            except Exception:
                obj_mag = 0.0  # Default bright
        else:
            raise ValueError(f"illegal planet number {body_id}.")

    # Check if object is below horizon
    if obj_alt < 0:
        # Apply HELFLAG options before returning
        use_dark_sky_early = bool(flags & SE_HELFLAG_VISLIM_DARK)
        exclude_moon_early = bool(flags & SE_HELFLAG_VISLIM_NOMOON)
        ret_sun_alt = -90.0 if use_dark_sky_early else sun_alt
        ret_moon_alt = -90.0 if exclude_moon_early else moon_alt
        dret = (
            0.0,
            obj_alt,
            obj_az,
            ret_sun_alt,
            sun_az,
            ret_moon_alt,
            moon_az,
            obj_mag,
        )
        return SE_HELFLAG_BELOW_HORIZON, dret

    # Apply HELFLAG options
    use_dark_sky = bool(flags & SE_HELFLAG_VISLIM_DARK)
    exclude_moon = bool(flags & SE_HELFLAG_VISLIM_NOMOON)

    if use_dark_sky:
        sun_alt = -90.0  # Assume Sun at nadir

    if exclude_moon:
        moon_alt = -90.0  # Assume Moon at nadir

    # Calculate atmospheric extinction coefficient
    # Simplified Schaefer model
    def calc_extinction_coeff(pressure_mbar, temp_c, humidity_pct, met_range_km):
        """Calculate total atmospheric extinction coefficient."""
        # Rayleigh scattering (molecular)
        k_rayleigh = 0.1451 * (pressure_mbar / 1013.25)

        # Aerosol scattering
        if met_range_km >= 1:
            # Meteorological range given
            k_aerosol = 3.912 / met_range_km - 0.106
        elif 0 < met_range_km < 1:
            # ktot given directly
            return met_range_km
        else:
            # Estimate from humidity
            # Higher humidity = more scattering
            k_aerosol = 0.1 + 0.2 * (humidity_pct / 100.0)

        # Ozone absorption (small contribution)
        k_ozone = 0.016

        return k_rayleigh + k_aerosol + k_ozone

    k_ext = calc_extinction_coeff(pressure, temperature, humidity_pct, met_range)

    # Calculate airmass
    def calc_airmass(altitude_deg):
        """Calculate airmass using Kasten & Young formula."""
        if altitude_deg <= 0:
            return 40.0  # Maximum airmass at/below horizon
        z = 90.0 - altitude_deg  # Zenith angle
        z_rad = math.radians(z)
        # Kasten & Young (1989)
        airmass = 1.0 / (math.cos(z_rad) + 0.50572 * (96.07995 - z) ** (-1.6364))
        return min(airmass, 40.0)

    airmass = calc_airmass(obj_alt)

    # Calculate sky brightness contribution from Sun
    def calc_sun_brightness_contribution(sun_altitude):
        """Calculate sky brightness factor from Sun position."""
        if sun_altitude >= 0:
            return 1e10  # Daylight - essentially infinite brightness
        elif sun_altitude >= -6:
            # Civil twilight
            return 10 ** (4.0 + sun_altitude * 0.4)
        elif sun_altitude >= -12:
            # Nautical twilight
            return 10 ** (1.6 + (sun_altitude + 6) * 0.2)
        elif sun_altitude >= -18:
            # Astronomical twilight
            return 10 ** (0.4 + (sun_altitude + 12) * 0.1)
        else:
            # Full night
            return 1.0

    # Calculate sky brightness contribution from Moon
    def calc_moon_brightness_contribution(moon_altitude, jd_ut):
        """Calculate sky brightness factor from Moon position and phase."""
        if moon_altitude <= 0:
            return 0.0

        # Get Moon phase
        try:
            moon_pheno, _ = swe_pheno_ut(jd_ut, SE_MOON, flags & 0xFF)
            phase_angle = moon_pheno[0]  # Phase angle in degrees
            illumination = (1 + math.cos(math.radians(phase_angle))) / 2.0
        except Exception:
            illumination = 0.5  # Assume half moon

        # Moon brightness increases with altitude and illumination
        moon_factor = illumination * math.sin(math.radians(moon_altitude))
        return moon_factor * 100.0

    sun_contrib = calc_sun_brightness_contribution(sun_alt)
    moon_contrib = (
        calc_moon_brightness_contribution(moon_alt, jd) if not exclude_moon else 0.0
    )

    # Total sky brightness (arbitrary units)
    total_sky_brightness = sun_contrib + moon_contrib + 1.0  # 1.0 = natural sky glow

    # Calculate limiting magnitude
    # Base limiting magnitude for perfect conditions: ~6.5 for naked eye
    base_limit_mag = 6.5

    # Adjust for observer characteristics
    # Age reduces sensitivity (approximately 0.05 mag per decade over 30)
    age_factor = max(0, (observer_age - 30) / 10.0) * 0.05

    # Snellen ratio adjustment (better eyes = higher limit)
    snellen_factor = -2.5 * math.log10(snellen_ratio) if snellen_ratio > 0 else 0

    # Airmass extinction
    airmass_extinction = k_ext * airmass

    # Sky brightness reduction of limiting magnitude
    if total_sky_brightness > 1:
        sky_reduction = 2.5 * math.log10(total_sky_brightness)
    else:
        sky_reduction = 0.0

    # Calculate final limiting magnitude
    limiting_mag = base_limit_mag - age_factor - snellen_factor - sky_reduction

    # Apply extinction to object magnitude (object appears fainter)
    apparent_obj_mag = obj_mag + airmass_extinction

    # Determine vision type based on sky brightness
    if sun_alt >= -6:
        # Photopic vision (daylight/bright twilight)
        vision_type = SE_HELFLAG_PHOTOPIC
    elif sun_alt >= -12:
        # Mixed/mesopic vision
        vision_type = SE_HELFLAG_MIXED
    else:
        # Scotopic vision (night)
        vision_type = SE_HELFLAG_SCOTOPIC

    # Build result tuple
    dret = (
        limiting_mag,  # 0: Limiting visual magnitude
        obj_alt,  # 1: Altitude of object
        obj_az,  # 2: Azimuth of object
        sun_alt,  # 3: Altitude of Sun
        sun_az,  # 4: Azimuth of Sun
        moon_alt,  # 5: Altitude of Moon
        moon_az,  # 6: Azimuth of Moon
        apparent_obj_mag,  # 7: Apparent magnitude of object (with extinction)
    )

    return vision_type, dret


# Alias for Swiss Ephemeris API compatibility
swe_vis_limit_mag = vis_limit_mag
