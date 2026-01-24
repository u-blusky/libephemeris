"""
Utility functions for libephemeris.

Provides helper functions compatible with pyswisseph API including
angular calculations and other mathematical utilities.
"""

import math
from typing import Optional, Tuple

# Azalt calculation method flags (compatible with pyswisseph)
SE_ECL2HOR: int = 0  # Ecliptic coordinates to horizontal
SE_EQU2HOR: int = 1  # Equatorial coordinates to horizontal

# Azalt_rev calculation method flags (compatible with pyswisseph)
SE_HOR2ECL: int = 0  # Horizontal to ecliptic coordinates
SE_HOR2EQU: int = 1  # Horizontal to equatorial coordinates

# Refraction calculation flags (compatible with pyswisseph)
SE_TRUE_TO_APP: int = 0  # True altitude to apparent altitude
SE_APP_TO_TRUE: int = 1  # Apparent altitude to true altitude


def cotrans_sp(
    coord: Tuple[float, float, float],
    speed: Tuple[float, float, float],
    obliquity: float,
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """
    Transform coordinates and velocities between ecliptic and equatorial systems.

    This function extends cotrans() to also transform velocity (speed) components,
    using analytical derivatives of the coordinate transformation equations.

    The direction of transformation depends on the sign of obliquity:
    - Positive obliquity: ecliptic (lon, lat) -> equatorial (RA, Dec)
    - Negative obliquity: equatorial (RA, Dec) -> ecliptic (lon, lat)

    Args:
        coord: Tuple of (longitude/RA, latitude/Dec, distance) in degrees
        speed: Tuple of (lon_speed, lat_speed, dist_speed) in degrees/day
        obliquity: Obliquity of the ecliptic in degrees.
                   Positive for ecliptic->equatorial, negative for equatorial->ecliptic.

    Returns:
        Tuple of (coord_transformed, speed_transformed) where:
        - coord_transformed: (new_lon/RA, new_lat/Dec, distance)
        - speed_transformed: (new_lon_speed, new_lat_speed, dist_speed)
        Distance and distance speed are unchanged by the transformation.

    Examples:
        >>> # Ecliptic to equatorial (positive obliquity)
        >>> coord, speed = cotrans_sp((90.0, 0.0, 1.0), (1.0, 0.0, 0.0), 23.4)
        >>> # Returns transformed position and velocity
    """
    lon = coord[0]
    lat = coord[1]
    dist = coord[2]
    lon_speed = speed[0]
    lat_speed = speed[1]
    dist_speed = speed[2]

    # Convert to radians
    lon_rad = math.radians(lon)
    lat_rad = math.radians(lat)
    eps_rad = math.radians(obliquity)

    # Precompute trig values
    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # Calculate the new latitude (Dec for ecl->eq, beta for eq->ecl)
    # sin(new_lat) = sin(lat) * cos(eps) - cos(lat) * sin(eps) * sin(lon)
    sin_new_lat = sin_lat * cos_eps - cos_lat * sin_eps * sin_lon
    # Clamp to avoid numerical issues with asin
    sin_new_lat = max(-1.0, min(1.0, sin_new_lat))
    new_lat_rad = math.asin(sin_new_lat)
    cos_new_lat = math.cos(new_lat_rad)

    # Calculate the new longitude (RA for ecl->eq, lambda for eq->ecl)
    # tan(new_lon) = (sin(lon) * cos(eps) + tan(lat) * sin(eps)) / cos(lon)
    tan_lat = math.tan(lat_rad)
    y = sin_lon * cos_eps + tan_lat * sin_eps
    x = cos_lon

    new_lon_rad = math.atan2(y, x)

    # Convert back to degrees
    new_lon = math.degrees(new_lon_rad)
    new_lat = math.degrees(new_lat_rad)

    # Normalize longitude to [0, 360)
    new_lon = new_lon % 360.0

    # --- Velocity transformation ---
    # Convert speed to radians/day for the calculation
    lon_speed_rad = math.radians(lon_speed)
    lat_speed_rad = math.radians(lat_speed)

    # Derivative of new latitude:
    # d/dt[sin(new_lat)] = cos(new_lat) * d(new_lat)/dt
    # d/dt[sin(lat)*cos(eps) - cos(lat)*sin(eps)*sin(lon)]
    #   = cos(lat)*cos(eps)*d(lat)/dt + sin(lat)*sin(eps)*sin(lon)*d(lat)/dt
    #     - cos(lat)*sin(eps)*cos(lon)*d(lon)/dt
    #
    # new_lat_speed = (1/cos(new_lat)) * [
    #   (cos(lat)*cos(eps) + sin(lat)*sin(eps)*sin(lon)) * lat_speed
    #   - cos(lat)*sin(eps)*cos(lon) * lon_speed
    # ]
    if abs(cos_new_lat) > 1e-10:
        new_lat_speed_rad = (
            (cos_lat * cos_eps + sin_lat * sin_eps * sin_lon) * lat_speed_rad
            - cos_lat * sin_eps * cos_lon * lon_speed_rad
        ) / cos_new_lat
    else:
        # At poles, latitude speed is undefined; use 0
        new_lat_speed_rad = 0.0

    # Derivative of new longitude:
    # new_lon = atan2(y, x) where y = sin(lon)*cos(eps) + tan(lat)*sin(eps), x = cos(lon)
    # d(new_lon)/dt = (x * dy/dt - y * dx/dt) / (x^2 + y^2)
    #
    # dx/dt = -sin(lon) * d(lon)/dt
    # dy/dt = cos(lon)*cos(eps)*d(lon)/dt + sec^2(lat)*sin(eps)*d(lat)/dt
    #       = cos(lon)*cos(eps)*d(lon)/dt + sin(eps)/(cos^2(lat))*d(lat)/dt
    dx_dt = -sin_lon * lon_speed_rad
    cos_lat_sq = cos_lat * cos_lat
    if abs(cos_lat_sq) > 1e-10:
        dy_dt = (
            cos_lon * cos_eps * lon_speed_rad + (sin_eps / cos_lat_sq) * lat_speed_rad
        )
    else:
        # At poles of input coordinates
        dy_dt = cos_lon * cos_eps * lon_speed_rad

    denom = x * x + y * y
    if abs(denom) > 1e-10:
        new_lon_speed_rad = (x * dy_dt - y * dx_dt) / denom
    else:
        new_lon_speed_rad = 0.0

    # Convert speeds back to degrees/day
    new_lon_speed = math.degrees(new_lon_speed_rad)
    new_lat_speed = math.degrees(new_lat_speed_rad)

    return (
        (new_lon, new_lat, dist),
        (new_lon_speed, new_lat_speed, dist_speed),
    )


def azalt(
    jd: float,
    calc_flag: int,
    lat: float,
    lon: float,
    altitude: float,
    pressure: float,
    temperature: float,
    coord: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """
    Convert equatorial or ecliptic coordinates to horizontal (azimuth/altitude).

    This function transforms celestial coordinates to horizontal coordinates
    (azimuth and altitude) for a given observer location and time.
    It accounts for atmospheric refraction.

    Compatible with pyswisseph's swe.azalt() function.

    Args:
        jd: Julian Day in Universal Time (UT1)
        calc_flag: Coordinate type flag:
            - SE_ECL2HOR (0): Input is ecliptic (longitude, latitude, distance)
            - SE_EQU2HOR (1): Input is equatorial (RA, Dec, distance)
        lat: Geographic latitude of observer in degrees (North positive)
        lon: Geographic longitude of observer in degrees (East positive)
        altitude: Observer altitude above sea level in meters
        pressure: Atmospheric pressure in mbar (hPa). Use 0 for no refraction.
        temperature: Atmospheric temperature in Celsius
        coord: Tuple of (longitude/RA, latitude/Dec, distance) in degrees

    Returns:
        Tuple of (azimuth, true_altitude, apparent_altitude) where:
            - azimuth: Degrees from South, westward (0=South, 90=West, 180=North, 270=East)
            - true_altitude: Geometric altitude without refraction (degrees)
            - apparent_altitude: Altitude with atmospheric refraction applied (degrees)

    Note:
        - Swiss Ephemeris convention: Azimuth is measured from South, westward.
          This differs from the common convention (from North, eastward).
          To convert: az_from_north = (azimuth + 180) % 360
        - Refraction is calculated using the Bennett formula, which is accurate
          to about 0.07 arcmin for altitudes above 15°.
        - When pressure=0, no refraction is applied (apparent_alt = true_alt)
        - For objects below the horizon (negative altitude), refraction is
          extrapolated but becomes less accurate.

    Example:
        >>> from libephemeris import azalt, SE_EQU2HOR, julday
        >>> jd = julday(2024, 6, 15, 12.0)
        >>> # RA=90°, Dec=23.5°, observer at lat=41.9°N, lon=12.5°E
        >>> az, alt_true, alt_app = azalt(jd, SE_EQU2HOR, 41.9, 12.5, 0, 1013.25, 15, (90, 23.5, 1))
    """
    from .time_utils import sidtime
    from skyfield.nutationlib import iau2000b_radians
    from .state import get_timescale

    # Get true obliquity (mean obliquity + nutation in obliquity)
    # This is needed for accurate ecliptic to equatorial conversion
    T = (jd - 2451545.0) / 36525.0

    # Mean obliquity (Laskar 1986 formula)
    eps0 = (
        84381.406 - 46.836769 * T - 0.0001831 * T * T + 0.00200340 * T * T * T
    ) / 3600.0  # Convert arcsec to degrees

    # Add nutation in obliquity for true obliquity
    ts = get_timescale()
    t = ts.ut1_jd(jd)
    dpsi_rad, deps_rad = iau2000b_radians(t)
    deps_deg = math.degrees(deps_rad)
    eps = eps0 + deps_deg  # True obliquity

    # Convert input coordinates to equatorial (RA, Dec) if ecliptic
    if calc_flag == SE_ECL2HOR:
        # Input is ecliptic: convert to equatorial
        ecl_lon = coord[0]
        ecl_lat = coord[1]
        dist = coord[2]

        # Convert ecliptic to equatorial
        # Note: cotrans uses convention where NEGATIVE obliquity converts
        # ecliptic to equatorial with correct sign (matching astronomical convention)
        eq_coord = cotrans((ecl_lon, ecl_lat, dist), -eps)
        ra = eq_coord[0]
        dec = eq_coord[1]
    else:
        # Input is already equatorial (RA, Dec)
        ra = coord[0]
        dec = coord[1]

    # Calculate Local Sidereal Time
    # We need nutation for proper sidereal time, but for azalt we can use mean
    nutation = 0.0  # Mean sidereal time approximation
    lst_hours = sidtime(jd, lon, eps, nutation)
    lst_deg = lst_hours * 15.0  # Convert hours to degrees

    # Calculate Hour Angle
    # H = LST - RA
    ha = (lst_deg - ra) % 360.0

    # Convert to radians for calculation
    ha_rad = math.radians(ha)
    dec_rad = math.radians(dec)
    lat_rad = math.radians(lat)

    # Calculate altitude (geometric/true)
    # sin(alt) = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(H)
    sin_alt = math.sin(lat_rad) * math.sin(dec_rad) + math.cos(lat_rad) * math.cos(
        dec_rad
    ) * math.cos(ha_rad)
    # Clamp to avoid numerical issues with asin
    sin_alt = max(-1.0, min(1.0, sin_alt))
    alt_true = math.degrees(math.asin(sin_alt))

    # Calculate azimuth
    # tan(A) = sin(H) / (cos(H) * sin(lat) - tan(dec) * cos(lat))
    # Using atan2 for proper quadrant handling
    sin_ha = math.sin(ha_rad)
    cos_ha = math.cos(ha_rad)
    cos_dec = math.cos(dec_rad)

    # Numerator and denominator for azimuth calculation
    # Swiss Ephemeris convention: azimuth from South, westward
    y = sin_ha * cos_dec
    x = cos_ha * cos_dec * math.sin(lat_rad) - math.sin(dec_rad) * math.cos(lat_rad)

    az_rad = math.atan2(y, x)
    azimuth = math.degrees(az_rad)

    # Normalize azimuth to 0-360 range
    azimuth = azimuth % 360.0
    if azimuth < 0:
        azimuth += 360.0

    # Calculate atmospheric refraction
    if pressure > 0 and alt_true > -2.0:
        # Bennett formula for atmospheric refraction
        # R = 1.02 / tan(h + 10.3/(h + 5.11)) [arcminutes]
        # Where h is apparent altitude in degrees
        # This is an iterative formula since R depends on apparent altitude

        # For true altitude, we use a modified approach:
        # Start with true altitude and calculate refraction
        # Refraction is applied as: apparent_alt = true_alt + R

        # Pressure and temperature correction factors
        # Standard conditions: 1010 mbar, 10°C
        pressure_factor = pressure / 1010.0
        temperature_factor = 283.0 / (273.0 + temperature)
        correction = pressure_factor * temperature_factor

        if alt_true > 15.0:
            # Simple formula for high altitudes
            # R = 58.1" * tan(z) - 0.07" * tan^3(z)  where z = zenith angle
            z_rad = math.radians(90.0 - alt_true)
            tan_z = math.tan(z_rad)
            refraction_arcsec = 58.1 * tan_z - 0.07 * tan_z**3
            refraction = refraction_arcsec / 3600.0 * correction
        elif alt_true > -1.0:
            # Bennett formula for lower altitudes
            # More accurate near the horizon
            h = alt_true
            if h < 0.01:
                h = 0.01  # Avoid division issues
            # Refraction in arcminutes
            r_arcmin = 1.02 / math.tan(math.radians(h + 10.3 / (h + 5.11)))
            refraction = r_arcmin / 60.0 * correction
        else:
            # Below horizon: extrapolate carefully
            # Use a simple linear extrapolation from horizon value
            h = -1.0
            r_arcmin = 1.02 / math.tan(math.radians(h + 10.3 / (h + 5.11)))
            refraction_at_horizon = r_arcmin / 60.0 * correction
            # Linear decrease below horizon
            refraction = refraction_at_horizon + (alt_true + 1.0) * 0.1

        alt_apparent = alt_true + refraction
    else:
        # No refraction correction
        alt_apparent = alt_true

    return (azimuth, alt_true, alt_apparent)


def azalt_rev(
    jd: float,
    calc_flag: int,
    lat: float,
    lon: float,
    altitude: float,
    azimuth: float,
    true_altitude: float,
) -> Tuple[float, float]:
    """
    Convert horizontal coordinates (azimuth/altitude) to equatorial or ecliptic.

    This function is the inverse of azalt(): it transforms horizontal coordinates
    (azimuth and true altitude) to celestial coordinates (equatorial or ecliptic)
    for a given observer location and time.

    Compatible with pyswisseph's swe.azalt_rev() function.

    Note: This function is not precisely the reverse of azalt(). If only an
    apparent altitude is available, the true altitude must first be computed
    using a refraction correction.

    Args:
        jd: Julian Day in Universal Time (UT1)
        calc_flag: Output coordinate type flag:
            - SE_HOR2EQU (1): Output is equatorial (RA, Dec)
            - SE_HOR2ECL (0): Output is ecliptic (longitude, latitude)
        lat: Geographic latitude of observer in degrees (North positive)
        lon: Geographic longitude of observer in degrees (East positive)
        altitude: Observer altitude above sea level in meters
        azimuth: Azimuth in degrees from South, westward
                 (0=South, 90=West, 180=North, 270=East)
        true_altitude: True (geometric) altitude above horizon in degrees

    Returns:
        Tuple of (x1, x2) where:
            - If calc_flag == SE_HOR2EQU: (Right Ascension, Declination) in degrees
            - If calc_flag == SE_HOR2ECL: (Ecliptic longitude, Ecliptic latitude) in degrees

    Example:
        >>> from libephemeris import azalt_rev, SE_HOR2EQU, julday
        >>> jd = julday(2024, 6, 15, 12.0)
        >>> # Object at azimuth 90° (West), altitude 45°, observer at Rome
        >>> ra, dec = azalt_rev(jd, SE_HOR2EQU, 41.9, 12.5, 0, 90.0, 45.0)
    """
    from .time_utils import sidtime
    from skyfield.nutationlib import iau2000b_radians
    from .state import get_timescale

    # Get true obliquity (mean obliquity + nutation in obliquity)
    # This is needed for accurate equatorial to ecliptic conversion
    T = (jd - 2451545.0) / 36525.0

    # Mean obliquity (Laskar 1986 formula)
    eps0 = (
        84381.406 - 46.836769 * T - 0.0001831 * T * T + 0.00200340 * T * T * T
    ) / 3600.0  # Convert arcsec to degrees

    # Add nutation in obliquity for true obliquity
    ts = get_timescale()
    t = ts.ut1_jd(jd)
    dpsi_rad, deps_rad = iau2000b_radians(t)
    deps_deg = math.degrees(deps_rad)
    eps = eps0 + deps_deg  # True obliquity

    # Convert to radians for calculation
    az_rad = math.radians(azimuth)
    alt_rad = math.radians(true_altitude)
    lat_rad = math.radians(lat)

    # Calculate declination from horizontal coordinates
    # For azimuth measured from South (westward), the formula is:
    # sin(dec) = sin(lat) * sin(alt) - cos(lat) * cos(alt) * cos(az)
    # Note: The minus sign is because azimuth is from South, not North
    sin_dec = math.sin(lat_rad) * math.sin(alt_rad) - math.cos(lat_rad) * math.cos(
        alt_rad
    ) * math.cos(az_rad)
    # Clamp to avoid numerical issues with asin
    sin_dec = max(-1.0, min(1.0, sin_dec))
    dec_rad = math.asin(sin_dec)
    dec = math.degrees(dec_rad)

    # Calculate hour angle from horizontal coordinates
    # From the forward azalt transform, we have:
    # y = sin(H) * cos(dec)  -> used in atan2 to get azimuth
    # x = cos(H) * cos(dec) * sin(lat) - sin(dec) * cos(lat)
    #
    # Therefore:
    # sin(H) = sin(az) * cos(alt) / cos(dec)
    # cos(H) = (cos(az) * cos(alt) + sin(dec) * cos(lat)) / (cos(dec) * sin(lat))
    cos_dec = math.cos(dec_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)

    if abs(cos_dec) < 1e-10 or abs(sin_lat) < 1e-10:
        # At poles or object at celestial pole, hour angle is undefined
        # Use 0 as default
        ha = 0.0
    else:
        sin_H = math.sin(az_rad) * math.cos(alt_rad) / cos_dec
        cos_H = (math.cos(az_rad) * math.cos(alt_rad) + sin_dec * cos_lat) / (
            cos_dec * sin_lat
        )

        ha_rad = math.atan2(sin_H, cos_H)
        ha = math.degrees(ha_rad)

    # Calculate Local Sidereal Time
    nutation = 0.0  # Mean sidereal time approximation (same as azalt)
    lst_hours = sidtime(jd, lon, eps, nutation)
    lst_deg = lst_hours * 15.0  # Convert hours to degrees

    # Calculate Right Ascension: RA = LST - H
    ra = (lst_deg - ha) % 360.0

    if calc_flag == SE_HOR2EQU:
        # Return equatorial coordinates (RA, Dec)
        return (ra, dec)
    else:
        # SE_HOR2ECL: Convert equatorial to ecliptic
        # cotrans with negative obliquity converts equatorial to ecliptic
        ecl_coord = cotrans((ra, dec, 1.0), eps)
        return (ecl_coord[0], ecl_coord[1])


def refrac(
    altitude: float,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    calc_flag: int = SE_TRUE_TO_APP,
) -> float:
    """
    Calculate true altitude from apparent altitude, or vice-versa.

    Atmospheric refraction makes celestial objects appear higher than their
    true (geometric) position. The effect is strongest near the horizon
    (about 34 arcminutes at 0 degrees) and negligible at high altitudes.

    This function converts between true (geometric) altitude and apparent
    (observed) altitude by adding or removing the refraction correction.

    Compatible with pyswisseph's swe.refrac() function.

    Args:
        altitude: Altitude in degrees. For SE_TRUE_TO_APP, this is the true
                  (geometric) altitude. For SE_APP_TO_TRUE, this is the apparent
                  (observed) altitude.
        pressure: Atmospheric pressure in mbar (hPa). Default is 1013.25 (sea level).
                  Use 0 to disable refraction correction (returns input altitude).
        temperature: Atmospheric temperature in Celsius. Default is 15.0.
        calc_flag: Direction of conversion:
            - SE_TRUE_TO_APP (0): Convert true altitude to apparent altitude
              (add refraction - object appears higher)
            - SE_APP_TO_TRUE (1): Convert apparent altitude to true altitude
              (subtract refraction - object's true position)

    Returns:
        The converted altitude in degrees:
        - For SE_TRUE_TO_APP: apparent altitude (= true altitude + refraction)
        - For SE_APP_TO_TRUE: true altitude (= apparent altitude - refraction)

    Notes:
        - At the horizon (0 degrees), refraction is approximately 34 arcminutes
          (0.567 degrees) under standard atmospheric conditions.
        - The formula used is the Saemund-Bennett formula, which is accurate
          to about 0.07 arcminutes for altitudes above 15 degrees.
        - For altitudes below -2 degrees, refraction is extrapolated linearly.
        - Pressure and temperature adjustments follow the standard formula:
          correction = (pressure / 1010) * (283 / (273 + temperature))

    Examples:
        >>> # True altitude at horizon -> apparent altitude is higher
        >>> refrac(0.0, 1013.25, 15.0, SE_TRUE_TO_APP)
        0.476...  # apparent altitude
        >>> # Apparent altitude at horizon -> true altitude (0° apparent = ~-0.5° true)
        >>> refrac(0.5, 1013.25, 15.0, SE_APP_TO_TRUE)
        0.0...  # approximately 0° true altitude
        >>> # No refraction when pressure is 0
        >>> refrac(10.0, 0, 15.0, SE_TRUE_TO_APP)
        10.0  # returns input altitude unchanged
    """
    # No refraction if pressure is zero or negative
    if pressure <= 0:
        return altitude

    # Pressure and temperature correction factors
    # Standard conditions: 1010 mbar, 10°C (283 K)
    pressure_factor = pressure / 1010.0
    temperature_factor = 283.0 / (273.0 + temperature)
    correction = pressure_factor * temperature_factor

    if calc_flag == SE_TRUE_TO_APP:
        # True altitude to apparent altitude (add refraction)
        alt = altitude

        if alt > 15.0:
            # Simple formula for high altitudes
            # R = 58.1" * tan(z) - 0.07" * tan^3(z) where z = zenith angle
            z_rad = math.radians(90.0 - alt)
            tan_z = math.tan(z_rad)
            refraction_arcsec = 58.1 * tan_z - 0.07 * tan_z**3
            refraction = refraction_arcsec / 3600.0 * correction
        elif alt > -2.0:
            # Bennett formula for lower altitudes (more accurate near horizon)
            # R = 1.02 / tan(h + 10.3/(h + 5.11)) arcminutes
            h = max(alt, 0.01)  # Avoid division issues
            r_arcmin = 1.02 / math.tan(math.radians(h + 10.3 / (h + 5.11)))
            refraction = r_arcmin / 60.0 * correction
        else:
            # Below -2 degrees: extrapolate linearly
            # Calculate refraction at -2 degrees and extrapolate
            h = -2.0
            r_arcmin = 1.02 / math.tan(math.radians(h + 10.3 / (h + 5.11)))
            refraction_at_minus2 = r_arcmin / 60.0 * correction
            # Linear decrease for more negative altitudes
            refraction = refraction_at_minus2 + (alt + 2.0) * 0.1

        return altitude + refraction

    else:
        # SE_APP_TO_TRUE: Apparent altitude to true altitude (subtract refraction)
        # We need to find what true altitude would give this apparent altitude

        # First, calculate the refraction at true altitude 0° (horizon)
        # This is the threshold below which APP_TO_TRUE returns input unchanged
        h = 0.01  # Small value near 0
        r_arcmin = 1.02 / math.tan(math.radians(h + 10.3 / (h + 5.11)))
        horizon_refraction = r_arcmin / 60.0 * correction

        # If apparent altitude is at or below the horizon refraction threshold,
        # return input unchanged (matching pyswisseph behavior)
        if altitude <= horizon_refraction:
            return altitude

        alt = altitude

        if alt > 15.0:
            # Simple formula for high altitudes
            z_rad = math.radians(90.0 - alt)
            tan_z = math.tan(z_rad)
            refraction_arcsec = 58.1 * tan_z - 0.07 * tan_z**3
            refraction = refraction_arcsec / 3600.0 * correction
        else:
            # Bennett formula (using apparent altitude directly is a good approximation)
            h = max(alt, 0.01)  # Avoid division issues
            r_arcmin = 1.02 / math.tan(math.radians(h + 10.3 / (h + 5.11)))
            refraction = r_arcmin / 60.0 * correction

        return altitude - refraction


def refrac_extended(
    altitude: float,
    altitude_geo: float,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    lapse_rate: Optional[float] = None,
    calc_flag: int = SE_TRUE_TO_APP,
) -> Tuple[float, Tuple[float, float, float, float]]:
    """
    Calculate true altitude from apparent altitude, or vice-versa (extended).

    This is an extended version of refrac() that includes:
    - Observer altitude above sea level
    - Atmospheric lapse rate (temperature variation with altitude)
    - Dip of the horizon calculation

    Compatible with pyswisseph's swe.refrac_extended() function.

    Args:
        altitude: Altitude of object above geometric horizon in degrees.
                  For SE_TRUE_TO_APP, this is the true (geometric) altitude.
                  For SE_APP_TO_TRUE, this is the apparent (observed) altitude.
        altitude_geo: Altitude of observer above sea level in meters.
        pressure: Atmospheric pressure in mbar (hPa). Default is 1013.25 (sea level).
                  Use 0 to disable refraction correction.
        temperature: Atmospheric temperature at observer in degrees Celsius.
                     Default is 15.0.
        lapse_rate: Temperature lapse rate dT/dh in degrees Kelvin per meter.
                    Default is None (uses global value from set_lapse_rate(),
                    or 0.0065 K/m if not set).
                    Typical values range from 0.0034 to 0.010 K/m.
        calc_flag: Direction of conversion:
            - SE_TRUE_TO_APP (0): Convert true altitude to apparent altitude
            - SE_APP_TO_TRUE (1): Convert apparent altitude to true altitude

    Returns:
        A tuple of (converted_altitude, details) where:
        - converted_altitude: The converted altitude in degrees
        - details: A tuple of 4 floats:
            - [0]: True altitude (input or computed)
            - [1]: Apparent altitude (input or computed)
            - [2]: Refraction amount in degrees
            - [3]: Dip of the horizon in degrees (negative, as horizon dips below
                   geometric horizontal for elevated observers)

    Notes:
        - The dip of the horizon accounts for both geometric effects and
          atmospheric refraction. An elevated observer sees the horizon
          below the geometric horizontal plane.
        - The lapse rate affects the dip calculation through atmospheric
          refraction near the horizon. Higher lapse rates result in less
          refraction of the horizon, yielding more negative dip values.
        - Standard atmospheric lapse rate is 0.0065 K/m (6.5°C per 1000m).

    Examples:
        >>> # At sea level, no dip of horizon
        >>> alt, (true, app, ref, dip) = refrac_extended(0.0, 0.0)
        >>> round(ref, 4)
        0.4721
        >>> round(dip, 4)
        0.0

        >>> # At 1000m elevation, horizon dips about 0.88 degrees
        >>> alt, (true, app, ref, dip) = refrac_extended(0.0, 1000.0)
        >>> round(dip, 2)
        -0.88
    """
    from .state import get_lapse_rate

    # Use global lapse rate if none provided
    if lapse_rate is None:
        lapse_rate = get_lapse_rate()

    # Earth's radius in meters
    EARTH_RADIUS = 6371000.0

    # Calculate refraction using the base refrac function
    if calc_flag == SE_TRUE_TO_APP:
        true_alt = altitude
        apparent_alt = refrac(altitude, pressure, temperature, SE_TRUE_TO_APP)
        refraction = apparent_alt - true_alt
    else:
        apparent_alt = altitude
        true_alt = refrac(altitude, pressure, temperature, SE_APP_TO_TRUE)
        refraction = apparent_alt - true_alt

    # Calculate dip of the horizon for elevated observers
    if altitude_geo <= 0:
        dip = 0.0
    else:
        # Geometric dip angle (without atmospheric refraction)
        # dip_geometric = arccos(R / (R + h)) ≈ sqrt(2h/R) for small h
        # Using the more accurate formula:
        ratio = EARTH_RADIUS / (EARTH_RADIUS + altitude_geo)
        if ratio >= 1.0:
            dip_geometric_rad = 0.0
        else:
            dip_geometric_rad = math.acos(ratio)

        dip_geometric = math.degrees(dip_geometric_rad)

        # Atmospheric refraction correction for the dip
        # The ray from the observer to the visible horizon is bent by
        # atmospheric refraction. This depends on the lapse rate.
        #
        # The refraction correction reduces the apparent dip (horizon appears
        # higher due to refraction).
        #
        # Based on empirical fitting to Swiss Ephemeris results, the
        # refraction coefficient k follows:
        # k = 0.1117 + 3.5516 * lapse_rate
        # where observed_dip = -dip_geometric * (1 - k)

        if lapse_rate > 0:
            # Calculate the refraction coefficient based on lapse rate
            refraction_coef = 0.1117 + 3.5516 * lapse_rate
        else:
            # If lapse_rate is 0 or negative, use geometric dip only
            refraction_coef = 0.0

        # Apply atmospheric correction: observed dip is less than geometric dip
        # due to refraction bending light downward
        dip = -dip_geometric * (1.0 - refraction_coef)

    # Return the converted altitude and detail tuple
    if calc_flag == SE_TRUE_TO_APP:
        return (apparent_alt, (true_alt, apparent_alt, refraction, dip))
    else:
        return (true_alt, (true_alt, apparent_alt, refraction, dip))


def cotrans(
    coord: Tuple[float, float, float], obliquity: float
) -> Tuple[float, float, float]:
    """
    Transform coordinates between ecliptic and equatorial systems.

    Compatible with pyswisseph's swe.cotrans() function.

    The direction of transformation depends on the sign of obliquity:
    - Positive obliquity: ecliptic (lon, lat) → equatorial (RA, Dec)
    - Negative obliquity: equatorial (RA, Dec) → ecliptic (lon, lat)

    Args:
        coord: Tuple of (longitude/RA, latitude/Dec, distance) in degrees
        obliquity: Obliquity of the ecliptic in degrees.
                   Positive for ecliptic→equatorial, negative for equatorial→ecliptic.

    Returns:
        Tuple of (transformed_lon/RA, transformed_lat/Dec, distance)
        Distance is unchanged by the transformation.

    Examples:
        >>> # Ecliptic to equatorial (positive obliquity)
        >>> cotrans((0.0, 0.0, 1.0), 23.4)
        (0.0, 0.0, 1.0)
        >>> # Equatorial to ecliptic (negative obliquity)
        >>> cotrans((0.0, 0.0, 1.0), -23.4)
        (0.0, 0.0, 1.0)
    """
    lon = coord[0]
    lat = coord[1]
    dist = coord[2]

    # Convert to radians
    lon_rad = math.radians(lon)
    lat_rad = math.radians(lat)
    eps_rad = math.radians(obliquity)

    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # Calculate the new latitude (Dec for ecl→eq, β for eq→ecl)
    # sin(new_lat) = sin(lat) * cos(eps) - cos(lat) * sin(eps) * sin(lon)
    sin_new_lat = sin_lat * cos_eps - cos_lat * sin_eps * sin_lon
    # Clamp to avoid numerical issues with asin
    sin_new_lat = max(-1.0, min(1.0, sin_new_lat))
    new_lat_rad = math.asin(sin_new_lat)

    # Calculate the new longitude (RA for ecl→eq, λ for eq→ecl)
    # tan(new_lon) = (sin(lon) * cos(eps) + tan(lat) * sin(eps)) / cos(lon)
    tan_lat = math.tan(lat_rad)
    y = sin_lon * cos_eps + tan_lat * sin_eps
    x = cos_lon

    new_lon_rad = math.atan2(y, x)

    # Convert back to degrees
    new_lon = math.degrees(new_lon_rad)
    new_lat = math.degrees(new_lat_rad)

    # Normalize longitude to [0, 360)
    new_lon = new_lon % 360.0

    return (new_lon, new_lat, dist)


def degnorm(angle: float) -> float:
    """
    Normalize an angle to the range [0, 360).

    Compatible with pyswisseph's swe.degnorm() function.
    Equivalent to `angle % 360` but correctly handles negative numbers.

    Args:
        angle: Angle in degrees (any value)

    Returns:
        Normalized angle in range [0, 360)

    Examples:
        >>> degnorm(45)
        45.0
        >>> degnorm(-45)
        315.0
        >>> degnorm(370)
        10.0
        >>> degnorm(-370)
        350.0
        >>> degnorm(360)
        0.0
        >>> degnorm(720)
        0.0
    """
    return angle % 360.0


TWO_PI = 2.0 * math.pi


def radnorm(angle: float) -> float:
    """
    Normalize an angle to the range [0, 2*pi).

    Compatible with pyswisseph's swe.radnorm() function.
    Equivalent to `angle % (2*pi)` but correctly handles negative numbers.

    Args:
        angle: Angle in radians (any value)

    Returns:
        Normalized angle in range [0, 2*pi)

    Examples:
        >>> import math
        >>> radnorm(math.pi / 4)  # 45 degrees
        0.7853981633974483
        >>> radnorm(-math.pi / 4)  # -45 degrees -> 315 degrees
        5.497787143782138
        >>> radnorm(3 * math.pi)  # 540 degrees -> 180 degrees
        3.141592653589793
        >>> radnorm(-3 * math.pi)  # -540 degrees -> 180 degrees
        3.141592653589793
        >>> radnorm(2 * math.pi)  # 360 degrees -> 0
        0.0
        >>> radnorm(4 * math.pi)  # 720 degrees -> 0
        0.0
    """
    return angle % TWO_PI


def difdeg2n(p1: float, p2: float) -> float:
    """
    Calculate distance in degrees p1 - p2 normalized to [-180;180].

    Compatible with pyswisseph's swe.difdeg2n() function.
    Computes the signed angular difference, handling 360° wrapping.

    Args:
        p1: First angle in degrees
        p2: Second angle in degrees

    Returns:
        Normalized difference in range [-180, 180]

    Examples:
        >>> difdeg2n(10, 20)
        -10.0
        >>> difdeg2n(350, 10)
        -20.0
        >>> difdeg2n(10, 350)
        20.0
        >>> difdeg2n(180, 0)
        180.0
    """
    diff = (p1 - p2) % 360.0
    if diff > 180.0:
        diff -= 360.0
    return diff


def difdegn(p1: float, p2: float) -> float:
    """
    Calculate distance in degrees p1 - p2 normalized to [0, 360).

    Compatible with pyswisseph's swe.difdegn() function.
    Computes the difference between two angles, always returning a positive
    value in the range [0, 360). Unlike difdeg2n() which returns [-180, 180],
    this function always returns a positive value.

    Args:
        p1: First angle in degrees
        p2: Second angle in degrees

    Returns:
        Normalized difference in range [0, 360)

    Examples:
        >>> difdegn(10, 20)
        350.0
        >>> difdegn(20, 10)
        10.0
        >>> difdegn(350, 10)
        340.0
        >>> difdegn(10, 350)
        20.0
        >>> difdegn(0, 0)
        0.0
    """
    return (p1 - p2) % 360.0


def difrad2n(p1: float, p2: float) -> float:
    """
    Calculate distance in radians p1 - p2 normalized to [-π, π].

    This is the radians equivalent of difdeg2n().
    Computes the signed angular difference, handling 2π wrapping.

    Args:
        p1: First angle in radians
        p2: Second angle in radians

    Returns:
        Normalized difference in range [-π, π]

    Examples:
        >>> import math
        >>> difrad2n(0.1, 0.2)
        -0.1
        >>> difrad2n(6.0, 0.2)  # Near 2π wraparound
        -0.48318530717958...
        >>> difrad2n(0.2, 6.0)  # Opposite direction
        0.48318530717958...
        >>> difrad2n(math.pi, 0)
        -3.141592653589793
    """
    diff = (p1 - p2) % TWO_PI
    if diff >= math.pi:
        diff -= TWO_PI
    return diff


# Constants for centiseconds calculations
# 1 centisecond = 1/100 arcsecond
# 360 degrees = 360 * 3600 * 100 = 129,600,000 centiseconds
CS360 = 360 * 3600 * 100  # 129600000 centiseconds in a full circle
CS180 = 180 * 3600 * 100  # 64800000 centiseconds in a half circle


def difcs2n(a: int, b: int) -> int:
    """
    Calculate distance in centiseconds a - b normalized to [-180°, +180°].

    This function computes the signed angular difference between two angles
    expressed in centiseconds (1/100 of an arcsecond), handling 360° wrapping.
    The result is normalized to the equivalent of [-180°, +180°] in centiseconds.

    Compatible with pyswisseph's swe.difcs2n() function.

    Args:
        a: First angle in centiseconds
        b: Second angle in centiseconds

    Returns:
        Normalized difference in range [-64800000, 64800000) centiseconds
        (equivalent to [-180°, +180°))

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - 360° = 129,600,000 centiseconds
        - 180° = 64,800,000 centiseconds
        - When the difference is exactly 180°, returns -180° (matching pyswisseph)

    Examples:
        >>> difcs2n(360000, 720000)  # 1° - 2° = -1° = -360000 cs
        -360000
        >>> difcs2n(129240000, 360000)  # 359° - 1° = -2° = -720000 cs (shorter path)
        -720000
        >>> difcs2n(360000, 129240000)  # 1° - 359° = 2° = 720000 cs (shorter path)
        720000
        >>> difcs2n(64800000, 0)  # 180° - 0° = -180° (pyswisseph convention)
        -64800000
    """
    diff = (a - b) % CS360
    if diff >= CS180:
        diff -= CS360
    return diff


def difcsn(a: int, b: int) -> int:
    """
    Calculate distance in centiseconds a - b normalized to [0, 360°).

    This function computes the angular difference between two angles
    expressed in centiseconds (1/100 of an arcsecond), handling 360° wrapping.
    The result is normalized to the equivalent of [0, 360°) in centiseconds,
    always returning a non-negative value.

    Compatible with pyswisseph's swe.difcsn() function.

    Args:
        a: First angle in centiseconds
        b: Second angle in centiseconds

    Returns:
        Normalized difference in range [0, 129600000) centiseconds
        (equivalent to [0°, 360°))

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - 360° = 129,600,000 centiseconds
        - Unlike difcs2n() which returns [-180°, +180°], this function
          always returns a positive value in [0°, 360°)

    Examples:
        >>> difcsn(720000, 360000)  # 2° - 1° = 1° = 360000 cs
        360000
        >>> difcsn(360000, 720000)  # 1° - 2° = 359° = 129240000 cs (positive)
        129240000
        >>> difcsn(360000, 129240000)  # 1° - 359° = 2° = 720000 cs
        720000
        >>> difcsn(0, 0)  # 0° - 0° = 0°
        0
    """
    return (a - b) % CS360


def csnorm(cs: int) -> int:
    """
    Normalize a value in centiseconds to the range [0, 360°).

    This function normalizes an angle expressed in centiseconds (1/100 of an
    arcsecond) to the equivalent of [0°, 360°) in centiseconds, i.e., the range
    [0, 129600000). Correctly handles negative numbers.

    Compatible with pyswisseph's swe.csnorm() function.

    Args:
        cs: Angle in centiseconds (any value)

    Returns:
        Normalized angle in range [0, 129600000) centiseconds
        (equivalent to [0°, 360°))

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - 360° = 129,600,000 centiseconds
        - This is the centisecond equivalent of degnorm() for degrees
          and radnorm() for radians

    Examples:
        >>> csnorm(360000)  # 1° stays as 1°
        360000
        >>> csnorm(-360000)  # -1° -> 359° = 129240000 cs
        129240000
        >>> csnorm(129960000)  # 361° -> 1° = 360000 cs
        360000
        >>> csnorm(-129960000)  # -361° -> 359° = 129240000 cs
        129240000
        >>> csnorm(129600000)  # 360° -> 0°
        0
        >>> csnorm(0)
        0
    """
    return cs % CS360


def csroundsec(cs: int) -> int:
    """
    Round a value in centiseconds to the nearest arcsecond.

    This function rounds a centisecond value to the nearest arcsecond,
    returning the result still in centiseconds (i.e., rounded to the
    nearest multiple of 100).

    Compatible with pyswisseph's swe.csroundsec() function.

    Args:
        cs: Angle in centiseconds (1/100 arcsecond)

    Returns:
        Centisecond value rounded to the nearest arcsecond (multiple of 100)

    Notes:
        - 1 centisecond = 1/100 arcsecond
        - Uses C-style integer division (truncation toward zero) with +50 offset
        - Special handling at degree boundaries (30° for positive, 90° for negative)

    Examples:
        >>> csroundsec(150)  # 1.50 arcseconds -> 2 arcseconds = 200 cs
        200
        >>> csroundsec(149)  # 1.49 arcseconds -> 1 arcsecond = 100 cs
        100
        >>> csroundsec(50)   # 0.50 arcseconds -> 1 arcsecond = 100 cs
        100
        >>> csroundsec(100)  # 1.00 arcsecond -> 1 arcsecond = 100 cs
        100
        >>> csroundsec(0)    # 0 centiseconds -> 0 arcseconds
        0
        >>> csroundsec(-150) # -1.50 arcseconds -> -1 arcsecond = -100 cs
        -100
    """
    # C-style integer division: truncates toward zero
    # For cs + 50, we need to use truncation toward zero, not floor division
    cs_plus_50 = cs + 50
    if cs_plus_50 >= 0:
        result = (cs_plus_50 // 100) * 100
    else:
        # C-style truncation toward zero for negative values
        result = -((-cs_plus_50) // 100) * 100

    # Special case for positive values at 30-degree boundaries (10800000 cs = 30°)
    # Values just below these boundaries round down instead of up
    if cs > 0 and result % 10800000 == 0 and result != 0 and cs < result:
        return result - 100

    # Special case for negative values at 90-degree boundaries (32400000 cs = 90°)
    if cs < 0:
        if result != 0 and result % 32400000 == 0 and cs <= result - 100:
            return result - 100
        # Values in (-100, 0) round to 0, values <= -100 round to -100 when formula gives 0
        if cs <= -100 and result == 0:
            return -100

    return result


def cs2degstr(cs: int) -> str:
    """
    Convert a value in centiseconds to a formatted degrees string.

    This function converts an angular measurement in centiseconds (1/100 of an
    arcsecond) to a human-readable string in the format "D°M'S.ss\"" where D is
    degrees, M is arcminutes, S is arcseconds, and ss is centiseconds (hundredths
    of arcsecond).

    Compatible with pyswisseph's swe.cs2degstr() function.

    Args:
        cs: Angle in centiseconds (any integer value)

    Returns:
        Formatted string representing the angle in degrees, minutes, seconds.
        Format: "DDD°MM'SS.ss\"" (e.g., "123°45'06.78\"")

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - 360° = 129,600,000 centiseconds
        - Negative values produce negative degree strings (e.g., "-1°00'00.00\"")
        - The seconds field includes two decimal places for centiseconds

    Examples:
        >>> cs2degstr(0)
        '  0° 0\\' 0.00"'
        >>> cs2degstr(360000)  # 1 degree
        '  1° 0\\' 0.00"'
        >>> cs2degstr(3723456)  # 10°20'34.56"
        ' 10°20\\'34.56"'
        >>> cs2degstr(-360000)  # -1 degree
        ' -1° 0\\' 0.00"'
    """
    # Handle sign
    if cs < 0:
        sign = -1
        cs = -cs
    else:
        sign = 1

    # Extract degrees, minutes, seconds, and centiseconds
    # 1 degree = 3600 * 100 = 360000 centiseconds
    # 1 minute = 60 * 100 = 6000 centiseconds
    # 1 second = 100 centiseconds
    degrees = cs // 360000
    remainder = cs % 360000
    minutes = remainder // 6000
    remainder = remainder % 6000
    seconds = remainder // 100
    centisecs = remainder % 100

    # Apply sign to degrees
    if sign < 0:
        degrees = -degrees

    # Format the string matching Swiss Ephemeris format
    # Format: "%3d°%2d'%2d.%02d"" with proper spacing
    return f"{degrees:3d}°{minutes:2d}'{seconds:2d}.{centisecs:02d}\""


def cs2lonlatstr(cs: int, plus_char: str, minus_char: str) -> str:
    """
    Convert a value in centiseconds to a formatted longitude/latitude string.

    This function converts an angular measurement in centiseconds (1/100 of an
    arcsecond) to a human-readable string with a directional character suffix.
    The format is "D°M'S\" X" where D is degrees, M is arcminutes, S is arcseconds,
    and X is the directional character (e.g., N/S for latitude, E/W for longitude).

    Compatible with pyswisseph's swe.cs2lonlatstr() function.

    Args:
        cs: Angle in centiseconds (any integer value)
        plus_char: Character to use for positive values (e.g., "N" or "E")
        minus_char: Character to use for negative values (e.g., "S" or "W")

    Returns:
        Formatted string representing the angle with directional character.
        Format: "DDD°MM'SS\" X" (e.g., "45°30'00\" N" or "122°15'30\" W")

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - Positive values use plus_char, negative values use minus_char
        - Seconds are rounded to whole numbers (no centiseconds shown)
        - Commonly used for geographic coordinates:
          - Latitude: plus_char="N", minus_char="S"
          - Longitude: plus_char="E", minus_char="W"

    Examples:
        >>> cs2lonlatstr(163800000, "N", "S")  # 45°30'00" North
        ' 45°30\\' 0" N'
        >>> cs2lonlatstr(-163800000, "N", "S")  # 45°30'00" South
        ' 45°30\\' 0" S'
        >>> cs2lonlatstr(440154000, "E", "W")  # 122°15'30" East
        '122°15\\'30" E'
        >>> cs2lonlatstr(-440154000, "E", "W")  # 122°15'30" West
        '122°15\\'30" W'
    """
    # Determine direction character based on sign
    if cs >= 0:
        direction = plus_char
    else:
        direction = minus_char
        cs = -cs

    # Extract degrees, minutes, and seconds
    # 1 degree = 3600 * 100 = 360000 centiseconds
    # 1 minute = 60 * 100 = 6000 centiseconds
    # 1 second = 100 centiseconds
    degrees = cs // 360000
    remainder = cs % 360000
    minutes = remainder // 6000
    remainder = remainder % 6000
    # Round centiseconds to nearest arcsecond
    seconds = (remainder + 50) // 100

    # Handle carry-over from rounding
    if seconds >= 60:
        seconds = 0
        minutes += 1
        if minutes >= 60:
            minutes = 0
            degrees += 1

    # Format the string matching Swiss Ephemeris format
    # Format: "%3d°%2d'%2d\" %c" with proper spacing
    return f"{degrees:3d}°{minutes:2d}'{seconds:2d}\" {direction}"


def cs2timestr(cs: int) -> str:
    """
    Convert a value in centiseconds to a formatted time string.

    This function converts a time measurement in centiseconds (1/100 of a second)
    to a human-readable string in the format "HH:MM:SS" where HH is hours, MM is
    minutes, and SS is seconds.

    Compatible with pyswisseph's swe.cs2timestr() function.

    Args:
        cs: Time in centiseconds (any integer value)

    Returns:
        Formatted string representing the time in hours, minutes, seconds.
        Format: "HH:MM:SS" (e.g., "12:34:56")

    Notes:
        - 1 centisecond = 1/100 second
        - 1 second = 100 centiseconds
        - 1 minute = 6000 centiseconds
        - 1 hour = 360000 centiseconds
        - Negative values produce negative hour strings (e.g., "-1:00:00")
        - Seconds are rounded to whole numbers (centiseconds are rounded)

    Examples:
        >>> cs2timestr(0)
        ' 0:00:00'
        >>> cs2timestr(360000)  # 1 hour
        ' 1:00:00'
        >>> cs2timestr(4526050)  # 12:34:21 (with rounding from .50)
        '12:34:21'
        >>> cs2timestr(-360000)  # -1 hour
        '-1:00:00'
    """
    # Handle sign
    if cs < 0:
        sign = -1
        cs = -cs
    else:
        sign = 1

    # Extract hours, minutes, and seconds
    # 1 hour = 60 * 60 * 100 = 360000 centiseconds
    # 1 minute = 60 * 100 = 6000 centiseconds
    # 1 second = 100 centiseconds
    hours = cs // 360000
    remainder = cs % 360000
    minutes = remainder // 6000
    remainder = remainder % 6000
    # Round centiseconds to nearest second
    seconds = (remainder + 50) // 100

    # Handle carry-over from rounding
    if seconds >= 60:
        seconds = 0
        minutes += 1
        if minutes >= 60:
            minutes = 0
            hours += 1

    # Apply sign to hours
    if sign < 0:
        hours = -hours

    # Format the string matching Swiss Ephemeris format
    # Format: "%2d:%02d:%02d" with proper spacing
    return f"{hours:2d}:{minutes:02d}:{seconds:02d}"


def deg_midp(a: float, b: float) -> float:
    """
    Calculate the midpoint between two angles in degrees.

    Handles wraparound at 360° correctly by finding the midpoint along the
    shorter arc between the two angles. For example, the midpoint between
    350° and 10° is 0° (or equivalently 360°), not 180°.

    Args:
        a: First angle in degrees (any value, will be normalized)
        b: Second angle in degrees (any value, will be normalized)

    Returns:
        Midpoint angle in range [0, 360)

    Examples:
        >>> deg_midp(0, 90)
        45.0
        >>> deg_midp(350, 10)
        0.0
        >>> deg_midp(10, 350)
        0.0
        >>> deg_midp(180, 0)
        90.0
        >>> deg_midp(170, 190)
        180.0
        >>> deg_midp(-10, 10)
        0.0
    """
    # Normalize both angles to [0, 360)
    a = a % 360.0
    b = b % 360.0

    # Calculate the difference
    diff = b - a

    # If the absolute difference is greater than 180, go the shorter way
    if diff > 180.0:
        diff -= 360.0
    elif diff < -180.0:
        diff += 360.0

    # Calculate midpoint along the shorter arc
    midp = a + diff / 2.0

    # Normalize result to [0, 360)
    return midp % 360.0


def rad_midp(a: float, b: float) -> float:
    """
    Calculate the midpoint between two angles in radians.

    Handles wraparound at 2*pi correctly by finding the midpoint along the
    shorter arc between the two angles. For example, the midpoint between
    5.5 radians and 0.5 radians is 0 (or equivalently 2*pi), not pi.

    Args:
        a: First angle in radians (any value, will be normalized)
        b: Second angle in radians (any value, will be normalized)

    Returns:
        Midpoint angle in range [0, 2*pi)

    Examples:
        >>> import math
        >>> rad_midp(0, math.pi / 2)  # 0 to 90 degrees -> 45 degrees
        0.7853981633974483
        >>> rad_midp(5.5, 0.5)  # Near 2*pi wraparound
        0.0  # approximately
        >>> rad_midp(math.pi, 0)  # 180 to 0 -> 90 degrees
        1.5707963267948966
    """
    # Normalize both angles to [0, 2*pi)
    a = a % TWO_PI
    b = b % TWO_PI

    # Calculate the difference
    diff = b - a

    # If the absolute difference is greater than pi, go the shorter way
    if diff > math.pi:
        diff -= TWO_PI
    elif diff < -math.pi:
        diff += TWO_PI

    # Calculate midpoint along the shorter arc
    midp = a + diff / 2.0

    # Normalize result to [0, 2*pi)
    return midp % TWO_PI


def d2l(value: float) -> int:
    """
    Convert a double (float) to a long integer with rounding.

    This function rounds a floating-point number to the nearest integer using
    "round half away from zero" semantics (also known as commercial rounding).
    This is the behavior used by the Swiss Ephemeris library internally.

    Compatible with pyswisseph's swe.d2l() function.

    Args:
        value: A floating-point number to convert.

    Returns:
        The nearest integer to value. For values exactly halfway between two
        integers, rounds away from zero (e.g., 0.5 -> 1, -0.5 -> -1).

    Notes:
        - This differs from Python's built-in round() function, which uses
          "round half to even" (banker's rounding) for Python 3.
        - Swiss Ephemeris uses this for internal conversions, but also exposes
          it publicly for consistency.

    Examples:
        >>> d2l(1.4)
        1
        >>> d2l(1.5)
        2
        >>> d2l(1.6)
        2
        >>> d2l(-1.4)
        -1
        >>> d2l(-1.5)
        -2
        >>> d2l(-1.6)
        -2
        >>> d2l(0.5)
        1
        >>> d2l(-0.5)
        -1
        >>> d2l(2.5)
        3
        >>> d2l(-2.5)
        -3
    """
    if value >= 0:
        return int(value + 0.5)
    else:
        return int(value - 0.5)


# Split degree flags (imported from constants for convenience)
SPLIT_DEG_ROUND_SEC: int = 1  # Round to seconds
SPLIT_DEG_ROUND_MIN: int = 2  # Round to minutes
SPLIT_DEG_ROUND_DEG: int = 4  # Round to degrees
SPLIT_DEG_ZODIACAL: int = 8  # Return zodiac sign number (0-11)
SPLIT_DEG_NAKSHATRA: int = 1024  # Return nakshatra number (0-26)
SPLIT_DEG_KEEP_SIGN: int = 16  # Don't round to next zodiac sign/nakshatra
SPLIT_DEG_KEEP_DEG: int = 32  # Don't round to next degree


def split_deg(degrees: float, roundflag: int = 0) -> Tuple[int, int, int, float, int]:
    """
    Split a degree value into sign/nakshatra, degrees, minutes, seconds, and fraction.

    This function decomposes an ecliptic longitude (or any angle in degrees) into
    its constituent parts: zodiac sign (or nakshatra), degrees within that sign,
    minutes, seconds, and fraction of second.

    Compatible with pyswisseph's swe.split_deg() function.

    Args:
        degrees: Position in decimal degrees (can be negative or > 360)
        roundflag: Bit flags combination indicating how to round and format:
            - 0: Don't round, return all components with full precision
            - SPLIT_DEG_ROUND_SEC (1): Round to nearest second
            - SPLIT_DEG_ROUND_MIN (2): Round to nearest minute
            - SPLIT_DEG_ROUND_DEG (4): Round to nearest degree
            - SPLIT_DEG_ZODIACAL (8): Return zodiac sign number (0-11, each 30 deg)
            - SPLIT_DEG_NAKSHATRA (1024): Return nakshatra number (0-26, each 13deg20')
            - SPLIT_DEG_KEEP_SIGN (16): Don't round to next zodiac sign/nakshatra
            - SPLIT_DEG_KEEP_DEG (32): Don't round to next degree

    Returns:
        Tuple of (deg, min, sec, secfr, sign) where:
        - deg: Degrees within sign (0-29 with ZODIACAL, 0-13 with NAKSHATRA,
               or total degrees without either flag)
        - min: Arc minutes (0-59)
        - sec: Arc seconds (0-59)
        - secfr: Fraction of arc second (0.0-0.999...), or the rounded seconds
                 value when rounding flags are used
        - sign: Zodiac sign (0-11), nakshatra (0-26), or +1/-1 for positive/negative

    Notes:
        - Without ZODIACAL or NAKSHATRA flag, sign returns +1 (positive) or -1 (negative)
        - With ZODIACAL flag: uses absolute value, sign is 0-11
          (0=Aries, 1=Taurus, ..., 11=Pisces)
        - With NAKSHATRA flag: uses absolute value, sign is 0-26
          (each nakshatra spans 13 deg 20 min = 13.333... degrees)
        - Rounding flags affect how values are truncated/rounded
        - KEEP_SIGN prevents rounding from advancing to the next sign
        - KEEP_DEG prevents rounding from advancing to the next degree

    Examples:
        >>> split_deg(45.5, 0)  # No flags: 45deg 30min
        (45, 30, 0, 0.0, 1)
        >>> split_deg(45.5, SPLIT_DEG_ZODIACAL)  # With zodiac: 15deg 30min Taurus
        (15, 30, 0, 0.0, 1)
        >>> split_deg(-30.5, 0)  # Negative: sign = -1
        (30, 30, 0, 0.0, -1)
        >>> split_deg(-30.5, SPLIT_DEG_ZODIACAL)  # Negative with zodiac: uses absolute
        (0, 30, 0, 0.0, 1)
    """
    # Handle sign/negative values
    is_negative = degrees < 0

    # Determine if we should use zodiacal or nakshatra division
    use_zodiacal = bool(roundflag & SPLIT_DEG_ZODIACAL)
    use_nakshatra = bool(roundflag & SPLIT_DEG_NAKSHATRA)
    round_sec = bool(roundflag & SPLIT_DEG_ROUND_SEC)
    round_min = bool(roundflag & SPLIT_DEG_ROUND_MIN)
    round_deg = bool(roundflag & SPLIT_DEG_ROUND_DEG)
    keep_sign = bool(roundflag & SPLIT_DEG_KEEP_SIGN)
    keep_deg = bool(roundflag & SPLIT_DEG_KEEP_DEG)

    # For zodiacal/nakshatra modes, pyswisseph uses absolute value
    if use_zodiacal or use_nakshatra:
        # Use absolute value (pyswisseph behavior for negative values)
        deg_abs = abs(degrees)

        if use_nakshatra:
            # Nakshatra: 27 divisions, each 13°20' = 13.333... degrees
            nakshatra_size = 360.0 / 27.0  # 13.333...
            sign_num = int(deg_abs / nakshatra_size)
            deg_in_sign = deg_abs - (sign_num * nakshatra_size)
            # Handle wraparound for values >= 360
            if sign_num >= 27:
                sign_num = sign_num % 27
        else:
            # Zodiacal: 12 signs, each 30 degrees
            sign_num = int(deg_abs / 30.0)
            deg_in_sign = deg_abs - (sign_num * 30.0)
            # Handle wraparound for values >= 360
            if sign_num >= 12:
                sign_num = sign_num % 12
    else:
        # No zodiacal/nakshatra: use absolute degrees with sign indicator
        deg_in_sign = abs(degrees)
        sign_num = 1 if not is_negative else -1

    # Split into degrees, minutes, seconds, fraction
    deg_part = int(deg_in_sign)
    remainder = (deg_in_sign - deg_part) * 60.0
    min_part = int(remainder)
    remainder = (remainder - min_part) * 60.0
    sec_part = int(remainder)
    secfr = remainder - sec_part

    # Store original values for ROUND_MIN and ROUND_DEG behavior
    orig_sec_part = sec_part
    orig_min_part = min_part

    # Apply rounding if requested
    if round_sec or round_min or round_deg:
        if round_sec:
            # Round to nearest second
            if secfr >= 0.5:
                sec_part += 1
            secfr = float(sec_part)  # pyswisseph copies sec value to secfr

            # Handle overflow: sec >= 60
            if sec_part >= 60:
                sec_part = 0
                min_part += 1
                # Handle overflow: min >= 60
                if min_part >= 60:
                    min_part = 0
                    deg_part += 1

        elif round_min:
            # Round to nearest minute
            # Check if seconds >= 30 to round up
            total_sec = orig_sec_part + secfr
            if total_sec >= 30.0:
                min_part += 1
            # pyswisseph: sec = original sec, secfr = original sec
            sec_part = orig_sec_part
            secfr = float(orig_sec_part)

            # Handle overflow: min >= 60
            if min_part >= 60:
                min_part = 0
                deg_part += 1

        elif round_deg:
            # Round to nearest degree
            # Check if minutes >= 30 to round up
            total_min = orig_min_part + (orig_sec_part + secfr) / 60.0
            if total_min >= 30.0:
                deg_part += 1
            # pyswisseph: min = 0 (not original), sec = original sec, secfr = original sec
            min_part = 0
            sec_part = orig_sec_part
            secfr = float(orig_sec_part)

        # Handle degree overflow for zodiacal/nakshatra
        if use_zodiacal or use_nakshatra:
            max_deg = 30 if use_zodiacal else 14
            max_signs = 12 if use_zodiacal else 27

            if deg_part >= max_deg:
                if keep_sign:
                    # Don't advance to next sign - cap at max
                    deg_part = max_deg - 1 if use_zodiacal else 13
                    min_part = 59
                    sec_part = 59
                    secfr = 59.0
                else:
                    # Advance to next sign
                    deg_part = 0
                    sign_num += 1
                    if sign_num >= max_signs:
                        sign_num = 0

    return (deg_part, min_part, sec_part, secfr, sign_num)


def swe_calc_angles(jd_ut: float, lat: float, lon: float):
    """
    Pre-calculate and cache astrological angles and planet positions
    for use with Arabic parts.

    Args:
        jd_ut: Julian Day (UT)
        lat: Latitude (degrees)
        lon: Longitude (degrees)

    Returns:
        Dictionary with calculated positions
    """
    from .state import set_angles_cache, set_topo
    from .angles import calc_angles
    from .planets import swe_calc_ut
    from .constants import SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS

    # Set observer location
    set_topo(lon, lat, 0)

    # Calculate angles
    angles_dict = calc_angles(jd_ut, lat, lon)

    # Calculate and add planet positions for Arabic parts
    sun_pos, _ = swe_calc_ut(jd_ut, SE_SUN, 0)
    moon_pos, _ = swe_calc_ut(jd_ut, SE_MOON, 0)
    mercury_pos, _ = swe_calc_ut(jd_ut, SE_MERCURY, 0)
    venus_pos, _ = swe_calc_ut(jd_ut, SE_VENUS, 0)

    angles_dict["Sun"] = sun_pos[0]
    angles_dict["Moon"] = moon_pos[0]
    angles_dict["Mercury"] = mercury_pos[0]
    angles_dict["Venus"] = venus_pos[0]

    # Cache for Arabic parts
    set_angles_cache(angles_dict)

    return angles_dict
