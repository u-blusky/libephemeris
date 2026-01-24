"""
Utility functions for libephemeris.

Provides helper functions compatible with pyswisseph API including
angular calculations and other mathematical utilities.
"""

import math
from typing import Tuple

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
    lapse_rate: float = 0.0065,
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
                    Default is 0.0065 K/m (standard atmosphere).
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
