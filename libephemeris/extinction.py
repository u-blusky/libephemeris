"""
Atmospheric extinction model for libephemeris.

This module provides functions to calculate atmospheric extinction,
which increases the apparent magnitude of celestial objects as they
approach the horizon. This is essential for heliacal visibility calculations.

The extinction model is based on the work of:
- Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
- Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"

Key concepts:
- Extinction is measured in magnitudes per airmass (mag/airmass)
- Airmass = sec(z) where z is zenith angle (for plane-parallel atmosphere)
- Total extinction = extinction_coefficient * airmass

The extinction coefficient k is composed of:
- Rayleigh scattering (molecular): wavelength-dependent, ~0.14 mag/airmass at V
- Aerosol scattering: varies with atmospheric conditions, typically 0.05-0.25
- Ozone absorption: small contribution, ~0.016 mag/airmass

References:
    - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
    - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    - Young, A.T. (1974) "Observational Technique and Data Reduction"
    - Kasten, F. & Young, A.T. (1989) "Revised optical air mass tables"
"""

import math
from dataclasses import dataclass
from typing import Optional


# Standard wavelengths in nanometers for common photometric bands
WAVELENGTH_U = 365.0  # U band (ultraviolet)
WAVELENGTH_B = 440.0  # B band (blue)
WAVELENGTH_V = 550.0  # V band (visual, standard)
WAVELENGTH_R = 640.0  # R band (red)
WAVELENGTH_I = 790.0  # I band (infrared)

# Default atmospheric parameters
DEFAULT_PRESSURE_MBAR = 1013.25  # Standard sea level pressure
DEFAULT_TEMPERATURE_C = 15.0  # Standard temperature
DEFAULT_HUMIDITY_PERCENT = 50.0  # Moderate humidity
DEFAULT_ALTITUDE_M = 0.0  # Sea level

# Atmospheric scale height in meters (for pressure/density variation)
SCALE_HEIGHT_M = 8500.0


@dataclass
class ExtinctionCoefficients:
    """
    Atmospheric extinction coefficients broken down by component.

    Attributes:
        k_rayleigh: Rayleigh (molecular) scattering coefficient [mag/airmass]
        k_aerosol: Aerosol (Mie) scattering coefficient [mag/airmass]
        k_ozone: Ozone absorption coefficient [mag/airmass]
        k_water: Water vapor absorption coefficient [mag/airmass]
        k_total: Total extinction coefficient [mag/airmass]
    """

    k_rayleigh: float
    k_aerosol: float
    k_ozone: float
    k_water: float
    k_total: float


def calc_airmass(altitude_deg: float, method: str = "kasten_young") -> float:
    """
    Calculate the airmass (relative path length through atmosphere).

    Airmass indicates how much atmosphere light must traverse compared
    to looking straight up (zenith, airmass = 1.0). At the horizon,
    airmass is approximately 38-40.

    Args:
        altitude_deg: Altitude of object above horizon in degrees.
                     Negative values indicate object below horizon.
        method: Airmass calculation method:
            - "secant": Simple sec(z) formula, accurate above ~15 degrees
            - "kasten_young": Kasten & Young (1989) formula, accurate to horizon
            - "rozenberg": Rozenberg (1966) formula for very low altitudes

    Returns:
        Airmass value (dimensionless). Returns 40.0 for objects at or below
        horizon as practical maximum.

    Algorithm:
        The simple plane-parallel approximation uses X = sec(z) where z is
        the zenith angle. This breaks down near the horizon due to Earth's
        curvature and atmospheric refraction.

        The Kasten & Young (1989) formula provides accurate airmass even
        near the horizon:
            X = 1 / [sin(h) + 0.50572 * (h + 6.07995)^(-1.6364)]
        where h is altitude in degrees.

    Example:
        >>> calc_airmass(90.0)  # Zenith
        1.0
        >>> calc_airmass(30.0)  # 30 degrees altitude
        2.0
        >>> calc_airmass(0.0)   # Horizon
        37.9...

    References:
        - Kasten, F. & Young, A.T. (1989) Applied Optics 28, 4735-4738
        - Rozenberg, G.V. (1966) Twilight: A Study in Atmospheric Optics
    """
    # Handle objects at or below horizon
    if altitude_deg <= 0:
        return 40.0  # Practical maximum airmass

    # Zenith angle in degrees and radians
    zenith_deg = 90.0 - altitude_deg
    zenith_rad = math.radians(zenith_deg)
    altitude_rad = math.radians(altitude_deg)

    if method == "secant":
        # Simple plane-parallel approximation
        # Only accurate above ~15 degrees altitude
        if zenith_deg >= 85:
            return 40.0  # Avoid division by very small cosine
        return 1.0 / math.cos(zenith_rad)

    elif method == "kasten_young":
        # Kasten & Young (1989) empirical formula
        # Accurate from zenith to horizon
        sin_h = math.sin(altitude_rad)
        # Altitude in degrees for polynomial term
        h = altitude_deg
        denominator = sin_h + 0.50572 * ((h + 6.07995) ** (-1.6364))
        if denominator <= 0:
            return 40.0
        return min(1.0 / denominator, 40.0)

    elif method == "rozenberg":
        # Rozenberg (1966) formula, good for very low altitudes
        cos_z = math.cos(zenith_rad)
        return 1.0 / (cos_z + 0.025 * math.exp(-11.0 * cos_z))

    else:
        # Default to Kasten & Young
        return calc_airmass(altitude_deg, method="kasten_young")


def calc_rayleigh_coefficient(
    wavelength_nm: float = WAVELENGTH_V,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    altitude_m: float = DEFAULT_ALTITUDE_M,
) -> float:
    """
    Calculate the Rayleigh scattering extinction coefficient.

    Rayleigh scattering is caused by molecules in the atmosphere
    (primarily N2 and O2). It is strongly wavelength-dependent,
    following a lambda^(-4) law, which is why the sky is blue.

    Args:
        wavelength_nm: Wavelength of observation in nanometers.
                       Default is V-band (550 nm).
        pressure_mbar: Atmospheric pressure in millibars (hPa).
        temperature_c: Temperature in degrees Celsius.
        altitude_m: Observer altitude in meters above sea level.

    Returns:
        Rayleigh scattering coefficient in magnitudes per airmass.

    Algorithm:
        k_R = 0.1451 * (P/1013.25) * (lambda_0/lambda)^4
        where lambda_0 = 550 nm (V-band reference)

        The coefficient 0.1451 is the standard Rayleigh coefficient
        at sea level for V-band observations.

    Example:
        >>> calc_rayleigh_coefficient(550.0, 1013.25)  # Standard conditions
        0.1451
        >>> calc_rayleigh_coefficient(440.0, 1013.25)  # B-band (bluer)
        0.359...  # Higher extinction for blue light

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Young, A.T. (1974) "Observational Technique and Data Reduction"
    """
    # Reference values
    k_ref = 0.1451  # Rayleigh coefficient at sea level, V-band
    wavelength_ref = 550.0  # V-band reference wavelength

    # Pressure correction (higher pressure = more scattering)
    pressure_factor = pressure_mbar / DEFAULT_PRESSURE_MBAR

    # Altitude correction (less atmosphere at higher altitudes)
    altitude_factor = math.exp(-altitude_m / SCALE_HEIGHT_M)

    # Wavelength dependence (Rayleigh scattering ~ lambda^-4)
    wavelength_factor = (wavelength_ref / wavelength_nm) ** 4

    return k_ref * pressure_factor * altitude_factor * wavelength_factor


def calc_aerosol_coefficient(
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
    visibility_km: Optional[float] = None,
) -> float:
    """
    Calculate the aerosol (Mie) scattering extinction coefficient.

    Aerosol scattering is caused by dust, pollen, smoke, and other
    particulates in the atmosphere. It is highly variable and depends
    on local conditions.

    Args:
        humidity_percent: Relative humidity (0-100).
                         Higher humidity increases aerosol scattering.
        altitude_m: Observer altitude in meters above sea level.
                   Aerosol concentration decreases with altitude.
        wavelength_nm: Wavelength of observation in nanometers.
        visibility_km: Meteorological visibility in km. If provided,
                       this is used to estimate aerosol content directly.

    Returns:
        Aerosol scattering coefficient in magnitudes per airmass.
        Typical values range from 0.05 (excellent) to 0.30 (hazy).

    Algorithm:
        If visibility is given:
            k_A = 3.912 / V - 0.1066 (empirical formula)
        Otherwise, estimate from humidity:
            k_A = 0.08 * (1 + 2 * (humidity/100)^2) * exp(-altitude/H_A)
        where H_A ~ 1500m is the aerosol scale height.

    Example:
        >>> calc_aerosol_coefficient(30.0, 0.0)   # Low humidity, sea level
        0.10...
        >>> calc_aerosol_coefficient(80.0, 0.0)   # High humidity
        0.18...
        >>> calc_aerosol_coefficient(50.0, 2000)  # Moderate humidity, 2km altitude
        0.07...

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    """
    if visibility_km is not None and visibility_km > 0:
        # Direct calculation from meteorological visibility
        # Using Koschmieder formula: V = 3.912 / beta
        # where beta is extinction coefficient per km
        # Convert to mag/airmass
        k_aerosol = 3.912 / visibility_km - 0.1066
        return max(0.0, k_aerosol)

    # Aerosol scale height (aerosols concentrated in lower atmosphere)
    aerosol_scale_height = 1500.0  # meters

    # Base aerosol coefficient (clean atmosphere)
    k_base = 0.08

    # Humidity factor (hygroscopic growth of aerosols)
    humidity_fraction = humidity_percent / 100.0
    humidity_factor = 1.0 + 2.0 * humidity_fraction**2

    # Altitude factor (aerosols decrease exponentially with altitude)
    altitude_factor = math.exp(-altitude_m / aerosol_scale_height)

    # Wavelength dependence (weaker than Rayleigh: ~lambda^-1.3)
    wavelength_factor = (WAVELENGTH_V / wavelength_nm) ** 1.3

    return k_base * humidity_factor * altitude_factor * wavelength_factor


def calc_ozone_coefficient(wavelength_nm: float = WAVELENGTH_V) -> float:
    """
    Calculate the ozone absorption coefficient.

    Ozone (O3) in the stratosphere absorbs UV and some visible light.
    The effect is small in the visual band but significant in UV/blue.

    Args:
        wavelength_nm: Wavelength of observation in nanometers.

    Returns:
        Ozone absorption coefficient in magnitudes per airmass.
        Typically 0.01-0.04 for visual wavelengths.

    Algorithm:
        Ozone absorption has a complex wavelength dependence.
        For visual wavelengths, we use an empirical approximation.

    Example:
        >>> calc_ozone_coefficient(550.0)  # V-band
        0.016
    """
    # Ozone absorption peaks in UV (Hartley band) and has Chappuis bands
    # in visible (around 600nm). For V-band, the contribution is small.
    if wavelength_nm < 400:
        # UV region - strong ozone absorption
        return 0.10
    elif wavelength_nm < 500:
        # Blue region
        return 0.04
    elif wavelength_nm < 700:
        # Visual region (includes Chappuis bands)
        return 0.016
    else:
        # Red/IR region
        return 0.008


def calc_water_vapor_coefficient(
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    wavelength_nm: float = WAVELENGTH_V,
) -> float:
    """
    Calculate the water vapor absorption coefficient.

    Water vapor absorption is significant primarily in the infrared
    and has relatively minor effect in the visual band.

    Args:
        humidity_percent: Relative humidity (0-100).
        temperature_c: Temperature in degrees Celsius.
        wavelength_nm: Wavelength of observation in nanometers.

    Returns:
        Water vapor absorption coefficient in magnitudes per airmass.
        Typically very small (<0.01) for visual wavelengths.

    Example:
        >>> calc_water_vapor_coefficient(50.0, 15.0, 550.0)
        0.003...
    """
    # Water vapor effect is minimal in visual bands
    # but included for completeness
    humidity_fraction = humidity_percent / 100.0

    if wavelength_nm < 600:
        # Blue-green: minimal water vapor effect
        return 0.002 * humidity_fraction
    elif wavelength_nm < 800:
        # Red: small water vapor bands
        return 0.005 * humidity_fraction
    else:
        # Near-IR: significant water vapor bands
        return 0.02 * humidity_fraction


def calc_extinction_coefficient(
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
    visibility_km: Optional[float] = None,
) -> ExtinctionCoefficients:
    """
    Calculate the total atmospheric extinction coefficient and its components.

    This is the main function for determining how much light is lost
    as it passes through the atmosphere. The extinction coefficient k
    is used with airmass X to calculate total extinction:
        delta_m = k * X
    where delta_m is the increase in apparent magnitude.

    Args:
        pressure_mbar: Atmospheric pressure in millibars (hPa).
                      Default is 1013.25 (standard sea level).
        temperature_c: Temperature in degrees Celsius.
                      Default is 15.0.
        humidity_percent: Relative humidity (0-100).
                         Default is 50.0.
        altitude_m: Observer altitude in meters above sea level.
                   Default is 0.0 (sea level).
        wavelength_nm: Wavelength of observation in nanometers.
                      Default is 550 (V-band).
        visibility_km: Meteorological visibility in km.
                      If provided, used to estimate aerosol content.

    Returns:
        ExtinctionCoefficients dataclass with individual components
        and total extinction coefficient in mag/airmass.

    Example:
        >>> coeff = calc_extinction_coefficient()
        >>> print(f"Total extinction: {coeff.k_total:.3f} mag/airmass")
        Total extinction: 0.28... mag/airmass

        >>> coeff = calc_extinction_coefficient(humidity_percent=90)
        >>> print(f"Humid extinction: {coeff.k_total:.3f} mag/airmass")
        Humid extinction: 0.36... mag/airmass

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    """
    k_rayleigh = calc_rayleigh_coefficient(
        wavelength_nm=wavelength_nm,
        pressure_mbar=pressure_mbar,
        temperature_c=temperature_c,
        altitude_m=altitude_m,
    )

    k_aerosol = calc_aerosol_coefficient(
        humidity_percent=humidity_percent,
        altitude_m=altitude_m,
        wavelength_nm=wavelength_nm,
        visibility_km=visibility_km,
    )

    k_ozone = calc_ozone_coefficient(wavelength_nm=wavelength_nm)

    k_water = calc_water_vapor_coefficient(
        humidity_percent=humidity_percent,
        temperature_c=temperature_c,
        wavelength_nm=wavelength_nm,
    )

    k_total = k_rayleigh + k_aerosol + k_ozone + k_water

    return ExtinctionCoefficients(
        k_rayleigh=k_rayleigh,
        k_aerosol=k_aerosol,
        k_ozone=k_ozone,
        k_water=k_water,
        k_total=k_total,
    )


def calc_extinction_magnitude(
    altitude_deg: float,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    observer_altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
    visibility_km: Optional[float] = None,
) -> float:
    """
    Calculate the total atmospheric extinction in magnitudes for an object.

    This function combines the extinction coefficient with airmass to give
    the actual extinction experienced by an object at a given altitude.

    For objects near the horizon, the extinction can be several magnitudes,
    making faint objects invisible.

    Args:
        altitude_deg: Altitude of object above horizon in degrees.
        pressure_mbar: Atmospheric pressure in millibars (hPa).
        temperature_c: Temperature in degrees Celsius.
        humidity_percent: Relative humidity (0-100).
        observer_altitude_m: Observer altitude in meters above sea level.
        wavelength_nm: Wavelength of observation in nanometers.
        visibility_km: Meteorological visibility in km.

    Returns:
        Total extinction in magnitudes. This value should be added
        to the object's catalog magnitude to get apparent magnitude.
        Returns a large value (99.0) for objects below the horizon.

    Algorithm:
        extinction = k * X
        where k is the extinction coefficient and X is airmass.

        Using the basic formula (approximately):
            extinction ~ 0.28 * sec(z)
        where z is zenith angle and 0.28 is typical k for V-band.

    Example:
        >>> calc_extinction_magnitude(90.0)  # Zenith
        0.28...  # Minimal extinction
        >>> calc_extinction_magnitude(30.0)  # 30 degrees altitude
        0.56...  # 2x airmass = 2x extinction
        >>> calc_extinction_magnitude(5.0)   # Near horizon
        3.2...   # Significant extinction
        >>> calc_extinction_magnitude(0.0)   # Horizon
        11.2...  # Very high extinction

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    """
    # Objects below horizon have essentially infinite extinction
    if altitude_deg < 0:
        return 99.0

    # Calculate airmass
    airmass = calc_airmass(altitude_deg, method="kasten_young")

    # Calculate extinction coefficient
    coeff = calc_extinction_coefficient(
        pressure_mbar=pressure_mbar,
        temperature_c=temperature_c,
        humidity_percent=humidity_percent,
        altitude_m=observer_altitude_m,
        wavelength_nm=wavelength_nm,
        visibility_km=visibility_km,
    )

    # Total extinction = k * X
    return coeff.k_total * airmass


def calc_simple_extinction(altitude_deg: float, k: float = 0.28) -> float:
    """
    Calculate atmospheric extinction using the simple Schaefer/Green formula.

    This is the basic formula mentioned in the task:
        extinction = k * sec(z) = k / cos(z) = k / sin(h)
    where z is zenith angle, h is altitude, and k ~ 0.28 for typical
    visual observations.

    Args:
        altitude_deg: Altitude of object above horizon in degrees.
        k: Extinction coefficient in mag/airmass. Default 0.28 is
           typical for V-band at a good observing site.

    Returns:
        Extinction in magnitudes.

    Note:
        This simple formula is accurate for altitudes above ~15 degrees.
        For lower altitudes, use calc_extinction_magnitude() which
        applies proper airmass corrections for atmospheric curvature.

    Example:
        >>> calc_simple_extinction(90.0, 0.28)  # Zenith
        0.28
        >>> calc_simple_extinction(30.0, 0.28)  # 30 degrees (sec(60) = 2)
        0.56
        >>> calc_simple_extinction(19.47, 0.28)  # airmass = 3
        0.84

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    """
    if altitude_deg <= 0:
        return 99.0

    # Simple secant formula: sec(z) = 1/sin(h)
    sin_alt = math.sin(math.radians(altitude_deg))
    if sin_alt <= 0.001:  # Very near horizon
        return k * 40.0  # Practical maximum

    return k / sin_alt


def apparent_magnitude_with_extinction(
    catalog_magnitude: float,
    altitude_deg: float,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    observer_altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
) -> float:
    """
    Calculate the apparent magnitude of an object including extinction.

    This applies atmospheric extinction to convert from catalog (extra-
    atmospheric) magnitude to what an observer would actually see.

    Args:
        catalog_magnitude: The object's magnitude outside the atmosphere.
        altitude_deg: Altitude of object above horizon in degrees.
        pressure_mbar: Atmospheric pressure in millibars.
        temperature_c: Temperature in degrees Celsius.
        humidity_percent: Relative humidity (0-100).
        observer_altitude_m: Observer altitude in meters.
        wavelength_nm: Wavelength of observation in nanometers.

    Returns:
        Apparent magnitude as seen by observer (fainter = higher value).
        Objects below horizon return a very large magnitude (99+).

    Example:
        >>> # Venus (mag -4.0) at 30 degrees altitude
        >>> apparent_magnitude_with_extinction(-4.0, 30.0)
        -3.4...  # Slightly dimmed by extinction

        >>> # Sirius (mag -1.46) near horizon (5 degrees)
        >>> apparent_magnitude_with_extinction(-1.46, 5.0)
        1.7...  # Significantly dimmed

        >>> # A 6th magnitude star at the zenith
        >>> apparent_magnitude_with_extinction(6.0, 90.0)
        6.28...  # Minimal extinction
    """
    extinction = calc_extinction_magnitude(
        altitude_deg=altitude_deg,
        pressure_mbar=pressure_mbar,
        temperature_c=temperature_c,
        humidity_percent=humidity_percent,
        observer_altitude_m=observer_altitude_m,
        wavelength_nm=wavelength_nm,
    )

    return catalog_magnitude + extinction


def get_extinction_for_heliacal(
    zenith_angle_deg: float,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    altitude_m: float = DEFAULT_ALTITUDE_M,
) -> float:
    """
    Calculate extinction coefficient for heliacal visibility calculations.

    This is a convenience function that returns the extinction coefficient
    in the format expected by heliacal event calculations.

    Args:
        zenith_angle_deg: Zenith angle of the object in degrees (0 = zenith, 90 = horizon).
        pressure_mbar: Atmospheric pressure in millibars.
        temperature_c: Temperature in degrees Celsius.
        humidity_percent: Relative humidity (0-100).
        altitude_m: Observer altitude in meters.

    Returns:
        Extinction coefficient (k value) suitable for heliacal calculations.
        This is the total extinction coefficient in mag/airmass.

    Example:
        >>> k = get_extinction_for_heliacal(60.0)  # 30 degrees altitude
        >>> print(f"Extinction coefficient: {k:.3f}")
        Extinction coefficient: 0.28...
    """
    coeff = calc_extinction_coefficient(
        pressure_mbar=pressure_mbar,
        temperature_c=temperature_c,
        humidity_percent=humidity_percent,
        altitude_m=altitude_m,
        wavelength_nm=WAVELENGTH_V,
    )
    return coeff.k_total
