"""
Atmospheric refraction via numerical ray-tracing through the ICAO Standard
Atmosphere (ISO 2533:1975).

Instead of relying on empirical curve-fits (Bennett 1982, Saemundsson 1986,
etc.), this module computes refraction from first principles by tracing a
light ray through a layered model atmosphere using Snell's law in spherical
geometry.

Physical model
--------------
The atmosphere is divided into concentric spherical shells.  Within each
shell the temperature profile follows the ICAO Standard Atmosphere:

    Layer             Altitude (km)    Lapse rate (K/km)    T_base (K)
    ─────────────     ──────────────   ─────────────────    ──────────
    Troposphere       0 – 11           -6.5                 288.15
    Tropopause        11 – 20           0.0                 216.65
    Stratosphere 1    20 – 32          +1.0                 216.65
    Stratosphere 2    32 – 47          +2.8                 228.65
    Upper             47 – 51           0.0                 270.65
    Mesosphere 1      51 – 71          -2.8                 270.65
    Mesosphere 2      71 – 84.852      -2.0                 214.65

Pressure at each altitude follows from the hydrostatic equation
integrated exactly for each layer type (isothermal or linear lapse).

The refractive index of dry air is computed from the gas-law
approximation (Barrell & Sears 1939):

    n(P, T) = 1 + 7.86e-4 * P / T

where P is in mbar and T in Kelvin.  This is accurate to ~0.5 ppm in
the visible band for dry air, sufficient for refraction work.

At each layer boundary, Snell's law in spherical geometry is applied:

    n_1 * r_1 * sin(z_1)  =  n_2 * r_2 * sin(z_2)

where r = R_earth + h is the geocentric radius and z is the local
zenith angle.  The total refraction is the accumulated angular
deviation of the ray.

For the SE_APP_TO_TRUE direction the function inverts the
SE_TRUE_TO_APP computation numerically via Newton-Raphson iteration
(typically 4-5 steps to full float64 convergence).

References
----------
- ISO 2533:1975 "Standard Atmosphere"
- ICAO Doc 7488/3 "Manual of the ICAO Standard Atmosphere"
- Barrell, H. & Sears, J.E. (1939), Phil. Trans. Roy. Soc. A, 238, 1
- Smart, W.M. (1977) "Textbook on Spherical Astronomy", Ch. VI
- Green, R.M. (1985) "Spherical Astronomy", Cambridge Univ. Press, Ch. 4
- Bomford, G. (1980) "Geodesy", 4th ed., Clarendon Press, §2.17-2.20
"""

from __future__ import annotations

import math
from typing import Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_R_EARTH: float = 6_371_000.0  # Mean Earth radius [m]
_G0: float = 9.80665  # Standard gravity [m/s^2]
_M_AIR: float = 0.028_964_4  # Molar mass of dry air [kg/mol]
_R_GAS: float = 8.314_462_618  # Universal gas constant [J/(mol*K)]

# ICAO atmosphere layers: (base_altitude_m, lapse_rate_K_per_m, T_base_K)
# Lapse rate is dT/dh (negative = temperature decreases with altitude).
_ICAO_LAYERS: list = [
    (0.0, -0.0065, 288.15),  # Troposphere
    (11_000.0, 0.0, 216.65),  # Tropopause (isothermal)
    (20_000.0, 0.001, 216.65),  # Stratosphere 1
    (32_000.0, 0.0028, 228.65),  # Stratosphere 2
    (47_000.0, 0.0, 270.65),  # Stratopause (isothermal)
    (51_000.0, -0.0028, 270.65),  # Mesosphere 1
    (71_000.0, -0.002, 214.65),  # Mesosphere 2
]
_ATMO_TOP: float = 84_852.0  # Top of modeled atmosphere [m]

# Standard sea-level conditions
_P0: float = 1013.25  # Standard pressure [mbar]
_T0: float = 288.15  # Standard temperature [K]

# Number of integration shells (adaptive: denser near ground)
_N_LAYERS_TROPO: int = 400  # 0 – 11 km (where refraction matters)
_N_LAYERS_UPPER: int = 60  # 11 – 85 km


# ---------------------------------------------------------------------------
# ICAO Standard Atmosphere: temperature & pressure profiles
# ---------------------------------------------------------------------------


def _icao_temperature(h: float) -> float:
    """Temperature at altitude *h* [m] according to the ICAO Standard
    Atmosphere.

    Returns temperature in Kelvin.
    """
    for i in range(len(_ICAO_LAYERS) - 1, -1, -1):
        h_base, lapse, T_base = _ICAO_LAYERS[i]
        if h >= h_base:
            return T_base + lapse * (h - h_base)
    return _ICAO_LAYERS[0][2]


def _icao_pressure(h: float) -> float:
    """Pressure at altitude *h* [m] according to the ICAO Standard
    Atmosphere.

    Uses the exact hydrostatic integral:
      - For isothermal layers:  P = P_base * exp(-g*M*(h-h_base) / (R*T))
      - For linear-lapse layers: P = P_base * (T/T_base)^(g*M / (R*lapse))

    Returns pressure in mbar (hPa).
    """
    # Walk up through layers, accumulating pressure from sea level.
    p = _P0
    exponent = _G0 * _M_AIR / _R_GAS  # g*M/R ≈ 0.03416 K/m

    for i, (h_base, lapse, T_base) in enumerate(_ICAO_LAYERS):
        # Top of this layer (= base of next, or atmosphere top)
        if i + 1 < len(_ICAO_LAYERS):
            h_top = _ICAO_LAYERS[i + 1][0]
        else:
            h_top = _ATMO_TOP

        if h <= h_top:
            # Target altitude is within this layer.
            dh = h - h_base
            if abs(lapse) < 1e-10:
                # Isothermal: P = P_base * exp(-exponent * dh / T_base)
                p *= math.exp(-exponent * dh / T_base)
            else:
                # Linear lapse: P = P_base * (T/T_base)^(exponent/lapse)
                T_at_h = T_base + lapse * dh
                p *= (T_at_h / T_base) ** (-exponent / lapse)
            return p

        # Move pressure to the top of this layer.
        dh = h_top - h_base
        if abs(lapse) < 1e-10:
            p *= math.exp(-exponent * dh / T_base)
        else:
            T_top = T_base + lapse * dh
            p *= (T_top / T_base) ** (-exponent / lapse)

    return p


def _refractive_index(p_mbar: float, T_kelvin: float) -> float:
    """Refractive index of dry air at pressure *p_mbar* and temperature
    *T_kelvin*.

    Derived from the Barrell & Sears (1939) group-refractive-index at
    standard conditions (n_s - 1 = 2.93e-4 at 0 degC, 1013.25 mbar)
    scaled to arbitrary P and T via the ideal gas law:

        n - 1 = (n_s - 1) * (P / P_s) * (T_s / T)
              = 2.93e-4 * P / 1013.25 * 273.15 / T
              = 7.905e-5 * P / T

    At sea level under ICAO standard conditions (1013.25 mbar, 288.15 K)
    this gives (n - 1) = 2.78e-4, in agreement with the Edlen (1966)
    value for the visible band to within 1%.  This is also consistent
    with the Smith & Weintraub (1953) radio/optical approximation
    (coefficient 7.76e-5).
    """
    if T_kelvin <= 0 or p_mbar <= 0:
        return 1.0
    return 1.0 + 7.905e-5 * p_mbar / T_kelvin


# ---------------------------------------------------------------------------
# Modified atmosphere: observer-supplied surface conditions
# ---------------------------------------------------------------------------


def _modified_temperature(
    h: float, obs_alt: float, surface_T: float, lapse_rate: float
) -> float:
    """Temperature at altitude *h*, using the observer's surface conditions
    in the troposphere and the ICAO profile above.

    For h < 11 km the temperature follows:
        T(h) = surface_T - lapse_rate * (h - obs_alt)

    Above 11 km the ICAO profile is used, anchored to whatever T exists
    at 11 km from the modified troposphere.
    """
    tropo_top = 11_000.0
    if h < tropo_top:
        return surface_T - lapse_rate * (h - obs_alt)
    # Above troposphere: use ICAO profile offset so that T is continuous
    # at the tropopause.
    T_at_tropo = surface_T - lapse_rate * (tropo_top - obs_alt)
    # Walk through upper ICAO layers
    for i in range(1, len(_ICAO_LAYERS)):
        h_base = _ICAO_LAYERS[i][0]
        lapse_i = _ICAO_LAYERS[i][1]
        if i + 1 < len(_ICAO_LAYERS):
            h_top = _ICAO_LAYERS[i + 1][0]
        else:
            h_top = _ATMO_TOP
        if h <= h_top:
            if i == 1:
                # Tropopause: use the anchored temperature
                return T_at_tropo + lapse_i * (h - h_base)
            else:
                # Higher layers: shift by the tropopause offset
                offset = T_at_tropo - _ICAO_LAYERS[1][2]
                T_base_shifted = _ICAO_LAYERS[i][2] + offset
                return T_base_shifted + lapse_i * (h - h_base)
    return 200.0  # Fallback for extreme altitudes


def _modified_pressure(
    h: float, obs_alt: float, obs_P: float, surface_T: float, lapse_rate: float
) -> float:
    """Pressure at altitude *h* using the modified atmosphere.

    Uses the exact hydrostatic integral from the observer's conditions.
    For a linear temperature lapse (troposphere):

        P(h) = P_obs * (T(h) / T_obs) ^ (g*M / (R*gamma))

    For isothermal layers:

        P(h) = P_base * exp(-g*M*(h - h_base) / (R*T))
    """
    exponent = _G0 * _M_AIR / _R_GAS  # g*M/R ≈ 0.03416 K/m
    tropo_top = 11_000.0

    if h <= tropo_top:
        # Troposphere: linear lapse rate from observer conditions.
        T_obs = surface_T
        T_h = T_obs - lapse_rate * (h - obs_alt)
        if abs(lapse_rate) < 1e-10:
            return obs_P * math.exp(-exponent * (h - obs_alt) / T_obs)
        return obs_P * (T_h / T_obs) ** (exponent / lapse_rate)

    # First, compute pressure at the tropopause using tropospheric formula.
    T_obs = surface_T
    T_tropo = T_obs - lapse_rate * (tropo_top - obs_alt)
    if abs(lapse_rate) < 1e-10:
        P_tropo = obs_P * math.exp(-exponent * (tropo_top - obs_alt) / T_obs)
    else:
        P_tropo = obs_P * (T_tropo / T_obs) ** (exponent / lapse_rate)

    # Above troposphere: use ICAO layer structure (isothermal or linear).
    p = P_tropo
    T_at_tropo = T_tropo  # Anchored temperature at 11 km

    for i in range(1, len(_ICAO_LAYERS)):
        h_base = _ICAO_LAYERS[i][0]
        lapse_i = _ICAO_LAYERS[i][1]
        if i + 1 < len(_ICAO_LAYERS):
            h_top_i = _ICAO_LAYERS[i + 1][0]
        else:
            h_top_i = _ATMO_TOP

        # Shift T_base to match the anchored tropopause temperature
        if i == 1:
            T_base_i = T_at_tropo
        else:
            offset = T_at_tropo - _ICAO_LAYERS[1][2]
            T_base_i = _ICAO_LAYERS[i][2] + offset

        if h <= h_top_i:
            dh = h - h_base
            if abs(lapse_i) < 1e-10:
                p *= math.exp(-exponent * dh / T_base_i)
            else:
                T_h = T_base_i + lapse_i * dh
                p *= (T_h / T_base_i) ** (-exponent / lapse_i)
            return p

        # Move pressure to top of this layer
        dh = h_top_i - h_base
        if abs(lapse_i) < 1e-10:
            p *= math.exp(-exponent * dh / T_base_i)
        else:
            T_top_i = T_base_i + lapse_i * dh
            p *= (T_top_i / T_base_i) ** (-exponent / lapse_i)

    return p


# ---------------------------------------------------------------------------
# Precomputed shell data cache
# ---------------------------------------------------------------------------

# Cache key: (obs_alt, obs_pressure, obs_T_K, lapse_rate)
# Cache value: list of (r_i, n_i) tuples
_shell_cache_key: tuple = ()
_shell_cache_data: list = []


def _get_shell_data(
    obs_alt: float, obs_pressure: float, obs_T_K: float, lapse_rate: float
) -> list:
    """Return a list of (r_i, n_i) tuples for all shell boundaries.

    Results are cached and reused when called with the same parameters
    (which is the common case: repeated calls from refrac/azalt with
    the same atmospheric conditions).
    """
    global _shell_cache_key, _shell_cache_data

    key = (obs_alt, obs_pressure, obs_T_K, lapse_rate)
    if key == _shell_cache_key and _shell_cache_data:
        return _shell_cache_data

    shells = _build_shell_altitudes(obs_alt)
    data = []
    for h in shells:
        r_i = _R_EARTH + h
        T_i = _modified_temperature(h, obs_alt, obs_T_K, lapse_rate)
        P_i = _modified_pressure(h, obs_alt, obs_pressure, obs_T_K, lapse_rate)
        n_i = _refractive_index(P_i, T_i)
        data.append((r_i, n_i))

    _shell_cache_key = key
    _shell_cache_data = data
    return data


# ---------------------------------------------------------------------------
# Core ray-tracing engine
# ---------------------------------------------------------------------------


def _build_shell_altitudes(obs_alt: float) -> list:
    """Build a list of shell boundary altitudes [m], denser near the
    ground where refraction gradients are steepest.

    Returns a sorted list of altitudes from *obs_alt* to the atmosphere
    top.
    """
    altitudes = []

    # Troposphere: dense sampling (every ~140 m)
    tropo_top = min(11_000.0, _ATMO_TOP)
    start = max(obs_alt, 0.0)
    if start < tropo_top:
        n = _N_LAYERS_TROPO
        dh = (tropo_top - start) / n
        for i in range(n + 1):
            altitudes.append(start + i * dh)

    # Upper atmosphere: coarser sampling (every ~1850 m)
    if tropo_top < _ATMO_TOP:
        n = _N_LAYERS_UPPER
        dh = (_ATMO_TOP - tropo_top) / n
        for i in range(1, n + 1):
            altitudes.append(tropo_top + i * dh)

    # Remove duplicates and sort
    altitudes = sorted(set(altitudes))

    # Ensure observer altitude is the first entry
    if altitudes and altitudes[0] > obs_alt:
        altitudes.insert(0, obs_alt)

    return altitudes


def _bending_at_observed_zenith(
    observed_zenith_deg: float,
    shell_data: list,
) -> float:
    """Total angular bending for a ray observed at *observed_zenith_deg*.

    The Bouguer (Snell) invariant C = n * r * sin(z) is conserved along
    the ray.  At each shell boundary the refractive index changes
    discretely, and the bending is:

        delta_i = z_above_i - z_below_i

    where *z_below* and *z_above* are the zenith angles in the medium
    below and above the boundary, both evaluated at the boundary
    radius *r_i*.  Going upward (dense -> less dense), delta_i > 0.

    The sum of all delta_i gives the total refraction R.

    Parameters
    ----------
    observed_zenith_deg : float
        Observed (apparent) zenith angle, in degrees.
    shell_data : list
        Precomputed list of (r_i, n_i) tuples from _get_shell_data().

    Returns refraction in degrees (>= 0).
    """
    z_rad = math.radians(observed_zenith_deg)

    # Observer is the first shell entry
    r_obs, n_obs = shell_data[0]

    # Bouguer invariant
    C = n_obs * r_obs * math.sin(z_rad)

    total_bending = 0.0
    n_prev = n_obs

    for i in range(1, len(shell_data)):
        r_i, n_i = shell_data[i]

        # Zenith angle just below the boundary (medium n_prev)
        sin_z_below = C / (n_prev * r_i)
        # Zenith angle just above the boundary (medium n_i)
        sin_z_above = C / (n_i * r_i)

        # Clamp to avoid domain errors at extreme zenith angles
        if sin_z_below > 1.0 or sin_z_above > 1.0:
            # Total internal reflection zone; stop tracing
            break

        z_below = math.asin(sin_z_below)
        z_above = math.asin(sin_z_above)

        # Bending at this interface (positive going upward)
        total_bending += z_above - z_below

        n_prev = n_i

    return max(0.0, math.degrees(total_bending))


def _trace_ray(
    true_zenith_deg: float,
    obs_alt: float = 0.0,
    obs_pressure: float = _P0,
    obs_temperature_C: float = 15.0,
    lapse_rate: float = 0.0065,
) -> float:
    """Compute refraction for a source at true zenith angle
    *true_zenith_deg*.

    The refraction R satisfies:

        z_true = z_apparent + R

    where z_apparent is the *observed* zenith angle and z_true is the
    geometric one.  Since the Bouguer invariant depends on z_apparent,
    we iterate:

        1. Start with z_obs = z_true  (R = 0 initial guess)
        2. Compute R = bending(z_obs)
        3. Update z_obs = z_true - R
        4. Repeat until convergence (typically 3-4 steps)

    Returns refraction in degrees (>= 0).
    """
    if true_zenith_deg <= 0.0:
        return 0.0
    if obs_pressure <= 0.0:
        return 0.0

    obs_T_K = obs_temperature_C + 273.15

    # Linear extrapolation for objects below the geometric horizon
    if true_zenith_deg > 91.0:
        r_91 = _trace_ray(91.0, obs_alt, obs_pressure, obs_temperature_C, lapse_rate)
        r_90 = _trace_ray(90.0, obs_alt, obs_pressure, obs_temperature_C, lapse_rate)
        slope = r_91 - r_90
        return max(0.0, r_91 + slope * (true_zenith_deg - 91.0))

    # Precompute shell (r, n) data (cached for repeated calls)
    shell_data = _get_shell_data(obs_alt, obs_pressure, obs_T_K, lapse_rate)
    if len(shell_data) < 2:
        return 0.0

    # Iterate: z_obs = z_true - R(z_obs)
    z_obs = true_zenith_deg
    R = 0.0
    for _ in range(6):
        R = _bending_at_observed_zenith(z_obs, shell_data)
        z_obs_new = true_zenith_deg - R
        if abs(z_obs_new - z_obs) < 1e-10:
            break
        z_obs = z_obs_new

    return max(0.0, R)


# ---------------------------------------------------------------------------
# Public functions (used internally by utils.py)
# ---------------------------------------------------------------------------


def calc_refraction_true_to_app(
    true_alt: float,
    pressure: float = _P0,
    temperature_C: float = 15.0,
    obs_alt: float = 0.0,
    lapse_rate: float = 0.0065,
) -> float:
    """Compute the refraction for a given true (geometric) altitude.

    Parameters
    ----------
    true_alt : float
        True altitude above the geometric horizon, in degrees.
    pressure : float
        Atmospheric pressure at the observer, in mbar.
    temperature_C : float
        Temperature at the observer, in degrees Celsius.
    obs_alt : float
        Observer altitude above sea level, in metres.
    lapse_rate : float
        Tropospheric lapse rate in K/m.

    Returns
    -------
    float
        Refraction in degrees (>= 0).
    """
    if pressure <= 0:
        return 0.0
    zenith = 90.0 - true_alt
    return _trace_ray(zenith, obs_alt, pressure, temperature_C, lapse_rate)


def calc_refraction_app_to_true(
    apparent_alt: float,
    pressure: float = _P0,
    temperature_C: float = 15.0,
    obs_alt: float = 0.0,
    lapse_rate: float = 0.0065,
) -> float:
    """Compute the refraction for a given apparent (observed) altitude.

    Uses Newton-Raphson iteration to invert the true-to-apparent
    mapping.

    Parameters
    ----------
    apparent_alt : float
        Apparent (observed) altitude, in degrees.
    pressure : float
        Atmospheric pressure at the observer, in mbar.
    temperature_C : float
        Temperature at the observer, in degrees Celsius.
    obs_alt : float
        Observer altitude above sea level, in metres.
    lapse_rate : float
        Tropospheric lapse rate in K/m.

    Returns
    -------
    float
        Refraction in degrees (>= 0).
    """
    if pressure <= 0:
        return 0.0

    # Initial estimate: true_alt ≈ apparent_alt
    true_est = apparent_alt

    for _ in range(8):
        r = calc_refraction_true_to_app(
            true_est, pressure, temperature_C, obs_alt, lapse_rate
        )
        # apparent = true + refraction  =>  residual = (true_est + r) - apparent
        residual = true_est + r - apparent_alt
        if abs(residual) < 1e-12:
            break

        # Numerical derivative dr/d(true_alt) via central difference
        dt = 0.001
        r_plus = calc_refraction_true_to_app(
            true_est + dt, pressure, temperature_C, obs_alt, lapse_rate
        )
        # d(apparent)/d(true) = 1 + dr/dtrue
        dapp_dtrue = 1.0 + (r_plus - r) / dt
        if abs(dapp_dtrue) < 1e-15:
            break
        true_est -= residual / dapp_dtrue

    # Refraction = apparent - true
    return max(0.0, apparent_alt - true_est)


def calc_dip(obs_alt: float, lapse_rate: float = 0.0065) -> float:
    """Dip of the horizon for an elevated observer.

    The geometric dip is reduced by terrestrial refraction whose
    strength depends on the temperature lapse rate.

    Parameters
    ----------
    obs_alt : float
        Observer altitude above sea level in metres.
    lapse_rate : float
        Tropospheric lapse rate in K/m.

    Returns
    -------
    float
        Dip angle in degrees (negative: horizon is below the
        geometric horizontal).

    References
    ----------
    Bomford, "Geodesy" (1980), 4th ed., §2.17-2.20
    Torge, "Geodesy" (2001), 3rd ed., §5.1.1
    """
    if obs_alt <= 0:
        return 0.0

    # Geometric dip: arccos(R / (R + h))
    ratio = _R_EARTH / (_R_EARTH + obs_alt)
    if ratio >= 1.0:
        return 0.0
    dip_geometric = math.degrees(math.acos(ratio))

    # Terrestrial refraction coefficient k as a function of lapse rate.
    # k = 0.1117 + 3.5516 * gamma  (linear model from geodetic literature).
    # At standard lapse rate (0.0065 K/m) this gives k ≈ 0.1348,
    # consistent with the geodetic standard k ≈ 1/7.4.
    if lapse_rate > 0:
        k = 0.1117 + 3.5516 * lapse_rate
    else:
        k = 0.0

    # Observed dip = geometric dip * (1 - k)
    return -dip_geometric * (1.0 - k)
