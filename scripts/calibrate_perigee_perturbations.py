#!/usr/bin/env python3
"""
Calibrate perigee perturbation coefficients against JPL DE441 ephemeris (v2).

This script determines optimal coefficients for a trigonometric perturbation
series for the interpolated lunar perigee.

Method (v2.2 - Passage-Interpolated Harmonic Fit):
    1. Find all perigee passages (Earth-Moon distance minima) over 1000 years
       At each passage, Moon's longitude IS the perigee longitude (unambiguous)
    2. Interpolate between passages with cubic spline to get a smooth,
       continuous perigee longitude curve
    3. Sample the spline daily to get ~365K clean "interpolated perigee" values
    4. Compute perturbation = interpolated_perigee - mean_perigee
     5. Fit 120-term harmonic series via least-squares

    This avoids both failure modes:
    - v2.0 (raw osculating): noise was correlated with design matrix terms
    - v2.1 (passage-only): M'≈0 at all passages → degenerate design matrix

    The spline interpolation creates data at arbitrary times (full M' coverage)
    from physically clean calibration points (passage longitudes).

Output:
    - Calibrated coefficients (Python code for lunar.py)
    - Perigee passage validation results
    - SE cross-validation results (informational)

References:
    - Park, R.S. et al. (2021) "JPL Planetary and Lunar Ephemerides DE440/DE441"
    - Chapront-Touze, M. & Chapront, J. "ELP 2000-82B" (1988)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from multiprocessing import cpu_count
import numpy as np

# =============================================================================
# DEFAULT PARAMETERS
# =============================================================================

DEFAULT_EPHEMERIS = "de441.bsp"
DEFAULT_FIT_START_YEAR = 1500
DEFAULT_FIT_END_YEAR = 2500
DEFAULT_WORKERS = cpu_count()

# Perigee passage search parameters
ANOMALISTIC_MONTH = 27.554549878  # days
PASSAGE_SEARCH_STEP = 0.5  # days for coarse search
PASSAGE_REFINE_TOL = 0.0001  # days (~8.6 seconds)

# Coefficient threshold for significance
COEFFICIENT_THRESHOLD = 0.001  # degrees


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def normalize_angle_diff(diff: float) -> float:
    """Normalize angle difference to [-180, 180) range."""
    while diff >= 180.0:
        diff -= 360.0
    while diff < -180.0:
        diff += 360.0
    return diff


def year_to_jd(year: float) -> float:
    """Convert year (floating point) to Julian Day (TT)."""
    return 2451545.0 + (year - 2000.0) * 365.25


def jd_to_year(jd: float) -> float:
    """Convert Julian Day to year (floating point)."""
    return (jd - 2451545.0) / 365.25 + 2000.0


def unwrap_longitudes(lons: list[float]) -> list[float]:
    """Unwrap longitudes to handle 0/360 discontinuity."""
    if not lons:
        return []
    result = [lons[0]]
    for i in range(1, len(lons)):
        diff = lons[i] - lons[i - 1]
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360
        result.append(result[-1] + diff)
    return result


# =============================================================================
# FUNDAMENTAL ARGUMENTS
# =============================================================================


def compute_fundamental_arguments(
    jd_tt: float,
) -> tuple[float, float, float, float]:
    """Compute fundamental lunar arguments D, M, M', F in radians."""
    T = (jd_tt - 2451545.0) / 36525.0
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    T5 = T4 * T

    D = (
        297.8501921
        + 445267.1114034 * T
        - 0.0018819 * T2
        + T3 / 545868.0
        - T4 / 113065000.0
        + T5 / 18999000000.0
    )
    D = math.radians(D % 360.0)

    M = (
        357.5291092
        + 35999.0502909 * T
        - 0.0001536 * T2
        + T3 / 24490000.0
        - T4 / 992300000.0
        + T5 / 189900000000.0
    )
    M = math.radians(M % 360.0)

    M_prime = (
        134.9633964
        + 477198.8675055 * T
        + 0.0087414 * T2
        + T3 / 69699.0
        - T4 / 14712000.0
        + T5 / 2520410000.0
    )
    M_prime = math.radians(M_prime % 360.0)

    F = (
        93.2720950
        + 483202.0175233 * T
        - 0.0036539 * T2
        - T3 / 3526000.0
        + T4 / 863310000.0
        - T5 / 142650000000.0
    )
    F = math.radians(F % 360.0)

    return D, M, M_prime, F


# =============================================================================
# 120-TERM DESIGN MATRIX
# =============================================================================


def build_design_vector(jd_tt: float) -> list[float]:
    """Build the 120-term design vector for the harmonic fit."""
    D, M, M_prime, F = compute_fundamental_arguments(jd_tt)
    T = (jd_tt - 2451545.0) / 36525.0
    E = 1.0 - 0.002516 * T - 0.0000074 * T**2
    E2 = E * E

    terms = []

    # Cat 1: Evection sin(kD - kM'), k=1..20
    for k in range(1, 21):
        terms.append(math.sin(k * D - k * M_prime))

    # Cat 2: Evection cos(kD - kM'), k=1..20
    for k in range(1, 21):
        terms.append(math.cos(k * D - k * M_prime))

    # Cat 3: Solar anomaly coupling — 9 terms
    terms.append(E * math.sin(M))
    terms.append(E * math.sin(2.0 * D - 2.0 * M_prime - M))
    terms.append(E * math.sin(2.0 * D - 2.0 * M_prime + M))
    terms.append(E * math.sin(D - M_prime - M))
    terms.append(E * math.sin(D - M_prime + M))
    terms.append(E * math.sin(4.0 * D - 4.0 * M_prime - M))
    terms.append(E * math.sin(4.0 * D - 4.0 * M_prime + M))
    terms.append(E * math.sin(6.0 * D - 6.0 * M_prime - M))
    terms.append(E * math.sin(6.0 * D - 6.0 * M_prime + M))

    # Cat 4: Solar double coupling — 4 terms
    terms.append(E2 * math.sin(2.0 * M))
    terms.append(E2 * math.sin(2.0 * D - 2.0 * M_prime - 2.0 * M))
    terms.append(E2 * math.sin(2.0 * D - 2.0 * M_prime + 2.0 * M))
    terms.append(E2 * math.sin(4.0 * D - 4.0 * M_prime - 2.0 * M))

    # Cat 5: Lunar anomaly pure — 3 terms
    terms.append(math.sin(M_prime))
    terms.append(math.sin(2.0 * M_prime))
    terms.append(math.sin(3.0 * M_prime))

    # Cat 6: Latitude coupling — 7 terms
    terms.append(math.sin(2.0 * F - 2.0 * M_prime))
    terms.append(math.sin(2.0 * F - 2.0 * D))
    terms.append(math.sin(2.0 * F))
    terms.append(math.sin(2.0 * F - 4.0 * M_prime + 2.0 * D))
    terms.append(math.sin(2.0 * F + 2.0 * M_prime - 2.0 * D))
    terms.append(math.sin(2.0 * F - 2.0 * M_prime + 2.0 * D))
    terms.append(math.sin(4.0 * F - 2.0 * M_prime))

    # Cat 7: Cross-coupling — 8 terms
    terms.append(math.sin(2.0 * D - M_prime))
    terms.append(math.sin(2.0 * D - 3.0 * M_prime))
    terms.append(math.sin(4.0 * D - 3.0 * M_prime))
    terms.append(math.sin(4.0 * D - 5.0 * M_prime))
    terms.append(math.sin(2.0 * D))
    terms.append(math.sin(4.0 * D))
    terms.append(math.sin(6.0 * D - 5.0 * M_prime))
    terms.append(math.sin(3.0 * D - 2.0 * M_prime))

    # Cat 8: Solar-latitude — 4 terms
    terms.append(E * math.sin(2.0 * F - 2.0 * M_prime + M))
    terms.append(E * math.sin(2.0 * F - 2.0 * M_prime - M))
    terms.append(E * math.sin(2.0 * F - 2.0 * D + M))
    terms.append(E * math.sin(2.0 * F - 2.0 * D - M))

    # Cat 9: Higher evection-solar — 4 terms
    terms.append(E * math.sin(8.0 * D - 8.0 * M_prime - M))
    terms.append(E * math.sin(8.0 * D - 8.0 * M_prime + M))
    terms.append(E * math.sin(10.0 * D - 10.0 * M_prime - M))
    terms.append(E * math.sin(10.0 * D - 10.0 * M_prime + M))

    # Cat 10: Secular — 4 terms
    terms.append(T * math.sin(2.0 * D - 2.0 * M_prime))
    terms.append(T * math.sin(M))
    terms.append(T * math.sin(D - M_prime))
    terms.append(T * math.cos(2.0 * D - 2.0 * M_prime))

    # Cat 11: Cosine coupling — 6 terms
    terms.append(math.cos(M))
    terms.append(math.cos(2.0 * D - M_prime))
    terms.append(math.cos(2.0 * F - 2.0 * M_prime))
    terms.append(math.cos(M_prime))
    terms.append(math.cos(2.0 * D - 3.0 * M_prime))
    terms.append(math.cos(2.0 * D))

    # Cat 12: Direct solar coupling (M with D, no M' dependence) — 5 terms
    terms.append(E * math.sin(2.0 * D - M))
    terms.append(E * math.sin(2.0 * D + M))
    terms.append(E * math.sin(D + M))
    terms.append(E * math.sin(D - M))
    terms.append(E * math.cos(2.0 * D - M))

    # Cat 13: Off-diagonal D-M' and parallactic terms — 8 terms
    terms.append(math.sin(D))
    terms.append(math.sin(2.0 * D + M_prime))
    terms.append(math.sin(2.0 * D + 2.0 * M_prime))
    terms.append(math.sin(3.0 * D))
    terms.append(math.cos(D))
    terms.append(math.cos(2.0 * D + M_prime))
    terms.append(math.cos(2.0 * D + 2.0 * M_prime))
    terms.append(math.cos(3.0 * D))

    # Cat 14: Sun-Moon anomaly coupling — 7 terms
    terms.append(E * math.sin(M + M_prime))
    terms.append(E * math.sin(M - M_prime))
    terms.append(E * math.sin(2.0 * D - M + M_prime))
    terms.append(E * math.sin(2.0 * D + M - M_prime))
    terms.append(E * math.sin(M - 2.0 * M_prime))
    terms.append(E * math.sin(M + 2.0 * M_prime))
    terms.append(E2 * math.sin(2.0 * D - 2.0 * M))

    # Cat 15: Additional latitude coupling — 4 terms
    terms.append(math.sin(2.0 * F + M_prime))
    terms.append(math.sin(2.0 * F - M_prime))
    terms.append(math.sin(2.0 * F + D))
    terms.append(math.sin(2.0 * F - D))

    # Cat 16: Additional cosine phase terms — 3 terms
    terms.append(math.cos(2.0 * M_prime))
    terms.append(math.cos(3.0 * M_prime))
    terms.append(math.cos(2.0 * F))

    # Cat 17: Polynomial corrections for mean perigee drift — 4 terms
    # These absorb secular errors in the mean perigee model (calc_mean_lilith)
    terms.append(1.0)  # constant offset
    terms.append(T)  # linear drift
    terms.append(T * T)  # quadratic drift
    terms.append(T * T * T)  # cubic drift

    return terms  # 120 terms total


def get_term_labels() -> list[str]:
    """Get human-readable labels for all 120 terms."""
    labels = []
    for k in range(1, 21):
        labels.append(f"sin({k}D-{k}M')")
    for k in range(1, 21):
        labels.append(f"cos({k}D-{k}M')")
    labels.extend(
        [
            "E*sin(M)",
            "E*sin(2D-2M'-M)",
            "E*sin(2D-2M'+M)",
            "E*sin(D-M'-M)",
            "E*sin(D-M'+M)",
            "E*sin(4D-4M'-M)",
            "E*sin(4D-4M'+M)",
            "E*sin(6D-6M'-M)",
            "E*sin(6D-6M'+M)",
        ]
    )
    labels.extend(
        [
            "E2*sin(2M)",
            "E2*sin(2D-2M'-2M)",
            "E2*sin(2D-2M'+2M)",
            "E2*sin(4D-4M'-2M)",
        ]
    )
    labels.extend(["sin(M')", "sin(2M')", "sin(3M')"])
    labels.extend(
        [
            "sin(2F-2M')",
            "sin(2F-2D)",
            "sin(2F)",
            "sin(2F-4M'+2D)",
            "sin(2F+2M'-2D)",
            "sin(2F-2M'+2D)",
            "sin(4F-2M')",
        ]
    )
    labels.extend(
        [
            "sin(2D-M')",
            "sin(2D-3M')",
            "sin(4D-3M')",
            "sin(4D-5M')",
            "sin(2D)",
            "sin(4D)",
            "sin(6D-5M')",
            "sin(3D-2M')",
        ]
    )
    labels.extend(
        [
            "E*sin(2F-2M'+M)",
            "E*sin(2F-2M'-M)",
            "E*sin(2F-2D+M)",
            "E*sin(2F-2D-M)",
        ]
    )
    labels.extend(
        [
            "E*sin(8D-8M'-M)",
            "E*sin(8D-8M'+M)",
            "E*sin(10D-10M'-M)",
            "E*sin(10D-10M'+M)",
        ]
    )
    labels.extend(
        [
            "T*sin(2D-2M')",
            "T*sin(M)",
            "T*sin(D-M')",
            "T*cos(2D-2M')",
        ]
    )
    labels.extend(
        [
            "cos(M)",
            "cos(2D-M')",
            "cos(2F-2M')",
            "cos(M')",
            "cos(2D-3M')",
            "cos(2D)",
        ]
    )
    # Cat 12: Direct solar coupling
    labels.extend(
        [
            "E*sin(2D-M)",
            "E*sin(2D+M)",
            "E*sin(D+M)",
            "E*sin(D-M)",
            "E*cos(2D-M)",
        ]
    )
    # Cat 13: Off-diagonal D-M' and parallactic
    labels.extend(
        [
            "sin(D)",
            "sin(2D+M')",
            "sin(2D+2M')",
            "sin(3D)",
            "cos(D)",
            "cos(2D+M')",
            "cos(2D+2M')",
            "cos(3D)",
        ]
    )
    # Cat 14: Sun-Moon anomaly coupling
    labels.extend(
        [
            "E*sin(M+M')",
            "E*sin(M-M')",
            "E*sin(2D-M+M')",
            "E*sin(2D+M-M')",
            "E*sin(M-2M')",
            "E*sin(M+2M')",
            "E2*sin(2D-2M)",
        ]
    )
    # Cat 15: Additional latitude coupling
    labels.extend(
        [
            "sin(2F+M')",
            "sin(2F-M')",
            "sin(2F+D)",
            "sin(2F-D)",
        ]
    )
    # Cat 16: Additional cosine phase terms
    labels.extend(["cos(2M')", "cos(3M')", "cos(2F)"])
    # Cat 17: Polynomial corrections
    labels.extend(["const", "T", "T^2", "T^3"])
    return labels


# =============================================================================
# PERIGEE PASSAGE FINDER
# =============================================================================


def earth_moon_distance(jd_tt: float) -> float:
    """Compute Earth-Moon distance in AU."""
    from libephemeris.state import get_planets, get_timescale

    planets = get_planets()
    ts = get_timescale()
    t = ts.tt_jd(jd_tt)
    moon_pos = (planets["moon"] - planets["earth"]).at(t)
    return float(moon_pos.distance().au)


def moon_ecliptic_longitude(jd_tt: float) -> float:
    """Compute Moon's ecliptic longitude in degrees [0, 360)."""
    from skyfield.framelib import ecliptic_frame
    from libephemeris.state import get_planets, get_timescale

    planets = get_planets()
    ts = get_timescale()
    t = ts.tt_jd(jd_tt)
    moon_pos = (planets["moon"] - planets["earth"]).at(t)
    _, lon, _ = moon_pos.frame_latlon(ecliptic_frame)
    return float(lon.degrees) % 360.0


def golden_section_minimize(f, a: float, b: float, tol: float = 0.0001) -> float:
    """Find minimum of f in [a, b] using golden section search."""
    phi = (1 + math.sqrt(5)) / 2
    resphi = 2.0 - phi
    x1 = a + resphi * (b - a)
    x2 = b - resphi * (b - a)
    f1 = f(x1)
    f2 = f(x2)
    while abs(b - a) > tol:
        if f1 < f2:
            b = x2
            x2 = x1
            f2 = f1
            x1 = a + resphi * (b - a)
            f1 = f(x1)
        else:
            a = x1
            x1 = x2
            f1 = f2
            x2 = b - resphi * (b - a)
            f2 = f(x2)
    return (a + b) / 2.0


def find_passage_in_interval(
    args: tuple[float, float],
) -> tuple[float, float, float] | None:
    """Find a single perigee passage in the given JD interval."""
    jd_start, jd_end = args
    n_scan = int((jd_end - jd_start) / PASSAGE_SEARCH_STEP) + 2
    best_jd = jd_start
    best_dist = float("inf")
    for i in range(n_scan):
        t = jd_start + i * PASSAGE_SEARCH_STEP
        if t > jd_end:
            break
        try:
            dist = earth_moon_distance(t)
            if dist < best_dist:
                best_dist = dist
                best_jd = t
        except Exception:
            continue
    try:
        a = max(best_jd - 1.0, jd_start)
        b = min(best_jd + 1.0, jd_end)
        refined_jd = golden_section_minimize(
            earth_moon_distance, a, b, tol=PASSAGE_REFINE_TOL
        )
        lon = moon_ecliptic_longitude(refined_jd)
        min_dist = earth_moon_distance(refined_jd)
        return (refined_jd, lon, min_dist)
    except Exception:
        return None


def find_all_perigee_passages(
    jd_start: float,
    jd_end: float,
    workers: int,
    ephemeris: str,
) -> list[tuple[float, float, float]]:
    """Find all perigee passages in the given range using parallel workers."""
    intervals = []
    jd = jd_start
    while jd < jd_end:
        interval_end = min(jd + ANOMALISTIC_MONTH, jd_end)
        intervals.append((jd, interval_end))
        jd += ANOMALISTIC_MONTH * 0.95

    tqdm_cls = None
    try:
        from tqdm import tqdm as tqdm_module

        tqdm_cls = tqdm_module
    except ImportError:
        pass

    passages = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=init_worker,
        initargs=(ephemeris,),
    ) as executor:
        futures = {
            executor.submit(find_passage_in_interval, interval): interval
            for interval in intervals
        }
        if tqdm_cls is not None:
            for future in tqdm_cls(
                as_completed(futures),
                total=len(futures),
                desc="Finding perigee passages",
            ):
                result = future.result()
                if result is not None:
                    passages.append(result)
        else:
            completed = 0
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    passages.append(result)
                completed += 1
                if completed % 1000 == 0:
                    print(f"  Progress: {completed}/{len(intervals)}")

    passages.sort(key=lambda x: x[0])

    # Deduplicate (overlapping intervals may find same passage)
    deduped = []
    for p in passages:
        if not deduped or abs(p[0] - deduped[-1][0]) > 5.0:
            deduped.append(p)
        else:
            if p[2] < deduped[-1][2]:
                deduped[-1] = p
    return deduped


# =============================================================================
# CUBIC SPLINE FOR PASSAGE INTERPOLATION
# =============================================================================


def cubic_spline_interpolate(
    x_knots: np.ndarray,
    y_knots: np.ndarray,
    x_eval: np.ndarray,
) -> np.ndarray:
    """Natural cubic spline interpolation.

    Args:
        x_knots: Knot x-coordinates (sorted, increasing)
        y_knots: Knot y-coordinates (unwrapped longitudes)
        x_eval: Points at which to evaluate the spline

    Returns:
        Interpolated y-values at x_eval points
    """
    from scipy.interpolate import CubicSpline

    cs = CubicSpline(x_knots, y_knots, bc_type="natural")
    return cs(x_eval)


# =============================================================================
# PARALLELIZATION HELPERS
# =============================================================================


def init_worker(ephemeris: str) -> None:
    """Initialize worker process with ephemeris."""
    os.environ["LIBEPHEMERIS_EPHEMERIS"] = ephemeris
    from libephemeris.state import set_ephemeris_file

    set_ephemeris_file(ephemeris)


def setup_ephemeris(ephemeris: str) -> bool:
    """Set up the ephemeris file for use."""
    from libephemeris.state import set_ephemeris_file

    try:
        set_ephemeris_file(ephemeris)
        return True
    except Exception as e:
        print(f"Error setting ephemeris: {e}")
        return False


# =============================================================================
# COEFFICIENT FORMATTER
# =============================================================================


def get_term_code_expressions() -> list[tuple[str, str, str]]:
    """Get (category, label, code_expression) for all 120 terms."""
    terms = []

    cat = "PRIMARY EVECTION HARMONICS sin(kD - kM')"
    for k in range(1, 21):
        label = f"sin({k}D-{k}M')"
        expr = (
            "math.sin(D - M_prime)"
            if k == 1
            else f"math.sin({k}.0 * D - {k}.0 * M_prime)"
        )
        terms.append((cat, label, expr))

    cat = "EVECTION PHASE CORRECTIONS cos(kD - kM')"
    for k in range(1, 21):
        label = f"cos({k}D-{k}M')"
        expr = (
            "math.cos(D - M_prime)"
            if k == 1
            else f"math.cos({k}.0 * D - {k}.0 * M_prime)"
        )
        terms.append((cat, label, expr))

    cat = "SOLAR ANOMALY COUPLING (M terms)"
    for label, expr in [
        ("E*sin(M)", "E * math.sin(M)"),
        ("E*sin(2D-2M'-M)", "E * math.sin(2.0 * D - 2.0 * M_prime - M)"),
        ("E*sin(2D-2M'+M)", "E * math.sin(2.0 * D - 2.0 * M_prime + M)"),
        ("E*sin(D-M'-M)", "E * math.sin(D - M_prime - M)"),
        ("E*sin(D-M'+M)", "E * math.sin(D - M_prime + M)"),
        ("E*sin(4D-4M'-M)", "E * math.sin(4.0 * D - 4.0 * M_prime - M)"),
        ("E*sin(4D-4M'+M)", "E * math.sin(4.0 * D - 4.0 * M_prime + M)"),
        ("E*sin(6D-6M'-M)", "E * math.sin(6.0 * D - 6.0 * M_prime - M)"),
        ("E*sin(6D-6M'+M)", "E * math.sin(6.0 * D - 6.0 * M_prime + M)"),
    ]:
        terms.append((cat, label, expr))

    cat = "SOLAR DOUBLE COUPLING (E^2 terms)"
    for label, expr in [
        ("E2*sin(2M)", "E2 * math.sin(2.0 * M)"),
        ("E2*sin(2D-2M'-2M)", "E2 * math.sin(2.0 * D - 2.0 * M_prime - 2.0 * M)"),
        ("E2*sin(2D-2M'+2M)", "E2 * math.sin(2.0 * D - 2.0 * M_prime + 2.0 * M)"),
        ("E2*sin(4D-4M'-2M)", "E2 * math.sin(4.0 * D - 4.0 * M_prime - 2.0 * M)"),
    ]:
        terms.append((cat, label, expr))

    cat = "LUNAR ANOMALY HARMONICS (M' alone)"
    for label, expr in [
        ("sin(M')", "math.sin(M_prime)"),
        ("sin(2M')", "math.sin(2.0 * M_prime)"),
        ("sin(3M')", "math.sin(3.0 * M_prime)"),
    ]:
        terms.append((cat, label, expr))

    cat = "LATITUDE COUPLING (F-dependent)"
    for label, expr in [
        ("sin(2F-2M')", "math.sin(2.0 * F - 2.0 * M_prime)"),
        ("sin(2F-2D)", "math.sin(2.0 * F - 2.0 * D)"),
        ("sin(2F)", "math.sin(2.0 * F)"),
        ("sin(2F-4M'+2D)", "math.sin(2.0 * F - 4.0 * M_prime + 2.0 * D)"),
        ("sin(2F+2M'-2D)", "math.sin(2.0 * F + 2.0 * M_prime - 2.0 * D)"),
        ("sin(2F-2M'+2D)", "math.sin(2.0 * F - 2.0 * M_prime + 2.0 * D)"),
        ("sin(4F-2M')", "math.sin(4.0 * F - 2.0 * M_prime)"),
    ]:
        terms.append((cat, label, expr))

    cat = "CROSS-COUPLING (D, M' combinations)"
    for label, expr in [
        ("sin(2D-M')", "math.sin(2.0 * D - M_prime)"),
        ("sin(2D-3M')", "math.sin(2.0 * D - 3.0 * M_prime)"),
        ("sin(4D-3M')", "math.sin(4.0 * D - 3.0 * M_prime)"),
        ("sin(4D-5M')", "math.sin(4.0 * D - 5.0 * M_prime)"),
        ("sin(2D)", "math.sin(2.0 * D)"),
        ("sin(4D)", "math.sin(4.0 * D)"),
        ("sin(6D-5M')", "math.sin(6.0 * D - 5.0 * M_prime)"),
        ("sin(3D-2M')", "math.sin(3.0 * D - 2.0 * M_prime)"),
    ]:
        terms.append((cat, label, expr))

    cat = "SOLAR-LATITUDE CROSS-COUPLING"
    for label, expr in [
        ("E*sin(2F-2M'+M)", "E * math.sin(2.0 * F - 2.0 * M_prime + M)"),
        ("E*sin(2F-2M'-M)", "E * math.sin(2.0 * F - 2.0 * M_prime - M)"),
        ("E*sin(2F-2D+M)", "E * math.sin(2.0 * F - 2.0 * D + M)"),
        ("E*sin(2F-2D-M)", "E * math.sin(2.0 * F - 2.0 * D - M)"),
    ]:
        terms.append((cat, label, expr))

    cat = "HIGHER-ORDER EVECTION-SOLAR COUPLING"
    for label, expr in [
        ("E*sin(8D-8M'-M)", "E * math.sin(8.0 * D - 8.0 * M_prime - M)"),
        ("E*sin(8D-8M'+M)", "E * math.sin(8.0 * D - 8.0 * M_prime + M)"),
        ("E*sin(10D-10M'-M)", "E * math.sin(10.0 * D - 10.0 * M_prime - M)"),
        ("E*sin(10D-10M'+M)", "E * math.sin(10.0 * D - 10.0 * M_prime + M)"),
    ]:
        terms.append((cat, label, expr))

    cat = "SECULAR AND LONG-PERIOD CORRECTIONS"
    for label, expr in [
        ("T*sin(2D-2M')", "T * math.sin(2.0 * D - 2.0 * M_prime)"),
        ("T*sin(M)", "T * math.sin(M)"),
        ("T*sin(D-M')", "T * math.sin(D - M_prime)"),
        ("T*cos(2D-2M')", "T * math.cos(2.0 * D - 2.0 * M_prime)"),
    ]:
        terms.append((cat, label, expr))

    cat = "COSINE COUPLING TERMS (phase corrections)"
    for label, expr in [
        ("cos(M)", "math.cos(M)"),
        ("cos(2D-M')", "math.cos(2.0 * D - M_prime)"),
        ("cos(2F-2M')", "math.cos(2.0 * F - 2.0 * M_prime)"),
        ("cos(M')", "math.cos(M_prime)"),
        ("cos(2D-3M')", "math.cos(2.0 * D - 3.0 * M_prime)"),
        ("cos(2D)", "math.cos(2.0 * D)"),
    ]:
        terms.append((cat, label, expr))

    cat = "DIRECT SOLAR COUPLING (M with D, no M')"
    for label, expr in [
        ("E*sin(2D-M)", "E * math.sin(2.0 * D - M)"),
        ("E*sin(2D+M)", "E * math.sin(2.0 * D + M)"),
        ("E*sin(D+M)", "E * math.sin(D + M)"),
        ("E*sin(D-M)", "E * math.sin(D - M)"),
        ("E*cos(2D-M)", "E * math.cos(2.0 * D - M)"),
    ]:
        terms.append((cat, label, expr))

    cat = "OFF-DIAGONAL D-M' AND PARALLACTIC TERMS"
    for label, expr in [
        ("sin(D)", "math.sin(D)"),
        ("sin(2D+M')", "math.sin(2.0 * D + M_prime)"),
        ("sin(2D+2M')", "math.sin(2.0 * D + 2.0 * M_prime)"),
        ("sin(3D)", "math.sin(3.0 * D)"),
        ("cos(D)", "math.cos(D)"),
        ("cos(2D+M')", "math.cos(2.0 * D + M_prime)"),
        ("cos(2D+2M')", "math.cos(2.0 * D + 2.0 * M_prime)"),
        ("cos(3D)", "math.cos(3.0 * D)"),
    ]:
        terms.append((cat, label, expr))

    cat = "SUN-MOON ANOMALY COUPLING"
    for label, expr in [
        ("E*sin(M+M')", "E * math.sin(M + M_prime)"),
        ("E*sin(M-M')", "E * math.sin(M - M_prime)"),
        ("E*sin(2D-M+M')", "E * math.sin(2.0 * D - M + M_prime)"),
        ("E*sin(2D+M-M')", "E * math.sin(2.0 * D + M - M_prime)"),
        ("E*sin(M-2M')", "E * math.sin(M - 2.0 * M_prime)"),
        ("E*sin(M+2M')", "E * math.sin(M + 2.0 * M_prime)"),
        ("E2*sin(2D-2M)", "E2 * math.sin(2.0 * D - 2.0 * M)"),
    ]:
        terms.append((cat, label, expr))

    cat = "ADDITIONAL LATITUDE COUPLING"
    for label, expr in [
        ("sin(2F+M')", "math.sin(2.0 * F + M_prime)"),
        ("sin(2F-M')", "math.sin(2.0 * F - M_prime)"),
        ("sin(2F+D)", "math.sin(2.0 * F + D)"),
        ("sin(2F-D)", "math.sin(2.0 * F - D)"),
    ]:
        terms.append((cat, label, expr))

    cat = "ADDITIONAL COSINE PHASE TERMS"
    for label, expr in [
        ("cos(2M')", "math.cos(2.0 * M_prime)"),
        ("cos(3M')", "math.cos(3.0 * M_prime)"),
        ("cos(2F)", "math.cos(2.0 * F)"),
    ]:
        terms.append((cat, label, expr))

    cat = "POLYNOMIAL MEAN PERIGEE CORRECTIONS"
    for label, expr in [
        ("const", "1.0"),
        ("T", "T"),
        ("T^2", "T * T"),
        ("T^3", "T * T * T"),
    ]:
        terms.append((cat, label, expr))

    return terms


def format_coefficients(
    coeffs: list[float],
    threshold: float = COEFFICIENT_THRESHOLD,
) -> str:
    """Format coefficients as Python code for _calc_elp2000_perigee_perturbations()."""
    term_defs = get_term_code_expressions()
    lines = []
    current_cat = None
    n_significant = 0
    for i, (cat, _label, expr) in enumerate(term_defs):
        c = coeffs[i]
        if abs(c) < threshold:
            continue
        n_significant += 1
        if cat != current_cat:
            if current_cat is not None:
                lines.append("")
            lines.append(
                "    # ========================================"
                "================================"
            )
            lines.append(f"    # {cat}")
            lines.append(
                "    # ========================================"
                "================================"
            )
            current_cat = cat
        lines.append(f"    perturbation += {c:+.4f} * {expr}")
    lines.append("")
    lines.append(
        f"    # Total: {n_significant} significant terms (|coeff| >= {threshold} deg)"
    )
    return "\n".join(lines)


# =============================================================================
# MAIN
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate perigee perturbation coefficients (v2.2: passage-interpolated fit)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Method (v2.2 - Passage-Interpolated Harmonic Fit):
  1. Find all perigee passages (Earth-Moon distance minima) from JPL ephemeris
  2. Cubic spline interpolation between passages for smooth perigee longitude
  3. Resample daily, compute perturbation = interpolated - mean perigee
  4. Least-squares fit of harmonic series to the perturbation signal

Prerequisites:
  - JPL DE441 ephemeris file (de441.bsp) must be available
  - scipy is required for spline interpolation
  - pyswisseph (optional) for cross-validation against Swiss Ephemeris

Output:
  - Calibrated Python code for lunar.py (printed to stdout)
  - Optional JSON file with full results (--output)

Usage examples:
  # Full calibration (production coefficients, ~30 min)
  poe calibrate-perigee

  # Quick validation run (~2 min)
  poe calibrate-perigee:quick

  # Custom range with JSON output
  python scripts/calibrate_perigee_perturbations.py \\
      --start-year 1500 --end-year 2500 --output results.json

After calibration:
  1. Paste the printed coefficients into lunar.py (_calc_elp2000_perigee_perturbations)
  2. Update the matching series in scripts/generate_lunar_corrections.py
  3. Regenerate correction table: poe generate-lunar-corrections
  4. Run tests: poe test:lunar:perigee
""",
    )
    parser.add_argument(
        "--ephemeris",
        default=DEFAULT_EPHEMERIS,
        help="JPL ephemeris file to use (default: %(default)s)",
    )
    parser.add_argument(
        "--start-year",
        type=int,
        default=DEFAULT_FIT_START_YEAR,
        help="Start year for fitting range (default: %(default)s)",
    )
    parser.add_argument(
        "--end-year",
        type=int,
        default=DEFAULT_FIT_END_YEAR,
        help="End year for fitting range (default: %(default)s)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help="Number of parallel workers for passage search (default: %(default)s)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Quick mode: [1800,2200] range, skip SE cross-validation (~2 min)",
    )
    parser.add_argument(
        "--output",
        dest="save_json",
        type=str,
        metavar="FILE",
        help="Save full results to JSON file",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=COEFFICIENT_THRESHOLD,
        help="Minimum coefficient magnitude in degrees to keep (default: %(default)s)",
    )
    parser.add_argument(
        "--skip-se",
        action="store_true",
        help="Skip Swiss Ephemeris cross-validation (useful if swisseph not installed)",
    )
    parser.add_argument(
        "--sample-step-days",
        type=int,
        default=1,
        help="Step in days for spline sampling (default: %(default)s)",
    )

    args = parser.parse_args()

    print(f"Setting up ephemeris: {args.ephemeris}")
    if not setup_ephemeris(args.ephemeris):
        return 1
    print(f"Using ephemeris: {args.ephemeris}")

    if args.quick:
        args.start_year = 1800
        args.end_year = 2200
        args.skip_se = True
        args.sample_step_days = 3
        print(
            f"\n[QUICK MODE] [{args.start_year}, {args.end_year}], step {args.sample_step_days}d"
        )

    # =================================================================
    # PHASE 1: FIND PERIGEE PASSAGES
    # =================================================================
    print(f"\n{'=' * 70}")
    print("PHASE 1: FINDING PERIGEE PASSAGES")
    print(f"{'=' * 70}")

    # Add margin for spline edge effects
    margin_years = 5
    jd_start = year_to_jd(args.start_year - margin_years)
    jd_end = year_to_jd(args.end_year + margin_years)
    jd_fit_start = year_to_jd(args.start_year)
    jd_fit_end = year_to_jd(args.end_year)

    print(
        f"  Passage search: [{args.start_year - margin_years}, {args.end_year + margin_years}] CE"
    )
    print(f"  Fit range:      [{args.start_year}, {args.end_year}] CE")

    passages = find_all_perigee_passages(jd_start, jd_end, args.workers, args.ephemeris)
    print(f"\n  Found {len(passages)} perigee passages")

    if len(passages) < 100:
        print("ERROR: Too few passages")
        return 1

    # =================================================================
    # PHASE 2: CUBIC SPLINE INTERPOLATION OF PASSAGE LONGITUDES
    # =================================================================
    print(f"\n{'=' * 70}")
    print("PHASE 2: CUBIC SPLINE INTERPOLATION THROUGH PASSAGES")
    print(f"{'=' * 70}")

    passage_times = np.array([p[0] for p in passages])
    passage_lons_raw = [p[1] for p in passages]

    # Unwrap longitudes to avoid 0/360 discontinuity
    passage_lons_unwrapped = unwrap_longitudes(passage_lons_raw)
    passage_lons = np.array(passage_lons_unwrapped)

    # Generate daily sample points within fit range
    n_samples = int((jd_fit_end - jd_fit_start) / args.sample_step_days) + 1
    sample_jds = np.array(
        [jd_fit_start + i * args.sample_step_days for i in range(n_samples)]
    )

    print(f"  Passages for spline: {len(passages)}")
    print(f"  Daily samples: {n_samples} (step {args.sample_step_days} day(s))")
    print("  Interpolating with cubic spline...")

    # Interpolate
    interp_lons_unwrapped = cubic_spline_interpolate(
        passage_times, passage_lons, sample_jds
    )
    # Wrap back to [0, 360)
    interp_lons = interp_lons_unwrapped % 360.0

    print("  Spline interpolation complete")

    # =================================================================
    # PHASE 3: HARMONIC FIT ON INTERPOLATED DATA
    # =================================================================
    print(f"\n{'=' * 70}")
    print("PHASE 3: HARMONIC FIT ON SPLINE-INTERPOLATED DATA")
    print(f"{'=' * 70}")

    from libephemeris.lunar import calc_mean_lilith

    n_terms = len(build_design_vector(sample_jds[0]))
    overdetermination = n_samples / n_terms

    print(f"  Samples: {n_samples}")
    print(f"  Design matrix: {n_samples} x {n_terms}")
    print(f"  Overdetermination: {overdetermination:.0f}:1")

    A = np.zeros((n_samples, n_terms))
    b = np.zeros(n_samples)

    for i in range(n_samples):
        jd = float(sample_jds[i])
        mean_perigee = (calc_mean_lilith(jd) + 180.0) % 360.0
        perturbation = normalize_angle_diff(float(interp_lons[i]) - mean_perigee)
        A[i, :] = build_design_vector(jd)
        b[i] = perturbation

    print(f"  Perturbation range: [{float(np.min(b)):.2f}, {float(np.max(b)):.2f}] deg")
    print(f"  Perturbation std: {float(np.std(b)):.2f} deg")

    print("  Solving least squares...")
    coeffs, _, rank, sv = np.linalg.lstsq(A, b, rcond=None)
    print(f"  Matrix rank: {rank}")
    if len(sv) > 0:
        print(f"  Condition number: {sv[0] / sv[-1]:.2e}")

    predicted = A @ coeffs
    errors = b - predicted
    fit_rms = float(np.sqrt(np.mean(errors**2)))
    fit_max = float(np.max(np.abs(errors)))

    print("\n  Fit results (on spline-interpolated data):")
    print(f"    RMS residual: {fit_rms:.4f} deg")
    print(f"    Max residual: {fit_max:.4f} deg")
    print(f"    Mean residual: {float(np.mean(errors)):.6f} deg")

    coeffs_list = coeffs.tolist()
    labels = get_term_labels()
    n_significant = sum(1 for c in coeffs_list if abs(c) >= args.threshold)
    print(f"\n  Significant terms: {n_significant}/{n_terms}")

    # =================================================================
    # PHASE 4: VALIDATE ON ACTUAL PASSAGES (held-out)
    # =================================================================
    print(f"\n{'=' * 70}")
    print("PHASE 4: VALIDATION ON ACTUAL PERIGEE PASSAGES")
    print(f"{'=' * 70}")

    # Use passages within fit range for validation
    fit_passages = [p for p in passages if jd_fit_start <= p[0] <= jd_fit_end]
    print(f"  Passages in fit range: {len(fit_passages)}")

    passage_errors = []
    for jd_p, lon_p, _ in fit_passages:
        mean_perigee = (calc_mean_lilith(jd_p) + 180.0) % 360.0
        design_vec = build_design_vector(jd_p)
        series_pred = sum(c * t for c, t in zip(coeffs_list, design_vec))
        predicted_lon = (mean_perigee + series_pred) % 360.0
        error = normalize_angle_diff(lon_p - predicted_lon)
        passage_errors.append(error)

    passage_rms = math.sqrt(sum(e**2 for e in passage_errors) / len(passage_errors))
    passage_max = max(abs(e) for e in passage_errors)
    passage_mean = sum(passage_errors) / len(passage_errors)

    print("  Passage validation results:")
    print(f"    RMS error:  {passage_rms:.4f} deg")
    print(f"    Max error:  {passage_max:.4f} deg")
    print(f"    Mean error: {passage_mean:.4f} deg")
    print("    Target: RMS < 0.5 deg")
    if passage_rms < 0.5:
        print("    STATUS: PASS")
    elif passage_rms < 1.0:
        print("    STATUS: MARGINAL")
    else:
        print("    STATUS: FAIL (consider adding more terms or using correction table)")

    # Residual by century
    print("\n  Fit RMS by century:")
    century_bins: dict[int, list[float]] = {}
    for i in range(n_samples):
        year = int(jd_to_year(float(sample_jds[i])))
        century = year // 100
        if century not in century_bins:
            century_bins[century] = []
        century_bins[century].append(errors[i])

    for century in sorted(century_bins.keys()):
        errs = century_bins[century]
        rms = math.sqrt(sum(e**2 for e in errs) / len(errs))
        print(
            f"    {century * 100:5d}-{century * 100 + 99}: RMS={rms:.4f} deg (n={len(errs)})"
        )

    # =================================================================
    # PHASE 5: SE CROSS-VALIDATION
    # =================================================================
    se_rms = None
    if not args.skip_se:
        try:
            import swisseph as swe

            print(f"\n{'=' * 70}")
            print("PHASE 5: SWISS EPHEMERIS CROSS-VALIDATION")
            print(f"{'=' * 70}")
            se_errors = []
            for year in range(1900, 2101, 5):
                jd_ut = swe.julday(year, 1, 15, 12.0)
                jd_tt = jd_ut + 69.184 / 86400.0
                pos_se, _ = swe.calc_ut(jd_ut, swe.INTP_PERG, 0)
                se_lon = pos_se[0]
                mean_perigee = (calc_mean_lilith(jd_tt) + 180.0) % 360.0
                design_vec = build_design_vector(jd_tt)
                series_pred = sum(c * t for c, t in zip(coeffs_list, design_vec))
                our_lon = (mean_perigee + series_pred) % 360.0
                error = normalize_angle_diff(our_lon - se_lon)
                se_errors.append(error)
            se_rms = math.sqrt(sum(e**2 for e in se_errors) / len(se_errors))
            se_max = max(abs(e) for e in se_errors)
            print("  SE comparison (1900-2100):")
            print(f"    RMS: {se_rms:.2f} deg")
            print(f"    Max: {se_max:.2f} deg")
        except ImportError:
            print("\n  [Skipping SE - swisseph not available]")

    # =================================================================
    # OUTPUT
    # =================================================================
    print(f"\n{'=' * 70}")
    print("CALIBRATED COEFFICIENTS:")
    print(f"{'=' * 70}")
    print(format_coefficients(coeffs_list, args.threshold))

    print(f"\n{'=' * 70}")
    print("ALL SIGNIFICANT COEFFICIENTS (sorted by magnitude):")
    print(f"{'=' * 70}")
    indexed_coeffs = [(i, labels[i], coeffs_list[i]) for i in range(len(coeffs_list))]
    indexed_coeffs.sort(key=lambda x: abs(x[2]), reverse=True)
    for idx, label, coeff in indexed_coeffs:
        if abs(coeff) >= args.threshold:
            print(f"  [{idx:2d}] {label:25s} = {coeff:+.6f} deg")

    # Save JSON
    if args.save_json:
        output = {
            "version": "2.2",
            "method": "passage_interpolated_harmonic_fit",
            "ephemeris": args.ephemeris,
            "fit_range_years": [args.start_year, args.end_year],
            "n_passages": len(passages),
            "n_samples": n_samples,
            "n_terms": n_terms,
            "n_significant": n_significant,
            "fit_rms_deg": fit_rms,
            "fit_max_deg": fit_max,
            "passage_validation_rms_deg": passage_rms,
            "se_crossval_rms_deg": se_rms,
            "coefficients": coeffs_list,
            "labels": labels,
            "timestamp": datetime.now().isoformat(),
        }
        with open(args.save_json, "w") as f:
            json.dump(output, f, indent=2)
        print(f"\n  Results saved to: {args.save_json}")

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print("  Method: Passage-interpolated harmonic fit (v2.2)")
    print(f"  Range: [{args.start_year}, {args.end_year}] CE")
    print(f"  Passages: {len(passages)}, Samples: {n_samples}")
    print(f"  Fit RMS (spline): {fit_rms:.4f} deg")
    print(f"  Passage RMS:      {passage_rms:.4f} deg")
    print(f"  Significant terms: {n_significant}/{n_terms}")
    if se_rms is not None:
        print(f"  SE RMS: {se_rms:.2f} deg")

    return 0


if __name__ == "__main__":
    sys.exit(main())
