#!/usr/bin/env python3
"""
Calibrate perigee perturbation coefficients against JPL DE441 ephemeris.

This script computes optimal coefficients for a 67-term trigonometric
perturbation series for the interpolated lunar perigee by fitting against
interpolated perigee positions derived from JPL DE441/DE441 ephemeris.

Method:
    The "interpolated perigee" is computed via quadratic polynomial regression
    on 9 samples of osculating perigee in a +/-28 day window. This smooths
    spurious short-period oscillations from the osculating elements.

    Perturbation = interpolated_perigee_JPL - mean_perigee

    A design matrix of 67 trigonometric terms is built and solved via
    least-squares to find optimal coefficients.

Output:
    - Calibrated coefficients (Python code for lunar.py)
    - Precision report (RMS/max error vs JPL)

References:
    - Park, R.S. et al. (2021) "The JPL Planetary and Lunar Ephemerides DE440 and DE441"
    - Chapront-Touze, M. & Chapront, J. "ELP 2000-82B" (1988)
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from multiprocessing import cpu_count
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    pass

DEFAULT_EPHEMERIS = "de441.bsp"
DEFAULT_FIT_START_YEAR = -2000
DEFAULT_FIT_END_YEAR = 4000
DEFAULT_FIT_STEP_DAYS = 14
DEFAULT_WORKERS = cpu_count()

HALF_WINDOW_DAYS = 28.0
N_SAMPLES = 9


def normalize_angle_diff(diff: float) -> float:
    while diff >= 180.0:
        diff -= 360.0
    while diff < -180.0:
        diff += 360.0
    return diff


def year_to_jd(year: float) -> float:
    return 2451545.0 + (year - 2000.0) * 365.25


def jd_to_year(jd: float) -> float:
    return (jd - 2451545.0) / 365.25 + 2000.0


def unwrap_longitudes(longitudes: list[float]) -> list[float]:
    if not longitudes:
        return []
    unwrapped = [longitudes[0]]
    for i in range(1, len(longitudes)):
        diff = longitudes[i] - longitudes[i - 1]
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360
        unwrapped.append(unwrapped[-1] + diff)
    return unwrapped


def compute_interpolated_perigee_jpl(jd_tt: float) -> float:
    from libephemeris.lunar import calc_osculating_perigee

    step = 2.0 * HALF_WINDOW_DAYS / (N_SAMPLES - 1)
    times = []
    lons = []
    for i in range(N_SAMPLES):
        t = jd_tt - HALF_WINDOW_DAYS + i * step
        lon, _, _ = calc_osculating_perigee(t)
        times.append(t - jd_tt)
        lons.append(lon)

    lons = unwrap_longitudes(lons)

    n = len(times)
    sum_t0 = float(n)
    sum_t1 = sum(times)
    sum_t2 = sum(t**2 for t in times)
    sum_t3 = sum(t**3 for t in times)
    sum_t4 = sum(t**4 for t in times)
    sum_y = sum(lons)
    sum_ty = sum(t * y for t, y in zip(times, lons))
    sum_t2y = sum(t**2 * y for t, y in zip(times, lons))

    A = [
        [sum_t4, sum_t3, sum_t2],
        [sum_t3, sum_t2, sum_t1],
        [sum_t2, sum_t1, sum_t0],
    ]
    rhs = [sum_t2y, sum_ty, sum_y]

    for col in range(3):
        max_row = col
        for row in range(col + 1, 3):
            if abs(A[row][col]) > abs(A[max_row][col]):
                max_row = row
        A[col], A[max_row] = A[max_row], A[col]
        rhs[col], rhs[max_row] = rhs[max_row], rhs[col]
        for row in range(col + 1, 3):
            if abs(A[col][col]) < 1e-30:
                continue
            factor = A[row][col] / A[col][col]
            for k in range(col, 3):
                A[row][k] -= factor * A[col][k]
            rhs[row] -= factor * rhs[col]

    x = [0.0, 0.0, 0.0]
    for i in range(2, -1, -1):
        x[i] = rhs[i]
        for j in range(i + 1, 3):
            x[i] -= A[i][j] * x[j]
        if abs(A[i][i]) > 1e-30:
            x[i] /= A[i][i]

    return x[2] % 360.0


def compute_mean_perigee_analytical(jd_tt: float) -> float:
    from libephemeris.lunar import calc_mean_lilith

    return (calc_mean_lilith(jd_tt) + 180.0) % 360.0


def compute_fundamental_arguments(jd_tt: float) -> tuple[float, float, float, float]:
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


def build_perturbation_vector(jd_tt: float) -> list[float]:
    D, M, M_prime, F = compute_fundamental_arguments(jd_tt)

    T = (jd_tt - 2451545.0) / 36525.0
    E = 1.0 - 0.002516 * T - 0.0000074 * T**2
    E2 = E * E

    terms = []

    for k in range(1, 19):
        terms.append(math.sin(k * D - k * M_prime))

    terms.append(E * math.sin(M))
    terms.append(E2 * math.sin(2.0 * M))
    terms.append(E * math.sin(2.0 * D - 2.0 * M_prime - M))
    terms.append(E * math.sin(2.0 * D - 2.0 * M_prime + M))
    terms.append(E * math.sin(D - M_prime - M))
    terms.append(E * math.sin(D - M_prime + M))
    terms.append(E * math.sin(4.0 * D - 4.0 * M_prime - M))
    terms.append(E * math.sin(4.0 * D - 4.0 * M_prime + M))
    terms.append(E * math.sin(6.0 * D - 6.0 * M_prime - M))
    terms.append(E2 * math.sin(2.0 * D - 2.0 * M_prime - 2.0 * M))
    terms.append(E2 * math.sin(2.0 * D - 2.0 * M_prime + 2.0 * M))

    terms.append(math.sin(M_prime))
    terms.append(math.sin(2.0 * M_prime))
    terms.append(math.sin(3.0 * M_prime))

    terms.append(math.sin(2.0 * F - 2.0 * M_prime))
    terms.append(math.sin(2.0 * F - 2.0 * D))
    terms.append(math.sin(2.0 * F))
    terms.append(math.sin(2.0 * F - 4.0 * M_prime + 2.0 * D))
    terms.append(math.sin(2.0 * F + 2.0 * M_prime - 2.0 * D))

    terms.append(math.sin(2.0 * D - M_prime))
    terms.append(math.sin(2.0 * D - 3.0 * M_prime))
    terms.append(math.sin(4.0 * D - 3.0 * M_prime))
    terms.append(math.sin(4.0 * D - 5.0 * M_prime))
    terms.append(math.sin(2.0 * D))
    terms.append(math.sin(4.0 * D))
    terms.append(E * math.sin(2.0 * F - 2.0 * M_prime + M))
    terms.append(E * math.sin(2.0 * F - 2.0 * M_prime - M))
    terms.append(E * math.sin(2.0 * F - 2.0 * D + M))
    terms.append(E * math.sin(2.0 * F - 2.0 * D - M))

    terms.append(T * math.sin(2.0 * D - 2.0 * M_prime))
    terms.append(T * math.sin(M))

    for k in range(1, 19):
        terms.append(math.cos(k * D - k * M_prime))

    return terms


def compute_sample_jd(jd_tt: float) -> tuple[float, float, list[float]] | None:
    try:
        interp_perigee = compute_interpolated_perigee_jpl(jd_tt)
    except Exception:
        return None

    mean_perigee = compute_mean_perigee_analytical(jd_tt)
    perturbation = normalize_angle_diff(interp_perigee - mean_perigee)
    pert_vector = build_perturbation_vector(jd_tt)

    return (jd_tt, perturbation, pert_vector)


def format_coefficients(coeffs: list[float]) -> str:
    lines = []
    idx = 0
    n_evection = 18
    threshold = 0.005

    lines.append(
        "    # ========================================================================"
    )
    lines.append("    # PRIMARY EVECTION HARMONICS (kD - kM')")
    lines.append(
        "    # ========================================================================"
    )
    for k in range(1, n_evection + 1):
        c = coeffs[idx]
        idx += 1
        if abs(c) > threshold:
            if k == 1:
                lines.append(f"    perturbation += {c:+.4f} * math.sin(D - M_prime)")
            else:
                lines.append(
                    f"    perturbation += {c:+.4f} * math.sin({k}.0 * D - {k}.0 * M_prime)"
                )

    lines.append("")
    lines.append(
        "    # ========================================================================"
    )
    lines.append("    # SOLAR ANOMALY COUPLING (M terms)")
    lines.append(
        "    # ========================================================================"
    )
    solar_exprs = [
        "E * math.sin(M)",
        "E2 * math.sin(2.0 * M)",
        "E * math.sin(2.0 * D - 2.0 * M_prime - M)",
        "E * math.sin(2.0 * D - 2.0 * M_prime + M)",
        "E * math.sin(D - M_prime - M)",
        "E * math.sin(D - M_prime + M)",
        "E * math.sin(4.0 * D - 4.0 * M_prime - M)",
        "E * math.sin(4.0 * D - 4.0 * M_prime + M)",
        "E * math.sin(6.0 * D - 6.0 * M_prime - M)",
        "E2 * math.sin(2.0 * D - 2.0 * M_prime - 2.0 * M)",
        "E2 * math.sin(2.0 * D - 2.0 * M_prime + 2.0 * M)",
    ]
    for expr in solar_exprs:
        c = coeffs[idx]
        idx += 1
        if abs(c) > threshold:
            lines.append(f"    perturbation += {c:+.4f} * {expr}")

    lines.append("")
    lines.append(
        "    # ========================================================================"
    )
    lines.append("    # LUNAR ANOMALY HARMONICS (M' alone)")
    lines.append(
        "    # ========================================================================"
    )
    lunar_exprs = [
        "math.sin(M_prime)",
        "math.sin(2.0 * M_prime)",
        "math.sin(3.0 * M_prime)",
    ]
    for expr in lunar_exprs:
        c = coeffs[idx]
        idx += 1
        if abs(c) > threshold:
            lines.append(f"    perturbation += {c:+.4f} * {expr}")

    lines.append("")
    lines.append(
        "    # ========================================================================"
    )
    lines.append("    # LATITUDE COUPLING TERMS (F-dependent)")
    lines.append(
        "    # ========================================================================"
    )
    lat_exprs = [
        "math.sin(2.0 * F - 2.0 * M_prime)",
        "math.sin(2.0 * F - 2.0 * D)",
        "math.sin(2.0 * F)",
        "math.sin(2.0 * F - 4.0 * M_prime + 2.0 * D)",
        "math.sin(2.0 * F + 2.0 * M_prime - 2.0 * D)",
    ]
    for expr in lat_exprs:
        c = coeffs[idx]
        idx += 1
        if abs(c) > threshold:
            lines.append(f"    perturbation += {c:+.4f} * {expr}")

    lines.append("")
    lines.append(
        "    # ========================================================================"
    )
    lines.append("    # CROSS-COUPLING TERMS (D, M, M', F combinations)")
    lines.append(
        "    # ========================================================================"
    )
    cross_exprs = [
        "math.sin(2.0 * D - M_prime)",
        "math.sin(2.0 * D - 3.0 * M_prime)",
        "math.sin(4.0 * D - 3.0 * M_prime)",
        "math.sin(4.0 * D - 5.0 * M_prime)",
        "math.sin(2.0 * D)",
        "math.sin(4.0 * D)",
        "E * math.sin(2.0 * F - 2.0 * M_prime + M)",
        "E * math.sin(2.0 * F - 2.0 * M_prime - M)",
        "E * math.sin(2.0 * F - 2.0 * D + M)",
        "E * math.sin(2.0 * F - 2.0 * D - M)",
    ]
    for expr in cross_exprs:
        c = coeffs[idx]
        idx += 1
        if abs(c) > threshold:
            lines.append(f"    perturbation += {c:+.4f} * {expr}")

    lines.append("")
    lines.append(
        "    # ========================================================================"
    )
    lines.append("    # SECULAR AND LONG-PERIOD CORRECTIONS")
    lines.append(
        "    # ========================================================================"
    )
    secular_exprs = [
        "T * math.sin(2.0 * D - 2.0 * M_prime)",
        "T * math.sin(M)",
    ]
    for expr in secular_exprs:
        c = coeffs[idx]
        idx += 1
        if abs(c) > threshold:
            lines.append(f"    perturbation += {c:+.4f} * {expr}")

    lines.append("")
    lines.append(
        "    # ========================================================================"
    )
    lines.append("    # COSINE TERMS (phase corrections) - cos(kD - kM')")
    lines.append(
        "    # ========================================================================"
    )
    for k in range(1, n_evection + 1):
        c = coeffs[idx]
        idx += 1
        if abs(c) > threshold:
            if k == 1:
                lines.append(f"    perturbation += {c:+.4f} * math.cos(D - M_prime)")
            else:
                lines.append(
                    f"    perturbation += {c:+.4f} * math.cos({k}.0 * D - {k}.0 * M_prime)"
                )

    return "\n".join(lines)


def init_worker(ephemeris: str) -> None:
    os.environ["LIBEPHEMERIS_EPHEMERIS"] = ephemeris
    from libephemeris.state import set_ephemeris_file

    set_ephemeris_file(ephemeris)


def setup_ephemeris(ephemeris: str) -> bool:
    from libephemeris.state import set_ephemeris_file

    try:
        set_ephemeris_file(ephemeris)
        return True
    except Exception as e:
        print(f"Error setting ephemeris: {e}")
        return False


def parallel_compute_samples(
    jd_list: list[float],
    workers: int,
    ephemeris: str,
    label: str = "samples",
) -> dict[float, tuple[float, list[float]]]:
    results = {}
    failed = 0

    tqdm_cls = None
    try:
        from tqdm import tqdm as tqdm_module

        tqdm_cls = tqdm_module
    except ImportError:
        pass

    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=init_worker,
        initargs=(ephemeris,),
    ) as executor:
        futures = {executor.submit(compute_sample_jd, jd): jd for jd in jd_list}

        if tqdm_cls is not None:
            for future in tqdm_cls(
                as_completed(futures), total=len(futures), desc=f"Computing {label}"
            ):
                result = future.result()
                if result is None:
                    failed += 1
                else:
                    jd, perturbation, pert_vector = result
                    results[jd] = (perturbation, pert_vector)
        else:
            completed = 0
            for future in as_completed(futures):
                result = future.result()
                if result is None:
                    failed += 1
                else:
                    jd, perturbation, pert_vector = result
                    results[jd] = (perturbation, pert_vector)
                completed += 1
                if completed % 1000 == 0:
                    print(
                        f"  Progress: {completed}/{len(jd_list)} ({100 * completed // len(jd_list)}%)"
                    )

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate perigee perturbation coefficients against JPL DE441",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--ephemeris",
        default=DEFAULT_EPHEMERIS,
        help=f"Ephemeris file to use (default: {DEFAULT_EPHEMERIS})",
    )
    parser.add_argument(
        "--start-year",
        type=int,
        default=DEFAULT_FIT_START_YEAR,
        help=f"Start year for fit (default: {DEFAULT_FIT_START_YEAR})",
    )
    parser.add_argument(
        "--end-year",
        type=int,
        default=DEFAULT_FIT_END_YEAR,
        help=f"End year for fit (default: {DEFAULT_FIT_END_YEAR})",
    )
    parser.add_argument(
        "--fit-step-days",
        type=int,
        default=DEFAULT_FIT_STEP_DAYS,
        help=f"Step in days for coefficient fitting (default: {DEFAULT_FIT_STEP_DAYS})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Number of parallel workers (default: {DEFAULT_WORKERS})",
    )
    parser.add_argument(
        "--test-single",
        type=float,
        help="Test computation for a single year (for debugging)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Quick run: fit on [1800,2200] with 30-day step",
    )

    args = parser.parse_args()

    print(f"Setting up ephemeris: {args.ephemeris}")
    if not setup_ephemeris(args.ephemeris):
        print(f"Make sure {args.ephemeris} is available in the ephemeris path")
        return 1
    print(f"Using ephemeris: {args.ephemeris}")

    if args.test_single is not None:
        print(f"\nTesting single year: {args.test_single}")
        jd = year_to_jd(args.test_single)
        print(f"  JD: {jd}")

        print("  Computing interpolated perigee (quadratic regression)...")
        interp_perigee = compute_interpolated_perigee_jpl(jd)
        print(f"    Interpolated: {interp_perigee:.6f} deg")

        print("  Computing mean perigee (analytical)...")
        mean_perigee = compute_mean_perigee_analytical(jd)
        print(f"    Mean: {mean_perigee:.6f} deg")

        perturbation = normalize_angle_diff(interp_perigee - mean_perigee)
        print(f"  Perturbation: {perturbation:.6f} deg")

        print("  Building perturbation vector...")
        pert_vector = build_perturbation_vector(jd)
        print(f"    Vector has {len(pert_vector)} terms")

        return 0

    if args.quick:
        args.start_year = 1800
        args.end_year = 2200
        args.fit_step_days = 30
        print(
            f"\n[QUICK MODE] Fit range [{args.start_year}, {args.end_year}], step {args.fit_step_days} days"
        )

    print(f"\n{'=' * 70}")
    print("COEFFICIENT FITTING")
    print(f"{'=' * 70}")
    print(f"  Fit range: [{args.start_year}, {args.end_year}] CE")
    print(f"  Fit step: {args.fit_step_days} days")
    print(
        f"  Smoothing: quadratic regression on +/-{HALF_WINDOW_DAYS:.0f} days ({N_SAMPLES} samples)"
    )

    jd_fit_start = year_to_jd(args.start_year)
    jd_fit_end = year_to_jd(args.end_year)
    n_fit_samples = int((jd_fit_end - jd_fit_start) / args.fit_step_days) + 1
    fit_jds = [jd_fit_start + i * args.fit_step_days for i in range(n_fit_samples)]

    print(f"  Planned fit samples: {n_fit_samples}")
    print(f"  Parallel workers: {args.workers}")

    fit_results = parallel_compute_samples(
        fit_jds, args.workers, args.ephemeris, label="fit samples"
    )

    n_successful = len(fit_results)
    n_failed = n_fit_samples - n_successful
    print(f"\n  Successful fit samples: {n_successful}/{n_fit_samples}")
    if n_failed > 0:
        print(f"  Failed fit samples: {n_failed}")

    sorted_fit_jds = sorted(fit_results.keys())
    n_terms = len(list(fit_results.values())[0][1])
    print(f"  Design matrix: {n_successful} x {n_terms}")

    A_fit = np.zeros((n_successful, n_terms))
    b_fit = np.zeros(n_successful)

    for i, jd in enumerate(sorted_fit_jds):
        perturbation, pert_vector = fit_results[jd]
        A_fit[i, :] = pert_vector
        b_fit[i] = perturbation

    print("  Solving least squares with numpy...")
    coeffs, residuals_arr, rank, sv = np.linalg.lstsq(A_fit, b_fit, rcond=None)
    print(f"  Matrix rank: {rank}")
    print(
        f"  Condition number: {sv[0] / sv[-1]:.2e}"
        if len(sv) > 0
        else "  No singular values"
    )

    predicted_fit = A_fit @ coeffs
    fit_errors = b_fit - predicted_fit
    fit_rms = float(np.sqrt(np.mean(fit_errors**2)))
    fit_max = float(np.max(np.abs(fit_errors)))
    print("\n  Fit range results:")
    print(f"    RMS residual: {fit_rms:.4f} deg")
    print(f"    Max residual: {fit_max:.4f} deg")

    print("\n" + "=" * 70)
    print("CALIBRATED COEFFICIENTS FOR _calc_elp2000_perigee_perturbations:")
    print("=" * 70)
    coeffs_list = coeffs.tolist()
    print(format_coefficients(coeffs_list))

    print("\n" + "=" * 70)
    print("RAW COEFFICIENT VALUES (for reference):")
    print("=" * 70)
    labels = [f"sin({k}D-{k}M')" for k in range(1, 19)]
    labels += [
        "E*sin(M)",
        "E2*sin(2M)",
        "E*sin(2D-2M'-M)",
        "E*sin(2D-2M'+M)",
        "E*sin(D-M'-M)",
        "E*sin(D-M'+M)",
        "E*sin(4D-4M'-M)",
        "E*sin(4D-4M'+M)",
        "E*sin(6D-6M'-M)",
        "E2*sin(2D-2M'-2M)",
        "E2*sin(2D-2M'+2M)",
    ]
    labels += ["sin(M')", "sin(2M')", "sin(3M')"]
    labels += [
        "sin(2F-2M')",
        "sin(2F-2D)",
        "sin(2F)",
        "sin(2F-4M'+2D)",
        "sin(2F+2M'-2D)",
    ]
    labels += [
        "sin(2D-M')",
        "sin(2D-3M')",
        "sin(4D-3M')",
        "sin(4D-5M')",
        "sin(2D)",
        "sin(4D)",
        "E*sin(2F-2M'+M)",
        "E*sin(2F-2M'-M)",
        "E*sin(2F-2D+M)",
        "E*sin(2F-2D-M)",
    ]
    labels += ["T*sin(2D-2M')", "T*sin(M)"]
    labels += [f"cos({k}D-{k}M')" for k in range(1, 19)]

    for i, (label, coeff) in enumerate(zip(labels, coeffs_list)):
        if abs(coeff) > 0.001:
            print(f"  [{i:2d}] {label:25s} = {coeff:+.6f}")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Ephemeris: {args.ephemeris}")
    print(
        f"Fit range: [{args.start_year}, {args.end_year}], step {args.fit_step_days} days"
    )
    print(f"Fit samples: {n_successful}")
    print(f"Fit RMS: {fit_rms:.4f} deg")
    print(f"Fit Max: {fit_max:.4f} deg")
    print(f"Trigonometric terms: {n_terms}")
    print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
