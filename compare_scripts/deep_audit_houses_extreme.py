"""Deep precision audit: house cusp calculations at extreme latitudes.

Compares libephemeris vs pyswisseph for house cusps across multiple house
systems, latitudes (including extreme/polar), and dates.
"""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import math
import traceback

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path("swisseph/ephe")

# ============================================================================
# Configuration
# ============================================================================

LATITUDES = [0.0, 23.44, 45.0, 60.0, 66.56, 70.0, 75.0, 80.0, 85.0, 89.0]

HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
    ("E", "Equal"),
    ("W", "Whole Sign"),
    ("O", "Porphyry"),
]

JULIAN_DAYS = [
    (2451545.0, "J2000 (2000-01-01 12:00)"),
    (2460400.0, "Modern (2024-04-03 12:00)"),
]

THRESHOLD_WARN = 0.01  # degrees (36 arcsec) — flag if exceeded

LON = 0.0  # Fixed longitude for all tests


def angular_diff(a: float, b: float) -> float:
    """Signed-magnitude angular difference accounting for 360° wrap."""
    d = abs(a - b)
    if d > 180.0:
        d = 360.0 - d
    return d


# ============================================================================
# Main audit
# ============================================================================

print("=" * 80)
print("DEEP AUDIT: House Cusps at Extreme Latitudes")
print("libephemeris (Skyfield mode) vs pyswisseph")
print("=" * 80)
print()
print(f"Latitudes tested : {LATITUDES}")
print(f"House systems    : {[s for s, _ in HOUSE_SYSTEMS]}")
print(f"Julian Days      : {[jd for jd, _ in JULIAN_DAYS]}")
print(f"Longitude        : {LON}°")
print(f"Warning threshold: {THRESHOLD_WARN}° ({THRESHOLD_WARN * 3600:.0f} arcsec)")
print()

# Collect results: key = (system, lat) -> list of max diffs across dates
results: dict[tuple[str, float], list[float]] = {}
flagged: list[tuple[str, str, float, float, str, str]] = []
errors: list[tuple[str, str, float, float, str, str]] = []

for sys_char, sys_name in HOUSE_SYSTEMS:
    print(f"--- {sys_name} ({sys_char}) ---")
    print(
        f"  {'Lat':>6s}  {'JD':>12s}  {'MaxCusp':>10s}  {'AscDiff':>10s}  "
        f"{'MCDiff':>10s}  {'Status'}"
    )

    for lat in LATITUDES:
        for jd, jd_label in JULIAN_DAYS:
            status = "OK"
            max_cusp_diff = 0.0
            asc_diff = 0.0
            mc_diff = 0.0
            detail = ""

            # --- pyswisseph ---
            try:
                cusps_swe, ascmc_swe = swe.houses(
                    jd, lat, LON, sys_char.encode("ascii")
                )
                swe_ok = True
            except Exception as e:
                swe_ok = False
                swe_err = str(e)

            # --- libephemeris ---
            try:
                cusps_lib, ascmc_lib = ephem.swe_houses(jd, lat, LON, sys_char)
                lib_ok = True
            except Exception as e:
                lib_ok = False
                lib_err = str(e)

            # --- Compare ---
            if swe_ok and lib_ok:
                # Compare 12 cusps
                cusp_diffs = []
                for i in range(min(len(cusps_swe), len(cusps_lib), 12)):
                    cusp_diffs.append(angular_diff(cusps_swe[i], cusps_lib[i]))
                max_cusp_diff = max(cusp_diffs) if cusp_diffs else 0.0

                # Compare Asc and MC
                asc_diff = angular_diff(ascmc_swe[0], ascmc_lib[0])
                mc_diff = angular_diff(ascmc_swe[1], ascmc_lib[1])

                overall_max = max(max_cusp_diff, asc_diff, mc_diff)

                if overall_max > THRESHOLD_WARN:
                    status = "*** FLAGGED ***"
                    # Find worst cusp
                    worst_cusp = cusp_diffs.index(max_cusp_diff) + 1
                    detail = (
                        f"worst=cusp{worst_cusp} "
                        f"swe={cusps_swe[worst_cusp - 1]:.6f} "
                        f"lib={cusps_lib[worst_cusp - 1]:.6f}"
                    )
                    flagged.append(
                        (
                            sys_char,
                            sys_name,
                            lat,
                            jd,
                            jd_label,
                            f"max_diff={overall_max:.6f}° {detail}",
                        )
                    )

                # Store result
                key = (sys_char, lat)
                if key not in results:
                    results[key] = []
                results[key].append(overall_max)

            elif swe_ok and not lib_ok:
                status = "LIB_ERROR"
                detail = lib_err[:60]
                errors.append(
                    (
                        sys_char,
                        sys_name,
                        lat,
                        jd,
                        jd_label,
                        f"libephemeris error: {lib_err[:80]}",
                    )
                )
            elif not swe_ok and lib_ok:
                status = "SWE_ERROR"
                detail = swe_err[:60]
                errors.append(
                    (
                        sys_char,
                        sys_name,
                        lat,
                        jd,
                        jd_label,
                        f"pyswisseph error: {swe_err[:80]}",
                    )
                )
            elif not swe_ok and not lib_ok:
                status = "BOTH_ERROR"
                detail = f"swe: {swe_err[:30]} | lib: {lib_err[:30]}"
                errors.append(
                    (
                        sys_char,
                        sys_name,
                        lat,
                        jd,
                        jd_label,
                        f"both failed: swe={swe_err[:40]} lib={lib_err[:40]}",
                    )
                )

            print(
                f"  {lat:6.2f}  {jd:12.1f}  {max_cusp_diff:10.6f}°  "
                f"{asc_diff:10.6f}°  {mc_diff:10.6f}°  {status}"
                + (f"  {detail}" if detail else "")
            )

    print()

# ============================================================================
# Summary table: max diff per (system, latitude) across all dates
# ============================================================================

print("=" * 80)
print("SUMMARY: Maximum cusp difference per system per latitude (degrees)")
print("=" * 80)
print()

# Header
header = f"{'System':>6s}"
for lat in LATITUDES:
    header += f"  {lat:7.2f}°"
print(header)
print("-" * len(header))

for sys_char, sys_name in HOUSE_SYSTEMS:
    row = f"{sys_char:>6s}"
    for lat in LATITUDES:
        key = (sys_char, lat)
        if key in results:
            max_diff = max(results[key])
            if max_diff > THRESHOLD_WARN:
                row += f"  {max_diff:7.4f}*"
            else:
                row += f"  {max_diff:8.5f}"
        else:
            row += f"  {'ERR':>8s}"
    print(row + f"  ({sys_name})")

print()
print(f"  * = exceeds {THRESHOLD_WARN}° ({THRESHOLD_WARN * 3600:.0f} arcsec) threshold")
print()

# ============================================================================
# Flagged cases detail
# ============================================================================

if flagged:
    print("=" * 80)
    print(f"FLAGGED CASES (diff > {THRESHOLD_WARN}°):")
    print("=" * 80)
    for sys_char, sys_name, lat, jd, jd_label, detail in flagged:
        print(f"  {sys_name} ({sys_char}) | lat={lat:6.2f}° | {jd_label} | {detail}")
    print()

# ============================================================================
# Error cases
# ============================================================================

if errors:
    print("=" * 80)
    print(f"ERROR CASES ({len(errors)} total):")
    print("=" * 80)
    for sys_char, sys_name, lat, jd, jd_label, detail in errors:
        print(f"  {sys_name} ({sys_char}) | lat={lat:6.2f}° | {jd_label} | {detail}")
    print()

# ============================================================================
# Overall verdict
# ============================================================================

print("=" * 80)
print("OVERALL VERDICT")
print("=" * 80)

total_comparisons = sum(len(v) for v in results.values())
total_flagged = len(flagged)
total_errors = len(errors)

# Compute global max per system
for sys_char, sys_name in HOUSE_SYSTEMS:
    sys_diffs = []
    for lat in LATITUDES:
        key = (sys_char, lat)
        if key in results:
            sys_diffs.extend(results[key])
    if sys_diffs:
        gmax = max(sys_diffs)
        gmean = sum(sys_diffs) / len(sys_diffs)
        mark = " *** CONCERNING ***" if gmax > THRESHOLD_WARN else ""
        print(
            f"  {sys_name:>15s} ({sys_char}): "
            f"max={gmax:.6f}°  mean={gmean:.6f}°  "
            f"({gmax * 3600:.1f} arcsec max){mark}"
        )

print()
print(f"Total comparisons : {total_comparisons}")
print(f"Flagged (>{THRESHOLD_WARN}°)  : {total_flagged}")
print(f"Errors            : {total_errors}")

if total_flagged == 0 and total_errors == 0:
    print("\nAll house systems match within tolerance at all latitudes.")
elif total_flagged > 0:
    print(
        f"\n{total_flagged} case(s) exceed the {THRESHOLD_WARN}° threshold — review flagged cases above."
    )
