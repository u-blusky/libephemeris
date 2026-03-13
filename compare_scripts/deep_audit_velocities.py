"""Deep precision audit: velocity/speed calculations vs pyswisseph."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import random
import statistics

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
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
    SEFLG_SWIEPH,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
    SEFLG_NONUT,
)

swe.set_ephe_path("swisseph/ephe")

random.seed(42)

# =============================================================================
# Configuration
# =============================================================================

BODIES = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
    (SE_URANUS, "Uranus"),
    (SE_NEPTUNE, "Neptune"),
    (SE_PLUTO, "Pluto"),
]

# JD range: 1900-01-01 to 2100-01-01 (approximate)
JD_1900 = 2415020.5
JD_2100 = 2488069.5

FLAG_COMBOS = [
    ("FLG_SPEED", SEFLG_SWIEPH | SEFLG_SPEED),
    ("FLG_SPEED|EQUATORIAL", SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL),
    ("FLG_SPEED|NONUT", SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_NONUT),
]

NUM_DATES = 20

# =============================================================================
# Helpers
# =============================================================================


def gen_random_jds(n: int) -> list[float]:
    """Generate n random JDs in 1900-2100."""
    return [random.uniform(JD_1900, JD_2100) for _ in range(n)]


def sign_str(v: float) -> str:
    return "+" if v >= 0 else "-"


# =============================================================================
# PART 1: Per-body speed precision at 20 random dates (FLG_SPEED)
# =============================================================================
print("=" * 100)
print("PART 1: VELOCITY PRECISION — FLG_SPEED — 10 bodies × 20 random dates")
print("=" * 100)
print()

jds = gen_random_jds(NUM_DATES)

# Collect all results for summary
all_issues: list[str] = []

flags = SEFLG_SWIEPH | SEFLG_SPEED

print(
    f"{'Body':>10s}  {'JD':>14s}  "
    f"{'Δ lon spd':>12s}  {'Δ lat spd':>12s}  {'Δ dist spd':>14s}  {'FLAG':>6s}"
)
print("-" * 90)

summary: dict[str, dict[str, list[float]]] = {}

for body_id, body_name in BODIES:
    summary[body_name] = {"lon_spd": [], "lat_spd": [], "dist_spd": []}
    for jd in jds:
        try:
            r_swe = swe.calc_ut(jd, body_id, flags)
            r_lib = ephem.swe_calc_ut(jd, body_id, flags)
            s, l = r_swe[0], r_lib[0]

            d_lon_spd = abs(s[3] - l[3])
            d_lat_spd = abs(s[4] - l[4])
            d_dist_spd = abs(s[5] - l[5])

            summary[body_name]["lon_spd"].append(d_lon_spd)
            summary[body_name]["lat_spd"].append(d_lat_spd)
            summary[body_name]["dist_spd"].append(d_dist_spd)

            flag = ""
            if d_lon_spd > 0.001:
                flag = "** LON"
                all_issues.append(
                    f"{body_name} JD={jd:.1f}: lon speed diff={d_lon_spd:.6f} °/day"
                )
            if d_lat_spd > 0.001:
                flag += " LAT" if flag else "** LAT"
            if d_dist_spd > 1e-6:
                flag += " DIST" if flag else "** DIST"

            if flag:
                print(
                    f"{body_name:>10s}  {jd:>14.2f}  "
                    f"{d_lon_spd:>12.8f}  {d_lat_spd:>12.8f}  {d_dist_spd:>14.10f}  {flag:>6s}"
                )
        except Exception as e:
            print(f"{body_name:>10s}  {jd:>14.2f}  ERROR: {e}")

print()
print("=" * 100)
print("PART 1 SUMMARY: Per-body max/mean speed differences (FLG_SPEED)")
print("=" * 100)
print()
print(
    f"{'Body':>10s}  "
    f"{'max Δlon':>12s}  {'mean Δlon':>12s}  "
    f"{'max Δlat':>12s}  {'mean Δlat':>12s}  "
    f"{'max Δdist':>14s}  {'mean Δdist':>14s}"
)
print("-" * 100)

for body_id, body_name in BODIES:
    d = summary[body_name]
    if d["lon_spd"]:
        print(
            f"{body_name:>10s}  "
            f"{max(d['lon_spd']):>12.8f}  {statistics.mean(d['lon_spd']):>12.8f}  "
            f"{max(d['lat_spd']):>12.8f}  {statistics.mean(d['lat_spd']):>12.8f}  "
            f"{max(d['dist_spd']):>14.10f}  {statistics.mean(d['dist_spd']):>14.10f}"
        )

# =============================================================================
# PART 2: Flag combinations (FLG_SPEED|EQUATORIAL, FLG_SPEED|NONUT)
# =============================================================================
print()
print("=" * 100)
print("PART 2: FLAG COMBINATIONS — speed differences per flag combo")
print("=" * 100)

for flag_name, flags_val in FLAG_COMBOS:
    print()
    print(f"--- {flag_name} (flags={flags_val}) ---")
    print(
        f"{'Body':>10s}  "
        f"{'max Δ[0]spd':>12s}  {'mean Δ[0]spd':>12s}  "
        f"{'max Δ[1]spd':>12s}  {'mean Δ[1]spd':>12s}  "
        f"{'max Δ[2]spd':>14s}  {'FLAG':>6s}"
    )
    print("-" * 90)

    for body_id, body_name in BODIES:
        errs_0: list[float] = []
        errs_1: list[float] = []
        errs_2: list[float] = []
        for jd in jds:
            try:
                r_swe = swe.calc_ut(jd, body_id, flags_val)
                r_lib = ephem.swe_calc_ut(jd, body_id, flags_val)
                s, l = r_swe[0], r_lib[0]
                errs_0.append(abs(s[3] - l[3]))
                errs_1.append(abs(s[4] - l[4]))
                errs_2.append(abs(s[5] - l[5]))
            except Exception:
                continue

        if errs_0:
            mx0 = max(errs_0)
            mx1 = max(errs_1)
            mx2 = max(errs_2)
            flag = ""
            if mx0 > 0.001:
                flag = "**"
                all_issues.append(
                    f"{body_name} [{flag_name}]: coord0 speed max diff={mx0:.6f}"
                )
            print(
                f"{body_name:>10s}  "
                f"{mx0:>12.8f}  {statistics.mean(errs_0):>12.8f}  "
                f"{mx1:>12.8f}  {statistics.mean(errs_1):>12.8f}  "
                f"{mx2:>14.10f}  {flag:>6s}"
            )

# =============================================================================
# PART 3: Retrograde precision — Mercury & Venus near stations
# =============================================================================
print()
print("=" * 100)
print("PART 3: RETROGRADE PRECISION — Mercury & Venus near station points")
print("=" * 100)
print()

# Mercury retrograde around JD 2460000 (~ early 2023)
# Venus retrograde around JD 2460500 (~ mid 2024)
# We scan a window around these dates to find where speed crosses zero.

RETRO_TESTS = [
    (SE_MERCURY, "Mercury", 2460000.0, 40.0),  # scan ±40 days
    (SE_VENUS, "Venus", 2460500.0, 60.0),  # scan ±60 days
]

flags = SEFLG_SWIEPH | SEFLG_SPEED

for body_id, body_name, jd_center, half_window in RETRO_TESTS:
    print(f"--- {body_name} near JD {jd_center:.0f} (±{half_window:.0f} days) ---")
    print(
        f"  {'JD':>14s}  "
        f"{'SWE lonSpd':>12s}  {'LIB lonSpd':>12s}  {'Δ lonSpd':>12s}  "
        f"{'SWE sign':>9s}  {'LIB sign':>9s}  {'SIGN_MISMATCH':>14s}"
    )
    print("  " + "-" * 95)

    # Sample at 0.5-day intervals
    sign_mismatches = 0
    max_diff = 0.0
    near_station_diffs: list[tuple[float, float, float, float]] = []

    steps = int(2 * half_window / 0.5) + 1
    for i in range(steps):
        jd = jd_center - half_window + i * 0.5
        try:
            r_swe = swe.calc_ut(jd, body_id, flags)
            r_lib = ephem.swe_calc_ut(jd, body_id, flags)
            s_spd = r_swe[0][3]
            l_spd = r_lib[0][3]
            diff = abs(s_spd - l_spd)
            max_diff = max(max_diff, diff)

            swe_sign = "retro" if s_spd < 0 else "direct"
            lib_sign = "retro" if l_spd < 0 else "direct"
            mismatch = "*** MISMATCH" if swe_sign != lib_sign else ""

            if mismatch:
                sign_mismatches += 1
                all_issues.append(
                    f"{body_name} JD={jd:.2f}: SIGN MISMATCH swe={s_spd:.6f} lib={l_spd:.6f}"
                )

            # Track near-station points (speed close to zero)
            if abs(s_spd) < 0.1 or abs(l_spd) < 0.1:
                near_station_diffs.append((jd, s_spd, l_spd, diff))

            # Only print interesting rows: near station or large diff
            if abs(s_spd) < 0.15 or abs(l_spd) < 0.15 or diff > 0.001 or mismatch:
                print(
                    f"  {jd:>14.2f}  "
                    f"{s_spd:>12.6f}  {l_spd:>12.6f}  {diff:>12.8f}  "
                    f"{swe_sign:>9s}  {lib_sign:>9s}  {mismatch:>14s}"
                )
        except Exception as e:
            print(f"  {jd:>14.2f}  ERROR: {e}")

    print()
    print(f"  {body_name} retrograde summary:")
    print(f"    Max lon speed diff in window: {max_diff:.8f} °/day")
    print(f"    Sign mismatches (retro vs direct): {sign_mismatches}")
    if near_station_diffs:
        station_max = max(d[3] for d in near_station_diffs)
        station_mean = statistics.mean(d[3] for d in near_station_diffs)
        print(
            f"    Near-station points (|speed| < 0.1): {len(near_station_diffs)} samples"
        )
        print(f"    Near-station max diff: {station_max:.8f} °/day")
        print(f"    Near-station mean diff: {station_mean:.8f} °/day")
    print()

# =============================================================================
# PART 4: Fine-grained station scan (0.01-day steps around station)
# =============================================================================
print("=" * 100)
print("PART 4: FINE-GRAINED STATION SCAN — 0.01-day steps")
print("=" * 100)
print()

# Find approximate station points from part 3, then zoom in
for body_id, body_name, jd_center, half_window in RETRO_TESTS:
    # Coarse scan to find sign changes
    prev_swe_spd = None
    station_jds: list[float] = []

    for i in range(int(2 * half_window)):
        jd = jd_center - half_window + i
        try:
            r_swe = swe.calc_ut(jd, body_id, flags)
            cur_spd = r_swe[0][3]
            if prev_swe_spd is not None and prev_swe_spd * cur_spd < 0:
                station_jds.append(jd - 0.5)
            prev_swe_spd = cur_spd
        except Exception:
            continue

    for station_jd in station_jds[:2]:  # up to 2 stations
        print(f"--- {body_name} station near JD {station_jd:.1f} ---")
        print(
            f"  {'JD':>14s}  "
            f"{'SWE lonSpd':>12s}  {'LIB lonSpd':>12s}  {'Δ':>12s}  "
            f"{'SWE sign':>9s}  {'LIB sign':>9s}  {'FLAG':>14s}"
        )
        print("  " + "-" * 85)

        for j in range(-50, 51):
            jd = station_jd + j * 0.01
            try:
                r_swe = swe.calc_ut(jd, body_id, flags)
                r_lib = ephem.swe_calc_ut(jd, body_id, flags)
                s_spd = r_swe[0][3]
                l_spd = r_lib[0][3]
                diff = abs(s_spd - l_spd)

                swe_sign = "retro" if s_spd < 0 else "direct"
                lib_sign = "retro" if l_spd < 0 else "direct"
                mismatch = "*** MISMATCH" if swe_sign != lib_sign else ""

                if abs(s_spd) < 0.05 or abs(l_spd) < 0.05 or mismatch:
                    print(
                        f"  {jd:>14.4f}  "
                        f"{s_spd:>12.8f}  {l_spd:>12.8f}  {diff:>12.8f}  "
                        f"{swe_sign:>9s}  {lib_sign:>9s}  {mismatch:>14s}"
                    )
            except Exception:
                continue
        print()

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("=" * 100)
print("FINAL SUMMARY")
print("=" * 100)
print()

# Overall per-body summary table
print("Per-body max lon speed difference (FLG_SPEED, °/day):")
print(
    f"{'Body':>10s}  {'Max Δ lon spd':>14s}  {'Mean Δ lon spd':>14s}  {'Status':>10s}"
)
print("-" * 55)
for body_id, body_name in BODIES:
    d = summary[body_name]
    if d["lon_spd"]:
        mx = max(d["lon_spd"])
        mn = statistics.mean(d["lon_spd"])
        status = "CONCERN" if mx > 0.001 else "OK"
        print(f"{body_name:>10s}  {mx:>14.8f}  {mn:>14.8f}  {status:>10s}")

print()
if all_issues:
    print(f"FLAGGED ISSUES ({len(all_issues)}):")
    for issue in all_issues:
        print(f"  !! {issue}")
else:
    print("No issues flagged. All velocity differences within tolerance.")

print()
print("Done.")
