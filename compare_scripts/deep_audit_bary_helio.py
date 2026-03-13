"""
Deep audit: libephemeris vs pyswisseph for barycentric, heliocentric,
equatorial, and XYZ coordinate modes.

Tests multiple flag combinations across 10 dates spanning 1950-2050.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SEFLG_SWIEPH,
    SEFLG_SPEED,
    SEFLG_HELCTR,
    SEFLG_BARYCTR,
    SEFLG_EQUATORIAL,
    SEFLG_XYZ,
)

swe.set_ephe_path("swisseph/ephe")

# ============================================================================
# CONSTANTS
# ============================================================================

BODY_NAMES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}

# 10 dates spanning 1950-2050 (JD values for Jan 1 12:00 TT of each year)
TEST_YEARS = [1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2040, 2050]


def year_to_jd(year: int) -> float:
    """Convert year to approximate JD (Jan 1, 12:00 UT)."""
    # J2000.0 = 2451545.0 (2000-Jan-1.5)
    return 2451545.0 + (year - 2000) * 365.25


TEST_JDS = [(y, year_to_jd(y)) for y in TEST_YEARS]

# Value labels per mode
LABELS_SPHERICAL = ["lon", "lat", "dist", "dlon/dt", "dlat/dt", "ddist/dt"]
LABELS_EQUATORIAL = ["RA", "Dec", "dist", "dRA/dt", "dDec/dt", "ddist/dt"]
LABELS_XYZ = ["X", "Y", "Z", "dX/dt", "dY/dt", "dZ/dt"]


# ============================================================================
# HELPERS
# ============================================================================


def angular_diff(a: float, b: float) -> float:
    """Angular difference accounting for 360-degree wrap."""
    d = abs(a - b) % 360
    return min(d, 360 - d)


def arcsec(deg: float) -> float:
    return deg * 3600.0


@dataclass
class CompareResult:
    mode: str
    year: int
    jd: float
    body_id: int
    body_name: str
    lib_vals: tuple | None = None
    swe_vals: tuple | None = None
    diffs: list[float] = field(default_factory=list)
    error: str | None = None
    flagged: bool = False
    flag_reasons: list[str] = field(default_factory=list)


def compare_one(
    mode: str,
    year: int,
    jd: float,
    body_id: int,
    body_name: str,
    flags_lib: int,
    flags_swe: int,
    labels: list[str],
    thresholds: list[float],
    use_angular_diff: list[bool],
) -> CompareResult:
    """Compare libephemeris vs pyswisseph for one body/date/flag combo."""
    res = CompareResult(
        mode=mode, year=year, jd=jd, body_id=body_id, body_name=body_name
    )

    # pyswisseph
    try:
        r_swe = swe.calc_ut(jd, body_id, flags_swe)
        res.swe_vals = tuple(r_swe[0]) if hasattr(r_swe[0], "__iter__") else r_swe[:6]
    except Exception as e:
        res.error = f"swe error: {e}"
        return res

    # libephemeris
    try:
        r_lib = ephem.swe_calc_ut(jd, body_id, flags_lib)
        res.lib_vals = tuple(r_lib[0])
    except Exception as e:
        res.error = f"lib error: {e}"
        return res

    # Compute diffs for all 6 values
    for i in range(6):
        sv = res.swe_vals[i]
        lv = res.lib_vals[i]
        if use_angular_diff[i]:
            d = angular_diff(sv, lv)
        else:
            d = abs(sv - lv)
        res.diffs.append(d)

        if d > thresholds[i]:
            res.flagged = True
            res.flag_reasons.append(f"{labels[i]}: {d:.8f} > {thresholds[i]}")

    return res


def print_mode_summary(mode: str, results: list[CompareResult], labels: list[str]):
    """Print summary table for a mode."""
    ok_results = [r for r in results if r.error is None and len(r.diffs) == 6]
    err_results = [r for r in results if r.error is not None]
    flagged_results = [r for r in results if r.flagged]

    print(f"\n{'=' * 100}")
    print(f"  MODE: {mode}")
    print(f"{'=' * 100}")
    print(
        f"  Total tests: {len(results)}  |  OK: {len(ok_results)}  |  "
        f"Errors: {len(err_results)}  |  Flagged: {len(flagged_results)}"
    )

    if err_results:
        print(f"\n  --- Errors ({len(err_results)}) ---")
        for r in err_results:
            print(f"    {r.body_name:10s} y={r.year}  {r.error}")

    if not ok_results:
        print("  No successful comparisons.")
        return

    # Aggregate stats per value
    print(
        f"\n  {'Value':<12s} {'Max Diff':>14s} {'Mean Diff':>14s} {'Median Diff':>14s}"
    )
    print(f"  {'-' * 12} {'-' * 14} {'-' * 14} {'-' * 14}")
    for i, label in enumerate(labels):
        dvals = sorted([r.diffs[i] for r in ok_results])
        mx = max(dvals)
        mn = sum(dvals) / len(dvals)
        mid = dvals[len(dvals) // 2]
        print(f"  {label:<12s} {mx:>14.8f} {mn:>14.8f} {mid:>14.8f}")

    # Per-body max diff for first 3 values (position)
    bodies_seen = sorted(set(r.body_id for r in ok_results))
    print(f"\n  Per-body max position diffs:")
    print(
        f"  {'Body':<10s} {'Max ' + labels[0]:>16s} {'Max ' + labels[1]:>16s} {'Max ' + labels[2]:>16s}"
    )
    print(f"  {'-' * 10} {'-' * 16} {'-' * 16} {'-' * 16}")
    for bid in bodies_seen:
        br = [r for r in ok_results if r.body_id == bid]
        bname = br[0].body_name
        m0 = max(r.diffs[0] for r in br)
        m1 = max(r.diffs[1] for r in br)
        m2 = max(r.diffs[2] for r in br)
        print(f"  {bname:<10s} {m0:>16.8f} {m1:>16.8f} {m2:>16.8f}")

    if flagged_results:
        print(f"\n  --- Flagged results ({len(flagged_results)}) ---")
        for r in flagged_results:
            reasons = "; ".join(r.flag_reasons)
            print(f"    {r.body_name:10s} y={r.year}  JD={r.jd:.1f}  {reasons}")
    else:
        print(f"\n  No flagged results — all within thresholds.")


# ============================================================================
# MODE 1: HELIOCENTRIC (bodies 2-9, Sun/Moon excluded)
# ============================================================================

print("=" * 100)
print("DEEP AUDIT: Barycentric, Heliocentric, Equatorial, XYZ modes")
print("libephemeris vs pyswisseph")
print("=" * 100)

helio_results: list[CompareResult] = []
helio_flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR
# Thresholds: lon>0.001°, lat>0.001°, dist>0.0001 AU, speed thresholds relaxed
HELIO_THRESHOLDS = [0.001, 0.001, 0.0001, 0.01, 0.01, 0.001]
HELIO_ANGULAR = [True, False, False, False, False, False]

print("\nRunning HELIOCENTRIC tests (bodies 2-9, 10 dates)...")
for year, jd in TEST_JDS:
    for body_id in range(2, 10):  # Mercury through Pluto
        r = compare_one(
            mode="HELIOCENTRIC",
            year=year,
            jd=jd,
            body_id=body_id,
            body_name=BODY_NAMES[body_id],
            flags_lib=helio_flags,
            flags_swe=helio_flags,
            labels=LABELS_SPHERICAL,
            thresholds=HELIO_THRESHOLDS,
            use_angular_diff=HELIO_ANGULAR,
        )
        helio_results.append(r)
        if r.flagged:
            print(f"  !! FLAGGED: {r.body_name} y={year}: {'; '.join(r.flag_reasons)}")

print_mode_summary("HELIOCENTRIC (SEFLG_HELCTR=8)", helio_results, LABELS_SPHERICAL)

# ============================================================================
# MODE 2: BARYCENTRIC (bodies 0-9)
# ============================================================================

bary_results: list[CompareResult] = []
bary_flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_BARYCTR
BARY_THRESHOLDS = [0.001, 0.001, 0.0001, 0.01, 0.01, 0.001]
BARY_ANGULAR = [True, False, False, False, False, False]

print("\nRunning BARYCENTRIC tests (bodies 0-9, 10 dates)...")
for year, jd in TEST_JDS:
    for body_id in range(0, 10):  # Sun through Pluto
        r = compare_one(
            mode="BARYCENTRIC",
            year=year,
            jd=jd,
            body_id=body_id,
            body_name=BODY_NAMES[body_id],
            flags_lib=bary_flags,
            flags_swe=bary_flags,
            labels=LABELS_SPHERICAL,
            thresholds=BARY_THRESHOLDS,
            use_angular_diff=BARY_ANGULAR,
        )
        bary_results.append(r)
        if r.flagged:
            print(f"  !! FLAGGED: {r.body_name} y={year}: {'; '.join(r.flag_reasons)}")

print_mode_summary("BARYCENTRIC (SEFLG_BARYCTR=16384)", bary_results, LABELS_SPHERICAL)

# ============================================================================
# MODE 3: EQUATORIAL (bodies 0,1,4,5 — Sun, Moon, Mars, Jupiter)
# ============================================================================

equat_results: list[CompareResult] = []
equat_flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL
# RA is in degrees (0-360) for both libs, so angular diff applies
EQUAT_THRESHOLDS = [0.001, 0.001, 0.0001, 0.01, 0.01, 0.001]
EQUAT_ANGULAR = [True, False, False, False, False, False]

print("\nRunning EQUATORIAL tests (Sun, Moon, Mars, Jupiter, 10 dates)...")
equat_bodies = [0, 1, 4, 5]
for year, jd in TEST_JDS:
    for body_id in equat_bodies:
        r = compare_one(
            mode="EQUATORIAL",
            year=year,
            jd=jd,
            body_id=body_id,
            body_name=BODY_NAMES[body_id],
            flags_lib=equat_flags,
            flags_swe=equat_flags,
            labels=LABELS_EQUATORIAL,
            thresholds=EQUAT_THRESHOLDS,
            use_angular_diff=EQUAT_ANGULAR,
        )
        equat_results.append(r)
        if r.flagged:
            print(f"  !! FLAGGED: {r.body_name} y={year}: {'; '.join(r.flag_reasons)}")

print_mode_summary(
    "EQUATORIAL (SEFLG_EQUATORIAL=2048)", equat_results, LABELS_EQUATORIAL
)

# ============================================================================
# MODE 4: XYZ CARTESIAN (bodies 0-9)
# ============================================================================

xyz_results: list[CompareResult] = []
xyz_flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_XYZ
# All values are linear (AU or AU/day), no angular wrapping
XYZ_THRESHOLDS = [0.00001, 0.00001, 0.00001, 0.001, 0.001, 0.001]
XYZ_ANGULAR = [False, False, False, False, False, False]

print("\nRunning XYZ CARTESIAN tests (bodies 0-9, 10 dates)...")
for year, jd in TEST_JDS:
    for body_id in range(0, 10):
        r = compare_one(
            mode="XYZ",
            year=year,
            jd=jd,
            body_id=body_id,
            body_name=BODY_NAMES[body_id],
            flags_lib=xyz_flags,
            flags_swe=xyz_flags,
            labels=LABELS_XYZ,
            thresholds=XYZ_THRESHOLDS,
            use_angular_diff=XYZ_ANGULAR,
        )
        xyz_results.append(r)
        if r.flagged:
            print(f"  !! FLAGGED: {r.body_name} y={year}: {'; '.join(r.flag_reasons)}")

print_mode_summary("XYZ CARTESIAN (SEFLG_XYZ=4096)", xyz_results, LABELS_XYZ)

# ============================================================================
# GRAND SUMMARY
# ============================================================================

all_modes = [
    ("HELIOCENTRIC", helio_results),
    ("BARYCENTRIC", bary_results),
    ("EQUATORIAL", equat_results),
    ("XYZ CARTESIAN", xyz_results),
]

print(f"\n{'=' * 100}")
print("  GRAND SUMMARY")
print(f"{'=' * 100}")
print(
    f"  {'Mode':<20s} {'Tests':>6s} {'OK':>6s} {'Errors':>7s} {'Flagged':>8s} {'Max pos diff':>16s}"
)
print(f"  {'-' * 20} {'-' * 6} {'-' * 6} {'-' * 7} {'-' * 8} {'-' * 16}")

any_concern = False
for mode_name, mode_results in all_modes:
    total = len(mode_results)
    ok = [r for r in mode_results if r.error is None and len(r.diffs) == 6]
    errs = [r for r in mode_results if r.error is not None]
    flagged = [r for r in mode_results if r.flagged]
    if ok:
        max_pos = max(max(r.diffs[0], r.diffs[1], r.diffs[2]) for r in ok)
    else:
        max_pos = float("nan")
    status = "CONCERN" if flagged else "OK"
    if flagged:
        any_concern = True
    print(
        f"  {mode_name:<20s} {total:>6d} {len(ok):>6d} {len(errs):>7d} "
        f"{len(flagged):>8d} {max_pos:>16.8f}  {status}"
    )

print()
if any_concern:
    print("  ** CONCERNS FOUND — see flagged results above for details. **")
else:
    print("  All modes PASSED within thresholds.")
print()
