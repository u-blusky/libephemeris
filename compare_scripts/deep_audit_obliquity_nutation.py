"""Deep precision audit: Obliquity, nutation, and precession vs pyswisseph."""

from __future__ import annotations

import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path("swisseph/ephe")

# =============================================================================
# Helpers
# =============================================================================


def year_to_jd(year: int) -> float:
    """Convert a calendar year to JD (Jan 1 noon)."""
    return swe.julday(year, 1, 1, 12.0)


def deg_to_arcsec(deg: float) -> float:
    """Convert degrees to arcseconds."""
    return deg * 3600.0


def fmt_deg(val: float) -> str:
    """Format degrees with high precision."""
    return f"{val:>16.10f}"


def fmt_arcsec(val: float) -> str:
    """Format arcseconds with appropriate precision."""
    return f"{val:>12.6f}"


def flag_diff(diff_arcsec: float, threshold: float = 0.01) -> str:
    """Return flag marker if difference exceeds threshold."""
    if abs(diff_arcsec) > threshold:
        return "  *** FAIL"
    return ""


# =============================================================================
# Discover constant values
# =============================================================================

print("=" * 90)
print("CONSTANT VALUES")
print("=" * 90)
print()

lib_ecl_nut = ephem.SE_ECL_NUT
swe_ecl_nut = swe.ECL_NUT
print(f"  ephem.SE_ECL_NUT = {lib_ecl_nut}")
print(f"  swe.ECL_NUT      = {swe_ecl_nut}")

lib_sun = ephem.SE_SUN
swe_sun = swe.SUN
print(f"  ephem.SE_SUN     = {lib_sun}")
print(f"  swe.SUN          = {swe_sun}")

lib_nonut = ephem.SEFLG_NONUT
swe_nonut = swe.FLG_NONUT
print(f"  ephem.SEFLG_NONUT = {lib_nonut}")
print(f"  swe.FLG_NONUT     = {swe_nonut}")
print()

# Track overall stats
all_diffs: dict[str, list[float]] = {
    "true_obl": [],
    "mean_obl": [],
    "nut_lon": [],
    "nut_obl": [],
    "derived_nut": [],
    "precession": [],
}
fail_count = 0

# =============================================================================
# PART 1 — OBLIQUITY TEST: True & mean obliquity + nutation at 20 dates
# =============================================================================

print("=" * 90)
print("PART 1: OBLIQUITY & NUTATION (ECL_NUT) — 20 dates spanning 1800-2200")
print("=" * 90)
print()

dates_obliquity = [
    ("1800", 1800),
    ("1820", 1820),
    ("1840", 1840),
    ("1860", 1860),
    ("1880", 1880),
    ("1900", 1900),
    ("1920", 1920),
    ("1940", 1940),
    ("1960", 1960),
    ("1980", 1980),
    ("1990", 1990),
    ("2000 (J2000)", 2000),
    ("2010", 2010),
    ("2020", 2020),
    ("2030", 2030),
    ("2050", 2050),
    ("2080", 2080),
    ("2100", 2100),
    ("2150", 2150),
    ("2200", 2200),
]

labels_ecl = [
    "True Obliquity",
    "Mean Obliquity",
    "Nutation Lon (dpsi)",
    "Nutation Obl (deps)",
]

for idx, (label, year) in enumerate(dates_obliquity):
    jd = year_to_jd(year)

    # libephemeris
    lib_pos, lib_flag = ephem.swe_calc_ut(jd, lib_ecl_nut, 0)
    # pyswisseph
    swe_pos, swe_flag = swe.calc_ut(jd, swe_ecl_nut, 0)

    print(f"--- {label}  (JD {jd:.1f}) ---")
    print(
        f"  {'Component':<22s}  {'SWE (deg)':>16s}  {'LIB (deg)':>16s}  "
        f"{'Diff (arcsec)':>12s}  FLAG"
    )

    for i, comp_label in enumerate(labels_ecl):
        s_val = swe_pos[i]
        l_val = lib_pos[i]
        diff_deg = l_val - s_val
        diff_as = deg_to_arcsec(diff_deg)

        key = ["true_obl", "mean_obl", "nut_lon", "nut_obl"][i]
        all_diffs[key].append(abs(diff_as))

        flag = flag_diff(diff_as)
        if flag:
            fail_count += 1

        print(
            f"  {comp_label:<22s}  {fmt_deg(s_val)}  {fmt_deg(l_val)}  "
            f"{fmt_arcsec(diff_as)}{flag}"
        )

    # Also verify index 4 and 5 are zero in both
    for z in (4, 5):
        if swe_pos[z] != 0.0 or lib_pos[z] != 0.0:
            print(f"  WARNING: index {z} not zero: swe={swe_pos[z]}, lib={lib_pos[z]}")

    print()


# Summary for Part 1
print("-" * 90)
print("PART 1 SUMMARY:")
for key, label in [
    ("true_obl", "True Obliquity"),
    ("mean_obl", "Mean Obliquity"),
    ("nut_lon", "Nutation Longitude"),
    ("nut_obl", "Nutation Obliquity"),
]:
    diffs = all_diffs[key]
    if diffs:
        print(
            f'  {label:<22s}  max={max(diffs):>10.6f}"  '
            f'mean={sum(diffs) / len(diffs):>10.6f}"  '
            f'min={min(diffs):>10.6f}"'
        )
print()

# =============================================================================
# PART 2 — NUTATION IMPACT TEST: Derived nutation from position differences
# =============================================================================

print("=" * 90)
print("PART 2: NUTATION IMPACT — Derived nutation vs ECL_NUT nutation")
print("=" * 90)
print()
print("Compute Sun longitude WITH and WITHOUT nutation. The difference is the")
print("nutation in longitude (dpsi). Compare this derived value to the ECL_NUT")
print("nutation_longitude value from each library.")
print()

# Use 10 dates for this test
dates_nut_impact = [
    ("1850", 1850),
    ("1900", 1900),
    ("1950", 1950),
    ("1980", 1980),
    ("2000", 2000),
    ("2010", 2010),
    ("2020", 2020),
    ("2030", 2030),
    ("2050", 2050),
    ("2100", 2100),
]

hdr_dpsi = 'Derived dpsi "'
hdr_ecl = 'ECL_NUT dpsi "'
hdr_diff = 'Diff "'
print(
    f"  {'Date':<10s}  {'Lib':<5s}  {'Sun(nut) deg':>16s}  {'Sun(nonut) deg':>16s}  "
    f"{hdr_dpsi:>14s}  {hdr_ecl:>14s}  {hdr_diff:>10s}  FLAG"
)
print("-" * 110)

for label, year in dates_nut_impact:
    jd = year_to_jd(year)

    for lib_name, calc_fn, ecl_nut_id, sun_id, nonut_flag in [
        ("SWE", swe.calc_ut, swe_ecl_nut, swe_sun, swe_nonut),
        ("LIB", ephem.swe_calc_ut, lib_ecl_nut, lib_sun, lib_nonut),
    ]:
        # Sun with nutation (default)
        pos_nut, _ = calc_fn(jd, sun_id, 0)
        # Sun without nutation
        pos_nonut, _ = calc_fn(jd, sun_id, nonut_flag)
        # Derived nutation in longitude
        derived_dpsi_deg = pos_nut[0] - pos_nonut[0]
        derived_dpsi_as = deg_to_arcsec(derived_dpsi_deg)

        # ECL_NUT nutation in longitude
        ecl_pos, _ = calc_fn(jd, ecl_nut_id, 0)
        ecl_dpsi_as = deg_to_arcsec(ecl_pos[2])  # index 2 = nutation longitude

        diff_as = derived_dpsi_as - ecl_dpsi_as

        if lib_name == "LIB":
            all_diffs["derived_nut"].append(abs(diff_as))

        flag = flag_diff(diff_as, threshold=0.01)
        if flag:
            fail_count += 1

        print(
            f"  {label:<10s}  {lib_name:<5s}  {fmt_deg(pos_nut[0])}  "
            f"{fmt_deg(pos_nonut[0])}  {fmt_arcsec(derived_dpsi_as)}  "
            f"{fmt_arcsec(ecl_dpsi_as)}  {fmt_arcsec(diff_as)}{flag}"
        )

    print()

# Cross-library comparison: does the derived nutation match between SWE and LIB?
print("-" * 110)
print("Cross-library derived nutation comparison (SWE vs LIB):")
hdr_swe_d = 'SWE derived "'
hdr_lib_d = 'LIB derived "'
hdr_diff2 = 'Diff "'
print(f"  {'Date':<10s}  {hdr_swe_d:>14s}  {hdr_lib_d:>14s}  {hdr_diff2:>10s}  FLAG")
print("-" * 80)

for label, year in dates_nut_impact:
    jd = year_to_jd(year)

    # SWE
    swe_nut, _ = swe.calc_ut(jd, swe_sun, 0)
    swe_nonut_pos, _ = swe.calc_ut(jd, swe_sun, swe_nonut)
    swe_derived = deg_to_arcsec(swe_nut[0] - swe_nonut_pos[0])

    # LIB
    lib_nut, _ = ephem.swe_calc_ut(jd, lib_sun, 0)
    lib_nonut_pos, _ = ephem.swe_calc_ut(jd, lib_sun, lib_nonut)
    lib_derived = deg_to_arcsec(lib_nut[0] - lib_nonut_pos[0])

    diff = lib_derived - swe_derived
    flag = flag_diff(diff, threshold=0.01)
    if flag:
        fail_count += 1

    print(
        f"  {label:<10s}  {fmt_arcsec(swe_derived)}  {fmt_arcsec(lib_derived)}  "
        f"{fmt_arcsec(diff)}{flag}"
    )

print()

# Summary for Part 2
print("PART 2 SUMMARY:")
diffs = all_diffs["derived_nut"]
if diffs:
    print(
        f'  Derived vs ECL_NUT (LIB)  max={max(diffs):>10.6f}"  '
        f'mean={sum(diffs) / len(diffs):>10.6f}"  '
        f'min={min(diffs):>10.6f}"'
    )
print()

# =============================================================================
# PART 3 — PRECESSION TEST: Sun longitude at dates far from J2000
# =============================================================================

print("=" * 90)
print("PART 3: PRECESSION — Sun longitude at dates far from J2000")
print("=" * 90)
print()
print("Precession shifts the ecliptic frame. Both libraries should agree closely.")
print()

dates_precession = [
    ("1800", 1800),
    ("1850", 1850),
    ("1900", 1900),
    ("1950", 1950),
    ("J2000", 2000),
    ("2010", 2010),
    ("2020", 2020),
    ("2050", 2050),
    ("2100", 2100),
    ("2150", 2150),
    ("2200", 2200),
]

print(
    f"  {'Date':<10s}  {'SWE Sun lon (deg)':>18s}  {'LIB Sun lon (deg)':>18s}  "
    f"{'Diff (arcsec)':>14s}  FLAG"
)
print("-" * 90)

for label, year in dates_precession:
    jd = year_to_jd(year)

    swe_pos, _ = swe.calc_ut(jd, swe_sun, 0)
    lib_pos, _ = ephem.swe_calc_ut(jd, lib_sun, 0)

    swe_lon = swe_pos[0]
    lib_lon = lib_pos[0]
    diff_as = deg_to_arcsec(lib_lon - swe_lon)

    all_diffs["precession"].append(abs(diff_as))
    flag = flag_diff(diff_as, threshold=0.01)
    if flag:
        fail_count += 1

    print(
        f"  {label:<10s}  {swe_lon:>18.10f}  {lib_lon:>18.10f}  "
        f"{fmt_arcsec(diff_as)}{flag}"
    )

print()

# Show that precession has accumulated over time
print("Precession accumulation check (longitude shift from J2000):")
jd_j2000 = year_to_jd(2000)
swe_j2000, _ = swe.calc_ut(jd_j2000, swe_sun, 0)
lib_j2000, _ = ephem.swe_calc_ut(jd_j2000, lib_sun, 0)

print(
    f"  {'Date':<10s}  {'SWE shift (deg)':>16s}  {'LIB shift (deg)':>16s}  "
    f"{'Shift diff (arcsec)':>20s}"
)
print("-" * 80)

for label, year in [("1900", 1900), ("2100", 2100)]:
    jd = year_to_jd(year)
    swe_pos, _ = swe.calc_ut(jd, swe_sun, 0)
    lib_pos, _ = ephem.swe_calc_ut(jd, lib_sun, 0)

    swe_shift = swe_pos[0] - swe_j2000[0]
    lib_shift = lib_pos[0] - lib_j2000[0]
    shift_diff = deg_to_arcsec(lib_shift - swe_shift)

    print(
        f"  {label:<10s}  {swe_shift:>16.8f}  {lib_shift:>16.8f}  "
        f"{fmt_arcsec(shift_diff)}"
    )

print()

# Summary for Part 3
print("PART 3 SUMMARY:")
diffs = all_diffs["precession"]
if diffs:
    print(
        f'  Sun longitude diff      max={max(diffs):>10.6f}"  '
        f'mean={sum(diffs) / len(diffs):>10.6f}"  '
        f'min={min(diffs):>10.6f}"'
    )
print()

# =============================================================================
# OVERALL SUMMARY
# =============================================================================

print("=" * 90)
print("OVERALL SUMMARY")
print("=" * 90)
print()
print(f'  Total FAIL flags (>0.01" threshold): {fail_count}')
print()

for key, label in [
    ("true_obl", "True Obliquity (ECL_NUT)"),
    ("mean_obl", "Mean Obliquity (ECL_NUT)"),
    ("nut_lon", "Nutation Longitude (ECL_NUT)"),
    ("nut_obl", "Nutation Obliquity (ECL_NUT)"),
    ("derived_nut", "Derived vs ECL_NUT Nutation (LIB)"),
    ("precession", "Sun Longitude (precession)"),
]:
    diffs = all_diffs[key]
    if diffs:
        max_d = max(diffs)
        status = "PASS" if max_d <= 0.01 else "FAIL"
        print(
            f'  [{status:>4s}]  {label:<40s}  max={max_d:>10.6f}"  '
            f'mean={sum(diffs) / len(diffs):>10.6f}"'
        )

print()

if fail_count == 0:
    print('ALL TESTS PASSED — all differences within 0.01" (sub-arcsecond precision)')
    sys.exit(0)
else:
    print(f'WARNING: {fail_count} value(s) exceeded 0.01" threshold')
    sys.exit(1)
