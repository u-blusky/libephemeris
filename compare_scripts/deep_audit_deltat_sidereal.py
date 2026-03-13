"""Deep precision audit: Delta-T and sidereal time calculations."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path("swisseph/ephe")

# =============================================================================
# Helpers
# =============================================================================


def year_to_jd(year: int) -> float:
    """Convert a calendar year to JD (Jan 1 noon).

    For negative years, use astronomical year numbering:
      500 BC  -> year -499
      1 BC    -> year 0
      1 AD    -> year 1
    """
    return swe.julday(year, 1, 1, 12.0)


def fmt_sec(val: float) -> str:
    """Format seconds with appropriate precision."""
    if abs(val) >= 100:
        return f"{val:>12.2f}"
    elif abs(val) >= 1:
        return f"{val:>12.4f}"
    else:
        return f"{val:>12.6f}"


# =============================================================================
# PART 1 — Delta-T audit
# =============================================================================

print("=" * 80)
print("PART 1: DELTA-T AUDIT")
print("=" * 80)
print()
print(
    f"{'Label':>20s}  {'Year':>8s}  {'JD':>14s}  "
    f"{'SWE (sec)':>12s}  {'LIB (sec)':>12s}  {'Diff (sec)':>12s}  {'FLAG':>6s}"
)
print("-" * 100)

deltat_dates = [
    ("Ancient: 500 BC", -499),
    ("Ancient: 0 AD", 0),
    ("Ancient: 500 AD", 500),
    ("Ancient: 1000 AD", 1000),
    ("Pre-tel: 1200", 1200),
    ("Pre-tel: 1400", 1400),
    ("Pre-tel: 1600", 1600),
    ("Historical: 1700", 1700),
    ("Historical: 1750", 1750),
    ("Historical: 1800", 1800),
    ("Historical: 1850", 1850),
    ("Modern: 1900", 1900),
    ("Modern: 1950", 1950),
    ("Modern: 1970", 1970),
    ("Modern: 1980", 1980),
    ("Modern: 1990", 1990),
    ("Modern: 2000", 2000),
    ("Modern: 2010", 2010),
    ("Modern: 2020", 2020),
    ("Modern: 2025", 2025),
    ("Future: 2050", 2050),
    ("Future: 2100", 2100),
    ("Future: 2200", 2200),
    ("Future: 2500", 2500),
]

deltat_flags = []
for label, year in deltat_dates:
    jd = year_to_jd(year)
    try:
        dt_swe = swe.deltat(jd) * 86400.0  # days -> seconds
        dt_lib = ephem.swe_deltat(jd) * 86400.0
        diff = dt_lib - dt_swe
        flag = "***" if abs(diff) > 1.0 else ""
        if flag:
            deltat_flags.append((label, year, diff))
        print(
            f"{label:>20s}  {year:>8d}  {jd:>14.2f}  "
            f"{fmt_sec(dt_swe)}  {fmt_sec(dt_lib)}  {fmt_sec(diff)}  {flag:>6s}"
        )
    except Exception as e:
        print(f"{label:>20s}  {year:>8d}  {jd:>14.2f}  ERROR: {e}")

print()
if deltat_flags:
    print(f"  *** {len(deltat_flags)} date(s) with Delta-T difference > 1 second:")
    for label, year, diff in deltat_flags:
        print(f"      {label}: diff = {diff:.4f} sec")
else:
    print("  All Delta-T differences within 1 second threshold.")


# =============================================================================
# PART 2 — Sidereal mode audit
# =============================================================================

print()
print("=" * 80)
print("PART 2: SIDEREAL POSITION AUDIT")
print("=" * 80)
print()

ayanamshas = [
    (0, "Fagan/Bradley"),
    (1, "Lahiri"),
    (3, "Raman"),
    (5, "Krishnamurti"),
]

bodies = [
    (0, "Sun"),
    (1, "Moon"),
]

test_jds = [
    (2451545.0, "J2000.0"),
    (2460400.0, "~2024-Apr"),
]

sidereal_flags = []

for sid_id, sid_name in ayanamshas:
    print(f"--- Ayanamsha: {sid_name} (mode={sid_id}) ---")
    arcsec_hdr = 'Diff(")'
    print(
        f"  {'Body':>6s}  {'Date':>12s}  {'JD':>14s}  "
        f"{'SWE lon(°)':>14s}  {'LIB lon(°)':>14s}  {'Diff(°)':>12s}  "
        f"{arcsec_hdr:>10s}  {'FLAG':>6s}"
    )
    print("  " + "-" * 94)

    for jd, date_label in test_jds:
        for body_id, body_name in bodies:
            try:
                # Set sidereal mode for both
                swe.set_sid_mode(sid_id)
                ephem.swe_set_sid_mode(sid_id)

                # Compute sidereal positions
                sp = swe.calc_ut(jd, body_id, swe.FLG_SIDEREAL)
                lp = ephem.swe_calc_ut(jd, body_id, ephem.FLG_SIDEREAL)

                swe_lon = (
                    float(sp[0][0])
                    if isinstance(sp[0], (list, tuple))
                    else float(sp[0])
                )
                lib_lon = (
                    float(lp[0][0])
                    if isinstance(lp[0], (list, tuple))
                    else float(lp[0])
                )

                diff_deg = lib_lon - swe_lon
                diff_arcsec = diff_deg * 3600.0
                flag = "***" if abs(diff_deg) > 0.001 else ""
                if flag:
                    sidereal_flags.append((sid_name, body_name, date_label, diff_deg))
                print(
                    f"  {body_name:>6s}  {date_label:>12s}  {jd:>14.1f}  "
                    f"{swe_lon:>14.6f}  {lib_lon:>14.6f}  {diff_deg:>12.6f}  "
                    f"{diff_arcsec:>10.2f}  {flag:>6s}"
                )
            except Exception as e:
                print(f"  {body_name:>6s}  {date_label:>12s}  {jd:>14.1f}  ERROR: {e}")

    print()

if sidereal_flags:
    print(
        f"  *** {len(sidereal_flags)} case(s) with sidereal longitude difference > 0.001°:"
    )
    for sid_name, body_name, date_label, diff in sidereal_flags:
        print(f"      {sid_name} / {body_name} @ {date_label}: diff = {diff:.6f}°")
else:
    print("  All sidereal longitude differences within 0.001° threshold.")


# =============================================================================
# PART 3 — Ayanamsha values audit
# =============================================================================

print()
print("=" * 80)
print("PART 3: AYANAMSHA VALUES AUDIT")
print("=" * 80)
print()

ayan_jds = [
    (year_to_jd(1900), "1900"),
    (year_to_jd(1950), "1950"),
    (2451545.0, "J2000.0"),
    (year_to_jd(2010), "2010"),
    (year_to_jd(2020), "2020"),
    (year_to_jd(2025), "2025"),
    (2460400.0, "~2024-Apr"),
    (year_to_jd(2050), "2050"),
]

ayan_flags = []

for sid_id, sid_name in ayanamshas:
    print(f"--- Ayanamsha: {sid_name} (mode={sid_id}) ---")
    arcsec_hdr = 'Diff(")'
    print(
        f"  {'Date':>12s}  {'JD':>14s}  "
        f"{'SWE (°)':>14s}  {'LIB (°)':>14s}  {'Diff(°)':>12s}  "
        f"{arcsec_hdr:>10s}  {'FLAG':>6s}"
    )
    print("  " + "-" * 80)

    for jd, date_label in ayan_jds:
        try:
            swe.set_sid_mode(sid_id)
            ephem.swe_set_sid_mode(sid_id)

            ayan_swe = float(swe.get_ayanamsa_ut(jd))
            ayan_lib = float(ephem.swe_get_ayanamsa_ut(jd))

            diff_deg = ayan_lib - ayan_swe
            diff_arcsec = diff_deg * 3600.0
            flag = "***" if abs(diff_deg) > 0.001 else ""
            if flag:
                ayan_flags.append((sid_name, date_label, diff_deg))
            print(
                f"  {date_label:>12s}  {jd:>14.2f}  "
                f"{ayan_swe:>14.6f}  {ayan_lib:>14.6f}  {diff_deg:>12.6f}  "
                f"{diff_arcsec:>10.2f}  {flag:>6s}"
            )
        except Exception as e:
            print(f"  {date_label:>12s}  {jd:>14.2f}  ERROR: {e}")

    print()

if ayan_flags:
    print(f"  *** {len(ayan_flags)} case(s) with ayanamsha difference > 0.001°:")
    for sid_name, date_label, diff in ayan_flags:
        print(f"      {sid_name} @ {date_label}: diff = {diff:.6f}°")
else:
    print("  All ayanamsha differences within 0.001° threshold.")


# =============================================================================
# SUMMARY
# =============================================================================

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
total_issues = len(deltat_flags) + len(sidereal_flags) + len(ayan_flags)
print(f"  Delta-T flags (>1 sec):            {len(deltat_flags)}")
print(f"  Sidereal longitude flags (>0.001°): {len(sidereal_flags)}")
print(f"  Ayanamsha flags (>0.001°):          {len(ayan_flags)}")
print(f"  Total flagged:                      {total_issues}")
if total_issues == 0:
    print("  RESULT: All checks passed within tolerances.")
else:
    print("  RESULT: Some checks exceeded tolerance thresholds — review above.")
print()
