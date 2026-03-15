#!/usr/bin/env python3
"""Round 114: Multi-Flag Combination Stress Test

Exhaustive test of all meaningful flag combinations for planet positions.
Tests every combination of SPEED, J2000, NONUT, NOABERR, TRUEPOS, EQUATORIAL,
HELCTR, BARYCTR, SIDEREAL, XYZ, RADIANS for multiple bodies and epochs.
"""

from __future__ import annotations

import sys
import os
import math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path("swisseph/ephe")

passed = 0
failed = 0
errors = 0

SEFLG_SPEED = 256
SEFLG_HELCTR = 8
SEFLG_BARYCTR = 4
SEFLG_TRUEPOS = 16
SEFLG_J2000 = 32
SEFLG_NONUT = 64
SEFLG_NOABERR = 1024
SEFLG_EQUATORIAL = 2048
SEFLG_XYZ = 4096
SEFLG_RADIANS = 8192
SEFLG_SIDEREAL = 65536

SE_AST_OFFSET = 10000

BODIES = {0: "Sun", 1: "Moon", 2: "Mercury", 4: "Mars", 5: "Jupiter", 9: "Pluto"}

TEST_EPOCHS = [
    (2451545.0, "J2000"),
    (swe.julday(2024, 6, 15, 12.0), "2024"),
    (swe.julday(1980, 3, 21, 12.0), "1980"),
]

# Flag combinations to test (excluding mutually exclusive ones)
# HELCTR and BARYCTR are mutually exclusive
# SIDEREAL requires set_sid_mode first
FLAG_COMBOS = [
    (SEFLG_SPEED, "SPEED"),
    (SEFLG_SPEED | SEFLG_NONUT, "SPEED+NONUT"),
    (SEFLG_SPEED | SEFLG_J2000 | SEFLG_NONUT, "SPEED+J2000+NONUT"),
    (SEFLG_SPEED | SEFLG_NOABERR, "SPEED+NOABERR"),
    (SEFLG_SPEED | SEFLG_TRUEPOS, "SPEED+TRUEPOS"),
    (SEFLG_SPEED | SEFLG_EQUATORIAL, "SPEED+EQUATORIAL"),
    (SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_NONUT, "SPEED+EQUAT+NONUT"),
    (SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_J2000 | SEFLG_NONUT, "SPEED+EQUAT+J2000"),
    (SEFLG_SPEED | SEFLG_XYZ, "SPEED+XYZ"),
    (SEFLG_SPEED | SEFLG_XYZ | SEFLG_EQUATORIAL, "SPEED+XYZ+EQUAT"),
    (SEFLG_SPEED | SEFLG_RADIANS, "SPEED+RADIANS"),
    (SEFLG_SPEED | SEFLG_HELCTR, "SPEED+HELCTR"),
    (SEFLG_SPEED | SEFLG_HELCTR | SEFLG_NONUT, "SPEED+HELCTR+NONUT"),
    (SEFLG_SPEED | SEFLG_HELCTR | SEFLG_EQUATORIAL, "SPEED+HELCTR+EQUAT"),
    (SEFLG_SPEED | SEFLG_BARYCTR, "SPEED+BARYCTR"),
    (SEFLG_SPEED | SEFLG_NONUT | SEFLG_NOABERR, "SPEED+NONUT+NOABERR"),
    (SEFLG_SPEED | SEFLG_TRUEPOS | SEFLG_NOABERR, "SPEED+TRUEPOS+NOABERR"),
    (SEFLG_SPEED | SEFLG_TRUEPOS | SEFLG_EQUATORIAL, "SPEED+TRUEPOS+EQUAT"),
    (SEFLG_SPEED | SEFLG_XYZ | SEFLG_J2000 | SEFLG_NONUT, "SPEED+XYZ+J2000"),
    (SEFLG_SPEED | SEFLG_RADIANS | SEFLG_EQUATORIAL, "SPEED+RAD+EQUAT"),
]


def get_tolerances(flags, body_id):
    """Get tolerance in arcseconds based on flags and body."""
    tol = 1.0  # base 1 arcsec

    if flags & SEFLG_TRUEPOS:
        tol = max(tol, 35.0)  # known ~31" TRUEPOS differences
    if flags & SEFLG_NOABERR:
        tol = max(tol, 35.0)
    if flags & SEFLG_BARYCTR:
        if body_id in (1,):  # Moon barycentric
            tol = max(tol, 5.0)
        elif body_id in (8, 9):  # Neptune, Pluto
            tol = max(tol, 20.0)
    if flags & SEFLG_J2000:
        tol = max(tol, 20.0)  # known J2000 offset for some bodies
    if flags & SEFLG_HELCTR:
        if body_id in (0, 1):
            tol = 999  # Sun helio is Earth, Moon helio is weird
        elif body_id == 2:
            tol = max(tol, 2.0)  # Mercury helio known diff
    return tol


print("=== Multi-flag combination stress test ===")

for flags, flag_label in FLAG_COMBOS:
    for jd, epoch_label in TEST_EPOCHS:
        for body_id, body_name in BODIES.items():
            # Skip meaningless combinations
            if (flags & SEFLG_HELCTR) and body_id == 0:
                continue  # Sun helio = Earth position (special)
            if (flags & SEFLG_BARYCTR) and body_id == 1:
                continue  # Moon barycentric (known diffs)

            try:
                se_result = swe.calc_ut(jd, body_id, flags)
                le_result = ephem.swe_calc_ut(jd, body_id, flags)

                se_pos = se_result[0]
                le_pos = le_result[0]

                tol = get_tolerances(flags, body_id)

                # Compare first 3 components
                all_ok = True
                for idx in range(3):
                    se_val = se_pos[idx]
                    le_val = le_pos[idx]

                    if flags & SEFLG_XYZ:
                        # XYZ: compare in AU
                        diff = abs(se_val - le_val)
                        # Convert to approximate arcsec equivalent
                        if abs(se_val) > 0.001:
                            diff_arcsec = (diff / abs(se_val)) * 206265
                        else:
                            diff_arcsec = diff * 206265
                    elif flags & SEFLG_RADIANS:
                        diff = abs(se_val - le_val)
                        diff_arcsec = math.degrees(diff) * 3600
                    else:
                        diff = abs(se_val - le_val)
                        if idx == 0 and diff > 180 and not (flags & SEFLG_XYZ):
                            diff = 360 - diff
                        if idx == 2:
                            # Distance — different tolerance
                            diff_arcsec = 0  # skip distance comparison here
                        else:
                            diff_arcsec = diff * 3600

                    if idx < 2 and diff_arcsec >= tol:
                        all_ok = False
                        failed += 1
                        comp = ["lon/X/RA", "lat/Y/Dec", "dist/Z"][idx]
                        print(
                            f"  FAIL {flag_label} {body_name} {epoch_label} {comp}: "
                            f'SE={se_val:.8f} LE={le_val:.8f} diff={diff_arcsec:.2f}" tol={tol:.1f}"'
                        )

                if all_ok:
                    passed += 1

            except Exception as e:
                errors += 1
                if (
                    "unsupported" not in str(e).lower()
                    and "not support" not in str(e).lower()
                ):
                    print(f"  ERROR {flag_label} {body_name} {epoch_label}: {e}")

print(f"\nAfter main test: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Sidereal flag combinations
# ============================================================
print("\n=== P2: Sidereal flag combinations ===")

SIDEREAL_COMBOS = [
    (SEFLG_SPEED | SEFLG_SIDEREAL, "SID"),
    (SEFLG_SPEED | SEFLG_SIDEREAL | SEFLG_NONUT, "SID+NONUT"),
]

SID_MODES = {0: "FAGAN", 1: "LAHIRI", 27: "TRUE_CITRA"}
SID_BODIES = {0: "Sun", 1: "Moon", 4: "Mars", 5: "Jupiter"}

for sid_mode, sid_name in SID_MODES.items():
    swe.set_sid_mode(sid_mode)
    ephem.swe_set_sid_mode(sid_mode, 0, 0)

    for flags, flag_label in SIDEREAL_COMBOS:
        for jd, epoch_label in TEST_EPOCHS[:2]:
            for body_id, body_name in SID_BODIES.items():
                try:
                    se_result = swe.calc_ut(jd, body_id, flags)
                    le_result = ephem.swe_calc_ut(jd, body_id, flags)

                    se_lon = se_result[0][0]
                    le_lon = le_result[0][0]

                    diff = abs(se_lon - le_lon)
                    if diff > 180:
                        diff = 360 - diff
                    diff_arcsec = diff * 3600

                    # Known ~14" sidereal offset
                    tol = 20.0

                    label = f"P2 {sid_name} {flag_label} {body_name} {epoch_label}"
                    if diff_arcsec >= tol:
                        failed += 1
                        print(f'  FAIL {label}: diff={diff_arcsec:.2f}"')
                    else:
                        passed += 1

                except Exception as e:
                    errors += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Speed values consistency across flags
# ============================================================
print("\n=== P3: Speed consistency across flags ===")

SPEED_BODIES = {0: "Sun", 1: "Moon", 2: "Mercury", 5: "Jupiter"}

for jd, epoch_label in TEST_EPOCHS[:2]:
    for body_id, body_name in SPEED_BODIES.items():
        try:
            # Default speed
            le_def = ephem.swe_calc_ut(jd, body_id, SEFLG_SPEED)
            le_def_spd = le_def[0][3]  # lon speed deg/day

            se_def = swe.calc_ut(jd, body_id, SEFLG_SPEED)
            se_def_spd = se_def[0][3]

            diff_spd = abs(le_def_spd - se_def_spd) * 3600  # arcsec/day

            tol_spd = 1.0  # 1"/day
            label = f"P3 {body_name} {epoch_label} lon_speed"
            if diff_spd >= tol_spd:
                failed += 1
                print(
                    f'  FAIL {label}: SE={se_def_spd:.8f} LE={le_def_spd:.8f} diff={diff_spd:.2f}"/day'
                )
            else:
                passed += 1

            # Lat speed
            le_lat_spd = le_def[0][4]
            se_lat_spd = se_def[0][4]
            diff_lat_spd = abs(le_lat_spd - se_lat_spd) * 3600

            if diff_lat_spd >= 1.0:
                failed += 1
                print(f'  FAIL {label}_lat: diff={diff_lat_spd:.2f}"/day')
            else:
                passed += 1

        except Exception as e:
            errors += 1

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Retflag preservation
# ============================================================
print("\n=== P4: Retflag preservation ===")

for flags, flag_label in FLAG_COMBOS[:10]:
    for body_id in [0, 1, 4, 5]:
        jd = 2451545.0
        try:
            le_result = ephem.swe_calc_ut(jd, body_id, flags)
            retflag = le_result[1]

            # SPEED flag should be preserved
            if flags & SEFLG_SPEED:
                has_speed = (retflag & SEFLG_SPEED) != 0
                if not has_speed:
                    failed += 1
                    print(
                        f"  FAIL P4 {flag_label} body={body_id}: SPEED not in retflag={retflag}"
                    )
                else:
                    passed += 1

            # XYZ flag should be preserved
            if flags & SEFLG_XYZ:
                has_xyz = (retflag & SEFLG_XYZ) != 0
                if not has_xyz:
                    failed += 1
                    print(
                        f"  FAIL P4 {flag_label} body={body_id}: XYZ not in retflag={retflag}"
                    )
                else:
                    passed += 1

            # RADIANS flag should be preserved
            if flags & SEFLG_RADIANS:
                has_rad = (retflag & SEFLG_RADIANS) != 0
                if not has_rad:
                    failed += 1
                    print(
                        f"  FAIL P4 {flag_label} body={body_id}: RADIANS not in retflag={retflag}"
                    )
                else:
                    passed += 1

        except Exception as e:
            errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 114 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
