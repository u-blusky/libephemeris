#!/usr/bin/env python3
"""
Round 23: cotrans / cotrans_sp Comprehensive Verification
==========================================================

Tests coordinate transformation functions (ecliptic <-> equatorial <-> horizon)
and their speed counterparts.

Parts:
  P1: Ecliptic → Equatorial (ecl2equ) at multiple obliquities
  P2: Equatorial → Ecliptic (equ2ecl) at multiple obliquities
  P3: Round-trip: ecl→equ→ecl identity test
  P4: cotrans_sp — speed transformation
  P5: Edge cases: poles, zero coordinates, 360° wrapping
  P6: Obliquity sweep: transform at obliquities from 22° to 24.5°
  P7: Full-sky grid: systematic coverage of coordinate space
"""

from __future__ import annotations

import os
import sys
import time
import traceback
import math

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as ephem
from libephemeris.constants import *

_EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
swe.set_ephe_path(_EPHE_PATH)


class R:
    def __init__(self, name):
        self.name = name
        self.passed = self.failed = self.skipped = 0
        self.failures = []
        self.max_diff = 0.0
        self.max_label = ""

    def ok(self, diff=0.0, label=""):
        self.passed += 1
        if diff > self.max_diff:
            self.max_diff = diff
            self.max_label = label

    def fail(self, msg):
        self.failed += 1
        self.failures.append(msg)

    def skip(self, msg=""):
        self.skipped += 1

    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'=' * 70}")
        print(f"  {self.name}: {self.passed}/{total} PASSED ({self.skipped} skip)")
        if self.max_diff > 0:
            print(f"  Max diff: {self.max_diff:.10f}° ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def angle_diff(a, b):
    """Signed angular difference handling wraparound."""
    d = a - b
    while d > 180:
        d -= 360
    while d < -180:
        d += 360
    return abs(d)


# Standard obliquity at J2000
EPS_J2000 = 23.4392911111  # IAU 2000 mean obliquity at J2000


# ============================================================
# PART 1: Ecliptic → Equatorial
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Ecliptic → Equatorial at J2000 obliquity")
    print("=" * 70)

    r = R("P1: Ecl→Equ")

    eps = EPS_J2000

    # Test points: (lon, lat) in ecliptic
    test_points = [
        (0.0, 0.0, "Vernal equinox"),
        (90.0, 0.0, "Summer solstice"),
        (180.0, 0.0, "Autumnal equinox"),
        (270.0, 0.0, "Winter solstice"),
        (45.0, 0.0, "45° ecl"),
        (135.0, 0.0, "135° ecl"),
        (225.0, 0.0, "225° ecl"),
        (315.0, 0.0, "315° ecl"),
        (0.0, 23.4, "North ecl pole area"),
        (0.0, -23.4, "South ecl pole area"),
        (90.0, 45.0, "High lat summer"),
        (180.0, -30.0, "Neg lat autumn"),
        (30.0, 5.0, "Typical planet"),
        (120.0, -3.0, "Typical planet 2"),
        (200.0, 1.5, "Typical planet 3"),
        (350.0, -0.5, "Near 360°"),
        (0.1, 0.0, "Near zero"),
        (359.9, 0.0, "Near 360°"),
        (90.0, 89.0, "Near north ecl pole"),
        (270.0, -89.0, "Near south ecl pole"),
    ]

    for lon, lat, desc in test_points:
        label = f"ecl({lon},{lat}) {desc}"

        # SE: cotrans with negative obliquity = ecl→equ
        se_result = swe.cotrans((lon, lat, 1.0), -eps)
        le_result = ephem.cotrans((lon, lat, 1.0), -eps)

        lon_diff = angle_diff(se_result[0], le_result[0])
        lat_diff = abs(se_result[1] - le_result[1])
        max_d = max(lon_diff, lat_diff)

        tol = 1e-8  # ~0.04 mas
        if max_d > tol:
            r.fail(f"{label}: lon_diff={lon_diff:.10f}° lat_diff={lat_diff:.10f}°")
        else:
            r.ok(max_d, label)

    print(f"  Tested {len(test_points)} points")
    return r.summary(), r


# ============================================================
# PART 2: Equatorial → Ecliptic
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Equatorial → Ecliptic at J2000 obliquity")
    print("=" * 70)

    r = R("P2: Equ→Ecl")

    eps = EPS_J2000

    # Test points: (RA, Dec) in equatorial
    test_points = [
        (0.0, 0.0, "RA=0 Dec=0"),
        (90.0, 0.0, "RA=6h Dec=0"),
        (180.0, 0.0, "RA=12h Dec=0"),
        (270.0, 0.0, "RA=18h Dec=0"),
        (0.0, 90.0, "North celestial pole"),
        (0.0, -90.0, "South celestial pole"),
        (45.0, 30.0, "Typical star"),
        (120.0, -15.0, "Southern star"),
        (200.0, 60.0, "High dec star"),
        (350.0, -45.0, "Southern star 2"),
        (83.633, 22.014, "Aldebaran area"),
        (101.287, -16.716, "Sirius area"),
        (213.915, 19.182, "Arcturus area"),
        (0.0, 23.44, "Near ecl pole dec"),
        (180.0, -23.44, "Near south ecl pole dec"),
    ]

    for ra, dec, desc in test_points:
        label = f"equ({ra},{dec}) {desc}"

        # SE: cotrans with positive obliquity = equ→ecl
        se_result = swe.cotrans((ra, dec, 1.0), eps)
        le_result = ephem.cotrans((ra, dec, 1.0), eps)

        lon_diff = angle_diff(se_result[0], le_result[0])
        lat_diff = abs(se_result[1] - le_result[1])
        max_d = max(lon_diff, lat_diff)

        tol = 1e-8
        if max_d > tol:
            r.fail(f"{label}: lon_diff={lon_diff:.10f}° lat_diff={lat_diff:.10f}°")
        else:
            r.ok(max_d, label)

    print(f"  Tested {len(test_points)} points")
    return r.summary(), r


# ============================================================
# PART 3: Round-trip identity: ecl→equ→ecl
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Round-trip ecl→equ→ecl identity")
    print("=" * 70)

    r = R("P3: Round-trip")

    eps = EPS_J2000

    # Generate test points on ecliptic grid
    lons = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    lats = [-80, -60, -30, -10, 0, 10, 30, 60, 80]

    for lon in lons:
        for lat in lats:
            label = f"ecl({lon},{lat})"

            # LE round-trip
            equ = ephem.cotrans((float(lon), float(lat), 1.0), -eps)
            ecl_back = ephem.cotrans(equ, eps)

            lon_diff = angle_diff(lon, ecl_back[0])
            lat_diff = abs(lat - ecl_back[1])
            max_d = max(lon_diff, lat_diff)

            tol = 1e-10  # Round-trip should be near-perfect
            if max_d > tol:
                r.fail(f"{label}: lon_err={lon_diff:.12f}° lat_err={lat_diff:.12f}°")
            else:
                r.ok(max_d, label)

    print(f"  Tested {len(lons) * len(lats)} round-trips")
    return r.summary(), r


# ============================================================
# PART 4: cotrans_sp — speed transformation
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: cotrans_sp — speed transformation")
    print("=" * 70)

    r = R("P4: cotrans_sp")

    eps = EPS_J2000

    # Test with position + speed tuples
    # Format: (lon, lat, dist, lon_speed, lat_speed, dist_speed)
    test_cases = [
        # (position, speed, description)
        ((30.0, 5.0, 1.0), (1.0, 0.1, 0.0), "Typical planet"),
        ((90.0, 0.0, 1.0), (0.5, 0.0, 0.0), "Solstice point"),
        ((180.0, -3.0, 1.0), (0.8, -0.05, 0.0), "Retrograde area"),
        ((270.0, 2.0, 1.0), (-0.3, 0.02, 0.0), "Retrograde"),
        ((0.0, 0.0, 1.0), (1.2, 0.0, 0.0), "Equinox"),
        ((45.0, 20.0, 1.0), (0.0, 0.5, 0.0), "Lat-only motion"),
        ((120.0, -45.0, 0.5), (0.3, -0.1, 0.01), "Close + south"),
        ((350.0, 0.0, 2.0), (1.5, 0.0, -0.001), "Far + fast"),
        # Moon-like speeds
        ((100.0, 5.0, 0.00257), (13.0, 0.5, 0.0), "Moon-like"),
        # Mercury-like speeds
        ((200.0, 2.0, 0.8), (1.8, 0.1, 0.0), "Mercury-like"),
    ]

    for pos, spd, desc in test_cases:
        label = f"sp {desc}"

        # SE cotrans_sp: takes 6-tuple (lon,lat,dist,lon_spd,lat_spd,dist_spd) + eps
        # Returns 6-tuple
        se_result = swe.cotrans_sp(pos + spd, -eps)

        # LE cotrans_sp: takes two 3-tuples + obliquity
        # Returns (pos_3tuple, speed_3tuple)
        le_pos, le_spd = ephem.cotrans_sp(pos, spd, -eps)
        le_result = le_pos + le_spd  # combine to 6-tuple for comparison

        # Compare all 6 components
        max_d = 0.0
        for i, (sv, lv) in enumerate(zip(se_result, le_result)):
            d = abs(sv - lv)
            if i == 0:  # lon/RA — handle wraparound
                d = angle_diff(sv, lv)
            if d > max_d:
                max_d = d

        tol = 1e-8
        if max_d > tol:
            r.fail(f"{label}: max_diff={max_d:.12f}")
        else:
            r.ok(max_d, label)

    # Also test equ→ecl direction
    for pos, spd, desc in test_cases[:5]:
        label = f"sp_rev {desc}"

        se_result = swe.cotrans_sp(pos + spd, eps)
        le_pos, le_spd = ephem.cotrans_sp(pos, spd, eps)
        le_result = le_pos + le_spd

        max_d = 0.0
        for i, (sv, lv) in enumerate(zip(se_result, le_result)):
            d = abs(sv - lv)
            if i == 0:
                d = angle_diff(sv, lv)
            if d > max_d:
                max_d = d

        tol = 1e-8
        if max_d > tol:
            r.fail(f"{label}: max_diff={max_d:.12f}")
        else:
            r.ok(max_d, label)

    print(f"  Tested {len(test_cases) + 5} speed transforms")
    return r.summary(), r


# ============================================================
# PART 5: Edge cases
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Edge cases — poles, wrapping, extreme values")
    print("=" * 70)

    r = R("P5: Edge Cases")

    eps = EPS_J2000

    edge_cases = [
        # Ecliptic poles
        (0.0, 90.0, "North ecliptic pole"),
        (0.0, -90.0, "South ecliptic pole"),
        (180.0, 90.0, "NP at lon=180"),
        # Near-pole
        (0.0, 89.999, "Near north pole"),
        (0.0, -89.999, "Near south pole"),
        # Wrapping
        (359.999, 0.0, "Near 360°"),
        (0.001, 0.0, "Near 0°"),
        (359.999, 0.001, "Near (360,0)"),
        # Zero
        (0.0, 0.0, "Origin"),
        # Exact quarters
        (90.0, 0.0, "90°"),
        (180.0, 0.0, "180°"),
        (270.0, 0.0, "270°"),
        (360.0, 0.0, "360° = 0°"),
    ]

    for lon, lat, desc in edge_cases:
        label = f"edge {desc}"

        try:
            se_result = swe.cotrans((lon, lat, 1.0), -eps)
            le_result = ephem.cotrans((lon, lat, 1.0), -eps)

            lon_diff = angle_diff(se_result[0], le_result[0])
            lat_diff = abs(se_result[1] - le_result[1])
            max_d = max(lon_diff, lat_diff)

            tol = 1e-6  # Slightly wider for edge cases
            if max_d > tol:
                r.fail(f"{label}: lon_diff={lon_diff:.10f}° lat_diff={lat_diff:.10f}°")
            else:
                r.ok(max_d, label)

        except Exception as e:
            r.fail(f"{label}: EXCEPTION {e}")

    print(f"  Tested {len(edge_cases)} edge cases")
    return r.summary(), r


# ============================================================
# PART 6: Obliquity sweep
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Obliquity sweep — transform at varying obliquities")
    print("=" * 70)

    r = R("P6: Obliquity Sweep")

    # Obliquity varies from ~22.1° to ~24.5° over 41000-year cycle
    obliquities = [22.0, 22.5, 23.0, 23.44, 23.5, 24.0, 24.5]

    # Fixed test point
    test_points = [
        (45.0, 5.0),
        (120.0, -10.0),
        (200.0, 30.0),
        (330.0, -45.0),
    ]

    for eps in obliquities:
        for lon, lat in test_points:
            label = f"eps={eps}° ({lon},{lat})"

            se_result = swe.cotrans((lon, lat, 1.0), -eps)
            le_result = ephem.cotrans((lon, lat, 1.0), -eps)

            lon_diff = angle_diff(se_result[0], le_result[0])
            lat_diff = abs(se_result[1] - le_result[1])
            max_d = max(lon_diff, lat_diff)

            tol = 1e-8
            if max_d > tol:
                r.fail(f"{label}: max_diff={max_d:.12f}°")
            else:
                r.ok(max_d, label)

        print(f"  eps={eps}°: tested {len(test_points)} points")

    return r.summary(), r


# ============================================================
# PART 7: Full-sky grid
# ============================================================
def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: Full-sky grid — systematic 15° coverage")
    print("=" * 70)

    r = R("P7: Full-Sky Grid")

    eps = EPS_J2000

    lons = range(0, 360, 15)
    lats = range(-75, 76, 15)

    max_diff_overall = 0.0

    for lon in lons:
        for lat in lats:
            label = f"grid({lon},{lat})"

            # ecl→equ
            se_equ = swe.cotrans((float(lon), float(lat), 1.0), -eps)
            le_equ = ephem.cotrans((float(lon), float(lat), 1.0), -eps)

            lon_d = angle_diff(se_equ[0], le_equ[0])
            lat_d = abs(se_equ[1] - le_equ[1])
            max_d = max(lon_d, lat_d)

            if max_d > max_diff_overall:
                max_diff_overall = max_d

            tol = 1e-8
            if max_d > tol:
                r.fail(f"{label}: diff={max_d:.12f}°")
            else:
                r.ok(max_d, label)

    n_total = len(list(lons)) * len(list(lats))
    print(
        f"  Tested {r.passed + r.failed} grid points, max_diff={max_diff_overall:.12f}°"
    )

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 23: cotrans / cotrans_sp Comprehensive Verification")
    print("=" * 70)

    start = time.time()
    all_ok = True
    all_results = []

    for pname, pfn in [
        ("P1", run_part1),
        ("P2", run_part2),
        ("P3", run_part3),
        ("P4", run_part4),
        ("P5", run_part5),
        ("P6", run_part6),
        ("P7", run_part7),
    ]:
        try:
            ok, res = pfn()
            all_results.append((pname, res))
            if not ok:
                all_ok = False
        except Exception as e:
            print(f"\n  {pname} CRASHED: {e}")
            traceback.print_exc()
            all_ok = False

    elapsed = time.time() - start

    print("\n" + "=" * 70)
    print("ROUND 23 FINAL SUMMARY")
    print("=" * 70)

    tp = tf = ts = 0
    for pn, res in all_results:
        st = "PASS" if res.failed == 0 else "FAIL"
        t = res.passed + res.failed
        print(f"  {pn} {res.name}: {res.passed}/{t} ({res.skipped} skip) [{st}]")
        tp += res.passed
        tf += res.failed
        ts += res.skipped

    print(f"\n  TOTAL: {tp}/{tp + tf} PASSED, {tf} FAILED, {ts} SKIPPED")
    print(f"  Time: {elapsed:.1f}s")
    if all_ok:
        print(f"\n  >>> ROUND 23: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 23: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
