"""Verify magnitude and phase angle fixes against Swiss Ephemeris."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path("swisseph/ephe")

bodies = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
]

# Test across multiple dates
test_jds = [
    2451545.0,  # J2000
    2460000.5,  # 2023
    2459580.5,  # 2022
    2458849.5,  # 2020
    2455197.5,  # 2010
    2452275.5,  # 2002
    2448000.5,  # 1990
    2444000.5,  # 1979
    2440000.5,  # 1968
    2460500.5,  # 2024
]

print("=" * 100)
print(
    f"{'Body':10s} | {'Mag Diff':>10s} | {'Phase Diff':>12s} | {'Elong Diff':>12s} | {'Disc Diff%':>10s} | {'SE Mag':>8s} | {'LIB Mag':>8s}"
)
print("=" * 100)

max_diffs = {}

for jd in test_jds:
    for bid, name in bodies:
        try:
            sp = swe.pheno_ut(jd, bid, 256)
            lp = ephem.swe_pheno_ut(jd, bid, 256)

            mag_diff = abs(float(sp[4]) - float(lp[0][4]))
            phase_diff_arcsec = abs(float(sp[0]) - float(lp[0][0])) * 3600
            elong_diff_arcsec = abs(float(sp[2]) - float(lp[0][2])) * 3600

            se_disc = float(sp[3])
            lib_disc = float(lp[0][3])
            disc_diff_pct = (
                abs(se_disc - lib_disc) / se_disc * 100 if se_disc > 0 else 0
            )

            key = name
            if key not in max_diffs:
                max_diffs[key] = {"mag": 0, "phase": 0, "elong": 0, "disc": 0}
            max_diffs[key]["mag"] = max(max_diffs[key]["mag"], mag_diff)
            max_diffs[key]["phase"] = max(max_diffs[key]["phase"], phase_diff_arcsec)
            max_diffs[key]["elong"] = max(max_diffs[key]["elong"], elong_diff_arcsec)
            max_diffs[key]["disc"] = max(max_diffs[key]["disc"], disc_diff_pct)

        except Exception as e:
            print(f"{name:10s} JD={jd}: ERROR: {e}")

print("\nMAX DIFFERENCES ACROSS ALL DATES:")
print("=" * 90)
print("Body       | Max Mag Diff | Max Phase(as) | Max Elong(as) | Max Disc%")
print("=" * 90)

for name in [b[1] for b in bodies]:
    if name in max_diffs:
        d = max_diffs[name]
        mag_ok = "ok" if d["mag"] < 0.05 else "FAIL"
        phase_ok = "ok" if d["phase"] < 30 else "FAIL"
        print(
            f"{name:10s} | {d['mag']:12.4f} {mag_ok:4s} | {d['phase']:11.2f} {phase_ok:4s} | {d['elong']:11.2f} | {d['disc']:9.2f}%"
        )
print("=" * 90)

for name in [b[1] for b in bodies]:
    if name in max_diffs:
        d = max_diffs[name]
        mag_ok = "✓" if d["mag"] < 0.05 else "✗"
        phase_ok = "✓" if d["phase"] < 30 else "✗"
        print(
            f'{name:10s} | {d["mag"]:12.4f} {mag_ok} | {d["phase"]:10.2f}" {phase_ok} | {d["elong"]:10.2f}" | {d["disc"]:9.2f}%'
        )
