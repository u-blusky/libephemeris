"""Find which dates cause Neptune magnitude discrepancy."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem
import math

swe.set_ephe_path("swisseph/ephe")

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

print("Neptune magnitude comparison per date:")
print("-" * 80)

for jd in test_jds:
    sp = swe.pheno_ut(jd, 8, 256)
    lp = ephem.swe_pheno_ut(jd, 8, 256)

    se_mag = float(sp[4])
    lib_mag = float(lp[0][4])
    diff = abs(se_mag - lib_mag)
    flag = " *** PROBLEM ***" if diff > 0.01 else ""

    print(f"  JD={jd}  SE={se_mag:.4f}  LIB={lib_mag:.4f}  diff={diff:.4f}{flag}")

# Also test Moon to find the crescent dates
print("\nMoon magnitude comparison per date:")
print("-" * 80)

for jd in test_jds:
    sp = swe.pheno_ut(jd, 1, 256)
    lp = ephem.swe_pheno_ut(jd, 1, 256)

    se_mag = float(sp[4])
    lib_mag = float(lp[0][4])
    se_phase = float(sp[0])
    diff = abs(se_mag - lib_mag)
    flag = " *** PROBLEM ***" if diff > 0.05 else ""

    print(
        f"  JD={jd}  SE={se_mag:.4f}  LIB={lib_mag:.4f}  diff={diff:.4f}  phase={se_phase:.2f}deg{flag}"
    )

# Now try Uranus too
print("\nUranus magnitude comparison per date:")
print("-" * 80)

for jd in test_jds:
    sp = swe.pheno_ut(jd, 7, 256)
    lp = ephem.swe_pheno_ut(jd, 7, 256)

    se_mag = float(sp[4])
    lib_mag = float(lp[0][4])
    diff = abs(se_mag - lib_mag)
    flag = " *** PROBLEM ***" if diff > 0.01 else ""

    print(f"  JD={jd}  SE={se_mag:.4f}  LIB={lib_mag:.4f}  diff={diff:.4f}{flag}")
