"""Investigate magnitude gaps: Sun, Moon thin crescent, Neptune, Uranus."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import math
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path("swisseph/ephe")

print("=== 1. SUN MAGNITUDE ===")
# Sun magnitude: constant 0.12 offset
for jd in [2451545.0, 2453000.0, 2455000.0, 2457000.0, 2459000.0]:
    sp = swe.pheno_ut(jd, 0, 256)
    lp = ephem.swe_pheno_ut(jd, 0, 256)
    se_mag = float(sp[4])
    lib_mag = float(lp[0][4])
    se_dist = swe.calc_ut(jd, 0, 256)[0][2]
    print(
        f"  JD={jd}: SE_mag={se_mag:.6f} LIB_mag={lib_mag:.6f} diff={se_mag - lib_mag:+.6f} dist={se_dist:.6f}AU"
    )

# SE Sun formula: V_sun = -26.74 + 5*log10(r)
# Our formula: V_sun = -26.74 + 5*log10(r) if r > 0 else -26.74
# Check if SE uses -26.86 or different V0
print("\n  Testing different V0 for Sun:")
jd = 2451545.0
sp = swe.pheno_ut(jd, 0, 256)
se_mag = float(sp[4])
se_dist = swe.calc_ut(jd, 0, 256)[0][2]
for v0 in [-26.74, -26.76, -26.77, -26.80, -26.86]:
    calc = v0 + 5.0 * math.log10(se_dist)
    print(f"    V0={v0}: calc={calc:.6f} SE={se_mag:.6f} diff={abs(calc - se_mag):.6f}")

print("\n=== 2. MOON MAGNITUDE AT THIN CRESCENT ===")
# The 5.6 mag worst case is at JD 2453676.6 - check the phase angle
jd = 2453676.6
sp = swe.pheno_ut(jd, 1, 256)
lp = ephem.swe_pheno_ut(jd, 1, 256)
print(f"  JD={jd}: SE_alpha={float(sp[0]):.2f} LIB_alpha={float(lp[0][0]):.2f}")
print(
    f"  SE_mag={float(sp[4]):.4f} LIB_mag={float(lp[0][4]):.4f} diff={abs(float(sp[4]) - float(lp[0][4])):.4f}"
)
# Check more thin crescent dates
print("\n  Moon magnitude at high phase angles:")
for i in range(2000):
    jd_test = 2451545.0 + i * 0.5
    sp = swe.pheno_ut(jd_test, 1, 256)
    alpha = float(sp[0])
    if alpha > 165:
        lp = ephem.swe_pheno_ut(jd_test, 1, 256)
        se_mag = float(sp[4])
        lib_mag = float(lp[0][4])
        diff = abs(se_mag - lib_mag)
        if diff > 0.5:
            print(
                f"    JD={jd_test:.1f} alpha={alpha:.2f} SE_mag={se_mag:.2f} LIB_mag={lib_mag:.2f} diff={diff:.2f}"
            )

print("\n=== 3. NEPTUNE MAGNITUDE ===")
# Constant 0.13 offset suggests wrong V(1,0) or formula
for jd in [2451545.0, 2453000.0, 2455000.0, 2457000.0, 2459000.0]:
    sp = swe.pheno_ut(jd, 8, 256)
    lp = ephem.swe_pheno_ut(jd, 8, 256)
    se_mag = float(sp[4])
    lib_mag = float(lp[0][4])
    print(
        f"  JD={jd}: SE_mag={se_mag:.4f} LIB_mag={lib_mag:.4f} diff={se_mag - lib_mag:+.4f}"
    )

print("\n=== 4. URANUS MAGNITUDE ===")
for jd in [2451545.0, 2453000.0, 2455000.0, 2457000.0, 2459000.0]:
    sp = swe.pheno_ut(jd, 7, 256)
    lp = ephem.swe_pheno_ut(jd, 7, 256)
    se_mag = float(sp[4])
    lib_mag = float(lp[0][4])
    print(
        f"  JD={jd}: SE_mag={se_mag:.4f} LIB_mag={lib_mag:.4f} diff={se_mag - lib_mag:+.4f}"
    )

print("\n=== 5. PLUTO MAGNITUDE ===")
for jd in [2451545.0, 2453000.0, 2455000.0, 2457000.0, 2459000.0]:
    sp = swe.pheno_ut(jd, 9, 256)
    lp = ephem.swe_pheno_ut(jd, 9, 256)
    se_mag = float(sp[4])
    lib_mag = float(lp[0][4])
    print(
        f"  JD={jd}: SE_mag={se_mag:.4f} LIB_mag={lib_mag:.4f} diff={se_mag - lib_mag:+.4f}"
    )
