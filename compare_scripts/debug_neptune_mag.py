"""Debug Neptune magnitude calculation in detail."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem
import math

swe.set_ephe_path("swisseph/ephe")

# Test multiple dates
test_jds = [2451545.0, 2460000.5, 2459580.5, 2458849.5, 2455197.5]

for jd in test_jds:
    sp = swe.pheno_ut(jd, 8, 256)  # Neptune = 8
    lp = ephem.swe_pheno_ut(jd, 8, 256)

    se_phase = float(sp[0])
    se_elong = float(sp[2])
    se_disc = float(sp[3])
    se_mag = float(sp[4])

    lib_phase = float(lp[0][0])
    lib_elong = float(lp[0][2])
    lib_disc = float(lp[0][3])
    lib_mag = float(lp[0][4])

    print(f"JD={jd}")
    print(
        f"  SE:  phase={se_phase:.6f}  elong={se_elong:.6f}  disc={se_disc:.8f}  mag={se_mag:.4f}"
    )
    print(
        f"  LIB: phase={lib_phase:.6f}  elong={lib_elong:.6f}  disc={lib_disc:.8f}  mag={lib_mag:.4f}"
    )
    print(f"  Diff: mag={se_mag - lib_mag:.4f}")

    # Manual calc: V = V0 + 5*log10(r*delta)
    # We need r (helio dist) and delta (geo dist)
    # From SE pheno: sp[1] = phase angle fraction (not useful here)
    # Let's get positions to compute distances
    pos_se = swe.calc_ut(jd, 8, 256)  # SEFLG_SPEED
    pos_lib = ephem.swe_calc_ut(jd, 8, 256)

    # SE calc returns ecliptic lon, lat, dist (AU from earth)
    se_delta = float(pos_se[0][2])  # geocentric distance
    lib_delta = float(pos_lib[0][2])

    # Get heliocentric distance
    pos_se_h = swe.calc_ut(jd, 8, 256 | 8)  # SEFLG_HELCTR = 8
    pos_lib_h = ephem.swe_calc_ut(jd, 8, 256 | 8)
    se_r = float(pos_se_h[0][2])
    lib_r = float(pos_lib_h[0][2])

    print(f"  SE:  r={se_r:.6f} AU  delta={se_delta:.6f} AU")
    print(f"  LIB: r={lib_r:.6f} AU  delta={lib_delta:.6f} AU")

    # Manual mag calculation: V = V0 + 5*log10(r*delta)
    se_manual_mag = -7.00 + 5 * math.log10(se_r * se_delta)
    lib_manual_mag = -7.00 + 5 * math.log10(lib_r * lib_delta)

    print(f"  Manual(V0=-7.00): SE={se_manual_mag:.4f}  LIB={lib_manual_mag:.4f}")
    print(f"  SE actual vs manual(SE): {se_mag - se_manual_mag:.4f}")
    print(f"  LIB actual vs manual(LIB): {lib_mag - lib_manual_mag:.4f}")
    print()
