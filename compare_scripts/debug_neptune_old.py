"""Debug Neptune magnitude at old dates in detail."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem
import math

swe.set_ephe_path("swisseph/ephe")

# Problem dates
test_jds = [
    (2460000.5, "2023 (ok)"),
    (2451545.0, "2000 (ok)"),
    (2448000.5, "1990 (bad)"),
    (2444000.5, "1979 (bad)"),
    (2440000.5, "1968 (bad)"),
    (2436000.5, "1957"),
    (2432000.5, "1946"),
]

for jd, label in test_jds:
    sp = swe.pheno_ut(jd, 8, 256)
    lp = ephem.swe_pheno_ut(jd, 8, 256)

    se_mag = float(sp[4])
    lib_mag = float(lp[0][4])

    # Get distances
    pos_geo = swe.calc_ut(jd, 8, 256)
    pos_geo_lib = ephem.swe_calc_ut(jd, 8, 256)
    se_delta = float(pos_geo[0][2])
    lib_delta = float(pos_geo_lib[0][2])

    pos_helio = swe.calc_ut(jd, 8, 256 | 8)
    pos_helio_lib = ephem.swe_calc_ut(jd, 8, 256 | 8)
    se_r = float(pos_helio[0][2])
    lib_r = float(pos_helio_lib[0][2])

    se_phase_deg = float(sp[0])
    lib_phase_deg = float(lp[0][0])

    # What V0 would SE need to give its magnitude?
    # V = V0 + 5*log10(r*delta) + phase_correction
    dist_term_se = 5 * math.log10(se_r * se_delta)
    implied_V0_se = se_mag - dist_term_se

    dist_term_lib = 5 * math.log10(lib_r * lib_delta)
    implied_V0_lib = lib_mag - dist_term_lib

    print(f"JD={jd} ({label})")
    print(
        f"  SE:  mag={se_mag:.4f}  r={se_r:.6f}  delta={se_delta:.6f}  phase={se_phase_deg:.4f}"
    )
    print(
        f"  LIB: mag={lib_mag:.4f}  r={lib_r:.6f}  delta={lib_delta:.6f}  phase={lib_phase_deg:.4f}"
    )
    print(f"  Implied V0: SE={implied_V0_se:.4f}  LIB={implied_V0_lib:.4f}")
    print(
        f"  Dist diff: r={abs(se_r - lib_r):.6f}  delta={abs(se_delta - lib_delta):.6f}"
    )
    print(f"  Mag diff: {se_mag - lib_mag:.4f}")
    print()
