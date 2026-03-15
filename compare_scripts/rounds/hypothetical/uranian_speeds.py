"""Round 143: Uranian/Hypothetical Body Speeds Deep.

Tests all 8 hypothetical bodies (Cupido through Poseidon) speed accuracy
across multiple epochs and flag combinations.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SEFLG_SPEED,
    SEFLG_SWIEPH,
    SEFLG_J2000,
    SEFLG_NONUT,
    SEFLG_HELCTR,
    SE_CUPIDO,
    SE_HADES,
    SE_ZEUS,
    SE_KRONOS,
    SE_APOLLON,
    SE_ADMETOS,
    SE_VULKANUS,
    SE_POSEIDON,
)

swe.set_ephe_path("swisseph/ephe")
ephem.swe_set_ephe_path("swisseph/ephe")

BODIES = {
    SE_CUPIDO: "Cupido",
    SE_HADES: "Hades",
    SE_ZEUS: "Zeus",
    SE_KRONOS: "Kronos",
    SE_APOLLON: "Apollon",
    SE_ADMETOS: "Admetos",
    SE_VULKANUS: "Vulkanus",
    SE_POSEIDON: "Poseidon",
}

TEST_JDS = [
    2433282.5,
    2440587.5,
    2444239.5,
    2447892.5,
    2451545.0,
    2453371.5,
    2455197.5,
    2457023.5,
    2458849.5,
    2460676.5,
    2462502.5,
    2466154.5,
    2469807.5,
]

FLAG_COMBOS = {
    "default": SEFLG_SWIEPH | SEFLG_SPEED,
    "J2000": SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_J2000,
    "NONUT": SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_NONUT,
    "HELCTR": SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR,
}

COMP_NAMES = ["lon", "lat", "dist", "lon_spd", "lat_spd", "dist_spd"]

passed = 0
failed = 0
errors = 0
total = 0

for jd in TEST_JDS:
    for body, bname in BODIES.items():
        for flag_name, flags in FLAG_COMBOS.items():
            try:
                se_result = swe.calc_ut(jd, body, flags)
                le_result = ephem.swe_calc_ut(jd, body, flags)

                for i, comp in enumerate(COMP_NAMES):
                    total += 1
                    se_val = se_result[0][i]
                    le_val = le_result[0][i]
                    diff = abs(le_val - se_val)

                    if comp in ("dist", "dist_spd"):
                        tol = 0.05  # AU - these are approximate
                        is_pass = diff < tol
                        diff_str = f"{diff:.6f} AU"
                    else:
                        diff_arcsec = diff * 3600.0
                        if "spd" in comp:
                            tol = 120.0  # 2 arcmin/day for speeds
                        else:
                            tol = 60.0  # 1 arcmin position
                        is_pass = diff_arcsec < tol
                        diff_str = f'{diff_arcsec:.3f}"'

                    if is_pass:
                        passed += 1
                    else:
                        failed += 1
                        print(
                            f"FAIL {bname:10s} {comp:8s} {flag_name:8s} JD={jd:.1f} SE={se_val:+14.6f} LE={le_val:+14.6f} diff={diff_str}"
                        )

            except Exception as e:
                errors += 1
                print(f"ERR  {bname:10s} {flag_name:8s} JD={jd:.1f}: {e}")

print(f"\n{'=' * 60}")
print(f"Round 143: Uranian/Hypothetical Body Speeds Deep")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
print(f"Errors:  {errors}")
