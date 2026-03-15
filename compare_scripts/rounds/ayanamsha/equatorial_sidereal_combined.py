"""Round 139: All Planets Equatorial+Sidereal Combined.

Tests the combination of SEFLG_EQUATORIAL + SEFLG_SIDEREAL for all planets.
This is a tricky flag combination where ayanamsha must be applied correctly
to RA/Dec coordinates.
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SEFLG_SPEED,
    SEFLG_SWIEPH,
    SEFLG_SIDEREAL,
    SEFLG_EQUATORIAL,
    SEFLG_NONUT,
    SEFLG_J2000,
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
    SE_CHIRON,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
)

swe.set_ephe_path("swisseph/ephe")
ephem.swe_set_ephe_path("swisseph/ephe")

BODIES = {
    SE_SUN: "Sun",
    SE_MOON: "Moon",
    SE_MERCURY: "Mercury",
    SE_VENUS: "Venus",
    SE_MARS: "Mars",
    SE_JUPITER: "Jupiter",
    SE_SATURN: "Saturn",
    SE_URANUS: "Uranus",
    SE_NEPTUNE: "Neptune",
    SE_PLUTO: "Pluto",
    SE_CHIRON: "Chiron",
    SE_CERES: "Ceres",
    SE_PALLAS: "Pallas",
    SE_JUNO: "Juno",
    SE_VESTA: "Vesta",
    SE_MEAN_NODE: "MeanNode",
    SE_TRUE_NODE: "TrueNode",
    SE_MEAN_APOG: "MeanApog",
    SE_OSCU_APOG: "OscuApog",
}

TEST_JDS = [
    2433282.5,
    2440587.5,
    2444239.5,
    2447892.5,
    2451545.0,
    2455197.5,
    2458849.5,
    2460676.5,
    2462502.5,
    2466154.5,
]

SIDEREAL_MODES = [0, 1, 3, 5, 27]  # Fagan, Lahiri, Raman, Krishnamurti, TrueCitra

FLAG_COMBOS = {
    "SID+EQ": SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_SIDEREAL | SEFLG_EQUATORIAL,
    "SID+EQ+NONUT": SEFLG_SWIEPH
    | SEFLG_SPEED
    | SEFLG_SIDEREAL
    | SEFLG_EQUATORIAL
    | SEFLG_NONUT,
}

COMP_NAMES = ["RA", "Dec", "dist", "RA_spd", "Dec_spd", "dist_spd"]

passed = 0
failed = 0
errors = 0
total = 0

for sid_mode in SIDEREAL_MODES:
    swe.set_sid_mode(sid_mode)
    ephem.swe_set_sid_mode(sid_mode)

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

                        if comp == "dist" or comp == "dist_spd":
                            tol = 0.001
                            is_pass = diff < tol
                        else:
                            diff_arcsec = diff * 3600.0
                            if body in (
                                SE_MEAN_NODE,
                                SE_TRUE_NODE,
                                SE_MEAN_APOG,
                                SE_OSCU_APOG,
                            ):
                                tol = 120.0  # 2 arcmin for analytical bodies
                            elif "spd" in comp and body == SE_OSCU_APOG:
                                tol = 200.0
                            else:
                                tol = 60.0  # 1 arcmin (includes ~14" sidereal offset)
                            is_pass = diff_arcsec < tol

                        if is_pass:
                            passed += 1
                        else:
                            diff_arcsec = diff * 3600.0
                            info = f'FAIL {bname:10s} {comp:8s} sid={sid_mode} JD={jd} {flag_name} SE={se_val:+14.6f} LE={le_val:+14.6f} diff={diff_arcsec:.3f}"'
                            print(info)
                            failed += 1

                except Exception as e:
                    errors += 1
                    if "Invalid Time" not in str(e):
                        print(
                            f"ERR  {bname:10s} sid={sid_mode} JD={jd} {flag_name}: {e}"
                        )

swe.set_sid_mode(0)
ephem.swe_set_sid_mode(0)

print(f"\n{'=' * 60}")
print(f"Round 139: All Planets Equatorial+Sidereal Combined")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
print(f"Errors:  {errors}")
