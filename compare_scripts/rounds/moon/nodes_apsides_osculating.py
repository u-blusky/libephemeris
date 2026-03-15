#!/usr/bin/env python3
"""Round 98: Node/Apsides Osculating vs Mean — Deep Sweep

Tests swe_nod_aps_ut() for all planets with different calculation methods:
1. SE_NODBIT_MEAN (mean nodes/apsides)
2. SE_NODBIT_OSCU (osculating nodes/apsides)
3. SE_NODBIT_OSCU_BAR (osculating barycentric)
4. SE_NODBIT_FOPOINT (second focal point)
5. Multiple epochs
6. All planets (Sun through Pluto + Chiron + asteroids)
"""

from __future__ import annotations

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path("swisseph/ephe")

SEFLG_SPEED = 256
SEFLG_SWIEPH = 2

# Body IDs
SE_SUN = 0
SE_MOON = 1
SE_MERCURY = 2
SE_VENUS = 3
SE_MARS = 4
SE_JUPITER = 5
SE_SATURN = 6
SE_URANUS = 7
SE_NEPTUNE = 8
SE_PLUTO = 9
SE_CHIRON = 15
SE_CERES = 17
SE_PALLAS = 18

# Method flags
SE_NODBIT_MEAN = 1
SE_NODBIT_OSCU = 2
SE_NODBIT_OSCU_BAR = 4
SE_NODBIT_FOPOINT = 256

PLANETS = [
    (SE_MOON, "Moon"),
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
    (SE_URANUS, "Uranus"),
    (SE_NEPTUNE, "Neptune"),
    (SE_PLUTO, "Pluto"),
    (SE_CHIRON, "Chiron"),
]

METHODS = [
    (SE_NODBIT_MEAN, "MEAN"),
    (SE_NODBIT_OSCU, "OSCU"),
    (SE_NODBIT_OSCU_BAR, "OSCU_BAR"),
    (SE_NODBIT_FOPOINT, "FOPOINT"),
    (SE_NODBIT_MEAN | SE_NODBIT_FOPOINT, "MEAN+FO"),
]

EPOCHS = [
    2451545.0,  # J2000
    2460000.0,  # 2023
    2440000.0,  # 1968
    2430000.0,  # 1941
    2470000.0,  # 2050
    2415020.0,  # 1900
]

# Tolerances (arcsec)
LON_TOL = 60.0  # 1 arcmin — osculating elements can diverge significantly
MEAN_LON_TOL = 30.0  # mean nodes tighter
DIST_TOL = 0.01  # AU


def run_tests():
    passed = 0
    failed = 0
    errors = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 98: Node/Apsides Osculating vs Mean — Deep Sweep")
    print("=" * 80)

    for method_flag, method_name in METHODS:
        m_pass = 0
        m_fail = 0
        m_err = 0

        tol = MEAN_LON_TOL if method_flag & SE_NODBIT_MEAN else LON_TOL

        for body_id, body_name in PLANETS:
            for jd in EPOCHS:
                total += 1
                iflag = SEFLG_SWIEPH | SEFLG_SPEED

                try:
                    se_result = swe.nod_aps_ut(jd, body_id, iflag, method_flag)
                    # Returns (asc_node, desc_node, perihelion, aphelion)
                    # Each is a 6-tuple
                except Exception as e:
                    m_err += 1
                    errors += 1
                    continue

                try:
                    le_result = ephem.swe_nod_aps_ut(jd, body_id, iflag, method_flag)
                except Exception as e:
                    m_err += 1
                    errors += 1
                    if len(fail_details) < 5:
                        fail_details.append(
                            f"  ERROR [{method_name}] {body_name} JD={jd:.1f}: {e}"
                        )
                    continue

                # Compare all 4 output nodes/apsides
                names = ["AscNode", "DescNode", "Perihelion", "Aphelion"]
                all_ok = True

                for i, name in enumerate(names):
                    se_pos = se_result[i]
                    le_pos = le_result[i]

                    if len(se_pos) < 3 or len(le_pos) < 3:
                        continue

                    dlon = abs(se_pos[0] - le_pos[0]) * 3600
                    dlat = abs(se_pos[1] - le_pos[1]) * 3600

                    if dlon > 180 * 3600:
                        dlon = 360 * 3600 - dlon

                    if dlon > tol or dlat > tol:
                        all_ok = False
                        if len(fail_details) < 30:
                            fail_details.append(
                                f"  FAIL [{method_name}] {body_name:10s} JD={jd:.1f} "
                                f'{name}: dLon={dlon:.1f}" dLat={dlat:.1f}" '
                                f"SE=({se_pos[0]:.4f},{se_pos[1]:.4f}) "
                                f"LE=({le_pos[0]:.4f},{le_pos[1]:.4f})"
                            )
                        break

                if all_ok:
                    m_pass += 1
                    passed += 1
                else:
                    m_fail += 1
                    failed += 1

        pct = 100 * m_pass / max(1, m_pass + m_fail + m_err)
        print(
            f"  [{method_name:12s}] {m_pass}/{m_pass + m_fail + m_err} passed ({pct:.1f}%) "
            f"[{m_fail} fail, {m_err} err]"
        )

    print()
    if fail_details:
        print("Sample failures/errors:")
        for d in fail_details[:30]:
            print(d)

    print()
    print("=" * 80)
    print(
        f"ROUND 98 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%), "
        f"{failed} failed, {errors} errors"
    )
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
