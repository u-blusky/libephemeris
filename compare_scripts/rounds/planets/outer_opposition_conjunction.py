#!/usr/bin/env python3
"""Round 123: Outer Planet Opposition/Conjunction Precision

Tests positions of Jupiter, Saturn, Uranus, Neptune at their opposition
and conjunction points — critical moments for astrological interpretation.
"""

from __future__ import annotations
import os, sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path("swisseph/ephe")

SEFLG_SPEED = 256
SE_SUN = 0
SE_JUPITER = 5
SE_SATURN = 6
SE_URANUS = 7
SE_NEPTUNE = 8


def find_oppositions_conjunctions(body_id, start_jd, count=15):
    """Find opposition/conjunction times by tracking Sun-planet elongation."""
    events = []
    jd = start_jd
    step = 1.0
    found = 0

    def elongation(jd_):
        sun = swe.calc_ut(jd_, SE_SUN, SEFLG_SPEED)[0][0]
        planet = swe.calc_ut(jd_, body_id, SEFLG_SPEED)[0][0]
        e = planet - sun
        while e > 180:
            e -= 360
        while e < -180:
            e += 360
        return e

    prev_elong = elongation(jd)
    jd += step

    for _ in range(50000):
        if found >= count:
            break
        curr_elong = elongation(jd)

        # Detect sign change (conjunction) or crossing ±180 (opposition)
        if prev_elong * curr_elong < 0 and abs(prev_elong) < 90:
            # Conjunction (elongation crosses 0)
            a, b = jd - step, jd
            for _ in range(50):
                mid = (a + b) / 2
                em = elongation(mid)
                if prev_elong > 0:
                    if em > 0:
                        a = mid
                    else:
                        b = mid
                else:
                    if em < 0:
                        a = mid
                    else:
                        b = mid
            events.append(("conjunction", (a + b) / 2))
            found += 1
        elif abs(prev_elong) > 170 and abs(curr_elong) > 170:
            # Check for opposition (elongation crosses ±180)
            if (prev_elong > 0 and curr_elong < 0) or (
                prev_elong < 0 and curr_elong > 0
            ):
                a, b = jd - step, jd
                for _ in range(50):
                    mid = (a + b) / 2
                    em = elongation(mid)
                    if abs(em) > 179.5:
                        if (em > 0) == (prev_elong > 0):
                            a = mid
                        else:
                            b = mid
                    else:
                        b = mid
                events.append(("opposition", (a + b) / 2))
                found += 1

        prev_elong = curr_elong
        jd += step

    return events


def main():
    print("=" * 80)
    print("ROUND 123: Outer Planet Opposition/Conjunction Precision")
    print("=" * 80)

    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    bodies = [
        (SE_JUPITER, "Jupiter", 20),
        (SE_SATURN, "Saturn", 15),
        (SE_URANUS, "Uranus", 10),
        (SE_NEPTUNE, "Neptune", 8),
    ]

    start_jd = 2451545.0

    for body_id, body_name, count in bodies:
        print(f"\n  Finding events for {body_name}...")
        events = find_oppositions_conjunctions(body_id, start_jd, count)
        print(f"  Found {len(events)} events")

        for etype, jd in events:
            # Test at event and nearby
            for offset in [-0.1, 0.0, 0.1]:
                test_jd = jd + offset
                try:
                    se_pos = swe.calc_ut(test_jd, body_id, SEFLG_SPEED)[0]
                    le_pos = ephem.swe_calc_ut(test_jd, body_id, SEFLG_SPEED)[0]
                except Exception:
                    continue

                for i, (label, tol) in enumerate(
                    [
                        ("lon", 2.0),
                        ("lat", 2.0),
                        ("dist", None),
                        ("lon_spd", 5.0),
                        ("lat_spd", 5.0),
                        ("dist_spd", None),
                    ]
                ):
                    if tol is None:
                        if abs(se_pos[i]) > 1e-10:
                            rel = abs(le_pos[i] - se_pos[i]) / abs(se_pos[i])
                            total_tests += 1
                            if rel < 1e-4:
                                total_pass += 1
                            else:
                                total_fail += 1
                                if len(failures) < 20:
                                    failures.append(
                                        f"  {body_name} {etype} {label}: rel={rel:.2e}"
                                    )
                        continue

                    diff = le_pos[i] - se_pos[i]
                    if i == 0:
                        if diff > 180:
                            diff -= 360
                        elif diff < -180:
                            diff += 360
                    diff_as = abs(diff) * 3600

                    total_tests += 1
                    if diff_as < tol:
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 30:
                            failures.append(
                                f'  {body_name} {etype} JD={test_jd:.2f} {label}: {diff_as:.4f}"'
                            )

    print(f"\n{'=' * 80}")
    pct = 100 * total_pass / total_tests if total_tests else 0
    print(f"ROUND 123 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)")
    print(f"  Failures: {total_fail}")
    print("=" * 80)
    if failures:
        print("\nSample failures:")
        for f in failures:
            print(f)
    if total_fail == 0:
        print("\nAll tests PASSED!")
    return total_fail


if __name__ == "__main__":
    sys.exit(main())
