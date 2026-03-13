"""Deep precision audit round 1: pheno_ut across 500+ dates for all bodies."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import math
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
labels = ["phase_angle", "phase", "elongation", "disc_diam", "magnitude"]

# Test across 500 dates spanning 20 years
dates = [2451545.0 + i * 14.6 for i in range(500)]

print("=== DEEP PRECISION AUDIT: pheno_ut across 500 dates ===\n")

for body_id, body_name in bodies:
    max_diffs = [0.0] * 5
    mean_diffs = [0.0] * 5
    count = 0
    worst_dates = [0.0] * 5

    for jd in dates:
        try:
            sp = swe.pheno_ut(jd, body_id, 256)
            lp = ephem.swe_pheno_ut(jd, body_id, 256)
        except Exception:
            continue

        count += 1
        for i in range(5):
            sv = float(sp[i])
            lv = float(lp[0][i])
            diff = abs(sv - lv)
            mean_diffs[i] += diff
            if diff > max_diffs[i]:
                max_diffs[i] = diff
                worst_dates[i] = jd

    if count > 0:
        mean_diffs = [d / count for d in mean_diffs]

    print(f"--- {body_name} ({count} dates) ---")
    for i, label in enumerate(labels):
        if label == "phase_angle":
            print(
                f'  {label:15s}: max={max_diffs[i] * 3600:8.1f}"  mean={mean_diffs[i] * 3600:8.2f}"  worst_jd={worst_dates[i]:.1f}'
            )
        elif label == "magnitude":
            print(
                f"  {label:15s}: max={max_diffs[i]:8.4f}mag  mean={mean_diffs[i]:8.5f}mag  worst_jd={worst_dates[i]:.1f}"
            )
        elif label == "disc_diam":
            # relative diff
            try:
                sp_worst = swe.pheno_ut(worst_dates[i], body_id, 256)
                rel = (
                    max_diffs[i] / float(sp_worst[i]) * 100
                    if float(sp_worst[i]) > 0
                    else 0
                )
            except Exception:
                rel = 0
            print(
                f"  {label:15s}: max={max_diffs[i]:.8f}deg  rel={rel:.2f}%  mean={mean_diffs[i]:.8f}deg"
            )
        else:
            print(f"  {label:15s}: max={max_diffs[i]:.8f}  mean={mean_diffs[i]:.8f}")
    print()

# Now test specific edge cases
print("\n=== EDGE CASES ===\n")

# Test at very early/late dates
edge_dates = [
    (2378496.5, "1800-01-01"),
    (2415020.5, "1900-01-01"),
    (2440587.5, "1970-01-01"),
    (2451545.0, "2000-01-01 (J2000)"),
    (2460310.5, "2024-01-01"),
    (2469807.5, "2050-01-01"),
    (2488069.5, "2100-01-01"),
]

for jd, label in edge_dates:
    print(f"--- JD={jd} ({label}) ---")
    for body_id, body_name in [(1, "Moon"), (2, "Mercury"), (5, "Jupiter")]:
        try:
            sp = swe.pheno_ut(jd, body_id, 256)
            lp = ephem.swe_pheno_ut(jd, body_id, 256)
            pa_diff = abs(float(sp[0]) - float(lp[0][0])) * 3600
            mag_diff = abs(float(sp[4]) - float(lp[0][4]))
            print(
                f'  {body_name:10s}: phase_angle diff={pa_diff:.1f}"  mag diff={mag_diff:.4f}'
            )
        except Exception as e:
            print(f"  {body_name:10s}: ERROR: {e}")
    print()
