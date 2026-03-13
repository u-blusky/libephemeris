"""Find Neptune's secular magnitude variation pattern in SE."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import math

swe.set_ephe_path("swisseph/ephe")

# Scan year by year from 1900 to 2025
print("Year       JD           SE_mag    r          delta      dist_term  implied_V0")
print("-" * 95)

for year in range(1900, 2026, 1):
    jd = 2451545.0 + (year - 2000) * 365.25  # approximate JD

    try:
        sp = swe.pheno_ut(jd, 8, 256)
        se_mag = float(sp[4])

        pos_geo = swe.calc_ut(jd, 8, 256)
        pos_helio = swe.calc_ut(jd, 8, 256 | 8)
        delta = float(pos_geo[0][2])
        r = float(pos_helio[0][2])

        dist_term = 5 * math.log10(r * delta)
        implied_V0 = se_mag - dist_term

        print(
            f"{year}  {jd:12.1f}  {se_mag:8.4f}  {r:9.6f}  {delta:9.6f}  {dist_term:9.4f}  {implied_V0:8.4f}"
        )
    except Exception as e:
        print(f"{year}  ERROR: {e}")
