"""Analyze SE Moon magnitude formula to reverse-engineer phase function."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import math
import numpy as np
import swisseph as swe

swe.set_ephe_path("swisseph/ephe")

_MEAN_MOON_DISTANCE_AU = 0.00257

# Collect SE data points
results = []
for i in range(2000):
    jd = 2451545.0 + i * 0.5
    sp = swe.pheno_ut(jd, 1, 256)
    se_alpha = sp[0]
    se_mag = sp[4]
    moon = swe.calc_ut(jd, 1, 256)
    dist = moon[0][2]
    results.append((se_alpha, se_mag, dist))

results.sort(key=lambda x: x[0])

d_mean = 0.00257
v0 = -12.73

print("=== High phase angle detail (>140 deg) ===")
print(f"{'alpha':>8} {'SE_mag':>8} {'dist':>10} {'f(a)':>8} {'AA_f':>8} {'resid':>8}")

for alpha, se_mag, dist in results:
    if alpha < 140:
        continue
    dist_corr = 5.0 * math.log10(dist / d_mean)
    f_alpha = se_mag - v0 - dist_corr
    aa_f = 0.026 * alpha + 4e-9 * alpha**4
    residual = f_alpha - aa_f
    print(
        f"{alpha:8.2f} {se_mag:8.4f} {dist:10.6f} {f_alpha:8.4f} {aa_f:8.4f} {residual:+8.4f}"
    )

# Now fit the full range properly
# Isolate f(alpha) = se_mag - v0 - dist_corr
alphas = []
f_values = []
for alpha, se_mag, dist in results:
    dist_corr = 5.0 * math.log10(dist / d_mean)
    f_alpha = se_mag - v0 - dist_corr
    alphas.append(alpha)
    f_values.append(f_alpha)

alphas = np.array(alphas)
f_values = np.array(f_values)

# Check what SE does at very high phase angles
# SE source uses: mag = -12.55 + 0.026*phase + 4e-9*phase^4
# OR: mag = -12.73 + 0.026*|i| + 4.0E-9*|i|^4  (Astronomical Almanac 2nd ed)
# Let's test: maybe SE uses no distance correction at all?

print("\n\n=== Testing: does SE apply distance correction? ===")
# If SE has NO distance correction, then SE_mag = v0 + f(alpha)
# So f(alpha) = SE_mag - v0
# Then we should see no correlation with distance at a given alpha

# Group by similar alpha values and check if distance affects magnitude
from collections import defaultdict

alpha_bins = defaultdict(list)
for alpha, se_mag, dist in results:
    bin_key = round(alpha / 5) * 5  # 5-degree bins
    alpha_bins[bin_key].append((alpha, se_mag, dist))

print(f"{'bin':>6} {'n':>4} {'mag_std':>8} {'dist_std':>10} {'corr':>8}")
for bin_key in sorted(alpha_bins.keys()):
    items = alpha_bins[bin_key]
    if len(items) < 3:
        continue
    mags = [x[1] for x in items]
    dists = [x[2] for x in items]
    mag_std = np.std(mags)
    dist_std = np.std(dists)
    if dist_std > 0 and mag_std > 0:
        corr = np.corrcoef(mags, dists)[0, 1]
    else:
        corr = 0
    print(f"{bin_key:6.0f} {len(items):4d} {mag_std:8.4f} {dist_std:10.6f} {corr:8.4f}")

# Test exact SE formula: V = -12.73 + 0.026*alpha + 4e-9*alpha^4
# with and without distance correction
print("\n\n=== Formula accuracy comparison ===")
for label, use_dist in [("With dist corr", True), ("No dist corr", False)]:
    errs_all = []
    errs_lt150 = []
    for alpha, se_mag, dist in results:
        dist_corr = 5.0 * math.log10(dist / d_mean) if use_dist else 0.0
        calc = -12.73 + dist_corr + 0.026 * alpha + 4e-9 * alpha**4
        err = abs(se_mag - calc)
        errs_all.append(err)
        if alpha < 150:
            errs_lt150.append(err)
    print(
        f"{label:20s}: all mean={np.mean(errs_all):.4f} max={np.max(errs_all):.4f} | <150 mean={np.mean(errs_lt150):.4f} max={np.max(errs_lt150):.4f}"
    )

# Try with -12.55 as V0 (some SE versions use this)
for v0_test in [-12.73, -12.55, -12.74, -12.70]:
    errs = []
    for alpha, se_mag, dist in results:
        dist_corr = 5.0 * math.log10(dist / d_mean)
        calc = v0_test + dist_corr + 0.026 * alpha + 4e-9 * alpha**4
        err = abs(se_mag - calc)
        errs.append(err)
    print(
        f"V0={v0_test:7.2f} with dist:  mean={np.mean(errs):.4f} max={np.max(errs):.4f}"
    )

# What if SE uses a different mean distance?
# Try the IAU conventional value: 384400 km = 0.002569 AU
for d_test in [0.002569, 0.002570, 0.002574, 0.00256955]:
    errs = []
    for alpha, se_mag, dist in results:
        dist_corr = 5.0 * math.log10(dist / d_test)
        calc = -12.73 + dist_corr + 0.026 * alpha + 4e-9 * alpha**4
        err = abs(se_mag - calc)
        errs.append(err)
    print(f"d_mean={d_test:.6f}: mean={np.mean(errs):.4f} max={np.max(errs):.4f}")
