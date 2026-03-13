"""
Deep audit: libephemeris vs pyswisseph at extreme dates.

Tests calc_ut positions for all 10 major bodies (Sun through Pluto) at
boundary and historical dates across the DE440s/DE440 ephemeris ranges.
"""

from __future__ import annotations

import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SEFLG_SWIEPH, SEFLG_SPEED

swe.set_ephe_path("swisseph/ephe")

# Default flags for calc_ut calls
DEFAULT_FLAGS = SEFLG_SWIEPH | SEFLG_SPEED

# Body constants
SUN = 0
MOON = 1
MERCURY = 2
VENUS = 3
MARS = 4
JUPITER = 5
SATURN = 6
URANUS = 7
NEPTUNE = 8
PLUTO = 9

BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
]

# Test dates organized by category
DATE_CATEGORIES = {
    "DE440s boundary": [1850, 2149],
    "DE440 boundary": [1551, 2649],
    "Historical": [1600, 1700, 1800],
    "Far future": [2200, 2500],
    "Modern reference": [2000, 2025],
}

# Thresholds
LON_THRESHOLD = 0.001  # degrees (= 3.6 arcsec)
LAT_THRESHOLD = 0.001  # degrees
DIST_THRESHOLD = 0.0001  # AU


def angular_diff(a: float, b: float) -> float:
    """Angular difference accounting for 360-degree wrap."""
    d = abs(a - b) % 360
    return min(d, 360 - d)


def arcsec(deg: float) -> float:
    return deg * 3600.0


# Storage for all results
results = []  # list of dicts
flags_list = []  # flagged results


def test_body_at_jd(body_id: int, body_name: str, jd: float, year: int, category: str):
    """Compare libephemeris vs pyswisseph for one body at one JD."""
    result = {
        "category": category,
        "year": year,
        "body": body_name,
        "jd": jd,
        "lon_lib": None,
        "lat_lib": None,
        "dist_lib": None,
        "lon_swe": None,
        "lat_swe": None,
        "dist_swe": None,
        "dlon": None,
        "dlat": None,
        "ddist": None,
        "error": None,
        "flagged": False,
        "flag_reasons": [],
    }

    try:
        r_swe = swe.calc_ut(jd, body_id, DEFAULT_FLAGS)
        s = r_swe[0]
        result["lon_swe"] = s[0]
        result["lat_swe"] = s[1]
        result["dist_swe"] = s[2]
    except Exception as e:
        result["error"] = f"swe: {e}"
        results.append(result)
        return

    try:
        r_lib = ephem.swe_calc_ut(jd, body_id, DEFAULT_FLAGS)
        l = r_lib[0]
        result["lon_lib"] = l[0]
        result["lat_lib"] = l[1]
        result["dist_lib"] = l[2]
    except Exception as e:
        result["error"] = f"lib: {e}"
        results.append(result)
        return

    dlon = angular_diff(s[0], l[0])
    dlat = abs(s[1] - l[1])
    ddist = abs(s[2] - l[2])

    result["dlon"] = dlon
    result["dlat"] = dlat
    result["ddist"] = ddist

    if dlon > LON_THRESHOLD:
        result["flagged"] = True
        result["flag_reasons"].append(f"lon {dlon:.6f}° > {LON_THRESHOLD}°")
    if dlat > LAT_THRESHOLD:
        result["flagged"] = True
        result["flag_reasons"].append(f"lat {dlat:.6f}° > {LAT_THRESHOLD}°")
    if ddist > DIST_THRESHOLD:
        result["flagged"] = True
        result["flag_reasons"].append(f"dist {ddist:.6f} AU > {DIST_THRESHOLD} AU")

    results.append(result)
    if result["flagged"]:
        flags_list.append(result)


# ============================================================================
# Run all tests
# ============================================================================
print("=" * 90)
print("DEEP AUDIT: libephemeris vs pyswisseph at extreme dates")
print("=" * 90)
print(
    f'Thresholds: lon > {LON_THRESHOLD}° ({arcsec(LON_THRESHOLD):.1f}")  '
    f'lat > {LAT_THRESHOLD}° ({arcsec(LAT_THRESHOLD):.1f}")  '
    f"dist > {DIST_THRESHOLD} AU"
)
print()

for category, years in DATE_CATEGORIES.items():
    print(f"--- {category} ---")
    for year in years:
        jd = swe.julday(year, 1, 1, 12.0)
        print(f"  Year {year}  (JD {jd:.1f})")
        for body_id, body_name in BODIES:
            test_body_at_jd(body_id, body_name, jd, year, category)
    print()

# ============================================================================
# Full summary table
# ============================================================================
print()
print("=" * 90)
print("FULL SUMMARY TABLE")
print("=" * 90)
print()

# Header
col_dlon_deg = "dLon(deg)"
col_dlon_as = "dLon(arcsec)"
col_dlat_deg = "dLat(deg)"
col_ddist = "dDist(AU)"
hdr = (
    f"{'Year':>6s}  {'Category':<17s}  {'Body':<9s}  "
    f"{col_dlon_deg:>12s}  {col_dlon_as:>12s}  "
    f"{col_dlat_deg:>12s}  {col_ddist:>12s}  "
    f"{'Status':>6s}"
)
print(hdr)
print("-" * len(hdr))

for r in results:
    if r["error"]:
        print(
            f"{r['year']:>6d}  {r['category']:<17s}  {r['body']:<9s}  "
            f"{'ERROR':>12s}  {'':>10s}  {'':>12s}  {'':>12s}  "
            f"  {r['error']}"
        )
        continue

    status = "FAIL" if r["flagged"] else "OK"
    dlon_as = arcsec(r["dlon"])
    print(
        f"{r['year']:>6d}  {r['category']:<17s}  {r['body']:<9s}  "
        f"{r['dlon']:>12.8f}  {dlon_as:>10.4f}  "
        f"{r['dlat']:>12.8f}  {r['ddist']:>12.8f}  "
        f"{status:>6s}"
    )

# ============================================================================
# Aggregated statistics per year
# ============================================================================
print()
print("=" * 90)
print("AGGREGATED STATISTICS PER YEAR")
print("=" * 90)
print()

from collections import defaultdict

year_stats = defaultdict(
    lambda: {
        "max_dlon": 0,
        "max_dlat": 0,
        "max_ddist": 0,
        "n": 0,
        "errors": 0,
        "flags": 0,
    }
)
for r in results:
    ys = year_stats[r["year"]]
    if r["error"]:
        ys["errors"] += 1
        continue
    ys["n"] += 1
    ys["max_dlon"] = max(ys["max_dlon"], r["dlon"])
    ys["max_dlat"] = max(ys["max_dlat"], r["dlat"])
    ys["max_ddist"] = max(ys["max_ddist"], r["ddist"])
    if r["flagged"]:
        ys["flags"] += 1

col_maxdlon_deg = "MaxDLon(deg)"
col_maxdlon_as = "MaxDLon(asec)"
col_maxdlat_deg = "MaxDLat(deg)"
col_maxddist = "MaxDDist(AU)"
hdr2 = (
    f"{'Year':>6s}  {'Bodies':>6s}  {'Errors':>6s}  {'Flags':>5s}  "
    f"{col_maxdlon_deg:>12s}  {col_maxdlon_as:>13s}  "
    f"{col_maxdlat_deg:>12s}  {col_maxddist:>13s}"
)
print(hdr2)
print("-" * len(hdr2))

for year in sorted(year_stats.keys()):
    ys = year_stats[year]
    print(
        f"{year:>6d}  {ys['n']:>6d}  {ys['errors']:>6d}  {ys['flags']:>5d}  "
        f"{ys['max_dlon']:>12.8f}  {arcsec(ys['max_dlon']):>11.4f}  "
        f"{ys['max_dlat']:>12.8f}  {ys['max_ddist']:>13.8f}"
    )

# ============================================================================
# Flagged results (concerning)
# ============================================================================
print()
print("=" * 90)
print("FLAGGED RESULTS (exceeding thresholds)")
print("=" * 90)
print()

if not flags_list:
    print(
        "  No results exceed thresholds. All body/date combinations are within tolerance."
    )
else:
    print(f"  Total flagged: {len(flags_list)} out of {len(results)} tests")
    print()
    for r in flags_list:
        reasons = "; ".join(r["flag_reasons"])
        dlon_as = arcsec(r["dlon"])
        print(
            f"  *** Year {r['year']} | {r['body']:<9s} | "
            f'dLon={r["dlon"]:.8f}° ({dlon_as:.4f}") | '
            f"dLat={r['dlat']:.8f}° | "
            f"dDist={r['ddist']:.8f} AU"
        )
        print(f"      Reason: {reasons}")
        print(
            f"      swe: lon={r['lon_swe']:.6f}  lat={r['lat_swe']:.6f}  dist={r['dist_swe']:.6f}"
        )
        print(
            f"      lib: lon={r['lon_lib']:.6f}  lat={r['lat_lib']:.6f}  dist={r['dist_lib']:.6f}"
        )
        print()

# ============================================================================
# Worst cases per body (across all dates)
# ============================================================================
print("=" * 90)
print("WORST CASE PER BODY (max longitude diff across all dates)")
print("=" * 90)
print()

body_worst = {}
for r in results:
    if r["error"]:
        continue
    name = r["body"]
    if name not in body_worst or r["dlon"] > body_worst[name]["dlon"]:
        body_worst[name] = r

hdr3 = (
    f"{'Body':<9s}  {'Year':>6s}  "
    f"{col_dlon_deg:>12s}  {col_dlon_as:>12s}  "
    f"{col_dlat_deg:>12s}  {col_ddist:>12s}  "
    f"{'Status':>6s}"
)
print(hdr3)
print("-" * len(hdr3))

for _, bname in BODIES:
    if bname not in body_worst:
        continue
    r = body_worst[bname]
    status = "FAIL" if r["flagged"] else "OK"
    print(
        f"{bname:<9s}  {r['year']:>6d}  "
        f"{r['dlon']:>12.8f}  {arcsec(r['dlon']):>10.4f}  "
        f"{r['dlat']:>12.8f}  {r['ddist']:>12.8f}  "
        f"{status:>6s}"
    )

# ============================================================================
# Final verdict
# ============================================================================
print()
print("=" * 90)
total = len(results)
errors = sum(1 for r in results if r["error"])
ok = sum(1 for r in results if not r["error"] and not r["flagged"])
flagged = len(flags_list)

print(f"FINAL VERDICT: {total} tests | {ok} OK | {flagged} FLAGGED | {errors} ERRORS")

if flagged == 0 and errors == 0:
    print(
        "ALL CLEAR: libephemeris matches pyswisseph within thresholds at all extreme dates."
    )
elif flagged > 0:
    print(
        f"ATTENTION: {flagged} test(s) exceed precision thresholds. See flagged results above."
    )
if errors > 0:
    print(f"NOTE: {errors} test(s) produced errors (likely ephemeris range issues).")

print("=" * 90)
