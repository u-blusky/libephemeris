"""
Debug script: Investigate barycentric Sun precision discrepancy.

CRITICAL FINDING: Previous audit used FLG_BARYCTR = 4, but that's actually
SEFLG_MOSEPH (Moshier ephemeris flag), NOT SEFLG_BARYCTR (16384).

This script tests BOTH values to confirm the root cause.
"""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SEFLG_BARYCTR,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_HELCTR,
    SE_SUN,
    SE_EARTH,
)

swe.set_ephe_path("swisseph/ephe")

# ============================================================================
# Flag value audit
# ============================================================================

print("=" * 100)
print("FLAG VALUE AUDIT - is the previous test using the right constant?")
print("=" * 100)
print(
    "  libephemeris SEFLG_BARYCTR = {} (0x{:04X})".format(SEFLG_BARYCTR, SEFLG_BARYCTR)
)
print("  libephemeris SEFLG_MOSEPH  = {} (0x{:04X})".format(SEFLG_MOSEPH, SEFLG_MOSEPH))

# Check pyswisseph's own constants
swe_baryctr = getattr(swe, "FLG_BARYCTR", None) or getattr(swe, "SEFLG_BARYCTR", None)
swe_moseph = getattr(swe, "FLG_MOSEPH", None) or getattr(swe, "SEFLG_MOSEPH", None)
print("  pyswisseph  FLG_BARYCTR   = {}".format(swe_baryctr))
print("  pyswisseph  FLG_MOSEPH    = {}".format(swe_moseph))
print()
print("  Previous audit used FLG_BARYCTR = 4  <--- THIS IS SEFLG_MOSEPH, NOT BARYCTR!")
print("  The 0.038 deg 'discrepancy' was likely Moshier vs JPL, not a barycentric bug.")


# ============================================================================
# Helpers
# ============================================================================


def year_to_jd(y: int, m: int = 1, d: int = 1) -> float:
    return swe.julday(y, m, d, 12.0)


def angular_diff(a: float, b: float) -> float:
    d = a - b
    while d > 180:
        d -= 360
    while d < -180:
        d += 360
    return d


# 20 test dates
TEST_DATES = [
    (1950, 1, 1),
    (1955, 6, 15),
    (1960, 3, 21),
    (1965, 9, 10),
    (1970, 1, 1),
    (1975, 7, 4),
    (1980, 12, 21),
    (1985, 4, 15),
    (1990, 1, 1),
    (1993, 8, 20),
    (1997, 3, 15),
    (2000, 1, 1),
    (2005, 6, 21),
    (2010, 1, 1),
    (2015, 9, 23),
    (2020, 1, 1),
    (2025, 3, 20),
    (2030, 12, 25),
    (2040, 6, 1),
    (2050, 1, 1),
]


# ============================================================================
# PART 1: Test with WRONG flag (4 = MOSEPH) - reproducing previous audit
# ============================================================================

print()
print("=" * 100)
print("PART 1: Using flag=4 (SEFLG_MOSEPH) - reproducing the previous audit's 'bug'")
print("  This is what the previous audit called 'FLG_BARYCTR = 4'")
print("=" * 100)
print(
    "{:>12s}  {:>12s} {:>12s} {:>10s}  {:>10s} {:>10s} {:>12s}".format(
        "Date", "LE lon", "SE lon", "dLon asec", "LE dist", "SE dist", "dDist"
    )
)
print("-" * 100)

WRONG_FLAG = 4  # This is MOSEPH, not BARYCTR
max_dlon_wrong = 0.0

for y, m, d in TEST_DATES:
    jd = year_to_jd(y, m, d)

    le = ephem.swe_calc_ut(jd, 0, WRONG_FLAG)[0]
    se = swe.calc_ut(jd, 0, WRONG_FLAG)[0]

    dlon = angular_diff(le[0], se[0]) * 3600  # arcsec
    ddist = le[2] - se[2]
    max_dlon_wrong = max(max_dlon_wrong, abs(dlon))

    print(
        "{:>4d}-{:02d}-{:02d}    {:>12.6f} {:>12.6f} {:>+10.3f}  "
        "{:>10.6f} {:>10.6f} {:>+12.8f}".format(
            y, m, d, le[0], se[0], dlon, le[2], se[2], ddist
        )
    )

print()
print(
    '  Max |dLon| with flag=4 (MOSEPH): {:.3f}" ({:.6f} deg)'.format(
        max_dlon_wrong, max_dlon_wrong / 3600
    )
)
print(
    "  pyswisseph uses Moshier analytical ephemeris; libephemeris ignores MOSEPH, always uses JPL."
)


# ============================================================================
# PART 2: Test with CORRECT flag (16384 = SEFLG_BARYCTR)
# ============================================================================

print()
print("=" * 100)
print(
    "PART 2: Using flag={} (SEFLG_BARYCTR) - the REAL barycentric test".format(
        SEFLG_BARYCTR
    )
)
print("=" * 100)
print(
    "{:>12s}  {:>12s} {:>12s} {:>10s}  {:>10s} {:>10s} {:>12s}".format(
        "Date", "LE lon", "SE lon", "dLon asec", "LE dist", "SE dist", "dDist"
    )
)
print("-" * 100)

BARY_FLAG = SEFLG_BARYCTR
max_dlon_correct = 0.0

for y, m, d in TEST_DATES:
    jd = year_to_jd(y, m, d)

    le = ephem.swe_calc_ut(jd, 0, BARY_FLAG)[0]
    se = swe.calc_ut(jd, 0, BARY_FLAG)[0]

    dlon = angular_diff(le[0], se[0]) * 3600  # arcsec
    ddist = le[2] - se[2]
    max_dlon_correct = max(max_dlon_correct, abs(dlon))

    print(
        "{:>4d}-{:02d}-{:02d}    {:>12.6f} {:>12.6f} {:>+10.3f}  "
        "{:>10.6f} {:>10.6f} {:>+12.8f}".format(
            y, m, d, le[0], se[0], dlon, le[2], se[2], ddist
        )
    )

print()
print(
    '  Max |dLon| with flag={} (BARYCTR): {:.3f}" ({:.6f} deg)'.format(
        SEFLG_BARYCTR, max_dlon_correct, max_dlon_correct / 3600
    )
)


# ============================================================================
# PART 3: Is libephemeris ignoring BARYCTR flag for Sun?
#         Compare geocentric vs barycentric from libephemeris
# ============================================================================

print()
print("=" * 100)
print("PART 3: Is libephemeris returning GEOCENTRIC when asked for BARYCENTRIC Sun?")
print("  Barycentric Sun distance should be small (~0.005-0.01 AU, SSB wobble)")
print("  Geocentric Sun distance should be ~1 AU")
print("=" * 100)

for y, m, d in [(2000, 1, 1), (2010, 1, 1), (2020, 1, 1)]:
    jd = year_to_jd(y, m, d)

    le_geo = ephem.swe_calc_ut(jd, 0, 0)[0]
    le_bary = ephem.swe_calc_ut(jd, 0, BARY_FLAG)[0]
    se_geo = swe.calc_ut(jd, 0, 0)[0]
    se_bary = swe.calc_ut(jd, 0, BARY_FLAG)[0]

    le_diff = angular_diff(le_geo[0], le_bary[0])
    se_diff = angular_diff(se_geo[0], se_bary[0])

    if abs(le_diff) < 0.001 and abs(le_geo[2] - le_bary[2]) < 0.001:
        assessment = "FLAG IGNORED!"
    elif le_bary[2] > 0.1:
        assessment = "DIST WRONG (~1 AU)"
    else:
        assessment = "OK"

    print()
    print("  {}-{:02d}-{:02d}  Assessment: {}".format(y, m, d, assessment))
    print("    lib  GEO:  lon={:12.6f}  dist={:12.8f}".format(le_geo[0], le_geo[2]))
    print("    lib  BARY: lon={:12.6f}  dist={:12.8f}".format(le_bary[0], le_bary[2]))
    print("    swe  GEO:  lon={:12.6f}  dist={:12.8f}".format(se_geo[0], se_geo[2]))
    print("    swe  BARY: lon={:12.6f}  dist={:12.8f}".format(se_bary[0], se_bary[2]))
    print("    lib GEO-BARY lon diff: {:+.6f} deg".format(le_diff))
    print("    swe GEO-BARY lon diff: {:+.6f} deg".format(se_diff))


# ============================================================================
# PART 4: Distance sanity check across all 20 dates
# ============================================================================

print()
print("=" * 100)
print("PART 4: Barycentric Sun distance sanity (should be ~0.005-0.01 AU)")
print("=" * 100)

for y, m, d in TEST_DATES:
    jd = year_to_jd(y, m, d)

    le_bary_dist = ephem.swe_calc_ut(jd, 0, BARY_FLAG)[0][2]
    se_bary_dist = swe.calc_ut(jd, 0, BARY_FLAG)[0][2]
    le_geo_dist = ephem.swe_calc_ut(jd, 0, 0)[0][2]

    flag = ""
    if le_bary_dist > 0.1:
        flag = " *** TOO LARGE - looks geocentric! ***"

    print(
        "  {:>4d}-{:02d}-{:02d}: lib_bary={:>12.8f}  swe_bary={:>12.8f}  "
        "lib_geo={:>10.6f}{}".format(
            y, m, d, le_bary_dist, se_bary_dist, le_geo_dist, flag
        )
    )


# ============================================================================
# PART 5: Return flag analysis
# ============================================================================

print()
print("=" * 100)
print("PART 5: Return flag analysis - is BARYCTR acknowledged?")
print("=" * 100)

jd = year_to_jd(2000)

for label, iflag in [
    ("Geocentric (flag=0)", 0),
    ("MOSEPH (flag=4)", SEFLG_MOSEPH),
    ("BARYCTR (flag={})".format(SEFLG_BARYCTR), SEFLG_BARYCTR),
    (
        "BARYCTR+SPEED (flag={})".format(SEFLG_BARYCTR | SEFLG_SPEED),
        SEFLG_BARYCTR | SEFLG_SPEED,
    ),
    ("HELCTR (flag={})".format(SEFLG_HELCTR), SEFLG_HELCTR),
]:
    le = ephem.swe_calc_ut(jd, 0, iflag)
    se = swe.calc_ut(jd, 0, iflag)

    le_rf = le[1]
    se_rf = se[1] if len(se) > 1 and isinstance(se[1], int) else "N/A"

    print()
    print("  {}:".format(label))
    print(
        "    lib retflag={} (0x{:05X})  BARYCTR={}  MOSEPH={}  HELCTR={}".format(
            le_rf,
            le_rf,
            "SET" if le_rf & SEFLG_BARYCTR else "no",
            "SET" if le_rf & SEFLG_MOSEPH else "no",
            "SET" if le_rf & SEFLG_HELCTR else "no",
        )
    )
    if isinstance(se_rf, int):
        print(
            "    swe retflag={} (0x{:05X})  BARYCTR={}  MOSEPH={}  HELCTR={}".format(
                se_rf,
                se_rf,
                "SET" if se_rf & SEFLG_BARYCTR else "no",
                "SET" if se_rf & SEFLG_MOSEPH else "no",
                "SET" if se_rf & SEFLG_HELCTR else "no",
            )
        )
    else:
        print("    swe retflag={}".format(se_rf))


# ============================================================================
# PART 6: Earth body (SE_EARTH=14) in barycentric mode
# ============================================================================

print()
print("=" * 100)
print("PART 6: Earth (body=14) in barycentric mode")
print("=" * 100)

for y, m, d in [(2000, 1, 1), (2010, 1, 1), (2020, 1, 1)]:
    jd = year_to_jd(y, m, d)
    print()
    print("  {}-{:02d}-{:02d}:".format(y, m, d))

    try:
        le = ephem.swe_calc_ut(jd, 14, BARY_FLAG)
        print(
            "    lib Earth bary: lon={:12.6f}  lat={:12.6f}  dist={:12.8f}  retflag={}".format(
                le[0][0], le[0][1], le[0][2], le[1]
            )
        )
    except Exception as e:
        print("    lib Earth bary: ERROR - {}".format(e))

    try:
        se = swe.calc_ut(jd, 14, BARY_FLAG)
        print(
            "    swe Earth bary: lon={:12.6f}  lat={:12.6f}  dist={:12.8f}  retflag={}".format(
                se[0][0], se[0][1], se[0][2], se[1]
            )
        )
    except Exception as e:
        print("    swe Earth bary: ERROR - {}".format(e))


# ============================================================================
# PART 7: Other bodies in barycentric - Sun-specific?
# ============================================================================

print()
print("=" * 100)
print(
    "PART 7: All major bodies in REAL barycentric mode (flag={}) at J2000".format(
        SEFLG_BARYCTR
    )
)
print("        Confirm whether discrepancy is Sun-specific")
print("=" * 100)

jd = year_to_jd(2000)
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

print(
    "  {:>10s}  {:>10s}  {:>10s}  {:>14s}  {:>12s}  {:>12s}  {:>10s}".format(
        "Body", "dLon asec", "dLat asec", "dDist AU", "LE dist", "SE dist", "Status"
    )
)
print("  " + "-" * 90)

for body_id, name in bodies:
    try:
        le = ephem.swe_calc_ut(jd, body_id, BARY_FLAG)[0]
        se = swe.calc_ut(jd, body_id, BARY_FLAG)[0]
        dlon = angular_diff(le[0], se[0]) * 3600
        dlat = (le[1] - se[1]) * 3600
        ddist = le[2] - se[2]
        status = "OK" if abs(dlon) < 3.6 else "PROBLEM"  # 3.6" = 0.001 deg
        print(
            "  {:>10s}  {:>+10.3f}  {:>+10.3f}  {:>+14.8f}  "
            "{:>12.6f}  {:>12.6f}  {:>10s}".format(
                name, dlon, dlat, ddist, le[2], se[2], status
            )
        )
    except Exception as e:
        print("  {:>10s}  ERROR: {}".format(name, e))


# ============================================================================
# PART 8: Also test with MOSEPH flag on other bodies
# ============================================================================

print()
print("=" * 100)
print("PART 8: Other bodies with flag=4 (MOSEPH) at J2000 - for comparison")
print("        If all bodies show similar discrepancy, confirms MOSEPH vs JPL issue")
print("=" * 100)

print(
    "  {:>10s}  {:>10s}  {:>14s}  {:>10s}".format(
        "Body", "dLon asec", "dDist AU", "Status"
    )
)
print("  " + "-" * 60)

for body_id, name in bodies:
    try:
        le = ephem.swe_calc_ut(jd, body_id, WRONG_FLAG)[0]
        se = swe.calc_ut(jd, body_id, WRONG_FLAG)[0]
        dlon = angular_diff(le[0], se[0]) * 3600
        ddist = le[2] - se[2]
        status = "OK" if abs(dlon) < 3.6 else "MOSEPH diff"
        print(
            "  {:>10s}  {:>+10.3f}  {:>+14.8f}  {:>10s}".format(
                name, dlon, ddist, status
            )
        )
    except Exception as e:
        print("  {:>10s}  ERROR: {}".format(name, e))


# ============================================================================
# SUMMARY
# ============================================================================

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)
print()
print("  Previous audit used FLG_BARYCTR = 4, which is actually SEFLG_MOSEPH.")
print(
    '  Max |dLon| with flag=4  (MOSEPH):  {:.3f}" ({:.6f} deg)'.format(
        max_dlon_wrong, max_dlon_wrong / 3600
    )
)
print(
    '  Max |dLon| with flag={} (BARYCTR): {:.3f}" ({:.6f} deg)'.format(
        SEFLG_BARYCTR, max_dlon_correct, max_dlon_correct / 3600
    )
)
print()
print('  If the MOSEPH discrepancy is large (~100+") and the BARYCTR discrepancy is')
print("  small (<5\"), then the 'bug' was simply using the wrong flag constant.")
print(
    "  pyswisseph with flag=4 uses its Moshier analytical ephemeris (lower precision),"
)
print("  while libephemeris ignores MOSEPH and always uses JPL DE440 (high precision).")
print()
