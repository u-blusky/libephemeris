"""
Deep comparative analysis round 2: libephemeris vs pyswisseph.

Focuses on areas not covered in round 1:
- NONUT + equatorial combinations
- NONUT for lunar bodies (MeanNode, MeanApog, TrueNode, OscuApog)
- Velocity precision for all bodies
- Planet-centric calculations (calc_pctr)
- SEFLG_NOGDEFL for individual planets
- Pholus/minor body precision deep-dive
- Eclipse location-specific functions
- Houses at extreme latitudes
- Sidereal + NONUT combinations
"""

from __future__ import annotations

import math
import os
import random
import statistics

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
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
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
    SE_INTP_APOG,
    SE_INTP_PERG,
    SE_CHIRON,
    SE_PHOLUS,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
    SEFLG_SWIEPH,
    SEFLG_SPEED,
    SEFLG_HELCTR,
    SEFLG_EQUATORIAL,
    SEFLG_J2000,
    SEFLG_NONUT,
    SEFLG_TRUEPOS,
    SEFLG_NOABERR,
    SEFLG_NOGDEFL,
    SEFLG_XYZ,
    SEFLG_RADIANS,
    SEFLG_SIDEREAL,
    SEFLG_BARYCTR,
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_LAHIRI,
)

swe.set_ephe_path("swisseph/ephe")
random.seed(42)

PLANETS = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
    (SE_URANUS, "Uranus"),
    (SE_NEPTUNE, "Neptune"),
    (SE_PLUTO, "Pluto"),
]

LUNAR_BODIES = [
    (SE_MEAN_NODE, "MeanNode"),
    (SE_TRUE_NODE, "TrueNode"),
    (SE_MEAN_APOG, "MeanApog"),
    (SE_OSCU_APOG, "OscuApog"),
    (SE_INTP_APOG, "IntpApog"),
    (SE_INTP_PERG, "IntpPerg"),
]

MINOR_BODIES = [
    (SE_CHIRON, "Chiron"),
    (SE_PHOLUS, "Pholus"),
    (SE_CERES, "Ceres"),
    (SE_PALLAS, "Pallas"),
    (SE_JUNO, "Juno"),
    (SE_VESTA, "Vesta"),
]


def angular_diff(a, b):
    d = abs(a - b) % 360
    return min(d, 360 - d)


def arcsec(deg):
    return deg * 3600


def gen_jds(n=100, start=2415020.5, end=2488069.5):
    return [random.uniform(start, end) for _ in range(n)]


issues = []


def report(cat, sev, desc, val, unit="arcsec"):
    issues.append({"cat": cat, "sev": sev, "desc": desc, "val": val, "unit": unit})


# ============================================================================
# 1. NONUT + EQUATORIAL for all planets
# ============================================================================
print("=" * 70)
print("1. NONUT + EQUATORIAL (200 dates x all planets)")
print("=" * 70)

jds = gen_jds(200)
flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_NONUT | SEFLG_EQUATORIAL

for pid, pname in PLANETS:
    lon_errs = []
    lat_errs = []
    for jd in jds:
        try:
            r_swe = swe.calc_ut(jd, pid, flags)
            r_lib = ephem.swe_calc_ut(jd, pid, flags)
            lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
            lat_errs.append(abs(r_swe[0][1] - r_lib[0][1]))
        except Exception:
            continue
    if lon_errs:
        mx_lon = max(lon_errs)
        mx_lat = max(lat_errs)
        if arcsec(mx_lon) > 0.5 or arcsec(mx_lat) > 0.5:
            print(
                f'  {pname:10s} RA: max={arcsec(mx_lon):8.3f}"  Dec: max={arcsec(mx_lat):8.3f}"'
            )
            if arcsec(mx_lon) > 5:
                report(
                    "nonut_eq",
                    "HIGH",
                    f'{pname} NONUT+EQ RA max={arcsec(mx_lon):.1f}"',
                    arcsec(mx_lon),
                )

# ============================================================================
# 2. NONUT for lunar bodies
# ============================================================================
print("\n" + "=" * 70)
print("2. NONUT for LUNAR BODIES (300 dates)")
print("=" * 70)

jds_l = gen_jds(300)
for flags_desc, flags_val in [
    ("NONUT ecliptic", SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_NONUT),
    ("NONUT+equatorial", SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_NONUT | SEFLG_EQUATORIAL),
]:
    print(f"  [{flags_desc}]")
    for bid, bname in LUNAR_BODIES:
        lon_errs = []
        for jd in jds_l:
            try:
                r_swe = swe.calc_ut(jd, bid, flags_val)
                r_lib = ephem.swe_calc_ut(jd, bid, flags_val)
                lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
            except Exception:
                continue
        if lon_errs:
            mx = max(lon_errs)
            mean = statistics.mean(lon_errs)
            print(f'    {bname:10s} max={arcsec(mx):8.3f}" mean={arcsec(mean):8.3f}"')
            if arcsec(mx) > 5 and bname not in ("IntpApog", "IntpPerg"):
                report(
                    "nonut_lunar",
                    "HIGH",
                    f'{bname} {flags_desc} max={arcsec(mx):.1f}"',
                    arcsec(mx),
                )

# ============================================================================
# 3. VELOCITY PRECISION for lunar bodies
# ============================================================================
print("\n" + "=" * 70)
print("3. VELOCITY PRECISION - lunar bodies (300 dates)")
print("=" * 70)

flags = SEFLG_SWIEPH | SEFLG_SPEED
for bid, bname in LUNAR_BODIES:
    speed_errs = []
    for jd in jds_l:
        try:
            r_swe = swe.calc_ut(jd, bid, flags)
            r_lib = ephem.swe_calc_ut(jd, bid, flags)
            if len(r_swe[0]) > 3 and len(r_lib[0]) > 3:
                speed_errs.append(abs(r_swe[0][3] - r_lib[0][3]))
        except Exception:
            continue
    if speed_errs:
        mx = max(speed_errs)
        mean = statistics.mean(speed_errs)
        if mx > 0.0001:
            print(f"  {bname:10s} lon_speed: max={mx:.8f} deg/day mean={mean:.8f}")
            if mx > 0.01:
                report(
                    "velocity_lunar",
                    "MEDIUM",
                    f"{bname} speed max={mx:.6f} deg/day",
                    mx,
                    "deg/day",
                )

# ============================================================================
# 4. BARYCENTRIC positions
# ============================================================================
print("\n" + "=" * 70)
print("4. BARYCENTRIC positions (100 dates)")
print("=" * 70)

jds_b = gen_jds(100)
flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_BARYCTR

for pid, pname in PLANETS:
    if pid == SE_MOON:
        continue  # barycentric Moon is special
    lon_errs = []
    for jd in jds_b:
        try:
            r_swe = swe.calc_ut(jd, pid, flags)
            r_lib = ephem.swe_calc_ut(jd, pid, flags)
            lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
        except Exception:
            continue
    if lon_errs:
        mx = max(lon_errs)
        if arcsec(mx) > 1.0:
            print(f'  {pname:10s} max={arcsec(mx):8.3f}"')
            if arcsec(mx) > 5:
                report(
                    "barycentric",
                    "MEDIUM",
                    f'{pname} bary max={arcsec(mx):.1f}"',
                    arcsec(mx),
                )

# ============================================================================
# 5. XYZ Cartesian output
# ============================================================================
print("\n" + "=" * 70)
print("5. XYZ CARTESIAN output (100 dates)")
print("=" * 70)

jds_x = gen_jds(100)
flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_XYZ

for pid, pname in PLANETS:
    xyz_errs = []
    for jd in jds_x:
        try:
            r_swe = swe.calc_ut(jd, pid, flags)
            r_lib = ephem.swe_calc_ut(jd, pid, flags)
            s, l = r_swe[0], r_lib[0]
            for i in range(3):
                if abs(s[i]) > 0.001:
                    xyz_errs.append(abs(s[i] - l[i]))
        except Exception:
            continue
    if xyz_errs:
        mx = max(xyz_errs)
        if mx > 1e-5:
            print(f"  {pname:10s} max_xyz_diff={mx:.10f} AU")
            if mx > 1e-3:
                report("xyz", "HIGH", f"{pname} XYZ max={mx:.6f} AU", mx, "AU")

# ============================================================================
# 6. RADIANS output
# ============================================================================
print("\n" + "=" * 70)
print("6. RADIANS output (50 dates)")
print("=" * 70)

jds_r = gen_jds(50)
flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_RADIANS

for pid, pname in [(SE_SUN, "Sun"), (SE_MOON, "Moon"), (SE_MARS, "Mars")]:
    rad_errs = []
    for jd in jds_r:
        try:
            r_swe = swe.calc_ut(jd, pid, flags)
            r_lib = ephem.swe_calc_ut(jd, pid, flags)
            s, l = r_swe[0], r_lib[0]
            # Compare radians (lon in radians)
            diff = abs(s[0] - l[0])
            if diff > math.pi:
                diff = 2 * math.pi - diff
            rad_errs.append(diff)
        except Exception:
            continue
    if rad_errs:
        mx = max(rad_errs)
        mx_arcsec = math.degrees(mx) * 3600
        if mx_arcsec > 0.5:
            print(f'  {pname:10s} max_rad_diff={mx:.10f} rad ({mx_arcsec:.3f}")')
            if mx_arcsec > 5:
                report(
                    "radians",
                    "MEDIUM",
                    f'{pname} radians max={mx_arcsec:.1f}"',
                    mx_arcsec,
                )

# ============================================================================
# 7. J2000 + EQUATORIAL combination
# ============================================================================
print("\n" + "=" * 70)
print("7. J2000 + EQUATORIAL (100 dates)")
print("=" * 70)

jds_j = gen_jds(100)
flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_J2000 | SEFLG_EQUATORIAL

for pid, pname in PLANETS:
    lon_errs = []
    for jd in jds_j:
        try:
            r_swe = swe.calc_ut(jd, pid, flags)
            r_lib = ephem.swe_calc_ut(jd, pid, flags)
            lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
        except Exception:
            continue
    if lon_errs:
        mx = max(lon_errs)
        if arcsec(mx) > 1.0:
            print(f'  {pname:10s} RA max={arcsec(mx):8.3f}"')
            if arcsec(mx) > 5:
                report(
                    "j2000_eq",
                    "MEDIUM",
                    f'{pname} J2000+EQ max={arcsec(mx):.1f}"',
                    arcsec(mx),
                )

# ============================================================================
# 8. NONUT + J2000 (should be equivalent to J2000)
# ============================================================================
print("\n" + "=" * 70)
print("8. NONUT + J2000 (100 dates, should equal J2000)")
print("=" * 70)

flags_nonut_j2k = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_NONUT | SEFLG_J2000

for pid, pname in [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
]:
    lon_errs = []
    for jd in jds_j:
        try:
            r_swe = swe.calc_ut(jd, pid, flags_nonut_j2k)
            r_lib = ephem.swe_calc_ut(jd, pid, flags_nonut_j2k)
            lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
        except Exception:
            continue
    if lon_errs:
        mx = max(lon_errs)
        if arcsec(mx) > 0.5:
            print(f'  {pname:10s} max={arcsec(mx):8.3f}"')

# ============================================================================
# 9. HELIOCENTRIC + various flags
# ============================================================================
print("\n" + "=" * 70)
print("9. HELIOCENTRIC + EQUATORIAL (100 dates)")
print("=" * 70)

jds_h = gen_jds(100)
flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR | SEFLG_EQUATORIAL

for pid, pname in PLANETS:
    if pid in (SE_SUN, SE_MOON):
        continue
    lon_errs = []
    for jd in jds_h:
        try:
            r_swe = swe.calc_ut(jd, pid, flags)
            r_lib = ephem.swe_calc_ut(jd, pid, flags)
            lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
        except Exception:
            continue
    if lon_errs:
        mx = max(lon_errs)
        if arcsec(mx) > 1.0:
            print(f'  {pname:10s} max={arcsec(mx):8.3f}"')
            if arcsec(mx) > 5:
                report(
                    "helio_eq",
                    "MEDIUM",
                    f'{pname} HELCTR+EQ max={arcsec(mx):.1f}"',
                    arcsec(mx),
                )

# ============================================================================
# 10. Minor body velocity deep-dive
# ============================================================================
print("\n" + "=" * 70)
print("10. MINOR BODY VELOCITY (200 dates)")
print("=" * 70)

jds_m = gen_jds(200)
flags = SEFLG_SWIEPH | SEFLG_SPEED
for bid, bname in MINOR_BODIES:
    speed_lon = []
    speed_lat = []
    speed_dist = []
    for jd in jds_m:
        try:
            r_swe = swe.calc_ut(jd, bid, flags)
            r_lib = ephem.swe_calc_ut(jd, bid, flags)
            s, l = r_swe[0], r_lib[0]
            if len(s) > 5 and len(l) > 5:
                speed_lon.append(abs(s[3] - l[3]))
                speed_lat.append(abs(s[4] - l[4]))
                speed_dist.append(abs(s[5] - l[5]))
        except Exception:
            continue
    if speed_lon:
        mx_lon = max(speed_lon)
        mx_lat = max(speed_lat)
        mx_dist = max(speed_dist)
        if mx_lon > 0.0001 or mx_lat > 0.0001:
            print(
                f"  {bname:10s} lon: {mx_lon:.8f}  lat: {mx_lat:.8f}  dist: {mx_dist:.10f}"
            )
            if mx_lon > 0.01:
                report(
                    "minor_vel",
                    "MEDIUM",
                    f"{bname} speed max={mx_lon:.6f}",
                    mx_lon,
                    "deg/day",
                )

# ============================================================================
# 11. Pholus position deep-dive (was 24\" in round 1)
# ============================================================================
print("\n" + "=" * 70)
print("11. PHOLUS POSITION DEEP-DIVE (500 dates)")
print("=" * 70)

jds_p = gen_jds(500)
flags = SEFLG_SWIEPH | SEFLG_SPEED
pholus_errs = []
worst_jd = None
worst_diff = 0
for jd in jds_p:
    try:
        r_swe = swe.calc_ut(jd, SE_PHOLUS, flags)
        r_lib = ephem.swe_calc_ut(jd, SE_PHOLUS, flags)
        diff = angular_diff(r_swe[0][0], r_lib[0][0])
        pholus_errs.append(diff)
        if diff > worst_diff:
            worst_diff = diff
            worst_jd = jd
    except Exception:
        continue

if pholus_errs:
    mx = max(pholus_errs)
    mean = statistics.mean(pholus_errs)
    p95 = sorted(pholus_errs)[int(len(pholus_errs) * 0.95)]
    print(
        f'  Pholus: max={arcsec(mx):.3f}"  mean={arcsec(mean):.3f}"  p95={arcsec(p95):.3f}"'
    )
    print(f"  Worst at JD={worst_jd:.1f}")

# ============================================================================
# 12. Sidereal + various combinations
# ============================================================================
print("\n" + "=" * 70)
print("12. SIDEREAL + NONUT (Lahiri, 100 dates)")
print("=" * 70)

jds_s = gen_jds(100)
swe.set_sid_mode(SE_SIDM_LAHIRI)
ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

for flags_desc, flags_val in [
    ("SID+NONUT", SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_SIDEREAL | SEFLG_NONUT),
    ("SID+EQ", SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_SIDEREAL | SEFLG_EQUATORIAL),
    (
        "SID+NONUT+EQ",
        SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_SIDEREAL | SEFLG_NONUT | SEFLG_EQUATORIAL,
    ),
]:
    print(f"  [{flags_desc}]")
    for pid, pname in [(SE_SUN, "Sun"), (SE_MOON, "Moon"), (SE_MARS, "Mars")]:
        lon_errs = []
        for jd in jds_s:
            try:
                r_swe = swe.calc_ut(jd, pid, flags_val)
                r_lib = ephem.swe_calc_ut(jd, pid, flags_val)
                lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
            except Exception:
                continue
        if lon_errs:
            mx = max(lon_errs)
            if arcsec(mx) > 0.5:
                print(f'    {pname:10s} max={arcsec(mx):8.3f}"')
                if arcsec(mx) > 10:
                    report(
                        "sid_combo",
                        "HIGH",
                        f'{pname} {flags_desc} max={arcsec(mx):.1f}"',
                        arcsec(mx),
                    )

swe.set_sid_mode(0)
ephem.swe_set_sid_mode(0)

# ============================================================================
# 13. calc_pctr (planet-centric)
# ============================================================================
print("\n" + "=" * 70)
print("13. PLANET-CENTRIC calc_pctr (50 dates)")
print("=" * 70)

jds_pc = gen_jds(50)
# Test: Moon as seen from Mars
for center_id, center_name in [(SE_MARS, "Mars"), (SE_JUPITER, "Jupiter")]:
    for target_id, target_name in [
        (SE_MOON, "Moon"),
        (SE_SUN, "Sun"),
        (SE_VENUS, "Venus"),
    ]:
        lon_errs = []
        for jd in jds_pc:
            try:
                r_swe = swe.calc_pctr(
                    jd, target_id, center_id, SEFLG_SWIEPH | SEFLG_SPEED
                )
                r_lib = ephem.swe_calc_pctr(
                    jd, target_id, center_id, SEFLG_SWIEPH | SEFLG_SPEED
                )
                lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
            except Exception:
                continue
        if lon_errs:
            mx = max(lon_errs)
            if arcsec(mx) > 1.0:
                print(
                    f'  {target_name:8s} from {center_name:8s}: max={arcsec(mx):8.3f}"'
                )
                if arcsec(mx) > 10:
                    report(
                        "pctr",
                        "MEDIUM",
                        f'{target_name} from {center_name} max={arcsec(mx):.1f}"',
                        arcsec(mx),
                    )

# ============================================================================
# 14. Solar eclipse WHERE (geographic)
# ============================================================================
print("\n" + "=" * 70)
print("14. SOLAR ECLIPSE WHERE (geographic details)")
print("=" * 70)

jd = 2451545.0
n_eclipses = 0
ecl_attr_errs = []
for _ in range(20):
    try:
        r_swe = swe.sol_eclipse_when_glob(jd, SEFLG_SWIEPH)
        r_lib = ephem.swe_sol_eclipse_when_glob(jd, SEFLG_SWIEPH)
        ecl_jd = r_swe[1][0]
        # Now test sol_eclipse_where at the max eclipse time
        try:
            w_swe = swe.sol_eclipse_where(ecl_jd, SEFLG_SWIEPH)
            w_lib = ephem.swe_sol_eclipse_where(ecl_jd, SEFLG_SWIEPH)
            # Compare geographic coords (geopos) and attributes
            if w_swe and w_lib:
                geo_swe = w_swe[1]  # (lon, lat, ...)
                geo_lib = w_lib[1]
                if len(geo_swe) >= 2 and len(geo_lib) >= 2:
                    geo_diff = math.sqrt(
                        (geo_swe[0] - geo_lib[0]) ** 2 + (geo_swe[1] - geo_lib[1]) ** 2
                    )
                    ecl_attr_errs.append(geo_diff)
                    n_eclipses += 1
        except Exception:
            pass
        jd = ecl_jd + 1
    except Exception:
        jd += 30

if ecl_attr_errs:
    mx = max(ecl_attr_errs)
    mean = statistics.mean(ecl_attr_errs)
    print(
        f"  sol_eclipse_where geo diff: max={mx:.4f}° mean={mean:.4f}° (n={n_eclipses})"
    )
    if mx > 1.0:
        report(
            "ecl_where", "MEDIUM", f"sol_eclipse_where geo max={mx:.2f}°", mx, "degrees"
        )

# ============================================================================
# 15. houses_ex (with SEFLG_SPEED)
# ============================================================================
print("\n" + "=" * 70)
print("15. HOUSES_EX cusp speeds (50 dates, key systems)")
print("=" * 70)

jds_he = gen_jds(50)
for sys_char in "PKOB":
    cusp_speed_errs = []
    for jd in jds_he:
        try:
            r_swe = swe.houses_ex(
                jd, 41.9, 12.5, ord(sys_char), SEFLG_SWIEPH | SEFLG_SPEED
            )
            r_lib = ephem.swe_houses_ex(
                jd, 41.9, 12.5, ord(sys_char), SEFLG_SWIEPH | SEFLG_SPEED
            )
            if r_swe and r_lib:
                cusps_swe, _, cusp_speeds_swe, _ = r_swe
                cusps_lib, _, cusp_speeds_lib, _ = r_lib
                for i in range(min(len(cusp_speeds_swe), len(cusp_speeds_lib))):
                    if cusp_speeds_swe[i] != 0 or cusp_speeds_lib[i] != 0:
                        cusp_speed_errs.append(
                            abs(cusp_speeds_swe[i] - cusp_speeds_lib[i])
                        )
        except Exception:
            continue
    if cusp_speed_errs:
        mx = max(cusp_speed_errs)
        mean = statistics.mean(cusp_speed_errs)
        if mx > 0.01:
            print(
                f"  System {sys_char}: cusp_speed max_diff={mx:.6f} deg/day mean={mean:.6f}"
            )

# ============================================================================
# FINAL REPORT
# ============================================================================
print("\n" + "=" * 70)
print("ROUND 2 FINAL REPORT: ACTIONABLE ISSUES")
print("=" * 70)

if not issues:
    print("  No actionable issues found!")
else:
    severity_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
    issues.sort(key=lambda x: (severity_order.get(x["sev"], 3), -x["val"]))
    for i in issues:
        print(
            f"  [{i['sev']:6s}] {i['cat']:20s} | {i['desc']} ({i['val']:.3f} {i['unit']})"
        )

print(
    f"\nTotal issues: {len(issues)} (HIGH={sum(1 for i in issues if i['sev'] == 'HIGH')}, MEDIUM={sum(1 for i in issues if i['sev'] == 'MEDIUM')})"
)
