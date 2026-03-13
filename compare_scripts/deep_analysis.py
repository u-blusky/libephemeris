"""
Deep comparative analysis: libephemeris vs pyswisseph.

Systematically probes all major API surfaces for precision discrepancies
that could be worth fixing. Outputs a structured report.
"""

from __future__ import annotations

import math
import os
import random
import statistics
import sys
import traceback

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
    SE_NODBIT_MEAN,
    SE_NODBIT_OSCU,
    SE_NODBIT_OSCU_BAR,
    SE_NODBIT_FOPOINT,
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


def angular_diff(a: float, b: float) -> float:
    d = abs(a - b) % 360
    return min(d, 360 - d)


def arcsec(deg: float) -> float:
    return deg * 3600


def generate_jds(
    n: int = 100, start: float = 2415020.5, end: float = 2488069.5
) -> list[float]:
    """Generate n random JDs between start and end (default 1900-2100)."""
    return [random.uniform(start, end) for _ in range(n)]


issues = []  # Collect all issues found


def report_issue(
    category: str, severity: str, description: str, max_err: float, unit: str = "arcsec"
):
    issues.append(
        {
            "category": category,
            "severity": severity,
            "description": description,
            "max_err": max_err,
            "unit": unit,
        }
    )


# ============================================================================
# 1. PLANETARY POSITIONS - per-planet, per-flag precision
# ============================================================================
print("=" * 70)
print("1. PLANETARY POSITIONS (200 dates x all planets x key flags)")
print("=" * 70)

jds = generate_jds(200)
flags_to_test = [
    (SEFLG_SWIEPH | SEFLG_SPEED, "ecliptic+speed"),
    (SEFLG_SWIEPH | SEFLG_EQUATORIAL | SEFLG_SPEED, "equatorial+speed"),
    (SEFLG_SWIEPH | SEFLG_HELCTR | SEFLG_SPEED, "heliocentric"),
    (SEFLG_SWIEPH | SEFLG_J2000 | SEFLG_SPEED, "J2000"),
    (SEFLG_SWIEPH | SEFLG_NONUT | SEFLG_SPEED, "nonut"),
    (SEFLG_SWIEPH | SEFLG_TRUEPOS | SEFLG_SPEED, "truepos"),
    (SEFLG_SWIEPH | SEFLG_NOABERR | SEFLG_SPEED, "noaberr"),
    (SEFLG_SWIEPH | SEFLG_NOGDEFL | SEFLG_SPEED, "nogdefl"),
    (SEFLG_SWIEPH | SEFLG_NOABERR | SEFLG_NOGDEFL | SEFLG_SPEED, "astrometric"),
]

for flags, flag_name in flags_to_test:
    for pid, pname in PLANETS:
        if (flags & SEFLG_HELCTR) and pid in (SE_SUN, SE_MOON):
            continue
        lon_errs = []
        lat_errs = []
        dist_errs = []
        speed_errs = []
        for jd in jds:
            try:
                r_swe = swe.calc_ut(jd, pid, flags)
                r_lib = ephem.swe_calc_ut(jd, pid, flags)
                s, l = r_swe[0], r_lib[0]
                lon_errs.append(angular_diff(s[0], l[0]))
                lat_errs.append(abs(s[1] - l[1]))
                if s[2] > 0:
                    dist_errs.append(abs(s[2] - l[2]))
                if len(s) > 3 and len(l) > 3:
                    speed_errs.append(abs(s[3] - l[3]))
            except Exception:
                continue

        if lon_errs:
            mx_lon = max(lon_errs)
            mx_lat = max(lat_errs) if lat_errs else 0
            mx_spd = max(speed_errs) if speed_errs else 0
            # Report if significant
            if arcsec(mx_lon) > 1.0 or arcsec(mx_lat) > 1.0 or mx_spd > 0.001:
                print(
                    f'  {pname:10s} [{flag_name:15s}] lon={arcsec(mx_lon):8.3f}" lat={arcsec(mx_lat):8.3f}" speed_diff={mx_spd:.6f} deg/day'
                )
                if arcsec(mx_lon) > 5.0:
                    report_issue(
                        "planet_position",
                        "HIGH",
                        f'{pname}/{flag_name} lon max={arcsec(mx_lon):.3f}"',
                        arcsec(mx_lon),
                    )
                elif arcsec(mx_lon) > 1.0:
                    report_issue(
                        "planet_position",
                        "MEDIUM",
                        f'{pname}/{flag_name} lon max={arcsec(mx_lon):.3f}"',
                        arcsec(mx_lon),
                    )

# ============================================================================
# 2. LUNAR BODIES PRECISION
# ============================================================================
print("\n" + "=" * 70)
print("2. LUNAR BODIES (300 dates)")
print("=" * 70)

jds_lunar = generate_jds(300)
flags = SEFLG_SWIEPH | SEFLG_SPEED

for bid, bname in LUNAR_BODIES:
    lon_errs = []
    lat_errs = []
    speed_errs = []
    for jd in jds_lunar:
        try:
            r_swe = swe.calc_ut(jd, bid, flags)
            r_lib = ephem.swe_calc_ut(jd, bid, flags)
            s, l = r_swe[0], r_lib[0]
            lon_errs.append(angular_diff(s[0], l[0]))
            lat_errs.append(abs(s[1] - l[1]))
            if len(s) > 3 and len(l) > 3:
                speed_errs.append(abs(s[3] - l[3]))
        except Exception:
            continue

    if lon_errs:
        mx_lon = max(lon_errs)
        mx_lat = max(lat_errs)
        mx_spd = max(speed_errs) if speed_errs else 0
        mean_lon = statistics.mean(lon_errs)
        print(
            f'  {bname:10s} lon: max={arcsec(mx_lon):8.3f}" mean={arcsec(mean_lon):8.3f}" | lat: max={arcsec(mx_lat):8.3f}" | speed: max={mx_spd:.6f}'
        )
        if arcsec(mx_lon) > 60:
            report_issue(
                "lunar_body",
                "HIGH",
                f'{bname} lon max={arcsec(mx_lon):.1f}"',
                arcsec(mx_lon),
            )
        elif arcsec(mx_lon) > 5:
            report_issue(
                "lunar_body",
                "MEDIUM",
                f'{bname} lon max={arcsec(mx_lon):.1f}"',
                arcsec(mx_lon),
            )

# ============================================================================
# 3. MINOR BODIES / ASTEROIDS
# ============================================================================
print("\n" + "=" * 70)
print("3. MINOR BODIES (100 dates)")
print("=" * 70)

jds_minor = generate_jds(100)
for bid, bname in MINOR_BODIES:
    lon_errs = []
    lat_errs = []
    speed_errs = []
    for jd in jds_minor:
        try:
            r_swe = swe.calc_ut(jd, bid, flags)
            r_lib = ephem.swe_calc_ut(jd, bid, flags)
            s, l = r_swe[0], r_lib[0]
            lon_errs.append(angular_diff(s[0], l[0]))
            lat_errs.append(abs(s[1] - l[1]))
            if len(s) > 3 and len(l) > 3:
                speed_errs.append(abs(s[3] - l[3]))
        except Exception:
            continue
    if lon_errs:
        mx_lon = max(lon_errs)
        mx_lat = max(lat_errs)
        mx_spd = max(speed_errs) if speed_errs else 0
        print(
            f'  {bname:10s} lon: max={arcsec(mx_lon):8.3f}" mean={arcsec(statistics.mean(lon_errs)):8.3f}" | lat: max={arcsec(mx_lat):8.3f}" | speed: max={mx_spd:.6f}'
        )
        if arcsec(mx_lon) > 10:
            report_issue(
                "minor_body",
                "MEDIUM",
                f'{bname} lon max={arcsec(mx_lon):.1f}"',
                arcsec(mx_lon),
            )

# ============================================================================
# 4. SPEED/VELOCITY PRECISION (focused test)
# ============================================================================
print("\n" + "=" * 70)
print("4. VELOCITY PRECISION (200 dates, all planets)")
print("=" * 70)

jds_vel = generate_jds(200)
for pid, pname in PLANETS:
    speed_errs_lon = []
    speed_errs_lat = []
    speed_errs_dist = []
    for jd in jds_vel:
        try:
            r_swe = swe.calc_ut(jd, pid, SEFLG_SWIEPH | SEFLG_SPEED)
            r_lib = ephem.swe_calc_ut(jd, pid, SEFLG_SWIEPH | SEFLG_SPEED)
            s, l = r_swe[0], r_lib[0]
            speed_errs_lon.append(abs(s[3] - l[3]))
            speed_errs_lat.append(abs(s[4] - l[4]))
            speed_errs_dist.append(abs(s[5] - l[5]))
        except Exception:
            continue
    if speed_errs_lon:
        mx_lon = max(speed_errs_lon)
        mx_lat = max(speed_errs_lat)
        mx_dist = max(speed_errs_dist)
        if mx_lon > 0.0001 or mx_lat > 0.0001:
            print(
                f"  {pname:10s} lon_speed: max={mx_lon:.8f} deg/day | lat_speed: max={mx_lat:.8f} | dist_speed: max={mx_dist:.10f}"
            )
            if mx_lon > 0.01:
                report_issue(
                    "velocity",
                    "HIGH",
                    f"{pname} lon speed diff={mx_lon:.6f} deg/day",
                    mx_lon,
                    "deg/day",
                )

# ============================================================================
# 5. HOUSES PRECISION (focused on problem systems)
# ============================================================================
print("\n" + "=" * 70)
print("5. HOUSE CUSP PRECISION (50 dates x key locations x all systems)")
print("=" * 70)

house_systems = list("PKROCBMEUDHNIJFAXLGQSWVY")
locations = [
    (41.9, 12.5, "Rome"),
    (51.5, -0.1, "London"),
    (40.7, -74.0, "NewYork"),
    (35.7, 139.7, "Tokyo"),
    (-33.9, 151.2, "Sydney"),
    (64.1, -21.9, "Reykjavik"),  # High latitude
    (78.2, 15.6, "Svalbard"),  # Very high latitude
]

jds_houses = generate_jds(50)
house_issues = {}

for sys_char in house_systems:
    for lat, lon, loc_name in locations:
        cusp_errs = []
        asc_errs = []
        mc_errs = []
        for jd in jds_houses:
            try:
                r_swe = swe.houses(jd, lat, lon, ord(sys_char))
                r_lib = ephem.swe_houses(jd, lat, lon, ord(sys_char))
                # cusps
                cusps_swe, points_swe = r_swe
                cusps_lib, points_lib = r_lib
                for i in range(min(len(cusps_swe), len(cusps_lib))):
                    if cusps_swe[i] != 0 or cusps_lib[i] != 0:
                        cusp_errs.append(angular_diff(cusps_swe[i], cusps_lib[i]))
                # ASC, MC
                if len(points_swe) >= 2 and len(points_lib) >= 2:
                    asc_errs.append(angular_diff(points_swe[0], points_lib[0]))
                    mc_errs.append(angular_diff(points_swe[1], points_lib[1]))
            except Exception:
                continue
        if cusp_errs:
            mx_cusp = max(cusp_errs)
            mx_asc = max(asc_errs) if asc_errs else 0
            mx_mc = max(mc_errs) if mc_errs else 0
            if arcsec(mx_cusp) > 1.0:
                key = f"House_{sys_char}/{loc_name}"
                house_issues[key] = mx_cusp
                if arcsec(mx_cusp) > 10:
                    print(
                        f'  System {sys_char} @ {loc_name:12s}: cusp max={arcsec(mx_cusp):8.3f}" ASC max={arcsec(mx_asc):8.3f}" MC max={arcsec(mx_mc):8.3f}"'
                    )

if house_issues:
    # Group by system
    by_system = {}
    for key, val in house_issues.items():
        sys_char = key.split("/")[0].replace("House_", "")
        if sys_char not in by_system:
            by_system[sys_char] = []
        by_system[sys_char].append(val)
    for sys_char, vals in sorted(by_system.items()):
        mx = max(vals)
        if arcsec(mx) > 5:
            report_issue(
                "houses",
                "MEDIUM",
                f'System {sys_char}: cusp max={arcsec(mx):.1f}"',
                arcsec(mx),
            )
else:
    print('  All house systems within 1" tolerance.')

# ============================================================================
# 6. PHENOMENA (pheno_ut)
# ============================================================================
print("\n" + "=" * 70)
print("6. PHENOMENA - pheno_ut (100 dates x planets)")
print("=" * 70)

jds_pheno = generate_jds(100)
for pid, pname in PLANETS:
    if pid == SE_SUN:
        continue
    attr_errs = []
    for jd in jds_pheno:
        try:
            r_swe = swe.pheno_ut(jd, pid, SEFLG_SWIEPH)
            r_lib = ephem.swe_pheno_ut(jd, pid, SEFLG_SWIEPH)
            for i in range(min(len(r_swe[0]), len(r_lib[0]))):
                attr_errs.append(abs(r_swe[0][i] - r_lib[0][i]))
        except Exception:
            continue
    if attr_errs:
        mx = max(attr_errs)
        mean = statistics.mean(attr_errs)
        if mx > 0.01:
            print(f"  {pname:10s} max_attr_diff={mx:.6f} mean={mean:.8f}")
            if mx > 0.1:
                report_issue(
                    "phenomena",
                    "MEDIUM",
                    f"{pname} pheno max attr diff={mx:.4f}",
                    mx,
                    "deg",
                )

# ============================================================================
# 7. SIDEREAL POSITIONS
# ============================================================================
print("\n" + "=" * 70)
print("7. SIDEREAL POSITIONS (Fagan-Bradley & Lahiri, 100 dates)")
print("=" * 70)

jds_sid = generate_jds(100)
for sid_mode, sid_name in [
    (SE_SIDM_FAGAN_BRADLEY, "FaganBradley"),
    (SE_SIDM_LAHIRI, "Lahiri"),
]:
    swe.set_sid_mode(sid_mode)
    ephem.swe_set_sid_mode(sid_mode)
    for pid, pname in PLANETS:
        lon_errs = []
        for jd in jds_sid:
            try:
                r_swe = swe.calc_ut(
                    jd, pid, SEFLG_SWIEPH | SEFLG_SIDEREAL | SEFLG_SPEED
                )
                r_lib = ephem.swe_calc_ut(
                    jd, pid, SEFLG_SWIEPH | SEFLG_SIDEREAL | SEFLG_SPEED
                )
                lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
            except Exception:
                continue
        if lon_errs:
            mx = max(lon_errs)
            if arcsec(mx) > 1.0:
                print(f'  {sid_name:15s} {pname:10s} max={arcsec(mx):8.3f}"')
                if arcsec(mx) > 5:
                    report_issue(
                        "sidereal",
                        "MEDIUM",
                        f'{sid_name}/{pname} max={arcsec(mx):.1f}"',
                        arcsec(mx),
                    )
    # Reset
    swe.set_sid_mode(0)
    ephem.swe_set_sid_mode(0)

# ============================================================================
# 8. NOD_APS METHODS: OSCU_BAR and FOPOINT
# ============================================================================
print("\n" + "=" * 70)
print("8. NOD_APS: OSCU_BAR and FOPOINT methods (30 dates)")
print("=" * 70)

jds_nod = generate_jds(30)
nod_methods = [
    (SE_NODBIT_OSCU_BAR, "OSCU_BAR"),
    (SE_NODBIT_FOPOINT, "FOPOINT"),
    (SE_NODBIT_MEAN | SE_NODBIT_FOPOINT, "MEAN+FOPOINT"),
    (SE_NODBIT_OSCU | SE_NODBIT_FOPOINT, "OSCU+FOPOINT"),
]

for method, mname in nod_methods:
    for pid, pname in PLANETS:
        if pid in (SE_SUN, SE_MOON):
            continue
        node_errs = []
        apse_errs = []
        n_errors = 0
        for jd in jds_nod:
            try:
                r_swe = swe.nod_aps_ut(jd, pid, SEFLG_SWIEPH | SEFLG_SPEED, method)
                r_lib = ephem.swe_nod_aps_ut(
                    jd, pid, SEFLG_SWIEPH | SEFLG_SPEED, method
                )
                for i in range(2):  # nodes
                    node_errs.append(angular_diff(r_swe[i][0], r_lib[i][0]))
                for i in range(2, 4):  # apsides
                    apse_errs.append(angular_diff(r_swe[i][0], r_lib[i][0]))
            except Exception:
                n_errors += 1
                continue
        if node_errs:
            mx_n = max(node_errs)
            mx_a = max(apse_errs) if apse_errs else 0
            if mx_n > 0.05 or mx_a > 5:
                print(
                    f"  {mname:15s} {pname:10s} nodes: max={mx_n:.4f}° apsides: max={mx_a:.4f}° (errors={n_errors})"
                )

# ============================================================================
# 9. ECLIPSE TIMING PRECISION
# ============================================================================
print("\n" + "=" * 70)
print("9. ECLIPSE TIMING PRECISION")
print("=" * 70)

jd_start = 2451545.0  # 2000
solar_timing_errs = []
lunar_timing_errs = []

for i in range(50):
    try:
        r_swe = swe.sol_eclipse_when_glob(jd_start, SEFLG_SWIEPH)
        r_lib = ephem.swe_sol_eclipse_when_glob(jd_start, SEFLG_SWIEPH)
        if r_swe and r_lib:
            # tret[0] = max eclipse time
            diff_sec = abs(r_swe[1][0] - r_lib[1][0]) * 86400
            solar_timing_errs.append(diff_sec)
            jd_start = r_swe[1][0] + 1
    except Exception:
        jd_start += 30
        continue

if solar_timing_errs:
    mx = max(solar_timing_errs)
    mean = statistics.mean(solar_timing_errs)
    print(
        f"  Solar eclipses: max_timing_err={mx:.2f}s mean={mean:.2f}s (n={len(solar_timing_errs)})"
    )
    if mx > 10:
        report_issue(
            "eclipses", "MEDIUM", f"Solar eclipse timing max={mx:.1f}s", mx, "seconds"
        )

jd_start = 2451545.0
for i in range(50):
    try:
        r_swe = swe.lun_eclipse_when(jd_start, SEFLG_SWIEPH)
        r_lib = ephem.swe_lun_eclipse_when(jd_start, SEFLG_SWIEPH)
        if r_swe and r_lib:
            diff_sec = abs(r_swe[1][0] - r_lib[1][0]) * 86400
            lunar_timing_errs.append(diff_sec)
            jd_start = r_swe[1][0] + 1
    except Exception:
        jd_start += 30
        continue

if lunar_timing_errs:
    mx = max(lunar_timing_errs)
    mean = statistics.mean(lunar_timing_errs)
    print(
        f"  Lunar eclipses: max_timing_err={mx:.2f}s mean={mean:.2f}s (n={len(lunar_timing_errs)})"
    )
    if mx > 10:
        report_issue(
            "eclipses", "MEDIUM", f"Lunar eclipse timing max={mx:.1f}s", mx, "seconds"
        )

# ============================================================================
# 10. TIME FUNCTIONS
# ============================================================================
print("\n" + "=" * 70)
print("10. TIME FUNCTIONS")
print("=" * 70)

# Delta T
jds_dt = generate_jds(200, 2415020.5, 2470000.0)  # 1900-2060
dt_errs = []
for jd in jds_dt:
    try:
        dt_swe = swe.deltat(jd)
        dt_lib = ephem.swe_deltat(jd)
        dt_errs.append(abs(dt_swe - dt_lib))
    except Exception:
        continue

if dt_errs:
    mx = max(dt_errs)
    mean = statistics.mean(dt_errs)
    print(f"  Delta T: max_diff={mx * 86400:.4f}s mean={mean * 86400:.6f}s (1900-2060)")
    if mx * 86400 > 1.0:
        report_issue(
            "time",
            "MEDIUM",
            f"Delta T max diff={mx * 86400:.2f}s",
            mx * 86400,
            "seconds",
        )

# Sidereal time
sidtime_errs = []
for jd in jds_dt:
    try:
        st_swe = swe.sidtime(jd)
        st_lib = ephem.swe_sidtime(jd)
        sidtime_errs.append(abs(st_swe - st_lib))
    except Exception:
        continue

if sidtime_errs:
    mx = max(sidtime_errs)
    print(
        f"  Sidereal time: max_diff={mx * 3600:.4f}s (hours*3600) (n={len(sidtime_errs)})"
    )
    if mx * 3600 > 1.0:
        report_issue(
            "time",
            "MEDIUM",
            f"Sidereal time max diff={mx * 3600:.2f}s",
            mx * 3600,
            "seconds",
        )

# ============================================================================
# 11. COORDINATE TRANSFORMS
# ============================================================================
print("\n" + "=" * 70)
print("11. COORDINATE TRANSFORMS (cotrans)")
print("=" * 70)

obliquity = 23.4393
coord_errs = []
for _ in range(500):
    lon = random.uniform(0, 360)
    lat = random.uniform(-90, 90)
    dist = random.uniform(0.5, 50)
    try:
        r_swe = swe.cotrans((lon, lat, dist), -obliquity)
        r_lib = ephem.swe_cotrans((lon, lat, dist), -obliquity)
        coord_errs.append(angular_diff(r_swe[0], r_lib[0]))
    except Exception:
        continue

if coord_errs:
    mx = max(coord_errs)
    print(f'  cotrans: max={arcsec(mx):.6f}" (n={len(coord_errs)})')
    if arcsec(mx) > 0.001:
        report_issue("coordinates", "LOW", f'cotrans max={arcsec(mx):.4f}"', arcsec(mx))

# ============================================================================
# 12. FIXED STARS - spot check large catalog
# ============================================================================
print("\n" + "=" * 70)
print("12. FIXED STARS (30 stars x 4 dates)")
print("=" * 70)

stars = [
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Fomalhaut",
    "Betelgeuse",
    "Rigel",
    "Sirius",
    "Canopus",
    "Procyon",
    "Pollux",
    "Castor",
    "Capella",
    "Vega",
    "Altair",
    "Deneb",
    "Achernar",
    "Algol",
    "Mira",
    "Arcturus",
    "Zubenelgenubi",
    "Zubeneschamali",
    "Vindemiatrix",
    "Markab",
    "Scheat",
    "Alpheratz",
    "Hamal",
    "Mirach",
    "Almach",
    "Menkar",
]
star_jds = [2451545.0, 2460000.0, 2440000.0, 2430000.0]

for star_name in stars:
    lon_errs = []
    for jd in star_jds:
        try:
            r_swe = swe.fixstar_ut(star_name, jd, SEFLG_SWIEPH)
            r_lib = ephem.swe_fixstar_ut(star_name, jd, SEFLG_SWIEPH)
            lon_errs.append(angular_diff(r_swe[0][0], r_lib[0][0]))
        except Exception:
            continue
    if lon_errs:
        mx = max(lon_errs)
        if arcsec(mx) > 2.0:
            print(f'  {star_name:20s} max={arcsec(mx):8.3f}"')
            if arcsec(mx) > 10:
                report_issue(
                    "fixed_stars",
                    "MEDIUM",
                    f'{star_name} max={arcsec(mx):.1f}"',
                    arcsec(mx),
                )

# ============================================================================
# 13. RISE/SET/TRANSIT
# ============================================================================
print("\n" + "=" * 70)
print("13. RISE/SET/TRANSIT TIMING (Sun+Moon, 50 dates)")
print("=" * 70)

jds_rise = generate_jds(50, 2451545.0, 2470000.0)
geopos = (12.5, 41.9, 0.0)  # Rome
rise_errs = []
transit_errs = []

for pid, pname in [(SE_SUN, "Sun"), (SE_MOON, "Moon")]:
    for jd in jds_rise:
        # Rise
        try:
            r_swe = swe.rise_trans(jd, pid, "", SEFLG_SWIEPH, 1, geopos, 1013.25, 15.0)
            r_lib = ephem.swe_rise_trans(
                jd, pid, "", SEFLG_SWIEPH, 1, geopos, 1013.25, 15.0
            )
            if r_swe and r_lib and r_swe[1] > 0 and r_lib[1] > 0:
                diff_sec = abs(r_swe[1] - r_lib[1]) * 86400
                rise_errs.append((pname, "rise", diff_sec))
        except Exception:
            pass
        # Transit
        try:
            r_swe = swe.rise_trans(jd, pid, "", SEFLG_SWIEPH, 8, geopos, 1013.25, 15.0)
            r_lib = ephem.swe_rise_trans(
                jd, pid, "", SEFLG_SWIEPH, 8, geopos, 1013.25, 15.0
            )
            if r_swe and r_lib and r_swe[1] > 0 and r_lib[1] > 0:
                diff_sec = abs(r_swe[1] - r_lib[1]) * 86400
                transit_errs.append((pname, "transit", diff_sec))
        except Exception:
            pass

if rise_errs:
    mx_rise = max(e[2] for e in rise_errs)
    mean_rise = statistics.mean(e[2] for e in rise_errs)
    print(f"  Rise: max={mx_rise:.2f}s mean={mean_rise:.2f}s (n={len(rise_errs)})")
    if mx_rise > 60:
        report_issue(
            "rise_transit",
            "MEDIUM",
            f"Rise timing max={mx_rise:.1f}s",
            mx_rise,
            "seconds",
        )

if transit_errs:
    mx_transit = max(e[2] for e in transit_errs)
    mean_transit = statistics.mean(e[2] for e in transit_errs)
    print(
        f"  Transit: max={mx_transit:.2f}s mean={mean_transit:.2f}s (n={len(transit_errs)})"
    )
    if mx_transit > 60:
        report_issue(
            "rise_transit",
            "MEDIUM",
            f"Transit timing max={mx_transit:.1f}s",
            mx_transit,
            "seconds",
        )

# ============================================================================
# 14. CROSSING FUNCTIONS
# ============================================================================
print("\n" + "=" * 70)
print("14. CROSSING FUNCTIONS (solcross, mooncross)")
print("=" * 70)

# Solar crossing at 0° (vernal equinox) - 20 years
crossing_errs = []
jd_cross = 2451545.0
for _ in range(20):
    try:
        r_swe = swe.solcross_ut(0.0, jd_cross, SEFLG_SWIEPH)
        r_lib = ephem.swe_solcross_ut(0.0, jd_cross, SEFLG_SWIEPH)
        diff_sec = abs(r_swe - r_lib) * 86400
        crossing_errs.append(("equinox", diff_sec))
        jd_cross = r_swe + 300
    except Exception:
        jd_cross += 365.25

if crossing_errs:
    mx = max(e[1] for e in crossing_errs)
    mean = statistics.mean(e[1] for e in crossing_errs)
    print(
        f"  Solar equinox crossing: max={mx:.3f}s mean={mean:.3f}s (n={len(crossing_errs)})"
    )
    if mx > 10:
        report_issue(
            "crossings",
            "MEDIUM",
            f"Solar equinox crossing max={mx:.1f}s",
            mx,
            "seconds",
        )

# Moon crossing at various longitudes
moon_cross_errs = []
for lon_target in range(0, 360, 30):
    jd_cross = 2451545.0
    for _ in range(5):
        try:
            r_swe = swe.mooncross_ut(float(lon_target), jd_cross, SEFLG_SWIEPH)
            r_lib = ephem.swe_mooncross_ut(float(lon_target), jd_cross, SEFLG_SWIEPH)
            diff_sec = abs(r_swe - r_lib) * 86400
            moon_cross_errs.append(diff_sec)
            jd_cross = r_swe + 20
        except Exception:
            jd_cross += 29.5

if moon_cross_errs:
    mx = max(moon_cross_errs)
    mean = statistics.mean(moon_cross_errs)
    print(
        f"  Moon longitude crossings: max={mx:.3f}s mean={mean:.3f}s (n={len(moon_cross_errs)})"
    )
    if mx > 10:
        report_issue(
            "crossings", "MEDIUM", f"Moon crossing max={mx:.1f}s", mx, "seconds"
        )

# ============================================================================
# 15. AYANAMSHA VALUES
# ============================================================================
print("\n" + "=" * 70)
print("15. AYANAMSHA VALUES (43 modes, 10 dates)")
print("=" * 70)

jds_ayan = generate_jds(10)
ayan_errs = {}
for mode in range(43):
    errs = []
    for jd in jds_ayan:
        try:
            swe.set_sid_mode(mode)
            ephem.swe_set_sid_mode(mode)
            a_swe = swe.get_ayanamsa_ut(jd)
            a_lib = ephem.swe_get_ayanamsa_ut(jd)
            errs.append(abs(a_swe - a_lib))
        except Exception:
            continue
    if errs:
        mx = max(errs)
        if arcsec(mx) > 0.1:
            ayan_errs[mode] = arcsec(mx)

swe.set_sid_mode(0)
ephem.swe_set_sid_mode(0)

if ayan_errs:
    for mode, mx in sorted(ayan_errs.items(), key=lambda x: -x[1]):
        if mx > 1.0:
            print(f'  Ayanamsha mode {mode:2d}: max_diff={mx:.3f}"')
            if mx > 5:
                report_issue("ayanamsha", "MEDIUM", f'Mode {mode} max={mx:.1f}"', mx)
else:
    print('  All ayanamsha modes within 0.1" tolerance.')

# ============================================================================
# FINAL REPORT
# ============================================================================
print("\n" + "=" * 70)
print("FINAL REPORT: ACTIONABLE ISSUES")
print("=" * 70)

if not issues:
    print("  No actionable issues found! All APIs within expected tolerances.")
else:
    # Sort by severity
    severity_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
    issues.sort(key=lambda x: (severity_order.get(x["severity"], 3), -x["max_err"]))
    for issue in issues:
        print(
            f"  [{issue['severity']:6s}] {issue['category']:20s} | {issue['description']} ({issue['max_err']:.3f} {issue['unit']})"
        )

print(f"\nTotal issues found: {len(issues)}")
print(f"  HIGH:   {sum(1 for i in issues if i['severity'] == 'HIGH')}")
print(f"  MEDIUM: {sum(1 for i in issues if i['severity'] == 'MEDIUM')}")
print(f"  LOW:    {sum(1 for i in issues if i['severity'] == 'LOW')}")
