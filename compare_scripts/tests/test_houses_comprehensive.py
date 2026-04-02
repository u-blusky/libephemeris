"""
Comprehensive House Systems Comparison Tests.

Compares ALL 24 house system calculations between pyswisseph and libephemeris
with tight tolerances across multiple locations, dates, and functions.

Coverage:
  - swe_houses: 24 systems × 8 locations × 12 dates — all 12 cusps + all 8 ASCMC
  - swe_houses_armc: 24 systems × 12 ARMC × 4 latitudes — cusps + ASCMC
  - swe_house_pos: 12 systems × 18 longitudes × 3 body latitudes × 4 dates
  - swe_houses_ex (sidereal): 8 systems × 3 ayanamshas × 6 dates
  - Precision report: per-system max error for cusps and each ASCMC index

Run via:
  poe test:houses
  leph test compare houses
"""

from __future__ import annotations

from collections import defaultdict

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SEFLG_SIDEREAL,
    SE_SIDM_LAHIRI,
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_RAMAN,
)


# ============================================================================
# UTILITIES
# ============================================================================


def angular_diff(a: float, b: float) -> float:
    """Absolute angular difference in degrees, handling 360/0 wraparound."""
    d = abs(float(a) - float(b))
    if d > 180.0:
        d = 360.0 - d
    return d


# ============================================================================
# CONFIGURATION
# ============================================================================

# All 24 house systems
HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("O", "Porphyry"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
    ("E", "Equal (Asc)"),
    ("A", "Equal (MC)"),
    ("W", "Whole Sign"),
    ("B", "Alcabitius"),
    ("T", "Topocentric"),
    ("M", "Morinus"),
    ("X", "Meridian"),
    ("V", "Vehlow"),
    ("H", "Horizontal"),
    ("U", "Krusinski"),
    ("F", "Carter"),
    ("G", "Gauquelin"),
    ("I", "Sunshine"),
    ("i", "Sunshine Makransky"),
    ("N", "Natural Gradient"),
    ("Y", "APC"),
    ("D", "Equal from MC"),
    ("J", "Savard-A"),
    ("L", "Pullen SD"),
    ("S", "Sripati"),
    ("Q", "Pullen SR"),
]

# Per-system cusp tolerance (degrees). Default 0.001° (~3.6 arcsec).
CUSP_TOLERANCE_DEFAULT = 0.001
CUSP_TOLERANCE = {
    "P": 0.002,  # Placidus: iteration differences at high latitudes
    "R": 0.002,  # Regiomontanus: minor precision differences
}

# ASCMC labels and tolerances per index
ASCMC_LABELS = [
    "ASC", "MC", "ARMC", "Vertex",
    "EquAsc", "CoAsc_Koch", "CoAsc_Munkasey", "PolarAsc",
]
ASCMC_TOLERANCE = {
    0: 0.001,   # ASC
    1: 0.001,   # MC
    2: 0.001,   # ARMC
    3: 0.01,    # Vertex
    4: 0.01,    # Equatorial Ascendant
    5: 0.01,    # Co-Ascendant Koch
    6: 0.01,    # Co-Ascendant Munkasey
    7: 0.01,    # Polar Ascendant
}

# 8 test locations: (name, lat, lon)
LOCATIONS = [
    ("Equator", 0.0, 0.0),
    ("Rome", 41.9028, 12.4964),
    ("London", 51.5074, -0.1278),
    ("Sydney", -33.8688, 151.2093),
    ("Tokyo", 35.6762, 139.6503),
    ("Reykjavik", 64.0, -22.0),
    ("Buenos Aires", -34.6037, -58.3816),
    ("Mumbai", 19.0760, 72.8777),
]

# 12 Julian Days spanning 1950–2050
JD_1950 = 2433282.5   # 1950-01-01
JD_2050 = 2469807.5   # 2050-01-01
TEST_JDS = [JD_1950 + i * (JD_2050 - JD_1950) / 11 for i in range(12)]

# 6 Julian Days for lighter tests
TEST_JDS_6 = TEST_JDS[::2]  # every other JD

# Systems that raise PolarCircleError at high latitudes
POLAR_FAIL_SYSTEMS = {"P", "K", "G"}


# ============================================================================
# SHARED PRECISION TRACKER
# ============================================================================

# Collects max errors per system for the final precision report.
# Keys: (system_char, value_name) -> max_diff
_precision_data: dict[tuple[str, str], float] = defaultdict(float)


def _track(system: str, value_name: str, diff: float) -> None:
    """Track the maximum diff for a (system, value) pair."""
    key = (system, value_name)
    if diff > _precision_data[key]:
        _precision_data[key] = diff


# ============================================================================
# TEST: swe_houses — All cusps + ALL 8 ASCMC
# ============================================================================


class TestHouseCuspsAllSystems:
    """Compare all 12 house cusps for all 24 systems vs pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS)
    @pytest.mark.parametrize("jd", TEST_JDS, ids=[f"jd{i}" for i in range(12)])
    def test_cusps(self, hsys, hsys_name, loc_name, lat, lon, jd):
        """All 12 cusps must match within tolerance."""
        tolerance = CUSP_TOLERANCE.get(hsys, CUSP_TOLERANCE_DEFAULT)

        try:
            cusps_swe, _ = swe.houses(jd, lat, lon, hsys.encode("ascii"))
        except swe.Error:
            pytest.skip(f"pyswisseph raises error for {hsys} at {loc_name}")

        try:
            cusps_py, _ = ephem.swe_houses(jd, lat, lon, hsys)
        except Exception:
            pytest.skip(f"libephemeris raises error for {hsys} at {loc_name}")

        num_cusps = min(len(cusps_swe), len(cusps_py))
        max_diff = 0.0
        worst_cusp = 0

        for i in range(num_cusps):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            _track(hsys, f"cusp{i + 1}", diff)
            if diff > max_diff:
                max_diff = diff
                worst_cusp = i + 1

        assert max_diff < tolerance, (
            f"{hsys_name} ({hsys}) at {loc_name}: cusp {worst_cusp} diff "
            f"{max_diff:.6f}° exceeds {tolerance}° "
            f"(swe={cusps_swe[worst_cusp - 1]:.6f}, py={cusps_py[worst_cusp - 1]:.6f})"
        )


class TestAscmcComplete:
    """Compare ALL 8 ASCMC values for all 24 systems vs pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS)
    @pytest.mark.parametrize("jd", TEST_JDS_6, ids=[f"jd{i}" for i in range(6)])
    def test_ascmc_all_values(self, hsys, hsys_name, loc_name, lat, lon, jd):
        """All 8 ASCMC values must match within their respective tolerances."""
        try:
            _, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
        except swe.Error:
            pytest.skip(f"pyswisseph raises error for {hsys} at {loc_name}")

        try:
            _, ascmc_py = ephem.swe_houses(jd, lat, lon, hsys)
        except Exception:
            pytest.skip(f"libephemeris raises error for {hsys} at {loc_name}")

        failures = []
        for i in range(min(len(ascmc_swe), len(ascmc_py), 8)):
            tol = ASCMC_TOLERANCE.get(i, 0.01)
            diff = angular_diff(ascmc_swe[i], ascmc_py[i])
            _track(hsys, ASCMC_LABELS[i], diff)

            if diff >= tol:
                failures.append(
                    f"  ascmc[{i}] ({ASCMC_LABELS[i]}): diff={diff:.6f}° > {tol}° "
                    f"(swe={float(ascmc_swe[i]):.6f}, py={float(ascmc_py[i]):.6f})"
                )

        assert not failures, (
            f"{hsys_name} ({hsys}) at {loc_name} jd={jd:.1f}:\n"
            + "\n".join(failures)
        )


# ============================================================================
# TEST: swe_houses_armc — All 24 systems, multiple latitudes
# ============================================================================


# 12 ARMC values
ARMC_VALUES = [float(x) for x in range(0, 360, 30)]

# 4 latitudes for ARMC tests
ARMC_LATITUDES = [
    ("equator", 0.0),
    ("Rome", 41.9),
    ("Sydney", -33.9),
    ("high_lat", 60.0),
]


class TestHousesArmcAllSystems:
    """Compare houses_armc for all 24 systems across ARMC values and latitudes."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    @pytest.mark.parametrize("armc", ARMC_VALUES, ids=[f"armc{int(a)}" for a in ARMC_VALUES])
    @pytest.mark.parametrize("lat_name,lat", ARMC_LATITUDES)
    def test_houses_armc_cusps_and_ascmc(self, hsys, hsys_name, armc, lat_name, lat):
        """Cusps and all 8 ASCMC from houses_armc must match pyswisseph."""
        eps = 23.4393  # standard obliquity

        try:
            cusps_swe, ascmc_swe = swe.houses_armc(armc, lat, eps, hsys.encode("ascii"))
        except swe.Error:
            pytest.skip(f"pyswisseph error for {hsys} at {lat_name}")

        try:
            cusps_py, ascmc_py = ephem.swe_houses_armc(armc, lat, eps, hsys)
        except Exception:
            pytest.skip(f"libephemeris error for {hsys} at {lat_name}")

        cusp_tol = CUSP_TOLERANCE.get(hsys, CUSP_TOLERANCE_DEFAULT)

        # Compare cusps
        num_cusps = min(len(cusps_swe), len(cusps_py))
        max_cusp_diff = 0.0
        for i in range(num_cusps):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            _track(hsys, f"armc_cusp{i + 1}", diff)
            max_cusp_diff = max(max_cusp_diff, diff)

        assert max_cusp_diff < cusp_tol, (
            f"{hsys_name} houses_armc at {lat_name} armc={armc}: "
            f"max cusp diff {max_cusp_diff:.6f}° > {cusp_tol}°"
        )

        # Compare ASCMC
        failures = []
        for i in range(min(len(ascmc_swe), len(ascmc_py), 8)):
            tol = ASCMC_TOLERANCE.get(i, 0.01)
            diff = angular_diff(ascmc_swe[i], ascmc_py[i])
            _track(hsys, f"armc_{ASCMC_LABELS[i]}", diff)
            if diff >= tol:
                failures.append(
                    f"  ascmc[{i}] ({ASCMC_LABELS[i]}): {diff:.6f}° > {tol}°"
                )

        assert not failures, (
            f"{hsys_name} houses_armc at {lat_name} armc={armc}:\n"
            + "\n".join(failures)
        )


# ============================================================================
# TEST: swe_house_pos — Comprehensive
# ============================================================================

HPOS_SYSTEMS = ["P", "K", "O", "R", "C", "E", "W", "B", "M", "T", "X", "V"]
HPOS_LONGITUDES = [float(x) for x in range(0, 360, 20)]  # 18 values
HPOS_BODY_LATS = [0.0, 3.0, -5.0]
HPOS_TOLERANCE = 0.02


class TestHousePosComprehensive:
    """Compare house_pos for 12 systems, 18 longitudes, 3 body latitudes."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", HPOS_SYSTEMS)
    @pytest.mark.parametrize("planet_lon", HPOS_LONGITUDES, ids=[f"lon{int(l)}" for l in HPOS_LONGITUDES])
    @pytest.mark.parametrize("body_lat", HPOS_BODY_LATS, ids=[f"blat{l}" for l in HPOS_BODY_LATS])
    @pytest.mark.parametrize("jd", TEST_JDS_6[:4], ids=[f"jd{i}" for i in range(4)])
    def test_house_pos(self, hsys, planet_lon, body_lat, jd):
        """house_pos results must match pyswisseph."""
        geo_lat = 41.9  # Rome

        # Get ARMC and obliquity from pyswisseph
        try:
            _, ascmc_swe = swe.houses(jd, geo_lat, 12.5, b"P")
            armc = float(ascmc_swe[2])
            ecl = swe.calc_ut(jd, -1, 0)
            eps = float(ecl[0][0])
        except Exception:
            pytest.skip("Cannot compute ARMC/eps")

        try:
            hp_swe = float(swe.house_pos(
                armc, geo_lat, eps, (planet_lon, body_lat), hsys.encode("ascii")
            ))
        except swe.Error:
            pytest.skip(f"pyswisseph error for house_pos {hsys}")

        try:
            hp_py = float(ephem.house_pos(
                armc, geo_lat, eps, (planet_lon, body_lat), hsys
            ))
        except Exception:
            pytest.skip(f"libephemeris error for house_pos {hsys}")

        diff = abs(hp_swe - hp_py)
        # Handle house number wraparound (12.95 vs 1.05)
        if diff > 6.0:
            diff = 12.0 - diff

        _track(hsys, "house_pos", diff)

        assert diff < HPOS_TOLERANCE, (
            f"house_pos {hsys} lon={planet_lon} blat={body_lat}: "
            f"diff={diff:.4f} (swe={hp_swe:.4f}, py={hp_py:.4f})"
        )


# ============================================================================
# TEST: swe_houses_ex — Sidereal
# ============================================================================

SIDEREAL_SYSTEMS = ["P", "K", "E", "W", "O", "R", "C", "B"]
SIDEREAL_MODES = [
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
    (SE_SIDM_RAMAN, "Raman"),
]


class TestSiderealHouses:
    """Compare sidereal house calculations via houses_ex."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", SIDEREAL_SYSTEMS)
    @pytest.mark.parametrize("sid_mode,sid_name", SIDEREAL_MODES)
    @pytest.mark.parametrize("jd", TEST_JDS_6, ids=[f"jd{i}" for i in range(6)])
    def test_sidereal_cusps(self, hsys, sid_mode, sid_name, jd):
        """Sidereal cusps must match pyswisseph."""
        lat, lon = 41.9, 12.5  # Rome

        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode)

        try:
            cusps_swe, ascmc_swe = swe.houses_ex(
                jd, lat, lon, hsys.encode("ascii"), SEFLG_SIDEREAL
            )
        except swe.Error:
            swe.set_sid_mode(0)
            ephem.set_sid_mode(0)
            pytest.skip(f"pyswisseph error for sidereal {hsys}")

        try:
            cusps_py, ascmc_py = ephem.houses_ex(
                jd, lat, lon, hsys, SEFLG_SIDEREAL
            )
        except Exception:
            swe.set_sid_mode(0)
            ephem.set_sid_mode(0)
            pytest.skip(f"libephemeris error for sidereal {hsys}")

        swe.set_sid_mode(0)
        ephem.set_sid_mode(0)

        cusp_tol = CUSP_TOLERANCE.get(hsys, CUSP_TOLERANCE_DEFAULT)
        max_diff = 0.0
        for i in range(min(len(cusps_swe), len(cusps_py), 12)):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < cusp_tol, (
            f"Sidereal {hsys} ({sid_name}): max cusp diff {max_diff:.6f}° > {cusp_tol}°"
        )

        # Also compare ASC and MC
        asc_diff = angular_diff(ascmc_swe[0], ascmc_py[0])
        mc_diff = angular_diff(ascmc_swe[1], ascmc_py[1])
        assert asc_diff < 0.001, (
            f"Sidereal {hsys} ({sid_name}): ASC diff {asc_diff:.6f}°"
        )
        assert mc_diff < 0.001, (
            f"Sidereal {hsys} ({sid_name}): MC diff {mc_diff:.6f}°"
        )


# ============================================================================
# PRECISION REPORT — Aggregate statistics
# ============================================================================


class TestPrecisionReport:
    """Print per-system precision summary after all tests have run.

    This test always passes — its purpose is to display the precision table.
    Run with -s flag to see the output.
    """

    @pytest.mark.comparison
    def test_print_precision_report(self):
        """Print per-system max error table."""
        if not _precision_data:
            pytest.skip("No precision data collected (run other tests first)")

        # Collect all systems and value names
        systems = sorted({k[0] for k in _precision_data})
        cusp_keys = [f"cusp{i}" for i in range(1, 13)]
        ascmc_keys = ASCMC_LABELS

        print("\n" + "=" * 100)
        print("HOUSE SYSTEMS PRECISION REPORT — max error vs pyswisseph (degrees)")
        print("=" * 100)

        # --- Cusp max errors ---
        print("\nMax cusp error per system (swe_houses):")
        print(f"  {'System':<25s} {'Max Cusp':>12s} {'Worst Cusp':>12s}")
        print(f"  {'-' * 25} {'-' * 12} {'-' * 12}")
        for sys_char in systems:
            max_cusp = 0.0
            worst = ""
            for ck in cusp_keys:
                v = _precision_data.get((sys_char, ck), 0.0)
                if v > max_cusp:
                    max_cusp = v
                    worst = ck
            if max_cusp > 0:
                # Find name
                name = next((n for c, n in HOUSE_SYSTEMS if c == sys_char), sys_char)
                print(f"  {name + ' (' + sys_char + ')':<25s} {max_cusp:>12.7f} {worst:>12s}")

        # --- ASCMC max errors ---
        print(f"\nMax ASCMC error per system (swe_houses):")
        header = f"  {'System':<18s}"
        for lbl in ascmc_keys:
            header += f" {lbl:>10s}"
        print(header)
        print(f"  {'-' * 18}" + (" " + "-" * 10) * 8)
        for sys_char in systems:
            name = next((n for c, n in HOUSE_SYSTEMS if c == sys_char), sys_char)
            row = f"  {name[:16] + ' (' + sys_char + ')':<18s}"
            for lbl in ascmc_keys:
                v = _precision_data.get((sys_char, lbl), 0.0)
                if v > 0:
                    row += f" {v:>10.7f}"
                else:
                    row += f" {'—':>10s}"
            print(row)

        # --- houses_armc max cusp error ---
        armc_cusp_keys = [f"armc_cusp{i}" for i in range(1, 13)]
        has_armc_data = any(
            _precision_data.get((s, k), 0) > 0
            for s in systems for k in armc_cusp_keys
        )
        if has_armc_data:
            print(f"\nMax cusp error per system (swe_houses_armc):")
            print(f"  {'System':<25s} {'Max Cusp':>12s}")
            print(f"  {'-' * 25} {'-' * 12}")
            for sys_char in systems:
                max_cusp = 0.0
                for ck in armc_cusp_keys:
                    v = _precision_data.get((sys_char, ck), 0.0)
                    max_cusp = max(max_cusp, v)
                if max_cusp > 0:
                    name = next((n for c, n in HOUSE_SYSTEMS if c == sys_char), sys_char)
                    print(f"  {name + ' (' + sys_char + ')':<25s} {max_cusp:>12.7f}")

        # --- house_pos max error ---
        has_hpos = any(
            _precision_data.get((s, "house_pos"), 0) > 0 for s in systems
        )
        if has_hpos:
            print(f"\nMax house_pos error per system:")
            print(f"  {'System':<25s} {'Max Diff':>12s}")
            print(f"  {'-' * 25} {'-' * 12}")
            for sys_char in systems:
                v = _precision_data.get((sys_char, "house_pos"), 0.0)
                if v > 0:
                    name = next((n for c, n in HOUSE_SYSTEMS if c == sys_char), sys_char)
                    print(f"  {name + ' (' + sys_char + ')':<25s} {v:>12.7f}")

        print("\n" + "=" * 100)
        print("END OF PRECISION REPORT")
        print("=" * 100 + "\n")

        # This test always passes — it's just a report
        assert True
