"""
Moshier Historical Deep-Date Comparison Tests (C vs Python).

Validates Moshier (SEFLG_MOSEPH) planetary calculations for deep historical
dates outside the DE440 range (1550-2650 CE), comparing pyswisseph (C library)
against libephemeris (Python reimplementation).

Historical dates are the PRIMARY use case for Moshier mode: archaeological
astronomy, astro-chronology, dating ancient eclipses, and biblical-era
planetary conjunctions all require reliable positions at dates far from J2000.

The VSOP87 (planets) and ELP2000-82B (Moon) series expansions degrade with
distance from J2000.0, and the C and Python implementations may diverge
differently due to:
- Different truncation of periodic series terms
- Different precession/nutation models
- Numerical precision differences in polynomial evaluation

This test produces a systematic [planet x epoch] precision matrix, documenting
which combinations are reliable for historical applications.

Tolerance bands (C-vs-Python longitude agreement):
    Planets (VSOP87: Sun, Mercury-Neptune):
        1000-1550 CE:  0.02 deg (72 arcsec)   - recent pre-DE440
        0-1000 CE:     0.10 deg (6 arcmin)     - classical/medieval
        -1000-0 CE:    0.50 deg (30 arcmin)    - ancient
        < -1000 CE:    1.00 deg (1 degree)     - deep ancient

    Moon (ELP2000-82B, diverges more than VSOP87):
        1000-1550 CE:  0.15 deg (9 arcmin)
        0-1000 CE:     0.50 deg (30 arcmin)
        -1000-0 CE:    1.00 deg (1 degree)
        < -1000 CE:    2.00 deg (2 degrees)

    Pluto (Chapront-Francou perturbation theory, very wide):
        1000-1550 CE:  2.0 deg
        0-1000 CE:     5.0 deg
        -1000-0 CE:   10.0 deg
        < -1000 CE:   20.0 deg

Total: 80 parametrized test cases (8 historical dates x 10 planets).
All individual tests marked @pytest.mark.slow.
"""

from __future__ import annotations

import pytest
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
    SEFLG_MOSEPH,
    SEFLG_SPEED,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# 8 deep historical dates outside DE440 range (1550-2650 CE).
# These are the primary use case for Moshier mode: archaeological astronomy,
# astro-chronology, dating ancient eclipses, biblical-era conjunctions.
# Year uses astronomical year numbering: year 0 = 1 BCE, -1 = 2 BCE, etc.
HISTORICAL_DATES = [
    (-3000, 6, 21, 12.0, "-3000 CE (pre-dynastic Egypt)"),
    (-2000, 1, 1, 12.0, "-2000 CE (early Mesopotamia)"),
    (-1000, 6, 15, 12.0, "-1000 CE (early Iron Age)"),
    (-500, 3, 21, 12.0, "-500 CE (classical Greece)"),
    (0, 1, 1, 12.0, "0 CE / 1 BCE"),
    (500, 3, 21, 12.0, "500 CE (early Middle Ages)"),
    (1000, 12, 25, 12.0, "1000 CE (high Middle Ages)"),
    (1200, 7, 4, 12.0, "1200 CE (late Middle Ages)"),
]

# All 10 major planets for comprehensive coverage.
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

# Era labels for the precision report (ordered by increasing distance
# from J2000, matching the tolerance graduation).
ERA_LABELS = {
    "deep_ancient": "< -1000 CE",
    "ancient": "-1000 to 0 CE",
    "medieval_early": "0 to 1000 CE",
    "medieval_late": "1000 to 1550 CE",
}


# ============================================================================
# TOLERANCE FUNCTIONS
# ============================================================================


def _tolerance_for_year(year: int, body_id: int) -> float:
    """Return graduated longitude tolerance in degrees based on era and body.

    VSOP87 (planets) and ELP2000-82B (Moon) are series expansions about
    J2000.0 that degrade with distance from that epoch. The C library
    (pyswisseph) and Python (libephemeris) use different truncations of
    these series, so differences between them grow with distance from J2000.

    Pluto uses Chapront-Francou perturbation theory with 2993 tabulated
    terms; C-vs-Python porting differences can reach several degrees at
    dates far from J2000.

    Args:
        year: Astronomical year (0 = 1 BCE, -1 = 2 BCE, etc.).
        body_id: Swiss Ephemeris body ID.

    Returns:
        Tolerance in degrees.
    """
    if body_id == SE_PLUTO:
        # Pluto Chapront-Francou theory has massive C-vs-Python divergence
        # at historical dates due to 2993-row perturbation table porting diffs.
        if year >= 1000:
            return 2.0
        elif year >= 0:
            return 5.0
        elif year >= -1000:
            return 10.0
        else:
            return 20.0

    if body_id == SE_MOON:
        # ELP2000-82B diverges more than VSOP87 between C and Python.
        # Empirically, Moon differences reach ~0.12 deg at 3000 CE
        # and grow further at historical dates.
        if year >= 1000:
            return 0.15
        elif year >= 0:
            return 0.5
        elif year >= -1000:
            return 1.0
        else:
            return 2.0

    # Planets (VSOP87 series): Sun, Mercury, Venus, Mars, Jupiter-Neptune
    if year >= 1000:
        return 0.02
    elif year >= 0:
        return 0.1
    elif year >= -1000:
        return 0.5
    else:
        return 1.0


def _distance_tolerance(year: int, body_id: int) -> float:
    """Return relative distance tolerance for a given body and era.

    Distance tolerances are graduated similarly to longitude, because
    orbital radius computation shares the same VSOP87/ELP series.

    Args:
        year: Astronomical year.
        body_id: Swiss Ephemeris body ID.

    Returns:
        Relative tolerance (dimensionless).
    """
    if body_id == SE_PLUTO:
        if year >= 0:
            return 0.05
        else:
            return 0.1
    if body_id == SE_MOON:
        if year >= 0:
            return 0.005
        else:
            return 0.01
    # Planets
    if year >= 0:
        return 0.002
    else:
        return 0.005


def _era_label(year: int) -> str:
    """Return era key for the precision report matrix."""
    if year >= 1000:
        return "medieval_late"
    elif year >= 0:
        return "medieval_early"
    elif year >= -1000:
        return "ancient"
    else:
        return "deep_ancient"


# ============================================================================
# HELPERS
# ============================================================================


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestHistoricalPositions:
    """Compare Moshier planetary positions at deep historical dates.

    80 parametrized test cases (8 dates x 10 planets) comparing C (pyswisseph)
    vs Python (libephemeris) with SEFLG_MOSEPH | SEFLG_SPEED.

    Each test validates longitude, latitude, and distance with graduated
    tolerances based on distance from J2000 and body type. Prints a
    degradation report visible with ``pytest -s`` for precision analysis.
    """

    @pytest.mark.slow
    @pytest.mark.comparison
    @pytest.mark.edge_case
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc",
        HISTORICAL_DATES,
        ids=[d[4] for d in HISTORICAL_DATES],
    )
    @pytest.mark.parametrize(
        "body_id,body_name",
        PLANETS,
        ids=[p[1] for p in PLANETS],
    )
    def test_historical_position(
        self, year, month, day, hour, date_desc, body_id, body_name
    ):
        """MOSEPH position at historical date should match between C and Python.

        Compares longitude, latitude, and distance for each planet/date
        combination using SEFLG_MOSEPH | SEFLG_SPEED with graduated
        tolerances based on distance from J2000 and body type.

        Prints degradation report for analysis with ``pytest -s``.
        """
        jd = swe.julday(year, month, day, hour)
        flag = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, body_id, flag)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        tol_lon = _tolerance_for_year(year, body_id)
        dist_from_j2000 = abs(year - 2000)

        # Distance comparison (relative)
        diff_dist_rel = 0.0
        dist_tol = _distance_tolerance(year, body_id)
        if pos_swe[2] > 0:
            diff_dist_rel = abs(pos_swe[2] - pos_py[2]) / pos_swe[2]

        # Degradation report visible with pytest -s
        print(
            f"\n  [HISTORICAL] {body_name} @ {date_desc} "
            f"(JD {jd:.1f}, {dist_from_j2000}y from J2000):\n"
            f"    lon: swe={pos_swe[0]:.6f} py={pos_py[0]:.6f} "
            f"diff={diff_lon:.6f} deg  tol={tol_lon:.3f} deg\n"
            f"    lat: swe={pos_swe[1]:.6f} py={pos_py[1]:.6f} "
            f"diff={diff_lat:.6f} deg  tol={tol_lon:.3f} deg\n"
            f"    dist: swe={pos_swe[2]:.8f} py={pos_py[2]:.8f} "
            f"rel_diff={diff_dist_rel:.6f}  tol={dist_tol:.4f}"
        )

        # Longitude assertion
        assert diff_lon < tol_lon, (
            f"{body_name} MOSEPH at {date_desc}: "
            f"longitude diff {diff_lon:.6f} deg exceeds "
            f"{tol_lon} deg tolerance\n"
            f"  pyswisseph:   {pos_swe[0]:.6f} deg\n"
            f"  libephemeris: {pos_py[0]:.6f} deg\n"
            f"  Distance from J2000: {dist_from_j2000} years"
        )

        # Latitude assertion (same tolerance band as longitude)
        assert diff_lat < tol_lon, (
            f"{body_name} MOSEPH at {date_desc}: "
            f"latitude diff {diff_lat:.6f} deg exceeds "
            f"{tol_lon} deg tolerance\n"
            f"  pyswisseph:   {pos_swe[1]:.6f} deg\n"
            f"  libephemeris: {pos_py[1]:.6f} deg\n"
            f"  Distance from J2000: {dist_from_j2000} years"
        )

        # Distance assertion (relative)
        if pos_swe[2] > 0:
            assert diff_dist_rel < dist_tol, (
                f"{body_name} MOSEPH at {date_desc}: "
                f"distance relative diff {diff_dist_rel:.6f} exceeds "
                f"{dist_tol} tolerance\n"
                f"  pyswisseph:   {pos_swe[2]:.8f} AU\n"
                f"  libephemeris: {pos_py[2]:.8f} AU\n"
                f"  Distance from J2000: {dist_from_j2000} years"
            )


class TestHistoricalPrecisionReport:
    """Produce a comprehensive precision matrix [planet x epoch].

    Runs all 80 combinations and collects max/mean errors per planet
    and per era, then prints a formatted report table. This test always
    passes (it is a reporting tool); actual assertions are in
    TestHistoricalPositions.
    """

    @pytest.mark.slow
    @pytest.mark.comparison
    def test_precision_matrix(self):
        """Generate and print the full [planet x epoch] precision matrix.

        Computes C-vs-Python differences for all 80 combinations and
        reports max/mean longitude error per planet per era, plus
        overall statistics per planet and per era.
        """
        flag = SEFLG_MOSEPH | SEFLG_SPEED

        # Collect results: {(planet_name, era_key): {"lon": [...], "lat": [...]}}
        results: dict[tuple[str, str], dict[str, list[float]]] = {}
        # Per-planet aggregates
        planet_errors: dict[str, list[float]] = {p[1]: [] for p in PLANETS}
        # Per-era aggregates
        era_errors: dict[str, list[float]] = {label: [] for label in ERA_LABELS}

        for body_id, body_name in PLANETS:
            for year, month, day, hour, _date_desc in HISTORICAL_DATES:
                jd = swe.julday(year, month, day, hour)

                pos_swe, _ = swe.calc_ut(jd, body_id, flag)
                pos_py, _ = ephem.swe_calc_ut(jd, body_id, flag)

                diff_lon = angular_diff(pos_swe[0], pos_py[0])
                diff_lat = abs(pos_swe[1] - pos_py[1])

                era = _era_label(year)
                key = (body_name, era)

                if key not in results:
                    results[key] = {"lon": [], "lat": []}

                results[key]["lon"].append(diff_lon)
                results[key]["lat"].append(diff_lat)

                planet_errors[body_name].append(diff_lon)
                era_errors[era].append(diff_lon)

        # Print precision matrix
        era_keys = list(ERA_LABELS.keys())

        print("\n" + "=" * 90)
        print(
            "MOSHIER HISTORICAL PRECISION MATRIX: "
            "C (pyswisseph) vs Python (libephemeris)"
        )
        print("=" * 90)

        # Header
        print(f"{'Planet':<12}", end="")
        for era in era_keys:
            print(f"  {ERA_LABELS[era]:>16}", end="")
        print(f"  {'Overall':>16}")
        print("-" * 90)

        # Body rows
        for _, body_name in PLANETS:
            print(f"{body_name:<12}", end="")
            for era in era_keys:
                key = (body_name, era)
                if key in results and results[key]["lon"]:
                    lons = results[key]["lon"]
                    max_err = max(lons)
                    mean_err = sum(lons) / len(lons)
                    print(f"  {max_err:7.4f}/{mean_err:7.4f}", end="")
                else:
                    print(f"  {'N/A':>16}", end="")

            # Overall for this planet
            all_errs = planet_errors[body_name]
            if all_errs:
                p_max = max(all_errs)
                p_mean = sum(all_errs) / len(all_errs)
                print(f"  {p_max:7.4f}/{p_mean:7.4f}")
            else:
                print(f"  {'N/A':>16}")

        print("-" * 90)

        # Era summary row
        print(f"{'Era max/avg':<12}", end="")
        for era in era_keys:
            errs = era_errors[era]
            if errs:
                print(
                    f"  {max(errs):7.4f}/{sum(errs) / len(errs):7.4f}",
                    end="",
                )
            else:
                print(f"  {'N/A':>16}", end="")
        print()

        print("=" * 90)
        print("Format: max_error / mean_error (degrees)")
        print(
            "Planet tolerance bands: "
            "<-1000: 1.0 deg, -1000-0: 0.5 deg, 0-1000: 0.1 deg, 1000+: 0.02 deg"
        )
        print("Moon: 2.0/1.0/0.5/0.15 deg  |  Pluto: 20.0/10.0/5.0/2.0 deg")
        print("=" * 90)
