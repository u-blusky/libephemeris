"""
Moshier Historical Birth Charts End-to-End Comparison Tests.

Validates complete horoscope calculations (10 planets + 4 lunar points +
12 house cusps + Ascendant + MC) for 5 famous people between pyswisseph
(C library) and libephemeris (Python reimplementation) using SEFLG_MOSEPH.

This is the gold standard for end-to-end Moshier validation: it replicates
exactly what an astrological software does when computing a natal chart.
If any component differs between C and Python (e.g. Ascendant in a different
sign, Moon in a different house), the user perceives the software as "wrong".

Famous people tested (from conftest.py famous_people fixture):
- Albert Einstein (1879-03-14, Ulm, Germany)
- Marilyn Monroe (1926-06-01, Los Angeles, USA)
- Mahatma Gandhi (1869-10-02, Porbandar, India)
- Pablo Picasso (1881-10-25, Malaga, Spain)
- Nelson Mandela (1918-07-18, Mvezo, South Africa)

Total validations: 5 persons x ~26 components = ~130 cross-library comparisons.
"""

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
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
)


# ============================================================================
# TOLERANCES
# ============================================================================

# Planet longitude tolerances (consistent with test_moshier_compare_planets.py)
# Birth dates span 1869-1926, i.e. 74-131 years from J2000, so we use the
# MOSHIER_EXTENDED_RANGE (0.03°) from test_moshier_compare_planets.py rather
# than the baseline 0.02° which is calibrated on dates near J2000.
PLANET_LON_TOL = 0.03  # degrees (~108 arcsec, extended range for historical dates)
PLUTO_LON_TOL = 2.0  # degrees (Chapront-Francou theory has large errors)

# Planet velocity tolerances
PLANET_SPEED_TOL = 0.02  # degrees/day
MOON_SPEED_TOL = 0.05  # degrees/day (Moon's high angular velocity)

# Lunar point tolerances (consistent with test_moshier_compare_lunar.py)
MEAN_NODE_TOL = 0.02  # degrees
TRUE_NODE_TOL = 0.2  # degrees
MEAN_LILITH_TOL = 0.02  # degrees
TRUE_LILITH_TOL = 0.2  # degrees

# House cusp tolerances (consistent with test_moshier_compare_houses.py)
CUSP_TOL = 0.02  # degrees (Placidus iterative algorithm amplifies diffs)
ASCMC_TOL = 0.02  # degrees


# ============================================================================
# CHART COMPONENTS
# ============================================================================

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

LUNAR_POINTS = [
    (SE_MEAN_NODE, "Mean Node", MEAN_NODE_TOL),
    (SE_TRUE_NODE, "True Node", TRUE_NODE_TOL),
    (SE_MEAN_APOG, "Mean Lilith", MEAN_LILITH_TOL),
    (SE_OSCU_APOG, "True Lilith", TRUE_LILITH_TOL),
]


# ============================================================================
# HELPERS
# ============================================================================


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def get_planet_lon_tol(planet_id: int) -> float:
    """Get longitude tolerance for a given planet."""
    if planet_id == SE_PLUTO:
        return PLUTO_LON_TOL
    return PLANET_LON_TOL


def get_planet_speed_tol(planet_id: int) -> float:
    """Get speed tolerance for a given planet."""
    if planet_id == SE_MOON:
        return MOON_SPEED_TOL
    return PLANET_SPEED_TOL


def find_house(lon: float, cusps: list) -> int:
    """Determine which house a given ecliptic longitude falls in.

    Args:
        lon: Ecliptic longitude in degrees [0, 360).
        cusps: List of 12 house cusp longitudes (0-indexed, cusps[0] = cusp 1).

    Returns:
        House number (1-12).
    """
    for i in range(12):
        cusp_start = cusps[i]
        cusp_end = cusps[(i + 1) % 12]

        if cusp_start < cusp_end:
            # Normal case: cusp range doesn't cross 0 degrees
            if cusp_start <= lon < cusp_end:
                return i + 1
        else:
            # Wrap-around: cusp range crosses 0 degrees (e.g., 350 to 10)
            if lon >= cusp_start or lon < cusp_end:
                return i + 1

    # Fallback (shouldn't happen with valid cusps)
    return 12


# ============================================================================
# TEST CLASS
# ============================================================================


class TestHistoricalBirthCharts:
    """End-to-end comparison of complete birth charts using Moshier.

    For each famous person: compute JD, then calculate all 10 planets +
    4 lunar points with SEFLG_MOSEPH | SEFLG_SPEED, then 12 Placidus cusps
    with swe_houses_ex(). Compare every component between C and Python.
    """

    @pytest.mark.comparison
    @pytest.mark.integration
    def test_full_chart_planets(self, famous_people):
        """Test all 10 planets for all 5 famous people.

        Validates longitude and speed for each planet in each birth chart.
        Total: 5 people x 10 planets = 50 planet comparisons.
        """
        flags_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flags_py = SEFLG_MOSEPH | SEFLG_SPEED

        failures = []

        for person in famous_people:
            name = person["name"]
            jd = swe.julday(
                person["year"], person["month"], person["day"], person["hour"]
            )

            for planet_id, planet_name in PLANETS:
                res_swe, _ = swe.calc_ut(jd, planet_id, flags_swe)
                res_py, _ = ephem.swe_calc_ut(jd, planet_id, flags_py)

                lon_diff = angular_diff(res_swe[0], res_py[0])
                lon_tol = get_planet_lon_tol(planet_id)

                if lon_diff >= lon_tol:
                    failures.append(
                        f"  {name} {planet_name}: lon diff {lon_diff:.6f}° "
                        f"(tol {lon_tol}°, swe={res_swe[0]:.4f}°, "
                        f"py={res_py[0]:.4f}°)"
                    )

                speed_diff = abs(res_swe[3] - res_py[3])
                speed_tol = get_planet_speed_tol(planet_id)

                if speed_diff >= speed_tol:
                    failures.append(
                        f"  {name} {planet_name}: speed diff {speed_diff:.6f}°/d "
                        f"(tol {speed_tol}°/d, swe={res_swe[3]:.6f}, "
                        f"py={res_py[3]:.6f})"
                    )

        assert len(failures) == 0, (
            f"Planet comparison failures ({len(failures)}):\n" + "\n".join(failures)
        )

    @pytest.mark.comparison
    @pytest.mark.integration
    def test_full_chart_lunar_points(self, famous_people):
        """Test all 4 lunar points for all 5 famous people.

        Validates longitude for Mean/True Node and Mean/True Lilith.
        Total: 5 people x 4 lunar points = 20 comparisons.
        """
        failures = []

        for person in famous_people:
            name = person["name"]
            jd = swe.julday(
                person["year"], person["month"], person["day"], person["hour"]
            )

            for point_id, point_name, tolerance in LUNAR_POINTS:
                pos_swe, _ = swe.calc_ut(jd, point_id, swe.FLG_MOSEPH)
                pos_py, _ = ephem.swe_calc_ut(jd, point_id, SEFLG_MOSEPH)

                diff = angular_diff(pos_swe[0], pos_py[0])

                if diff >= tolerance:
                    failures.append(
                        f"  {name} {point_name}: lon diff {diff:.6f}° "
                        f"(tol {tolerance}°, swe={pos_swe[0]:.4f}°, "
                        f"py={pos_py[0]:.4f}°)"
                    )

        assert len(failures) == 0, (
            f"Lunar point comparison failures ({len(failures)}):\n"
            + "\n".join(failures)
        )

    @pytest.mark.comparison
    @pytest.mark.integration
    def test_full_chart_house_cusps(self, famous_people):
        """Test all 12 Placidus house cusps for all 5 famous people.

        Validates cusp longitudes using Placidus (the most common system).
        Total: 5 people x 12 cusps = 60 cusp comparisons.
        """
        failures = []

        for person in famous_people:
            name = person["name"]
            jd = swe.julday(
                person["year"], person["month"], person["day"], person["hour"]
            )
            lat = person["lat"]
            lon = person["lon"]

            cusps_swe, _ = swe.houses_ex(jd, lat, lon, b"P", swe.FLG_MOSEPH)
            cusps_py, _ = ephem.swe_houses_ex(jd, lat, lon, "P", SEFLG_MOSEPH)

            for i in range(12):
                diff = angular_diff(cusps_swe[i], cusps_py[i])

                if diff >= CUSP_TOL:
                    failures.append(
                        f"  {name} cusp {i + 1}: diff {diff:.6f}° "
                        f"(tol {CUSP_TOL}°, swe={cusps_swe[i]:.4f}°, "
                        f"py={cusps_py[i]:.4f}°)"
                    )

        assert len(failures) == 0, (
            f"House cusp comparison failures ({len(failures)}):\n" + "\n".join(failures)
        )

    @pytest.mark.comparison
    @pytest.mark.integration
    def test_full_chart_ascendant_mc(self, famous_people):
        """Test Ascendant and MC for all 5 famous people.

        These are the two most critical angles in any birth chart.
        Total: 5 people x 2 angles = 10 comparisons.
        """
        failures = []

        for person in famous_people:
            name = person["name"]
            jd = swe.julday(
                person["year"], person["month"], person["day"], person["hour"]
            )
            lat = person["lat"]
            lon = person["lon"]

            _, ascmc_swe = swe.houses_ex(jd, lat, lon, b"P", swe.FLG_MOSEPH)
            _, ascmc_py = ephem.swe_houses_ex(jd, lat, lon, "P", SEFLG_MOSEPH)

            asc_diff = angular_diff(ascmc_swe[0], ascmc_py[0])
            if asc_diff >= ASCMC_TOL:
                failures.append(
                    f"  {name} ASC: diff {asc_diff:.6f}° "
                    f"(tol {ASCMC_TOL}°, swe={ascmc_swe[0]:.4f}°, "
                    f"py={ascmc_py[0]:.4f}°)"
                )

            mc_diff = angular_diff(ascmc_swe[1], ascmc_py[1])
            if mc_diff >= ASCMC_TOL:
                failures.append(
                    f"  {name} MC: diff {mc_diff:.6f}° "
                    f"(tol {ASCMC_TOL}°, swe={ascmc_swe[1]:.4f}°, "
                    f"py={ascmc_py[1]:.4f}°)"
                )

        assert len(failures) == 0, (
            f"ASC/MC comparison failures ({len(failures)}):\n" + "\n".join(failures)
        )

    @pytest.mark.comparison
    @pytest.mark.integration
    def test_full_chart_house_placement(self, famous_people):
        """Test that all planets fall in the same house in both implementations.

        This is the ultimate user-facing validation: if the C library places
        the Moon in house 4 but the Python library places it in house 5,
        the user sees a completely different chart. We verify that every
        planet (using C-library positions and cusps as reference) is placed
        in the same house by both implementations.

        Total: 5 people x 10 planets = 50 house placement checks.
        """
        flags_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flags_py = SEFLG_MOSEPH | SEFLG_SPEED

        failures = []

        for person in famous_people:
            name = person["name"]
            jd = swe.julday(
                person["year"], person["month"], person["day"], person["hour"]
            )
            lat = person["lat"]
            lon = person["lon"]

            # Get house cusps from both implementations
            cusps_swe, _ = swe.houses_ex(jd, lat, lon, b"P", swe.FLG_MOSEPH)
            cusps_py, _ = ephem.swe_houses_ex(jd, lat, lon, "P", SEFLG_MOSEPH)

            for planet_id, planet_name in PLANETS:
                res_swe, _ = swe.calc_ut(jd, planet_id, flags_swe)
                res_py, _ = ephem.swe_calc_ut(jd, planet_id, flags_py)

                # Determine house placement using each implementation's
                # own cusps and planet positions
                house_swe = find_house(res_swe[0], list(cusps_swe[:12]))
                house_py = find_house(res_py[0], list(cusps_py[:12]))

                if house_swe != house_py:
                    failures.append(
                        f"  {name} {planet_name}: C=house {house_swe} vs "
                        f"Py=house {house_py} "
                        f"(lon_swe={res_swe[0]:.4f}°, lon_py={res_py[0]:.4f}°)"
                    )

        assert len(failures) == 0, (
            f"House placement mismatches ({len(failures)}):\n" + "\n".join(failures)
        )


# ============================================================================
# PARAMETRIZED PER-PERSON TESTS
# ============================================================================

# Birth data for parametrize (mirrors conftest.py famous_people fixture)
FAMOUS_PEOPLE = [
    ("Einstein", 1879, 3, 14, 11.5, 48.4010, 9.9876),
    ("Monroe", 1926, 6, 1, 9.5, 34.0522, -118.2437),
    ("Gandhi", 1869, 10, 2, 7.2, 21.6417, 69.6293),
    ("Picasso", 1881, 10, 25, 23.25, 36.7213, -4.4214),
    ("Mandela", 1918, 7, 18, 14.0, -31.5875, 28.7833),
]


class TestPerPersonPlanets:
    """Parametrized per-person, per-planet tests for clear failure reporting."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize(
        "person,year,month,day,hour,lat,lon",
        FAMOUS_PEOPLE,
        ids=[p[0] for p in FAMOUS_PEOPLE],
    )
    def test_planet_longitude(
        self, person, year, month, day, hour, lat, lon, planet_id, planet_name
    ):
        """Test Moshier planet longitude for a famous person's birth chart."""
        jd = swe.julday(year, month, day, hour)
        flags_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flags_py = SEFLG_MOSEPH | SEFLG_SPEED

        res_swe, _ = swe.calc_ut(jd, planet_id, flags_swe)
        res_py, _ = ephem.swe_calc_ut(jd, planet_id, flags_py)

        diff = angular_diff(res_swe[0], res_py[0])
        tol = get_planet_lon_tol(planet_id)

        assert diff < tol, (
            f"{person} {planet_name}: lon diff {diff:.6f}° "
            f"(tol {tol}°, swe={res_swe[0]:.4f}°, py={res_py[0]:.4f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize(
        "person,year,month,day,hour,lat,lon",
        FAMOUS_PEOPLE,
        ids=[p[0] for p in FAMOUS_PEOPLE],
    )
    def test_planet_speed(
        self, person, year, month, day, hour, lat, lon, planet_id, planet_name
    ):
        """Test Moshier planet speed for a famous person's birth chart."""
        jd = swe.julday(year, month, day, hour)
        flags_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flags_py = SEFLG_MOSEPH | SEFLG_SPEED

        res_swe, _ = swe.calc_ut(jd, planet_id, flags_swe)
        res_py, _ = ephem.swe_calc_ut(jd, planet_id, flags_py)

        diff = abs(res_swe[3] - res_py[3])
        tol = get_planet_speed_tol(planet_id)

        assert diff < tol, (
            f"{person} {planet_name}: speed diff {diff:.6f}°/d "
            f"(tol {tol}°/d, swe={res_swe[3]:.6f}, py={res_py[3]:.6f})"
        )


class TestPerPersonLunarPoints:
    """Parametrized per-person, per-lunar-point tests."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "point_id,point_name,tolerance",
        LUNAR_POINTS,
        ids=[lp[1] for lp in LUNAR_POINTS],
    )
    @pytest.mark.parametrize(
        "person,year,month,day,hour,lat,lon",
        FAMOUS_PEOPLE,
        ids=[p[0] for p in FAMOUS_PEOPLE],
    )
    def test_lunar_point_longitude(
        self, person, year, month, day, hour, lat, lon, point_id, point_name, tolerance
    ):
        """Test Moshier lunar point longitude for a famous person's birth chart."""
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, point_id, swe.FLG_MOSEPH)
        pos_py, _ = ephem.swe_calc_ut(jd, point_id, SEFLG_MOSEPH)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tolerance, (
            f"{person} {point_name}: lon diff {diff:.6f}° "
            f"(tol {tolerance}°, swe={pos_swe[0]:.4f}°, py={pos_py[0]:.4f}°)"
        )


class TestPerPersonHouseCusps:
    """Parametrized per-person house cusp and angle tests."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "person,year,month,day,hour,lat,lon",
        FAMOUS_PEOPLE,
        ids=[p[0] for p in FAMOUS_PEOPLE],
    )
    def test_placidus_cusps(self, person, year, month, day, hour, lat, lon):
        """Test all 12 Placidus cusps for a famous person's birth chart."""
        jd = swe.julday(year, month, day, hour)

        cusps_swe, _ = swe.houses_ex(jd, lat, lon, b"P", swe.FLG_MOSEPH)
        cusps_py, _ = ephem.swe_houses_ex(jd, lat, lon, "P", SEFLG_MOSEPH)

        max_diff = 0.0
        worst_cusp = 0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            if diff > max_diff:
                max_diff = diff
                worst_cusp = i + 1

        assert max_diff < CUSP_TOL, (
            f"{person} Placidus: max cusp diff {max_diff:.6f}° at cusp {worst_cusp} "
            f"(tol {CUSP_TOL}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "person,year,month,day,hour,lat,lon",
        FAMOUS_PEOPLE,
        ids=[p[0] for p in FAMOUS_PEOPLE],
    )
    def test_ascendant(self, person, year, month, day, hour, lat, lon):
        """Test Ascendant for a famous person's birth chart."""
        jd = swe.julday(year, month, day, hour)

        _, ascmc_swe = swe.houses_ex(jd, lat, lon, b"P", swe.FLG_MOSEPH)
        _, ascmc_py = ephem.swe_houses_ex(jd, lat, lon, "P", SEFLG_MOSEPH)

        diff = angular_diff(ascmc_swe[0], ascmc_py[0])

        assert diff < ASCMC_TOL, (
            f"{person} ASC: diff {diff:.6f}° "
            f"(tol {ASCMC_TOL}°, swe={ascmc_swe[0]:.4f}°, py={ascmc_py[0]:.4f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "person,year,month,day,hour,lat,lon",
        FAMOUS_PEOPLE,
        ids=[p[0] for p in FAMOUS_PEOPLE],
    )
    def test_mc(self, person, year, month, day, hour, lat, lon):
        """Test MC (Medium Coeli) for a famous person's birth chart."""
        jd = swe.julday(year, month, day, hour)

        _, ascmc_swe = swe.houses_ex(jd, lat, lon, b"P", swe.FLG_MOSEPH)
        _, ascmc_py = ephem.swe_houses_ex(jd, lat, lon, "P", SEFLG_MOSEPH)

        diff = angular_diff(ascmc_swe[1], ascmc_py[1])

        assert diff < ASCMC_TOL, (
            f"{person} MC: diff {diff:.6f}° "
            f"(tol {ASCMC_TOL}°, swe={ascmc_swe[1]:.4f}°, py={ascmc_py[1]:.4f}°)"
        )
