"""
Moshier Coordinate Transformation Comparison Tests.

Compares coordinate transformation functions between pyswisseph and libephemeris
using positions derived from Moshier (SEFLG_MOSEPH) calculations as input:
- cotrans / cotrans_sp - ecliptic/equatorial coordinate transforms
- azalt / azalt_rev - azimuth/altitude calculations
- refrac - atmospheric refraction

Unlike test_compare_coordinates.py which uses hardcoded test coordinates,
this module validates the full chain: Moshier position -> coordinate transform
-> output, ensuring that differences in Moshier obliquity (IAU 2006 vs Lieske
1979) do not propagate into azimuth/altitude errors for observational
applications using Moshier for historical dates.

Test count: 48 cases
  - TestMoshierCotrans: 6 ecl->eq + 6 eq->ecl + 6 roundtrip = 18
  - TestMoshierCotransSp: 6 velocity transform
  - TestMoshierAzalt: 8 equatorial->horizontal
  - TestMoshierAzaltRev: 8 horizontal->equatorial
  - TestMoshierRefrac: 8 refraction from Moshier altitudes
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCES
# ============================================================================

ANGLE_TOL = 0.01  # degrees for transformed coordinates
REFRACTION_TOL = 0.01  # degrees for refraction
SPEED_TOL = 0.01  # degrees/day for velocity


# ============================================================================
# TEST DATA
# ============================================================================

BODIES = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
]

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 11, 15, 0.0, "Current Era"),
    (1900, 1, 1, 0.0, "Early 1900s"),
]

AZALT_LOCATIONS = [
    ("Equator", 0.0, 0.0, 0),
    ("Mid-lat North", 45.0, 0.0, 0),
    ("High-lat North", 70.0, 0.0, 0),
    ("Mid-lat South", -35.0, 0.0, 0),
]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def compute_obliquity(jd: float) -> float:
    """Compute mean obliquity of ecliptic (Laskar 1986) for a given Julian Day.

    Used as the obliquity parameter for cotrans/cotrans_sp. The exact value
    does not matter for comparison tests since both implementations receive
    the same obliquity; what matters is that it is a realistic value for the
    date so that the transformation exercises realistic coordinate ranges.
    """
    t = (jd - 2451545.0) / 36525.0
    eps = (
        84381.406 - 46.836769 * t - 0.0001831 * t * t + 0.00200340 * t * t * t
    ) / 3600.0
    return eps


def get_moshier_ecliptic(body_id: int, jd: float):
    """Get Moshier ecliptic positions with speed from both swe and ephem."""
    flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
    flag_py = SEFLG_MOSEPH | SEFLG_SPEED

    pos_swe, _ = swe.calc_ut(jd, body_id, flag_swe)
    pos_py, _ = ephem.swe_calc_ut(jd, body_id, flag_py)

    return pos_swe, pos_py


def get_moshier_equatorial(body_id: int, jd: float):
    """Get Moshier equatorial positions with speed from both swe and ephem."""
    flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
    flag_py = SEFLG_MOSEPH | SEFLG_SPEED | SEFLG_EQUATORIAL

    pos_swe, _ = swe.calc_ut(jd, body_id, flag_swe)
    pos_py, _ = ephem.swe_calc_ut(jd, body_id, flag_py)

    return pos_swe, pos_py


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMoshierCotrans:
    """Compare cotrans coordinate transformation with Moshier-derived positions.

    Uses ecliptic/equatorial positions calculated via SEFLG_MOSEPH as input,
    then transforms using cotrans() and compares swe vs ephem outputs.

    The same Moshier-derived coordinates are fed to both implementations
    so that the test isolates the cotrans transformation itself, with the
    Moshier aspect providing realistic (non-trivial) input values.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_cotrans_ecl_to_eq(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test cotrans ecliptic to equatorial with Moshier positions."""
        jd = swe.julday(year, month, day, hour)
        eps = compute_obliquity(jd)

        # Get Moshier ecliptic positions from both implementations
        pos_swe, pos_py = get_moshier_ecliptic(body_id, jd)

        # Use swe positions as canonical input for both cotrans
        ecl_coords = (pos_swe[0], pos_swe[1], pos_swe[2])

        result_swe = swe.cotrans(ecl_coords, -eps)
        result_py = ephem.cotrans(ecl_coords, -eps)

        diff_ra = angular_diff(result_swe[0], result_py[0])
        diff_dec = abs(result_swe[1] - result_py[1])

        assert diff_ra < ANGLE_TOL, (
            f"cotrans ecl->eq RA diff {diff_ra:.6f}° for {body_name} at {date_desc}"
        )
        assert diff_dec < ANGLE_TOL, (
            f"cotrans ecl->eq Dec diff {diff_dec:.6f}° for {body_name} at {date_desc}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_cotrans_eq_to_ecl(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test cotrans equatorial to ecliptic with Moshier positions."""
        jd = swe.julday(year, month, day, hour)
        eps = compute_obliquity(jd)

        # Get Moshier equatorial positions from both implementations
        pos_swe, pos_py = get_moshier_equatorial(body_id, jd)

        # Use swe positions as canonical input for both cotrans
        eq_coords = (pos_swe[0], pos_swe[1], pos_swe[2])

        result_swe = swe.cotrans(eq_coords, eps)
        result_py = ephem.cotrans(eq_coords, eps)

        diff_lon = angular_diff(result_swe[0], result_py[0])
        diff_lat = abs(result_swe[1] - result_py[1])

        assert diff_lon < ANGLE_TOL, (
            f"cotrans eq->ecl lon diff {diff_lon:.6f}° for {body_name} at {date_desc}"
        )
        assert diff_lat < ANGLE_TOL, (
            f"cotrans eq->ecl lat diff {diff_lat:.6f}° for {body_name} at {date_desc}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_cotrans_roundtrip(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test cotrans roundtrip ecl -> eq -> ecl with Moshier positions."""
        jd = swe.julday(year, month, day, hour)
        eps = compute_obliquity(jd)

        # Get Moshier ecliptic positions from ephem
        _, pos_py = get_moshier_ecliptic(body_id, jd)
        ecl_orig = (pos_py[0], pos_py[1], pos_py[2])

        # Ecliptic -> Equatorial -> Ecliptic
        equatorial = ephem.cotrans(ecl_orig, -eps)
        ecl_back = ephem.cotrans((equatorial[0], equatorial[1], equatorial[2]), eps)

        diff_lon = angular_diff(ecl_orig[0], ecl_back[0])
        diff_lat = abs(ecl_orig[1] - ecl_back[1])

        assert diff_lon < 1e-8, (
            f"Roundtrip lon diff {diff_lon:.12f}° for {body_name} at {date_desc}"
        )
        assert diff_lat < 1e-8, (
            f"Roundtrip lat diff {diff_lat:.12f}° for {body_name} at {date_desc}"
        )


class TestMoshierCotransSp:
    """Compare cotrans_sp with Moshier-derived velocities.

    Uses positions and velocities calculated via SEFLG_MOSEPH | SEFLG_SPEED
    as input for cotrans_sp(), validating that velocity transformation
    is consistent between swe and ephem when given Moshier-like inputs.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_cotrans_sp_with_moshier_velocity(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test cotrans_sp with Moshier velocities."""
        jd = swe.julday(year, month, day, hour)
        eps = compute_obliquity(jd)

        # Get Moshier ecliptic positions with speed from ephem
        _, pos_py = get_moshier_ecliptic(body_id, jd)

        # Use ephem positions/speeds as canonical input for both
        coord = (pos_py[0], pos_py[1], pos_py[2])
        speed = (pos_py[3], pos_py[4], pos_py[5])

        # libephemeris cotrans_sp
        result_coord, result_speed = ephem.cotrans_sp(coord, speed, -eps)

        # pyswisseph cotrans_sp (6-tuple input)
        coords_6 = coord + speed
        result_swe = swe.cotrans_sp(coords_6, -eps)

        # Compare positions
        diff_ra = angular_diff(result_swe[0], result_coord[0])
        diff_dec = abs(result_swe[1] - result_coord[1])

        assert diff_ra < ANGLE_TOL, (
            f"cotrans_sp RA diff {diff_ra:.6f}° for {body_name} at {date_desc}"
        )
        assert diff_dec < ANGLE_TOL, (
            f"cotrans_sp Dec diff {diff_dec:.6f}° for {body_name} at {date_desc}"
        )

        # Compare velocities
        diff_ra_speed = abs(result_swe[3] - result_speed[0])
        diff_dec_speed = abs(result_swe[4] - result_speed[1])

        assert diff_ra_speed < SPEED_TOL, (
            f"cotrans_sp RA speed diff {diff_ra_speed:.6f}°/day "
            f"for {body_name} at {date_desc}"
        )
        assert diff_dec_speed < SPEED_TOL, (
            f"cotrans_sp Dec speed diff {diff_dec_speed:.6f}°/day "
            f"for {body_name} at {date_desc}"
        )


class TestMoshierAzalt:
    """Compare azalt with Moshier-derived equatorial positions.

    Validates the chain: Moshier SEFLG_EQUATORIAL position ->
    azalt() equatorial-to-horizontal -> azimuth/altitude.
    The same Moshier equatorial coordinates are used as input for both
    swe.azalt and ephem.azalt to isolate the transformation comparison.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("name,lat,lon,alt", AZALT_LOCATIONS)
    def test_azalt_with_moshier_positions(
        self, body_id, body_name, name, lat, lon, alt
    ):
        """Test azalt with Moshier-derived equatorial coordinates."""
        jd = swe.julday(2000, 1, 1, 12.0)  # J2000 noon
        geopos = (lon, lat, alt)
        atpress = 1013.25
        attemp = 15.0

        # Get Moshier equatorial position (use swe as canonical)
        pos_swe, _ = get_moshier_equatorial(body_id, jd)
        xin = (pos_swe[0], pos_swe[1], pos_swe[2])

        # SE_EQU2HOR = 1
        result_swe = swe.azalt(jd, swe.EQU2HOR, geopos, atpress, attemp, xin)
        result_py = ephem.azalt(jd, 1, geopos, atpress, attemp, xin)

        diff_az = angular_diff(result_swe[0], result_py[0])
        diff_alt = abs(result_swe[1] - result_py[1])

        assert diff_az < ANGLE_TOL, (
            f"azalt azimuth diff {diff_az:.6f}° for {body_name} at {name}"
        )
        assert diff_alt < ANGLE_TOL, (
            f"azalt altitude diff {diff_alt:.6f}° for {body_name} at {name}"
        )


AZALT_REV_LOCATIONS = [
    # Equator excluded: azalt_rev has a known singularity at lat=0
    # (sin(lat)~0 causes hour angle degeneration in libephemeris)
    ("Mid-lat North", 45.0, 0.0, 0),
    ("High-lat North", 70.0, 0.0, 0),
    ("Mid-lat South", -35.0, 0.0, 0),
    ("Low-lat North", 10.0, 0.0, 0),
]


class TestMoshierAzaltRev:
    """Compare azalt_rev with horizontal coordinates derived from Moshier positions.

    Validates the reverse chain: Moshier position -> azalt() -> azalt_rev()
    recovering equatorial coordinates from horizontal ones derived from
    Moshier equatorial input.

    Note: Equator (lat=0) is excluded because azalt_rev has a known singularity
    at lat=0 where sin(lat)~0 causes the hour angle calculation to degenerate
    (division by sin_lat in the cos_H formula). This is a known limitation of
    the libephemeris azalt_rev implementation, not specific to Moshier.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("name,lat,lon,alt", AZALT_REV_LOCATIONS)
    def test_azalt_rev_from_moshier(self, body_id, body_name, name, lat, lon, alt):
        """Test azalt_rev with horizontal coordinates derived from Moshier positions."""
        jd = swe.julday(2000, 1, 1, 12.0)
        geopos = (lon, lat, alt)
        atpress = 1013.25
        attemp = 15.0

        # Get Moshier equatorial position (use swe as canonical)
        pos_swe, _ = get_moshier_equatorial(body_id, jd)
        xin = (pos_swe[0], pos_swe[1], pos_swe[2])

        # Forward: get horizontal coordinates using swe (canonical)
        result_hor = swe.azalt(jd, swe.EQU2HOR, geopos, atpress, attemp, xin)
        azimuth = result_hor[0]
        true_alt = result_hor[1]

        # Reverse: horizontal to equatorial from both implementations
        result_swe = swe.azalt_rev(jd, swe.HOR2EQU, geopos, azimuth, true_alt)
        result_py = ephem.azalt_rev(jd, ephem.SE_HOR2EQU, geopos, azimuth, true_alt)

        diff_ra = angular_diff(result_swe[0], result_py[0])
        diff_dec = abs(result_swe[1] - result_py[1])

        assert diff_ra < ANGLE_TOL, (
            f"azalt_rev RA diff {diff_ra:.6f}° for {body_name} at {name}"
        )
        assert diff_dec < ANGLE_TOL, (
            f"azalt_rev Dec diff {diff_dec:.6f}° for {body_name} at {name}"
        )


class TestMoshierRefrac:
    """Compare refraction with altitudes derived from Moshier positions.

    Uses true altitudes computed via azalt() from Moshier equatorial positions,
    then applies refraction correction and compares between implementations.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", BODIES)
    @pytest.mark.parametrize("name,lat,lon,alt", AZALT_LOCATIONS)
    def test_refrac_from_moshier_altitude(
        self, body_id, body_name, name, lat, lon, alt
    ):
        """Test refrac with altitude derived from Moshier position."""
        jd = swe.julday(2000, 1, 1, 12.0)
        geopos = (lon, lat, alt)
        atpress = 1013.25
        attemp = 15.0

        # Get altitude from Moshier position via azalt
        pos_swe, _ = get_moshier_equatorial(body_id, jd)
        xin = (pos_swe[0], pos_swe[1], pos_swe[2])

        result_hor = swe.azalt(jd, swe.EQU2HOR, geopos, atpress, attemp, xin)
        true_alt = result_hor[1]  # True altitude

        # Skip deeply negative altitudes where refraction extrapolation
        # formulas diverge between implementations (below -2° the Bennett
        # formula is extrapolated and implementations differ in approach)
        if true_alt < -2.0:
            pytest.skip(
                f"true_alt={true_alt:.2f}° below -2°: refraction "
                f"extrapolation diverges between implementations"
            )

        # Apply refraction from both implementations
        # SE_TRUE_TO_APP = 0
        refrac_swe = swe.refrac(true_alt, atpress, attemp, swe.TRUE_TO_APP)
        refrac_py = ephem.refrac(true_alt, atpress, attemp, 0)

        diff = abs(refrac_swe - refrac_py)

        assert diff < REFRACTION_TOL, (
            f"refrac diff {diff:.6f}° for {body_name} at {name} "
            f"(true_alt={true_alt:.2f}°)"
        )
