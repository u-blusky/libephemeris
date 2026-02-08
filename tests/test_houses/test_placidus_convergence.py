"""
Tests for Placidus house convergence accuracy.

These tests verify that the Placidus house system iteration converges
to match Swiss Ephemeris within 0.0001° (0.36 arcseconds).
"""

import pytest
import swisseph as swe
import libephemeris as ephem


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# High precision tolerance: 0.0001° = 0.36 arcseconds
HIGH_PRECISION_TOL = 0.0001


class TestPlacidusConvergence:
    """Tests for Placidus house cusp convergence to Swiss Ephemeris precision."""

    @pytest.mark.parametrize(
        "name,jd,lat,lon",
        [
            # john_lennon (1940-10-09 07:00:00, Liverpool)
            ("john_lennon", 2428867.7916666665, 53.4084, -2.9916),
            # paul_mccartney (1942-06-18 14:30:00, Liverpool)
            ("paul_mccartney", 2431259.4375, 53.4084, -2.9916),
            # johnny_depp (1963-06-09 03:45:00, San Francisco)
            ("johnny_depp", 2439606.15625, 37.7749, -122.4194),
        ],
    )
    def test_placidus_high_precision_match(self, name, jd, lat, lon):
        """Test that Placidus cusps match Swiss Ephemeris within 0.0001°."""
        cusps_swe, _ = swe.houses(jd, lat, lon, b"P")
        cusps_py, _ = ephem.swe_houses(jd, lat, lon, ord("P"))

        max_diff = 0.0
        worst_cusp = 0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            if diff > max_diff:
                max_diff = diff
                worst_cusp = i + 1

        assert max_diff < HIGH_PRECISION_TOL, (
            f"{name}: Placidus cusp {worst_cusp} differs by {max_diff:.6f}° "
            f"(target: <{HIGH_PRECISION_TOL}°)"
        )

    @pytest.mark.parametrize(
        "name,lat,lon",
        [
            ("Rome", 41.9028, 12.4964),
            ("New York", 40.7128, -74.0060),
            ("Sydney", -33.8688, 151.2093),
            ("Equator", 0.0, 0.0),
            ("London", 51.5074, -0.1278),
            ("Tokyo", 35.6762, 139.6503),
        ],
    )
    @pytest.mark.parametrize(
        "year,month,day,hour",
        [
            (2000, 1, 1, 12.0),  # J2000
            (2024, 6, 15, 0.0),  # Summer 2024
            (1980, 5, 20, 14.5),  # Past
            (1940, 10, 9, 7.0),  # john_lennon era
        ],
    )
    def test_placidus_convergence_across_dates(
        self, name, lat, lon, year, month, day, hour
    ):
        """Test Placidus convergence across various dates and locations."""
        jd = swe.julday(year, month, day, hour)

        cusps_swe, _ = swe.houses(jd, lat, lon, b"P")
        cusps_py, _ = ephem.swe_houses(jd, lat, lon, ord("P"))

        max_diff = 0.0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < HIGH_PRECISION_TOL, (
            f"Placidus at {name} ({year}-{month:02d}-{day:02d}): "
            f"max diff {max_diff:.6f}° exceeds {HIGH_PRECISION_TOL}°"
        )


class TestPlacidusIntermediateCusps:
    """
    Test intermediate cusps (2, 3, 5, 6, 8, 9, 11, 12) specifically.

    These cusps require iterative calculation and are most affected by
    the convergence algorithm's precision.
    """

    @pytest.mark.parametrize(
        "name,jd,lat,lon",
        [
            ("john_lennon", 2428867.7916666665, 53.4084, -2.9916),
            ("paul_mccartney", 2431259.4375, 53.4084, -2.9916),
        ],
    )
    def test_intermediate_cusps_precision(self, name, jd, lat, lon):
        """Test that intermediate cusps (calculated iteratively) are precise."""
        cusps_swe, _ = swe.houses(jd, lat, lon, b"P")
        cusps_py, _ = ephem.swe_houses(jd, lat, lon, ord("P"))

        # Intermediate cusps (not Asc/IC/Dsc/MC)
        intermediate_cusps = [2, 3, 5, 6, 8, 9, 11, 12]

        for cusp_num in intermediate_cusps:
            i = cusp_num - 1  # 0-indexed
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            assert diff < HIGH_PRECISION_TOL, (
                f"{name}: Cusp {cusp_num} diff {diff:.6f}° exceeds {HIGH_PRECISION_TOL}°"
            )
