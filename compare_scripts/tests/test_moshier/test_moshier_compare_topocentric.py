"""
Moshier Topocentric Cross-Library Comparison Tests.

Compares topocentric (SEFLG_TOPOCTR) planetary calculations between pyswisseph
(C library) and libephemeris (Python reimplementation) in Moshier mode.

PROBLEM: At planets.py:1045-1049, the topocentric correction in Moshier mode
has a comment `# Moon parallax is significant (~1 degree max)` followed by
`pass  # TODO: Implement proper topocentric correction if needed`. This means
SEFLG_TOPOCTR is silently ignored in Moshier mode. The C library applies the
parallax correction; for the Moon, the difference is up to ~1 degree.

IMPACT: The Moon parallax of ~1 degree renders the Moon position unusable for
local eclipses, occultations, rise/set from the Earth's surface, and any
application requiring topocentric position. Applications migrating from
pyswisseph with FLG_TOPOCTR | FLG_MOSEPH will get Moon positions wrong by
up to 1 degree WITHOUT ANY WARNING. This is the most severe violation of
1:1 compatibility.

These tests document the discrepancy with xfail markers, serving as a roadmap
for implementation and as a warning for users.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_TOPOCTR,
)


# ============================================================================
# CONSTANTS
# ============================================================================

# Observer position: Rome, Italy (lat=41.9, lon=12.5, alt=0m)
TOPO_LON = 12.5
TOPO_LAT = 41.9
TOPO_ALT = 0.0

# Three test dates spanning different eras
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 6, 15, 20.0, "Summer 2024 evening"),
    (1985, 3, 21, 6.0, "Vernal equinox 1985"),
]

# Bodies to test
BODIES = [
    (SE_MOON, "Moon"),
    (SE_SUN, "Sun"),
]

# Moshier C-vs-Python tolerances (same as test_moshier_compare_planets.py)
MOSHIER_LONGITUDE = 0.02  # degrees (~72 arcsec)
MOSHIER_LATITUDE = 0.02  # degrees
MOSHIER_DISTANCE_REL = 0.001  # relative (dimensionless)


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
# MODULE-LEVEL XFAIL MARKER
# ============================================================================

# All tests in this module are expected to fail because SEFLG_TOPOCTR is
# silently ignored in libephemeris Moshier mode (planets.py:1049 TODO).
# strict=False allows tests to pass (xpass) if the correction is later
# implemented, without breaking the test suite.
pytestmark = pytest.mark.xfail(
    reason=(
        "SEFLG_TOPOCTR silently ignored in Moshier mode "
        "(planets.py:1049 TODO: topocentric correction not implemented)"
    ),
    strict=False,
)


# ============================================================================
# TEST CLASS
# ============================================================================


class TestMoshierTopocentric:
    """Compare Moshier topocentric calculations between pyswisseph and libephemeris.

    The C Swiss Ephemeris applies topocentric (parallax) correction in Moshier
    mode via SEFLG_TOPOCTR. Libephemeris silently ignores this flag
    (planets.py:1045-1049), returning geocentric positions instead.

    For the Moon, the parallax can reach ~1 degree, making this the most
    significant Moshier compatibility gap. For the Sun, the parallax is
    ~8.8 arcsec (negligible at our tolerance level).

    Each test sets the observer position via set_topo() in both libraries
    before computing topocentric positions with SEFLG_MOSEPH | SEFLG_TOPOCTR.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_topocentric_moon(self, year, month, day, hour, date_desc):
        """Test Moshier topocentric Moon position matches pyswisseph.

        Expected to fail: Moon parallax (~1 degree max) is not applied by
        libephemeris in Moshier mode. The C library computes proper
        topocentric correction using the observer's geographic position.

        This discrepancy makes the Moon position unusable for:
        - Local eclipse predictions
        - Occultation calculations
        - Rise/set times from Earth's surface
        - Any application requiring apparent position from a specific location
        """
        jd = swe.julday(year, month, day, hour)
        flags_swe = swe.FLG_MOSEPH | swe.FLG_TOPOCTR | swe.FLG_SPEED
        flags_py = SEFLG_MOSEPH | SEFLG_TOPOCTR | SEFLG_SPEED

        swe.set_topo(TOPO_LON, TOPO_LAT, TOPO_ALT)
        ephem.swe_set_topo(TOPO_LON, TOPO_LAT, TOPO_ALT)

        pos_swe, _ = swe.calc_ut(jd, SE_MOON, flags_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MOON, flags_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        assert diff_lon < MOSHIER_LONGITUDE, (
            f"Moshier Topocentric Moon @ {date_desc}: "
            f"longitude diff {diff_lon:.6f}° exceeds tolerance {MOSHIER_LONGITUDE}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°). "
            f"This is caused by missing parallax correction in planets.py:1049."
        )
        assert diff_lat < MOSHIER_LATITUDE, (
            f"Moshier Topocentric Moon @ {date_desc}: "
            f"latitude diff {diff_lat:.6f}° exceeds tolerance {MOSHIER_LATITUDE}° "
            f"(swe={pos_swe[1]:.6f}°, lib={pos_py[1]:.6f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_topocentric_sun(self, year, month, day, hour, date_desc):
        """Test Moshier topocentric Sun position matches pyswisseph.

        The Sun's parallax is very small (~8.8 arcsec max), well within
        our Moshier C-vs-Python tolerance of 0.02° (~72 arcsec). This test
        may pass (xpass) even without topocentric implementation, since the
        Sun parallax is below the Moshier comparison tolerance.
        """
        jd = swe.julday(year, month, day, hour)
        flags_swe = swe.FLG_MOSEPH | swe.FLG_TOPOCTR | swe.FLG_SPEED
        flags_py = SEFLG_MOSEPH | SEFLG_TOPOCTR | SEFLG_SPEED

        swe.set_topo(TOPO_LON, TOPO_LAT, TOPO_ALT)
        ephem.swe_set_topo(TOPO_LON, TOPO_LAT, TOPO_ALT)

        pos_swe, _ = swe.calc_ut(jd, SE_SUN, flags_swe)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_SUN, flags_py)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        assert diff_lon < MOSHIER_LONGITUDE, (
            f"Moshier Topocentric Sun @ {date_desc}: "
            f"longitude diff {diff_lon:.6f}° exceeds tolerance {MOSHIER_LONGITUDE}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_py[0]:.6f}°)"
        )
        assert diff_lat < MOSHIER_LATITUDE, (
            f"Moshier Topocentric Sun @ {date_desc}: "
            f"latitude diff {diff_lat:.6f}° exceeds tolerance {MOSHIER_LATITUDE}° "
            f"(swe={pos_swe[1]:.6f}°, lib={pos_py[1]:.6f}°)"
        )
