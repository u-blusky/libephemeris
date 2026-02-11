"""
Moshier Return Flag Cross-Library Comparison Tests.

Compares the return flag (second element of the tuple) from swe.calc_ut() (C
pyswisseph) and ephem.swe_calc_ut() (Python libephemeris) with SEFLG_MOSEPH
for 8 flag combinations x 2 planets = 16 test cases.

PROBLEM: libephemeris returns iflag unchanged as return flag (planets.py:1087);
the C Swiss Ephemeris modifies the return flag to indicate which flags were
actually applied:
  - SEFLG_SIDEREAL: C adds SEFLG_NONUT (64) because nutation is bypassed
  - SEFLG_HELCTR: C adds SEFLG_NOGDEFL (512) + SEFLG_NOABERR (1024) because
    aberration and gravitational deflection are geocentric corrections
  - SEFLG_J2000: C adds SEFLG_NONUT (64) because J2000 frame has no nutation

IMPACT: The return flag is the only runtime mechanism to determine which
ephemeris was actually used; applications implementing automatic fallback
(try SEFLG_SWIEPH, if fails try SEFLG_MOSEPH, check retflag) depend on
return flag correctness. Discrepancies can cause infinite loops or incorrect
fallback behavior.

Each test documents the exact bit-level discrepancy between C and Python.
Tests where libephemeris returns a different retflag than pyswisseph are
marked xfail(strict=True) to document the known incompatibility.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_VENUS,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
    SEFLG_SIDEREAL,
    SEFLG_HELCTR,
    SEFLG_TOPOCTR,
    SEFLG_ICRS,
    SEFLG_J2000,
    SEFLG_NONUT,
    SEFLG_NOABERR,
    SEFLG_NOGDEFL,
)


# ============================================================================
# CONSTANTS
# ============================================================================

# Two representative planets for each flag combination:
# - Sun (SE_SUN=0): the default body, special in heliocentric mode
# - Venus (SE_VENUS=3): a standard planet with no special-case handling
PLANETS = [
    (SE_SUN, "Sun"),
    (SE_VENUS, "Venus"),
]

# Julian Day at J2000.0 epoch
JD_J2000 = 2451545.0

# Flag combination descriptors
# Each entry: (label, extra_flags, needs_setup, expected_match)
# needs_setup: "sidereal" | "topocentric" | None
# expected_match: True if C and Python retflags should be identical
FLAG_COMBOS = [
    ("MOSEPH", 0, None, True),
    ("MOSEPH|SPEED", SEFLG_SPEED, None, True),
    ("MOSEPH|EQUATORIAL", SEFLG_EQUATORIAL, None, True),
    ("MOSEPH|SIDEREAL", SEFLG_SIDEREAL, "sidereal", False),
    ("MOSEPH|HELCTR", SEFLG_HELCTR, None, False),
    ("MOSEPH|TOPOCTR", SEFLG_TOPOCTR, "topocentric", True),
    ("MOSEPH|ICRS", SEFLG_ICRS, None, True),
    ("MOSEPH|J2000", SEFLG_J2000, None, False),
]


# ============================================================================
# HELPERS
# ============================================================================


def _flag_bits_description(flag: int) -> str:
    """Return a human-readable description of set flag bits."""
    known_bits = [
        (0x00001, "SEFLG_JPLEPH"),
        (0x00002, "SEFLG_SWIEPH"),
        (0x00004, "SEFLG_MOSEPH"),
        (0x00008, "SEFLG_HELCTR"),
        (0x00010, "SEFLG_TRUEPOS"),
        (0x00020, "SEFLG_J2000"),
        (0x00040, "SEFLG_NONUT"),
        (0x00080, "SEFLG_SPEED3"),
        (0x00100, "SEFLG_SPEED"),
        (0x00200, "SEFLG_NOGDEFL"),
        (0x00400, "SEFLG_NOABERR"),
        (0x00800, "SEFLG_EQUATORIAL"),
        (0x01000, "SEFLG_XYZ"),
        (0x02000, "SEFLG_RADIANS"),
        (0x04000, "SEFLG_BARYCTR"),
        (0x08000, "SEFLG_TOPOCTR"),
        (0x10000, "SEFLG_SIDEREAL"),
        (0x20000, "SEFLG_ICRS"),
    ]
    parts = []
    for bit, name in known_bits:
        if flag & bit:
            parts.append(name)
    return " | ".join(parts) if parts else "0"


def _xor_bits_description(xor: int) -> str:
    """Describe which flag bits differ (XOR between C and Python retflags)."""
    if xor == 0:
        return "no differences"
    return _flag_bits_description(xor)


def _setup_for_combo(needs_setup: str | None) -> None:
    """Configure sidereal mode or topocentric position if needed."""
    if needs_setup == "sidereal":
        swe.set_sid_mode(0)  # Fagan/Bradley
        ephem.swe_set_sid_mode(0)
    elif needs_setup == "topocentric":
        swe.set_topo(12.4964, 41.9028, 0)  # Rome
        ephem.swe_set_topo(12.4964, 41.9028, 0)


# ============================================================================
# TEST CLASS
# ============================================================================


class TestMoshierReturnFlag:
    """Compare return flag (retflag) between pyswisseph C and libephemeris Python.

    The C Swiss Ephemeris modifies the return flag to communicate which
    calculation adjustments were actually applied. Libephemeris currently
    returns iflag unchanged (planets.py:1087), which creates bit-level
    discrepancies for SIDEREAL, HELCTR, and J2000 modes.

    Each test case compares retflag bit-by-bit and documents the exact
    discrepancy. Cases where a mismatch is expected are marked as
    xfail(strict=True) to ensure the test suite passes while documenting
    the known incompatibility.
    """

    # ------------------------------------------------------------------
    # MOSEPH only (no extra flags) — retflag should match
    # ------------------------------------------------------------------

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_retflag_moseph_only(self, planet_id, planet_name):
        """MOSEPH only: retflag should match between C and Python.

        The C library returns SEFLG_MOSEPH (4) unchanged when no other
        modifier flags are set. Libephemeris also returns iflag=4 unchanged.
        """
        flag = SEFLG_MOSEPH
        _, retflag_c = swe.calc_ut(JD_J2000, planet_id, flag)
        _, retflag_py = ephem.swe_calc_ut(JD_J2000, planet_id, flag)

        assert retflag_c == retflag_py, (
            f"{planet_name} MOSEPH retflag mismatch: "
            f"C={retflag_c} ({_flag_bits_description(retflag_c)}), "
            f"Py={retflag_py} ({_flag_bits_description(retflag_py)}), "
            f"XOR={retflag_c ^ retflag_py} ({_xor_bits_description(retflag_c ^ retflag_py)})"
        )

    # ------------------------------------------------------------------
    # MOSEPH|SPEED — retflag should match
    # ------------------------------------------------------------------

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_retflag_moseph_speed(self, planet_id, planet_name):
        """MOSEPH|SPEED: retflag should match between C and Python.

        SEFLG_SPEED (256) is always applied; the C library does not modify
        the return flag for this combination.
        """
        flag = SEFLG_MOSEPH | SEFLG_SPEED
        _, retflag_c = swe.calc_ut(JD_J2000, planet_id, flag)
        _, retflag_py = ephem.swe_calc_ut(JD_J2000, planet_id, flag)

        assert retflag_c == retflag_py, (
            f"{planet_name} MOSEPH|SPEED retflag mismatch: "
            f"C={retflag_c} ({_flag_bits_description(retflag_c)}), "
            f"Py={retflag_py} ({_flag_bits_description(retflag_py)}), "
            f"XOR={retflag_c ^ retflag_py} ({_xor_bits_description(retflag_c ^ retflag_py)})"
        )

    # ------------------------------------------------------------------
    # MOSEPH|EQUATORIAL — retflag should match
    # ------------------------------------------------------------------

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_retflag_moseph_equatorial(self, planet_id, planet_name):
        """MOSEPH|EQUATORIAL: retflag should match between C and Python.

        SEFLG_EQUATORIAL (2048) triggers coordinate transformation from
        ecliptic to equatorial (RA/Dec). The C library does not modify the
        return flag for this transformation.
        """
        flag = SEFLG_MOSEPH | SEFLG_EQUATORIAL
        _, retflag_c = swe.calc_ut(JD_J2000, planet_id, flag)
        _, retflag_py = ephem.swe_calc_ut(JD_J2000, planet_id, flag)

        assert retflag_c == retflag_py, (
            f"{planet_name} MOSEPH|EQUATORIAL retflag mismatch: "
            f"C={retflag_c} ({_flag_bits_description(retflag_c)}), "
            f"Py={retflag_py} ({_flag_bits_description(retflag_py)}), "
            f"XOR={retflag_c ^ retflag_py} ({_xor_bits_description(retflag_c ^ retflag_py)})"
        )

    # ------------------------------------------------------------------
    # MOSEPH|SIDEREAL — retflag MISMATCH: C adds SEFLG_NONUT (64)
    # ------------------------------------------------------------------

    @pytest.mark.comparison
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "libephemeris returns iflag unchanged (planets.py:1087); "
            "C library adds SEFLG_NONUT (64) to retflag in sidereal mode "
            "because nutation correction is bypassed"
        ),
    )
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_retflag_moseph_sidereal(self, planet_id, planet_name):
        """MOSEPH|SIDEREAL: retflag MISMATCH — C adds SEFLG_NONUT (64).

        In sidereal mode, the C library sets SEFLG_NONUT in the return flag
        to indicate that nutation correction was not applied (sidereal
        coordinates bypass the nutation step). Libephemeris returns iflag
        unchanged.

        Expected discrepancy:
            C retflag:  MOSEPH | SIDEREAL | NONUT  = 65604
            Py retflag: MOSEPH | SIDEREAL          = 65540
            XOR: SEFLG_NONUT (64)
        """
        _setup_for_combo("sidereal")
        flag = SEFLG_MOSEPH | SEFLG_SIDEREAL
        _, retflag_c = swe.calc_ut(JD_J2000, planet_id, flag)
        _, retflag_py = ephem.swe_calc_ut(JD_J2000, planet_id, flag)

        assert retflag_c == retflag_py, (
            f"{planet_name} MOSEPH|SIDEREAL retflag mismatch: "
            f"C={retflag_c} ({_flag_bits_description(retflag_c)}), "
            f"Py={retflag_py} ({_flag_bits_description(retflag_py)}), "
            f"XOR={retflag_c ^ retflag_py} ({_xor_bits_description(retflag_c ^ retflag_py)})"
        )

    # ------------------------------------------------------------------
    # MOSEPH|HELCTR — retflag MISMATCH: C adds NOGDEFL (512) + NOABERR (1024)
    # ------------------------------------------------------------------

    @pytest.mark.comparison
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "libephemeris returns iflag unchanged (planets.py:1087); "
            "C library adds SEFLG_NOGDEFL (512) + SEFLG_NOABERR (1024) to "
            "retflag in heliocentric mode because aberration and gravitational "
            "deflection are geocentric-only corrections"
        ),
    )
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_retflag_moseph_helctr(self, planet_id, planet_name):
        """MOSEPH|HELCTR: retflag MISMATCH — C adds NOGDEFL + NOABERR.

        In heliocentric mode, the C library sets SEFLG_NOGDEFL (512) and
        SEFLG_NOABERR (1024) in the return flag because aberration and
        gravitational light deflection are geocentric corrections that do
        not apply to heliocentric positions. Libephemeris returns iflag
        unchanged.

        Expected discrepancy:
            C retflag:  MOSEPH | HELCTR | NOGDEFL | NOABERR  = 1548
            Py retflag: MOSEPH | HELCTR                      = 12
            XOR: SEFLG_NOGDEFL (512) | SEFLG_NOABERR (1024) = 1536
        """
        flag = SEFLG_MOSEPH | SEFLG_HELCTR
        # Use Venus for heliocentric; Sun heliocentric is a degenerate case
        # but both planets exhibit the same retflag discrepancy
        _, retflag_c = swe.calc_ut(JD_J2000, planet_id, flag)
        _, retflag_py = ephem.swe_calc_ut(JD_J2000, planet_id, flag)

        assert retflag_c == retflag_py, (
            f"{planet_name} MOSEPH|HELCTR retflag mismatch: "
            f"C={retflag_c} ({_flag_bits_description(retflag_c)}), "
            f"Py={retflag_py} ({_flag_bits_description(retflag_py)}), "
            f"XOR={retflag_c ^ retflag_py} ({_xor_bits_description(retflag_c ^ retflag_py)})"
        )

    # ------------------------------------------------------------------
    # MOSEPH|TOPOCTR — retflag should match
    # ------------------------------------------------------------------

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_retflag_moseph_topoctr(self, planet_id, planet_name):
        """MOSEPH|TOPOCTR: retflag should match between C and Python.

        SEFLG_TOPOCTR (32768) applies topocentric correction. The C library
        does not remove or modify this flag in the return value when
        Moshier mode is used (unlike some documentation that suggests
        TOPOCTR may be stripped in certain modes).
        """
        _setup_for_combo("topocentric")
        flag = SEFLG_MOSEPH | SEFLG_TOPOCTR
        _, retflag_c = swe.calc_ut(JD_J2000, planet_id, flag)
        _, retflag_py = ephem.swe_calc_ut(JD_J2000, planet_id, flag)

        assert retflag_c == retflag_py, (
            f"{planet_name} MOSEPH|TOPOCTR retflag mismatch: "
            f"C={retflag_c} ({_flag_bits_description(retflag_c)}), "
            f"Py={retflag_py} ({_flag_bits_description(retflag_py)}), "
            f"XOR={retflag_c ^ retflag_py} ({_xor_bits_description(retflag_c ^ retflag_py)})"
        )

    # ------------------------------------------------------------------
    # MOSEPH|ICRS — retflag should match
    # ------------------------------------------------------------------

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_retflag_moseph_icrs(self, planet_id, planet_name):
        """MOSEPH|ICRS: retflag should match between C and Python.

        SEFLG_ICRS (131072) selects the ICRS reference frame. The C library
        returns the flag unchanged in the return value.
        """
        flag = SEFLG_MOSEPH | SEFLG_ICRS
        _, retflag_c = swe.calc_ut(JD_J2000, planet_id, flag)
        _, retflag_py = ephem.swe_calc_ut(JD_J2000, planet_id, flag)

        assert retflag_c == retflag_py, (
            f"{planet_name} MOSEPH|ICRS retflag mismatch: "
            f"C={retflag_c} ({_flag_bits_description(retflag_c)}), "
            f"Py={retflag_py} ({_flag_bits_description(retflag_py)}), "
            f"XOR={retflag_c ^ retflag_py} ({_xor_bits_description(retflag_c ^ retflag_py)})"
        )

    # ------------------------------------------------------------------
    # MOSEPH|J2000 — retflag MISMATCH: C adds SEFLG_NONUT (64)
    # ------------------------------------------------------------------

    @pytest.mark.comparison
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "libephemeris returns iflag unchanged (planets.py:1087); "
            "C library adds SEFLG_NONUT (64) to retflag in J2000 mode "
            "because the J2000 reference frame does not apply nutation"
        ),
    )
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_retflag_moseph_j2000(self, planet_id, planet_name):
        """MOSEPH|J2000: retflag MISMATCH — C adds SEFLG_NONUT (64).

        In J2000 mode, the C library sets SEFLG_NONUT in the return flag
        because the J2000 ecliptic frame does not include nutation
        correction. Libephemeris returns iflag unchanged.

        Expected discrepancy:
            C retflag:  MOSEPH | J2000 | NONUT  = 100
            Py retflag: MOSEPH | J2000          = 36
            XOR: SEFLG_NONUT (64)
        """
        flag = SEFLG_MOSEPH | SEFLG_J2000
        _, retflag_c = swe.calc_ut(JD_J2000, planet_id, flag)
        _, retflag_py = ephem.swe_calc_ut(JD_J2000, planet_id, flag)

        assert retflag_c == retflag_py, (
            f"{planet_name} MOSEPH|J2000 retflag mismatch: "
            f"C={retflag_c} ({_flag_bits_description(retflag_c)}), "
            f"Py={retflag_py} ({_flag_bits_description(retflag_py)}), "
            f"XOR={retflag_c ^ retflag_py} ({_xor_bits_description(retflag_c ^ retflag_py)})"
        )
