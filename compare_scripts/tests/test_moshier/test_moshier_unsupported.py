"""
Moshier Ephemeris Unsupported Bodies Tests.

Validates that CalculationError is raised when attempting to calculate
positions for bodies not supported by the Moshier ephemeris:
- Chiron and other centaurs
- Main belt asteroids (Ceres, Pallas, Juno, Vesta)
- Trans-Neptunian Objects (TNOs)
- Hypothetical planets (Uranians, etc.)

These bodies require SPK kernel files and are not available in the
semi-analytical Moshier ephemeris.
"""

from __future__ import annotations

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_PLUTO,
    SE_CHIRON,
    SE_PHOLUS,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
    SE_AST_OFFSET,
    SE_FICT_OFFSET,
    SE_CUPIDO,
    SE_HADES,
    SE_ZEUS,
    SE_KRONOS,
    SE_APOLLON,
    SE_ADMETOS,
    SE_VULKANUS,
    SE_POSEIDON,
    SE_ISIS,
    SE_ERIS,
    SE_SEDNA,
    SE_INTP_APOG,
    SE_INTP_PERG,
    SEFLG_MOSEPH,
    SEFLG_SWIEPH,
)
from libephemeris.exceptions import CalculationError


# =============================================================================
# TEST CONFIGURATIONS
# =============================================================================

# Centaurs - not supported in Moshier
CENTAURS = [
    (SE_CHIRON, "Chiron"),
    (SE_PHOLUS, "Pholus"),
]

# Main belt asteroids - not supported in Moshier
MAIN_BELT_ASTEROIDS = [
    (SE_CERES, "Ceres"),
    (SE_PALLAS, "Pallas"),
    (SE_JUNO, "Juno"),
    (SE_VESTA, "Vesta"),
]

# Hypothetical Uranian planets - not supported in Moshier
URANIAN_PLANETS = [
    (SE_CUPIDO, "Cupido"),
    (SE_HADES, "Hades"),
    (SE_ZEUS, "Zeus"),
    (SE_KRONOS, "Kronos"),
    (SE_APOLLON, "Apollon"),
    (SE_ADMETOS, "Admetos"),
    (SE_VULKANUS, "Vulkanus"),
    (SE_POSEIDON, "Poseidon"),
]

# Other hypothetical bodies - not supported in Moshier
OTHER_HYPOTHETICAL = [
    (SE_ISIS, "Transpluto/Isis"),
]

# Trans-Neptunian Objects (sample) - not supported in Moshier
TNOS = [
    (SE_ERIS, "Eris"),
    (SE_SEDNA, "Sedna"),
]

# Interpolated lunar apsides - not supported in Moshier
INTERPOLATED_LUNAR = [
    (SE_INTP_APOG, "Interpolated Apogee"),
    (SE_INTP_PERG, "Interpolated Perigee"),
]

# Numbered asteroids via SE_AST_OFFSET - not supported in Moshier
NUMBERED_ASTEROIDS = [
    (SE_AST_OFFSET + 1, "1 Ceres (numbered)"),
    (SE_AST_OFFSET + 2, "2 Pallas (numbered)"),
    (SE_AST_OFFSET + 4, "4 Vesta (numbered)"),
    (SE_AST_OFFSET + 433, "433 Eros"),
    (SE_AST_OFFSET + 2060, "2060 Chiron (numbered)"),
]

# All unsupported bodies for comprehensive testing
ALL_UNSUPPORTED = (
    CENTAURS
    + MAIN_BELT_ASTEROIDS
    + URANIAN_PLANETS
    + OTHER_HYPOTHETICAL
    + INTERPOLATED_LUNAR
    + NUMBERED_ASTEROIDS
)

# Standard test date
STANDARD_JD = 2451545.0  # J2000.0


# =============================================================================
# TEST CLASSES
# =============================================================================


class TestCentaursRaiseCalculationError:
    """Verify centaurs raise CalculationError in Moshier mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", CENTAURS)
    def test_centaur_raises_error(self, body_id, body_name):
        """Centaurs should raise CalculationError with SEFLG_MOSEPH."""
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)

        err = exc_info.value
        # Error message should mention the body name
        assert body_name in str(err) or f"ID {body_id}" in str(err)
        # Should mention Moshier or lack of support
        assert "Moshier" in str(err) or "not available" in str(err)
        # Should suggest using SEFLG_SWIEPH
        assert "SEFLG_SWIEPH" in str(err)


class TestMainBeltAsteroidsRaiseCalculationError:
    """Verify main belt asteroids raise CalculationError in Moshier mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", MAIN_BELT_ASTEROIDS)
    def test_main_belt_asteroid_raises_error(self, body_id, body_name):
        """Main belt asteroids should raise CalculationError with SEFLG_MOSEPH."""
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)

        err = exc_info.value
        assert body_name in str(err) or f"ID {body_id}" in str(err)
        assert "Moshier" in str(err) or "not available" in str(err)


class TestUranianPlanetsRaiseCalculationError:
    """Verify hypothetical Uranian planets raise CalculationError in Moshier mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", URANIAN_PLANETS)
    def test_uranian_planet_raises_error(self, body_id, body_name):
        """Uranian planets should raise CalculationError with SEFLG_MOSEPH."""
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)

        err = exc_info.value
        # Should mention the body
        assert f"ID {body_id}" in str(err) or body_name in str(err)
        # Should mention Moshier limitation
        assert "Moshier" in str(err) or "not available" in str(err)


class TestTNOsRaiseCalculationError:
    """Verify Trans-Neptunian Objects raise CalculationError in Moshier mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", TNOS)
    def test_tno_raises_error(self, body_id, body_name):
        """TNOs should raise CalculationError with SEFLG_MOSEPH."""
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)

        err = exc_info.value
        assert "Moshier" in str(err) or "not available" in str(err)


class TestNumberedAsteroidsRaiseCalculationError:
    """Verify numbered asteroids (via SE_AST_OFFSET) raise CalculationError."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", NUMBERED_ASTEROIDS)
    def test_numbered_asteroid_raises_error(self, body_id, body_name):
        """Numbered asteroids should raise CalculationError with SEFLG_MOSEPH."""
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)

        err = exc_info.value
        assert "Moshier" in str(err) or "not available" in str(err)


class TestInterpolatedLunarApsidesRaiseCalculationError:
    """Verify interpolated lunar apsides raise CalculationError in Moshier mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", INTERPOLATED_LUNAR)
    def test_interpolated_apside_raises_error(self, body_id, body_name):
        """Interpolated lunar apsides should raise CalculationError with SEFLG_MOSEPH."""
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)

        err = exc_info.value
        assert "Moshier" in str(err) or "not available" in str(err)


# Note: TestUnsupportedBodiesWorkWithJPL was removed because it depends on
# having SPK kernels installed, which is not guaranteed in all environments.
# The main purpose of this test module is to verify that unsupported bodies
# raise CalculationError in Moshier mode, not to test JPL availability.


class TestErrorMessageQuality:
    """Verify error messages are helpful and informative."""

    @pytest.mark.comparison
    def test_chiron_error_message_format(self):
        """Chiron error message should be clear and helpful."""
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc_ut(STANDARD_JD, SE_CHIRON, SEFLG_MOSEPH)

        err_str = str(exc_info.value)

        # Should contain body name
        assert "Chiron" in err_str
        # Should contain body ID
        assert "15" in err_str or "ID 15" in err_str
        # Should suggest the fix
        assert "SEFLG_SWIEPH" in err_str
        # Should mention Moshier
        assert "Moshier" in err_str

    @pytest.mark.comparison
    def test_ceres_error_message_format(self):
        """Ceres error message should be clear and helpful."""
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc_ut(STANDARD_JD, SE_CERES, SEFLG_MOSEPH)

        err_str = str(exc_info.value)

        # Should contain body name
        assert "Ceres" in err_str
        # Should contain body ID
        assert "17" in err_str or "ID 17" in err_str


class TestSupportedBodiesDoNotRaiseError:
    """Verify that supported bodies do NOT raise CalculationError in Moshier mode."""

    SUPPORTED_BODIES = [
        (SE_SUN, "Sun"),
        (SE_MOON, "Moon"),
        (SE_MERCURY, "Mercury"),
        (SE_PLUTO, "Pluto"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", SUPPORTED_BODIES)
    def test_supported_body_does_not_raise(self, body_id, body_name):
        """Supported bodies should work without raising CalculationError."""
        # Should not raise
        pos, flag = ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)

        assert len(pos) == 6
        assert 0 <= pos[0] < 360


class TestSweCalcUnsupportedBodies:
    """Test unsupported bodies with swe_calc (TT time) in addition to swe_calc_ut."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", CENTAURS)
    def test_swe_calc_centaur_raises_error(self, body_id, body_name):
        """swe_calc should also raise CalculationError for unsupported bodies."""
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc(STANDARD_JD, body_id, SEFLG_MOSEPH)

        err = exc_info.value
        assert "Moshier" in str(err) or "not available" in str(err)


class TestComprehensiveUnsupportedBodyCoverage:
    """Comprehensive test covering all unsupported bodies."""

    @pytest.mark.comparison
    @pytest.mark.slow
    def test_all_unsupported_bodies_summary(self):
        """Summary test for all unsupported bodies."""
        passed = []
        failed = []

        for body_id, body_name in ALL_UNSUPPORTED:
            try:
                ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)
                # If we get here, it didn't raise - that's a failure
                failed.append(
                    f"{body_name} (ID {body_id}): did not raise CalculationError"
                )
            except CalculationError:
                # This is expected
                passed.append(f"{body_name} (ID {body_id})")
            except Exception as e:
                # Unexpected exception type
                failed.append(
                    f"{body_name} (ID {body_id}): raised {type(e).__name__}: {e}"
                )

        total = len(ALL_UNSUPPORTED)
        pass_count = len(passed)

        assert pass_count == total, (
            f"Unsupported body test: {pass_count}/{total} correctly raised CalculationError.\n"
            f"Failures:\n  " + "\n  ".join(failed)
        )
