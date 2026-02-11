"""
Moshier Unsupported Bodies Cross-Library Error Comparison Tests.

Compares the error behavior for unsupported bodies in Moshier mode between
pyswisseph (C library) and libephemeris (Python reimplementation).

While test_moshier_unsupported.py verifies that libephemeris raises
CalculationError for unsupported bodies, this module verifies that pyswisseph
ALSO rejects the same bodies, and documents the exact error type and message
from each library for cross-library migration compatibility.

This is critical for applications migrating from pyswisseph to libephemeris:
code with try/except for swisseph.Error needs to know that libephemeris raises
a compatible error (CalculationError inherits from Error) for the same bodies.

Body categories tested (8 key body IDs):
- Centaur: SE_CHIRON (15)
- Main belt asteroids: SE_CERES (17), SE_PALLAS (18), SE_JUNO (19), SE_VESTA (20)
- Hypothetical Uranian: SE_CUPIDO (40)
- TNO: SE_ERIS (146199)
- Numbered asteroid: SE_AST_OFFSET + 433 (10433, Eros)

KNOWN BEHAVIORAL DIVERGENCE:
    SE_CUPIDO (40) - pyswisseph ACCEPTS this body in Moshier mode because
    hypothetical Uranian planets use analytical Keplerian orbital elements
    (not SPK kernels). libephemeris rejects it. This is a documented
    incompatibility for migration.

Reference: libephemeris/planets.py:830-839 (Moshier unsupported body check)
"""

from __future__ import annotations

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_CHIRON,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
    SE_CUPIDO,
    SE_ERIS,
    SE_AST_OFFSET,
    SEFLG_MOSEPH,
)
from libephemeris.exceptions import CalculationError, Error


# =============================================================================
# TEST CONFIGURATIONS
# =============================================================================

# Standard test date
STANDARD_JD = 2451545.0  # J2000.0

# Bodies rejected by BOTH pyswisseph and libephemeris in Moshier mode.
# These require SPK kernels that are not available in Moshier.
BODIES_REJECTED_BY_BOTH = [
    (SE_CHIRON, "Chiron", "centaur"),
    (SE_CERES, "Ceres", "main_belt_asteroid"),
    (SE_PALLAS, "Pallas", "main_belt_asteroid"),
    (SE_JUNO, "Juno", "main_belt_asteroid"),
    (SE_VESTA, "Vesta", "main_belt_asteroid"),
    (SE_ERIS, "Eris", "tno"),
    (SE_AST_OFFSET + 433, "433 Eros", "numbered_asteroid"),
]

# Bodies with DIVERGENT behavior between the two libraries.
# pyswisseph ACCEPTS these in Moshier mode (analytical Keplerian elements),
# but libephemeris REJECTS them.
BODIES_DIVERGENT = [
    (SE_CUPIDO, "Cupido", "hypothetical_uranian"),
]

# All 8 key body IDs for comprehensive coverage
ALL_BODIES = BODIES_REJECTED_BY_BOTH + BODIES_DIVERGENT


def _call_pyswisseph(jd: float, body_id: int, flags: int) -> dict:
    """Call pyswisseph calc_ut and capture result or error.

    Returns:
        Dict with keys:
            - 'raised': bool - whether an exception was raised
            - 'exc_type': str - exception class name (if raised)
            - 'exc_message': str - exception message (if raised)
            - 'result': tuple - (position, flag) if no exception
    """
    try:
        result = swe.calc_ut(jd, body_id, flags)
        return {
            "raised": False,
            "exc_type": None,
            "exc_message": None,
            "result": result,
        }
    except Exception as e:
        return {
            "raised": True,
            "exc_type": type(e).__name__,
            "exc_message": str(e),
            "result": None,
        }


def _call_libephemeris(jd: float, body_id: int, flags: int) -> dict:
    """Call libephemeris swe_calc_ut and capture result or error.

    Returns:
        Dict with keys:
            - 'raised': bool - whether an exception was raised
            - 'exc_type': str - exception class name (if raised)
            - 'exc_message': str - exception message (if raised)
            - 'result': tuple - (position, flag) if no exception
    """
    try:
        result = ephem.swe_calc_ut(jd, body_id, flags)
        return {
            "raised": False,
            "exc_type": None,
            "exc_message": None,
            "result": result,
        }
    except Exception as e:
        return {
            "raised": True,
            "exc_type": type(e).__name__,
            "exc_message": str(e),
            "result": None,
        }


# =============================================================================
# TEST CLASS: Bodies rejected by both libraries
# =============================================================================


class TestMoshierUnsupportedBothReject:
    """Compare error behavior for bodies rejected by both libraries.

    For each of 7 body IDs (centaurs, asteroids, TNOs), verifies that:
    1. Both pyswisseph and libephemeris reject the body in Moshier mode
    2. Both raise an exception (not a silent partial result)
    3. The libephemeris error is catchable as Error (pyswisseph-compatible)
    4. Documents the exact exception type and message from each library
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        BODIES_REJECTED_BY_BOTH,
        ids=[f"{name}_ID{bid}" for bid, name, _ in BODIES_REJECTED_BY_BOTH],
    )
    def test_both_libraries_reject_body(self, body_id, body_name, category):
        """Both libraries should reject unsupported body in Moshier mode.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (centaur, main_belt_asteroid, etc.).
        """
        swe_result = _call_pyswisseph(STANDARD_JD, body_id, SEFLG_MOSEPH)
        lib_result = _call_libephemeris(STANDARD_JD, body_id, SEFLG_MOSEPH)

        # Both must raise an exception
        assert swe_result["raised"], (
            f"pyswisseph did NOT raise an error for {body_name} (ID {body_id}) "
            f"with SEFLG_MOSEPH. Got result: {swe_result['result']}"
        )
        assert lib_result["raised"], (
            f"libephemeris did NOT raise an error for {body_name} (ID {body_id}) "
            f"with SEFLG_MOSEPH. Got result: {lib_result['result']}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        BODIES_REJECTED_BY_BOTH,
        ids=[f"{name}_ID{bid}" for bid, name, _ in BODIES_REJECTED_BY_BOTH],
    )
    def test_pyswisseph_error_type(self, body_id, body_name, category):
        """pyswisseph should raise an exception catchable via except Exception.

        Documents the exact exception type pyswisseph raises for each body.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (centaur, main_belt_asteroid, etc.).
        """
        swe_result = _call_pyswisseph(STANDARD_JD, body_id, SEFLG_MOSEPH)

        assert swe_result["raised"], (
            f"pyswisseph did NOT raise for {body_name} (ID {body_id})"
        )
        # Document: pyswisseph typically raises swisseph.Error
        assert swe_result["exc_type"] in ("Error", "SwissephError"), (
            f"pyswisseph raised unexpected {swe_result['exc_type']} "
            f"for {body_name} (ID {body_id}): {swe_result['exc_message']}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        BODIES_REJECTED_BY_BOTH,
        ids=[f"{name}_ID{bid}" for bid, name, _ in BODIES_REJECTED_BY_BOTH],
    )
    def test_pyswisseph_error_is_not_silent(self, body_id, body_name, category):
        """pyswisseph must not return a partial result with no error.

        Some C libraries return partial results (e.g., zeros) instead of
        raising errors. This verifies pyswisseph properly raises.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (centaur, main_belt_asteroid, etc.).
        """
        swe_result = _call_pyswisseph(STANDARD_JD, body_id, SEFLG_MOSEPH)

        assert swe_result["raised"], (
            f"pyswisseph returned a result for {body_name} (ID {body_id}) "
            f"with SEFLG_MOSEPH instead of raising an error. "
            f"Result: {swe_result['result']}. "
            f"This means pyswisseph silently accepts this body in Moshier mode, "
            f"which creates a behavioral incompatibility with libephemeris."
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        BODIES_REJECTED_BY_BOTH,
        ids=[f"{name}_ID{bid}" for bid, name, _ in BODIES_REJECTED_BY_BOTH],
    )
    def test_cross_library_error_documentation(self, body_id, body_name, category):
        """Document and compare error details from both libraries.

        This test captures and compares the full error profile from both
        libraries for each body, ensuring consistent rejection behavior.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (centaur, main_belt_asteroid, etc.).
        """
        swe_result = _call_pyswisseph(STANDARD_JD, body_id, SEFLG_MOSEPH)
        lib_result = _call_libephemeris(STANDARD_JD, body_id, SEFLG_MOSEPH)

        # Both must reject the body
        assert swe_result["raised"] and lib_result["raised"], (
            f"Rejection mismatch for {body_name} (ID {body_id}, {category}): "
            f"pyswisseph raised={swe_result['raised']}, "
            f"libephemeris raised={lib_result['raised']}"
        )

        # Document the error profile for each library
        # pyswisseph: typically raises swisseph.Error with C-level message
        # libephemeris: raises CalculationError (subclass of Error) with
        # Python-level message mentioning body name, ID, and Moshier
        #
        # The messages will differ (C vs Python), but both reject the body.
        # This is the expected and correct behavior for migration.
        assert swe_result["exc_type"] is not None
        assert lib_result["exc_type"] == "CalculationError"


# =============================================================================
# TEST CLASS: libephemeris behavior for ALL 8 bodies
# =============================================================================


class TestMoshierUnsupportedLibephemeris:
    """Verify libephemeris rejects all 8 bodies and errors are well-formed.

    Tests that apply to all 8 bodies regardless of pyswisseph behavior:
    - libephemeris raises CalculationError
    - Error is catchable as Error (base class)
    - Error message mentions the body ID/name
    - Error message mentions Moshier
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        ALL_BODIES,
        ids=[f"{name}_ID{bid}" for bid, name, _ in ALL_BODIES],
    )
    def test_libephemeris_raises_calculation_error(self, body_id, body_name, category):
        """libephemeris should raise CalculationError for unsupported bodies.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (centaur, main_belt_asteroid, etc.).
        """
        with pytest.raises(CalculationError) as exc_info:
            ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)

        err = exc_info.value
        # CalculationError must be catchable as Error (pyswisseph compat)
        assert isinstance(err, Error), (
            f"CalculationError for {body_name} is not an instance of Error. "
            f"Migration code using 'except Error' would miss this."
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        ALL_BODIES,
        ids=[f"{name}_ID{bid}" for bid, name, _ in ALL_BODIES],
    )
    def test_error_messages_mention_body(self, body_id, body_name, category):
        """Error messages from libephemeris should reference the body.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (centaur, main_belt_asteroid, etc.).
        """
        lib_result = _call_libephemeris(STANDARD_JD, body_id, SEFLG_MOSEPH)

        assert lib_result["raised"]
        # libephemeris error must mention body ID or name
        msg = lib_result["exc_message"]
        assert str(body_id) in msg or body_name in msg, (
            f"libephemeris error message does not mention body {body_name} "
            f"(ID {body_id}). Message: {msg}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        ALL_BODIES,
        ids=[f"{name}_ID{bid}" for bid, name, _ in ALL_BODIES],
    )
    def test_libephemeris_error_mentions_moshier(self, body_id, body_name, category):
        """libephemeris error should mention Moshier to explain the limitation.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (centaur, main_belt_asteroid, etc.).
        """
        lib_result = _call_libephemeris(STANDARD_JD, body_id, SEFLG_MOSEPH)

        assert lib_result["raised"]
        msg = lib_result["exc_message"]
        assert "Moshier" in msg or "not available" in msg, (
            f"libephemeris error for {body_name} does not mention 'Moshier' or "
            f"'not available'. Message: {msg}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        ALL_BODIES,
        ids=[f"{name}_ID{bid}" for bid, name, _ in ALL_BODIES],
    )
    def test_libephemeris_error_catchable_as_base_error(
        self, body_id, body_name, category
    ):
        """libephemeris error must be catchable by 'except Error'.

        This is the critical migration compatibility check: applications that
        use 'except swe.Error' must also catch libephemeris errors when
        switching to 'except ephem.Error'.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (centaur, main_belt_asteroid, etc.).
        """
        caught = False
        try:
            ephem.swe_calc_ut(STANDARD_JD, body_id, SEFLG_MOSEPH)
        except Error:
            caught = True
        except Exception as e:
            pytest.fail(
                f"libephemeris raised {type(e).__name__} for {body_name} "
                f"(ID {body_id}) which is NOT catchable by 'except Error'. "
                f"Message: {e}"
            )

        assert caught, (
            f"libephemeris did not raise any error for {body_name} (ID {body_id})"
        )


# =============================================================================
# TEST CLASS: Behavioral divergence documentation (SE_CUPIDO)
# =============================================================================


class TestMoshierUnsupportedDivergence:
    """Document behavioral divergences between pyswisseph and libephemeris.

    SE_CUPIDO (40) and other hypothetical Uranian planets are ACCEPTED by
    pyswisseph in Moshier mode because the C Swiss Ephemeris implements them
    using analytical Keplerian orbital elements (not SPK kernels).

    libephemeris rejects these bodies because _MOSHIER_SUPPORTED_BODIES does
    not include hypothetical planets.

    This is a known incompatibility: applications using pyswisseph to compute
    hypothetical Uranian positions in Moshier mode will get a CalculationError
    when switching to libephemeris.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        BODIES_DIVERGENT,
        ids=[f"{name}_ID{bid}" for bid, name, _ in BODIES_DIVERGENT],
    )
    def test_pyswisseph_accepts_body(self, body_id, body_name, category):
        """pyswisseph accepts this body in Moshier mode (analytical elements).

        Documents that pyswisseph returns a valid result (not an error) for
        hypothetical Uranian planets in Moshier mode, using Keplerian orbital
        elements built into the C library.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (hypothetical_uranian).
        """
        swe_result = _call_pyswisseph(STANDARD_JD, body_id, SEFLG_MOSEPH)

        assert not swe_result["raised"], (
            f"pyswisseph unexpectedly raised for {body_name} (ID {body_id}) "
            f"with SEFLG_MOSEPH: {swe_result['exc_type']}: "
            f"{swe_result['exc_message']}. "
            f"Expected: pyswisseph accepts hypothetical Uranians in Moshier mode."
        )

        # Verify it returned a valid-looking result
        pos, flag = swe_result["result"]
        assert len(pos) == 6, (
            f"pyswisseph returned {len(pos)} values for {body_name}, expected 6"
        )
        assert 0 <= pos[0] < 360, (
            f"pyswisseph returned invalid longitude {pos[0]} for {body_name}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        BODIES_DIVERGENT,
        ids=[f"{name}_ID{bid}" for bid, name, _ in BODIES_DIVERGENT],
    )
    def test_libephemeris_rejects_body(self, body_id, body_name, category):
        """libephemeris rejects this body in Moshier mode.

        Documents that libephemeris raises CalculationError for hypothetical
        Uranian planets in Moshier mode, because _MOSHIER_SUPPORTED_BODIES
        does not include them.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (hypothetical_uranian).
        """
        lib_result = _call_libephemeris(STANDARD_JD, body_id, SEFLG_MOSEPH)

        assert lib_result["raised"], (
            f"libephemeris did NOT raise for {body_name} (ID {body_id}) "
            f"with SEFLG_MOSEPH. Got result: {lib_result['result']}"
        )
        assert lib_result["exc_type"] == "CalculationError", (
            f"libephemeris raised {lib_result['exc_type']} instead of "
            f"CalculationError for {body_name}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name,category",
        BODIES_DIVERGENT,
        ids=[f"{name}_ID{bid}" for bid, name, _ in BODIES_DIVERGENT],
    )
    def test_divergence_documented(self, body_id, body_name, category):
        """Document the full divergence profile for migration awareness.

        Captures both libraries' responses to highlight the incompatibility:
        - pyswisseph: returns a valid result (Keplerian analytical elements)
        - libephemeris: raises CalculationError

        Applications migrating from pyswisseph that use hypothetical Uranian
        planets in Moshier mode must add try/except handling for this case.

        Args:
            body_id: Swiss Ephemeris body constant.
            body_name: Human-readable body name.
            category: Body category (hypothetical_uranian).
        """
        swe_result = _call_pyswisseph(STANDARD_JD, body_id, SEFLG_MOSEPH)
        lib_result = _call_libephemeris(STANDARD_JD, body_id, SEFLG_MOSEPH)

        # Confirm the divergence: pyswisseph accepts, libephemeris rejects
        assert not swe_result["raised"], (
            f"Divergence changed: pyswisseph now rejects {body_name}"
        )
        assert lib_result["raised"], (
            f"Divergence changed: libephemeris now accepts {body_name}"
        )

        # Document: pyswisseph returns valid position data
        pos, _ = swe_result["result"]
        assert 0 <= pos[0] < 360

        # Document: libephemeris raises CalculationError with Moshier mention
        assert lib_result["exc_type"] == "CalculationError"
        assert (
            "Moshier" in lib_result["exc_message"]
            or "not available" in lib_result["exc_message"]
        )
