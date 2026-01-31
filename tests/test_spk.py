"""
Tests for SPK kernel support in libephemeris.

These tests verify:
- SPK body registration and unregistration
- NAIF ID deduction
- Fallback to Keplerian when SPK not available
- Coverage checking
- Error handling (SPKNotFoundError)
- jplephem version compatibility

Note: Integration tests that download actual SPK files from Horizons
are skipped by default. Set LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 to run them.
"""

import os
import pytest
from unittest.mock import patch, MagicMock
from packaging import version

import jplephem

import libephemeris as eph
from libephemeris import spk, state
from libephemeris.constants import (
    SE_CHIRON,
    SE_CERES,
    SE_ERIS,
    NAIF_CHIRON,
    NAIF_CERES,
    NAIF_ERIS,
    NAIF_ASTEROID_OFFSET,
    SEFLG_SPEED,
)
from libephemeris.exceptions import SPKNotFoundError


class TestJplephemVersion:
    """Test jplephem version requirements for optimal SPK reading performance."""

    # Minimum required version for all features (from skyfield dependency)
    MINIMUM_VERSION = "2.13"
    # Recommended version with NumPy compatibility fixes
    RECOMMENDED_VERSION = "2.24"

    def test_jplephem_installed(self):
        """Verify jplephem is installed."""
        assert jplephem is not None
        assert hasattr(jplephem, "__version__")

    def test_jplephem_minimum_version(self):
        """Verify jplephem meets minimum version requirement.

        Version 2.13+ is required by skyfield for:
        - OutOfRangeError with array attribute for date errors
        - Essential SPK reading capabilities
        """
        installed = version.parse(jplephem.__version__)
        minimum = version.parse(self.MINIMUM_VERSION)
        assert installed >= minimum, (
            f"jplephem {jplephem.__version__} is below minimum version "
            f"{self.MINIMUM_VERSION} required for SPK support"
        )

    def test_jplephem_recommended_version(self):
        """Verify jplephem meets recommended version for optimal performance.

        Version 2.24+ includes:
        - NumPy deprecation warning fixes (.shape -> .reshape())
        - Better compatibility with modern NumPy versions
        """
        installed = version.parse(jplephem.__version__)
        recommended = version.parse(self.RECOMMENDED_VERSION)
        assert installed >= recommended, (
            f"jplephem {jplephem.__version__} is below recommended version "
            f"{self.RECOMMENDED_VERSION}. Consider upgrading for NumPy compatibility."
        )

    def test_jplephem_spk_module_available(self):
        """Verify jplephem.spk module is available for SPK reading."""
        from jplephem.spk import SPK

        assert SPK is not None
        # Verify SPK has essential methods
        assert hasattr(SPK, "open")

    def test_jplephem_daf_module_available(self):
        """Verify jplephem.daf module is available for DAF file reading."""
        from jplephem.daf import DAF

        assert DAF is not None


class TestNaifIdDeduction:
    """Test NAIF ID deduction from body identifiers."""

    def test_extract_number_from_digits(self):
        """Extract asteroid number from pure digit string."""
        assert spk._extract_asteroid_number("2060") == 2060
        assert spk._extract_asteroid_number("136199") == 136199

    def test_extract_number_with_name(self):
        """Extract asteroid number from 'number name' format."""
        assert spk._extract_asteroid_number("2060 Chiron") == 2060
        assert spk._extract_asteroid_number("136199 Eris") == 136199

    def test_extract_number_with_parentheses(self):
        """Extract asteroid number from '(number)' format."""
        assert spk._extract_asteroid_number("(2060)") == 2060
        assert spk._extract_asteroid_number("(136199) Eris") == 136199

    def test_extract_number_name_only(self):
        """Return None for name-only identifiers."""
        assert spk._extract_asteroid_number("Chiron") is None
        assert spk._extract_asteroid_number("Eris") is None

    def test_deduce_naif_id(self):
        """Deduce NAIF ID from body identifier."""
        assert spk._deduce_naif_id("2060") == 2060 + NAIF_ASTEROID_OFFSET
        assert spk._deduce_naif_id("2060 Chiron") == 2060 + NAIF_ASTEROID_OFFSET
        assert spk._deduce_naif_id("Chiron") is None


class TestBuildHorizonsUrl:
    """Test Horizons API URL building."""

    def test_basic_url(self):
        """Build basic Horizons URL."""
        url = spk._build_horizons_url("Chiron", "2000-01-01", "2100-01-01")
        assert "COMMAND='Chiron'" in url
        assert "EPHEM_TYPE=SPK" in url
        assert "CENTER='500@0'" in url
        assert "START_TIME='2000-01-01'" in url
        assert "STOP_TIME='2100-01-01'" in url

    def test_custom_center(self):
        """Build URL with custom center."""
        url = spk._build_horizons_url(
            "Chiron", "2000-01-01", "2100-01-01", center="@sun"
        )
        assert "CENTER='@sun'" in url


class TestSanitizeFilename:
    """Test filename sanitization."""

    def test_simple_name(self):
        """Sanitize simple name."""
        assert spk._sanitize_filename("Chiron") == "chiron"

    def test_name_with_spaces(self):
        """Sanitize name with spaces."""
        assert spk._sanitize_filename("2060 Chiron") == "2060_chiron"

    def test_name_with_parentheses(self):
        """Sanitize name with special characters."""
        assert spk._sanitize_filename("(136199) Eris") == "136199_eris"


class TestSpkRegistration:
    """Test SPK body registration."""

    def setup_method(self):
        """Clear state before each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        """Clear state after each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_list_spk_bodies_empty(self):
        """List SPK bodies when none registered."""
        assert eph.list_spk_bodies() == {}

    def test_get_spk_body_info_not_registered(self):
        """Get info for unregistered body."""
        assert eph.get_spk_body_info(SE_CHIRON) is None

    def test_unregister_nonexistent(self):
        """Unregister non-existent body (should not error)."""
        eph.unregister_spk_body(SE_CHIRON)  # Should not raise

    def test_register_file_not_found(self):
        """Register with non-existent file raises SPKNotFoundError."""
        with pytest.raises(SPKNotFoundError) as exc_info:
            eph.register_spk_body(SE_CHIRON, "/nonexistent/file.bsp", NAIF_CHIRON)

        # Verify the error contains helpful information
        error = exc_info.value
        assert "/nonexistent/file.bsp" in str(error)
        assert error.filepath == "/nonexistent/file.bsp"

    def test_register_file_not_found_with_body_name(self):
        """SPKNotFoundError includes body name when available."""
        with pytest.raises(SPKNotFoundError) as exc_info:
            eph.register_spk_body(SE_CHIRON, "/path/to/missing.bsp", NAIF_CHIRON)

        error = exc_info.value
        # Should have Chiron as body name
        assert error.body_name == "Chiron"
        # Should include Horizons ID
        assert error.body_id == "2060"
        # Error message should include instructions
        assert "download_spk" in str(error)

    def test_register_file_not_found_with_eris(self):
        """SPKNotFoundError works for SE_ERIS (TNO with high ID offset)."""
        with pytest.raises(SPKNotFoundError) as exc_info:
            eph.register_spk_body(SE_ERIS, "/path/to/missing_eris.bsp", NAIF_ERIS)

        error = exc_info.value
        assert error.body_name == "Eris"
        assert error.body_id == "136199"

    def test_spk_not_found_error_helpful_message(self):
        """SPKNotFoundError message includes multiple options to obtain SPK."""
        with pytest.raises(SPKNotFoundError) as exc_info:
            eph.register_spk_body(SE_CERES, "/missing/ceres.bsp", NAIF_CERES)

        message = str(exc_info.value)
        # Check that key instructions are included
        assert "download_spk" in message
        assert "download_and_register_spk" in message
        assert "set_auto_spk_download" in message
        assert "libephemeris.scripts.download_spk" in message


class TestCalcWithoutSpk:
    """Test calculations without SPK (Keplerian fallback)."""

    def setup_method(self):
        """Clear SPK state before each test."""
        state._SPK_BODY_MAP.clear()

    def teardown_method(self):
        """Clear SPK state after each test."""
        state._SPK_BODY_MAP.clear()

    def test_chiron_keplerian(self):
        """Calculate Chiron position using Keplerian model."""
        # J2000.0 epoch
        jd = 2451545.0
        pos, flags = eph.calc_ut(jd, SE_CHIRON, SEFLG_SPEED)

        # Should return valid position
        assert 0 <= pos[0] < 360  # Longitude in range
        assert -90 <= pos[1] <= 90  # Latitude in range
        assert pos[2] > 0  # Distance positive

    def test_ceres_keplerian(self):
        """Calculate Ceres position using Keplerian model."""
        jd = 2451545.0
        pos, flags = eph.calc_ut(jd, SE_CERES, SEFLG_SPEED)

        assert 0 <= pos[0] < 360
        assert -90 <= pos[1] <= 90
        assert pos[2] > 0


class TestSpkCalcIntegration:
    """Integration tests for SPK calculations.

    These tests mock the SPK kernel to avoid network dependency.
    """

    def setup_method(self):
        """Clear state before each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        """Clear state after each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_spk_calc_returns_none_when_not_registered(self):
        """calc_spk_body_position returns None for unregistered bodies."""
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        result = spk.calc_spk_body_position(t, SE_CHIRON, 0)
        assert result is None


class TestNaifConstants:
    """Test NAIF ID constants."""

    def test_naif_offset(self):
        """NAIF offset is 2000000."""
        assert NAIF_ASTEROID_OFFSET == 2000000

    def test_naif_chiron(self):
        """Chiron NAIF ID is correct."""
        assert NAIF_CHIRON == 2002060

    def test_naif_ceres(self):
        """Ceres NAIF ID is correct."""
        assert NAIF_CERES == 2000001

    def test_naif_eris(self):
        """Eris NAIF ID is correct."""
        assert NAIF_ERIS == 2136199


class TestEphemerisContextSpk:
    """Test SPK support in EphemerisContext."""

    def setup_method(self):
        """Clear state before each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        """Clear state after each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_context_has_spk_methods(self):
        """EphemerisContext has SPK registration methods."""
        ctx = eph.EphemerisContext()
        assert hasattr(ctx, "register_spk_body")
        assert hasattr(ctx, "unregister_spk_body")
        assert hasattr(ctx, "get_spk_body_info")
        assert hasattr(ctx, "list_spk_bodies")

    def test_context_list_spk_bodies_empty(self):
        """Context lists empty SPK bodies when none registered."""
        ctx = eph.EphemerisContext()
        assert ctx.list_spk_bodies() == {}

    def test_context_get_spk_body_info_not_registered(self):
        """Context returns None for unregistered body."""
        ctx = eph.EphemerisContext()
        assert ctx.get_spk_body_info(SE_CHIRON) is None


# =============================================================================
# NETWORK-DEPENDENT TESTS (Skipped by default)
# =============================================================================


@pytest.mark.skipif(
    os.environ.get("LIBEPHEMERIS_TEST_SPK_DOWNLOAD") != "1",
    reason="Set LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 to run network tests",
)
class TestSpkDownloadIntegration:
    """Integration tests that download actual SPK files.

    These tests require network access and may be slow.
    Skipped by default; set LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 to run.
    """

    def setup_method(self):
        """Clear state before each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        """Clear state after each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_download_chiron_spk(self, tmp_path):
        """Download Chiron SPK from Horizons."""
        path = eph.download_spk(
            body="2060",
            start="2020-01-01",
            end="2025-01-01",
            path=str(tmp_path),
        )

        assert os.path.exists(path)
        assert path.endswith(".bsp")

    def test_download_and_register(self, tmp_path):
        """Download and register Chiron SPK."""
        path = eph.download_and_register_spk(
            body="2060",
            ipl=SE_CHIRON,
            start="2020-01-01",
            end="2025-01-01",
            path=str(tmp_path),
        )

        assert os.path.exists(path)
        assert eph.get_spk_body_info(SE_CHIRON) is not None

        # Calculate position using SPK
        jd = 2459215.5  # 2021-01-01
        pos, _ = eph.calc_ut(jd, SE_CHIRON, SEFLG_SPEED)

        assert 0 <= pos[0] < 360
        assert pos[2] > 0
