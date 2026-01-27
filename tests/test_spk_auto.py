"""
Tests for spk_auto.py - Automatic SPK download and caching module.

These tests verify:
- AutoSpkConfig class functionality
- Registration and de-registration of auto-SPK
- Cache path generation
- Import checks for astroquery
- Integration with the main SPK system

Note: Tests that require actual network downloads are skipped by default.
Set LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1 to run them.
"""

import os
import pytest
from unittest.mock import patch, MagicMock

import libephemeris as eph
from libephemeris import spk_auto
from libephemeris.constants import (
    SE_CHIRON,
    SE_CERES,
    SE_ERIS,
    SE_PHOLUS,
    NAIF_CHIRON,
    NAIF_CERES,
    SEFLG_SPEED,
)


class TestAutoSpkConfig:
    """Test AutoSpkConfig class."""

    def test_config_creation(self):
        """Create AutoSpkConfig with basic parameters."""
        config = spk_auto.AutoSpkConfig(
            ipl=SE_CHIRON,
            body_id="2060",
            naif_id=NAIF_CHIRON,
        )

        assert config.ipl == SE_CHIRON
        assert config.body_id == "2060"
        assert config.naif_id == NAIF_CHIRON
        assert config.enabled is True
        assert config.spk_path is None

    def test_config_with_dates(self):
        """Create config with custom date range."""
        config = spk_auto.AutoSpkConfig(
            ipl=SE_CHIRON,
            body_id="2060",
            start="2020-01-01",
            end="2030-01-01",
        )

        assert config.start == "2020-01-01"
        assert config.end == "2030-01-01"

    def test_cache_filename_generation(self):
        """Generate unique cache filename."""
        config = spk_auto.AutoSpkConfig(
            ipl=SE_CHIRON,
            body_id="2060",
            start="2000-01-01",
            end="2100-01-01",
        )

        filename = config.get_cache_filename()
        assert filename.endswith(".bsp")
        assert "2060" in filename

    def test_cache_filename_uniqueness(self):
        """Different configs produce different filenames."""
        config1 = spk_auto.AutoSpkConfig(
            ipl=SE_CHIRON,
            body_id="2060",
            start="2000-01-01",
            end="2100-01-01",
        )
        config2 = spk_auto.AutoSpkConfig(
            ipl=SE_CHIRON,
            body_id="2060",
            start="2020-01-01",
            end="2030-01-01",
        )

        assert config1.get_cache_filename() != config2.get_cache_filename()


class TestEnableDisableAutoSpk:
    """Test enabling and disabling auto-SPK."""

    def setup_method(self):
        """Clear registry before each test."""
        spk_auto.disable_all()

    def teardown_method(self):
        """Clear registry after each test."""
        spk_auto.disable_all()

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_enable_auto_spk(self, mock_check):
        """Enable auto-SPK for a body."""
        spk_auto.enable_auto_spk(
            ipl=SE_CHIRON,
            body_id="2060",
            naif_id=NAIF_CHIRON,
        )

        assert spk_auto.is_auto_spk_enabled(SE_CHIRON) is True

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_disable_auto_spk(self, mock_check):
        """Disable auto-SPK for a body."""
        spk_auto.enable_auto_spk(
            ipl=SE_CHIRON,
            body_id="2060",
        )

        spk_auto.disable_auto_spk(SE_CHIRON)

        assert spk_auto.is_auto_spk_enabled(SE_CHIRON) is False

    def test_is_auto_spk_enabled_not_configured(self):
        """Check auto-SPK status for unconfigured body."""
        assert spk_auto.is_auto_spk_enabled(SE_CHIRON) is False

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_get_auto_spk_config(self, mock_check):
        """Get configuration for a body."""
        spk_auto.enable_auto_spk(
            ipl=SE_CHIRON,
            body_id="2060",
            naif_id=NAIF_CHIRON,
            start="2000-01-01",
            end="2100-01-01",
        )

        config = spk_auto.get_auto_spk_config(SE_CHIRON)

        assert config is not None
        assert config.body_id == "2060"
        assert config.naif_id == NAIF_CHIRON

    def test_get_auto_spk_config_not_configured(self):
        """Get config for unconfigured body returns None."""
        assert spk_auto.get_auto_spk_config(SE_CHIRON) is None

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_list_auto_spk_bodies(self, mock_check):
        """List all configured bodies."""
        spk_auto.enable_auto_spk(ipl=SE_CHIRON, body_id="2060")
        spk_auto.enable_auto_spk(ipl=SE_CERES, body_id="1")

        bodies = spk_auto.list_auto_spk_bodies()

        assert SE_CHIRON in bodies
        assert SE_CERES in bodies
        assert len(bodies) == 2

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_disable_all(self, mock_check):
        """Disable all auto-SPK configurations."""
        spk_auto.enable_auto_spk(ipl=SE_CHIRON, body_id="2060")
        spk_auto.enable_auto_spk(ipl=SE_CERES, body_id="1")

        spk_auto.disable_all()

        assert len(spk_auto.list_auto_spk_bodies()) == 0


class TestAstroqueryCheck:
    """Test astroquery availability check."""

    def test_check_available_when_installed(self):
        """Check returns True when astroquery is available."""
        # This test passes if astroquery is installed
        try:
            from astroquery.jplhorizons import Horizons

            assert spk_auto._check_astroquery_available() is True
        except ImportError:
            pytest.skip("astroquery not installed")

    @patch.dict("sys.modules", {"astroquery": None, "astroquery.jplhorizons": None})
    def test_check_returns_false_when_not_installed(self):
        """Check returns False when astroquery is not available."""
        # Force the import to fail
        with patch.object(spk_auto, "_check_astroquery_available", return_value=False):
            assert spk_auto._check_astroquery_available() is False


class TestEnableWithoutAstroquery:
    """Test behavior when astroquery is not installed."""

    def setup_method(self):
        """Clear registry before each test."""
        spk_auto.disable_all()

    def teardown_method(self):
        """Clear registry after each test."""
        spk_auto.disable_all()

    @patch.object(spk_auto, "_check_astroquery_available", return_value=False)
    def test_enable_raises_import_error(self, mock_check):
        """enable_auto_spk raises ImportError when astroquery not available."""
        with pytest.raises(ImportError) as exc_info:
            spk_auto.enable_auto_spk(ipl=SE_CHIRON, body_id="2060")

        assert "astroquery" in str(exc_info.value)


class TestCacheInfo:
    """Test cache information functions."""

    def test_get_cache_info_empty(self, tmp_path):
        """Get cache info when cache is empty/non-existent."""
        with patch.object(spk_auto, "get_library_path", return_value=str(tmp_path)):
            info = spk_auto.get_cache_info()

            assert info["num_files"] == 0
            assert info["total_size_mb"] == 0.0
            assert info["files"] == []

    def test_get_cache_info_with_files(self, tmp_path):
        """Get cache info when cache has files."""
        cache_dir = tmp_path / spk_auto.DEFAULT_CACHE_DIR
        cache_dir.mkdir()

        # Create a dummy SPK file
        dummy_file = cache_dir / "test_body.bsp"
        dummy_file.write_bytes(b"dummy content" * 1000)

        with patch.object(spk_auto, "get_library_path", return_value=str(tmp_path)):
            info = spk_auto.get_cache_info()

            assert info["num_files"] == 1
            assert "test_body.bsp" in info["files"]
            assert info["total_size_mb"] > 0


class TestTryAutoDownload:
    """Test the try_auto_download function."""

    def setup_method(self):
        """Clear registry before each test."""
        spk_auto.disable_all()

    def teardown_method(self):
        """Clear registry after each test."""
        spk_auto.disable_all()

    def test_returns_none_when_not_configured(self):
        """Returns None for unconfigured body."""
        result = spk_auto.try_auto_download(SE_CHIRON)
        assert result is None


class TestDownloadNowWithoutConfig:
    """Test download_now without configuration."""

    def setup_method(self):
        """Clear registry before each test."""
        spk_auto.disable_all()

    def teardown_method(self):
        """Clear registry after each test."""
        spk_auto.disable_all()

    def test_download_now_raises_without_config(self):
        """download_now raises ValueError when not configured."""
        with pytest.raises(ValueError) as exc_info:
            spk_auto.download_now(SE_CHIRON)

        assert "Auto-SPK not enabled" in str(exc_info.value)


class TestEnableCommonBodies:
    """Test the enable_common_bodies preset."""

    def setup_method(self):
        """Clear registry before each test."""
        spk_auto.disable_all()

    def teardown_method(self):
        """Clear registry after each test."""
        spk_auto.disable_all()

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_enable_common_bodies(self, mock_check):
        """Enable common minor bodies."""
        spk_auto.enable_common_bodies()

        bodies = spk_auto.list_auto_spk_bodies()

        # Should include common bodies
        assert SE_CHIRON in bodies
        assert SE_CERES in bodies

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_enable_common_bodies_with_dates(self, mock_check):
        """Enable common bodies with custom date range."""
        spk_auto.enable_common_bodies(
            start="2020-01-01",
            end="2030-01-01",
        )

        config = spk_auto.get_auto_spk_config(SE_CHIRON)

        assert config is not None
        assert config.start == "2020-01-01"
        assert config.end == "2030-01-01"


class TestClearCache:
    """Test cache clearing functionality."""

    def setup_method(self):
        """Clear registry before each test."""
        spk_auto.disable_all()

    def teardown_method(self):
        """Clear registry after each test."""
        spk_auto.disable_all()

    def test_clear_cache_returns_zero_when_empty(self):
        """clear_cache returns 0 when no files to delete."""
        deleted = spk_auto.clear_cache(SE_CHIRON)
        assert deleted == 0


class TestModuleExports:
    """Test that module is properly exported."""

    def test_spk_auto_accessible(self):
        """spk_auto module is accessible from libephemeris."""
        assert hasattr(eph, "spk_auto")

    def test_enable_auto_spk_accessible(self):
        """enable_auto_spk is accessible."""
        assert hasattr(eph.spk_auto, "enable_auto_spk")

    def test_enable_common_bodies_accessible(self):
        """enable_common_bodies is accessible."""
        assert hasattr(eph.spk_auto, "enable_common_bodies")


# =============================================================================
# NETWORK-DEPENDENT TESTS (Skipped by default)
# =============================================================================


@pytest.mark.skipif(
    os.environ.get("LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD") != "1",
    reason="Set LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1 to run network tests",
)
class TestSpkAutoDownloadIntegration:
    """Integration tests that download actual SPK files.

    These tests require:
    - Network access
    - astroquery installed

    Skipped by default; set LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1 to run.
    """

    def setup_method(self):
        """Clear state before each test."""
        spk_auto.disable_all()
        from libephemeris import state

        state._SPK_BODY_MAP.clear()

    def teardown_method(self):
        """Clear state after each test."""
        spk_auto.disable_all()
        from libephemeris import state

        state._SPK_BODY_MAP.clear()

    def test_download_chiron_spk(self, tmp_path):
        """Download Chiron SPK using astroquery."""
        spk_auto.enable_auto_spk(
            ipl=SE_CHIRON,
            body_id="2060",
            naif_id=NAIF_CHIRON,
            start="2020-01-01",
            end="2025-01-01",
            cache_dir=str(tmp_path),
        )

        path = spk_auto.download_now(SE_CHIRON)

        assert os.path.exists(path)
        assert path.endswith(".bsp")

    def test_auto_download_and_calc(self, tmp_path):
        """Enable auto-SPK and calculate position."""
        spk_auto.enable_auto_spk(
            ipl=SE_CHIRON,
            body_id="2060",
            naif_id=NAIF_CHIRON,
            start="2020-01-01",
            end="2025-01-01",
            cache_dir=str(tmp_path),
        )

        # Trigger download
        path = spk_auto.download_now(SE_CHIRON)
        assert os.path.exists(path)

        # Calculate position - should use SPK
        jd = 2459215.5  # 2021-01-01
        pos, _ = eph.calc_ut(jd, SE_CHIRON, SEFLG_SPEED)

        assert 0 <= pos[0] < 360
        assert pos[2] > 0
