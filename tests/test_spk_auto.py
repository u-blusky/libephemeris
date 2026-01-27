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


class TestEnsureCacheDir:
    """Test ensure_cache_dir function."""

    def test_creates_default_cache_dir(self, tmp_path):
        """Create default cache directory if it doesn't exist."""
        with patch.object(spk_auto, "get_library_path", return_value=str(tmp_path)):
            expected_dir = tmp_path / spk_auto.DEFAULT_CACHE_DIR
            assert not expected_dir.exists()

            result = spk_auto.ensure_cache_dir()

            assert expected_dir.exists()
            assert result == str(expected_dir)

    def test_creates_custom_cache_dir(self, tmp_path):
        """Create custom cache directory if it doesn't exist."""
        custom_dir = tmp_path / "my_custom_cache"
        assert not custom_dir.exists()

        result = spk_auto.ensure_cache_dir(str(custom_dir))

        assert custom_dir.exists()
        assert result == str(custom_dir)

    def test_returns_existing_dir(self, tmp_path):
        """Return existing cache directory path."""
        existing_dir = tmp_path / "existing_cache"
        existing_dir.mkdir()

        result = spk_auto.ensure_cache_dir(str(existing_dir))

        assert result == str(existing_dir)

    def test_returns_absolute_path(self, tmp_path, monkeypatch):
        """Always returns absolute path."""
        monkeypatch.chdir(tmp_path)
        relative_path = "relative_cache"

        result = spk_auto.ensure_cache_dir(relative_path)

        assert os.path.isabs(result)
        assert result == os.path.abspath(relative_path)


class TestGetCachePath:
    """Test get_cache_path function."""

    def test_numeric_body_id_int(self, tmp_path):
        """Get cache path for numeric body ID (int)."""
        with patch.object(spk_auto, "get_library_path", return_value=str(tmp_path)):
            result = spk_auto.get_cache_path(2060)

            assert result.endswith("2060.bsp")
            assert spk_auto.DEFAULT_CACHE_DIR in result

    def test_numeric_body_id_str(self, tmp_path):
        """Get cache path for numeric body ID (string)."""
        with patch.object(spk_auto, "get_library_path", return_value=str(tmp_path)):
            result = spk_auto.get_cache_path("136199")

            assert result.endswith("136199.bsp")

    def test_named_body_id(self, tmp_path):
        """Get cache path for named body."""
        with patch.object(spk_auto, "get_library_path", return_value=str(tmp_path)):
            result = spk_auto.get_cache_path("Chiron")

            assert result.endswith("chiron.bsp")

    def test_body_id_with_spaces(self, tmp_path):
        """Sanitize body ID with spaces."""
        with patch.object(spk_auto, "get_library_path", return_value=str(tmp_path)):
            result = spk_auto.get_cache_path("2060 Chiron")

            assert " " not in result
            assert result.endswith(".bsp")
            assert "2060_chiron" in result

    def test_custom_cache_dir(self, tmp_path):
        """Get cache path with custom cache directory."""
        custom_dir = tmp_path / "custom"
        result = spk_auto.get_cache_path("2060", str(custom_dir))

        assert str(custom_dir) in result
        assert result.endswith("2060.bsp")

    def test_returns_absolute_path(self, tmp_path):
        """Always returns absolute path."""
        result = spk_auto.get_cache_path("2060", str(tmp_path))

        assert os.path.isabs(result)


class TestCacheInfoStandalone:
    """Test cache_info standalone function."""

    def test_empty_cache(self, tmp_path):
        """Return empty info when cache doesn't exist."""
        non_existent = tmp_path / "nonexistent"
        info = spk_auto.cache_info(str(non_existent))

        assert info["cache_dir"] == str(non_existent)
        assert info["num_files"] == 0
        assert info["total_size_mb"] == 0.0
        assert info["files"] == []

    def test_with_files(self, tmp_path):
        """Return correct info when cache has files."""
        # Create some dummy SPK files (large enough to register non-zero MB)
        (tmp_path / "body1.bsp").write_bytes(b"x" * (1024 * 1024))  # 1 MB
        (tmp_path / "body2.bsp").write_bytes(b"y" * (1024 * 512))  # 0.5 MB
        (tmp_path / "not_spk.txt").write_bytes(b"ignored")

        info = spk_auto.cache_info(str(tmp_path))

        assert info["cache_dir"] == str(tmp_path)
        assert info["num_files"] == 2
        assert "body1.bsp" in info["files"]
        assert "body2.bsp" in info["files"]
        assert "not_spk.txt" not in info["files"]
        assert info["total_size_mb"] >= 1.0  # At least 1 MB

    def test_default_cache_dir(self, tmp_path):
        """Use default cache directory when none specified."""
        with patch.object(spk_auto, "get_library_path", return_value=str(tmp_path)):
            info = spk_auto.cache_info()

            expected_dir = os.path.join(str(tmp_path), spk_auto.DEFAULT_CACHE_DIR)
            assert info["cache_dir"] == expected_dir

    def test_returns_absolute_path(self, tmp_path, monkeypatch):
        """Always returns absolute path in cache_dir."""
        monkeypatch.chdir(tmp_path)
        relative_path = "relative_cache"
        os.makedirs(relative_path, exist_ok=True)

        info = spk_auto.cache_info(relative_path)

        assert os.path.isabs(info["cache_dir"])


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


# =============================================================================
# TESTS FOR auto_get_spk FUNCTION
# =============================================================================


class TestJdToIsoDate:
    """Test Julian Day to ISO date conversion."""

    def test_j2000_epoch(self):
        """Test J2000.0 epoch conversion."""
        # J2000.0 = 2451545.0 = 2000-01-01 12:00 TT
        # The function converts to date (ignoring time)
        result = spk_auto._jd_to_iso_date(2451545.0)
        assert result == "2000-01-01"

    def test_2020_epoch(self):
        """Test a date in 2020."""
        # 2020-01-01 = JD 2458849.5
        result = spk_auto._jd_to_iso_date(2458849.5)
        assert result == "2020-01-01"

    def test_2030_epoch(self):
        """Test a date in 2030."""
        # 2030-01-01 = JD 2462502.5
        result = spk_auto._jd_to_iso_date(2462502.5)
        assert result == "2030-01-01"


class TestGenerateSpkCacheFilename:
    """Test SPK cache filename generation."""

    def test_basic_filename(self):
        """Test basic filename generation."""
        filename = spk_auto._generate_spk_cache_filename("2060", 2458849.5, 2462502.5)
        assert filename == "2060_2458849_2462502.bsp"

    def test_filename_with_name(self):
        """Test filename with body name."""
        filename = spk_auto._generate_spk_cache_filename("Chiron", 2458849.5, 2462502.5)
        assert filename == "chiron_2458849_2462502.bsp"

    def test_filename_sanitization(self):
        """Test that special characters are sanitized."""
        filename = spk_auto._generate_spk_cache_filename(
            "2060 Chiron", 2458849.5, 2462502.5
        )
        assert " " not in filename
        assert filename.endswith(".bsp")


class TestFindCoveringSpk:
    """Test finding existing SPK files that cover a date range."""

    def test_no_cache_dir(self, tmp_path):
        """Returns None if cache directory doesn't exist."""
        result = spk_auto._find_covering_spk(
            "2060", 2458849.5, 2462502.5, str(tmp_path / "nonexistent")
        )
        assert result is None

    def test_empty_cache_dir(self, tmp_path):
        """Returns None if cache directory is empty."""
        result = spk_auto._find_covering_spk(
            "2060", 2458849.5, 2462502.5, str(tmp_path)
        )
        assert result is None

    def test_finds_exact_match(self, tmp_path):
        """Finds an SPK file with exact matching range."""
        # Create a dummy SPK file with matching name
        spk_file = tmp_path / "2060_2458849_2462502.bsp"
        spk_file.write_bytes(b"dummy")

        result = spk_auto._find_covering_spk(
            "2060", 2458849.5, 2462502.5, str(tmp_path)
        )
        assert result == str(spk_file)

    def test_finds_covering_file(self, tmp_path):
        """Finds an SPK file that covers the requested range."""
        # Create a file with wider range
        spk_file = tmp_path / "2060_2450000_2470000.bsp"
        spk_file.write_bytes(b"dummy")

        result = spk_auto._find_covering_spk(
            "2060", 2458849.5, 2462502.5, str(tmp_path)
        )
        assert result == str(spk_file)

    def test_ignores_non_covering_file(self, tmp_path):
        """Ignores an SPK file that doesn't cover the range."""
        # Create a file with narrower range
        spk_file = tmp_path / "2060_2459000_2460000.bsp"
        spk_file.write_bytes(b"dummy")

        result = spk_auto._find_covering_spk(
            "2060", 2458849.5, 2462502.5, str(tmp_path)
        )
        assert result is None

    def test_ignores_different_body(self, tmp_path):
        """Ignores SPK files for different bodies."""
        # Create a file for different body
        spk_file = tmp_path / "5145_2458849_2462502.bsp"
        spk_file.write_bytes(b"dummy")

        result = spk_auto._find_covering_spk(
            "2060", 2458849.5, 2462502.5, str(tmp_path)
        )
        assert result is None


class TestAutoGetSpkValidation:
    """Test auto_get_spk input validation."""

    @patch.object(spk_auto, "_check_astroquery_available", return_value=False)
    def test_raises_import_error_without_astroquery(self, mock_check, tmp_path):
        """Raises ImportError when astroquery is not available."""
        with pytest.raises(ImportError) as exc_info:
            spk_auto.auto_get_spk("2060", 2458849.5, 2462502.5, str(tmp_path))

        assert "astroquery" in str(exc_info.value)

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_raises_value_error_for_invalid_range(self, mock_check, tmp_path):
        """Raises ValueError when jd_end <= jd_start."""
        with pytest.raises(ValueError) as exc_info:
            spk_auto.auto_get_spk("2060", 2462502.5, 2458849.5, str(tmp_path))

        assert "must be greater than" in str(exc_info.value)


class TestAutoGetSpkCaching:
    """Test auto_get_spk caching behavior."""

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_returns_cached_file(self, mock_check, tmp_path):
        """Returns existing cached file without downloading."""
        # Create a cached SPK file
        spk_file = tmp_path / "2060_2458849_2462502.bsp"
        spk_file.write_bytes(b"cached SPK data")

        # Should return cached file without calling download
        with patch.object(spk_auto, "_download_spk_astroquery") as mock_download:
            result = spk_auto.auto_get_spk("2060", 2458849.5, 2462502.5, str(tmp_path))

            assert result == str(spk_file)
            mock_download.assert_not_called()

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    def test_uses_covering_file(self, mock_check, tmp_path):
        """Uses a file that covers the requested range."""
        # Create a file with wider range
        wide_file = tmp_path / "2060_2450000_2470000.bsp"
        wide_file.write_bytes(b"wide range SPK")

        with patch.object(spk_auto, "_download_spk_astroquery") as mock_download:
            result = spk_auto.auto_get_spk("2060", 2458849.5, 2462502.5, str(tmp_path))

            assert result == str(wide_file)
            mock_download.assert_not_called()

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    @patch.object(spk_auto, "_download_spk_astroquery")
    def test_downloads_when_not_cached(self, mock_download, mock_check, tmp_path):
        """Downloads SPK when not in cache."""
        result = spk_auto.auto_get_spk("2060", 2458849.5, 2462502.5, str(tmp_path))

        # Should have called download
        mock_download.assert_called_once()

        # Check that the correct dates were passed
        call_args = mock_download.call_args
        assert call_args[1]["body_id"] == "2060"
        assert call_args[1]["start"] == "2020-01-01"
        assert call_args[1]["end"] == "2030-01-01"

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    @patch.object(spk_auto, "_download_spk_astroquery")
    def test_creates_cache_directory(self, mock_download, mock_check, tmp_path):
        """Creates cache directory if it doesn't exist."""
        cache_dir = tmp_path / "new_cache"
        assert not cache_dir.exists()

        spk_auto.auto_get_spk("2060", 2458849.5, 2462502.5, str(cache_dir))

        assert cache_dir.exists()


class TestAutoGetSpkDefaultDirectory:
    """Test auto_get_spk default cache directory."""

    def test_default_cache_dir_constant(self):
        """Default cache directory is set correctly."""
        expected = os.path.join(os.path.expanduser("~"), ".libephemeris", "spk")
        assert spk_auto.DEFAULT_AUTO_SPK_DIR == expected

    @patch.object(spk_auto, "_check_astroquery_available", return_value=True)
    @patch.object(spk_auto, "_download_spk_astroquery")
    def test_uses_default_directory(
        self, mock_download, mock_check, tmp_path, monkeypatch
    ):
        """Uses default directory when none specified."""
        # Monkeypatch the default directory
        test_default = str(tmp_path / "default_cache")
        monkeypatch.setattr(spk_auto, "DEFAULT_AUTO_SPK_DIR", test_default)

        spk_auto.auto_get_spk("2060", 2458849.5, 2462502.5)

        # Check that the output path is in the default directory
        call_args = mock_download.call_args
        assert test_default in call_args[1]["output_path"]


# Network-dependent tests for auto_get_spk
@pytest.mark.skipif(
    os.environ.get("LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD") != "1",
    reason="Set LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1 to run network tests",
)
class TestAutoGetSpkIntegration:
    """Integration tests for auto_get_spk that require network access."""

    def test_auto_get_spk_chiron(self, tmp_path):
        """Download Chiron SPK via auto_get_spk."""
        jd_start = 2459215.5  # 2021-01-01
        jd_end = 2460676.5  # 2025-01-01

        spk_path = spk_auto.auto_get_spk("2060", jd_start, jd_end, str(tmp_path))

        assert os.path.exists(spk_path)
        assert spk_path.endswith(".bsp")
        assert os.path.getsize(spk_path) > 0

    def test_auto_get_spk_cached(self, tmp_path):
        """Second call returns cached file."""
        jd_start = 2459215.5
        jd_end = 2460676.5

        # First call downloads
        path1 = spk_auto.auto_get_spk("2060", jd_start, jd_end, str(tmp_path))

        # Second call should return cached file
        path2 = spk_auto.auto_get_spk("2060", jd_start, jd_end, str(tmp_path))

        assert path1 == path2


# =============================================================================
# TESTS FOR is_spk_cached FUNCTION
# =============================================================================


class TestIsSpkCached:
    """Test is_spk_cached function."""

    def test_returns_false_when_cache_dir_not_exists(self, tmp_path):
        """Returns False if cache directory doesn't exist."""
        result = spk_auto.is_spk_cached(
            "2060", 2458849.5, 2462502.5, str(tmp_path / "nonexistent")
        )
        assert result is False

    def test_returns_false_when_cache_empty(self, tmp_path):
        """Returns False if cache directory is empty."""
        result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5, str(tmp_path))
        assert result is False

    def test_returns_false_when_no_matching_body(self, tmp_path):
        """Returns False if no SPK file matches the body."""
        # Create a file for different body
        spk_file = tmp_path / "5145_2450000_2470000.bsp"
        spk_file.write_bytes(b"dummy")

        result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5, str(tmp_path))
        assert result is False

    def test_returns_false_when_file_invalid(self, tmp_path):
        """Returns False if SPK file cannot be parsed."""
        # Create a file with matching name but invalid content
        spk_file = tmp_path / "2060_2450000_2470000.bsp"
        spk_file.write_bytes(b"not a valid SPK file")

        result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5, str(tmp_path))
        assert result is False

    def test_uses_default_cache_dir(self, tmp_path, monkeypatch):
        """Uses default cache directory when none specified."""
        # Monkeypatch the default directory
        test_default = str(tmp_path / "default_cache")
        monkeypatch.setattr(spk_auto, "DEFAULT_AUTO_SPK_DIR", test_default)

        # Should return False since directory doesn't exist
        result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5)
        assert result is False

    def test_sanitizes_body_id_with_spaces(self, tmp_path):
        """Sanitizes body ID with spaces for filename matching."""
        # Create a file with sanitized name
        spk_file = tmp_path / "2060_chiron_2450000_2470000.bsp"
        spk_file.write_bytes(b"dummy")

        # Should attempt to match (will fail because file is invalid, but shouldn't crash)
        result = spk_auto.is_spk_cached(
            "2060 Chiron", 2458849.5, 2462502.5, str(tmp_path)
        )
        assert result is False  # False because file content is invalid

    def test_sanitizes_body_id_lowercase(self, tmp_path):
        """Body ID is lowercased for filename matching."""
        # The function should match "chiron" regardless of input case
        spk_file = tmp_path / "chiron_2450000_2470000.bsp"
        spk_file.write_bytes(b"dummy")

        # Should attempt to match (will fail because file is invalid, but shouldn't crash)
        result = spk_auto.is_spk_cached("Chiron", 2458849.5, 2462502.5, str(tmp_path))
        assert result is False

    def test_ignores_non_bsp_files(self, tmp_path):
        """Ignores non-.bsp files in cache directory."""
        # Create a file with wrong extension
        txt_file = tmp_path / "2060_2450000_2470000.txt"
        txt_file.write_bytes(b"not an SPK file")

        result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5, str(tmp_path))
        assert result is False


class TestIsSpkCachedWithMockedCoverage:
    """Test is_spk_cached with mocked get_spk_coverage."""

    def test_returns_true_when_spk_covers_range(self, tmp_path):
        """Returns True when SPK file covers the requested range."""
        # Create a dummy SPK file
        spk_file = tmp_path / "2060_2450000_2470000.bsp"
        spk_file.write_bytes(b"dummy SPK")

        # Mock get_spk_coverage in the spk module (where it's imported from)
        with patch("libephemeris.spk.get_spk_coverage") as mock_coverage:
            mock_coverage.return_value = (2450000.0, 2470000.0)

            result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5, str(tmp_path))
            assert result is True
            mock_coverage.assert_called()

    def test_returns_false_when_spk_range_too_narrow(self, tmp_path):
        """Returns False when SPK file doesn't cover the full range."""
        # Create a dummy SPK file
        spk_file = tmp_path / "2060_2459000_2460000.bsp"
        spk_file.write_bytes(b"dummy SPK")

        # Mock get_spk_coverage to return narrower range
        with patch("libephemeris.spk.get_spk_coverage") as mock_coverage:
            mock_coverage.return_value = (2459000.0, 2460000.0)

            result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5, str(tmp_path))
            assert result is False

    def test_returns_false_when_coverage_returns_none(self, tmp_path):
        """Returns False when get_spk_coverage returns None."""
        # Create a dummy SPK file
        spk_file = tmp_path / "2060_2450000_2470000.bsp"
        spk_file.write_bytes(b"dummy SPK")

        # Mock get_spk_coverage to return None
        with patch("libephemeris.spk.get_spk_coverage") as mock_coverage:
            mock_coverage.return_value = None

            result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5, str(tmp_path))
            assert result is False

    def test_tries_multiple_files(self, tmp_path):
        """Tries multiple SPK files until one covers the range."""
        # Create multiple SPK files
        spk_file1 = tmp_path / "2060_2459000_2460000.bsp"
        spk_file1.write_bytes(b"narrow range")
        spk_file2 = tmp_path / "2060_2450000_2470000.bsp"
        spk_file2.write_bytes(b"wide range")

        call_count = 0
        files_checked = []

        def mock_coverage_side_effect(path):
            nonlocal call_count
            call_count += 1
            files_checked.append(path)
            if "2459000" in path:
                return (2459000.0, 2460000.0)  # Too narrow
            elif "2450000" in path:
                return (2450000.0, 2470000.0)  # Covers range
            return None

        with patch(
            "libephemeris.spk.get_spk_coverage",
            side_effect=mock_coverage_side_effect,
        ):
            result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5, str(tmp_path))
            assert result is True
            # Should have checked at least one file
            assert call_count >= 1

    def test_handles_coverage_exception(self, tmp_path):
        """Handles exceptions from get_spk_coverage gracefully."""
        # Create a dummy SPK file
        spk_file = tmp_path / "2060_2450000_2470000.bsp"
        spk_file.write_bytes(b"dummy SPK")

        # Mock get_spk_coverage to raise exception
        with patch("libephemeris.spk.get_spk_coverage") as mock_coverage:
            mock_coverage.side_effect = Exception("Failed to parse SPK")

            # Should not raise, just return False
            result = spk_auto.is_spk_cached("2060", 2458849.5, 2462502.5, str(tmp_path))
            assert result is False


class TestIsSpkCachedModuleExport:
    """Test that is_spk_cached is properly exported."""

    def test_is_spk_cached_accessible(self):
        """is_spk_cached is accessible from spk_auto module."""
        assert hasattr(spk_auto, "is_spk_cached")
        assert callable(spk_auto.is_spk_cached)

    def test_is_spk_cached_via_libephemeris(self):
        """is_spk_cached is accessible via libephemeris.spk_auto."""
        assert hasattr(eph.spk_auto, "is_spk_cached")
        assert callable(eph.spk_auto.is_spk_cached)
