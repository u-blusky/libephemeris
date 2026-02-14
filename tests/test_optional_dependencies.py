"""
Tests for optional dependencies documentation and availability checking.

This module tests that:
1. Optional dependencies are correctly documented in README.md
2. Optional dependencies are defined in pyproject.toml extras
3. The library gracefully handles missing optional dependencies
4. Optional features provide clear error messages when dependencies are missing
"""

import re
import sys
from pathlib import Path
from unittest.mock import patch

import pytest


class TestOptionalDependenciesDocumentation:
    """Tests for optional dependencies documentation in README.md."""

    @pytest.fixture
    def readme_content(self) -> str:
        """Load README.md content."""
        readme_path = Path(__file__).parent.parent / "README.md"
        return readme_path.read_text()

    def test_pyerfa_documented_in_readme(self, readme_content: str):
        """pyerfa should be documented in the README."""
        assert "pyerfa" in readme_content
        assert "IAU 2006" in readme_content or "required" in readme_content.lower()

    def test_astroquery_documented_in_readme(self, readme_content: str):
        """astroquery should be documented as an optional dependency."""
        assert "astroquery" in readme_content
        assert "SPK" in readme_content or "Horizons" in readme_content

    def test_astropy_documented_in_readme(self, readme_content: str):
        """astropy should be documented as an optional dependency."""
        assert "astropy" in readme_content
        assert (
            "star catalog" in readme_content.lower() or "star_catalog" in readme_content
        )

    def test_optional_dependencies_section_exists(self, readme_content: str):
        """README should have an Optional Dependencies section."""
        assert "Optional Dependencies" in readme_content

    def test_installation_commands_documented(self, readme_content: str):
        """README should include installation commands for optional dependencies."""
        # Check for pip install extras syntax
        assert "pip install libephemeris[" in readme_content


class TestPyprojectOptionalDependencies:
    """Tests for optional dependencies in pyproject.toml."""

    @pytest.fixture
    def pyproject_content(self) -> str:
        """Load pyproject.toml content."""
        pyproject_path = Path(__file__).parent.parent / "pyproject.toml"
        return pyproject_path.read_text()

    def test_pyerfa_is_required_dependency(self, pyproject_content: str):
        """pyproject.toml should list pyerfa as a required dependency."""
        # pyerfa was promoted from optional [precision] to required dependency
        assert "pyerfa" in pyproject_content

    def test_astroquery_is_required_dependency(self, pyproject_content: str):
        """pyproject.toml should list astroquery as a required dependency."""
        assert "astroquery" in pyproject_content

    def test_stars_extra_defined(self, pyproject_content: str):
        """pyproject.toml should define a 'stars' extra with astropy."""
        assert "stars" in pyproject_content
        assert "astropy" in pyproject_content

    def test_all_extra_defined(self, pyproject_content: str):
        """pyproject.toml should define an 'all' extra."""
        # Check for 'all' extra that includes optional dependencies
        lines = pyproject_content.split("\n")
        in_all_section = False
        all_deps = []
        for line in lines:
            if line.strip().startswith("all = ["):
                in_all_section = True
            elif in_all_section:
                if line.strip() == "]":
                    break
                all_deps.append(line.strip())

        # 'all' extra should include rebound (N-body integration)
        all_deps_str = " ".join(all_deps)
        assert "rebound" in all_deps_str


class TestPyerfaRequiredImport:
    """Tests for pyerfa required dependency (promoted from optional)."""

    def test_lunar_module_imports(self):
        """lunar module should import successfully."""
        from libephemeris import lunar

        assert lunar is not None

    def test_true_node_works(self):
        """True node calculation should work with pyerfa."""
        from libephemeris import calc_ut, SE_TRUE_NODE

        jd = 2451545.0  # J2000.0
        result, _ = calc_ut(jd, SE_TRUE_NODE, 0)

        # Should return a valid longitude
        assert 0 <= result[0] < 360

    def test_erfa_import_exists_in_core_modules(self):
        """Core modules should import erfa directly (required dependency)."""
        from libephemeris import cache

        cache_source = Path(cache.__file__).read_text()
        assert "import erfa" in cache_source


class TestAstroqueryOptionalImport:
    """Tests for astroquery optional dependency handling."""

    def test_spk_auto_module_imports_without_astroquery(self):
        """spk_auto module should import even without astroquery."""
        from libephemeris import spk_auto

        assert spk_auto is not None

    def test_astroquery_availability_check_function_exists(self):
        """spk_auto should have a function to check astroquery availability."""
        from libephemeris import spk_auto

        assert hasattr(spk_auto, "_check_astroquery_available")
        # The function should return a boolean
        result = spk_auto._check_astroquery_available()
        assert isinstance(result, bool)

    def test_enable_auto_spk_raises_without_astroquery(self):
        """enable_auto_spk should raise ImportError when astroquery is not available."""
        from libephemeris import spk_auto
        from libephemeris.constants import SE_CHIRON

        # Mock astroquery as unavailable
        with patch.object(spk_auto, "_check_astroquery_available", return_value=False):
            with pytest.raises(ImportError) as exc_info:
                spk_auto.enable_auto_spk(
                    ipl=SE_CHIRON,
                    body_id="2060",
                )
            assert "astroquery" in str(exc_info.value)


class TestOptionalDependencyErrorMessages:
    """Tests for helpful error messages when optional dependencies are missing."""

    def test_spk_download_error_message_mentions_astroquery(self):
        """SPK download functions should mention astroquery in error messages."""
        from libephemeris import spk_auto

        with patch.object(spk_auto, "_check_astroquery_available", return_value=False):
            with pytest.raises(ImportError) as exc_info:
                spk_auto.auto_get_spk(
                    body_id="2060",
                    jd_start=2451545.0,
                    jd_end=2488069.5,
                )
            error_msg = str(exc_info.value).lower()
            assert "astroquery" in error_msg
            assert "pip install" in error_msg

    def test_download_spk_from_horizons_error_message(self):
        """download_spk_from_horizons should provide helpful error message."""
        from libephemeris import spk_auto

        with patch.object(spk_auto, "_check_astroquery_available", return_value=False):
            with pytest.raises(ImportError) as exc_info:
                spk_auto.download_spk_from_horizons(
                    body_id="2060",
                    jd_start=2451545.0,
                    jd_end=2488069.5,
                    output_path="/tmp/test.bsp",
                )
            error_msg = str(exc_info.value).lower()
            assert "astroquery" in error_msg


class TestOptionalDependencyFeatures:
    """Tests that optional features work correctly when dependencies are available."""

    def test_pyerfa_enhances_lunar_node_precision(self):
        """When pyerfa is available, it should be used for lunar node calculations."""
        try:
            import erfa

            has_erfa = True
        except ImportError:
            has_erfa = False

        if not has_erfa:
            pytest.skip("pyerfa not installed")

        # If pyerfa is available, lunar calculations should still work
        from libephemeris import calc_ut, SE_TRUE_NODE

        jd = 2451545.0
        result, _ = calc_ut(jd, SE_TRUE_NODE, 0)
        assert 0 <= result[0] < 360

    def test_astroquery_check_returns_true_when_installed(self):
        """_check_astroquery_available should return True when astroquery is installed."""
        try:
            from astroquery.jplhorizons import Horizons

            has_astroquery = True
        except ImportError:
            has_astroquery = False

        if not has_astroquery:
            pytest.skip("astroquery not installed")

        from libephemeris import spk_auto

        assert spk_auto._check_astroquery_available() is True
