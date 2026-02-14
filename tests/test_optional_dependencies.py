"""
Tests for dependency documentation and availability checking.

This module tests that:
1. Required and optional dependencies are correctly documented in README.md
2. Optional dependencies are defined in pyproject.toml extras
3. Required dependencies (pyerfa, astroquery) are available and functional
"""

from pathlib import Path

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


class TestAstroqueryRequiredImport:
    """Tests for astroquery required dependency (promoted from optional)."""

    def test_spk_auto_module_imports(self):
        """spk_auto module should import successfully."""
        from libephemeris import spk_auto

        assert spk_auto is not None

    def test_astroquery_availability_check_returns_true(self):
        """astroquery availability check should always return True (required dep)."""
        from libephemeris import spk_auto

        assert hasattr(spk_auto, "_check_astroquery_available")
        result = spk_auto._check_astroquery_available()
        assert result is True

    def test_astroquery_import_works(self):
        """astroquery should be importable (required dependency)."""
        from astroquery.jplhorizons import Horizons  # noqa: F401


class TestDependencyFeatureIntegration:
    """Tests that required dependencies are properly integrated."""

    def test_spk_auto_functions_exist(self):
        """Key spk_auto functions should be available."""
        from libephemeris import spk_auto

        assert hasattr(spk_auto, "enable_auto_spk")
        assert hasattr(spk_auto, "auto_get_spk")
        assert hasattr(spk_auto, "download_spk_from_horizons")


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
