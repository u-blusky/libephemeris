"""Tests for package metadata and distribution.

This module verifies that the package is correctly configured for PyPI distribution.
"""

import importlib.metadata
import pathlib

import pytest

import libephemeris


class TestPackageMetadata:
    """Test package metadata from pyproject.toml."""

    def test_version_is_string(self):
        """Version should be a string."""
        assert isinstance(libephemeris.__version__, str)

    def test_version_format(self):
        """Version should follow semantic versioning pattern."""
        version = libephemeris.__version__
        parts = version.split(".")
        assert len(parts) >= 2, "Version should have at least major.minor"
        # Check major and minor are numeric
        assert parts[0].isdigit(), "Major version should be numeric"
        assert parts[1].isdigit(), "Minor version should be numeric"

    def test_version_matches_pyproject(self):
        """Package __version__ should match pyproject.toml."""
        metadata = importlib.metadata.metadata("libephemeris")
        assert libephemeris.__version__ == metadata["Version"]

    def test_author_is_set(self):
        """Author should be defined."""
        assert hasattr(libephemeris, "__author__")
        assert libephemeris.__author__ == "Giacomo Battaglia"

    def test_license_is_set(self):
        """License should be defined."""
        assert hasattr(libephemeris, "__license__")
        assert libephemeris.__license__ == "LGPL-3.0"


class TestPackageDistribution:
    """Test package distribution configuration."""

    def test_package_is_installed(self):
        """Package should be installed and discoverable."""
        metadata = importlib.metadata.metadata("libephemeris")
        assert metadata is not None

    def test_package_name(self):
        """Package name should be libephemeris."""
        metadata = importlib.metadata.metadata("libephemeris")
        assert metadata["Name"] == "libephemeris"

    def test_package_description(self):
        """Package should have a description."""
        metadata = importlib.metadata.metadata("libephemeris")
        summary = metadata["Summary"]
        assert summary is not None
        assert len(summary) > 0
        assert "ephemeris" in summary.lower() or "astronomical" in summary.lower()

    def test_package_requires_python(self):
        """Package should specify Python version requirement."""
        metadata = importlib.metadata.metadata("libephemeris")
        requires_python = metadata.get("Requires-Python")
        assert requires_python is not None
        assert "3.10" in requires_python or ">=" in requires_python

    def test_package_has_dependencies(self):
        """Package should have required dependencies."""
        requires = importlib.metadata.requires("libephemeris")
        assert requires is not None
        # Filter out optional dependencies (those with ; extra ==)
        core_deps = [r for r in requires if "extra" not in r]
        dep_names = [r.split()[0].lower() for r in core_deps]
        assert "skyfield" in dep_names or any("skyfield" in d for d in dep_names)

    def test_package_homepage_url(self):
        """Package should have homepage URL."""
        metadata = importlib.metadata.metadata("libephemeris")
        # Check for Homepage or Project-URL with Homepage
        project_urls = metadata.get_all("Project-URL") or []
        has_homepage = any("Homepage" in url for url in project_urls)
        has_home_page = metadata.get("Home-page") is not None
        assert has_homepage or has_home_page, "Package should have a homepage URL"


class TestPackageFiles:
    """Test that required package files exist."""

    @pytest.fixture
    def project_root(self):
        """Get the project root directory."""
        return pathlib.Path(__file__).parent.parent

    def test_readme_exists(self, project_root):
        """README.md should exist in project root."""
        readme = project_root / "README.md"
        assert readme.exists(), "README.md not found"
        assert readme.stat().st_size > 0, "README.md is empty"

    def test_license_exists(self, project_root):
        """LICENSE file should exist in project root."""
        license_file = project_root / "LICENSE"
        assert license_file.exists(), "LICENSE not found"
        content = license_file.read_text()
        assert "GNU LESSER GENERAL PUBLIC LICENSE" in content

    def test_pyproject_toml_exists(self, project_root):
        """pyproject.toml should exist in project root."""
        pyproject = project_root / "pyproject.toml"
        assert pyproject.exists(), "pyproject.toml not found"

    def test_py_typed_marker_exists(self, project_root):
        """py.typed marker should exist for type hints."""
        py_typed = project_root / "libephemeris" / "py.typed"
        assert py_typed.exists(), "py.typed marker not found"


class TestPackageImports:
    """Test that public API is correctly exported."""

    def test_core_functions_importable(self):
        """Core pyswisseph-compatible functions should be importable."""
        # Time functions
        assert hasattr(libephemeris, "swe_julday")
        assert hasattr(libephemeris, "julday")
        assert hasattr(libephemeris, "swe_revjul")
        assert hasattr(libephemeris, "revjul")

        # Calculation functions
        assert hasattr(libephemeris, "swe_calc_ut")
        assert hasattr(libephemeris, "calc_ut")
        assert hasattr(libephemeris, "swe_calc")
        assert hasattr(libephemeris, "calc")

        # House functions
        assert hasattr(libephemeris, "swe_houses")
        assert hasattr(libephemeris, "houses")

        # Ayanamsa functions
        assert hasattr(libephemeris, "swe_get_ayanamsa_ut")
        assert hasattr(libephemeris, "get_ayanamsa_ut")
        assert hasattr(libephemeris, "swe_set_sid_mode")
        assert hasattr(libephemeris, "set_sid_mode")

    def test_context_api_importable(self):
        """Thread-safe EphemerisContext should be importable."""
        assert hasattr(libephemeris, "EphemerisContext")
        from libephemeris import EphemerisContext

        assert EphemerisContext is not None

    def test_constants_importable(self):
        """Constants should be importable from libephemeris.constants."""
        from libephemeris.constants import SE_SUN, SE_MOON, SEFLG_SWIEPH, SEFLG_SPEED

        assert SE_SUN == 0
        assert SE_MOON == 1
        assert isinstance(SEFLG_SWIEPH, int)
        assert isinstance(SEFLG_SPEED, int)

    def test_exception_importable(self):
        """Error exception should be importable."""
        assert hasattr(libephemeris, "Error")
        from libephemeris import Error

        assert issubclass(Error, Exception)


class TestPackageAllExports:
    """Test __all__ exports."""

    def test_all_is_defined(self):
        """__all__ should be defined."""
        assert hasattr(libephemeris, "__all__")
        assert isinstance(libephemeris.__all__, list)
        assert len(libephemeris.__all__) > 0

    def test_all_items_exist(self):
        """All items in __all__ should exist in the module."""
        for name in libephemeris.__all__:
            assert hasattr(libephemeris, name), (
                f"{name} in __all__ but not in libephemeris"
            )
