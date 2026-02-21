"""
Tests for TRUE_LILITH_METHODS.md documentation consistency.

Verifies that the precision claims in TRUE_LILITH_METHODS.md match the
current precision values documented in PRECISION.md, ensuring users
receive accurate information for method selection decisions.
"""

import pathlib
import re

import pytest


class TestTrueLilithMethodsDocumentation:
    """Verify TRUE_LILITH_METHODS.md documentation is accurate and consistent."""

    @pytest.fixture
    def true_lilith_methods_doc(self) -> str:
        """Load TRUE_LILITH_METHODS.md content."""
        doc_path = (
            pathlib.Path(__file__).parent.parent.parent
            / "docs"
            / "TRUE_LILITH_METHODS.md"
        )
        assert doc_path.exists(), "TRUE_LILITH_METHODS.md not found"
        return doc_path.read_text()

    @pytest.fixture
    def precision_doc(self) -> str:
        """Load PRECISION.md content."""
        doc_path = pathlib.Path(__file__).parent.parent.parent / "docs" / "PRECISION.md"
        assert doc_path.exists(), "PRECISION.md not found"
        return doc_path.read_text()

    @pytest.mark.precision
    def test_true_lilith_precision_values_present(self, true_lilith_methods_doc):
        """Verify TRUE_LILITH_METHODS.md contains current precision measurements."""
        content = true_lilith_methods_doc

        # Check for documented precision values (~52 arcsec mean, ~235 arcsec max)
        assert "52 arcsec" in content or "~52" in content, (
            "Missing True Lilith mean difference (~52 arcsec)"
        )
        assert "235 arcsec" in content or "~235" in content, (
            "Missing True Lilith max difference (~235 arcsec)"
        )
        assert "0.015" in content, (
            "Missing True Lilith mean difference in degrees (~0.015)"
        )
        assert "0.065" in content, (
            "Missing True Lilith max difference in degrees (~0.065)"
        )

    @pytest.mark.precision
    def test_comparison_table_exists(self, true_lilith_methods_doc):
        """Verify comparison table with True Lilith vs Mean Lilith vs Interpolated exists."""
        content = true_lilith_methods_doc

        # Check for comparison table content
        assert "True Lilith" in content and "Mean Lilith" in content, (
            "Missing True Lilith vs Mean Lilith comparison"
        )
        assert "Interpolated" in content, "Missing Interpolated Apogee comparison"
        assert "1.1" in content or "1.1°" in content, (
            "Missing Interpolated Apogee precision (~1.1 degrees)"
        )

    @pytest.mark.precision
    def test_method_selection_guide_exists(self, true_lilith_methods_doc):
        """Verify method selection guide exists with recommendations."""
        content = true_lilith_methods_doc

        # Check for selection guide section
        assert "Recommendation" in content or "recommended" in content.lower(), (
            "Missing method recommendation section"
        )
        assert "osculating" in content.lower(), (
            "Missing osculating apogee use case description"
        )

    @pytest.mark.precision
    def test_sub_arcminute_precision_documented(self, true_lilith_methods_doc):
        """Verify sub-arcminute precision claim is documented."""
        content = true_lilith_methods_doc

        # Check for sub-arcminute precision claim
        assert (
            "sub-arcminute" in content.lower() or "sub arcminute" in content.lower()
        ), "Missing sub-arcminute precision claim for True Lilith"

    @pytest.mark.precision
    def test_precision_consistency_with_precision_md(
        self, true_lilith_methods_doc, precision_doc
    ):
        """Verify precision values are consistent between TRUE_LILITH_METHODS.md and PRECISION.md."""
        # Extract True Lilith precision from PRECISION.md
        precision_content = precision_doc

        # Check that both documents mention the same precision values
        # Mean: ~52 arcsec
        assert "52 arcsec" in precision_content, (
            "PRECISION.md missing 52 arcsec value for True Lilith"
        )
        assert (
            "52 arcsec" in true_lilith_methods_doc or "~52" in true_lilith_methods_doc
        ), "TRUE_LILITH_METHODS.md missing 52 arcsec value for True Lilith"

        # Max: ~235 arcsec
        assert "235 arcsec" in precision_content, (
            "PRECISION.md missing 235 arcsec value for True Lilith"
        )
        assert (
            "235 arcsec" in true_lilith_methods_doc or "~235" in true_lilith_methods_doc
        ), "TRUE_LILITH_METHODS.md missing 235 arcsec value for True Lilith"

        # Degree values
        assert "0.015" in precision_content and "0.015" in true_lilith_methods_doc, (
            "Inconsistent mean precision in degrees (~0.015)"
        )
        assert "0.065" in precision_content and "0.065" in true_lilith_methods_doc, (
            "Inconsistent max precision in degrees (~0.065)"
        )

    @pytest.mark.precision
    def test_mean_lilith_precision_documented(self, true_lilith_methods_doc):
        """Verify Mean Lilith precision is documented for comparison."""
        content = true_lilith_methods_doc

        assert "12 arcsec" in content or "0.003" in content, (
            "Missing Mean Lilith precision (~12 arcsec or ~0.003 deg)"
        )

    @pytest.mark.precision
    def test_eccentricity_vector_method_documented(self, true_lilith_methods_doc):
        """Verify eccentricity vector method is documented."""
        content = true_lilith_methods_doc

        assert "eccentricity vector" in content.lower(), (
            "Missing eccentricity vector method description"
        )
        assert "calc_true_lilith" in content, (
            "Missing calc_true_lilith function reference"
        )

    @pytest.mark.precision
    def test_500_date_sampling_mentioned(self, true_lilith_methods_doc):
        """Verify 500-date random sampling methodology is mentioned."""
        content = true_lilith_methods_doc

        assert "500" in content and (
            "date" in content.lower() or "sample" in content.lower()
        ), "Missing 500-date random sampling methodology description"

    @pytest.mark.precision
    def test_pyswisseph_compatibility_documented(self, true_lilith_methods_doc):
        """Verify pyswisseph compatibility is documented."""
        content = true_lilith_methods_doc

        assert "pyswisseph" in content or "SE_OSCU_APOG" in content, (
            "Missing pyswisseph compatibility documentation"
        )
        assert "SE_OSCU_APOG" in content, (
            "Missing SE_OSCU_APOG constant reference for compatibility"
        )
