"""
Unit tests for swe_fixstar_mag and swe_fixstar2_mag functions.

Tests the magnitude lookup functions that return only the visual magnitude
of a fixed star without calculating position.
"""

import pytest
import libephemeris as ephem


@pytest.mark.unit
class TestFixstarMag:
    """Tests for swe_fixstar_mag()."""

    def test_fixstar_mag_regulus(self):
        """Test magnitude lookup for Regulus."""
        mag, err = ephem.swe_fixstar_mag("Regulus")

        assert err == "", f"Unexpected error: {err}"
        assert mag == pytest.approx(1.40, abs=0.01), f"Regulus magnitude: {mag}"

    def test_fixstar_mag_spica(self):
        """Test magnitude lookup for Spica."""
        mag, err = ephem.swe_fixstar_mag("Spica")

        assert err == "", f"Unexpected error: {err}"
        assert mag == pytest.approx(1.04, abs=0.01), f"Spica magnitude: {mag}"

    def test_fixstar_mag_case_insensitive(self):
        """Test that star name lookup is case insensitive."""
        mag_upper, err1 = ephem.swe_fixstar_mag("REGULUS")
        mag_lower, err2 = ephem.swe_fixstar_mag("regulus")
        mag_mixed, err3 = ephem.swe_fixstar_mag("ReGuLuS")

        assert err1 == "" and err2 == "" and err3 == ""
        assert mag_upper == mag_lower == mag_mixed

    def test_fixstar_mag_unknown_star(self):
        """Test handling of unknown star name."""
        mag, err = ephem.swe_fixstar_mag("UnknownStar")

        assert "not found" in err, "Should report star not found"
        assert mag == 0.0, "Magnitude should be 0.0 on error"

    def test_fixstar_mag_with_comma_in_name(self):
        """Test star name with comma (common in catalogs)."""
        mag, err = ephem.swe_fixstar_mag("Regulus,alLeo")

        assert err == "", "Should handle comma in star name"
        assert mag == pytest.approx(1.40, abs=0.01)

    def test_fixstar_mag_whitespace_in_name(self):
        """Test star name with extra whitespace."""
        mag1, _ = ephem.swe_fixstar_mag("Regulus")
        mag2, _ = ephem.swe_fixstar_mag("  Regulus  ")

        assert mag1 == mag2, "Should handle whitespace in name"

    def test_fixstar_mag_return_structure(self):
        """Test the return structure of swe_fixstar_mag."""
        result = ephem.swe_fixstar_mag("Regulus")

        # Should return a 2-tuple
        assert len(result) == 2, "Should return 2-tuple (magnitude, error)"

        mag, err = result

        assert isinstance(mag, float), "Magnitude should be float"
        assert isinstance(err, str), "Error should be string"

    def test_fixstar_mag_alias(self):
        """Test that fixstar_mag is an alias for swe_fixstar_mag."""
        mag1, err1 = ephem.swe_fixstar_mag("Regulus")
        mag2, err2 = ephem.fixstar_mag("Regulus")

        assert mag1 == mag2
        assert err1 == err2


@pytest.mark.unit
class TestFixstar2Mag:
    """Tests for swe_fixstar2_mag()."""

    def test_fixstar2_mag_exact_name_regulus(self):
        """Test exact name lookup for Regulus magnitude."""
        name, mag, err = ephem.swe_fixstar2_mag("Regulus")

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert mag == pytest.approx(1.40, abs=0.01), f"Regulus magnitude: {mag}"

    def test_fixstar2_mag_exact_name_spica(self):
        """Test exact name lookup for Spica magnitude."""
        name, mag, err = ephem.swe_fixstar2_mag("Spica")

        assert err == "", f"Unexpected error: {err}"
        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"
        assert mag == pytest.approx(1.04, abs=0.01), f"Spica magnitude: {mag}"

    def test_fixstar2_mag_case_insensitive(self):
        """Test that star name lookup is case insensitive."""
        name_upper, mag_upper, _ = ephem.swe_fixstar2_mag("REGULUS")
        name_lower, mag_lower, _ = ephem.swe_fixstar2_mag("regulus")
        name_mixed, mag_mixed, _ = ephem.swe_fixstar2_mag("ReGuLuS")

        assert name_upper == name_lower == name_mixed == "Regulus,alLeo"
        assert mag_upper == mag_lower == mag_mixed

    def test_fixstar2_mag_hip_number(self):
        """Test lookup by Hipparcos catalog number."""
        name, mag, err = ephem.swe_fixstar2_mag("49669")

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert mag == pytest.approx(1.40, abs=0.01)

    def test_fixstar2_mag_hip_number_with_comma(self):
        """Test lookup by HIP number with leading comma (Swiss Ephemeris format)."""
        name, mag, err = ephem.swe_fixstar2_mag(",65474")

        assert err == "", f"Unexpected error: {err}"
        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"
        assert mag == pytest.approx(1.04, abs=0.01)

    def test_fixstar2_mag_nomenclature(self):
        """Test lookup by Bayer/Flamsteed nomenclature."""
        name, mag, err = ephem.swe_fixstar2_mag("alLeo")

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert mag == pytest.approx(1.40, abs=0.01)

    def test_fixstar2_mag_nomenclature_case_insensitive(self):
        """Test nomenclature lookup is case insensitive."""
        name1, mag1, err1 = ephem.swe_fixstar2_mag("alVir")
        name2, mag2, err2 = ephem.swe_fixstar2_mag("ALVIR")
        name3, mag3, err3 = ephem.swe_fixstar2_mag("AlViR")

        assert err1 == "" and err2 == "" and err3 == ""
        assert name1 == name2 == name3 == "Spica,alVir"
        assert mag1 == mag2 == mag3

    def test_fixstar2_mag_partial_name(self):
        """Test partial name lookup (prefix search)."""
        name, mag, err = ephem.swe_fixstar2_mag("Reg")

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

    def test_fixstar2_mag_partial_name_spica(self):
        """Test partial name lookup for Spica."""
        name, mag, err = ephem.swe_fixstar2_mag("Spi")

        assert err == "", f"Unexpected error: {err}"
        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"

    def test_fixstar2_mag_unknown_star(self):
        """Test handling of unknown star name."""
        name, mag, err = ephem.swe_fixstar2_mag("UnknownStar")

        assert "not found" in err, "Should report star not found"
        assert name == "", "Name should be empty on error"
        assert mag == 0.0, "Magnitude should be 0.0 on error"

    def test_fixstar2_mag_unknown_hip_number(self):
        """Test handling of unknown HIP number."""
        name, mag, err = ephem.swe_fixstar2_mag("99999")

        assert "not found" in err, "Should report star not found"
        assert name == "", "Name should be empty on error"

    def test_fixstar2_mag_empty_string(self):
        """Test handling of empty star name."""
        name, mag, err = ephem.swe_fixstar2_mag("")

        assert "Empty" in err, "Should report empty name"
        assert name == "", "Name should be empty on error"

    def test_fixstar2_mag_whitespace_handling(self):
        """Test star name with extra whitespace."""
        name1, mag1, _ = ephem.swe_fixstar2_mag("Regulus")
        name2, mag2, _ = ephem.swe_fixstar2_mag("  Regulus  ")

        assert name1 == name2, "Should handle whitespace in name"
        assert mag1 == mag2, "Magnitudes should match"

    def test_fixstar2_mag_comma_format_input(self):
        """Test input with comma (catalog format)."""
        name, mag, err = ephem.swe_fixstar2_mag("Regulus,alLeo")

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

    def test_fixstar2_mag_return_structure(self):
        """Test the return structure of swe_fixstar2_mag."""
        result = ephem.swe_fixstar2_mag("Regulus")

        # Should return a 3-tuple (name, magnitude, error)
        assert len(result) == 3, "Should return 3-tuple (name, magnitude, error)"

        name, mag, err = result

        # Name should be string
        assert isinstance(name, str), "Name should be string"
        assert name != "", "Name should not be empty on success"

        # Magnitude should be float
        assert isinstance(mag, float), "Magnitude should be float"

        # err should be string
        assert isinstance(err, str), "error should be string"
        assert err == "", "Error should be empty on success"

    def test_fixstar2_mag_vs_fixstar_mag_consistency(self):
        """Test that fixstar2_mag and fixstar_mag give same magnitude for same star."""
        name, mag2, err2 = ephem.swe_fixstar2_mag("Regulus")
        mag1, err1 = ephem.swe_fixstar_mag("Regulus")

        assert err1 == "" and err2 == ""
        assert mag1 == mag2, "Magnitude should match"

    def test_fixstar2_mag_alias(self):
        """Test that fixstar2_mag is an alias for swe_fixstar2_mag."""
        name1, mag1, err1 = ephem.swe_fixstar2_mag("Regulus")
        name2, mag2, err2 = ephem.fixstar2_mag("Regulus")

        assert name1 == name2
        assert mag1 == mag2
        assert err1 == err2


@pytest.mark.unit
class TestMagnitudeValues:
    """Tests to validate magnitude values are astronomically correct."""

    def test_regulus_magnitude_value(self):
        """Regulus is a 1st magnitude star (V=1.40)."""
        mag, _ = ephem.swe_fixstar_mag("Regulus")
        # Regulus V magnitude is approximately 1.40
        assert 1.35 < mag < 1.45, f"Regulus magnitude {mag} not in expected range"

    def test_spica_magnitude_value(self):
        """Spica is a 1st magnitude star (V=1.04)."""
        mag, _ = ephem.swe_fixstar_mag("Spica")
        # Spica V magnitude is approximately 1.04
        assert 0.99 < mag < 1.09, f"Spica magnitude {mag} not in expected range"

    def test_spica_brighter_than_regulus(self):
        """Spica (V=1.04) should be brighter than Regulus (V=1.40)."""
        mag_spica, _ = ephem.swe_fixstar_mag("Spica")
        mag_regulus, _ = ephem.swe_fixstar_mag("Regulus")

        # Lower magnitude means brighter
        assert mag_spica < mag_regulus, "Spica should be brighter than Regulus"
