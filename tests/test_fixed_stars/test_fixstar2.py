"""
Unit tests for swe_fixstar2 and swe_fixstar2_ut functions.

Tests the enhanced fixed star lookup functions that support:
- Exact star name lookup
- Hipparcos catalog number lookup
- Partial name search
- Bayer/Flamsteed nomenclature lookup
- Returns the full star name with the position
"""

import pytest
import libephemeris as ephem


@pytest.mark.unit
class TestFixstar2TT:
    """Tests for swe_fixstar2() with Terrestrial Time."""

    def test_fixstar2_exact_name_regulus(self, standard_jd):
        """Test exact name lookup for Regulus."""
        name, pos, retflag, err = ephem.swe_fixstar2("Regulus", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.2f} out of range"
        assert -1 < pos[1] < 2, f"Regulus lat: {pos[1]:.2f} out of range"
        assert pos[2] > 1000, "Fixed stars should be very distant"

    def test_fixstar2_exact_name_spica(self, standard_jd):
        """Test exact name lookup for Spica."""
        name, pos, retflag, err = ephem.swe_fixstar2("Spica", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"
        assert 203 < pos[0] < 205, f"Spica lon: {pos[0]:.2f} out of range"
        assert -3 < pos[1] < -1, f"Spica lat: {pos[1]:.2f} out of range"

    def test_fixstar2_case_insensitive(self, standard_jd):
        """Test that star name lookup is case insensitive."""
        name_upper, pos_upper, _, _ = ephem.swe_fixstar2("REGULUS", standard_jd, 0)
        name_lower, pos_lower, _, _ = ephem.swe_fixstar2("regulus", standard_jd, 0)
        name_mixed, pos_mixed, _, _ = ephem.swe_fixstar2("ReGuLuS", standard_jd, 0)

        assert name_upper == name_lower == name_mixed == "Regulus,alLeo"
        assert pos_upper[0] == pos_lower[0] == pos_mixed[0]

    def test_fixstar2_hip_number(self, standard_jd):
        """Test lookup by Hipparcos catalog number."""
        # Regulus is HIP 49669
        name, pos, retflag, err = ephem.swe_fixstar2("49669", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.2f} out of range"

    def test_fixstar2_hip_number_with_comma(self, standard_jd):
        """Test lookup by HIP number with leading comma (Swiss Ephemeris format)."""
        # Spica is HIP 65474
        name, pos, retflag, err = ephem.swe_fixstar2(",65474", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"
        assert 203 < pos[0] < 205, f"Spica lon: {pos[0]:.2f} out of range"

    def test_fixstar2_nomenclature(self, standard_jd):
        """Test lookup by Bayer/Flamsteed nomenclature."""
        name, pos, retflag, err = ephem.swe_fixstar2("alLeo", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

    def test_fixstar2_nomenclature_case_insensitive(self, standard_jd):
        """Test nomenclature lookup is case insensitive."""
        name1, _, _, err1 = ephem.swe_fixstar2("alVir", standard_jd, 0)
        name2, _, _, err2 = ephem.swe_fixstar2("ALVIR", standard_jd, 0)
        name3, _, _, err3 = ephem.swe_fixstar2("AlViR", standard_jd, 0)

        assert err1 == "" and err2 == "" and err3 == ""
        assert name1 == name2 == name3 == "Spica,alVir"

    def test_fixstar2_partial_name(self, standard_jd):
        """Test partial name lookup (prefix search)."""
        name, pos, retflag, err = ephem.swe_fixstar2("Reg", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

    def test_fixstar2_partial_name_spica(self, standard_jd):
        """Test partial name lookup for Spica."""
        name, pos, retflag, err = ephem.swe_fixstar2("Spi", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"

    def test_fixstar2_unknown_star(self, standard_jd):
        """Test handling of unknown star name."""
        name, pos, retflag, err = ephem.swe_fixstar2("UnknownStar", standard_jd, 0)

        assert "not found" in err, "Should report star not found"
        assert name == "", "Name should be empty on error"
        assert pos == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), "Position should be zero tuple"

    def test_fixstar2_unknown_hip_number(self, standard_jd):
        """Test handling of unknown HIP number."""
        name, pos, retflag, err = ephem.swe_fixstar2("99999", standard_jd, 0)

        assert "not found" in err, "Should report star not found"
        assert name == "", "Name should be empty on error"

    def test_fixstar2_empty_string(self, standard_jd):
        """Test handling of empty star name."""
        name, pos, retflag, err = ephem.swe_fixstar2("", standard_jd, 0)

        assert "Empty" in err, "Should report empty name"
        assert name == "", "Name should be empty on error"

    def test_fixstar2_whitespace_handling(self, standard_jd):
        """Test star name with extra whitespace."""
        name1, pos1, _, _ = ephem.swe_fixstar2("Regulus", standard_jd, 0)
        name2, pos2, _, _ = ephem.swe_fixstar2("  Regulus  ", standard_jd, 0)

        assert name1 == name2, "Should handle whitespace in name"
        assert pos1[0] == pos2[0], "Positions should match"

    def test_fixstar2_comma_format_input(self, standard_jd):
        """Test input with comma (catalog format)."""
        name, pos, retflag, err = ephem.swe_fixstar2("Regulus,alLeo", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

    def test_fixstar2_return_structure(self, standard_jd):
        """Test the return structure of swe_fixstar2."""
        result = ephem.swe_fixstar2("Regulus", standard_jd, 0)

        # Should return a 4-tuple (name, pos, iflag, error)
        assert len(result) == 4, "Should return 4-tuple (name, pos, iflag, error)"

        name, pos, iflag, err = result

        # Name should be string
        assert isinstance(name, str), "Name should be string"
        assert name != "", "Name should not be empty on success"

        # Position should be 6-tuple
        assert len(pos) == 6, "Position should be 6-tuple"
        assert all(isinstance(p, float) for p in pos), "All positions should be floats"

        # iflag should be int
        assert isinstance(iflag, int), "iflag should be int"

        # err should be string
        assert isinstance(err, str), "error should be string"
        assert err == "", "Error should be empty on success"

    def test_fixstar2_vs_fixstar_consistency(self, standard_jd):
        """Test that fixstar2 and fixstar give same position for same star."""
        name2, pos2, _, err2 = ephem.swe_fixstar2("Regulus", standard_jd, 0)
        pos1, _, err1 = ephem.swe_fixstar("Regulus", standard_jd, 0)

        assert err1 == "" and err2 == ""
        assert pos1[0] == pos2[0], "Longitude should match"
        assert pos1[1] == pos2[1], "Latitude should match"
        assert pos1[2] == pos2[2], "Distance should match"

    def test_fixstar2_alias(self, standard_jd):
        """Test that fixstar2 is an alias for swe_fixstar2."""
        name1, pos1, flag1, err1 = ephem.swe_fixstar2("Regulus", standard_jd, 0)
        name2, pos2, flag2, err2 = ephem.fixstar2("Regulus", standard_jd, 0)

        assert name1 == name2
        assert pos1 == pos2
        assert flag1 == flag2
        assert err1 == err2


@pytest.mark.unit
class TestFixstar2UT:
    """Tests for swe_fixstar2_ut() with Universal Time."""

    def test_fixstar2_ut_basic(self, standard_jd):
        """Test basic fixstar2_ut functionality."""
        name, pos, retflag, err = ephem.swe_fixstar2_ut("Regulus", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.2f} out of range"

    def test_fixstar2_ut_hip_number(self, standard_jd):
        """Test fixstar2_ut with HIP number lookup."""
        name, pos, retflag, err = ephem.swe_fixstar2_ut("65474", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Spica,alVir", f"Expected 'Spica,alVir', got '{name}'"
        assert 203 < pos[0] < 205, f"Spica lon: {pos[0]:.2f} out of range"

    def test_fixstar2_ut_partial_name(self, standard_jd):
        """Test fixstar2_ut with partial name lookup."""
        name, pos, retflag, err = ephem.swe_fixstar2_ut("Spi", standard_jd, 0)

        assert err == "", f"Unexpected error: {err}"
        assert name == "Spica,alVir"

    def test_fixstar2_ut_unknown_star(self, standard_jd):
        """Test fixstar2_ut handling of unknown star."""
        name, pos, retflag, err = ephem.swe_fixstar2_ut("UnknownStar", standard_jd, 0)

        assert "not found" in err
        assert name == ""
        assert pos == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    def test_fixstar2_ut_alias(self, standard_jd):
        """Test that fixstar2_ut is an alias for swe_fixstar2_ut."""
        name1, pos1, flag1, err1 = ephem.swe_fixstar2_ut("Spica", standard_jd, 0)
        name2, pos2, flag2, err2 = ephem.fixstar2_ut("Spica", standard_jd, 0)

        assert name1 == name2
        assert pos1 == pos2
        assert flag1 == flag2
        assert err1 == err2

    def test_fixstar2_ut_vs_fixstar_ut_consistency(self, standard_jd):
        """Test that fixstar2_ut and fixstar_ut give same position for same star."""
        name2, pos2, _, err2 = ephem.swe_fixstar2_ut("Spica", standard_jd, 0)
        pos1, _, err1 = ephem.swe_fixstar_ut("Spica", standard_jd, 0)

        assert err1 == "" and err2 == ""
        assert pos1[0] == pos2[0], "Longitude should match"
        assert pos1[1] == pos2[1], "Latitude should match"

    def test_fixstar2_ut_vs_tt_relationship(self, standard_jd):
        """
        Test the relationship between fixstar2 (TT) and fixstar2_ut (UT).

        When called with the same JD value, they should give slightly different
        results because of Delta T.
        """
        jd = standard_jd

        name_tt, pos_tt, _, _ = ephem.swe_fixstar2("Regulus", jd, 0)
        name_ut, pos_ut, _, _ = ephem.swe_fixstar2_ut("Regulus", jd, 0)

        # Names should be the same
        assert name_tt == name_ut == "Regulus,alLeo"

        # Positions should be slightly different due to Delta T
        diff = abs(pos_tt[0] - pos_ut[0])
        assert diff < 0.01, f"Difference too large: {diff}"


@pytest.mark.unit
class TestFixstar2AndFixstar2UtConsistency:
    """Tests for consistency between fixstar2 and fixstar2_ut."""

    def test_deltat_effect(self):
        """
        Test that the Delta T effect is correctly applied.

        If we calculate a JD in UT, then convert to TT by adding Delta T,
        calling fixstar2(jd_tt) should give the same result as fixstar2_ut(jd_ut).
        """
        jd_ut = ephem.swe_julday(2000, 1, 1, 12.0)

        # Get Delta T for this date
        delta_t = ephem.swe_deltat(jd_ut)  # Returns Delta T in days
        jd_tt = jd_ut + delta_t

        # Calculate positions
        name_ut, pos_from_ut, _, _ = ephem.swe_fixstar2_ut("Regulus", jd_ut, 0)
        name_tt, pos_from_tt, _, _ = ephem.swe_fixstar2("Regulus", jd_tt, 0)

        # Names should match
        assert name_ut == name_tt == "Regulus,alLeo"

        # Positions should be essentially identical (allowing for floating point)
        diff = abs(pos_from_ut[0] - pos_from_tt[0])
        assert diff < 1e-6, f"Positions should match: diff={diff}"

    def test_both_functions_return_names(self, standard_jd):
        """Test that both functions return proper star names."""
        name_tt, pos_tt, _, err_tt = ephem.swe_fixstar2("Spica", standard_jd, 0)
        name_ut, pos_ut, _, err_ut = ephem.swe_fixstar2_ut("Spica", standard_jd, 0)

        assert err_tt == "" and err_ut == ""
        assert name_tt == name_ut == "Spica,alVir"
        assert 203 < pos_tt[0] < 205
        assert 203 < pos_ut[0] < 205

    def test_multiple_dates_proper_motion(self):
        """Test that positions change with proper motion over time."""
        jd_2000 = ephem.swe_julday(2000, 1, 1, 12.0)
        jd_2050 = ephem.swe_julday(2050, 1, 1, 12.0)

        name_2000, pos_2000, _, err1 = ephem.swe_fixstar2("Regulus", jd_2000, 0)
        name_2050, pos_2050, _, err2 = ephem.swe_fixstar2("Regulus", jd_2050, 0)

        assert err1 == "" and err2 == ""
        assert name_2000 == name_2050 == "Regulus,alLeo"

        # Star should have moved due to proper motion + precession
        diff = abs(pos_2050[0] - pos_2000[0])
        assert 0.001 < diff < 5.0, f"Proper motion over 50 years: {diff}"
