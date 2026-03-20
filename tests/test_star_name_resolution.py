"""
Unit tests for star name resolution with aliases and wildcards.

Tests the pyswisseph-compatible star name resolution system including:
- STAR_ALIASES dictionary with alternative names
- resolve_star_name function with prefix/fuzzy matching
- swe_fixstar_ut with canonical star name return
"""

import pytest
import libephemeris as ephem
from libephemeris.exceptions import Error
from libephemeris.fixed_stars import (
    STAR_ALIASES,
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
    swe_fixstar_ut,
    swe_fixstar,
)
from libephemeris.constants import (
    SE_REGULUS,
    SE_SPICA_STAR,
    SE_ALGOL,
    SE_SIRIUS,
    SE_ALDEBARAN,
    SE_ANTARES,
    SE_VEGA,
    SE_POLARIS,
    SE_FOMALHAUT,
    SE_BETELGEUSE,
    SE_RIGEL,
    SE_PROCYON,
    SE_ARCTURUS,
    SE_DENEB,
    SE_POLLUX,
    SE_CASTOR,
    SE_ALTAIR,
)


@pytest.mark.unit
class TestStarAliasesDictionary:
    """Tests for the STAR_ALIASES dictionary structure."""

    def test_star_aliases_has_minimum_entries(self):
        """Verify STAR_ALIASES has at least 300 entries as required."""
        assert len(STAR_ALIASES) >= 300, (
            f"STAR_ALIASES should have at least 300 entries, "
            f"but has {len(STAR_ALIASES)}"
        )

    def test_star_aliases_values_are_valid_star_ids(self):
        """Verify all alias values map to valid star IDs in catalog."""
        valid_star_ids = {entry.id for entry in STAR_CATALOG}
        for alias, star_id in STAR_ALIASES.items():
            assert star_id in valid_star_ids, (
                f"Alias '{alias}' maps to invalid star ID {star_id}"
            )

    def test_star_aliases_keys_are_normalized(self):
        """Verify all alias keys are normalized (uppercase ASCII, as-is for Unicode)."""
        for alias in STAR_ALIASES:
            # Check that ASCII characters are uppercase
            ascii_chars = "".join(c for c in alias if c.isascii())
            assert ascii_chars == ascii_chars.upper(), (
                f"ASCII chars in alias '{alias}' should be uppercase"
            )


@pytest.mark.unit
class TestResolveStarName:
    """Tests for the resolve_star_name function."""

    def test_empty_name_returns_none(self):
        """Empty string should return None."""
        assert resolve_star_name("") is None
        assert resolve_star_name("   ") is None

    def test_exact_canonical_name_match(self):
        """Test exact match against canonical names."""
        assert resolve_star_name("Regulus") == SE_REGULUS
        assert resolve_star_name("Spica") == SE_SPICA_STAR
        assert resolve_star_name("Algol") == SE_ALGOL
        assert resolve_star_name("Sirius") == SE_SIRIUS

    def test_case_insensitive_matching(self):
        """Test that matching is case-insensitive."""
        assert resolve_star_name("REGULUS") == SE_REGULUS
        assert resolve_star_name("regulus") == SE_REGULUS
        assert resolve_star_name("ReGuLuS") == SE_REGULUS
        assert resolve_star_name("SIRIUS") == SE_SIRIUS
        assert resolve_star_name("sirius") == SE_SIRIUS

    def test_whitespace_stripping(self):
        """Test that whitespace is stripped."""
        assert resolve_star_name("  Regulus  ") == SE_REGULUS
        assert resolve_star_name("\tSirius\n") == SE_SIRIUS

    def test_alias_matching(self):
        """Test matching against STAR_ALIASES."""
        # Regulus aliases
        assert resolve_star_name("Cor Leonis") == SE_REGULUS
        assert resolve_star_name("Alpha Leo") == SE_REGULUS
        assert resolve_star_name("Alpha Leonis") == SE_REGULUS

        # Algol aliases
        assert resolve_star_name("Demon Star") == SE_ALGOL
        assert resolve_star_name("Beta Persei") == SE_ALGOL

        # Sirius aliases
        assert resolve_star_name("Dog Star") == SE_SIRIUS
        assert resolve_star_name("Alpha CMa") == SE_SIRIUS

    def test_comma_prefix_partial_match(self):
        """Test comma-prefix partial matching (pyswisseph convention)."""
        # ,alg should find Algol
        result = resolve_star_name(",alg")
        assert result == SE_ALGOL, f"Expected Algol, got {result}"

        # ,reg should find Regulus
        result = resolve_star_name(",reg")
        assert result == SE_REGULUS, f"Expected Regulus, got {result}"

        # ,sir should find Sirius
        result = resolve_star_name(",sir")
        assert result == SE_SIRIUS, f"Expected Sirius, got {result}"

        # ,spi should find Spica
        result = resolve_star_name(",spi")
        assert result == SE_SPICA_STAR, f"Expected Spica, got {result}"

    def test_comma_prefix_with_whitespace(self):
        """Test comma-prefix with trailing whitespace."""
        assert resolve_star_name(", alg") == SE_ALGOL
        assert resolve_star_name(",  reg  ") == SE_REGULUS

    def test_comma_prefix_empty_returns_none(self):
        """Test that comma with no prefix returns None."""
        assert resolve_star_name(",") is None
        assert resolve_star_name(",  ") is None

    def test_trailing_comma_format(self):
        """Test pyswisseph format with trailing comma (star_name,nomenclature)."""
        assert resolve_star_name("Regulus,alLeo") == SE_REGULUS
        assert resolve_star_name("Spica,alVir") == SE_SPICA_STAR

    def test_fuzzy_matching_for_short_inputs(self):
        """Test fuzzy matching (substring) for short inputs."""
        # "SIR" should match SIRIUS via substring
        result = resolve_star_name("SIR")
        assert result == SE_SIRIUS, f"Expected Sirius for 'SIR', got {result}"

    def test_unknown_name_returns_none(self):
        """Test that unknown star names return None."""
        assert resolve_star_name("UnknownStar") is None
        assert resolve_star_name("NotARealStar") is None
        assert resolve_star_name("xyz123") is None


@pytest.mark.unit
class TestGetCanonicalStarName:
    """Tests for the get_canonical_star_name function."""

    def test_valid_star_ids_return_names(self):
        """Test that valid star IDs return canonical names."""
        assert get_canonical_star_name(SE_REGULUS) == "Regulus"
        assert get_canonical_star_name(SE_SPICA_STAR) == "Spica"
        assert get_canonical_star_name(SE_ALGOL) == "Algol"
        assert get_canonical_star_name(SE_SIRIUS) == "Sirius"

    def test_invalid_star_id_returns_none(self):
        """Test that invalid star IDs return None."""
        assert get_canonical_star_name(-1) is None
        assert get_canonical_star_name(999999) is None


@pytest.mark.unit
class TestSweFixstarUtWithAliases:
    """Tests for swe_fixstar_ut with star name resolution."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch as Julian Day."""
        return 2451545.0

    def test_exact_name_works(self, standard_jd):
        """Test swe_fixstar_ut with exact canonical names."""
        pos, name, retflag = swe_fixstar_ut("Regulus", standard_jd, 0)
        assert name == "Regulus"
        assert 149 < pos[0] < 151  # Regulus longitude

    def test_returns_canonical_name(self, standard_jd):
        """Test that swe_fixstar_ut returns canonical star name."""
        pos, name, retflag = swe_fixstar_ut("Alpha Leo", standard_jd, 0)
        assert name == "Regulus", f"Expected 'Regulus', got '{name}'"

        pos, name, retflag = swe_fixstar_ut("Dog Star", standard_jd, 0)
        assert name == "Sirius", f"Expected 'Sirius', got '{name}'"

    def test_case_insensitive(self, standard_jd):
        """Test case-insensitive star name lookup."""
        pos1, name1, _ = swe_fixstar_ut("SIRIUS", standard_jd, 0)
        pos2, name2, _ = swe_fixstar_ut("sirius", standard_jd, 0)
        pos3, name3, _ = swe_fixstar_ut("Sirius", standard_jd, 0)

        assert name1 == name2 == name3 == "Sirius"
        assert abs(pos1[0] - pos2[0]) < 0.0001
        assert abs(pos1[0] - pos3[0]) < 0.0001

    def test_comma_prefix_search(self, standard_jd):
        """Test comma-prefix partial search."""
        pos, name, retflag = swe_fixstar_ut(",alg", standard_jd, 0)
        assert name == "Algol", f"Expected 'Algol', got '{name}'"

    def test_bayer_designation(self, standard_jd):
        """Test Bayer designation lookup."""
        pos, name, _ = swe_fixstar_ut("Alpha Leo", standard_jd, 0)
        assert name == "Regulus"

        pos, name, _ = swe_fixstar_ut("Alpha CMa", standard_jd, 0)
        assert name == "Sirius"

    def test_identical_results_for_aliases(self, standard_jd):
        """Test that all aliases for a star return identical positions."""
        # Get Sirius position using different aliases
        pos_sirius, _, _ = swe_fixstar_ut("Sirius", standard_jd, 0)
        pos_dog_star, _, _ = swe_fixstar_ut("Dog Star", standard_jd, 0)
        pos_alpha_cma, _, _ = swe_fixstar_ut("Alpha CMa", standard_jd, 0)

        # All should return the same position
        assert abs(pos_sirius[0] - pos_dog_star[0]) < 0.0001
        assert abs(pos_sirius[0] - pos_alpha_cma[0]) < 0.0001

    def test_unknown_star_returns_error(self, standard_jd):
        """Test that unknown star raises an error."""
        with pytest.raises(Exception):
            swe_fixstar_ut("UnknownStar", standard_jd, 0)

    def test_backward_compatibility_regulus(self, standard_jd):
        """Test backward compatibility - Regulus still works."""
        pos, name, retflag = swe_fixstar_ut("Regulus", standard_jd, 0)
        assert name == "Regulus"
        assert 149 < pos[0] < 151
        assert -1 < pos[1] < 2

    def test_backward_compatibility_spica(self, standard_jd):
        """Test backward compatibility - Spica still works."""
        pos, name, retflag = swe_fixstar_ut("Spica", standard_jd, 0)
        assert name == "Spica"
        assert 203 < pos[0] < 205
        assert -3 < pos[1] < -1


@pytest.mark.unit
class TestSweFixstarWithAliases:
    """Tests for swe_fixstar (TT version) with star name resolution."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch as Julian Day."""
        return 2451545.0

    def test_exact_name_works(self, standard_jd):
        """Test swe_fixstar with exact canonical names."""
        pos, name, retflag = swe_fixstar("Regulus", standard_jd, 0)
        assert name == "Regulus"

    def test_alias_works(self, standard_jd):
        """Test swe_fixstar with aliases."""
        pos, name, retflag = swe_fixstar("Alpha Leo", standard_jd, 0)
        assert name == "Regulus"

    def test_comma_prefix_works(self, standard_jd):
        """Test swe_fixstar with comma-prefix search."""
        pos, name, retflag = swe_fixstar(",sir", standard_jd, 0)
        assert name == "Sirius"


@pytest.mark.unit
class TestAllStarsHaveAliases:
    """Tests to verify all catalog stars have proper aliases."""

    def test_all_catalog_stars_are_resolvable(self):
        """Test that all stars in catalog can be resolved by name."""
        for entry in STAR_CATALOG:
            result = resolve_star_name(entry.name)
            assert result == entry.id, (
                f"Star '{entry.name}' should resolve to {entry.id}, got {result}"
            )

    def test_bayer_designations_work(self):
        """Test various Bayer designation formats."""
        bayer_tests = [
            ("Alpha Leo", SE_REGULUS),
            ("Alpha Leonis", SE_REGULUS),
            ("Alpha Vir", SE_SPICA_STAR),
            ("Alpha Virginis", SE_SPICA_STAR),
            ("Beta Per", SE_ALGOL),
            ("Beta Persei", SE_ALGOL),
            ("Alpha Tau", SE_ALDEBARAN),
            ("Alpha Sco", SE_ANTARES),
            ("Alpha Lyr", SE_VEGA),
            ("Alpha Ori", SE_BETELGEUSE),
            ("Beta Ori", SE_RIGEL),
            ("Alpha CMi", SE_PROCYON),
            ("Alpha Boo", SE_ARCTURUS),
            ("Alpha Cyg", SE_DENEB),
            ("Beta Gem", SE_POLLUX),
            ("Alpha Gem", SE_CASTOR),
            ("Alpha Aql", SE_ALTAIR),
        ]
        for designation, expected_id in bayer_tests:
            result = resolve_star_name(designation)
            assert result == expected_id, (
                f"Bayer designation '{designation}' should resolve to "
                f"{expected_id}, got {result}"
            )

    def test_common_names_work(self):
        """Test various common/traditional names."""
        common_name_tests = [
            ("Dog Star", SE_SIRIUS),
            ("Demon Star", SE_ALGOL),
            ("North Star", SE_POLARIS),
            ("Pole Star", SE_POLARIS),
        ]
        for name, expected_id in common_name_tests:
            result = resolve_star_name(name)
            assert result == expected_id, (
                f"Common name '{name}' should resolve to {expected_id}, got {result}"
            )


@pytest.mark.integration
class TestPyswissephCompatibility:
    """Integration tests comparing against pyswisseph behavior."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch as Julian Day."""
        return 2451545.0

    def test_comma_prefix_finds_algol(self, standard_jd):
        """Test that ,alg finds Algol (pyswisseph behavior)."""
        pos, name, retflag = swe_fixstar_ut(",alg", standard_jd, 0)
        assert name == "Algol"
        # Algol is around 26 Taurus = 56 degrees
        assert 55 < pos[0] < 57

    def test_alpha_leo_finds_regulus(self, standard_jd):
        """Test that Alpha Leo finds Regulus (pyswisseph behavior)."""
        pos, name, retflag = swe_fixstar_ut("Alpha Leo", standard_jd, 0)
        assert name == "Regulus"

    def test_sirius_case_variations_identical(self, standard_jd):
        """Test SIRIUS and sirius return identical results."""
        pos_upper, name_upper, _ = swe_fixstar_ut("SIRIUS", standard_jd, 0)
        pos_lower, name_lower, _ = swe_fixstar_ut("sirius", standard_jd, 0)
        pos_mixed, name_mixed, _ = swe_fixstar_ut("Sirius", standard_jd, 0)

        # All should return Sirius as canonical name
        assert name_upper == name_lower == name_mixed == "Sirius"

        # All should return identical positions
        assert abs(pos_upper[0] - pos_lower[0]) < 0.0001
        assert abs(pos_upper[0] - pos_mixed[0]) < 0.0001

    def test_alpha_cma_finds_sirius(self, standard_jd):
        """Test that Alpha CMa finds Sirius (pyswisseph behavior)."""
        pos, name, retflag = swe_fixstar_ut("Alpha CMa", standard_jd, 0)
        assert name == "Sirius"


@pytest.mark.unit
class TestReturnTypeStructure:
    """Tests for correct return type structure."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch as Julian Day."""
        return 2451545.0

    def test_swe_fixstar_ut_return_structure(self, standard_jd):
        """Test swe_fixstar_ut returns correct tuple structure."""
        result = swe_fixstar_ut("Regulus", standard_jd, 0)

        assert isinstance(result, tuple)
        assert len(result) == 3

        pos, name, retflag = result

        # Position tuple should have 6 elements
        assert isinstance(pos, tuple)
        assert len(pos) == 6

        # name should be a string (canonical star name)
        assert isinstance(name, str)
        assert name == "Regulus"

        # retflag should be an int
        assert isinstance(retflag, int)

    def test_swe_fixstar_return_structure(self, standard_jd):
        """Test swe_fixstar returns correct tuple structure."""
        result = swe_fixstar("Spica", standard_jd, 0)

        assert isinstance(result, tuple)
        assert len(result) == 3

        pos, name, retflag = result

        assert isinstance(pos, tuple)
        assert len(pos) == 6
        assert isinstance(name, str)
        assert name == "Spica"
        assert isinstance(retflag, int)
