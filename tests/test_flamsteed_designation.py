"""
Unit tests for Flamsteed designation star search.

Tests the ability to search for stars by their Flamsteed designations using
number + constellation genitive form, e.g., "32 Leonis", "87 Virginis".

This feature is required for swe_fixstar2 compatibility with designation search.
"""

import pytest

from libephemeris.fixed_stars import (
    _parse_flamsteed_designation,
    resolve_star_name,
    _resolve_star2,
    swe_fixstar2_ut,
    swe_fixstar2,
    CONSTELLATION_ABBREV,
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
    SE_ARCTURUS,
    SE_DENEB,
    SE_CASTOR,
    SE_POLLUX,
    SE_ALCYONE,
    SE_ASTEROPE,
    SE_CELAENO,
    SE_ELECTRA,
    SE_MAIA,
    SE_MEROPE,
    SE_TAYGETA,
    SE_ATLAS,
    SE_PLEIONE,
)


@pytest.mark.unit
class TestFlamsteedDesignationParsing:
    """Tests for _parse_flamsteed_designation function."""

    def test_parse_32_leonis(self):
        """32 Leonis should parse to '32 LEO'."""
        result = _parse_flamsteed_designation("32 Leonis")
        assert result == "32 LEO"

    def test_parse_87_leonis(self):
        """87 Leonis (Regulus) should parse to '87 LEO'."""
        result = _parse_flamsteed_designation("87 Leonis")
        assert result == "87 LEO"

    def test_parse_67_virginis(self):
        """67 Virginis (Spica) should parse to '67 VIR'."""
        result = _parse_flamsteed_designation("67 Virginis")
        assert result == "67 VIR"

    def test_parse_26_persei(self):
        """26 Persei (Algol) should parse to '26 PER'."""
        result = _parse_flamsteed_designation("26 Persei")
        assert result == "26 PER"

    def test_parse_9_canis_majoris(self):
        """9 Canis Majoris (Sirius) should parse to '9 CMA'."""
        result = _parse_flamsteed_designation("9 Canis Majoris")
        assert result == "9 CMA"

    def test_parse_87_tauri(self):
        """87 Tauri (Aldebaran) should parse to '87 TAU'."""
        result = _parse_flamsteed_designation("87 Tauri")
        assert result == "87 TAU"

    def test_parse_21_scorpii(self):
        """21 Scorpii (Antares) should parse to '21 SCO'."""
        result = _parse_flamsteed_designation("21 Scorpii")
        assert result == "21 SCO"

    def test_parse_3_lyrae(self):
        """3 Lyrae (Vega) should parse to '3 LYR'."""
        result = _parse_flamsteed_designation("3 Lyrae")
        assert result == "3 LYR"

    def test_parse_case_insensitive(self):
        """Parsing should be case-insensitive."""
        result1 = _parse_flamsteed_designation("87 LEONIS")
        result2 = _parse_flamsteed_designation("87 leonis")
        result3 = _parse_flamsteed_designation("87 Leonis")
        assert result1 == result2 == result3 == "87 LEO"

    def test_parse_nominative_form(self):
        """Should support nominative constellation names like '87 Leo'."""
        result = _parse_flamsteed_designation("87 Leo")
        assert result == "87 LEO"

    def test_parse_invalid_no_number(self):
        """Non-numeric first part should return None."""
        result = _parse_flamsteed_designation("Alpha Leonis")
        assert result is None

    def test_parse_invalid_constellation(self):
        """Invalid constellation should return None."""
        result = _parse_flamsteed_designation("87 Xyzzy")
        assert result is None

    def test_parse_empty_string(self):
        """Empty string should return None."""
        result = _parse_flamsteed_designation("")
        assert result is None

    def test_parse_number_only(self):
        """Number only (missing constellation) should return None."""
        result = _parse_flamsteed_designation("87")
        assert result is None

    def test_parse_multi_word_constellation(self):
        """Should support multi-word constellations like 'Ursa Major'."""
        result = _parse_flamsteed_designation("50 Ursae Majoris")
        assert result == "50 UMA"

    def test_parse_with_whitespace(self):
        """Should handle extra whitespace."""
        result = _parse_flamsteed_designation("  87   Leonis  ")
        assert result == "87 LEO"

    # Pleiades cluster stars (Flamsteed numbers in Taurus)
    def test_parse_21_tauri(self):
        """21 Tauri (Asterope) should parse to '21 TAU'."""
        result = _parse_flamsteed_designation("21 Tauri")
        assert result == "21 TAU"

    def test_parse_16_tauri(self):
        """16 Tauri (Celaeno) should parse to '16 TAU'."""
        result = _parse_flamsteed_designation("16 Tauri")
        assert result == "16 TAU"

    def test_parse_17_tauri(self):
        """17 Tauri (Electra) should parse to '17 TAU'."""
        result = _parse_flamsteed_designation("17 Tauri")
        assert result == "17 TAU"

    def test_parse_20_tauri(self):
        """20 Tauri (Maia) should parse to '20 TAU'."""
        result = _parse_flamsteed_designation("20 Tauri")
        assert result == "20 TAU"

    def test_parse_23_tauri(self):
        """23 Tauri (Merope) should parse to '23 TAU'."""
        result = _parse_flamsteed_designation("23 Tauri")
        assert result == "23 TAU"


@pytest.mark.unit
class TestConstellationCoverage:
    """Test that major constellations are covered for Flamsteed parsing."""

    def test_zodiacal_constellations(self):
        """All zodiacal constellations should be parseable."""
        zodiacal = [
            ("1 Arietis", "ARI"),
            ("1 Tauri", "TAU"),
            ("1 Geminorum", "GEM"),
            ("1 Cancri", "CNC"),
            ("1 Leonis", "LEO"),
            ("1 Virginis", "VIR"),
            ("1 Librae", "LIB"),
            ("1 Scorpii", "SCO"),
            ("1 Sagittarii", "SGR"),
            ("1 Capricorni", "CAP"),
            ("1 Aquarii", "AQR"),
            ("1 Piscium", "PSC"),
        ]
        for designation, expected_abbrev in zodiacal:
            result = _parse_flamsteed_designation(designation)
            assert result == f"1 {expected_abbrev}", f"Failed for {designation}"

    def test_prominent_constellations(self):
        """Prominent constellations should be parseable."""
        constellations = [
            ("1 Orionis", "ORI"),
            ("1 Cygni", "CYG"),
            ("1 Lyrae", "LYR"),
            ("1 Aquilae", "AQL"),
            ("1 Persei", "PER"),
            ("1 Ursae Majoris", "UMA"),
            ("1 Ursae Minoris", "UMI"),
            ("1 Draconis", "DRA"),
            ("1 Centauri", "CEN"),
            ("1 Crucis", "CRU"),
        ]
        for designation, expected_abbrev in constellations:
            result = _parse_flamsteed_designation(designation)
            assert result == f"1 {expected_abbrev}", f"Failed for {designation}"


@pytest.mark.unit
class TestResolveStarNameWithFlamsteed:
    """Tests for resolve_star_name with Flamsteed designations."""

    def test_resolve_87_leonis(self):
        """87 Leonis should resolve to Regulus (SE_REGULUS)."""
        result = resolve_star_name("87 Leonis")
        assert result == SE_REGULUS

    def test_resolve_87_leo(self):
        """87 Leo should also resolve to Regulus."""
        result = resolve_star_name("87 Leo")
        assert result == SE_REGULUS

    def test_resolve_67_virginis(self):
        """67 Virginis should resolve to Spica (SE_SPICA_STAR)."""
        result = resolve_star_name("67 Virginis")
        assert result == SE_SPICA_STAR

    def test_resolve_26_persei(self):
        """26 Persei should resolve to Algol (SE_ALGOL)."""
        result = resolve_star_name("26 Persei")
        assert result == SE_ALGOL

    def test_resolve_9_canis_majoris(self):
        """9 Canis Majoris should resolve to Sirius (SE_SIRIUS)."""
        result = resolve_star_name("9 Canis Majoris")
        assert result == SE_SIRIUS

    def test_resolve_87_tauri(self):
        """87 Tauri should resolve to Aldebaran (SE_ALDEBARAN)."""
        result = resolve_star_name("87 Tauri")
        assert result == SE_ALDEBARAN

    def test_resolve_21_scorpii(self):
        """21 Scorpii should resolve to Antares (SE_ANTARES)."""
        result = resolve_star_name("21 Scorpii")
        assert result == SE_ANTARES

    def test_resolve_3_lyrae(self):
        """3 Lyrae should resolve to Vega (SE_VEGA)."""
        result = resolve_star_name("3 Lyrae")
        assert result == SE_VEGA

    def test_resolve_1_ursae_minoris(self):
        """1 Ursae Minoris should resolve to Polaris (SE_POLARIS)."""
        result = resolve_star_name("1 Ursae Minoris")
        assert result == SE_POLARIS

    def test_resolve_16_bootis(self):
        """16 Bootis should resolve to Arcturus (SE_ARCTURUS)."""
        result = resolve_star_name("16 Bootis")
        assert result == SE_ARCTURUS

    def test_resolve_50_cygni(self):
        """50 Cygni should resolve to Deneb (SE_DENEB)."""
        result = resolve_star_name("50 Cygni")
        assert result == SE_DENEB

    def test_resolve_66_geminorum(self):
        """66 Geminorum should resolve to Castor (SE_CASTOR)."""
        result = resolve_star_name("66 Geminorum")
        assert result == SE_CASTOR

    def test_resolve_78_geminorum(self):
        """78 Geminorum should resolve to Pollux (SE_POLLUX)."""
        result = resolve_star_name("78 Geminorum")
        assert result == SE_POLLUX

    # Pleiades cluster tests
    def test_resolve_25_tauri(self):
        """25 Tauri should resolve to Alcyone (SE_ALCYONE)."""
        result = resolve_star_name("25 Tauri")
        assert result == SE_ALCYONE

    def test_resolve_21_tauri(self):
        """21 Tauri should resolve to Asterope (SE_ASTEROPE)."""
        result = resolve_star_name("21 Tauri")
        assert result == SE_ASTEROPE

    def test_resolve_case_insensitive(self):
        """Resolution should be case-insensitive."""
        result1 = resolve_star_name("87 LEONIS")
        result2 = resolve_star_name("87 leonis")
        result3 = resolve_star_name("87 Leonis")
        assert result1 == result2 == result3 == SE_REGULUS


@pytest.mark.unit
class TestResolve2StarWithFlamsteed:
    """Tests for _resolve_star2 with Flamsteed designations."""

    def test_resolve2_87_leonis(self):
        """87 Leonis should resolve to Regulus catalog entry."""
        entry, err = _resolve_star2("87 Leonis")
        assert err is None
        assert entry is not None
        assert entry.id == SE_REGULUS
        assert entry.name == "Regulus"

    def test_resolve2_67_virginis(self):
        """67 Virginis should resolve to Spica catalog entry."""
        entry, err = _resolve_star2("67 Virginis")
        assert err is None
        assert entry is not None
        assert entry.id == SE_SPICA_STAR
        assert entry.name == "Spica"

    def test_resolve2_26_persei(self):
        """26 Persei should resolve to Algol catalog entry."""
        entry, err = _resolve_star2("26 Persei")
        assert err is None
        assert entry is not None
        assert entry.id == SE_ALGOL
        assert entry.name == "Algol"

    def test_resolve2_star_not_in_catalog(self):
        """Flamsteed designation for star not in catalog should return error."""
        entry, err = _resolve_star2("99 Leonis")  # Not in catalog
        assert entry is None
        assert err is not None
        assert "could not find" in err.lower()


@pytest.mark.unit
class TestSweFixstar2WithFlamsteed:
    """Tests for swe_fixstar2 and swe_fixstar2_ut with Flamsteed designations."""

    def test_swe_fixstar2_ut_87_leonis(self):
        """swe_fixstar2_ut should find Regulus via '87 Leonis'."""
        jd = 2451545.0  # J2000.0
        name, pos, flag, err = swe_fixstar2_ut("87 Leonis", jd, 0)
        assert err == "", f"Unexpected error: {err}"
        assert "Regulus" in name
        assert 110 < pos[0] < 160  # Rough longitude check

    def test_swe_fixstar2_ut_67_virginis(self):
        """swe_fixstar2_ut should find Spica via '67 Virginis'."""
        jd = 2451545.0
        name, pos, flag, err = swe_fixstar2_ut("67 Virginis", jd, 0)
        assert err == "", f"Unexpected error: {err}"
        assert "Spica" in name
        assert 200 < pos[0] < 210  # Rough longitude check

    def test_swe_fixstar2_26_persei(self):
        """swe_fixstar2 (TT) should find Algol via '26 Persei'."""
        jd = 2451545.0
        name, pos, flag, err = swe_fixstar2("26 Persei", jd, 0)
        assert err == "", f"Unexpected error: {err}"
        assert "Algol" in name

    def test_swe_fixstar2_21_tauri(self):
        """swe_fixstar2 should find Asterope via '21 Tauri'."""
        jd = 2451545.0
        name, pos, flag, err = swe_fixstar2("21 Tauri", jd, 0)
        assert err == "", f"Unexpected error: {err}"
        assert "Asterope" in name

    def test_swe_fixstar2_ut_star_not_found(self):
        """swe_fixstar2_ut should return error for star not in catalog."""
        jd = 2451545.0
        name, pos, flag, err = swe_fixstar2_ut("99 Leonis", jd, 0)
        assert err != ""
        assert "could not find" in err.lower()


@pytest.mark.unit
class TestFlamsteedEdgeCases:
    """Edge case tests for Flamsteed designation handling."""

    def test_leading_zeros_not_supported(self):
        """Numbers with leading zeros are still valid Flamsteed numbers."""
        # "087 Leonis" should work the same as "87 Leonis"
        result = _parse_flamsteed_designation("087 Leonis")
        # Python's isdigit() accepts "087" as digits
        assert result == "087 LEO"

    def test_very_large_number(self):
        """Large Flamsteed numbers should parse correctly."""
        result = _parse_flamsteed_designation("999 Leonis")
        assert result == "999 LEO"

    def test_single_digit_number(self):
        """Single digit Flamsteed numbers should parse correctly."""
        result = _parse_flamsteed_designation("1 Leonis")
        assert result == "1 LEO"

    def test_mixed_case_constellation(self):
        """Mixed case constellation names should work."""
        result = _parse_flamsteed_designation("87 lEoNiS")
        assert result == "87 LEO"

    def test_extra_spaces(self):
        """Extra spaces should be handled correctly."""
        result = _parse_flamsteed_designation("  87   Leonis  ")
        assert result == "87 LEO"
