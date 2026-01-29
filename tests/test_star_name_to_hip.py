"""
Unit tests for the STAR_NAME_TO_HIP mapping and get_hip_from_star_name function.

Tests the mapping from common star names, Bayer designations, and Flamsteed
numbers to Hipparcos (HIP) catalog numbers.

Data sourced from IAU Working Group on Star Names (WGSN) catalog.
"""

import pytest
from libephemeris.fixed_stars import (
    STAR_NAME_TO_HIP,
    get_hip_from_star_name,
    STAR_CATALOG,
)


@pytest.mark.unit
class TestStarNameToHipMapping:
    """Tests for the STAR_NAME_TO_HIP dictionary structure."""

    def test_mapping_has_minimum_entries(self):
        """Verify STAR_NAME_TO_HIP has a substantial number of entries."""
        # Should have at least 400 entries (common names + Bayer + Flamsteed)
        assert len(STAR_NAME_TO_HIP) >= 400, (
            f"STAR_NAME_TO_HIP should have at least 400 entries, "
            f"but has {len(STAR_NAME_TO_HIP)}"
        )

    def test_mapping_keys_are_uppercase(self):
        """Verify all mapping keys are uppercase for consistent lookup."""
        for key in STAR_NAME_TO_HIP:
            assert key == key.upper(), f"Key '{key}' should be uppercase"

    def test_mapping_values_are_integers(self):
        """Verify all mapping values are integers (HIP numbers)."""
        for name, hip in STAR_NAME_TO_HIP.items():
            assert isinstance(hip, int), (
                f"HIP number for '{name}' should be int, got {type(hip)}"
            )

    def test_no_zero_hip_numbers(self):
        """Verify there are no zero HIP numbers (invalid)."""
        for name, hip in STAR_NAME_TO_HIP.items():
            # HIP numbers should be positive (or -1 for no HIP)
            assert hip != 0, f"HIP number for '{name}' should not be 0"


@pytest.mark.unit
class TestGetHipFromStarName:
    """Tests for the get_hip_from_star_name function."""

    def test_empty_name_returns_none(self):
        """Empty string should return None."""
        assert get_hip_from_star_name("") is None
        assert get_hip_from_star_name("   ") is None

    def test_none_input_returns_none(self):
        """None input should return None."""
        assert get_hip_from_star_name(None) is None

    def test_unknown_name_returns_none(self):
        """Unknown star name should return None."""
        assert get_hip_from_star_name("UnknownStar") is None
        assert get_hip_from_star_name("NotARealStar") is None
        assert get_hip_from_star_name("xyz123") is None

    def test_case_insensitive_lookup(self):
        """Lookup should be case-insensitive."""
        # Test various case combinations
        assert get_hip_from_star_name("Regulus") == 49669
        assert get_hip_from_star_name("REGULUS") == 49669
        assert get_hip_from_star_name("regulus") == 49669
        assert get_hip_from_star_name("ReGuLuS") == 49669

        assert get_hip_from_star_name("Sirius") == 32349
        assert get_hip_from_star_name("SIRIUS") == 32349
        assert get_hip_from_star_name("sirius") == 32349

    def test_whitespace_stripping(self):
        """Whitespace should be stripped from input."""
        assert get_hip_from_star_name("  Regulus  ") == 49669
        assert get_hip_from_star_name("\tSirius\n") == 32349
        assert get_hip_from_star_name("  Alpha Leo  ") == 49669


@pytest.mark.unit
class TestCommonStarNames:
    """Tests for common/proper star names (IAU-approved)."""

    def test_royal_stars_of_persia(self):
        """Test the four Royal Stars of Persia."""
        # Aldebaran - Watcher of the East
        assert get_hip_from_star_name("Aldebaran") == 21421
        # Regulus - Watcher of the North
        assert get_hip_from_star_name("Regulus") == 49669
        # Antares - Watcher of the West
        assert get_hip_from_star_name("Antares") == 80763
        # Fomalhaut - Watcher of the South
        assert get_hip_from_star_name("Fomalhaut") == 113368

    def test_brightest_stars(self):
        """Test the brightest stars in the night sky."""
        assert get_hip_from_star_name("Sirius") == 32349  # Brightest
        assert get_hip_from_star_name("Canopus") == 30438
        assert get_hip_from_star_name("Arcturus") == 69673
        assert get_hip_from_star_name("Vega") == 91262
        assert get_hip_from_star_name("Capella") == 24608
        assert get_hip_from_star_name("Rigel") == 24436
        assert get_hip_from_star_name("Procyon") == 37279
        assert get_hip_from_star_name("Betelgeuse") == 27989
        assert get_hip_from_star_name("Achernar") == 7588
        assert get_hip_from_star_name("Altair") == 97649
        assert get_hip_from_star_name("Deneb") == 102098

    def test_navigation_stars(self):
        """Test important navigation stars."""
        assert get_hip_from_star_name("Polaris") == 11767  # North Star
        assert get_hip_from_star_name("Acrux") == 60718
        assert get_hip_from_star_name("Canopus") == 30438
        assert get_hip_from_star_name("Alpheratz") == 677

    def test_pleiades_stars(self):
        """Test the named Pleiades stars."""
        assert get_hip_from_star_name("Alcyone") == 17702
        assert get_hip_from_star_name("Atlas") == 17847
        assert get_hip_from_star_name("Electra") == 17499
        assert get_hip_from_star_name("Maia") == 17573
        assert get_hip_from_star_name("Merope") == 17608
        assert get_hip_from_star_name("Taygeta") == 17531
        assert get_hip_from_star_name("Pleione") == 17851
        assert get_hip_from_star_name("Celaeno") == 17489
        assert get_hip_from_star_name("Asterope") == 17579

    def test_big_dipper_stars(self):
        """Test the Big Dipper (Ursa Major) stars."""
        assert get_hip_from_star_name("Dubhe") == 54061
        assert get_hip_from_star_name("Merak") == 53910
        assert get_hip_from_star_name("Phecda") == 58001
        assert get_hip_from_star_name("Megrez") == 59774
        assert get_hip_from_star_name("Alioth") == 62956
        assert get_hip_from_star_name("Mizar") == 65378
        assert get_hip_from_star_name("Alkaid") == 67301
        assert get_hip_from_star_name("Alcor") == 65477

    def test_orion_belt_stars(self):
        """Test Orion's belt stars."""
        assert get_hip_from_star_name("Alnitak") == 26727
        assert get_hip_from_star_name("Alnilam") == 26311
        assert get_hip_from_star_name("Mintaka") == 25930

    def test_astrologically_significant_stars(self):
        """Test astrologically significant fixed stars."""
        assert get_hip_from_star_name("Spica") == 65474
        assert get_hip_from_star_name("Algol") == 14576
        assert get_hip_from_star_name("Alphecca") == 76267
        assert get_hip_from_star_name("Zubenelgenubi") == 72622
        assert get_hip_from_star_name("Zubeneschamali") == 74785


@pytest.mark.unit
class TestBayerDesignations:
    """Tests for Bayer designation lookups."""

    def test_alpha_stars_full_name(self):
        """Test Alpha stars with full constellation names."""
        assert get_hip_from_star_name("Alpha Leonis") == 49669  # Regulus
        assert get_hip_from_star_name("Alpha Virginis") == 65474  # Spica
        assert get_hip_from_star_name("Alpha Scorpii") == 80763  # Antares
        assert get_hip_from_star_name("Alpha Tauri") == 21421  # Aldebaran
        assert get_hip_from_star_name("Alpha Lyrae") == 91262  # Vega
        assert get_hip_from_star_name("Alpha Orionis") == 27989  # Betelgeuse
        assert get_hip_from_star_name("Alpha Canis Majoris") == 32349  # Sirius

    def test_alpha_stars_abbreviated(self):
        """Test Alpha stars with abbreviated constellation names."""
        assert get_hip_from_star_name("Alpha Leo") == 49669  # Regulus
        assert get_hip_from_star_name("Alpha Vir") == 65474  # Spica
        assert get_hip_from_star_name("Alpha Sco") == 80763  # Antares
        assert get_hip_from_star_name("Alpha Tau") == 21421  # Aldebaran
        assert get_hip_from_star_name("Alpha Lyr") == 91262  # Vega
        assert get_hip_from_star_name("Alpha Ori") == 27989  # Betelgeuse
        assert get_hip_from_star_name("Alpha CMA") == 32349  # Sirius

    def test_beta_stars(self):
        """Test Beta star designations."""
        assert get_hip_from_star_name("Beta Orionis") == 24436  # Rigel
        assert get_hip_from_star_name("Beta Ori") == 24436
        assert get_hip_from_star_name("Beta Persei") == 14576  # Algol
        assert get_hip_from_star_name("Beta Per") == 14576
        assert get_hip_from_star_name("Beta Geminorum") == 37826  # Pollux
        assert get_hip_from_star_name("Beta Gem") == 37826
        assert get_hip_from_star_name("Beta Leonis") == 57632  # Denebola
        assert get_hip_from_star_name("Beta Leo") == 57632

    def test_gamma_stars(self):
        """Test Gamma star designations."""
        assert get_hip_from_star_name("Gamma Orionis") == 25336  # Bellatrix
        assert get_hip_from_star_name("Gamma Ori") == 25336
        assert get_hip_from_star_name("Gamma Leonis") == 50583  # Algieba
        assert get_hip_from_star_name("Gamma Leo") == 50583
        assert get_hip_from_star_name("Gamma Crucis") == 61084  # Gacrux
        assert get_hip_from_star_name("Gamma Cru") == 61084

    def test_other_greek_letters(self):
        """Test other Greek letter designations."""
        # Delta stars
        assert get_hip_from_star_name("Delta Orionis") == 25930  # Mintaka
        assert get_hip_from_star_name("Delta Leonis") == 54872  # Zosma
        # Epsilon stars
        assert get_hip_from_star_name("Epsilon Orionis") == 26311  # Alnilam
        assert get_hip_from_star_name("Epsilon Sagittarii") == 90185  # Kaus Australis
        # Zeta stars
        assert get_hip_from_star_name("Zeta Orionis") == 26727  # Alnitak
        # Lambda stars
        assert get_hip_from_star_name("Lambda Scorpii") == 85927  # Shaula


@pytest.mark.unit
class TestFlamsteedDesignations:
    """Tests for Flamsteed number lookups."""

    def test_pleiades_flamsteed_numbers(self):
        """Test Pleiades stars by Flamsteed number."""
        assert get_hip_from_star_name("16 Tauri") == 17489  # Celaeno
        assert get_hip_from_star_name("17 Tauri") == 17499  # Electra
        assert get_hip_from_star_name("19 Tauri") == 17531  # Taygeta
        assert get_hip_from_star_name("20 Tauri") == 17573  # Maia
        assert get_hip_from_star_name("21 Tauri") == 17579  # Asterope
        assert get_hip_from_star_name("23 Tauri") == 17608  # Merope
        assert get_hip_from_star_name("27 Tauri") == 17847  # Atlas
        assert get_hip_from_star_name("28 Tauri") == 17851  # Pleione

    def test_notable_flamsteed_stars(self):
        """Test other notable stars with Flamsteed numbers."""
        assert get_hip_from_star_name("51 Pegasi") == 113357  # Helvetios
        assert get_hip_from_star_name("55 Cancri") == 43587  # Copernicus
        assert get_hip_from_star_name("47 Ursae Majoris") == 53721  # Chalawan
        assert get_hip_from_star_name("80 Ursae Majoris") == 65477  # Alcor


@pytest.mark.unit
class TestAlternativeNames:
    """Tests for alternative names and historical spellings."""

    def test_common_alternative_names(self):
        """Test commonly used alternative names."""
        assert get_hip_from_star_name("Dog Star") == 32349  # Sirius
        assert get_hip_from_star_name("North Star") == 11767  # Polaris
        assert get_hip_from_star_name("Pole Star") == 11767  # Polaris
        assert get_hip_from_star_name("Cor Leonis") == 49669  # Regulus

    def test_historical_spellings(self):
        """Test historical spelling variants."""
        assert get_hip_from_star_name("Wega") == 91262  # Vega
        assert get_hip_from_star_name("Benetnash") == 67301  # Alkaid
        assert get_hip_from_star_name("Agena") == 68702  # Hadar

    def test_latin_names(self):
        """Test Latin constellation component names."""
        assert get_hip_from_star_name("Becrux") == 62434  # Beta Crucis = Mimosa


@pytest.mark.unit
class TestStarsWithoutHipNumbers:
    """Tests for stars without HIP numbers (returns -1)."""

    def test_exoplanet_host_stars(self):
        """Test exoplanet host stars discovered by transit (no HIP)."""
        # These stars don't have HIP numbers, return -1
        assert get_hip_from_star_name("Absolutno") == -1  # XO-5
        assert get_hip_from_star_name("Amansinaya") == -1  # WASP-34
        assert get_hip_from_star_name("Anadolu") == -1  # WASP-52


@pytest.mark.unit
class TestConsistencyWithStarCatalog:
    """Tests to verify consistency with existing STAR_CATALOG."""

    def test_catalog_stars_have_matching_hip(self):
        """Verify stars in STAR_CATALOG have matching HIP in STAR_NAME_TO_HIP."""
        for entry in STAR_CATALOG:
            # Get HIP from our new mapping
            hip = get_hip_from_star_name(entry.name)
            if hip is not None and hip > 0:
                assert hip == entry.hip_number, (
                    f"HIP mismatch for {entry.name}: "
                    f"STAR_NAME_TO_HIP has {hip}, STAR_CATALOG has {entry.hip_number}"
                )

    def test_selected_catalog_stars(self):
        """Test selected stars are consistent between mappings."""
        test_stars = [
            ("Regulus", 49669),
            ("Spica", 65474),
            ("Aldebaran", 21421),
            ("Antares", 80763),
            ("Sirius", 32349),
            ("Vega", 91262),
            ("Arcturus", 69673),
            ("Capella", 24608),
            ("Procyon", 37279),
            ("Pollux", 37826),
        ]
        for name, expected_hip in test_stars:
            assert get_hip_from_star_name(name) == expected_hip


@pytest.mark.unit
class TestSpecificIauStars:
    """Tests for specific IAU-approved star names."""

    def test_recently_approved_names(self):
        """Test stars with recently IAU-approved names."""
        # These names were approved by IAU WGSN
        assert get_hip_from_star_name("Cervantes") == 86796  # Mu Arae
        assert get_hip_from_star_name("Copernicus") == 43587  # 55 Cancri
        assert get_hip_from_star_name("Helvetios") == 113357  # 51 Pegasi
        assert get_hip_from_star_name("Ran") == 16537  # Epsilon Eridani
        assert get_hip_from_star_name("Fomalhaut") == 113368

    def test_southern_hemisphere_stars(self):
        """Test Southern Cross and southern navigation stars."""
        assert get_hip_from_star_name("Acrux") == 60718  # Alpha Crucis
        assert get_hip_from_star_name("Mimosa") == 62434  # Beta Crucis
        assert get_hip_from_star_name("Gacrux") == 61084  # Gamma Crucis
        assert get_hip_from_star_name("Imai") == 59747  # Delta Crucis
        assert get_hip_from_star_name("Hadar") == 68702  # Beta Centauri
        assert get_hip_from_star_name("Rigil Kentaurus") == 71683  # Alpha Centauri
