"""
Unit tests for zodiacal constellation bright stars.

This module tests the bright stars from each zodiacal constellation used in
astrological interpretation. The 12 zodiacal constellations lie along the
ecliptic and are associated with the astrological signs.

Zodiacal constellations and their key stars:
- Aries: Hamal, Sheratan, Mesarthim
- Taurus: Aldebaran, Elnath (+ Pleiades, Hyades clusters)
- Gemini: Pollux, Castor
- Cancer: Acubens, Tarf, Asellus Borealis, Asellus Australis
- Leo: Regulus, Denebola, Algieba, Zosma
- Virgo: Spica, Vindemiatrix
- Libra: Zubenelgenubi, Zubeneschamali
- Scorpius: Antares, Shaula, Sargas, Dschubba, Graffias, Lesath
- Sagittarius: Kaus Australis, Nunki, Kaus Media, Kaus Borealis, Ascella
- Capricornus: Deneb Algedi, Algedi, Dabih, Nashira
- Aquarius: Sadalsuud, Sadalmelik, Skat
- Pisces: Eta Piscium, Alrescha
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    # Aries
    SE_HAMAL,
    SE_SHERATAN,
    SE_MESARTHIM,
    # Cancer
    SE_ACUBENS,
    SE_TARF,
    SE_ASELLUS_BOREALIS,
    SE_ASELLUS_AUSTRALIS,
    # Sagittarius
    SE_KAUS_AUSTRALIS,
    SE_NUNKI,
    SE_KAUS_MEDIA,
    SE_KAUS_BOREALIS,
    SE_ASCELLA,
    # Capricornus
    SE_DENEB_ALGEDI,
    SE_ALGEDI,
    SE_DABIH,
    SE_NASHIRA,
    # Aquarius
    SE_SADALSUUD,
    SE_SADALMELIK,
    SE_SKAT,
    # Pisces
    SE_ETA_PISCIUM,
    SE_ALRESCHA,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# ======== ARIES CONSTELLATION STARS ========
ARIES_STARS = [
    (SE_HAMAL, "Hamal", 9884, 2.00),  # Alpha Ari - brightest
    (SE_SHERATAN, "Sheratan", 8903, 2.64),  # Beta Ari
    (SE_MESARTHIM, "Mesarthim", 8832, 3.88),  # Gamma Ari
]


# ======== CANCER CONSTELLATION STARS ========
CANCER_STARS = [
    (SE_TARF, "Tarf", 42911, 3.52),  # Beta Cnc - brightest
    (SE_ASELLUS_AUSTRALIS, "Asellus Australis", 42911, 3.94),  # Delta Cnc
    (SE_ACUBENS, "Acubens", 44066, 4.25),  # Alpha Cnc
    (SE_ASELLUS_BOREALIS, "Asellus Borealis", 42806, 4.66),  # Gamma Cnc
]


# ======== SAGITTARIUS CONSTELLATION STARS ========
SAGITTARIUS_STARS = [
    (SE_KAUS_AUSTRALIS, "Kaus Australis", 90185, 1.85),  # Epsilon Sgr - brightest
    (SE_NUNKI, "Nunki", 92855, 2.02),  # Sigma Sgr
    (SE_ASCELLA, "Ascella", 93506, 2.59),  # Zeta Sgr
    (SE_KAUS_MEDIA, "Kaus Media", 89931, 2.70),  # Delta Sgr
    (SE_KAUS_BOREALIS, "Kaus Borealis", 90496, 2.81),  # Lambda Sgr
]


# ======== CAPRICORNUS CONSTELLATION STARS ========
CAPRICORNUS_STARS = [
    (SE_DENEB_ALGEDI, "Deneb Algedi", 107556, 2.81),  # Delta Cap - brightest
    (SE_DABIH, "Dabih", 100345, 3.08),  # Beta Cap
    (SE_ALGEDI, "Algedi", 100064, 3.57),  # Alpha Cap
    (SE_NASHIRA, "Nashira", 106985, 3.68),  # Gamma Cap
]


# ======== AQUARIUS CONSTELLATION STARS ========
AQUARIUS_STARS = [
    (SE_SADALSUUD, "Sadalsuud", 106278, 2.87),  # Beta Aqr - brightest
    (SE_SADALMELIK, "Sadalmelik", 109074, 2.96),  # Alpha Aqr
    (SE_SKAT, "Skat", 113136, 3.27),  # Delta Aqr
]


# ======== PISCES CONSTELLATION STARS ========
PISCES_STARS = [
    (SE_ETA_PISCIUM, "Eta Piscium", 5742, 3.62),  # Eta Psc - brightest
    (SE_ALRESCHA, "Alrescha", 7097, 3.82),  # Alpha Psc
]


# Combine all new zodiacal stars for parametrized tests
ALL_NEW_ZODIACAL_STARS = (
    ARIES_STARS
    + CANCER_STARS
    + SAGITTARIUS_STARS
    + CAPRICORNUS_STARS
    + AQUARIUS_STARS
    + PISCES_STARS
)


@pytest.mark.unit
class TestZodiacalStarsCatalog:
    """Test that all zodiacal constellation stars are in the catalog."""

    def test_all_new_zodiacal_stars_present(self):
        """Verify all new zodiacal stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in ALL_NEW_ZODIACAL_STARS:
            assert star_id in catalog_ids, (
                f"Zodiacal star {name} (ID={star_id}) not in catalog"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", ALL_NEW_ZODIACAL_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each zodiacal star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Stars should have proper motion data
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", ALL_NEW_ZODIACAL_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each zodiacal star has correct Hipparcos number."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.hip_number == hip, (
            f"Star {name} should have HIP {hip}, got {entry.hip_number}"
        )


@pytest.mark.unit
class TestZodiacalStarsCalculation:
    """Test position calculations for zodiacal constellation stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", ALL_NEW_ZODIACAL_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each zodiacal star returns a reasonable position."""
        pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_aries_stars_in_aries_region(self, standard_jd):
        """Test that Aries stars are in the Aries/Taurus ecliptic region."""
        for star_id, name, _, _ in ARIES_STARS:
            pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)
            # Aries stars should be between ~0 and ~45 degrees ecliptic
            assert 0 < pos[0] < 50, (
                f"{name} should be in Aries region, got {pos[0]:.1f} degrees"
            )

    def test_sagittarius_stars_in_sagittarius_region(self, standard_jd):
        """Test that Sagittarius stars are in the Sagittarius ecliptic region."""
        for star_id, name, _, _ in SAGITTARIUS_STARS:
            pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)
            # Sagittarius stars should be between ~260 and ~300 degrees ecliptic
            assert 255 < pos[0] < 310, (
                f"{name} should be in Sagittarius region, got {pos[0]:.1f} degrees"
            )

    def test_kaus_australis_is_brightest_sagittarius(self, standard_jd):
        """Test that Kaus Australis is the brightest Sagittarius star."""
        kaus_entry = None
        for e in STAR_CATALOG:
            if e.id == SE_KAUS_AUSTRALIS:
                kaus_entry = e
                break

        assert kaus_entry is not None, "Kaus Australis not found"

        for star_id, name, _, mag in SAGITTARIUS_STARS:
            if star_id != SE_KAUS_AUSTRALIS:
                assert kaus_entry.magnitude < mag, (
                    f"Kaus Australis ({kaus_entry.magnitude}) "
                    f"should be brighter than {name} ({mag})"
                )

    def test_hamal_is_brightest_aries(self, standard_jd):
        """Test that Hamal is the brightest Aries star."""
        hamal_entry = None
        for e in STAR_CATALOG:
            if e.id == SE_HAMAL:
                hamal_entry = e
                break

        assert hamal_entry is not None, "Hamal not found"

        for star_id, name, _, mag in ARIES_STARS:
            if star_id != SE_HAMAL:
                assert hamal_entry.magnitude < mag, (
                    f"Hamal ({hamal_entry.magnitude}) "
                    f"should be brighter than {name} ({mag})"
                )


@pytest.mark.unit
class TestZodiacalStarsNameResolution:
    """Test name resolution for zodiacal constellation stars."""

    # Aries
    def test_resolve_hamal(self):
        """Test Hamal name resolution."""
        assert resolve_star_name("Hamal") == SE_HAMAL
        assert resolve_star_name("Alpha Arietis") == SE_HAMAL
        assert resolve_star_name("Alpha Ari") == SE_HAMAL

    def test_resolve_sheratan(self):
        """Test Sheratan name resolution."""
        assert resolve_star_name("Sheratan") == SE_SHERATAN
        assert resolve_star_name("Beta Arietis") == SE_SHERATAN
        assert resolve_star_name("Beta Ari") == SE_SHERATAN

    # Cancer
    def test_resolve_tarf(self):
        """Test Tarf name resolution."""
        assert resolve_star_name("Tarf") == SE_TARF
        assert resolve_star_name("Beta Cancri") == SE_TARF
        assert resolve_star_name("Beta Cnc") == SE_TARF

    def test_resolve_asellus_stars(self):
        """Test Asellus (donkey) stars name resolution."""
        assert resolve_star_name("Asellus Borealis") == SE_ASELLUS_BOREALIS
        assert resolve_star_name("Northern Donkey") == SE_ASELLUS_BOREALIS
        assert resolve_star_name("Asellus Australis") == SE_ASELLUS_AUSTRALIS
        assert resolve_star_name("Southern Donkey") == SE_ASELLUS_AUSTRALIS

    # Sagittarius
    def test_resolve_kaus_australis(self):
        """Test Kaus Australis name resolution."""
        assert resolve_star_name("Kaus Australis") == SE_KAUS_AUSTRALIS
        assert resolve_star_name("Epsilon Sagittarii") == SE_KAUS_AUSTRALIS
        assert resolve_star_name("Epsilon Sgr") == SE_KAUS_AUSTRALIS

    def test_resolve_nunki(self):
        """Test Nunki name resolution."""
        assert resolve_star_name("Nunki") == SE_NUNKI
        assert resolve_star_name("Sigma Sagittarii") == SE_NUNKI
        assert resolve_star_name("Sigma Sgr") == SE_NUNKI

    # Capricornus
    def test_resolve_algedi(self):
        """Test Algedi name resolution."""
        assert resolve_star_name("Algedi") == SE_ALGEDI
        assert resolve_star_name("Alpha Capricorni") == SE_ALGEDI
        assert resolve_star_name("Alpha Cap") == SE_ALGEDI
        assert resolve_star_name("Giedi") == SE_ALGEDI

    def test_resolve_nashira(self):
        """Test Nashira name resolution."""
        assert resolve_star_name("Nashira") == SE_NASHIRA
        assert resolve_star_name("Gamma Capricorni") == SE_NASHIRA
        assert resolve_star_name("Fortunate One") == SE_NASHIRA

    # Aquarius
    def test_resolve_sadalsuud(self):
        """Test Sadalsuud name resolution."""
        assert resolve_star_name("Sadalsuud") == SE_SADALSUUD
        assert resolve_star_name("Beta Aquarii") == SE_SADALSUUD
        assert resolve_star_name("Beta Aqr") == SE_SADALSUUD

    def test_resolve_sadalmelik(self):
        """Test Sadalmelik name resolution."""
        assert resolve_star_name("Sadalmelik") == SE_SADALMELIK
        assert resolve_star_name("Alpha Aquarii") == SE_SADALMELIK
        assert resolve_star_name("Alpha Aqr") == SE_SADALMELIK

    # Pisces
    def test_resolve_alrescha(self):
        """Test Alrescha name resolution."""
        assert resolve_star_name("Alrescha") == SE_ALRESCHA
        assert resolve_star_name("Alpha Piscium") == SE_ALRESCHA
        assert resolve_star_name("The Knot") == SE_ALRESCHA

    def test_resolve_eta_piscium(self):
        """Test Eta Piscium name resolution."""
        assert resolve_star_name("Eta Piscium") == SE_ETA_PISCIUM
        assert resolve_star_name("Eta Psc") == SE_ETA_PISCIUM

    # Canonical names
    def test_canonical_names(self):
        """Test canonical name retrieval for new zodiacal stars."""
        assert get_canonical_star_name(SE_HAMAL) == "Hamal"
        assert get_canonical_star_name(SE_TARF) == "Tarf"
        assert get_canonical_star_name(SE_KAUS_AUSTRALIS) == "Kaus Australis"
        assert get_canonical_star_name(SE_SADALSUUD) == "Sadalsuud"
        assert get_canonical_star_name(SE_ALRESCHA) == "Alrescha"


@pytest.mark.unit
class TestZodiacalConstellationCoverage:
    """Test that all 12 zodiacal constellations have stars."""

    def test_all_zodiacal_constellations_have_stars(self):
        """Verify each zodiacal constellation has at least one star."""
        # Map constellation abbreviations to expected nomenclature patterns
        zodiacal_constellations = {
            "Aries": "Ari",
            "Taurus": "Tau",
            "Gemini": "Gem",
            "Cancer": "Cnc",
            "Leo": "Leo",
            "Virgo": "Vir",
            "Libra": "Lib",
            "Scorpius": "Sco",
            "Sagittarius": "Sgr",
            "Capricornus": "Cap",
            "Aquarius": "Aqr",
            "Pisces": "Psc",
        }

        for constellation, abbrev in zodiacal_constellations.items():
            found = False
            for entry in STAR_CATALOG:
                if abbrev in entry.nomenclature:
                    found = True
                    break
            assert found, (
                f"Zodiacal constellation {constellation} ({abbrev}) "
                "has no stars in catalog"
            )

    def test_each_zodiacal_constellation_has_minimum_stars(self):
        """Verify each zodiacal constellation has at least 2 bright stars."""
        zodiacal_constellations = {
            "Aries": ("Ari", 3),  # Hamal, Sheratan, Mesarthim
            "Taurus": ("Tau", 2),  # Aldebaran, Elnath (+ clusters)
            "Gemini": ("Gem", 2),  # Pollux, Castor
            "Cancer": ("Cnc", 4),  # Tarf, Acubens, Asellus Borealis/Australis
            "Leo": ("Leo", 4),  # Regulus, Denebola, Algieba, Zosma
            "Virgo": ("Vir", 2),  # Spica, Vindemiatrix
            "Libra": ("Lib", 2),  # Zubenelgenubi, Zubeneschamali
            "Scorpius": ("Sco", 6),  # Antares, Shaula, etc.
            "Sagittarius": ("Sgr", 5),  # Kaus Australis, Nunki, etc.
            "Capricornus": ("Cap", 4),  # Deneb Algedi, Algedi, Dabih, Nashira
            "Aquarius": ("Aqr", 3),  # Sadalsuud, Sadalmelik, Skat
            "Pisces": ("Psc", 2),  # Eta Piscium, Alrescha
        }

        for constellation, (abbrev, min_count) in zodiacal_constellations.items():
            count = sum(1 for entry in STAR_CATALOG if abbrev in entry.nomenclature)
            assert count >= min_count, (
                f"Zodiacal constellation {constellation} should have at least "
                f"{min_count} stars, found {count}"
            )
