"""
Unit tests for Behenian fixed stars.

The 15 Behenian fixed stars are from medieval astrology (Heinrich Cornelius Agrippa).
These stars were considered particularly powerful in magical and astrological workings.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SE_ALGOL,
    SE_ALCYONE,
    SE_ALDEBARAN,
    SE_REGULUS,
    SE_ALKAID,
    SE_ALGORAB,
    SE_SPICA_STAR,
    SE_ARCTURUS,
    SE_ALPHECCA,
    SE_ANTARES,
    SE_VEGA,
    SE_DENEB_ALGEDI,
    SE_FOMALHAUT,
    SE_DENEB,
    SE_MARKAB,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The 15 Behenian fixed stars as defined in medieval astrological tradition
BEHENIAN_STARS = [
    # (constant, name, approximate_lon_j2000 in ecliptic degrees)
    (SE_ALGOL, "Algol", 56.0),  # ~26° Taurus = 56° ecliptic
    (SE_ALCYONE, "Alcyone", 60.0),  # ~0° Gemini = 60° ecliptic (Pleiades)
    (SE_ALDEBARAN, "Aldebaran", 70.0),  # ~10° Gemini = 70° ecliptic
    (SE_REGULUS, "Regulus", 149.0),  # ~29° Leo = 149° ecliptic
    (SE_ALKAID, "Alkaid", 177.0),  # ~27° Virgo = 177° ecliptic (Eta Ursae Majoris)
    (SE_ALGORAB, "Algorab", 193.0),  # ~13° Libra = 193° ecliptic (Delta Corvi)
    (SE_SPICA_STAR, "Spica", 203.0),  # ~23° Libra = 203° ecliptic
    (SE_ARCTURUS, "Arcturus", 203.0),  # ~24° Libra = 203-204° ecliptic
    (
        SE_ALPHECCA,
        "Alphecca",
        222.0,
    ),  # ~12° Scorpio = 222° ecliptic (Alpha Coronae Borealis)
    (SE_ANTARES, "Antares", 249.0),  # ~9° Sagittarius = 249° ecliptic
    (SE_VEGA, "Vega", 275.0),  # ~15° Capricorn = 275° ecliptic
    (
        SE_DENEB_ALGEDI,
        "Deneb Algedi",
        323.0,
    ),  # ~23° Aquarius = 323° ecliptic (Delta Capricorni)
    (SE_FOMALHAUT, "Fomalhaut", 333.0),  # ~3° Pisces = 333° ecliptic
    (SE_DENEB, "Deneb", 335.0),  # ~5° Pisces = 335° ecliptic (Alpha Cygni)
    (SE_MARKAB, "Markab", 353.0),  # ~23° Pisces = 353° ecliptic (Alpha Pegasi)
]


@pytest.mark.unit
class TestBehenianStarsCatalog:
    """Test that all 15 Behenian stars are in the catalog."""

    def test_all_behenian_stars_present(self):
        """Verify all 15 Behenian stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _ in BEHENIAN_STARS:
            assert star_id in catalog_ids, (
                f"Behenian star {name} (ID={star_id}) not in catalog"
            )

    def test_behenian_count(self):
        """Verify we have exactly 15 Behenian stars defined."""
        assert len(BEHENIAN_STARS) == 15, "Should have exactly 15 Behenian stars"

    @pytest.mark.parametrize("star_id,name,_", BEHENIAN_STARS)
    def test_star_has_proper_motion(self, star_id, name, _):
        """Verify each Behenian star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,_", BEHENIAN_STARS)
    def test_star_has_valid_magnitude(self, star_id, name, _):
        """Verify each Behenian star has a valid magnitude."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert -2 < entry.magnitude < 5, (
            f"Star {name} should have reasonable magnitude, got {entry.magnitude}"
        )


@pytest.mark.unit
class TestBehenianStarsCalculation:
    """Test position calculations for Behenian stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,approx_lon", BEHENIAN_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, approx_lon):
        """Test each Behenian star returns a reasonable position."""
        pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}° out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}° out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    @pytest.mark.parametrize("star_id,name,approx_lon", BEHENIAN_STARS)
    def test_star_position_in_expected_range(
        self, standard_jd, star_id, name, approx_lon
    ):
        """Test each Behenian star is near its expected longitude."""
        pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)

        # Allow 10° tolerance for approximate positions (precession moves stars)
        diff = abs(pos[0] - approx_lon)
        if diff > 180:
            diff = 360 - diff

        assert diff < 15, f"{name} expected ~{approx_lon}° but got {pos[0]:.1f}°"


@pytest.mark.unit
class TestBehenianStarsNameResolution:
    """Test name resolution for Behenian stars."""

    def test_resolve_algol(self):
        """Test Algol name resolution."""
        assert resolve_star_name("Algol") == SE_ALGOL
        assert resolve_star_name("Demon Star") == SE_ALGOL
        assert resolve_star_name("Beta Persei") == SE_ALGOL

    def test_resolve_alcyone(self):
        """Test Alcyone (Pleiades) name resolution."""
        assert resolve_star_name("Alcyone") == SE_ALCYONE
        assert resolve_star_name("Pleiades") == SE_ALCYONE
        assert resolve_star_name("Eta Tauri") == SE_ALCYONE
        assert resolve_star_name("Seven Sisters") == SE_ALCYONE

    def test_resolve_aldebaran(self):
        """Test Aldebaran name resolution."""
        assert resolve_star_name("Aldebaran") == SE_ALDEBARAN
        assert resolve_star_name("Alpha Tauri") == SE_ALDEBARAN
        assert resolve_star_name("Eye of Taurus") == SE_ALDEBARAN

    def test_resolve_regulus(self):
        """Test Regulus name resolution."""
        assert resolve_star_name("Regulus") == SE_REGULUS
        assert resolve_star_name("Alpha Leonis") == SE_REGULUS
        assert resolve_star_name("Cor Leonis") == SE_REGULUS

    def test_resolve_alkaid(self):
        """Test Alkaid name resolution."""
        assert resolve_star_name("Alkaid") == SE_ALKAID
        assert resolve_star_name("Eta Ursae Majoris") == SE_ALKAID
        assert resolve_star_name("Benetnash") == SE_ALKAID

    def test_resolve_algorab(self):
        """Test Algorab name resolution."""
        assert resolve_star_name("Algorab") == SE_ALGORAB
        assert resolve_star_name("Delta Corvi") == SE_ALGORAB
        assert resolve_star_name("Crow's Wing") == SE_ALGORAB

    def test_resolve_spica(self):
        """Test Spica name resolution."""
        assert resolve_star_name("Spica") == SE_SPICA_STAR
        assert resolve_star_name("Alpha Virginis") == SE_SPICA_STAR

    def test_resolve_arcturus(self):
        """Test Arcturus name resolution."""
        assert resolve_star_name("Arcturus") == SE_ARCTURUS
        assert resolve_star_name("Alpha Bootis") == SE_ARCTURUS
        assert resolve_star_name("Bear Guard") == SE_ARCTURUS

    def test_resolve_alphecca(self):
        """Test Alphecca name resolution."""
        assert resolve_star_name("Alphecca") == SE_ALPHECCA
        assert resolve_star_name("Alpha Coronae Borealis") == SE_ALPHECCA
        assert resolve_star_name("Gemma") == SE_ALPHECCA

    def test_resolve_antares(self):
        """Test Antares name resolution."""
        assert resolve_star_name("Antares") == SE_ANTARES
        assert resolve_star_name("Alpha Scorpii") == SE_ANTARES
        assert resolve_star_name("Rival of Mars") == SE_ANTARES

    def test_resolve_vega(self):
        """Test Vega name resolution."""
        assert resolve_star_name("Vega") == SE_VEGA
        assert resolve_star_name("Alpha Lyrae") == SE_VEGA
        assert resolve_star_name("Harp Star") == SE_VEGA

    def test_resolve_deneb_algedi(self):
        """Test Deneb Algedi name resolution."""
        assert resolve_star_name("Deneb Algedi") == SE_DENEB_ALGEDI
        assert resolve_star_name("Delta Capricorni") == SE_DENEB_ALGEDI
        assert resolve_star_name("Tail of the Goat") == SE_DENEB_ALGEDI

    def test_resolve_fomalhaut(self):
        """Test Fomalhaut name resolution."""
        assert resolve_star_name("Fomalhaut") == SE_FOMALHAUT
        assert resolve_star_name("Alpha Piscis Austrini") == SE_FOMALHAUT
        assert resolve_star_name("Fish's Mouth") == SE_FOMALHAUT

    def test_resolve_deneb(self):
        """Test Deneb name resolution."""
        assert resolve_star_name("Deneb") == SE_DENEB
        assert resolve_star_name("Alpha Cygni") == SE_DENEB
        assert resolve_star_name("Tail of Hen") == SE_DENEB

    def test_resolve_markab(self):
        """Test Markab name resolution."""
        assert resolve_star_name("Markab") == SE_MARKAB
        assert resolve_star_name("Alpha Pegasi") == SE_MARKAB
        assert resolve_star_name("Saddle") == SE_MARKAB

    @pytest.mark.parametrize("star_id,name,_", BEHENIAN_STARS)
    def test_canonical_name_retrieval(self, star_id, name, _):
        """Test canonical name retrieval for all Behenian stars."""
        canonical = get_canonical_star_name(star_id)
        assert canonical == name, f"Expected {name}, got {canonical}"


@pytest.mark.unit
class TestBehenianStarsProperMotion:
    """Test proper motion effects for Behenian stars over time."""

    def test_proper_motion_over_50_years(self):
        """Test that Behenian stars move over 50 years."""
        jd1 = ephem.swe_julday(2000, 1, 1, 12.0)
        jd2 = ephem.swe_julday(2050, 1, 1, 12.0)

        # Test a few high proper motion stars
        test_stars = [
            (SE_ARCTURUS, "Arcturus"),  # High proper motion
            (SE_VEGA, "Vega"),  # Moderate proper motion
            (SE_REGULUS, "Regulus"),  # Low proper motion
        ]

        for star_id, name in test_stars:
            pos1, _ = ephem.swe_calc_ut(jd1, star_id, 0)
            pos2, _ = ephem.swe_calc_ut(jd2, star_id, 0)

            diff = abs(pos2[0] - pos1[0])
            if diff > 180:
                diff = 360 - diff

            # Should move at least 0.01° over 50 years (precession alone is ~0.7°)
            assert diff > 0.01, f"{name} should show movement over 50 years"


@pytest.mark.unit
class TestNewBehenianStarsData:
    """Test the newly added Behenian stars have correct data."""

    def test_alcyone_data(self):
        """Test Alcyone (Pleiades) catalog entry."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == SE_ALCYONE:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Alcyone"
        assert entry.nomenclature == "etTau"
        assert entry.hip_number == 17702
        assert 56 < entry.data.ra_j2000 < 57  # ~56.87°
        assert 24 < entry.data.dec_j2000 < 25  # ~24.1°
        assert 2.5 < entry.magnitude < 3.0  # ~2.87

    def test_algorab_data(self):
        """Test Algorab (Delta Corvi) catalog entry."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == SE_ALGORAB:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Algorab"
        assert entry.nomenclature == "deCrv"
        assert entry.hip_number == 60965
        assert 187 < entry.data.ra_j2000 < 188  # ~187.47°
        assert -17 < entry.data.dec_j2000 < -16  # ~-16.52°
        assert 2.5 < entry.magnitude < 3.5  # ~2.95

    def test_alphecca_data(self):
        """Test Alphecca (Alpha Coronae Borealis) catalog entry."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == SE_ALPHECCA:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Alphecca"
        assert entry.nomenclature == "alCrB"
        assert entry.hip_number == 76267
        assert 233 < entry.data.ra_j2000 < 234  # ~233.67°
        assert 26 < entry.data.dec_j2000 < 27  # ~26.71°
        assert 2.0 < entry.magnitude < 2.5  # ~2.23

    def test_deneb_algedi_data(self):
        """Test Deneb Algedi (Delta Capricorni) catalog entry."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == SE_DENEB_ALGEDI:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Deneb Algedi"
        assert entry.nomenclature == "deCap"
        assert entry.hip_number == 107556
        assert 326 < entry.data.ra_j2000 < 327  # ~326.76°
        assert -17 < entry.data.dec_j2000 < -16  # ~-16.13°
        assert 2.5 < entry.magnitude < 3.0  # ~2.81
