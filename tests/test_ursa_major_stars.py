"""
Unit tests for Ursa Major (Big Dipper) constellation stars.

The Big Dipper is one of the most recognizable asterisms in the northern sky,
consisting of 7 bright stars that form a distinctive ladle/plough shape:

Bowl stars:
- Dubhe (Alpha UMa) - pointer star, upper lip of bowl
- Merak (Beta UMa) - pointer star, lower lip of bowl
- Phecda (Gamma UMa) - bottom corner of bowl
- Megrez (Delta UMa) - connects bowl to handle

Handle stars:
- Alioth (Epsilon UMa) - brightest of the seven
- Mizar (Zeta UMa) - famous double star, "horse"
- Alcor (80 UMa) - optical double with Mizar, "rider"
- Alkaid (Eta UMa) - end of handle

The two "pointer" stars (Dubhe and Merak) point to Polaris, the North Star.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SE_DUBHE,
    SE_MERAK,
    SE_PHECDA,
    SE_MEGREZ,
    SE_ALIOTH,
    SE_MIZAR,
    SE_ALCOR,
    SE_ALKAID,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The 7 Big Dipper stars plus Alcor (visual companion to Mizar)
# (constant, name, hip_number, magnitude)
BIG_DIPPER_STARS = [
    (SE_DUBHE, "Dubhe", 54061, 1.79),  # Alpha UMa - pointer, upper
    (SE_MERAK, "Merak", 53910, 2.37),  # Beta UMa - pointer, lower
    (SE_PHECDA, "Phecda", 58001, 2.44),  # Gamma UMa - bowl corner
    (SE_MEGREZ, "Megrez", 59774, 3.31),  # Delta UMa - bowl-handle junction
    (SE_ALIOTH, "Alioth", 62956, 1.77),  # Epsilon UMa - brightest
    (SE_MIZAR, "Mizar", 65378, 2.23),  # Zeta UMa - double star
    (SE_ALCOR, "Alcor", 65477, 3.99),  # 80 UMa - companion to Mizar
    (SE_ALKAID, "Alkaid", 67301, 1.86),  # Eta UMa - end of handle
]

# The bowl stars (without handle)
BOWL_STARS = [
    (SE_DUBHE, "Dubhe"),
    (SE_MERAK, "Merak"),
    (SE_PHECDA, "Phecda"),
    (SE_MEGREZ, "Megrez"),
]

# The handle stars
HANDLE_STARS = [
    (SE_MEGREZ, "Megrez"),  # junction point
    (SE_ALIOTH, "Alioth"),
    (SE_MIZAR, "Mizar"),
    (SE_ALKAID, "Alkaid"),
]

# Pointer stars (point to Polaris)
POINTER_STARS = [
    (SE_MERAK, "Merak"),
    (SE_DUBHE, "Dubhe"),
]


@pytest.mark.unit
class TestBigDipperStarsCatalog:
    """Test that all Big Dipper stars are in the catalog."""

    def test_all_big_dipper_stars_present(self):
        """Verify all 8 Big Dipper stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in BIG_DIPPER_STARS:
            assert star_id in catalog_ids, (
                f"Big Dipper star {name} (ID={star_id}) not in catalog"
            )

    def test_big_dipper_star_count(self):
        """Verify we have exactly 8 Big Dipper stars defined (7 + Alcor)."""
        assert len(BIG_DIPPER_STARS) == 8, (
            "Should have exactly 8 Big Dipper stars (including Alcor)"
        )

    def test_alioth_is_brightest(self):
        """Verify Alioth is the brightest Big Dipper star (lowest magnitude)."""
        alioth_mag = None
        other_mags = []

        for star_id, name, _, mag in BIG_DIPPER_STARS:
            if star_id == SE_ALIOTH:
                alioth_mag = mag
            else:
                other_mags.append((name, mag))

        assert alioth_mag is not None, "Alioth not found"
        for name, mag in other_mags:
            assert alioth_mag < mag, (
                f"Alioth ({alioth_mag}) should be brighter than {name} ({mag})"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", BIG_DIPPER_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each Big Dipper star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Big Dipper stars should have proper motion data
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", BIG_DIPPER_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each Big Dipper star has correct Hipparcos number."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.hip_number == hip, (
            f"Star {name} should have HIP {hip}, got {entry.hip_number}"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", BIG_DIPPER_STARS)
    def test_star_nomenclature_is_ursa_major(self, star_id, name, hip, mag):
        """Verify each Big Dipper star has UMa in nomenclature."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert "UMa" in entry.nomenclature, (
            f"Star {name} nomenclature should contain 'UMa', got {entry.nomenclature}"
        )


@pytest.mark.unit
class TestBigDipperStarsCalculation:
    """Test position calculations for Big Dipper stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", BIG_DIPPER_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each Big Dipper star returns a reasonable position."""
        pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_bowl_stars_form_quadrilateral(self, standard_jd):
        """Test that the bowl stars form a rough quadrilateral."""
        positions = []
        for star_id, name in BOWL_STARS:
            pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)
            positions.append((name, pos[0], pos[1]))

        # Bowl stars should be within a reasonable region
        lons = [p[1] for p in positions]
        lon_range = max(lons) - min(lons)
        assert lon_range < 20, f"Bowl longitude spread {lon_range}deg too large"

        lats = [p[2] for p in positions]
        lat_range = max(lats) - min(lats)
        assert lat_range < 15, f"Bowl latitude spread {lat_range}deg too large"

    def test_handle_stars_roughly_aligned(self, standard_jd):
        """Test that handle stars form a rough line."""
        positions = []
        for star_id, name in HANDLE_STARS:
            pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)
            positions.append((name, pos[0], pos[1]))

        # Handle stars should span along an arc
        lons = [p[1] for p in positions]
        lon_range = max(lons) - min(lons)
        # Handle spans from Megrez to Alkaid, about 20-25 degrees
        assert 10 < lon_range < 35, (
            f"Handle longitude spread {lon_range}deg not in expected range"
        )

    def test_mizar_alcor_close_together(self, standard_jd):
        """Test that Mizar and Alcor are very close (famous naked-eye double)."""
        mizar_pos, _ = ephem.swe_calc_ut(standard_jd, SE_MIZAR, 0)
        alcor_pos, _ = ephem.swe_calc_ut(standard_jd, SE_ALCOR, 0)

        # Calculate angular separation (simplified for small angles)
        lon_diff = abs(mizar_pos[0] - alcor_pos[0])
        lat_diff = abs(mizar_pos[1] - alcor_pos[1])
        separation = (lon_diff**2 + lat_diff**2) ** 0.5

        # Mizar and Alcor are about 12 arcminutes apart = 0.2 degrees
        assert separation < 0.5, (
            f"Mizar-Alcor separation {separation:.3f}deg too large, should be ~0.2deg"
        )

    def test_big_dipper_in_leo_virgo_region(self, standard_jd):
        """Test that Big Dipper stars are in the Leo/Virgo ecliptic region (~120-180 deg)."""
        # Dubhe is a good reference for the Big Dipper's position
        pos, _ = ephem.swe_calc_ut(standard_jd, SE_DUBHE, 0)

        # Should be in the Leo/Virgo region
        assert 130 < pos[0] < 180, (
            f"Dubhe should be in Leo/Virgo region, got {pos[0]:.1f} degrees"
        )

    def test_phecda_position(self, standard_jd):
        """Test that Phecda (new star) is positioned correctly in the bowl."""
        pos, _ = ephem.swe_calc_ut(standard_jd, SE_PHECDA, 0)

        # Phecda should be near 147 degrees ecliptic longitude (Leo region)
        assert 140 < pos[0] < 160, (
            f"Phecda longitude {pos[0]:.1f} out of expected range"
        )

        # Phecda has positive ecliptic latitude (north of ecliptic)
        assert 40 < pos[1] < 60, f"Phecda latitude {pos[1]:.2f} out of expected range"

    def test_megrez_position(self, standard_jd):
        """Test that Megrez (new star) is positioned correctly at bowl-handle junction."""
        pos, _ = ephem.swe_calc_ut(standard_jd, SE_MEGREZ, 0)

        # Megrez should be near 152 degrees ecliptic longitude
        assert 145 < pos[0] < 165, (
            f"Megrez longitude {pos[0]:.1f} out of expected range"
        )

        # Megrez has positive ecliptic latitude
        assert 45 < pos[1] < 65, f"Megrez latitude {pos[1]:.2f} out of expected range"


@pytest.mark.unit
class TestBigDipperStarsNameResolution:
    """Test name resolution for Big Dipper stars."""

    def test_resolve_dubhe(self):
        """Test Dubhe name resolution."""
        assert resolve_star_name("Dubhe") == SE_DUBHE
        assert resolve_star_name("Alpha Ursae Majoris") == SE_DUBHE
        assert resolve_star_name("Alpha UMa") == SE_DUBHE

    def test_resolve_merak(self):
        """Test Merak name resolution."""
        assert resolve_star_name("Merak") == SE_MERAK
        assert resolve_star_name("Beta Ursae Majoris") == SE_MERAK
        assert resolve_star_name("Beta UMa") == SE_MERAK

    def test_resolve_phecda(self):
        """Test Phecda name resolution."""
        assert resolve_star_name("Phecda") == SE_PHECDA
        assert resolve_star_name("Gamma Ursae Majoris") == SE_PHECDA
        assert resolve_star_name("Gamma UMa") == SE_PHECDA
        assert resolve_star_name("Phad") == SE_PHECDA
        assert resolve_star_name("Phekda") == SE_PHECDA

    def test_resolve_megrez(self):
        """Test Megrez name resolution."""
        assert resolve_star_name("Megrez") == SE_MEGREZ
        assert resolve_star_name("Delta Ursae Majoris") == SE_MEGREZ
        assert resolve_star_name("Delta UMa") == SE_MEGREZ
        assert resolve_star_name("Kaffa") == SE_MEGREZ

    def test_resolve_alioth(self):
        """Test Alioth name resolution."""
        assert resolve_star_name("Alioth") == SE_ALIOTH
        assert resolve_star_name("Epsilon Ursae Majoris") == SE_ALIOTH
        assert resolve_star_name("Epsilon UMa") == SE_ALIOTH

    def test_resolve_mizar(self):
        """Test Mizar name resolution."""
        assert resolve_star_name("Mizar") == SE_MIZAR
        assert resolve_star_name("Zeta Ursae Majoris") == SE_MIZAR
        assert resolve_star_name("Zeta UMa") == SE_MIZAR
        assert resolve_star_name("Horse and Rider") == SE_MIZAR

    def test_resolve_alcor(self):
        """Test Alcor name resolution."""
        assert resolve_star_name("Alcor") == SE_ALCOR
        assert resolve_star_name("80 Ursae Majoris") == SE_ALCOR
        assert resolve_star_name("80 UMa") == SE_ALCOR
        assert resolve_star_name("Suha") == SE_ALCOR

    def test_resolve_alkaid(self):
        """Test Alkaid name resolution."""
        assert resolve_star_name("Alkaid") == SE_ALKAID
        assert resolve_star_name("Eta Ursae Majoris") == SE_ALKAID
        assert resolve_star_name("Eta UMa") == SE_ALKAID
        assert resolve_star_name("Benetnash") == SE_ALKAID

    def test_canonical_names(self):
        """Test canonical name retrieval for Big Dipper stars."""
        assert get_canonical_star_name(SE_DUBHE) == "Dubhe"
        assert get_canonical_star_name(SE_MERAK) == "Merak"
        assert get_canonical_star_name(SE_PHECDA) == "Phecda"
        assert get_canonical_star_name(SE_MEGREZ) == "Megrez"
        assert get_canonical_star_name(SE_ALIOTH) == "Alioth"
        assert get_canonical_star_name(SE_MIZAR) == "Mizar"
        assert get_canonical_star_name(SE_ALCOR) == "Alcor"
        assert get_canonical_star_name(SE_ALKAID) == "Alkaid"
