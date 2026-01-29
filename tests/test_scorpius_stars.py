"""
Unit tests for Scorpius constellation stars.

Scorpius is one of the zodiac constellations, featuring these 6 major stars:
- Antares (Alpha Sco) - Red supergiant, heart of the scorpion (Royal Star)
- Shaula (Lambda Sco) - Scorpion's sting
- Sargas (Theta Sco) - Scorpion's tail
- Dschubba (Delta Sco) - Scorpion's forehead
- Graffias (Beta Sco) - Scorpion's claw
- Lesath (Upsilon Sco) - Stinger
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SE_ANTARES,
    SE_SHAULA,
    SE_SARGAS,
    SE_DSCHUBBA,
    SE_GRAFFIAS,
    SE_LESATH,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The 6 major Scorpius stars with their properties
# (constant, name, hip_number, magnitude)
SCORPIUS_STARS = [
    (SE_ANTARES, "Antares", 80763, 1.06),  # Alpha Sco - red supergiant, heart
    (SE_SHAULA, "Shaula", 85927, 1.62),  # Lambda Sco - sting
    (SE_SARGAS, "Sargas", 86228, 1.87),  # Theta Sco - tail
    (SE_DSCHUBBA, "Dschubba", 78401, 2.32),  # Delta Sco - forehead
    (SE_GRAFFIAS, "Graffias", 78820, 2.56),  # Beta Sco - claw
    (SE_LESATH, "Lesath", 85696, 2.70),  # Upsilon Sco - stinger
]

# The sting stars (Shaula and Lesath form a close pair)
SCORPIUS_STING_STARS = [
    (SE_SHAULA, "Shaula"),
    (SE_LESATH, "Lesath"),
]


@pytest.mark.unit
class TestScorpiusStarsCatalog:
    """Test that all 6 major Scorpius stars are in the catalog."""

    def test_all_scorpius_stars_present(self):
        """Verify all 6 Scorpius stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in SCORPIUS_STARS:
            assert star_id in catalog_ids, (
                f"Scorpius star {name} (ID={star_id}) not in catalog"
            )

    def test_scorpius_star_count(self):
        """Verify we have exactly 6 major Scorpius stars defined."""
        assert len(SCORPIUS_STARS) == 6, "Should have exactly 6 major Scorpius stars"

    def test_antares_is_brightest(self):
        """Verify Antares is the brightest Scorpius star (lowest magnitude)."""
        antares_mag = None
        other_mags = []

        for star_id, name, _, mag in SCORPIUS_STARS:
            if star_id == SE_ANTARES:
                antares_mag = mag
            else:
                other_mags.append((name, mag))

        assert antares_mag is not None, "Antares not found"
        for name, mag in other_mags:
            assert antares_mag < mag, (
                f"Antares ({antares_mag}) should be brighter than {name} ({mag})"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", SCORPIUS_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each Scorpius star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Scorpius stars should have proper motion data
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", SCORPIUS_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each Scorpius star has correct Hipparcos number."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.hip_number == hip, (
            f"Star {name} should have HIP {hip}, got {entry.hip_number}"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", SCORPIUS_STARS)
    def test_star_nomenclature_is_scorpius(self, star_id, name, hip, mag):
        """Verify each Scorpius star has Sco in nomenclature."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert "Sco" in entry.nomenclature, (
            f"Star {name} nomenclature should contain 'Sco', got {entry.nomenclature}"
        )


@pytest.mark.unit
class TestScorpiusStarsCalculation:
    """Test position calculations for Scorpius stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", SCORPIUS_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each Scorpius star returns a reasonable position."""
        pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_scorpius_sting_stars_close(self, standard_jd):
        """Test that Shaula and Lesath (sting stars) are close together."""
        positions = []
        for star_id, name in SCORPIUS_STING_STARS:
            pos, _ = ephem.swe_calc_ut(standard_jd, star_id, 0)
            positions.append((name, pos[0], pos[1]))

        # Sting stars should be within ~2 degrees of each other
        lons = [p[1] for p in positions]
        lon_diff = abs(lons[0] - lons[1])
        assert lon_diff < 3, f"Sting stars longitude spread {lon_diff}deg too large"

        lats = [p[2] for p in positions]
        lat_diff = abs(lats[0] - lats[1])
        assert lat_diff < 3, f"Sting stars latitude spread {lat_diff}deg too large"

    def test_scorpius_in_sagittarius_region(self, standard_jd):
        """Test that Scorpius stars are in the Sagittarius ecliptic region (~240-270 deg)."""
        # Antares is a good reference for Scorpius position
        pos, _ = ephem.swe_calc_ut(standard_jd, SE_ANTARES, 0)

        # Should be near 249 degrees ecliptic (Sagittarius)
        assert 240 < pos[0] < 260, (
            f"Antares should be near 249 deg ecliptic, got {pos[0]:.1f}"
        )

    def test_antares_below_ecliptic(self, standard_jd):
        """Test that Antares is located below the ecliptic (southern latitude)."""
        pos, _ = ephem.swe_calc_ut(standard_jd, SE_ANTARES, 0)

        # Antares should have negative ecliptic latitude (below ecliptic)
        assert pos[1] < 0, f"Antares should be below ecliptic, got lat {pos[1]:.2f}"

        # Should be around -4 to -5 degrees
        assert -10 < pos[1] < 0, f"Antares latitude {pos[1]:.2f} out of expected range"


@pytest.mark.unit
class TestScorpiusStarsNameResolution:
    """Test name resolution for Scorpius stars."""

    def test_resolve_antares(self):
        """Test Antares name resolution."""
        assert resolve_star_name("Antares") == SE_ANTARES
        assert resolve_star_name("Alpha Scorpii") == SE_ANTARES
        assert resolve_star_name("Alpha Sco") == SE_ANTARES
        assert resolve_star_name("Rival of Mars") == SE_ANTARES
        assert resolve_star_name("Heart of Scorpion") == SE_ANTARES

    def test_resolve_shaula(self):
        """Test Shaula name resolution."""
        assert resolve_star_name("Shaula") == SE_SHAULA
        assert resolve_star_name("Lambda Scorpii") == SE_SHAULA
        assert resolve_star_name("Lambda Sco") == SE_SHAULA
        assert resolve_star_name("Scorpion's Sting") == SE_SHAULA

    def test_resolve_sargas(self):
        """Test Sargas name resolution."""
        assert resolve_star_name("Sargas") == SE_SARGAS
        assert resolve_star_name("Theta Scorpii") == SE_SARGAS
        assert resolve_star_name("Theta Sco") == SE_SARGAS
        assert resolve_star_name("Girtab") == SE_SARGAS
        assert resolve_star_name("Scorpion's Tail") == SE_SARGAS

    def test_resolve_dschubba(self):
        """Test Dschubba name resolution."""
        assert resolve_star_name("Dschubba") == SE_DSCHUBBA
        assert resolve_star_name("Delta Scorpii") == SE_DSCHUBBA
        assert resolve_star_name("Delta Sco") == SE_DSCHUBBA
        assert resolve_star_name("Scorpion's Forehead") == SE_DSCHUBBA

    def test_resolve_graffias(self):
        """Test Graffias name resolution."""
        assert resolve_star_name("Graffias") == SE_GRAFFIAS
        assert resolve_star_name("Beta Scorpii") == SE_GRAFFIAS
        assert resolve_star_name("Beta Sco") == SE_GRAFFIAS
        assert resolve_star_name("Acrab") == SE_GRAFFIAS
        assert resolve_star_name("Akrab") == SE_GRAFFIAS

    def test_resolve_lesath(self):
        """Test Lesath name resolution."""
        assert resolve_star_name("Lesath") == SE_LESATH
        assert resolve_star_name("Upsilon Scorpii") == SE_LESATH
        assert resolve_star_name("Upsilon Sco") == SE_LESATH
        assert resolve_star_name("Stinger") == SE_LESATH

    def test_canonical_names(self):
        """Test canonical name retrieval for Scorpius stars."""
        assert get_canonical_star_name(SE_ANTARES) == "Antares"
        assert get_canonical_star_name(SE_SHAULA) == "Shaula"
        assert get_canonical_star_name(SE_SARGAS) == "Sargas"
        assert get_canonical_star_name(SE_DSCHUBBA) == "Dschubba"
        assert get_canonical_star_name(SE_GRAFFIAS) == "Graffias"
        assert get_canonical_star_name(SE_LESATH) == "Lesath"
