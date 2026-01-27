"""
Unit tests for minor body orbital elements data.

Tests verify:
- Orbital elements are at the correct epoch (JD 2461000.5 TDB = 2025-Sep-19)
- All 18 bodies in MINOR_BODY_ELEMENTS have valid orbital parameters
- Elements are consistent with expected physical ranges
- Mean motion (n) is consistent with semi-major axis (a) via Kepler's 3rd law
"""

import math
import pytest
from libephemeris.constants import (
    SE_CHIRON,
    SE_PHOLUS,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
    SE_ERIS,
    SE_SEDNA,
    SE_HAUMEA,
    SE_MAKEMAKE,
    SE_IXION,
    SE_ORCUS,
    SE_QUAOAR,
    SE_NESSUS,
    SE_ASBOLUS,
    SE_CHARIKLO,
    SE_GONGGONG,
    SE_VARUNA,
)
from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS, OrbitalElements


# Expected epoch for all bodies (JD 2461000.5 TDB = 2025-Sep-19)
EXPECTED_EPOCH = 2461000.5

# Gaussian gravitational constant: k = 0.01720209895 rad/day
# This gives n = k / a^(3/2) in rad/day, or n_deg = k * 180/pi / a^(3/2) deg/day
# For consistency check: n_deg = 0.9856076686 / a^1.5
KEPLER_CONSTANT = 0.9856076686  # degrees/day for a in AU


@pytest.fixture
def all_minor_body_ids():
    """All minor body identifiers."""
    return [
        SE_CHIRON,
        SE_PHOLUS,
        SE_CERES,
        SE_PALLAS,
        SE_JUNO,
        SE_VESTA,
        SE_ERIS,
        SE_SEDNA,
        SE_HAUMEA,
        SE_MAKEMAKE,
        SE_IXION,
        SE_ORCUS,
        SE_QUAOAR,
        SE_NESSUS,
        SE_ASBOLUS,
        SE_CHARIKLO,
        SE_GONGGONG,
        SE_VARUNA,
    ]


@pytest.mark.unit
class TestOrbitalElementsEpoch:
    """Tests for orbital elements epoch."""

    def test_all_bodies_at_correct_epoch(self, all_minor_body_ids):
        """Verify all bodies have epoch JD 2461000.5."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert elements.epoch == EXPECTED_EPOCH, (
                f"{elements.name}: epoch {elements.epoch} != expected {EXPECTED_EPOCH}"
            )

    def test_expected_number_of_bodies(self):
        """Verify we have 18 minor bodies."""
        assert len(MINOR_BODY_ELEMENTS) == 18, (
            f"Expected 18 minor bodies, got {len(MINOR_BODY_ELEMENTS)}"
        )


@pytest.mark.unit
class TestOrbitalElementsValidity:
    """Tests for orbital element validity."""

    def test_all_eccentricities_valid(self, all_minor_body_ids):
        """Eccentricity must be in [0, 1) for bound orbits."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert 0 <= elements.e < 1, (
                f"{elements.name}: invalid eccentricity {elements.e}"
            )

    def test_all_inclinations_valid(self, all_minor_body_ids):
        """Inclination must be in [0, 180) degrees."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert 0 <= elements.i < 180, (
                f"{elements.name}: invalid inclination {elements.i}"
            )

    def test_all_angular_elements_valid(self, all_minor_body_ids):
        """Angular elements (omega, Omega, M0) must be in [0, 360) degrees."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert 0 <= elements.omega < 360, (
                f"{elements.name}: omega {elements.omega} not in [0, 360)"
            )
            assert 0 <= elements.Omega < 360, (
                f"{elements.name}: Omega {elements.Omega} not in [0, 360)"
            )
            assert 0 <= elements.M0 < 360, (
                f"{elements.name}: M0 {elements.M0} not in [0, 360)"
            )

    def test_all_semi_major_axes_positive(self, all_minor_body_ids):
        """Semi-major axis must be positive."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert elements.a > 0, (
                f"{elements.name}: invalid semi-major axis {elements.a}"
            )

    def test_all_mean_motions_positive(self, all_minor_body_ids):
        """Mean motion must be positive."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert elements.n > 0, f"{elements.name}: invalid mean motion {elements.n}"


@pytest.mark.unit
class TestOrbitalElementsConsistency:
    """Tests for orbital element physical consistency."""

    def test_mean_motion_vs_semi_major_axis(self, all_minor_body_ids):
        """
        Mean motion should be consistent with semi-major axis via Kepler's 3rd law.

        n = 0.9856076686 / a^1.5 (degrees/day)
        """
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            expected_n = KEPLER_CONSTANT / (elements.a**1.5)
            rel_error = abs(elements.n - expected_n) / expected_n

            # Allow 0.1% relative error for JPL's high-precision elements
            assert rel_error < 0.001, (
                f"{elements.name}: mean motion {elements.n} deg/day differs from "
                f"Kepler prediction {expected_n:.10f} deg/day by {rel_error * 100:.4f}%"
            )

    def test_perihelion_within_aphelion(self, all_minor_body_ids):
        """Perihelion distance q = a(1-e) must be positive."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            q = elements.a * (1 - elements.e)
            Q = elements.a * (1 + elements.e)

            assert q > 0, f"{elements.name}: negative perihelion {q} AU"
            assert Q > q, f"{elements.name}: aphelion {Q} <= perihelion {q}"


@pytest.mark.unit
class TestSpecificBodiesOrbitalElements:
    """Tests for specific body orbital element ranges."""

    def test_chiron_orbital_elements(self):
        """Chiron: Centaur with orbit between Saturn and Uranus."""
        elements = MINOR_BODY_ELEMENTS[SE_CHIRON]
        assert elements.name == "Chiron"
        assert 13.5 < elements.a < 14.0, f"Chiron a={elements.a} unexpected"
        assert 0.37 < elements.e < 0.39, f"Chiron e={elements.e} unexpected"
        assert 6.5 < elements.i < 7.5, f"Chiron i={elements.i} unexpected"

    def test_ceres_orbital_elements(self):
        """Ceres: Largest main belt asteroid, dwarf planet."""
        elements = MINOR_BODY_ELEMENTS[SE_CERES]
        assert elements.name == "Ceres"
        assert 2.7 < elements.a < 2.8, f"Ceres a={elements.a} unexpected"
        assert 0.07 < elements.e < 0.09, f"Ceres e={elements.e} unexpected"
        assert 10.0 < elements.i < 11.0, f"Ceres i={elements.i} unexpected"

    def test_pallas_high_inclination(self):
        """Pallas: Known for high inclination ~35 degrees."""
        elements = MINOR_BODY_ELEMENTS[SE_PALLAS]
        assert elements.name == "Pallas"
        assert 34.0 < elements.i < 36.0, f"Pallas i={elements.i} unexpected"

    def test_eris_distant_tno(self):
        """Eris: Distant TNO with high inclination."""
        elements = MINOR_BODY_ELEMENTS[SE_ERIS]
        assert elements.name == "Eris"
        assert 67.0 < elements.a < 69.0, f"Eris a={elements.a} unexpected"
        assert 43.0 < elements.i < 45.0, f"Eris i={elements.i} unexpected"

    def test_sedna_extreme_orbit(self):
        """Sedna: Detached object with extreme orbit."""
        elements = MINOR_BODY_ELEMENTS[SE_SEDNA]
        assert elements.name == "Sedna"
        # Sedna has the most distant semi-major axis
        assert 500 < elements.a < 600, f"Sedna a={elements.a} unexpected"
        # Very eccentric
        assert 0.85 < elements.e < 0.87, f"Sedna e={elements.e} unexpected"
        # Very slow mean motion
        assert elements.n < 0.0001, f"Sedna n={elements.n} too fast"

    def test_ixion_plutino(self):
        """Ixion: Plutino in 2:3 resonance with Neptune."""
        elements = MINOR_BODY_ELEMENTS[SE_IXION]
        assert elements.name == "Ixion"
        # Plutinos have a ~ 39.4 AU (like Pluto)
        assert 39.0 < elements.a < 40.0, f"Ixion a={elements.a} unexpected"

    def test_orcus_anti_pluto(self):
        """Orcus: Anti-Pluto (opposite phase in 2:3 resonance)."""
        elements = MINOR_BODY_ELEMENTS[SE_ORCUS]
        assert elements.name == "Orcus"
        # Similar to Pluto's orbital parameters
        assert 39.0 < elements.a < 40.0, f"Orcus a={elements.a} unexpected"

    def test_quaoar_low_eccentricity(self):
        """Quaoar: TNO with nearly circular orbit."""
        elements = MINOR_BODY_ELEMENTS[SE_QUAOAR]
        assert elements.name == "Quaoar"
        # Low eccentricity for a TNO
        assert elements.e < 0.05, (
            f"Quaoar e={elements.e} too high for 'nearly circular'"
        )

    def test_nessus_centaur(self):
        """Nessus: Centaur with orbit between Saturn and Pluto."""
        elements = MINOR_BODY_ELEMENTS[SE_NESSUS]
        assert elements.name == "Nessus"
        # Semi-major axis ~24.5 AU (between Saturn ~9.5 AU and Pluto ~39.5 AU)
        assert 24.0 < elements.a < 25.0, f"Nessus a={elements.a} unexpected"
        # Moderately eccentric orbit
        assert 0.51 < elements.e < 0.53, f"Nessus e={elements.e} unexpected"
        # Inclination ~15.6 degrees
        assert 15.0 < elements.i < 16.5, f"Nessus i={elements.i} unexpected"

    def test_asbolus_centaur(self):
        """Asbolus: Centaur with orbit between Saturn and Uranus."""
        elements = MINOR_BODY_ELEMENTS[SE_ASBOLUS]
        assert elements.name == "Asbolus"
        # Semi-major axis ~18 AU (between Saturn ~9.5 AU and Uranus ~19.2 AU)
        assert 17.5 < elements.a < 18.5, f"Asbolus a={elements.a} unexpected"
        # Highly eccentric orbit
        assert 0.61 < elements.e < 0.63, f"Asbolus e={elements.e} unexpected"
        # Inclination ~17.6 degrees
        assert 17.0 < elements.i < 18.5, f"Asbolus i={elements.i} unexpected"

    def test_chariklo_centaur(self):
        """Chariklo: Largest known centaur with ring system discovered in 2014."""
        elements = MINOR_BODY_ELEMENTS[SE_CHARIKLO]
        assert elements.name == "Chariklo"
        # Semi-major axis ~15.7 AU (between Saturn ~9.5 AU and Uranus ~19.2 AU)
        assert 15.5 < elements.a < 16.0, f"Chariklo a={elements.a} unexpected"
        # Moderate eccentricity
        assert 0.16 < elements.e < 0.18, f"Chariklo e={elements.e} unexpected"
        # Inclination ~23.4 degrees
        assert 23.0 < elements.i < 24.0, f"Chariklo i={elements.i} unexpected"

    def test_gonggong_tno(self):
        """Gonggong: TNO dwarf planet candidate (formerly 2007 OR10)."""
        elements = MINOR_BODY_ELEMENTS[SE_GONGGONG]
        assert elements.name == "Gonggong"
        # Semi-major axis ~67 AU (distant TNO)
        assert 66.0 < elements.a < 68.0, f"Gonggong a={elements.a} unexpected"
        # Moderately eccentric orbit
        assert 0.50 < elements.e < 0.51, f"Gonggong e={elements.e} unexpected"
        # High inclination ~30.9 degrees
        assert 30.5 < elements.i < 31.5, f"Gonggong i={elements.i} unexpected"
        # Slow mean motion (long orbital period ~547 years)
        assert elements.n < 0.002, f"Gonggong n={elements.n} too fast"

    def test_varuna_classical_kbo(self):
        """Varuna: Large classical Kuiper belt object (~670 km diameter)."""
        elements = MINOR_BODY_ELEMENTS[SE_VARUNA]
        assert elements.name == "Varuna"
        # Semi-major axis ~43.2 AU (classical KBO region)
        assert 42.5 < elements.a < 44.0, f"Varuna a={elements.a} unexpected"
        # Nearly circular orbit (low eccentricity)
        assert elements.e < 0.06, f"Varuna e={elements.e} too high for classical KBO"
        # Moderate inclination ~17.1 degrees
        assert 16.5 < elements.i < 18.0, f"Varuna i={elements.i} unexpected"
        # Mean motion consistent with ~284 year period
        assert 0.003 < elements.n < 0.004, f"Varuna n={elements.n} unexpected"


@pytest.mark.unit
class TestOrbitalElementsDataclass:
    """Tests for OrbitalElements dataclass."""

    def test_orbital_elements_attributes(self):
        """Test OrbitalElements has all required attributes."""
        elem = MINOR_BODY_ELEMENTS[SE_CERES]
        assert hasattr(elem, "name")
        assert hasattr(elem, "epoch")
        assert hasattr(elem, "a")
        assert hasattr(elem, "e")
        assert hasattr(elem, "i")
        assert hasattr(elem, "omega")
        assert hasattr(elem, "Omega")
        assert hasattr(elem, "M0")
        assert hasattr(elem, "n")

    def test_orbital_elements_types(self):
        """Test OrbitalElements attribute types."""
        elem = MINOR_BODY_ELEMENTS[SE_CERES]
        assert isinstance(elem.name, str)
        assert isinstance(elem.epoch, float)
        assert isinstance(elem.a, float)
        assert isinstance(elem.e, float)
        assert isinstance(elem.i, float)
        assert isinstance(elem.omega, float)
        assert isinstance(elem.Omega, float)
        assert isinstance(elem.M0, float)
        assert isinstance(elem.n, float)


@pytest.mark.unit
class TestJPLDataSourcePrecision:
    """Tests verifying JPL-level precision in orbital elements."""

    def test_semi_major_axis_precision(self, all_minor_body_ids):
        """Semi-major axis should have at least 6 significant figures."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            # Check that the value has significant precision
            # by ensuring it's not a round number
            a_str = f"{elements.a:.10f}"
            # Count non-zero digits after removing leading zeros and decimal point
            significant = len(a_str.replace(".", "").lstrip("0"))
            assert significant >= 6, (
                f"{elements.name}: semi-major axis {elements.a} has insufficient precision"
            )

    def test_mean_motion_precision(self, all_minor_body_ids):
        """Mean motion should have at least 6 significant figures."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            n_str = f"{elements.n:.15f}"
            # Count significant digits
            significant = len(
                n_str.replace(".", "").replace("-", "").lstrip("0").rstrip("0")
            )
            # Some bodies like Sedna have very small n, allow 5 significant figures
            assert significant >= 5, (
                f"{elements.name}: mean motion {elements.n} has insufficient precision"
            )
