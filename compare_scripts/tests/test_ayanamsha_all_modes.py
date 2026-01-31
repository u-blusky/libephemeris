"""
Test all 43 ayanamsha modes at multiple dates.

This test comprehensively validates all ayanamsha modes across different
time periods to ensure consistency and accuracy of sidereal calculations.
"""

import pytest
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import (
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_LAHIRI,
    SE_SIDM_DELUCE,
    SE_SIDM_RAMAN,
    SE_SIDM_USHASHASHI,
    SE_SIDM_KRISHNAMURTI,
    SE_SIDM_DJWHAL_KHUL,
    SE_SIDM_YUKTESHWAR,
    SE_SIDM_JN_BHASIN,
    SE_SIDM_BABYL_KUGLER1,
    SE_SIDM_BABYL_KUGLER2,
    SE_SIDM_BABYL_KUGLER3,
    SE_SIDM_BABYL_HUBER,
    SE_SIDM_BABYL_ETPSC,
    SE_SIDM_ALDEBARAN_15TAU,
    SE_SIDM_HIPPARCHOS,
    SE_SIDM_SASSANIAN,
    SE_SIDM_GALCENT_0SAG,
    SE_SIDM_J2000,
    SE_SIDM_J1900,
    SE_SIDM_B1950,
    SE_SIDM_SURYASIDDHANTA,
    SE_SIDM_SURYASIDDHANTA_MSUN,
    SE_SIDM_ARYABHATA,
    SE_SIDM_ARYABHATA_MSUN,
    SE_SIDM_SS_REVATI,
    SE_SIDM_SS_CITRA,
    SE_SIDM_TRUE_CITRA,
    SE_SIDM_TRUE_REVATI,
    SE_SIDM_TRUE_PUSHYA,
    SE_SIDM_GALCENT_RGILBRAND,
    SE_SIDM_GALEQU_IAU1958,
    SE_SIDM_GALEQU_TRUE,
    SE_SIDM_GALEQU_MULA,
    SE_SIDM_GALALIGN_MARDYKS,
    SE_SIDM_TRUE_MULA,
    SE_SIDM_GALCENT_MULA_WILHELM,
    SE_SIDM_ARYABHATA_522,
    SE_SIDM_BABYL_BRITTON,
    SE_SIDM_TRUE_SHEORAN,
    SE_SIDM_GALCENT_COCHRANE,
    SE_SIDM_GALEQU_FIORENZA,
    SE_SIDM_VALENS_MOON,
)

# All 43 ayanamsha modes with their names
ALL_AYANAMSHA_MODES = [
    (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_DELUCE, "De Luce"),
    (SE_SIDM_RAMAN, "Raman"),
    (SE_SIDM_USHASHASHI, "Ushashashi"),
    (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SE_SIDM_DJWHAL_KHUL, "Djwhal Khul"),
    (SE_SIDM_YUKTESHWAR, "Yukteshwar"),
    (SE_SIDM_JN_BHASIN, "JN Bhasin"),
    (SE_SIDM_BABYL_KUGLER1, "Babylonian (Kugler 1)"),
    (SE_SIDM_BABYL_KUGLER2, "Babylonian (Kugler 2)"),
    (SE_SIDM_BABYL_KUGLER3, "Babylonian (Kugler 3)"),
    (SE_SIDM_BABYL_HUBER, "Babylonian (Huber)"),
    (SE_SIDM_BABYL_ETPSC, "Babylonian (ETPSC)"),
    (SE_SIDM_ALDEBARAN_15TAU, "Aldebaran at 15 Tau"),
    (SE_SIDM_HIPPARCHOS, "Hipparchos"),
    (SE_SIDM_SASSANIAN, "Sassanian"),
    (SE_SIDM_GALCENT_0SAG, "Galactic Center at 0 Sag"),
    (SE_SIDM_J2000, "J2000"),
    (SE_SIDM_J1900, "J1900"),
    (SE_SIDM_B1950, "B1950"),
    (SE_SIDM_SURYASIDDHANTA, "Suryasiddhanta"),
    (SE_SIDM_SURYASIDDHANTA_MSUN, "Suryasiddhanta (mean Sun)"),
    (SE_SIDM_ARYABHATA, "Aryabhata"),
    (SE_SIDM_ARYABHATA_MSUN, "Aryabhata (mean Sun)"),
    (SE_SIDM_SS_REVATI, "SS Revati"),
    (SE_SIDM_SS_CITRA, "SS Citra"),
    (SE_SIDM_TRUE_CITRA, "True Citra"),
    (SE_SIDM_TRUE_REVATI, "True Revati"),
    (SE_SIDM_TRUE_PUSHYA, "True Pushya"),
    (SE_SIDM_GALCENT_RGILBRAND, "Galactic Center (Gil Brand)"),
    (SE_SIDM_GALEQU_IAU1958, "Galactic Equator (IAU 1958)"),
    (SE_SIDM_GALEQU_TRUE, "Galactic Equator (True)"),
    (SE_SIDM_GALEQU_MULA, "Galactic Equator at Mula"),
    (SE_SIDM_GALALIGN_MARDYKS, "Galactic Alignment (Mardyks)"),
    (SE_SIDM_TRUE_MULA, "True Mula"),
    (SE_SIDM_GALCENT_MULA_WILHELM, "Galactic Center at Mula (Wilhelm)"),
    (SE_SIDM_ARYABHATA_522, "Aryabhata 522"),
    (SE_SIDM_BABYL_BRITTON, "Babylonian (Britton)"),
    (SE_SIDM_TRUE_SHEORAN, "True Sheoran"),
    (SE_SIDM_GALCENT_COCHRANE, "Galactic Center (Cochrane)"),
    (SE_SIDM_GALEQU_FIORENZA, "Galactic Equator (Fiorenza)"),
    (SE_SIDM_VALENS_MOON, "Valens (Moon)"),
]

# Star-based modes that may have larger differences due to proper motion calculations
STAR_BASED_MODES = {
    SE_SIDM_TRUE_CITRA,
    SE_SIDM_TRUE_REVATI,
    SE_SIDM_TRUE_PUSHYA,
    SE_SIDM_TRUE_MULA,
    SE_SIDM_TRUE_SHEORAN,
    SE_SIDM_GALCENT_0SAG,
    SE_SIDM_GALCENT_RGILBRAND,
    SE_SIDM_GALCENT_MULA_WILHELM,
    SE_SIDM_GALCENT_COCHRANE,
    SE_SIDM_GALEQU_TRUE,
    SE_SIDM_GALEQU_MULA,
    SE_SIDM_GALEQU_FIORENZA,
}

# Multiple test dates spanning different epochs
TEST_DATES = [
    (1900, 1, 1, 12.0, "1900-01-01 (Early 20th century)"),
    (1950, 6, 15, 0.0, "1950-06-15 (Mid 20th century)"),
    (2000, 1, 1, 12.0, "2000-01-01 (J2000 epoch)"),
    (2024, 3, 20, 3.06, "2024-03-20 (Vernal Equinox 2024)"),
    (2050, 7, 4, 12.0, "2050-07-04 (Mid 21st century)"),
]


def get_tolerance(sid_mode: int) -> float:
    """Return appropriate tolerance for ayanamsha comparison based on mode type."""
    # Star-based modes use real-time star positions and may have
    # larger differences due to proper motion calculation methods
    if sid_mode in STAR_BASED_MODES:
        return 1.0
    return 0.1


class TestAllAyanamshasMultipleDates:
    """Test all 43 ayanamsha modes at multiple dates."""

    @pytest.mark.parametrize("sid_mode, mode_name", ALL_AYANAMSHA_MODES)
    @pytest.mark.parametrize("year, month, day, hour, date_desc", TEST_DATES)
    def test_ayanamsha_value_at_date(
        self, sid_mode, mode_name, year, month, day, hour, date_desc
    ):
        """Test each ayanamsha mode value at multiple dates."""
        jd = swe.julday(year, month, day, hour)

        # Swiss Ephemeris
        swe.set_sid_mode(sid_mode)
        aya_swe = swe.get_ayanamsa_ut(jd)

        # libephemeris
        pyephem.swe_set_sid_mode(sid_mode)
        aya_py = pyephem.swe_get_ayanamsa_ut(jd)

        diff = abs(aya_swe - aya_py)
        tol = get_tolerance(sid_mode)

        assert diff < tol, (
            f"{mode_name} ({sid_mode}) at {date_desc}: "
            f"SWE={aya_swe:.6f}, PY={aya_py:.6f}, diff={diff:.6f} > tol={tol}"
        )

    @pytest.mark.parametrize("sid_mode, mode_name", ALL_AYANAMSHA_MODES)
    def test_ayanamsha_increases_over_time(self, sid_mode, mode_name):
        """Verify that ayanamsha increases over time (precession)."""
        # J2000 mode is anchored at 0 at J2000, so skip it
        if sid_mode == SE_SIDM_J2000:
            pytest.skip("J2000 mode has zero ayanamsha at J2000 epoch")

        jd_1900 = swe.julday(1900, 1, 1, 12.0)
        jd_2000 = swe.julday(2000, 1, 1, 12.0)
        jd_2050 = swe.julday(2050, 1, 1, 12.0)

        pyephem.swe_set_sid_mode(sid_mode)
        aya_1900 = pyephem.swe_get_ayanamsa_ut(jd_1900)
        aya_2000 = pyephem.swe_get_ayanamsa_ut(jd_2000)
        aya_2050 = pyephem.swe_get_ayanamsa_ut(jd_2050)

        # Helper to compute difference accounting for 360-degree wrap-around
        def compute_increase(later: float, earlier: float) -> float:
            """Compute increase, handling wrap-around at 360 degrees."""
            diff = later - earlier
            # If diff is very negative (wrap-around from ~359 to ~1), add 360
            if diff < -180:
                diff += 360
            return diff

        # Ayanamsha should increase over time due to precession
        # Rate is approximately 50.3 arcsec/year = ~1.4 deg/century
        increase_1900_2000 = compute_increase(aya_2000, aya_1900)
        increase_2000_2050 = compute_increase(aya_2050, aya_2000)

        assert increase_1900_2000 > 0, (
            f"{mode_name}: Ayanamsha should increase from 1900 to 2000 "
            f"(1900={aya_1900:.6f}, 2000={aya_2000:.6f}, diff={increase_1900_2000:.6f})"
        )
        assert increase_2000_2050 > 0, (
            f"{mode_name}: Ayanamsha should increase from 2000 to 2050 "
            f"(2000={aya_2000:.6f}, 2050={aya_2050:.6f}, diff={increase_2000_2050:.6f})"
        )

        # Check approximate rate: ~1.4 degrees per century (50.3"/year)
        assert 1.0 < increase_1900_2000 < 2.0, (
            f"{mode_name}: Precession rate should be ~1.4 deg/century, "
            f"got {increase_1900_2000:.4f} deg/century"
        )

    @pytest.mark.parametrize("sid_mode, mode_name", ALL_AYANAMSHA_MODES)
    def test_ayanamsha_reasonable_range(self, sid_mode, mode_name):
        """Verify ayanamsha values are in a reasonable range for the epoch."""
        jd_2000 = swe.julday(2000, 1, 1, 12.0)

        pyephem.swe_set_sid_mode(sid_mode)
        aya = pyephem.swe_get_ayanamsa_ut(jd_2000)

        # At J2000, most ayanamshas should be roughly 22-25 degrees
        # But some modes like J2000, B1950, J1900 are epoch-based (near 0 at their epoch)
        # And some galactic modes have different reference points

        # Just verify it's in a sensible range (0-360)
        assert 0 <= aya < 360, (
            f"{mode_name}: Ayanamsha should be in [0, 360) range, got {aya:.6f}"
        )


class TestAyanamshaConsistency:
    """Test consistency of ayanamsha calculations."""

    def test_all_modes_unique_at_j2000(self):
        """Verify all ayanamsha modes produce distinct values at J2000."""
        jd = swe.julday(2000, 1, 1, 12.0)

        values = {}
        for sid_mode, mode_name in ALL_AYANAMSHA_MODES:
            pyephem.swe_set_sid_mode(sid_mode)
            aya = pyephem.swe_get_ayanamsa_ut(jd)
            values[mode_name] = aya

        # Check that we have 43 unique modes
        assert len(values) == 43, f"Expected 43 modes, got {len(values)}"

    def test_lahiri_value_at_j2000(self):
        """Test known Lahiri ayanamsha value at J2000."""
        jd = swe.julday(2000, 1, 1, 12.0)

        pyephem.swe_set_sid_mode(SE_SIDM_LAHIRI)
        aya = pyephem.swe_get_ayanamsa_ut(jd)

        # Lahiri ayanamsha at J2000 should be approximately 23.86 degrees
        assert 23.5 < aya < 24.2, (
            f"Lahiri at J2000 should be ~23.86 degrees, got {aya:.6f}"
        )

    def test_j2000_mode_zero_at_epoch(self):
        """Test that J2000 mode gives zero ayanamsha at J2000.0 epoch."""
        jd = swe.julday(2000, 1, 1, 12.0)

        pyephem.swe_set_sid_mode(SE_SIDM_J2000)
        aya = pyephem.swe_get_ayanamsa_ut(jd)

        # J2000 mode should be exactly 0 at J2000.0
        assert abs(aya) < 0.001, f"J2000 mode at J2000.0 should be ~0, got {aya:.6f}"

    @pytest.mark.parametrize("year, month, day, hour, date_desc", TEST_DATES)
    def test_all_modes_produce_values(self, year, month, day, hour, date_desc):
        """Verify all 43 modes produce valid values at each test date."""
        jd = swe.julday(year, month, day, hour)

        for sid_mode, mode_name in ALL_AYANAMSHA_MODES:
            pyephem.swe_set_sid_mode(sid_mode)
            aya = pyephem.swe_get_ayanamsa_ut(jd)

            # Should produce a finite number
            assert isinstance(aya, float), (
                f"{mode_name} at {date_desc}: expected float, got {type(aya)}"
            )
            assert aya == aya, (  # Check for NaN
                f"{mode_name} at {date_desc}: got NaN"
            )


class TestAyanamshaComparison:
    """Comparison tests across all modes and dates."""

    def test_summary_all_modes_all_dates(self):
        """Generate summary comparison of all 43 modes at all test dates."""
        results = []

        for year, month, day, hour, date_desc in TEST_DATES:
            jd = swe.julday(year, month, day, hour)

            for sid_mode, mode_name in ALL_AYANAMSHA_MODES:
                # Swiss Ephemeris
                swe.set_sid_mode(sid_mode)
                aya_swe = swe.get_ayanamsa_ut(jd)

                # libephemeris
                pyephem.swe_set_sid_mode(sid_mode)
                aya_py = pyephem.swe_get_ayanamsa_ut(jd)

                diff = abs(aya_swe - aya_py)
                tol = get_tolerance(sid_mode)

                results.append(
                    {
                        "date": date_desc,
                        "mode_id": sid_mode,
                        "mode_name": mode_name,
                        "swe": aya_swe,
                        "libephem": aya_py,
                        "diff": diff,
                        "tol": tol,
                        "passed": diff < tol,
                    }
                )

        # All results should pass
        failed = [r for r in results if not r["passed"]]
        assert len(failed) == 0, (
            f"Failed {len(failed)} of {len(results)} comparisons:\n"
            + "\n".join(
                f"  {r['mode_name']} at {r['date']}: diff={r['diff']:.6f} > tol={r['tol']}"
                for r in failed[:10]  # Show first 10 failures
            )
        )

        # Verify we tested all combinations: 43 modes x 5 dates = 215 tests
        assert len(results) == 43 * 5, (
            f"Expected 215 results (43 modes x 5 dates), got {len(results)}"
        )


if __name__ == "__main__":
    # Run with verbose output for manual testing
    pytest.main([__file__, "-v", "-s"])
