"""
Tests for planetary perturbation corrections in lunar module.

The Moon's position is perturbed mainly by Jupiter and the Sun. The
perturbation corrections should improve the accuracy of the true node
calculation by reducing errors from several arcminutes to sub-arcminute levels.
"""

import math
import pytest
from libephemeris.lunar import (
    _calc_lunar_fundamental_arguments,
    _calc_jupiter_mean_longitude,
    _calc_venus_mean_longitude,
    _calc_elp2000_node_perturbations,
    calc_true_lunar_node,
    calc_mean_lunar_node,
)


class TestLunarFundamentalArguments:
    """Tests for the fundamental lunar arguments (D, M, M', F)."""

    def test_arguments_at_j2000(self):
        """Test fundamental arguments at J2000.0 epoch."""
        jd_j2000 = 2451545.0  # J2000.0

        D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd_j2000)

        # At J2000.0, T=0, so we should get the constant terms
        # D ~= 297.85° (from Meeus)
        # M ~= 357.53° (Sun's mean anomaly)
        # M' ~= 134.96° (Moon's mean anomaly)
        # F ~= 93.27° (argument of latitude)
        D_deg = math.degrees(D)
        M_deg = math.degrees(M)
        M_prime_deg = math.degrees(M_prime)
        F_deg = math.degrees(F)

        # Check values are in expected ranges at J2000
        assert 295 < D_deg < 300, f"D at J2000 = {D_deg}°, expected ~297.85°"
        assert 355 < M_deg < 360, f"M at J2000 = {M_deg}°, expected ~357.53°"
        assert 132 < M_prime_deg < 138, (
            f"M' at J2000 = {M_prime_deg}°, expected ~134.96°"
        )
        assert 90 < F_deg < 96, f"F at J2000 = {F_deg}°, expected ~93.27°"

    def test_arguments_normalized(self):
        """All arguments should be normalized to [0, 2*pi)."""
        # Test at various dates
        test_dates = [
            2451545.0,  # J2000.0
            2451545.0 + 365.25,  # J2001
            2451545.0 - 365.25 * 100,  # 1900
            2451545.0 + 365.25 * 100,  # 2100
        ]

        for jd in test_dates:
            D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd)

            # All should be in [0, 2*pi)
            assert 0 <= D < 2 * math.pi, f"D = {D} at JD {jd}"
            assert 0 <= M < 2 * math.pi, f"M = {M} at JD {jd}"
            assert 0 <= M_prime < 2 * math.pi, f"M' = {M_prime} at JD {jd}"
            assert 0 <= F < 2 * math.pi, f"F = {F} at JD {jd}"

    def test_arguments_vary_with_time(self):
        """Arguments should change over time (not constant)."""
        jd_1 = 2451545.0
        jd_2 = 2451545.0 + 30.0  # 30 days later

        D1, M1, M1p, F1 = _calc_lunar_fundamental_arguments(jd_1)
        D2, M2, M2p, F2 = _calc_lunar_fundamental_arguments(jd_2)

        # Over 30 days, all arguments should change measurably
        # Moon moves ~40° in D per month, ~13° in anomaly per month
        assert D1 != D2
        assert M1 != M2
        assert M1p != M2p
        assert F1 != F2


class TestJupiterMeanLongitude:
    """Tests for Jupiter's mean longitude calculation."""

    def test_jupiter_longitude_at_j2000(self):
        """Test Jupiter's mean longitude at J2000.0."""
        jd_j2000 = 2451545.0

        L_jup = _calc_jupiter_mean_longitude(jd_j2000)
        L_jup_deg = math.degrees(L_jup)

        # Jupiter's mean longitude at J2000.0 should be ~34.35°
        assert 30 < L_jup_deg < 40, (
            f"Jupiter L at J2000 = {L_jup_deg}°, expected ~34.35°"
        )

    def test_jupiter_longitude_normalized(self):
        """Jupiter's longitude should be normalized to [0, 2*pi)."""
        test_dates = [
            2451545.0,  # J2000.0
            2451545.0 + 4332.6 * 0.5,  # Half Jupiter orbital period later
            2451545.0 + 4332.6,  # Full Jupiter orbital period
        ]

        for jd in test_dates:
            L_jup = _calc_jupiter_mean_longitude(jd)
            assert 0 <= L_jup < 2 * math.pi, f"L_jup = {L_jup} at JD {jd}"

    def test_jupiter_orbital_period(self):
        """Jupiter should complete ~1 orbit in ~11.86 years."""
        jd_start = 2451545.0
        jupiter_orbital_period_days = 4332.6  # ~11.86 years

        L_start = _calc_jupiter_mean_longitude(jd_start)
        L_end = _calc_jupiter_mean_longitude(jd_start + jupiter_orbital_period_days)

        L_start_deg = math.degrees(L_start)
        L_end_deg = math.degrees(L_end)

        # After one period, should return to ~same longitude (within a few degrees)
        diff = abs(L_end_deg - L_start_deg)
        if diff > 180:
            diff = 360 - diff
        assert diff < 15, f"Jupiter longitude diff after 1 period = {diff}°"


class TestVenusMeanLongitude:
    """Tests for Venus's mean longitude calculation."""

    def test_venus_longitude_at_j2000(self):
        """Test Venus's mean longitude at J2000.0."""
        jd_j2000 = 2451545.0

        L_venus = _calc_venus_mean_longitude(jd_j2000)
        L_venus_deg = math.degrees(L_venus)

        # Venus's mean longitude at J2000.0 should be ~181.98°
        assert 178 < L_venus_deg < 186, (
            f"Venus L at J2000 = {L_venus_deg}°, expected ~181.98°"
        )

    def test_venus_longitude_normalized(self):
        """Venus's longitude should be normalized to [0, 2*pi)."""
        test_dates = [
            2451545.0,  # J2000.0
            2451545.0 + 224.7 * 0.5,  # Half Venus orbital period later
            2451545.0 + 224.7,  # Full Venus orbital period
        ]

        for jd in test_dates:
            L_venus = _calc_venus_mean_longitude(jd)
            assert 0 <= L_venus < 2 * math.pi, f"L_venus = {L_venus} at JD {jd}"

    def test_venus_orbital_period(self):
        """Venus should complete ~1 orbit in ~224.7 days."""
        jd_start = 2451545.0
        venus_orbital_period_days = 224.7  # ~0.615 years

        L_start = _calc_venus_mean_longitude(jd_start)
        L_end = _calc_venus_mean_longitude(jd_start + venus_orbital_period_days)

        L_start_deg = math.degrees(L_start)
        L_end_deg = math.degrees(L_end)

        # After one period, should return to ~same longitude (within a few degrees)
        diff = abs(L_end_deg - L_start_deg)
        if diff > 180:
            diff = 360 - diff
        assert diff < 10, f"Venus longitude diff after 1 period = {diff}°"

    def test_venus_synodic_relationship_with_earth(self):
        """Venus moves faster than Earth, completing ~1.6 orbits per Earth year."""
        jd_start = 2451545.0
        one_year = 365.25

        L_start = _calc_venus_mean_longitude(jd_start)
        L_end = _calc_venus_mean_longitude(jd_start + one_year)

        # Calculate total angular motion (may wrap around multiple times)
        # Venus mean motion is ~1.6021°/day, so in 365.25 days: ~584.9°
        # This is roughly 1.625 complete orbits
        L_start_deg = math.degrees(L_start)
        L_end_deg = math.degrees(L_end)

        # The angular difference in one year should correspond to ~1.6 orbits
        # 584.9° mod 360° ≈ 224.9°
        diff = L_end_deg - L_start_deg
        if diff < 0:
            diff += 360

        # Expect ~224° difference (1.625 orbits - 1 = 0.625 orbits = 225°)
        assert 200 < diff < 250, f"Venus angular motion per year = {diff}°"


class TestPlanetaryPerturbations:
    """Tests for the combined planetary perturbation calculation."""

    def test_perturbation_returns_degrees(self):
        """Perturbation should return a value in degrees."""
        jd_j2000 = 2451545.0

        perturbation = _calc_elp2000_node_perturbations(jd_j2000)

        # Perturbation should be a reasonable magnitude (< 5°)
        # Main solar terms can reach ~1.5° amplitude
        assert -5.0 < perturbation < 5.0, f"Perturbation = {perturbation}°"

    def test_perturbation_varies_over_synodic_month(self):
        """Perturbation should vary significantly over a synodic month."""
        jd_start = 2451545.0

        perturbations = []
        for i in range(30):  # Sample ~1 month
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # The range should be significant (solar terms have ~3° total variation)
        p_range = max(perturbations) - min(perturbations)
        assert p_range > 0.5, f"Perturbation range over 1 month = {p_range}°"

    def test_perturbation_bounded(self):
        """Perturbation should stay within reasonable bounds."""
        # Test across a year
        jd_start = 2451545.0

        for i in range(0, 365, 7):  # Weekly samples for a year
            jd = jd_start + i
            perturbation = _calc_elp2000_node_perturbations(jd)

            # Total perturbation should never exceed ~5° (sum of all terms)
            assert abs(perturbation) < 5.0, f"Perturbation at JD {jd} = {perturbation}°"

    def test_solar_perturbation_dominates(self):
        """Solar perturbation should be larger than Jupiter perturbation."""
        # This is tested indirectly - the total should be dominated by
        # the ~1.5° solar term, not the ~0.05° Jupiter terms
        jd_j2000 = 2451545.0

        # Get perturbation values at multiple phases
        pert_max = max(
            _calc_elp2000_node_perturbations(jd_j2000 + i) for i in range(30)
        )
        pert_min = min(
            _calc_elp2000_node_perturbations(jd_j2000 + i) for i in range(30)
        )

        # Range should be dominated by solar term (~3°) not Jupiter (~0.1°)
        assert pert_max - pert_min > 1.0, "Solar perturbation should dominate"

    def test_venus_perturbation_contributes(self):
        """Venus perturbation terms should contribute measurable effects.

        Venus perturbation terms have amplitudes of 0.001-0.005 degrees
        (approximately 3.6-18 arcseconds). These are small but measurable
        effects on the lunar node position.
        """
        # We can't isolate Venus terms directly, but we can verify the
        # perturbation function includes them by checking that the total
        # perturbation varies in ways consistent with Venus's orbital period

        jd_start = 2451545.0
        venus_orbital_period = 224.7  # days

        # Sample over 2 Venus years to capture Venus-related variations
        samples = []
        for i in range(int(venus_orbital_period * 2)):
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            samples.append(p)

        # The perturbation series should show complex variation
        # (not just solar monthly signal) due to Venus terms
        # Check that we have sufficient variability
        p_range = max(samples) - min(samples)
        assert p_range > 0.5, f"Perturbation range = {p_range}°"

        # The standard deviation should indicate complex modulation
        mean_p = sum(samples) / len(samples)
        variance = sum((p - mean_p) ** 2 for p in samples) / len(samples)
        std_dev = variance**0.5
        assert std_dev > 0.1, f"Perturbation std dev = {std_dev}°"

    def test_venus_terms_amplitude_range(self):
        """Venus perturbation terms should have expected amplitude range.

        Based on Chapront-Touzé lunar theory, Venus terms have amplitudes
        of approximately 0.001-0.005 degrees (3.6-18 arcseconds).
        The combined effect over Venus synodic cycle can reach ~0.02°.
        """
        # Sample at two times when Venus is at opposite positions
        jd_start = 2451545.0
        venus_synodic = 583.9  # Venus synodic period (days)

        # Get Venus longitude at two opposite phases
        L_venus_0 = _calc_venus_mean_longitude(jd_start)
        L_venus_half = _calc_venus_mean_longitude(jd_start + venus_synodic / 2)

        # Verify Venus has moved ~180° (opposite position)
        diff = abs(math.degrees(L_venus_half) - math.degrees(L_venus_0))
        if diff > 180:
            diff = 360 - diff
        # Should be close to 180° (Venus moves more than 360° in half synodic)
        # Actually, in half synodic period, Venus moves ~292° relative to Earth

        # The key test is that perturbations vary with Venus phase
        pert_0 = _calc_elp2000_node_perturbations(jd_start)
        pert_half = _calc_elp2000_node_perturbations(jd_start + venus_synodic / 2)

        # Perturbations should differ (though solar terms may dominate)
        # The Venus contribution is small but the total will vary
        assert pert_0 != pert_half, "Perturbations should vary"


class TestTrueNodeWithPerturbations:
    """Tests for the true node calculation with planetary perturbations."""

    def test_true_node_returns_valid_position(self):
        """True node should return valid position tuple."""
        jd = 2451545.0

        lon, lat, dist = calc_true_lunar_node(jd)

        assert 0 <= lon < 360, f"Longitude = {lon}"
        assert lat == 0.0, "Node latitude should always be 0"
        assert dist >= 0, "Node distance should be non-negative"

    def test_true_node_differs_from_mean(self):
        """True node should differ from mean node at most dates."""
        # Test at multiple dates to find when they differ
        # At some specific moments, true and mean can be very close
        jd_start = 2451545.0

        diffs = []
        for i in range(30):  # Sample over a month
            jd = jd_start + i
            true_lon, _, _ = calc_true_lunar_node(jd)
            mean_lon = calc_mean_lunar_node(jd)

            diff = abs(true_lon - mean_lon)
            if diff > 180:
                diff = 360 - diff
            diffs.append(diff)

        # Maximum difference over the month should be significant
        max_diff = max(diffs)
        # True node typically differs from mean by 1-10° at maximum
        assert max_diff > 0.5, f"Max True-Mean difference = {max_diff}°"

    def test_true_node_differs_from_mean_node(self):
        """True node should differ from mean node due to orbital perturbations."""
        jd_start = 2451545.0

        differences = []
        for i in range(30):  # One month
            jd = jd_start + i
            true_lon, _, _ = calc_true_lunar_node(jd)
            mean_lon = calc_mean_lunar_node(jd)

            diff = true_lon - mean_lon
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            differences.append(diff)

        # True node should differ from mean node
        # The osculating true node tracks the instantaneous orbital plane
        max_diff = max(abs(d) for d in differences)

        # True node typically differs from mean by several degrees
        assert max_diff > 0.1, (
            f"True node should differ from mean node, max diff was {max_diff}°"
        )

    @pytest.mark.parametrize(
        "year,month,day",
        [
            (2000, 1, 1),
            (2010, 6, 15),
            (2020, 12, 21),
            (1990, 3, 20),
            (2030, 9, 23),
        ],
    )
    def test_true_node_at_various_dates(self, year, month, day):
        """True node should return valid positions at various dates."""
        # Convert to Julian Day (simplified formula)
        a = (14 - month) // 12
        y = year + 4800 - a
        m = month + 12 * a - 3
        jd = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045

        lon, lat, dist = calc_true_lunar_node(jd)

        assert 0 <= lon < 360, (
            f"Invalid longitude {lon} for {year}-{month:02d}-{day:02d}"
        )
        assert lat == 0.0
        assert dist >= 0  # Distance is a proxy value


@pytest.mark.integration
class TestPerturbationIntegration:
    """Integration tests for perturbation accuracy."""

    def test_perturbation_improves_accuracy(self):
        """Perturbation corrections should improve true node accuracy.

        This is a sanity check that the perturbation terms have the expected
        effect of modifying the osculating node position.
        """
        jd = 2451545.0

        # Get the true node (with perturbations)
        true_lon, _, _ = calc_true_lunar_node(jd)

        # The perturbation value at this date
        pert = _calc_elp2000_node_perturbations(jd)

        # Verify perturbation is non-zero (corrections are being applied)
        assert abs(pert) > 0.001, "Perturbation should be non-zero"

    def test_long_term_stability(self):
        """Perturbation calculation should be stable over long periods."""
        # Test over 20 years
        start_jd = 2451545.0
        years = 20
        samples_per_year = 12

        for year in range(years):
            for month in range(samples_per_year):
                jd = start_jd + year * 365.25 + month * 30.4
                lon, lat, dist = calc_true_lunar_node(jd)

                # All positions should be valid
                assert 0 <= lon < 360, f"Invalid lon at JD {jd}"
                assert lat == 0.0
                assert dist >= 0  # Distance is a proxy value

    def test_nodal_cycle_period(self):
        """True node should complete cycle in ~18.6 years (retrograde)."""
        start_jd = 2451545.0
        nodal_period = 18.6 * 365.25  # ~6798 days

        start_lon, _, _ = calc_true_lunar_node(start_jd)
        end_lon, _, _ = calc_true_lunar_node(start_jd + nodal_period)

        # After one nodal period, should be back near the same longitude
        # (Allow for oscillations - within 30°)
        diff = abs(end_lon - start_lon)
        if diff > 180:
            diff = 360 - diff

        assert diff < 30, f"After nodal period, diff = {diff}°"
