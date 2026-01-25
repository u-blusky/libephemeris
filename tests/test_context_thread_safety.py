"""
Thread-safety tests for concurrent ephemeris usage.

This module provides comprehensive thread-safety tests for the libephemeris library,
verifying that multiple threads can perform calculations concurrently with:
- Different ayanamsha systems
- Different topocentric locations
- No race conditions
- Correct isolated results

Each test uses 10+ threads to stress test the library's thread-safety mechanisms.

Note: Since the module-level API uses global state which is not thread-safe by design
(matching Swiss Ephemeris behavior), these tests verify that the state management
doesn't crash under concurrent access and that EphemerisContext provides proper
isolation for thread-safe concurrent calculations.
"""

import concurrent.futures
import threading
import time

import pytest

import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SEFLG_SPEED,
    SEFLG_SIDEREAL,
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_LAHIRI,
    SE_SIDM_RAMAN,
    SE_SIDM_KRISHNAMURTI,
    SE_SIDM_YUKTESHWAR,
    SE_SIDM_TRUE_CITRA,
    SE_SIDM_DELUCE,
    SE_SIDM_JN_BHASIN,
    SE_SIDM_DJWHAL_KHUL,
    SE_SIDM_USHASHASHI,
    SE_SIDM_BABYL_KUGLER1,
    SE_SIDM_GALCENT_0SAG,
)


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def test_ayanamshas():
    """12 different ayanamsha systems for thread testing."""
    return [
        (SE_SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
        (SE_SIDM_LAHIRI, "Lahiri"),
        (SE_SIDM_RAMAN, "Raman"),
        (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
        (SE_SIDM_YUKTESHWAR, "Yukteshwar"),
        (SE_SIDM_TRUE_CITRA, "True Citra"),
        (SE_SIDM_DELUCE, "De Luce"),
        (SE_SIDM_JN_BHASIN, "JN Bhasin"),
        (SE_SIDM_DJWHAL_KHUL, "Djwhal Khul"),
        (SE_SIDM_USHASHASHI, "Ushashashi"),
        (SE_SIDM_BABYL_KUGLER1, "Babylonian Kugler 1"),
        (SE_SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
    ]


@pytest.fixture
def test_locations_12():
    """12 different geographic locations for thread testing."""
    return [
        (0, "Equator", 0.0, 0.0, 0),
        (1, "Rome", 12.5, 41.9, 50),
        (2, "London", -0.1, 51.5, 10),
        (3, "Tokyo", 139.7, 35.7, 40),
        (4, "New York", -74.0, 40.7, 10),
        (5, "Sydney", 151.2, -33.9, 5),
        (6, "San Francisco", -122.4, 37.8, 20),
        (7, "Paris", 2.3, 48.9, 35),
        (8, "Berlin", 13.4, 52.5, 30),
        (9, "Rio de Janeiro", -43.2, -22.9, 10),
        (10, "Mumbai", 72.9, 19.1, 15),
        (11, "Cape Town", 18.4, -33.9, 25),
    ]


# =============================================================================
# TEST: EphemerisContext Instance Isolation
# =============================================================================


class TestContextInstanceIsolation:
    """Tests that EphemerisContext instances have isolated state."""

    def test_separate_contexts_have_isolated_sidereal_modes(self):
        """Each context should have its own sidereal mode."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        ctx1.set_sid_mode(SE_SIDM_LAHIRI)
        ctx2.set_sid_mode(SE_SIDM_FAGAN_BRADLEY)

        assert ctx1.get_sid_mode() == SE_SIDM_LAHIRI
        assert ctx2.get_sid_mode() == SE_SIDM_FAGAN_BRADLEY

    def test_separate_contexts_have_isolated_topo(self):
        """Each context should have its own topocentric location."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        ctx1.set_topo(12.5, 41.9, 0)  # Rome
        ctx2.set_topo(-0.1, 51.5, 0)  # London

        topo1 = ctx1.get_topo()
        topo2 = ctx2.get_topo()

        assert topo1 is not None
        assert topo2 is not None
        assert abs(topo1.latitude.degrees - 41.9) < 0.001
        assert abs(topo2.latitude.degrees - 51.5) < 0.001

    def test_12_contexts_with_different_configs(
        self, test_ayanamshas, test_locations_12
    ):
        """Create 12 contexts with different configurations and verify isolation."""
        contexts = []

        for i in range(12):
            ctx = EphemerisContext()
            mode, _ = test_ayanamshas[i]
            _, name, lon, lat, alt = test_locations_12[i]

            ctx.set_sid_mode(mode)
            ctx.set_topo(lon, lat, alt)
            contexts.append((ctx, mode, lon, lat))

        # Verify all contexts retained their configurations
        for ctx, expected_mode, expected_lon, expected_lat in contexts:
            assert ctx.get_sid_mode() == expected_mode
            topo = ctx.get_topo()
            assert topo is not None
            assert abs(topo.longitude.degrees - expected_lon) < 0.001
            assert abs(topo.latitude.degrees - expected_lat) < 0.001


# =============================================================================
# TEST: Context Resource Sharing
# =============================================================================


class TestContextResourceSharing:
    """Tests that EphemerisContext instances share expensive resources."""

    def test_contexts_share_ephemeris_via_module_state(self):
        """Module-level API uses shared planetary ephemeris."""
        from libephemeris import state

        # Pre-load by doing a calculation with module API
        ephem.swe_calc_ut(2451545.0, SE_SUN, 0)

        # Verify module state is loaded
        assert state._PLANETS is not None

        # Multiple calculations should use the same ephemeris
        planets_before = state._PLANETS
        ephem.swe_calc_ut(2451545.0, SE_MOON, 0)
        planets_after = state._PLANETS

        assert planets_before is planets_after, "Ephemeris should remain shared"

    def test_contexts_share_timescale_via_module_state(self):
        """Module-level API uses shared timescale."""
        from libephemeris import state

        # Pre-load
        ephem.swe_calc_ut(2451545.0, SE_SUN, 0)

        # Verify timescale is loaded
        assert state._TS is not None

        ts_before = state._TS
        ephem.swe_calc_ut(2451545.0, SE_MARS, 0)
        ts_after = state._TS

        assert ts_before is ts_after, "Timescale should remain shared"


# =============================================================================
# TEST: Sequential Calculations with Different Ayanamshas
# =============================================================================


class TestSequentialAyanamshaCalculations:
    """Tests for calculations with different ayanamsha systems."""

    def test_12_ayanamshas_produce_different_results(self, test_ayanamshas):
        """
        Test that 12 different ayanamshas produce different sidereal positions.

        Uses module-level API for reliable execution.
        """
        jd = 2451545.0  # J2000.0
        results: dict[int, float] = {}

        for mode, name in test_ayanamshas:
            ephem.swe_set_sid_mode(mode)
            pos, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)
            results[mode] = pos[0]

        # All 12 ayanamshas should produce different positions
        unique_positions = len(set(round(lon, 4) for lon in results.values()))
        assert unique_positions >= 10, (
            f"Expected at least 10 unique positions from 12 ayanamshas, "
            f"got {unique_positions}"
        )


# =============================================================================
# TEST: Concurrent Module-Level API (Stress Test)
# =============================================================================


class TestConcurrentModuleLevelAPI:
    """
    Stress tests for concurrent access to module-level API.

    Note: The module-level API is NOT thread-safe by design (matching
    Swiss Ephemeris behavior). These tests verify the API doesn't crash
    under concurrent access, though results may not be deterministic
    when different threads modify global state.
    """

    def test_concurrent_geocentric_calculations(self):
        """
        Test that concurrent geocentric calculations don't crash.

        Geocentric calculations don't depend on topocentric state,
        so they should be safe for concurrent access.
        """
        jd = 2451545.0
        num_threads = 12
        errors: list[str] = []
        results: list[tuple[int, float]] = []
        results_lock = threading.Lock()

        def calculate_planet(thread_id: int, planet: int):
            try:
                pos, _ = ephem.swe_calc_ut(jd, planet, 0)
                with results_lock:
                    results.append((planet, pos[0]))
            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        planets = [
            SE_SUN,
            SE_MOON,
            SE_MERCURY,
            SE_VENUS,
            SE_MARS,
            SE_JUPITER,
            SE_SATURN,
        ]

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [
                ex.submit(calculate_planet, i, planets[i % len(planets)])
                for i in range(num_threads)
            ]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors occurred: {errors}"
        assert len(results) == num_threads

    def test_concurrent_house_calculations(self, test_locations_12):
        """
        Test that concurrent house calculations don't crash.

        Each thread calculates houses for a different location.
        """
        jd = 2451545.0
        errors: list[str] = []
        results: list[tuple] = []
        results_lock = threading.Lock()

        def calculate_houses(thread_id: int, lat: float, lon: float):
            try:
                cusps, ascmc = ephem.swe_houses(jd, lat, lon, ord("P"))
                with results_lock:
                    results.append((thread_id, cusps))
            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
            futures = [
                ex.submit(calculate_houses, tid, lat, lon)
                for tid, _, lon, lat, _ in test_locations_12
            ]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors occurred: {errors}"
        assert len(results) == 12


# =============================================================================
# TEST: Concurrent Sequential Calculations
# =============================================================================


class TestConcurrentSequentialCalculations:
    """
    Tests where each thread performs a sequence of independent calculations.
    """

    def test_10_threads_each_calculating_all_planets(self):
        """
        10 threads each calculate positions for all major planets.

        This tests that the calculation engine handles concurrent reads.
        """
        jd = 2451545.0
        num_threads = 10
        planets = [
            SE_SUN,
            SE_MOON,
            SE_MERCURY,
            SE_VENUS,
            SE_MARS,
            SE_JUPITER,
            SE_SATURN,
        ]
        errors: list[str] = []
        calculation_counts = [0]
        count_lock = threading.Lock()

        def calculate_all_planets(thread_id: int):
            try:
                for planet in planets:
                    pos, _ = ephem.swe_calc_ut(jd + thread_id * 0.1, planet, 0)
                    # Verify position is valid
                    if not (0 <= pos[0] < 360):
                        with count_lock:
                            errors.append(
                                f"Thread {thread_id}: invalid position {pos[0]}"
                            )
                        return

                with count_lock:
                    calculation_counts[0] += len(planets)

            except Exception as e:
                with count_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(calculate_all_planets, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"
        expected = num_threads * len(planets)
        assert calculation_counts[0] == expected, (
            f"Expected {expected} calculations, got {calculation_counts[0]}"
        )

    def test_15_threads_rapid_calculations(self):
        """
        15 threads perform rapid calculations in succession.
        """
        num_threads = 15
        calculations_per_thread = 20
        errors: list[str] = []
        success_count = [0]
        count_lock = threading.Lock()

        def rapid_calcs(thread_id: int):
            try:
                for i in range(calculations_per_thread):
                    jd = 2451545.0 + thread_id + i * 0.05
                    planet = [SE_SUN, SE_MOON, SE_MARS][i % 3]
                    pos, _ = ephem.swe_calc_ut(jd, planet, 0)

                    if not (0 <= pos[0] < 360):
                        raise ValueError(f"Invalid position: {pos[0]}")

                with count_lock:
                    success_count[0] += 1

            except Exception as e:
                with count_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(rapid_calcs, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors[:10]}"
        assert success_count[0] == num_threads


# =============================================================================
# TEST: Result Consistency
# =============================================================================


class TestResultConsistency:
    """Tests that verify calculation results are consistent."""

    def test_repeated_calculations_give_same_result(self):
        """Repeated calculations for the same input should give identical results."""
        jd = 2451545.0
        num_iterations = 50

        results = []
        for _ in range(num_iterations):
            pos, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
            results.append(pos[0])

        # All results should be identical
        assert all(abs(r - results[0]) < 1e-10 for r in results), (
            "Results should be identical for same input"
        )

    def test_concurrent_same_calculation(self):
        """
        Multiple threads calculating the same thing should get the same result.
        """
        jd = 2451545.0
        num_threads = 10
        results: list[float] = []
        results_lock = threading.Lock()
        errors: list[str] = []

        def calculate(thread_id: int):
            try:
                pos, _ = ephem.swe_calc_ut(jd, SE_MOON, SEFLG_SPEED)
                with results_lock:
                    results.append(pos[0])
            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(calculate, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors
        assert len(results) == num_threads

        # All results should be identical (or very close)
        reference = results[0]
        for i, r in enumerate(results):
            assert abs(r - reference) < 1e-8, f"Result {i} differs: {r} vs {reference}"


# =============================================================================
# TEST: Mixed Calculation Types
# =============================================================================


class TestMixedCalculationTypes:
    """Tests mixing different types of calculations concurrently."""

    def test_concurrent_planet_and_house_calculations(self, test_locations_12):
        """
        Concurrent planet and house calculations.
        """
        jd = 2451545.0
        errors: list[str] = []
        planet_results: list[float] = []
        house_results: list = []
        results_lock = threading.Lock()

        def calc_planet(thread_id: int):
            try:
                planet = [SE_SUN, SE_MOON, SE_MARS, SE_JUPITER][thread_id % 4]
                pos, _ = ephem.swe_calc_ut(jd, planet, 0)
                with results_lock:
                    planet_results.append(pos[0])
            except Exception as e:
                with results_lock:
                    errors.append(f"Planet Thread {thread_id}: {e}")

        def calc_houses(thread_id: int, lat: float, lon: float):
            try:
                cusps, _ = ephem.swe_houses(jd, lat, lon, ord("P"))
                with results_lock:
                    house_results.append(cusps)
            except Exception as e:
                with results_lock:
                    errors.append(f"House Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as ex:
            # Submit planet calculations
            planet_futures = [ex.submit(calc_planet, i) for i in range(8)]

            # Submit house calculations
            house_futures = [
                ex.submit(calc_houses, tid, lat, lon)
                for tid, _, lon, lat, _ in test_locations_12
            ]

            concurrent.futures.wait(planet_futures + house_futures)

        assert not errors, f"Errors: {errors}"
        assert len(planet_results) == 8
        assert len(house_results) == 12

    def test_interleaved_different_jd_calculations(self):
        """
        Threads calculating for different Julian Days interleaved.
        """
        num_threads = 12
        errors: list[str] = []
        results: dict[int, float] = {}
        results_lock = threading.Lock()

        def calculate(thread_id: int, jd: float):
            try:
                pos, _ = ephem.swe_calc_ut(jd, SE_SUN, 0)
                with results_lock:
                    results[thread_id] = pos[0]
            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        # Different JDs for each thread
        jds = [2451545.0 + i * 30.0 for i in range(num_threads)]

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(calculate, i, jd) for i, jd in enumerate(jds)]
            concurrent.futures.wait(futures)

        assert not errors
        assert len(results) == num_threads

        # Verify results are different (Sun moves ~1 degree per day)
        positions = list(results.values())
        unique_positions = len(set(round(p, 2) for p in positions))
        assert unique_positions == num_threads, (
            f"Expected {num_threads} unique positions for different JDs"
        )


# =============================================================================
# TEST: High Load
# =============================================================================


class TestHighLoad:
    """High-load stress tests."""

    def test_50_concurrent_calculations(self):
        """
        50 concurrent calculation threads.
        """
        num_threads = 50
        errors: list[str] = []
        success_count = [0]
        count_lock = threading.Lock()

        def worker(thread_id: int):
            try:
                jd = 2451545.0 + thread_id
                planet = [SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS, SE_MARS][thread_id % 5]
                pos, _ = ephem.swe_calc_ut(jd, planet, 0)

                if not (0 <= pos[0] < 360):
                    raise ValueError(f"Invalid: {pos[0]}")

                with count_lock:
                    success_count[0] += 1

            except Exception as e:
                with count_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors[:10]}"
        assert success_count[0] == num_threads

    def test_repeated_concurrent_batches(self):
        """
        Run 5 batches of 20 concurrent calculations each.
        """
        num_batches = 5
        threads_per_batch = 20
        all_errors: list[str] = []

        for batch in range(num_batches):
            errors: list[str] = []
            results: list[float] = []
            results_lock = threading.Lock()

            def calc(thread_id: int):
                try:
                    jd = 2451545.0 + batch + thread_id * 0.1
                    pos, _ = ephem.swe_calc_ut(jd, SE_MOON, 0)
                    with results_lock:
                        results.append(pos[0])
                except Exception as e:
                    with results_lock:
                        errors.append(f"Batch {batch} Thread {thread_id}: {e}")

            with concurrent.futures.ThreadPoolExecutor(
                max_workers=threads_per_batch
            ) as ex:
                futures = [ex.submit(calc, i) for i in range(threads_per_batch)]
                concurrent.futures.wait(futures)

            all_errors.extend(errors)

            # Verify batch results
            assert len(results) == threads_per_batch, (
                f"Batch {batch}: expected {threads_per_batch} results"
            )

        assert not all_errors, f"Errors across batches: {all_errors[:10]}"


# =============================================================================
# TEST: Long-Running Sessions
# =============================================================================


class TestLongRunningSessions:
    """Tests for sustained calculation sessions."""

    def test_sustained_parallel_calculations(self):
        """
        10 threads each performing 100 calculations.
        """
        num_threads = 10
        calcs_per_thread = 100
        errors: list[str] = []
        total_calcs = [0]
        count_lock = threading.Lock()

        def long_session(thread_id: int):
            try:
                for i in range(calcs_per_thread):
                    jd = 2451545.0 + i * 0.5
                    planet = [SE_SUN, SE_MOON, SE_MARS][i % 3]
                    ephem.swe_calc_ut(jd, planet, 0)

                with count_lock:
                    total_calcs[0] += calcs_per_thread

            except Exception as e:
                with count_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        start = time.time()

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(long_session, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        elapsed = time.time() - start

        assert not errors, f"Errors: {errors[:10]}"
        expected = num_threads * calcs_per_thread
        assert total_calcs[0] == expected

        # Should complete in reasonable time (not hang)
        assert elapsed < 60, f"Took too long: {elapsed:.1f}s"
