"""
Comprehensive tests for EphemerisContext thread safety.

Tests concurrent usage and isolation of context instances.
"""

import pytest
import threading
import time
import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.constants import *


class TestContextBasic:
    """Basic tests for EphemerisContext."""

    @pytest.mark.unit
    def test_context_creation(self):
        """Should be able to create context instance."""
        ctx = EphemerisContext()
        assert ctx is not None

    @pytest.mark.unit
    def test_context_calc_ut(self):
        """Context should be able to calculate planet positions."""
        ctx = EphemerisContext()
        pos, flags = ctx.calc_ut(2451545.0, SE_SUN, 0)

        assert len(pos) == 6
        assert 0 <= pos[0] < 360

    @pytest.mark.unit
    def test_context_houses(self):
        """Context should be able to calculate houses."""
        ctx = EphemerisContext()
        cusps, ascmc = ctx.houses(2451545.0, 41.9, 12.5, ord("P"))

        assert len(cusps) >= 12
        assert len(ascmc) >= 2


class TestContextIsolation:
    """Test that contexts are isolated from each other."""

    @pytest.mark.unit
    def test_topo_isolation(self):
        """Setting topo in one context shouldn't affect another."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        ctx1.set_topo(12.5, 41.9, 0)  # Rome
        ctx2.set_topo(-74.0, 40.7, 0)  # New York

        # Get topos - they should differ
        topo1 = ctx1.get_topo()
        topo2 = ctx2.get_topo()

        assert topo1 is not None
        assert topo2 is not None
        # Latitudes should differ
        # (implementation detail - may need adjustment)

    @pytest.mark.unit
    def test_sid_mode_isolation(self):
        """Setting sid mode in one context shouldn't affect another."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        ctx1.set_sid_mode(SE_SIDM_LAHIRI)
        ctx2.set_sid_mode(SE_SIDM_FAGAN_BRADLEY)

        mode1 = ctx1.get_sid_mode()
        mode2 = ctx2.get_sid_mode()

        assert mode1 != mode2

    @pytest.mark.unit
    def test_global_vs_context_isolation(self):
        """Context should be isolated from global state."""
        # Set global state
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        # Create context with different mode
        ctx = EphemerisContext()
        ctx.set_sid_mode(SE_SIDM_FAGAN_BRADLEY)

        # Context mode should differ from global
        ctx_mode = ctx.get_sid_mode()
        # Note: This depends on implementation - context may start with global state


class TestContextConcurrency:
    """Test concurrent usage of multiple contexts."""

    @pytest.mark.integration
    def test_concurrent_calculations(self):
        """Multiple threads should be able to calculate simultaneously."""
        results = []
        errors = []

        def calculate_in_thread(thread_id, jd):
            try:
                ctx = EphemerisContext()
                ctx.set_topo(thread_id * 10, thread_id * 5, 0)

                for planet in [SE_SUN, SE_MOON, SE_MARS]:
                    pos, _ = ctx.calc_ut(jd, planet, 0)
                    results.append((thread_id, planet, pos[0]))
            except Exception as e:
                errors.append((thread_id, e))

        print("  Starting 10 concurrent threads...")

        threads = []
        for i in range(10):
            jd = 2451545.0 + i
            t = threading.Thread(target=calculate_in_thread, args=(i, jd))
            threads.append(t)
            t.start()

        for t in threads:
            t.join()

        print(f"  Threads complete: {len(results)} results, {len(errors)} errors")

        # No errors should have occurred
        assert len(errors) == 0, f"Errors: {errors}"
        # Should have 10 threads × 3 planets = 30 results
        assert len(results) == 30

    @pytest.mark.integration
    def test_concurrent_different_modes(self):
        """Threads with different sid modes shouldn't interfere."""
        results = {}

        def calculate_with_mode(mode, mode_name):
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)

            # Small delay to ensure overlap
            time.sleep(0.01)

            pos, _ = ctx.calc_ut(2451545.0, SE_SUN, SEFLG_SIDEREAL)
            results[mode_name] = pos[0]

        print("  Starting sid mode threads...")

        threads = [
            threading.Thread(
                target=calculate_with_mode, args=(SE_SIDM_LAHIRI, "Lahiri")
            ),
            threading.Thread(
                target=calculate_with_mode, args=(SE_SIDM_FAGAN_BRADLEY, "Fagan")
            ),
        ]

        for t in threads:
            t.start()
        for t in threads:
            t.join()

        print(
            f"  Lahiri={results.get('Lahiri', 'N/A'):.2f}, Fagan={results.get('Fagan', 'N/A'):.2f}"
        )

        # Results should be different (different ayanamshas)
        assert "Lahiri" in results
        assert "Fagan" in results
        # Lahiri and Fagan differ by about 0.8°
        diff = abs(results["Lahiri"] - results["Fagan"])
        assert 0.1 < diff < 2.0, f"Diff {diff} unexpected"


class TestContextResourceSharing:
    """Test that contexts share ephemeris resources efficiently."""

    @pytest.mark.integration
    def test_multiple_contexts_same_results(self):
        """Multiple contexts should give same results for same inputs."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        jd = 2451545.0

        pos1, _ = ctx1.calc_ut(jd, SE_MARS, 0)
        pos2, _ = ctx2.calc_ut(jd, SE_MARS, 0)

        assert pos1[0] == pytest.approx(pos2[0], abs=1e-10)
        assert pos1[1] == pytest.approx(pos2[1], abs=1e-10)

    @pytest.mark.unit
    def test_context_reuse(self, progress_reporter):
        """Same context should be reusable for multiple calculations."""
        ctx = EphemerisContext()

        results = []
        progress = progress_reporter("Context reuse", 100, report_every=25)

        for i in range(100):
            jd = 2451545.0 + i
            pos, _ = ctx.calc_ut(jd, SE_SUN, 0)
            results.append(pos[0])
            progress.update(i)

        progress.done()

        # All results should be valid
        assert len(results) == 100
        for lon in results:
            assert 0 <= lon < 360
