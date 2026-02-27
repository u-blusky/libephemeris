"""
Tests for EphemerisContext LEB integration.

Validates that EphemerisContext.calc_ut() / calc() correctly use the LEB
fast path when a .leb file is configured, and falls back to Skyfield
when needed.
"""

from __future__ import annotations


import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SE_MEAN_NODE,
    SE_MOON,
    SE_SUN,
    SEFLG_EQUATORIAL,
    SEFLG_SPEED,
    SEFLG_TOPOCTR,
)
from libephemeris.context import EphemerisContext


class TestContextLEBFastPath:
    """Verify that EphemerisContext uses the LEB fast path."""

    @pytest.fixture
    def ctx_with_leb(self, test_leb_file):
        """Create an EphemerisContext with .leb file configured."""
        ctx = EphemerisContext()
        ctx.set_leb_file(test_leb_file)
        yield ctx
        ctx.set_leb_file(None)

    @pytest.fixture
    def jd_mid(self, test_leb_file):
        """Return the midpoint JD of the test LEB file."""
        from libephemeris.leb_reader import LEBReader

        reader = LEBReader(test_leb_file)
        jd_start, jd_end = reader.jd_range
        reader.close()
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_calc_ut_with_leb(self, ctx_with_leb, jd_mid):
        """calc_ut() should produce valid results when .leb is set."""
        result, retflag = ctx_with_leb.calc_ut(jd_mid, SE_SUN, SEFLG_SPEED)

        assert len(result) == 6
        assert 0.0 <= result[0] < 360.0
        assert -90.0 <= result[1] <= 90.0
        # Sun speed should be ~1 deg/day
        assert abs(result[3]) > 0.5

    @pytest.mark.integration
    def test_calc_ut_matches_global(self, ctx_with_leb, jd_mid, test_leb_file):
        """Context calc_ut() with LEB should match global swe_calc_ut() with LEB."""
        # Set up global LEB too
        old_leb = None
        try:
            from libephemeris import state

            old_leb = state._LEB_FILE
            ephem.set_leb_file(test_leb_file)

            ctx_result, _ = ctx_with_leb.calc_ut(jd_mid, SE_SUN, SEFLG_SPEED)
            global_result, _ = ephem.swe_calc_ut(jd_mid, SE_SUN, SEFLG_SPEED)

            # Should match exactly (same LEB reader, same code path)
            for i in range(6):
                assert abs(ctx_result[i] - global_result[i]) < 1e-10, (
                    f"Component {i}: ctx={ctx_result[i]}, global={global_result[i]}"
                )
        finally:
            ephem.set_leb_file(old_leb)

    @pytest.mark.integration
    def test_calc_tt_with_leb(self, ctx_with_leb, jd_mid):
        """calc() (TT) should produce valid results when .leb is set."""
        result, retflag = ctx_with_leb.calc(jd_mid, SE_SUN, SEFLG_SPEED)

        assert len(result) == 6
        assert 0.0 <= result[0] < 360.0

    @pytest.mark.integration
    def test_calc_ut_moon(self, ctx_with_leb, jd_mid):
        """Moon position via context + LEB should be valid."""
        result, _ = ctx_with_leb.calc_ut(jd_mid, SE_MOON, SEFLG_SPEED)

        assert 0.0 <= result[0] < 360.0
        assert -10.0 <= result[1] <= 10.0  # Moon latitude within ~5 deg

    @pytest.mark.integration
    def test_calc_ut_mean_node(self, ctx_with_leb, jd_mid):
        """Mean node via context + LEB (ecliptic direct pipeline)."""
        result, _ = ctx_with_leb.calc_ut(jd_mid, SE_MEAN_NODE, SEFLG_SPEED)

        assert 0.0 <= result[0] < 360.0
        # Mean node retrograde speed: ~-0.053 deg/day
        assert result[3] < 0.0, f"Mean node speed should be negative: {result[3]}"

    @pytest.mark.integration
    def test_calc_ut_equatorial(self, ctx_with_leb, jd_mid):
        """Equatorial coordinates via context + LEB."""
        result, _ = ctx_with_leb.calc_ut(jd_mid, SE_SUN, SEFLG_EQUATORIAL | SEFLG_SPEED)

        assert 0.0 <= result[0] < 360.0
        assert -90.0 <= result[1] <= 90.0


class TestContextLEBFallback:
    """Verify that context falls back correctly when LEB can't serve a request."""

    @pytest.fixture
    def ctx_with_leb(self, test_leb_file):
        ctx = EphemerisContext()
        ctx.set_leb_file(test_leb_file)
        yield ctx
        ctx.set_leb_file(None)

    @pytest.fixture
    def jd_mid(self, test_leb_file):
        from libephemeris.leb_reader import LEBReader

        reader = LEBReader(test_leb_file)
        jd_start, jd_end = reader.jd_range
        reader.close()
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_fallback_for_topoctr(self, ctx_with_leb, jd_mid):
        """SEFLG_TOPOCTR should fall back to Skyfield (not crash)."""
        ctx_with_leb.set_topo(12.5, 41.9, 0)
        result, _ = ctx_with_leb.calc_ut(jd_mid, SE_SUN, SEFLG_TOPOCTR | SEFLG_SPEED)

        # Should get valid results from Skyfield fallback
        assert 0.0 <= result[0] < 360.0

    @pytest.mark.integration
    def test_fallback_for_unknown_body(self, ctx_with_leb, jd_mid):
        """Body not in .leb should fall back to Skyfield."""
        # SE_MARS (4) is in the test LEB, but try a body that isn't
        # Jupiter (5) is not in the test fixture
        result, _ = ctx_with_leb.calc_ut(jd_mid, 5, SEFLG_SPEED)

        # Should get valid results from Skyfield fallback
        assert 0.0 <= result[0] < 360.0


class TestContextLEBGracefulError:
    """Test graceful error handling for invalid .leb files."""

    @pytest.mark.integration
    def test_invalid_leb_path_falls_back(self):
        """Context with invalid .leb path should fall back to Skyfield."""
        ctx = EphemerisContext()
        ctx.set_leb_file("/nonexistent/path/fake.leb")

        # get_leb_reader() should return None (not raise)
        reader = ctx.get_leb_reader()
        assert reader is None

        # calc_ut should still work via Skyfield fallback
        jd = 2460000.0  # ~2023
        result, _ = ctx.calc_ut(jd, SE_SUN, SEFLG_SPEED)
        assert 0.0 <= result[0] < 360.0

    @pytest.mark.integration
    def test_global_invalid_leb_path_falls_back(self):
        """Global set_leb_file with invalid path should fall back."""
        from libephemeris import state

        old_leb = state._LEB_FILE
        try:
            ephem.set_leb_file("/nonexistent/path/fake.leb")
            reader = state.get_leb_reader()
            assert reader is None

            # calc_ut should still work via Skyfield fallback
            jd = 2460000.0
            result, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
            assert 0.0 <= result[0] < 360.0
        finally:
            ephem.set_leb_file(old_leb)


class TestContextLEBGlobalFallthrough:
    """Test that context picks up the global LEB reader when no context-local one."""

    @pytest.mark.integration
    def test_context_uses_global_leb(self, test_leb_file):
        """Context without own .leb should use global reader if set."""
        from libephemeris import state
        from libephemeris.leb_reader import LEBReader

        old_leb = state._LEB_FILE
        try:
            # Set global LEB
            ephem.set_leb_file(test_leb_file)

            # Create context WITHOUT setting .leb
            ctx = EphemerisContext()
            assert ctx._leb_file is None

            reader_tmp = LEBReader(test_leb_file)
            jd_start, jd_end = reader_tmp.jd_range
            jd_mid = (jd_start + jd_end) / 2.0
            reader_tmp.close()

            # calc_ut should use the global LEB reader
            result, _ = ctx.calc_ut(jd_mid, SE_SUN, SEFLG_SPEED)
            assert 0.0 <= result[0] < 360.0
            assert abs(result[3]) > 0.5  # Sun speed confirms it worked
        finally:
            ephem.set_leb_file(old_leb)
