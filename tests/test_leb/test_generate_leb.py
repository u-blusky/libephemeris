"""
Tests for the LEB generator script (generate_leb.py).

Tests that the generator produces valid .leb files and that round-trip
(generate → read → evaluate) produces correct results.
"""

from __future__ import annotations

import os
import sys

import pytest

# Ensure scripts directory is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))

from libephemeris.constants import SE_EARTH, SE_MOON, SE_SUN
from libephemeris.leb_format import (
    BODY_PARAMS,
    MAGIC,
)
from libephemeris.leb_reader import LEBReader
from libephemeris.time_utils import swe_julday


class TestAssembleLeb:
    """Test the full LEB assembly pipeline."""

    @pytest.mark.integration
    def test_generates_valid_file(self, test_leb_file):
        """Generated .leb file should exist and be non-empty."""
        assert os.path.exists(test_leb_file)
        assert os.path.getsize(test_leb_file) > 0

    @pytest.mark.integration
    def test_file_has_correct_magic(self, test_leb_file):
        """Generated file should start with LEB1 magic bytes."""
        with open(test_leb_file, "rb") as f:
            magic = f.read(4)
        assert magic == MAGIC

    @pytest.mark.integration
    def test_reader_opens_generated_file(self, test_leb_file):
        """LEBReader should successfully parse the generated file."""
        reader = LEBReader(test_leb_file)
        assert reader is not None
        jd_start, jd_end = reader.jd_range
        assert jd_start < jd_end
        reader.close()

    @pytest.mark.integration
    def test_generated_file_has_expected_bodies(self, test_leb_file):
        """Generated file should contain the requested bodies."""
        with LEBReader(test_leb_file) as reader:
            # conftest requests: SE_SUN(0), SE_MOON(1), SE_MARS(4),
            # SE_EARTH(14), SE_MEAN_NODE(10)
            assert reader.has_body(SE_SUN)
            assert reader.has_body(SE_MOON)
            assert reader.has_body(SE_EARTH)

    @pytest.mark.integration
    def test_generated_file_has_nutation(self, test_leb_file):
        """Generated file should include nutation data."""
        with LEBReader(test_leb_file) as reader:
            assert reader._nutation is not None

    @pytest.mark.integration
    def test_generated_file_has_delta_t(self, test_leb_file):
        """Generated file should include Delta-T data."""
        with LEBReader(test_leb_file) as reader:
            assert len(reader._delta_t_jds) > 0

    @pytest.mark.integration
    def test_body_coord_types_correct(self, test_leb_file):
        """Body coordinate types should match BODY_PARAMS."""
        with LEBReader(test_leb_file) as reader:
            for body_id, body in reader._bodies.items():
                if body_id in BODY_PARAMS:
                    expected_coord = BODY_PARAMS[body_id][2]
                    assert body.coord_type == expected_coord, (
                        f"Body {body_id}: coord_type={body.coord_type}, "
                        f"expected={expected_coord}"
                    )


class TestGenerateMinimal:
    """Test minimal file generation (Sun + Earth only, 1 year)."""

    @pytest.mark.integration
    def test_minimal_file_works(self, test_leb_file_minimal):
        """Minimal .leb file should work correctly."""
        with LEBReader(test_leb_file_minimal) as reader:
            assert reader.has_body(SE_SUN)
            assert reader.has_body(SE_EARTH)
            # Should not have other bodies
            assert not reader.has_body(SE_MOON)

    @pytest.mark.integration
    def test_minimal_roundtrip(self, test_leb_file_minimal):
        """Sun position should survive generate→read roundtrip."""
        with LEBReader(test_leb_file_minimal) as reader:
            jd_start, jd_end = reader.jd_range
            jd_mid = (jd_start + jd_end) / 2.0

            pos, vel = reader.eval_body(SE_SUN, jd_mid)
            # Sun should be near SSB (< 0.02 AU typically)
            import math

            dist = math.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
            assert dist < 0.1, f"Sun at {dist} AU from SSB"


class TestGenerateSingleBody:
    """Test individual body generation."""

    @pytest.mark.integration
    @pytest.mark.slow
    def test_generate_sun(self, tmp_path):
        """Generate Sun-only LEB file and verify."""
        from scripts.generate_leb import assemble_leb

        path = str(tmp_path / "sun_only.leb")
        jd_start = swe_julday(2024, 1, 1, 0.0)
        jd_end = swe_julday(2024, 7, 1, 0.0)

        assemble_leb(
            output=path,
            jd_start=jd_start,
            jd_end=jd_end,
            bodies=[0],  # Sun only
            workers=1,
            verbose=False,
        )

        with LEBReader(path) as reader:
            assert reader.has_body(SE_SUN)
            jd_mid = (jd_start + jd_end) / 2.0
            pos, vel = reader.eval_body(SE_SUN, jd_mid)
            assert len(pos) == 3
            assert len(vel) == 3


class TestChebyshevFitting:
    """Test Chebyshev polynomial fitting quality."""

    @pytest.mark.integration
    def test_sun_fit_accuracy(self, test_leb_file):
        """Sun Chebyshev fit should reproduce Skyfield within 0.01 arcsec."""
        from libephemeris.planets import get_planet_target
        from libephemeris.state import get_planets, get_timescale

        with LEBReader(test_leb_file) as reader:
            jd_start, jd_end = reader.jd_range
            jd_mid = (jd_start + jd_end) / 2.0

            pos_leb, _ = reader.eval_body(SE_SUN, jd_mid)

            planets = get_planets()
            ts = get_timescale()
            target = get_planet_target(planets, "sun")
            t = ts.tt_jd(jd_mid)
            ref_pos = target.at(t).position.au

            max_err_arcsec = 0.0
            for c in range(3):
                err_au = abs(pos_leb[c] - float(ref_pos[c]))
                err_arcsec = err_au * 206265.0
                max_err_arcsec = max(max_err_arcsec, err_arcsec)

            assert max_err_arcsec < 0.01, f"Sun fit error = {max_err_arcsec:.4f} arcsec"

    @pytest.mark.integration
    def test_moon_fit_accuracy(self, test_leb_file):
        """Moon Chebyshev fit should reproduce Skyfield within 0.01 arcsec."""
        from libephemeris.planets import get_planet_target
        from libephemeris.state import get_planets, get_timescale

        with LEBReader(test_leb_file) as reader:
            jd_start, jd_end = reader.jd_range
            jd_mid = (jd_start + jd_end) / 2.0

            pos_leb, _ = reader.eval_body(SE_MOON, jd_mid)

            planets = get_planets()
            ts = get_timescale()
            # Moon position from Skyfield (geocentric ICRS)
            target = get_planet_target(planets, "moon")
            t = ts.tt_jd(jd_mid)
            ref_pos = target.at(t).position.au

            max_err_arcsec = 0.0
            for c in range(3):
                err_au = abs(pos_leb[c] - float(ref_pos[c]))
                err_arcsec = err_au * 206265.0
                max_err_arcsec = max(max_err_arcsec, err_arcsec)

            assert max_err_arcsec < 0.01, (
                f"Moon fit error = {max_err_arcsec:.4f} arcsec"
            )
