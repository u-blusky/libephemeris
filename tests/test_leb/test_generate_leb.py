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
            # ICRS barycentric: Sun is near SSB, typically <0.01 AU
            import math

            dist = math.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
            assert dist < 0.02, f"Sun barycentric distance = {dist} AU (too large)"


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
        """Sun through full pipeline should match Skyfield within 5 arcsec.

        Note: Tolerance is generous (5") for on-the-fly test files with default
        segment parameters.  Production accuracy validated by compare/ tests.
        """
        import libephemeris as ephem
        from libephemeris.constants import SEFLG_SPEED
        from libephemeris.fast_calc import fast_calc_ut

        with LEBReader(test_leb_file) as reader:
            jd_start, jd_end = reader.jd_range
            jd_mid = (jd_start + jd_end) / 2.0

            # Full pipeline: LEB → fast_calc → geocentric ecliptic
            fast_result, _ = fast_calc_ut(reader, jd_mid, SE_SUN, SEFLG_SPEED)

            # Skyfield reference (geocentric ecliptic of date)
            ref, _ = ephem.swe_calc_ut(jd_mid, SE_SUN, SEFLG_SPEED)

            lon_err = abs(fast_result[0] - ref[0])
            if lon_err > 180.0:
                lon_err = 360.0 - lon_err
            lon_err_arcsec = lon_err * 3600.0

            lat_err_arcsec = abs(fast_result[1] - ref[1]) * 3600.0

            max_err_arcsec = max(lon_err_arcsec, lat_err_arcsec)
            assert max_err_arcsec < 5.0, f"Sun fit error = {max_err_arcsec:.4f} arcsec"

    @pytest.mark.integration
    def test_moon_fit_accuracy(self, test_leb_file):
        """Moon through full pipeline should match Skyfield within 5 arcsec.

        Note: Tolerance is generous (5") for on-the-fly test files with default
        segment parameters.  Production accuracy validated by compare/ tests.
        """
        import libephemeris as ephem
        from libephemeris.constants import SEFLG_SPEED
        from libephemeris.fast_calc import fast_calc_ut

        with LEBReader(test_leb_file) as reader:
            jd_start, jd_end = reader.jd_range
            jd_mid = (jd_start + jd_end) / 2.0

            # Full pipeline: LEB → fast_calc → geocentric ecliptic
            fast_result, _ = fast_calc_ut(reader, jd_mid, SE_MOON, SEFLG_SPEED)

            # Skyfield reference (geocentric ecliptic of date)
            ref, _ = ephem.swe_calc_ut(jd_mid, SE_MOON, SEFLG_SPEED)

            lon_err = abs(fast_result[0] - ref[0])
            if lon_err > 180.0:
                lon_err = 360.0 - lon_err
            lon_err_arcsec = lon_err * 3600.0

            lat_err_arcsec = abs(fast_result[1] - ref[1]) * 3600.0

            max_err_arcsec = max(lon_err_arcsec, lat_err_arcsec)
            assert max_err_arcsec < 5.0, f"Moon fit error = {max_err_arcsec:.4f} arcsec"
