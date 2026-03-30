"""Tests for TAI time functions: UTC↔TAI, TT↔TAI round-trips."""

from __future__ import annotations

import math
import pytest
import libephemeris as swe

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestUtcToTaiRoundTrip:
    """Test UTC -> TAI -> UTC round-trips."""

    @pytest.mark.parametrize(
        "year,month,day,hour,minute,second",
        [
            (2000, 1, 1, 12, 0, 0.0),
            (2020, 6, 15, 6, 30, 45.5),
            (2010, 3, 1, 0, 0, 0.0),
            (1990, 12, 31, 23, 59, 30.0),
            (2024, 7, 4, 18, 15, 22.7),
        ],
    )
    def test_utc_tai_utc_round_trip(self, year, month, day, hour, minute, second):
        """UTC -> TAI -> UTC should recover original date components."""
        jd_tai = swe.utc_to_tai_jd(year, month, day, hour, minute, second)
        assert math.isfinite(jd_tai), f"TAI JD is not finite: {jd_tai}"

        y2, m2, d2, h2, min2, sec2 = swe.tai_jd_to_utc(jd_tai)
        assert y2 == year
        assert m2 == month
        assert d2 == day
        assert h2 == hour
        assert min2 == minute
        assert abs(sec2 - second) < 0.001, f"Second mismatch: {sec2} vs {second}"

    def test_tai_ahead_of_utc(self):
        """TAI JD should be ahead of UTC JD (by leap seconds)."""
        jd_tai = swe.utc_to_tai_jd(2020, 1, 1, 0, 0, 0.0)
        jd_utc = swe.julday(2020, 1, 1, 0.0)
        # TAI = UTC + leap_seconds. In 2020, TAI-UTC = 37s
        diff_seconds = (jd_tai - jd_utc) * 86400.0
        assert 30.0 < diff_seconds < 45.0, f"TAI-UTC = {diff_seconds}s, expected ~37s"


@pytest.mark.unit
class TestTtToTaiRoundTrip:
    """Test TT -> TAI -> TT round-trips."""

    @pytest.mark.parametrize(
        "jd_tt",
        [
            JD_J2000,
            JD_J2000 + 3652.5,  # ~2010
            JD_J2000 - 3652.5,  # ~1990
            JD_J2000 + 7305.0,  # ~2020
            JD_J2000 + 18262.5,  # ~2050
        ],
    )
    def test_tt_tai_tt_round_trip(self, jd_tt):
        """TT -> TAI -> TT should recover the original JD exactly."""
        jd_tai = swe.tt_to_tai_jd(jd_tt)
        jd_tt_out = swe.tai_to_tt_jd(jd_tai)
        assert abs(jd_tt_out - jd_tt) < 1e-12, (
            f"TT round-trip error: {abs(jd_tt_out - jd_tt)} days"
        )

    def test_tt_tai_offset_is_32184ms(self):
        """TT - TAI should be exactly 32.184 seconds."""
        jd_tt = JD_J2000
        jd_tai = swe.tt_to_tai_jd(jd_tt)
        diff_seconds = (jd_tt - jd_tai) * 86400.0
        assert abs(diff_seconds - 32.184) < 0.001, (
            f"TT-TAI = {diff_seconds}s, expected 32.184s"
        )


@pytest.mark.unit
class TestTaiUtcJdConversion:
    """Test tai_to_utc_jd and utc_jd_to_tai conversions."""

    @pytest.mark.parametrize(
        "jd_utc_approx",
        [
            JD_J2000,
            JD_J2000 + 3652.5,
            JD_J2000 + 7305.0,
            JD_J2000 - 7305.0,
            JD_J2000 + 10957.5,
        ],
    )
    def test_utc_tai_utc_jd_round_trip(self, jd_utc_approx):
        """UTC JD -> TAI JD -> UTC JD should recover within ~1 microsecond."""
        jd_tai = swe.utc_jd_to_tai(jd_utc_approx)
        jd_utc_out = swe.tai_to_utc_jd(jd_tai)
        error_us = abs(jd_utc_out - jd_utc_approx) * 86400e6
        assert error_us < 1.0, f"UTC-TAI-UTC round-trip error: {error_us} µs"

    def test_tai_utc_jd_offset_positive(self):
        """TAI JD should be ahead of UTC JD by a positive number of seconds."""
        jd_utc = JD_J2000
        jd_tai = swe.utc_jd_to_tai(jd_utc)
        assert jd_tai > jd_utc, "TAI should be ahead of UTC"


@pytest.mark.unit
class TestGetTaiUtc:
    """Test get_tai_utc_for_jd leap second lookup."""

    def test_j2000_leap_seconds(self):
        """At J2000 (2000-01-01), TAI-UTC should be 32s."""
        leap = swe.get_tai_utc_for_jd(JD_J2000)
        assert abs(leap - 32.0) < 0.5, f"TAI-UTC at J2000: {leap}s, expected 32"

    def test_2017_leap_seconds(self):
        """After 2017-01-01, TAI-UTC should be 37s."""
        jd_2020 = swe.julday(2020, 1, 1, 0.0)
        leap = swe.get_tai_utc_for_jd(jd_2020)
        assert abs(leap - 37.0) < 0.5, f"TAI-UTC in 2020: {leap}s, expected 37"

    def test_leap_seconds_increase_over_time(self):
        """TAI-UTC should generally increase over decades."""
        jd_1980 = swe.julday(1980, 1, 1, 0.0)
        jd_2020 = swe.julday(2020, 1, 1, 0.0)
        leap_1980 = swe.get_tai_utc_for_jd(jd_1980)
        leap_2020 = swe.get_tai_utc_for_jd(jd_2020)
        assert leap_2020 > leap_1980, (
            f"Leap seconds should increase: {leap_1980} (1980) vs {leap_2020} (2020)"
        )

    def test_leap_seconds_always_positive(self):
        """TAI-UTC should be positive for modern dates."""
        for year in [1975, 1985, 1995, 2005, 2015, 2025]:
            jd = swe.julday(year, 6, 1, 0.0)
            leap = swe.get_tai_utc_for_jd(jd)
            assert leap > 0, f"TAI-UTC at {year}: {leap}s (expected positive)"
