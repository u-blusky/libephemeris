"""
Tests for auto-download SPK in strict precision mode.

Tests verify:
- When auto_spk_download is enabled and strict_precision is True,
  auto-download is attempted before raising SPKRequiredError
- If auto-download succeeds, calculation works
- If auto-download fails, SPKRequiredError is still raised
- Logging messages are produced during auto-download
"""

import pytest
import logging
from unittest.mock import patch, MagicMock
import libephemeris as eph
from libephemeris import state
from libephemeris.constants import (
    SE_CHIRON,
    SE_CERES,
    SEFLG_SPEED,
)


@pytest.fixture(autouse=True)
def cleanup():
    """Reset state before and after each test."""
    # Clear any SPK registrations
    state._SPK_BODY_MAP.clear()
    # Reset settings
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)
    yield
    # Cleanup after test
    state._SPK_BODY_MAP.clear()
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)


class TestAutoDownloadInStrictMode:
    """Test auto-download behavior when strict precision is enabled."""

    def test_auto_download_attempted_before_error(self):
        """Auto-download should be attempted before raising SPKRequiredError."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch(
            "libephemeris.minor_bodies.ensure_major_asteroid_spk"
        ) as mock_ensure:
            mock_ensure.return_value = False  # Simulate download failure

            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)

            # Verify that ensure_major_asteroid_spk was called
            mock_ensure.assert_called_once()
            # First arg should be SE_CHIRON
            assert mock_ensure.call_args[0][0] == SE_CHIRON

    def test_auto_download_success_returns_result(self):
        """If auto-download succeeds, should return calculation result."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        # Mock ensure_major_asteroid_spk to return True (success)
        # and spk.calc_spk_body_position to return a valid result
        mock_position = (100.0, 5.0, 15.0, 0.1, 0.01, 0.001)

        with patch(
            "libephemeris.minor_bodies.ensure_major_asteroid_spk"
        ) as mock_ensure:
            mock_ensure.return_value = True

            # We need to mock spk.calc_spk_body_position to return a value
            # on the second call (after ensure_major_asteroid_spk succeeds)
            call_count = [0]

            def mock_calc_spk(t, ipl, iflag):
                call_count[0] += 1
                if call_count[0] == 1:
                    # First call returns None (SPK not registered yet)
                    return None
                else:
                    # Second call after ensure returns the position
                    return mock_position

            with patch(
                "libephemeris.spk.calc_spk_body_position", side_effect=mock_calc_spk
            ):
                pos, flags = eph.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)

                # Should have returned our mock position
                assert pos[0] == pytest.approx(100.0)
                assert pos[1] == pytest.approx(5.0)
                assert pos[2] == pytest.approx(15.0)

    def test_no_auto_download_when_disabled(self):
        """Auto-download should not be attempted when disabled."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(False)

        with patch(
            "libephemeris.minor_bodies.ensure_major_asteroid_spk"
        ) as mock_ensure:
            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)

            # ensure_major_asteroid_spk should NOT be called
            mock_ensure.assert_not_called()

    def test_auto_download_logs_info_message(self):
        """Auto-download should log an info message."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch(
            "libephemeris.minor_bodies.ensure_major_asteroid_spk"
        ) as mock_ensure:
            mock_ensure.return_value = False

            # Mock get_logger from logging_config (where it's defined)
            with patch("libephemeris.logging_config.get_logger") as mock_get_logger:
                mock_logger = MagicMock()
                mock_get_logger.return_value = mock_logger

                with pytest.raises(eph.SPKRequiredError):
                    eph.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)

                # Check that info logging was called with the body name
                mock_logger.info.assert_called_once()
                log_call_args = mock_logger.info.call_args[0]
                assert "Auto-downloading" in log_call_args[0]
                assert "Chiron" in str(log_call_args)

    def test_download_exception_still_raises_spk_error(self):
        """If download raises exception, SPKRequiredError should still be raised."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch(
            "libephemeris.minor_bodies.ensure_major_asteroid_spk"
        ) as mock_ensure:
            mock_ensure.side_effect = RuntimeError("Network error")

            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)

    def test_multiple_bodies_auto_download(self):
        """Auto-download should work for multiple major asteroids."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        bodies_called = []

        def track_calls(body_id, jd=None):
            bodies_called.append(body_id)
            return False  # Fail download

        with patch(
            "libephemeris.minor_bodies.ensure_major_asteroid_spk",
            side_effect=track_calls,
        ):
            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)

            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, SE_CERES, SEFLG_SPEED)

        assert SE_CHIRON in bodies_called
        assert SE_CERES in bodies_called


class TestAutoDownloadWithStrictDisabled:
    """Test that auto-download path is used with strict mode disabled."""

    def test_keplerian_fallback_when_strict_disabled(self):
        """When strict is disabled, should fall back to Keplerian even if download fails."""
        eph.set_strict_precision(False)
        eph.set_auto_spk_download(True)

        # Even without SPK, should get a result via Keplerian fallback
        pos, flags = eph.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)

        # Should return valid position
        assert 0 <= pos[0] < 360  # Longitude in range
        assert -90 <= pos[1] <= 90  # Latitude in range
        assert pos[2] > 0  # Distance positive


class TestAutoDownloadWithFirstTryPath:
    """Test that the first try_auto_spk_download path still works."""

    def test_first_try_path_still_attempted(self):
        """The initial _try_auto_spk_download should still be called."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        # Track if _try_auto_spk_download is called
        with patch("libephemeris.planets._try_auto_spk_download") as mock_try:
            mock_try.return_value = None  # Return None to continue to strict check

            with patch(
                "libephemeris.minor_bodies.ensure_major_asteroid_spk"
            ) as mock_ensure:
                mock_ensure.return_value = False

                with pytest.raises(eph.SPKRequiredError):
                    eph.calc_ut(2451545.0, SE_CHIRON, SEFLG_SPEED)

                # Both paths should have been attempted
                mock_try.assert_called_once()
                mock_ensure.assert_called_once()


class TestAutoDownloadJDPassthrough:
    """Test that Julian Day is passed correctly to ensure_major_asteroid_spk."""

    def test_jd_passed_to_ensure(self):
        """The Julian Day should be passed to ensure_major_asteroid_spk."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        test_jd = 2460000.5  # A specific JD

        with patch(
            "libephemeris.minor_bodies.ensure_major_asteroid_spk"
        ) as mock_ensure:
            mock_ensure.return_value = False

            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(test_jd, SE_CHIRON, SEFLG_SPEED)

            # Check that jd was passed (second argument)
            assert mock_ensure.call_args[0][0] == SE_CHIRON
            # The jd is passed as t.tt which will be close to test_jd
            # (there's a small UT->TT conversion difference)
            passed_jd = mock_ensure.call_args[0][1]
            assert abs(passed_jd - test_jd) < 0.001  # Within ~1 minute
