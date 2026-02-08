"""
Tests for libephemeris centralized logging configuration.

This module tests the logging_config module to ensure:
- Logger is properly configured
- Log messages have the correct format
- Logging can be enabled/disabled
- Log levels can be changed
"""

from __future__ import annotations

import io
import logging
import sys
from unittest.mock import patch

import pytest


class TestGetLogger:
    """Tests for get_logger() function."""

    def test_get_logger_returns_logger(self):
        """get_logger should return a logging.Logger instance."""
        from libephemeris.logging_config import get_logger, _logger_configured

        # Reset the logger configured flag for testing
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()
        assert isinstance(logger, logging.Logger)
        assert logger.name == "libephemeris"

    def test_get_logger_has_handler(self):
        """Logger should have at least one handler configured."""
        from libephemeris.logging_config import get_logger
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()
        assert len(logger.handlers) >= 1

    def test_get_logger_handler_is_stderr(self):
        """Logger handler should output to stderr."""
        from libephemeris.logging_config import get_logger
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()
        # At least one handler should be a StreamHandler to stderr
        stderr_handlers = [
            h
            for h in logger.handlers
            if isinstance(h, logging.StreamHandler) and h.stream == sys.stderr
        ]
        assert len(stderr_handlers) >= 1

    def test_get_logger_is_idempotent(self):
        """Calling get_logger multiple times should return the same logger."""
        from libephemeris.logging_config import get_logger
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger1 = get_logger()
        handler_count = len(logger1.handlers)

        logger2 = get_logger()
        # Should not add duplicate handlers
        assert len(logger2.handlers) == handler_count
        assert logger1 is logger2


class TestLogFormat:
    """Tests for log message formatting."""

    def test_log_format_contains_libephemeris_prefix(self):
        """Log messages should be prefixed with [libephemeris]."""
        from libephemeris.logging_config import get_logger
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        # Capture stderr output
        captured = io.StringIO()
        with patch.object(sys, "stderr", captured):
            # Create a new handler that writes to our captured stream
            test_handler = logging.StreamHandler(captured)
            test_handler.setLevel(logging.INFO)
            from libephemeris.logging_config import LibephemerisFormatter

            test_handler.setFormatter(LibephemerisFormatter())

            # Temporarily add our test handler
            logger.addHandler(test_handler)
            try:
                logger.info("Test message")
            finally:
                logger.removeHandler(test_handler)

        output = captured.getvalue()
        assert "[libephemeris]" in output

    def test_log_format_contains_level(self):
        """Log messages should contain the log level."""
        from libephemeris.logging_config import LibephemerisFormatter

        formatter = LibephemerisFormatter()
        record = logging.LogRecord(
            name="libephemeris",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="Test message",
            args=(),
            exc_info=None,
        )

        formatted = formatter.format(record)
        assert "INFO" in formatted
        assert "Test message" in formatted


class TestSetLogLevel:
    """Tests for set_log_level() function."""

    def test_set_log_level_changes_logger_level(self):
        """set_log_level should change the logger level."""
        from libephemeris.logging_config import get_logger, set_log_level
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        set_log_level(logging.DEBUG)
        assert logger.level == logging.DEBUG

        set_log_level(logging.WARNING)
        assert logger.level == logging.WARNING

        # Reset to default
        set_log_level(logging.INFO)


class TestDisableLogging:
    """Tests for disable_logging() function."""

    def test_disable_logging_silences_logger(self):
        """disable_logging should prevent all log output."""
        from libephemeris.logging_config import (
            get_logger,
            disable_logging,
            enable_logging,
        )
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        # Disable logging
        disable_logging()

        # The level should be above CRITICAL
        assert logger.level > logging.CRITICAL

        # Re-enable for other tests
        enable_logging()


class TestEnableLogging:
    """Tests for enable_logging() function."""

    def test_enable_logging_restores_default_level(self):
        """enable_logging should restore the default INFO level."""
        from libephemeris.logging_config import (
            get_logger,
            disable_logging,
            enable_logging,
            DEFAULT_LEVEL,
        )
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        disable_logging()
        enable_logging()

        assert logger.level == DEFAULT_LEVEL

    def test_enable_logging_with_custom_level(self):
        """enable_logging should accept a custom level."""
        from libephemeris.logging_config import get_logger, enable_logging
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        enable_logging(logging.DEBUG)
        assert logger.level == logging.DEBUG

        # Reset to default
        enable_logging()


class TestFormatFileSize:
    """Tests for format_file_size() function."""

    def test_format_bytes(self):
        """format_file_size should format bytes correctly."""
        from libephemeris.logging_config import format_file_size

        assert format_file_size(0) == "0.0 B"
        assert format_file_size(512) == "512.0 B"

    def test_format_kilobytes(self):
        """format_file_size should format kilobytes correctly."""
        from libephemeris.logging_config import format_file_size

        assert format_file_size(1024) == "1.0 KB"
        assert format_file_size(2048) == "2.0 KB"

    def test_format_megabytes(self):
        """format_file_size should format megabytes correctly."""
        from libephemeris.logging_config import format_file_size

        assert format_file_size(1024 * 1024) == "1.0 MB"
        # DE440 is about 114 MB
        result = format_file_size(119537664)
        assert "MB" in result
        assert float(result.split()[0]) > 100

    def test_format_gigabytes(self):
        """format_file_size should format gigabytes correctly."""
        from libephemeris.logging_config import format_file_size

        assert format_file_size(1024 * 1024 * 1024) == "1.0 GB"


class TestLoggingExportedFromInit:
    """Tests that logging functions are exported from libephemeris.__init__."""

    def test_get_logger_exported(self):
        """get_logger should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "get_logger")

    def test_set_log_level_exported(self):
        """set_log_level should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "set_log_level")

    def test_disable_logging_exported(self):
        """disable_logging should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "disable_logging")

    def test_enable_logging_exported(self):
        """enable_logging should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "enable_logging")

    def test_format_file_size_exported(self):
        """format_file_size should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "format_file_size")


class TestLoggerNameConstant:
    """Tests for the LOGGER_NAME constant."""

    def test_logger_name_is_libephemeris(self):
        """LOGGER_NAME should be 'libephemeris'."""
        from libephemeris.logging_config import LOGGER_NAME

        assert LOGGER_NAME == "libephemeris"
