"""
Tests for progress reporting in run_all_compare.py.

Tests the run_module function's progress reporting feature that provides
feedback every 30 seconds during subprocess execution.
"""

import sys
import os
import time
from unittest.mock import patch, MagicMock
import pytest

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from run_all_compare import run_module, run_all_modules, TQDM_AVAILABLE, ModuleResult


class TestRunModuleProgressReporting:
    """Tests for run_module progress reporting."""

    def test_run_module_returns_tuple(self):
        """run_module should return a tuple of (passed, total, output)."""
        # Mock subprocess to simulate a quick script execution
        mock_process = MagicMock()
        mock_process.stdout = MagicMock()
        mock_process.stderr = MagicMock()
        mock_process.stdout.readline = MagicMock(
            side_effect=["Total tests:   10\n", "Passed:        8 ✓\n", ""]
        )
        mock_process.stderr.readline = MagicMock(return_value="")
        mock_process.poll = MagicMock(side_effect=[None, None, None, 0])

        with patch("run_all_compare.subprocess.Popen", return_value=mock_process):
            with patch("run_all_compare.time.sleep"):
                passed, total, output = run_module("test_script.py", verbose=False)

        assert isinstance(passed, int)
        assert isinstance(total, int)
        assert isinstance(output, str)

    def test_run_module_accepts_progress_interval(self):
        """run_module should accept progress_interval parameter."""
        mock_process = MagicMock()
        mock_process.stdout = MagicMock()
        mock_process.stderr = MagicMock()
        mock_process.stdout.readline = MagicMock(return_value="")
        mock_process.stderr.readline = MagicMock(return_value="")
        mock_process.poll = MagicMock(return_value=0)

        with patch("run_all_compare.subprocess.Popen", return_value=mock_process):
            with patch("run_all_compare.time.sleep"):
                # Should not raise any exceptions
                passed, total, output = run_module(
                    "test_script.py", verbose=False, progress_interval=60
                )

        assert isinstance(passed, int)

    def test_run_module_parses_output_correctly(self):
        """run_module should correctly parse passed/total from output."""
        mock_process = MagicMock()
        mock_process.stdout = MagicMock()
        mock_process.stderr = MagicMock()
        mock_process.stdout.readline = MagicMock(
            side_effect=[
                "Running tests...\n",
                "Total tests:   39\n",
                "Passed:        24 ✓\n",
                "",
            ]
        )
        mock_process.stderr.readline = MagicMock(return_value="")
        mock_process.poll = MagicMock(side_effect=[None, None, None, None, 0])

        with patch("run_all_compare.subprocess.Popen", return_value=mock_process):
            with patch("run_all_compare.time.sleep"):
                passed, total, output = run_module("test_script.py", verbose=False)

        assert total == 39, "Should have parsed total tests"
        assert passed == 24, "Should have parsed passed tests"
        assert passed <= total, "Passed should not exceed total"

    def test_run_module_handles_verbose_flag(self):
        """run_module should pass verbose flag to subprocess."""
        captured_cmds = []

        def capture_popen(cmd, **kwargs):
            captured_cmds.append(cmd)
            mock_process = MagicMock()
            mock_process.stdout = MagicMock()
            mock_process.stderr = MagicMock()
            mock_process.stdout.readline = MagicMock(return_value="")
            mock_process.stderr.readline = MagicMock(return_value="")
            mock_process.poll = MagicMock(return_value=0)
            return mock_process

        with patch("run_all_compare.subprocess.Popen", side_effect=capture_popen):
            with patch("run_all_compare.time.sleep"):
                run_module("test_script.py", verbose=False)
                run_module("test_script.py", verbose=True)

        # Verbose run should have --verbose flag
        assert "--verbose" not in captured_cmds[0]
        assert "--verbose" in captured_cmds[1]


class TestProgressOutput:
    """Tests for progress output during slow operations."""

    def test_progress_prints_on_interval(self, capsys):
        """Progress should be printed at specified intervals during long runs."""
        mock_process = MagicMock()
        mock_process.stdout = MagicMock()
        mock_process.stderr = MagicMock()
        mock_process.stdout.readline = MagicMock(return_value="")
        mock_process.stderr.readline = MagicMock(return_value="")
        mock_process.kill = MagicMock()

        poll_count = [0]

        def poll_side_effect():
            poll_count[0] += 1
            if poll_count[0] < 20:
                return None
            return 0

        mock_process.poll = poll_side_effect

        # Mock time to progress 35 seconds (past the 30 second interval)
        time_values = [0.0]

        def mock_time():
            result = time_values[0]
            time_values[0] += 2.0  # 2 seconds per call
            return result

        with patch("run_all_compare.subprocess.Popen", return_value=mock_process):
            with patch("run_all_compare.time.time", side_effect=mock_time):
                with patch("run_all_compare.time.sleep"):
                    run_module("test_script.py", verbose=False, progress_interval=30)

        captured = capsys.readouterr()
        # Should have printed progress with elapsed/remaining info
        assert "Running test_script.py" in captured.out
        assert "elapsed" in captured.out

    def test_progress_shows_spinner(self, capsys):
        """Progress should show spinner characters."""
        mock_process = MagicMock()
        mock_process.stdout = MagicMock()
        mock_process.stderr = MagicMock()
        mock_process.stdout.readline = MagicMock(return_value="")
        mock_process.stderr.readline = MagicMock(return_value="")

        poll_count = [0]

        def poll_side_effect():
            poll_count[0] += 1
            if poll_count[0] < 10:
                return None
            return 0

        mock_process.poll = poll_side_effect

        # Mock time to trigger multiple progress reports
        time_values = [0.0]

        def mock_time():
            result = time_values[0]
            time_values[0] += 35.0  # Each call jumps 35 seconds
            return result

        with patch("run_all_compare.subprocess.Popen", return_value=mock_process):
            with patch("run_all_compare.time.time", side_effect=mock_time):
                with patch("run_all_compare.time.sleep"):
                    run_module("test_script.py", verbose=False, progress_interval=30)

        captured = capsys.readouterr()
        # Should have spinner character (Braille pattern)
        assert any(c in captured.out for c in "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")


class TestRunAllModulesProgress:
    """Tests for run_all_modules with progress bar."""

    def test_run_all_modules_shows_module_count(self, capsys):
        """run_all_modules should show [1/N] style progress."""

        def mock_run_module(script, verbose=False, progress_interval=30):
            return (10, 10, "Total tests:   10\nPassed:        10 ✓\n")

        with patch("run_all_compare.run_module", side_effect=mock_run_module):
            test_modules = [
                ("test", "test_script.py", "Test module", True),
            ]
            results = run_all_modules(test_modules, verbose=False)

        captured = capsys.readouterr()
        assert "[1/1]" in captured.out

    def test_run_all_modules_shows_progress_for_multiple(self, capsys):
        """run_all_modules should show progress for multiple modules."""

        def mock_run_module(script, verbose=False, progress_interval=30):
            return (5, 5, "Total tests:   5\nPassed:        5 ✓\n")

        with patch("run_all_compare.run_module", side_effect=mock_run_module):
            test_modules = [
                ("mod1", "mod1.py", "Module 1", True),
                ("mod2", "mod2.py", "Module 2", True),
                ("mod3", "mod3.py", "Module 3", True),
            ]
            results = run_all_modules(test_modules, verbose=False)

        captured = capsys.readouterr()
        assert "[1/3]" in captured.out
        assert "[2/3]" in captured.out
        assert "[3/3]" in captured.out

    def test_run_all_modules_returns_module_results(self):
        """run_all_modules should return list of ModuleResult objects."""

        def mock_run_module(script, verbose=False, progress_interval=30):
            return (8, 10, "Total tests:   10\nPassed:        8 ✓\n")

        with patch("run_all_compare.run_module", side_effect=mock_run_module):
            test_modules = [
                ("test", "test_script.py", "Test module", True),
            ]
            results = run_all_modules(test_modules, verbose=False)

        assert len(results) == 1
        assert isinstance(results[0], ModuleResult)
        assert results[0].name == "test"
        assert results[0].passed == 8
        assert results[0].total == 10

    def test_tqdm_available_flag_exists(self):
        """TQDM_AVAILABLE flag should exist."""
        assert isinstance(TQDM_AVAILABLE, bool)


class TestTimeoutHandling:
    """Tests for timeout handling with progress reporting."""

    def test_run_module_handles_timeout(self):
        """run_module should handle timeout gracefully."""
        mock_process = MagicMock()
        mock_process.stdout = MagicMock()
        mock_process.stderr = MagicMock()
        mock_process.stdout.readline = MagicMock(return_value="")
        mock_process.stderr.readline = MagicMock(return_value="")
        mock_process.poll = MagicMock(return_value=None)  # Never completes
        mock_process.kill = MagicMock()

        # Mock time to simulate timeout (jump 400 seconds)
        time_values = [0]

        def mock_time():
            result = time_values[0]
            time_values[0] += 100  # Jump 100 seconds each call
            return result

        with patch("run_all_compare.subprocess.Popen", return_value=mock_process):
            with patch("run_all_compare.time.time", side_effect=mock_time):
                with patch("run_all_compare.time.sleep"):
                    passed, total, output = run_module("test.py", verbose=False)

        assert passed == 0
        assert total == 0
        assert "TIMEOUT" in output
        mock_process.kill.assert_called_once()

    def test_run_module_timeout_message_includes_limit(self):
        """Timeout message should mention the 5 minute limit."""
        mock_process = MagicMock()
        mock_process.stdout = MagicMock()
        mock_process.stderr = MagicMock()
        mock_process.stdout.readline = MagicMock(return_value="")
        mock_process.stderr.readline = MagicMock(return_value="")
        mock_process.poll = MagicMock(return_value=None)
        mock_process.kill = MagicMock()

        time_values = [0]

        def mock_time():
            result = time_values[0]
            time_values[0] += 400
            return result

        with patch("run_all_compare.subprocess.Popen", return_value=mock_process):
            with patch("run_all_compare.time.time", side_effect=mock_time):
                with patch("run_all_compare.time.sleep"):
                    _, _, output = run_module("test.py", verbose=False)

        assert "5 minute" in output


class TestModuleResultDataclass:
    """Tests for ModuleResult dataclass."""

    def test_module_result_pass_rate(self):
        """ModuleResult should calculate pass rate correctly."""
        result = ModuleResult(name="test", passed=80, total=100, duration=1.0)
        assert result.pass_rate == 80.0

    def test_module_result_pass_rate_zero_total(self):
        """ModuleResult should handle zero total gracefully."""
        result = ModuleResult(name="test", passed=0, total=0, duration=1.0)
        assert result.pass_rate == 0.0

    def test_module_result_status_pass(self):
        """ModuleResult status should be PASS when all tests pass."""
        result = ModuleResult(name="test", passed=10, total=10, duration=1.0)
        assert result.status == "PASS"

    def test_module_result_status_fail(self):
        """ModuleResult status should be FAIL when some tests fail."""
        result = ModuleResult(name="test", passed=8, total=10, duration=1.0)
        assert result.status == "FAIL"

    def test_module_result_status_error(self):
        """ModuleResult status should be ERROR when error is set."""
        result = ModuleResult(
            name="test", passed=0, total=0, duration=1.0, error="Failed to run"
        )
        assert result.status == "ERROR"


class TestExceptionHandling:
    """Tests for exception handling in run_module."""

    def test_run_module_handles_popen_exception(self):
        """run_module should handle Popen exceptions gracefully."""
        with patch(
            "run_all_compare.subprocess.Popen",
            side_effect=OSError("Failed to start process"),
        ):
            passed, total, output = run_module("test.py", verbose=False)

        assert passed == 0
        assert total == 0
        assert "ERROR" in output

    def test_run_module_handles_generic_exception(self):
        """run_module should handle generic exceptions gracefully."""
        with patch(
            "run_all_compare.subprocess.Popen",
            side_effect=RuntimeError("Unexpected error"),
        ):
            passed, total, output = run_module("test.py", verbose=False)

        assert passed == 0
        assert total == 0
        assert "ERROR" in output
        assert "Unexpected error" in output
