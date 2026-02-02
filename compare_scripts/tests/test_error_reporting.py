"""
Tests for improved error reporting in run_all_compare.py.

Tests the analyze_error and format_error_report functions that provide
detailed error categorization, suggestions, and context for module failures.
"""

import sys
import os
from unittest.mock import patch, MagicMock
import pytest

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from run_all_compare import (
    analyze_error,
    format_error_report,
    ErrorInfo,
    ModuleResult,
    run_all_modules,
)


class TestAnalyzeError:
    """Tests for the analyze_error function."""

    def test_analyze_error_returns_none_for_empty_output(self):
        """analyze_error should return None for empty output."""
        assert analyze_error("") is None

    def test_analyze_error_returns_none_for_no_error(self):
        """analyze_error should return None when no error is detected."""
        output = """
        Running tests...
        Total tests:   10
        Passed:        10 ✓
        All tests passed!
        """
        assert analyze_error(output) is None

    def test_analyze_error_detects_import_error(self):
        """analyze_error should detect and categorize import errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 5, in <module>
            import missing_module
        ModuleNotFoundError: No module named 'missing_module'
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "import_error"
        assert "ModuleNotFoundError" in result.message
        assert "pip install" in result.suggestion

    def test_analyze_error_detects_general_import_error(self):
        """analyze_error should detect general ImportError."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 5, in <module>
            from package import something
        ImportError: cannot import name 'something' from 'package'
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "import_error"
        assert "PYTHONPATH" in result.suggestion

    def test_analyze_error_detects_syntax_error(self):
        """analyze_error should detect syntax errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 10
            print("hello"
                        ^
        SyntaxError: unexpected EOF while parsing
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "syntax_error"
        assert "SyntaxError" in result.message
        assert "syntax" in result.suggestion.lower()

    def test_analyze_error_detects_indentation_error(self):
        """analyze_error should detect indentation errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 15
            return result
            ^
        IndentationError: unexpected indent
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "syntax_error"
        assert "IndentationError" in result.message

    def test_analyze_error_detects_assertion_error(self):
        """analyze_error should detect assertion errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 20, in test_something
            assert result == expected
        AssertionError: Values differ: 1.23 != 1.24
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "assertion_error"
        assert "assertion" in result.suggestion.lower()

    def test_analyze_error_detects_file_not_found_error(self):
        """analyze_error should detect file not found errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 8, in <module>
            with open("data/ephemeris.dat") as f:
        FileNotFoundError: [Errno 2] No such file or directory: 'data/ephemeris.dat'
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "file_error"
        assert "ephemeris" in result.suggestion.lower()

    def test_analyze_error_detects_timeout(self):
        """analyze_error should detect timeout errors."""
        output = "TIMEOUT: Module exceeded 5 minute limit"
        result = analyze_error(output)
        assert result is not None
        assert result.category == "timeout_error"
        assert "time limit" in result.suggestion.lower()

    def test_analyze_error_detects_attribute_error(self):
        """analyze_error should detect attribute errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 12, in <module>
            obj.missing_method()
        AttributeError: 'MyClass' object has no attribute 'missing_method'
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "attribute_error"
        assert "API" in result.suggestion

    def test_analyze_error_detects_type_error(self):
        """analyze_error should detect type errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 15, in <module>
            func(1, 2, 3)
        TypeError: func() takes 2 positional arguments but 3 were given
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "type_error"
        assert "type" in result.suggestion.lower()

    def test_analyze_error_detects_value_error(self):
        """analyze_error should detect value errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 18, in <module>
            int("not_a_number")
        ValueError: invalid literal for int() with base 10: 'not_a_number'
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "value_error"

    def test_analyze_error_detects_key_error(self):
        """analyze_error should detect key errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 20, in <module>
            data["missing_key"]
        KeyError: 'missing_key'
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "key_error"
        assert "Dictionary" in result.suggestion

    def test_analyze_error_detects_memory_error(self):
        """analyze_error should detect memory errors."""
        output = """
        Traceback (most recent call last):
          File "compare_test.py", line 25, in <module>
            large_list = [0] * 10**12
        MemoryError
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "resource_error"
        assert "memory" in result.suggestion.lower()

    def test_analyze_error_extracts_last_10_lines(self):
        """analyze_error should extract the last 10 non-empty lines."""
        lines = [f"Line {i}" for i in range(20)]
        lines.append("Traceback (most recent call last):")
        lines.append("ERROR: Something went wrong")
        output = "\n".join(lines)

        result = analyze_error(output)
        assert result is not None
        assert len(result.last_lines) <= 10
        assert "ERROR" in result.last_lines[-1]

    def test_analyze_error_handles_unknown_error(self):
        """analyze_error should handle unknown error types."""
        output = """
        ERROR: Something unexpected happened
        This is an unusual error message
        """
        result = analyze_error(output)
        assert result is not None
        assert result.category == "unknown_error"
        assert "Review" in result.suggestion


class TestFormatErrorReport:
    """Tests for the format_error_report function."""

    def test_format_error_report_includes_category(self):
        """format_error_report should include the error category."""
        error_info = ErrorInfo(
            category="import_error",
            message="ModuleNotFoundError: No module named 'missing'",
            suggestion="Install the missing module",
            last_lines=["line 1", "line 2"],
        )
        result = format_error_report(error_info)
        assert "Import Error" in result  # Title case conversion

    def test_format_error_report_includes_message(self):
        """format_error_report should include the error message."""
        error_info = ErrorInfo(
            category="syntax_error",
            message="SyntaxError: invalid syntax",
            suggestion="Fix the syntax",
            last_lines=["error line"],
        )
        result = format_error_report(error_info)
        assert "SyntaxError: invalid syntax" in result

    def test_format_error_report_includes_suggestion(self):
        """format_error_report should include the suggestion."""
        error_info = ErrorInfo(
            category="file_error",
            message="FileNotFoundError: file.txt",
            suggestion="Check file paths and permissions",
            last_lines=["error line"],
        )
        result = format_error_report(error_info)
        assert "Check file paths and permissions" in result

    def test_format_error_report_includes_last_lines(self):
        """format_error_report should include the last lines of output."""
        last_lines = [
            "Running test...",
            "Test started",
            "ERROR: Test failed",
        ]
        error_info = ErrorInfo(
            category="assertion_error",
            message="AssertionError",
            suggestion="Check assertions",
            last_lines=last_lines,
        )
        result = format_error_report(error_info)
        for line in last_lines:
            assert line in result

    def test_format_error_report_structure(self):
        """format_error_report should have proper structure."""
        error_info = ErrorInfo(
            category="type_error",
            message="TypeError: wrong type",
            suggestion="Check argument types",
            last_lines=["line 1", "line 2"],
        )
        result = format_error_report(error_info)

        # Check structure
        assert "Error Type:" in result
        assert "Message:" in result
        assert "Suggestion:" in result
        assert "Last 10 lines" in result
        assert "-" * 40 in result


class TestErrorInfoDataclass:
    """Tests for the ErrorInfo dataclass."""

    def test_error_info_creation(self):
        """ErrorInfo should be creatable with all fields."""
        error_info = ErrorInfo(
            category="test_error",
            message="Test message",
            suggestion="Test suggestion",
            last_lines=["line 1", "line 2"],
        )
        assert error_info.category == "test_error"
        assert error_info.message == "Test message"
        assert error_info.suggestion == "Test suggestion"
        assert error_info.last_lines == ["line 1", "line 2"]


class TestModuleResultWithErrorInfo:
    """Tests for ModuleResult with ErrorInfo integration."""

    def test_module_result_with_error_info(self):
        """ModuleResult should accept error_info parameter."""
        error_info = ErrorInfo(
            category="import_error",
            message="Import failed",
            suggestion="Install module",
            last_lines=["error"],
        )
        result = ModuleResult(
            name="test",
            passed=0,
            total=0,
            duration=1.0,
            error="Import failed",
            error_info=error_info,
        )
        assert result.error_info is not None
        assert result.error_info.category == "import_error"

    def test_module_result_without_error_info(self):
        """ModuleResult should work without error_info."""
        result = ModuleResult(
            name="test",
            passed=10,
            total=10,
            duration=1.0,
        )
        assert result.error_info is None
        assert result.status == "PASS"


class TestRunAllModulesErrorReporting:
    """Tests for error reporting in run_all_modules."""

    def test_run_all_modules_captures_error_info(self, capsys):
        """run_all_modules should capture and display error info."""

        def mock_run_module(script, verbose=False, progress_interval=30):
            return (
                0,
                0,
                """
                Traceback (most recent call last):
                  File "test.py", line 1
                    import missing_module
                ModuleNotFoundError: No module named 'missing_module'
                """,
            )

        with patch("run_all_compare.run_module", side_effect=mock_run_module):
            test_modules = [
                ("test", "test_script.py", "Test module", True),
            ]
            results = run_all_modules(test_modules, verbose=False)

        assert len(results) == 1
        result = results[0]
        assert result.error is not None
        assert result.error_info is not None
        assert result.error_info.category == "import_error"

        captured = capsys.readouterr()
        assert "Import Error" in captured.out
        assert "pip install" in captured.out

    def test_run_all_modules_no_error_info_for_success(self):
        """run_all_modules should not create error_info for successful runs."""

        def mock_run_module(script, verbose=False, progress_interval=30):
            return (10, 10, "Total tests:   10\nPassed:        10 ✓\n")

        with patch("run_all_compare.run_module", side_effect=mock_run_module):
            test_modules = [
                ("test", "test_script.py", "Test module", True),
            ]
            results = run_all_modules(test_modules, verbose=False)

        assert len(results) == 1
        result = results[0]
        assert result.error is None
        assert result.error_info is None

    def test_run_all_modules_displays_last_lines(self, capsys):
        """run_all_modules should display last 10 lines on error."""

        output_lines = [f"Log line {i}" for i in range(15)]
        output_lines.append("Traceback (most recent call last):")
        output_lines.append("  File 'test.py', line 10")
        output_lines.append("AssertionError: Test failed")
        full_output = "\n".join(output_lines)

        def mock_run_module(script, verbose=False, progress_interval=30):
            return (0, 0, full_output)

        with patch("run_all_compare.run_module", side_effect=mock_run_module):
            test_modules = [
                ("test", "test_script.py", "Test module", True),
            ]
            results = run_all_modules(test_modules, verbose=False)

        captured = capsys.readouterr()
        # Should show last lines
        assert "Last 10 lines" in captured.out
        assert "AssertionError" in captured.out


class TestErrorPatternCoverage:
    """Tests to ensure all error patterns are correctly matched."""

    @pytest.mark.parametrize(
        "error_text,expected_category",
        [
            ("ModuleNotFoundError: No module", "import_error"),
            ("ImportError: cannot import", "import_error"),
            ("SyntaxError: invalid syntax", "syntax_error"),
            ("IndentationError: unexpected", "syntax_error"),
            ("AssertionError: test failed", "assertion_error"),
            ("FileNotFoundError: file.txt", "file_error"),
            ("PermissionError: access denied", "file_error"),
            ("TimeoutError: operation timed out", "timeout_error"),
            ("TIMEOUT: Module exceeded", "timeout_error"),
            ("MemoryError", "resource_error"),
            ("RecursionError: maximum recursion", "resource_error"),
            ("AttributeError: 'obj' has no", "attribute_error"),
            ("TypeError: unsupported operand", "type_error"),
            ("ValueError: invalid literal", "value_error"),
            ("KeyError: 'missing_key'", "key_error"),
            ("ZeroDivisionError: division by zero", "math_error"),
            ("ConnectionError: connection refused", "network_error"),
        ],
    )
    def test_error_pattern_categorization(self, error_text, expected_category):
        """Test that each error pattern is correctly categorized."""
        output = f"Traceback (most recent call last):\n{error_text}"
        result = analyze_error(output)
        assert result is not None
        assert result.category == expected_category
