"""
Tests for the Error exception class.

This ensures that libephemeris.Error is compatible with swisseph.Error (swe.Error).
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.exceptions import Error


class TestErrorClass:
    """Test the Error exception class."""

    @pytest.mark.unit
    def test_error_is_exception(self):
        """Error should inherit from Exception."""
        assert issubclass(Error, Exception)

    @pytest.mark.unit
    def test_error_is_type(self):
        """Error should be a class (type)."""
        assert isinstance(Error, type)

    @pytest.mark.unit
    def test_error_can_be_raised(self):
        """Error should be raisable with a message."""
        with pytest.raises(Error) as exc_info:
            raise Error("test error message")
        assert str(exc_info.value) == "test error message"

    @pytest.mark.unit
    def test_error_can_be_raised_without_message(self):
        """Error should be raisable without a message."""
        with pytest.raises(Error):
            raise Error()

    @pytest.mark.unit
    def test_error_can_be_caught_as_exception(self):
        """Error should be catchable as Exception."""
        try:
            raise Error("test")
        except Exception as e:
            assert isinstance(e, Error)
            assert str(e) == "test"

    @pytest.mark.unit
    def test_error_exported_from_libephemeris(self):
        """Error should be accessible via libephemeris.Error."""
        assert hasattr(ephem, "Error")
        assert ephem.Error is Error

    @pytest.mark.unit
    def test_error_in_all(self):
        """Error should be in __all__ of libephemeris."""
        assert "Error" in ephem.__all__


class TestErrorCompatibilityWithPySwisseph:
    """Test that Error has the same structure as swe.Error."""

    @pytest.mark.unit
    def test_same_base_class(self):
        """Both Error classes should inherit from Exception."""
        assert issubclass(Error, Exception)
        assert issubclass(swe.Error, Exception)

    @pytest.mark.unit
    def test_same_type(self):
        """Both Error classes should be types."""
        assert isinstance(Error, type)
        assert isinstance(swe.Error, type)

    @pytest.mark.unit
    def test_error_instances_are_exception(self):
        """Instances of both Error classes should be Exception instances."""
        ephem_error = Error("test")
        swe_error = swe.Error("test")
        assert isinstance(ephem_error, Exception)
        assert isinstance(swe_error, Exception)

    @pytest.mark.unit
    def test_error_message_handling(self):
        """Both Error classes should handle messages the same way."""
        ephem_error = Error("ephemeris not found")
        swe_error = swe.Error("ephemeris not found")
        assert str(ephem_error) == str(swe_error)


class TestErrorUsagePatterns:
    """Test common usage patterns for Error."""

    @pytest.mark.unit
    def test_try_except_pattern(self):
        """Standard try/except pattern should work."""
        caught = False
        try:
            raise Error("calculation failed")
        except Error:
            caught = True
        assert caught

    @pytest.mark.unit
    def test_error_with_formatted_message(self):
        """Error should work with formatted strings."""
        planet_id = 999
        with pytest.raises(Error) as exc_info:
            raise Error(f"Unknown planet ID: {planet_id}")
        assert "999" in str(exc_info.value)

    @pytest.mark.unit
    def test_error_chaining(self):
        """Error should support exception chaining."""
        try:
            try:
                raise ValueError("original error")
            except ValueError as e:
                raise Error("wrapped error") from e
        except Error as e:
            assert str(e) == "wrapped error"
            assert isinstance(e.__cause__, ValueError)
            assert str(e.__cause__) == "original error"

    @pytest.mark.unit
    def test_error_args(self):
        """Error should store args like standard exceptions."""
        err = Error("test", "arg2", 123)
        assert err.args == ("test", "arg2", 123)


class TestSPKNotFoundError:
    """Test the SPKNotFoundError exception class."""

    @pytest.mark.unit
    def test_spk_not_found_error_is_error_subclass(self):
        """SPKNotFoundError should inherit from Error."""
        from libephemeris.exceptions import SPKNotFoundError

        assert issubclass(SPKNotFoundError, Error)

    @pytest.mark.unit
    def test_spk_not_found_error_exported(self):
        """SPKNotFoundError should be accessible via libephemeris.SPKNotFoundError."""
        import libephemeris as eph
        from libephemeris.exceptions import SPKNotFoundError

        assert hasattr(eph, "SPKNotFoundError")
        assert eph.SPKNotFoundError is SPKNotFoundError

    @pytest.mark.unit
    def test_spk_not_found_error_in_all(self):
        """SPKNotFoundError should be in __all__ of libephemeris."""
        import libephemeris as eph

        assert "SPKNotFoundError" in eph.__all__

    @pytest.mark.unit
    def test_spk_not_found_error_basic(self):
        """SPKNotFoundError can be created with a message."""
        from libephemeris.exceptions import SPKNotFoundError

        err = SPKNotFoundError("test message")
        assert str(err) == "test message"

    @pytest.mark.unit
    def test_spk_not_found_error_with_attributes(self):
        """SPKNotFoundError stores filepath, body_name, and body_id."""
        from libephemeris.exceptions import SPKNotFoundError

        err = SPKNotFoundError(
            "test message",
            filepath="/path/to/file.bsp",
            body_name="Chiron",
            body_id="2060",
        )
        assert err.filepath == "/path/to/file.bsp"
        assert err.body_name == "Chiron"
        assert err.body_id == "2060"

    @pytest.mark.unit
    def test_spk_not_found_error_from_filepath(self):
        """SPKNotFoundError.from_filepath creates helpful error message."""
        from libephemeris.exceptions import SPKNotFoundError

        err = SPKNotFoundError.from_filepath(
            filepath="/path/to/chiron.bsp",
            body_name="Chiron",
            body_id="2060",
        )
        message = str(err)

        # Check that key information is in the message
        assert "/path/to/chiron.bsp" in message
        assert "download_spk" in message
        assert "download_and_register_spk" in message
        assert "set_auto_spk_download" in message
        assert "2060" in message  # body_id should be in examples
        assert "CHIRON" in message  # body_name should be in examples

    @pytest.mark.unit
    def test_spk_not_found_error_from_filepath_no_body_info(self):
        """SPKNotFoundError.from_filepath works without body info."""
        from libephemeris.exceptions import SPKNotFoundError

        err = SPKNotFoundError.from_filepath(
            filepath="/path/to/unknown.bsp",
        )
        message = str(err)

        # Should still have instructions even without body info
        assert "/path/to/unknown.bsp" in message
        assert "download_spk" in message

    @pytest.mark.unit
    def test_spk_not_found_error_repr(self):
        """SPKNotFoundError has informative repr."""
        from libephemeris.exceptions import SPKNotFoundError

        err = SPKNotFoundError(
            "test",
            filepath="/test.bsp",
            body_name="Ceres",
            body_id="1",
        )
        r = repr(err)
        assert "SPKNotFoundError" in r
        assert "/test.bsp" in r
        assert "Ceres" in r
