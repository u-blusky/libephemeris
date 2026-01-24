"""Exception classes for libephemeris.

This module provides pyswisseph-compatible exception classes.
"""


class Error(Exception):
    """Swiss Ephemeris error.

    This exception is raised for ephemeris-related errors such as:
    - Ephemeris files not found
    - Unsupported planet or body ID
    - Dates out of range for available ephemeris data
    - Fixed star not found
    - Calculation failures

    This class is designed to be compatible with swisseph.Error (swe.Error)
    from pyswisseph, allowing client code that catches swe.Error to work
    unchanged with libephemeris.Error.

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.calc_ut(jd, invalid_planet_id)
        ... except ephem.Error as e:
        ...     print(f"Ephemeris error: {e}")
    """

    pass
