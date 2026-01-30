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


class PolarCircleError(Error):
    """Error raised when house calculation fails at polar latitudes.

    Certain house systems (Placidus, Koch, Gauquelin) cannot be calculated
    at latitudes beyond the polar circle (~66.5°) because some ecliptic
    points never rise or set (circumpolar behavior).

    This exception provides detailed information about the polar circle
    condition to help users understand why the calculation failed and
    what alternatives are available.

    Attributes:
        latitude: The geographic latitude that triggered the error
        threshold: The polar circle threshold latitude for the given obliquity
        obliquity: The true obliquity of the ecliptic used in the calculation
        house_system: The requested house system that failed (e.g., 'P', 'K', 'G')
        message: Human-readable error message

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.swe_houses(jd, 70.0, 0.0, ord('P'))
        ... except ephem.PolarCircleError as e:
        ...     print(f"Polar error at {e.latitude}° (threshold: {e.threshold}°)")
        ...     # Use fallback house system
        ...     cusps, ascmc = ephem.swe_houses(jd, 70.0, 0.0, ord('O'))

    See Also:
        swe_houses_with_fallback: Automatically falls back to Porphyry
        get_polar_latitude_threshold: Returns the threshold for a given obliquity
    """

    def __init__(
        self,
        message: str,
        latitude: float | None = None,
        threshold: float | None = None,
        obliquity: float | None = None,
        house_system: str | None = None,
    ):
        super().__init__(message)
        self.latitude = latitude
        self.threshold = threshold
        self.obliquity = obliquity
        self.house_system = house_system
        self.message = message

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"PolarCircleError({self.message!r}, latitude={self.latitude}, "
            f"threshold={self.threshold}, obliquity={self.obliquity}, "
            f"house_system={self.house_system!r})"
        )


class EphemerisRangeError(Error):
    """Error raised when a calculation date is outside the ephemeris range.

    This exception is raised when the requested Julian Day falls outside the
    date range covered by the loaded ephemeris file (e.g., DE440 covers
    1550-2650).

    This exception provides detailed information about the supported date range
    to help users understand why the calculation failed and what dates are valid.

    Attributes:
        requested_jd: The Julian Day that was requested
        start_jd: The start of the supported Julian Day range
        end_jd: The end of the supported Julian Day range
        start_date: The start date in calendar format (YYYY-MM-DD)
        end_date: The end date in calendar format (YYYY-MM-DD)
        body_id: The body ID that was being calculated (if available)
        body_name: Human-readable name of the body (if available)
        ephemeris_file: The ephemeris file in use (if available)
        message: Human-readable error message

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.calc_ut(5000000.0, ephem.SE_SUN, 0)  # Year ~8827 AD
        ... except ephem.EphemerisRangeError as e:
        ...     print(f"Date JD {e.requested_jd} is outside range")
        ...     print(f"Supported: JD {e.start_jd} to {e.end_jd}")
        ...     print(f"           ({e.start_date} to {e.end_date})")

    See Also:
        get_current_file_data: Returns information about the loaded ephemeris
        set_ephemeris_file: Load a different ephemeris with different date range
    """

    def __init__(
        self,
        message: str,
        requested_jd: float | None = None,
        start_jd: float | None = None,
        end_jd: float | None = None,
        start_date: str | None = None,
        end_date: str | None = None,
        body_id: int | None = None,
        body_name: str | None = None,
        ephemeris_file: str | None = None,
    ):
        super().__init__(message)
        self.requested_jd = requested_jd
        self.start_jd = start_jd
        self.end_jd = end_jd
        self.start_date = start_date
        self.end_date = end_date
        self.body_id = body_id
        self.body_name = body_name
        self.ephemeris_file = ephemeris_file
        self.message = message

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"EphemerisRangeError({self.message!r}, "
            f"requested_jd={self.requested_jd}, "
            f"start_jd={self.start_jd}, end_jd={self.end_jd}, "
            f"start_date={self.start_date!r}, end_date={self.end_date!r}, "
            f"body_id={self.body_id}, body_name={self.body_name!r}, "
            f"ephemeris_file={self.ephemeris_file!r})"
        )
