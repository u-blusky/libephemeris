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
