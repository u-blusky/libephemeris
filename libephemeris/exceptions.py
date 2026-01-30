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


class SPKNotFoundError(Error):
    """Error raised when an SPK file is requested but not available.

    This exception is raised when attempting to register or load an SPK file
    that does not exist. It provides helpful instructions for obtaining the
    required SPK file.

    SPK (SPICE kernel) files contain high-precision ephemeris data for minor
    bodies (asteroids, TNOs, comets) and are downloaded from JPL Horizons.

    Attributes:
        filepath: Path to the SPK file that was not found
        body_name: Name of the celestial body (if known)
        body_id: JPL Horizons identifier (if known)
        message: Human-readable error message with instructions

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.register_spk_body(ephem.SE_CHIRON, "/path/to/chiron.bsp", 2002060)
        ... except ephem.SPKNotFoundError as e:
        ...     print(f"Missing SPK: {e.filepath}")
        ...     print("Instructions:", e.message)

    See Also:
        download_spk: Download SPK file from JPL Horizons
        download_and_register_spk: Download and register in one step
        enable_auto_spk: Enable automatic SPK downloading
    """

    def __init__(
        self,
        message: str,
        filepath: str | None = None,
        body_name: str | None = None,
        body_id: str | None = None,
    ):
        super().__init__(message)
        self.filepath = filepath
        self.body_name = body_name
        self.body_id = body_id
        self.message = message

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"SPKNotFoundError({self.message!r}, "
            f"filepath={self.filepath!r}, body_name={self.body_name!r}, "
            f"body_id={self.body_id!r})"
        )

    @classmethod
    def from_filepath(
        cls,
        filepath: str,
        body_name: str | None = None,
        body_id: str | None = None,
    ) -> "SPKNotFoundError":
        """Create an SPKNotFoundError with a helpful message.

        Args:
            filepath: Path to the missing SPK file
            body_name: Optional name of the celestial body
            body_id: Optional JPL Horizons identifier

        Returns:
            SPKNotFoundError with formatted instructions
        """
        lines = [f"SPK file not found: {filepath}"]
        lines.append("")
        lines.append("To obtain this SPK file, you can:")
        lines.append("")
        lines.append("1. Download manually using libephemeris.download_spk():")
        if body_id:
            lines.append("   >>> import libephemeris as eph")
            lines.append(f'   >>> eph.download_spk("{body_id}")')
        else:
            lines.append("   >>> import libephemeris as eph")
            lines.append('   >>> eph.download_spk("<body_id_or_name>")')
        lines.append("")
        lines.append("2. Download and register in one step:")
        if body_id and body_name:
            lines.append(
                f'   >>> eph.download_and_register_spk("{body_id}", eph.SE_{body_name.upper()})'
            )
        else:
            lines.append(
                '   >>> eph.download_and_register_spk("<body_id>", <libephemeris_body_id>)'
            )
        lines.append("")
        lines.append("3. Enable automatic downloading (requires astroquery):")
        lines.append("   >>> eph.set_auto_spk_download(True)")
        lines.append("")
        lines.append("4. Use the CLI download script:")
        lines.append(
            "   $ python -m libephemeris.scripts.download_spk --bodies <body_name>"
        )
        lines.append("")
        lines.append("For more information, see the libephemeris documentation.")

        message = "\n".join(lines)
        return cls(
            message=message,
            filepath=filepath,
            body_name=body_name,
            body_id=body_id,
        )


class CoordinateError(Error):
    """Error raised when geographic coordinates are invalid.

    This exception is raised when latitude or longitude values are outside
    their valid ranges:
    - Latitude: must be in [-90, 90] degrees
    - Longitude: must be in [-180, 180] degrees

    Attributes:
        message: Human-readable error message
        coordinate_name: Name of the coordinate ("latitude" or "longitude")
        value: The invalid coordinate value
        min_value: Minimum allowed value for this coordinate
        max_value: Maximum allowed value for this coordinate

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.set_topo(0.0, 91.0, 0.0)  # Invalid latitude
        ... except ephem.CoordinateError as e:
        ...     print(f"Invalid {e.coordinate_name}: {e.value}")
        ...     print(f"Valid range: [{e.min_value}, {e.max_value}]")
        Invalid latitude: 91.0
        Valid range: [-90, 90]

    See Also:
        set_topo: Sets observer location (validates coordinates)
        swe_houses: House calculation (validates latitude)
    """

    def __init__(
        self,
        message: str,
        coordinate_name: str | None = None,
        value: float | None = None,
        min_value: float | None = None,
        max_value: float | None = None,
    ):
        super().__init__(message)
        self.message = message
        self.coordinate_name = coordinate_name
        self.value = value
        self.min_value = min_value
        self.max_value = max_value

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"CoordinateError({self.message!r}, "
            f"coordinate_name={self.coordinate_name!r}, value={self.value}, "
            f"min_value={self.min_value}, max_value={self.max_value})"
        )


def validate_latitude(lat: float, func_name: str = "") -> None:
    """Validate that latitude is within valid range [-90, 90].

    Args:
        lat: Geographic latitude in degrees
        func_name: Optional function name for error message context

    Raises:
        CoordinateError: If latitude is outside [-90, 90] range

    Example:
        >>> validate_latitude(45.0)  # OK
        >>> validate_latitude(91.0)  # Raises CoordinateError
    """
    if lat < -90.0 or lat > 90.0:
        prefix = f"{func_name}: " if func_name else ""
        message = (
            f"{prefix}latitude {lat} is out of valid range. "
            f"Latitude must be between -90 and 90 degrees."
        )
        raise CoordinateError(
            message=message,
            coordinate_name="latitude",
            value=lat,
            min_value=-90.0,
            max_value=90.0,
        )


def validate_longitude(lon: float, func_name: str = "") -> None:
    """Validate that longitude is within valid range [-180, 180].

    Args:
        lon: Geographic longitude in degrees
        func_name: Optional function name for error message context

    Raises:
        CoordinateError: If longitude is outside [-180, 180] range

    Example:
        >>> validate_longitude(12.5)  # OK
        >>> validate_longitude(200.0)  # Raises CoordinateError
    """
    if lon < -180.0 or lon > 180.0:
        prefix = f"{func_name}: " if func_name else ""
        message = (
            f"{prefix}longitude {lon} is out of valid range. "
            f"Longitude must be between -180 and 180 degrees."
        )
        raise CoordinateError(
            message=message,
            coordinate_name="longitude",
            value=lon,
            min_value=-180.0,
            max_value=180.0,
        )


def validate_coordinates(lat: float, lon: float, func_name: str = "") -> None:
    """Validate both latitude and longitude coordinates.

    Args:
        lat: Geographic latitude in degrees
        lon: Geographic longitude in degrees
        func_name: Optional function name for error message context

    Raises:
        CoordinateError: If latitude is outside [-90, 90] or
                        longitude is outside [-180, 180]

    Example:
        >>> validate_coordinates(41.9, 12.5)  # OK (Rome)
        >>> validate_coordinates(91.0, 0.0)  # Raises CoordinateError
    """
    validate_latitude(lat, func_name)
    validate_longitude(lon, func_name)


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
