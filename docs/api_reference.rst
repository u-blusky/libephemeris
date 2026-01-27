libephemeris API Reference
==========================

.. module:: libephemeris
   :synopsis: A pure-Python Swiss Ephemeris-compatible astronomical library.

This is the complete API reference for libephemeris, a pure-Python implementation
of astronomical ephemeris calculations compatible with the Swiss Ephemeris.

.. contents:: Table of Contents
   :local:
   :depth: 2


Exceptions
----------

.. autoexception:: Error

   Swiss Ephemeris-compatible exception class.

   This exception is raised for ephemeris-related errors such as:

   - Ephemeris files not found
   - Unsupported planet or body ID
   - Dates out of range for available ephemeris data
   - Fixed star not found
   - Calculation failures

   Example:
       >>> import libephemeris as ephem
       >>> try:
       ...     ephem.calc_ut(jd, invalid_planet_id)
       ... except ephem.Error as e:
       ...     print(f"Ephemeris error: {e}")


Thread-Safe Context API
-----------------------

EphemerisContext
~~~~~~~~~~~~~~~~

.. autoclass:: EphemerisContext

   Thread-safe context for ephemeris calculations.

   Each instance maintains independent calculation state while sharing
   expensive resources (ephemeris files, timescale) across instances.

   **Attributes:**

   - ``topo`` (Optional[Topos]): Observer topocentric location (or None for geocentric)
   - ``sidereal_mode`` (int): Active sidereal mode ID (1 = Lahiri by default)
   - ``sidereal_t0`` (float): Reference epoch for custom ayanamsha (JD)
   - ``sidereal_ayan_t0`` (float): Ayanamsha value at reference epoch (degrees)

   **Example:**

   >>> ctx = EphemerisContext()
   >>> ctx.set_topo(12.5, 41.9, 0)  # Rome coordinates
   >>> ctx.set_sid_mode(1)  # Lahiri ayanamsha
   >>> pos, flags = ctx.calc_ut(2451545.0, SE_MARS, SEFLG_SIDEREAL)
   >>> print(f"Mars at {pos[0]:.2f} sidereal longitude")

   **Thread Safety:**

   Each context instance is independent. Multiple threads can create
   and use their own contexts without interference.

   **Methods:**

   .. method:: __init__(ephe_path=None, ephe_file="de421.bsp")

      Initialize a new ephemeris context.

      :param ephe_path: Optional path to directory containing ephemeris files.
                        If None, uses default workspace directory.
      :type ephe_path: str, optional
      :param ephe_file: Ephemeris file to use (default: "de421.bsp")
      :type ephe_file: str

   .. method:: set_topo(lon, lat, alt)

      Set observer's topocentric location for calculations.

      :param lon: Geographic longitude in degrees (East positive)
      :type lon: float
      :param lat: Geographic latitude in degrees (North positive)
      :type lat: float
      :param alt: Elevation above sea level in meters
      :type alt: float

   .. method:: get_topo()

      Get the observer's topocentric location.

      :returns: Current observer location or None if not set
      :rtype: Optional[Topos]

   .. method:: set_sid_mode(mode, t0=0.0, ayan_t0=0.0)

      Set the sidereal mode (ayanamsha system) for calculations.

      :param mode: Sidereal mode ID (SE_SIDM_*) or 255 for custom
      :type mode: int
      :param t0: Reference epoch (Julian Day) for custom ayanamsha (default: J2000.0)
      :type t0: float
      :param ayan_t0: Ayanamsha value at t0 in degrees (for custom mode)
      :type ayan_t0: float

   .. method:: get_sid_mode(full=False)

      Get the current sidereal mode configuration.

      :param full: If True, return (mode, t0, ayan_t0); if False, return only mode ID
      :type full: bool
      :returns: Sidereal mode ID, or full configuration tuple
      :rtype: int or tuple

   .. method:: calc_ut(tjd_ut, ipl, iflag)

      Calculate planetary position for Universal Time.

      :param tjd_ut: Julian Day in Universal Time (UT1)
      :type tjd_ut: float
      :param ipl: Planet/body ID (SE_SUN, SE_MOON, etc.)
      :type ipl: int
      :param iflag: Calculation flags (SEFLG_SPEED, SEFLG_HELCTR, etc.)
      :type iflag: int
      :returns: Tuple of (position_tuple, return_flag) where position_tuple is
                (longitude, latitude, distance, speed_lon, speed_lat, speed_dist)
      :rtype: tuple

   .. method:: calc(tjd, ipl, iflag)

      Calculate planetary position for Ephemeris Time (ET/TT).

      :param tjd: Julian Day in Terrestrial Time (TT/ET)
      :type tjd: float
      :param ipl: Planet/body ID (SE_SUN, SE_MOON, etc.)
      :type ipl: int
      :param iflag: Calculation flags (SEFLG_SPEED, SEFLG_HELCTR, etc.)
      :type iflag: int
      :returns: Tuple of (position_tuple, return_flag)
      :rtype: tuple

   .. method:: houses(tjd_ut, lat, lon, hsys)

      Calculate house cusps and angles.

      :param tjd_ut: Julian Day in Universal Time
      :type tjd_ut: float
      :param lat: Geographic latitude in degrees
      :type lat: float
      :param lon: Geographic longitude in degrees
      :type lon: float
      :param hsys: House system character code (ord('P') for Placidus, etc.)
      :type hsys: int
      :returns: Tuple of (cusps, ascmc) where cusps is list of 12 house cusp
                longitudes and ascmc is list of angles [ASC, MC, ARMC, Vertex, ...]
      :rtype: tuple

   .. classmethod:: close()

      Close all shared ephemeris resources and release file handles.

      This class method closes the shared SPK kernel file handles and resets
      all shared resources. Call this when you want to:

      - Free memory and file handles in long-running applications
      - Switch to a different ephemeris file
      - Ensure clean state in test suites


Time Functions
--------------

swe_julday / julday
~~~~~~~~~~~~~~~~~~~

.. function:: swe_julday(year, month, day, hour, gregflag=SE_GREG_CAL)

   Convert calendar date to Julian Day number.

   :param year: Calendar year (negative for BCE)
   :type year: int
   :param month: Month (1-12)
   :type month: int
   :param day: Day of month (1-31)
   :type day: int
   :param hour: Decimal hour (0.0-23.999...)
   :type hour: float
   :param gregflag: SE_GREG_CAL (1) for Gregorian, SE_JUL_CAL (0) for Julian
   :type gregflag: int
   :returns: Julian Day number (days since JD 0.0 = noon Jan 1, 4713 BCE)
   :rtype: float

   **Example:**

   >>> jd = swe_julday(2000, 1, 1, 12.0)  # J2000.0 epoch
   >>> print(jd)
   2451545.0

   **Note:**

   Transition date: Oct 15, 1582 (Gregorian) = Oct 5, 1582 (Julian)
   JD 2451545.0 = Jan 1, 2000 12:00 TT (J2000.0 epoch)


swe_revjul / revjul
~~~~~~~~~~~~~~~~~~~

.. function:: swe_revjul(jd, gregflag=SE_GREG_CAL)

   Convert Julian Day number to calendar date.

   :param jd: Julian Day number
   :type jd: float
   :param gregflag: SE_GREG_CAL (1) for Gregorian, SE_JUL_CAL (0) for Julian
   :type gregflag: int
   :returns: Tuple of (year, month, day, hour)
   :rtype: tuple[int, int, int, float]

   **Example:**

   >>> year, month, day, hour = swe_revjul(2451545.0)
   >>> print(f"{year}-{month:02d}-{day:02d} {hour:.1f}h")
   2000-01-01 12.0h


swe_deltat / deltat
~~~~~~~~~~~~~~~~~~~

.. function:: swe_deltat(tjd)

   Calculate Delta T (TT - UT1) for a given Julian Day.

   :param tjd: Julian Day number in UT1
   :type tjd: float
   :returns: Delta T in days (TT - UT1)
   :rtype: float

   **Note:**

   Delta T accounts for Earth's irregular rotation and is required
   for high-precision planetary calculations. Values are obtained
   from IERS (International Earth Rotation Service) data.

   - For modern dates: ~0.0008 days (~69 seconds as of 2024)
   - For historical dates: Calculated from polynomial models


swe_deltat_ex / deltat_ex
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_deltat_ex(tjd, ephe_flag=SEFLG_SWIEPH)

   Calculate Delta T with explicit ephemeris source specification.

   :param tjd: Julian Day number in UT1
   :type tjd: float
   :param ephe_flag: Ephemeris selection flag (SEFLG_SWIEPH, SEFLG_JPLEPH, SEFLG_MOSEPH)
   :type ephe_flag: int
   :returns: Tuple of (delta_t, error_message)
   :rtype: tuple[float, str]

   **Example:**

   >>> dt, err = swe_deltat_ex(2451545.0, SEFLG_SWIEPH)
   >>> print(f"Delta T: {dt * 86400:.2f} seconds")
   Delta T: 63.83 seconds


date_conversion
~~~~~~~~~~~~~~~

.. function:: date_conversion(year, month, day, hour, calendar)

   Convert a date between Julian and Gregorian calendars.

   :param year: Calendar year
   :type year: int
   :param month: Month (1-12)
   :type month: int
   :param day: Day of month (1-31)
   :type day: int
   :param hour: Decimal hour (0.0-23.999...)
   :type hour: float
   :param calendar: Target calendar - 'j' for Julian or 'g' for Gregorian
   :type calendar: str
   :returns: Tuple of (year, month, day, hour) in the requested calendar
   :rtype: tuple[int, int, int, float]
   :raises ValueError: If calendar is not 'j' or 'g'

   **Example:**

   >>> # Convert first Gregorian date to Julian
   >>> date_conversion(1582, 10, 15, 12.0, 'j')
   (1582, 10, 5, 12.0)


day_of_week
~~~~~~~~~~~

.. function:: day_of_week(jd)

   Calculate day of week for a Julian Day.

   :param jd: Julian Day number
   :type jd: float
   :returns: Day of week (0=Monday, 6=Sunday)
   :rtype: int


utc_to_jd
~~~~~~~~~

.. function:: utc_to_jd(year, month, day, hour, minute, second, gregflag=SE_GREG_CAL)

   Convert UTC date/time to Julian Day numbers.

   :param year: Year
   :type year: int
   :param month: Month (1-12)
   :type month: int
   :param day: Day of month
   :type day: int
   :param hour: Hour (0-23)
   :type hour: int
   :param minute: Minute (0-59)
   :type minute: int
   :param second: Second (0-59, can include fractional part)
   :type second: float
   :param gregflag: Calendar flag (SE_GREG_CAL or SE_JUL_CAL)
   :type gregflag: int
   :returns: Tuple of (jd_et, jd_ut) - Julian Days in ET and UT
   :rtype: tuple[float, float]


sidtime / sidtime0
~~~~~~~~~~~~~~~~~~

.. function:: sidtime(jd_ut)

   Calculate local apparent sidereal time.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :returns: Sidereal time in hours (0-24)
   :rtype: float

.. function:: sidtime0(jd_ut, eps, nut)

   Calculate sidereal time with given nutation values.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param eps: Mean obliquity in degrees
   :type eps: float
   :param nut: Nutation in longitude in degrees
   :type nut: float
   :returns: Sidereal time in hours (0-24)
   :rtype: float


time_equ
~~~~~~~~

.. function:: time_equ(jd_ut)

   Calculate equation of time.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :returns: Equation of time in decimal days
   :rtype: float


Planet Calculation Functions
----------------------------

swe_calc_ut / calc_ut
~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_calc_ut(tjd_ut, ipl, iflag)

   Calculate planetary position for Universal Time.

   This is the primary function for computing planetary positions.

   :param tjd_ut: Julian Day in Universal Time (UT1)
   :type tjd_ut: float
   :param ipl: Planet/body ID (SE_SUN=0, SE_MOON=1, SE_MERCURY=2, etc.)
   :type ipl: int
   :param iflag: Calculation flags (bitwise OR of SEFLG_* constants)
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple[tuple[float, float, float, float, float, float], int]

   **Position tuple contents:**

   - ``[0]`` longitude: Ecliptic longitude in degrees (0-360)
   - ``[1]`` latitude: Ecliptic latitude in degrees
   - ``[2]`` distance: Distance in AU
   - ``[3]`` speed_lon: Daily motion in longitude (degrees/day)
   - ``[4]`` speed_lat: Daily motion in latitude (degrees/day)
   - ``[5]`` speed_dist: Daily motion in distance (AU/day)

   **Calculation Flags:**

   - ``SEFLG_SPEED``: Include velocity (always calculated)
   - ``SEFLG_HELCTR``: Heliocentric instead of geocentric
   - ``SEFLG_TOPOCTR``: Topocentric (requires swe_set_topo)
   - ``SEFLG_SIDEREAL``: Sidereal zodiac (requires swe_set_sid_mode)
   - ``SEFLG_EQUATORIAL``: Return RA/Dec instead of lon/lat
   - ``SEFLG_J2000``: J2000.0 reference frame
   - ``SEFLG_TRUEPOS``: True geometric position (no light time)

   **Example:**

   >>> pos, retflag = swe_calc_ut(2451545.0, SE_MARS, SEFLG_SPEED)
   >>> lon, lat, dist = pos[0], pos[1], pos[2]
   >>> print(f"Mars: {lon:.4f} lon, {lat:.4f} lat, {dist:.6f} AU")


swe_calc / calc
~~~~~~~~~~~~~~~

.. function:: swe_calc(tjd, ipl, iflag)

   Calculate planetary position for Ephemeris Time (ET/TT).

   Similar to swe_calc_ut() but takes Terrestrial Time instead of Universal Time.

   :param tjd: Julian Day in Terrestrial Time (TT/ET)
   :type tjd: float
   :param ipl: Planet/body ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple

   **Note:**

   TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
   which varies from ~32 seconds (year 2000) to minutes (historical times).
   For most astrological applications, use swe_calc_ut() instead.


swe_calc_pctr / calc_pctr
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_calc_pctr(tjd_ut, ipl, iplctr, iflag)

   Calculate planetary position as seen from another planet (planet-centric).

   :param tjd_ut: Julian Day in Universal Time (UT1)
   :type tjd_ut: float
   :param ipl: Target planet/body ID
   :type ipl: int
   :param iplctr: Observer/center planet ID (the body from which to observe)
   :type iplctr: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple

   **Example:**

   >>> # Position of Moon as seen from Mars
   >>> pos, retflag = swe_calc_pctr(2451545.0, SE_MOON, SE_MARS, SEFLG_SPEED)
   >>> print(f"Moon longitude from Mars: {pos[0]:.2f} degrees")


get_planet_name
~~~~~~~~~~~~~~~

.. function:: get_planet_name(planet_id)

   Get the human-readable name of a planet given its ID.

   :param planet_id: Planet/body ID (SE_SUN, SE_MOON, etc.)
   :type planet_id: int
   :returns: Human-readable planet name, or "Unknown (ID)" for unrecognized IDs
   :rtype: str

   **Example:**

   >>> get_planet_name(0)
   'Sun'
   >>> get_planet_name(4)
   'Mars'


Orbital Elements
----------------

swe_nod_aps / nod_aps
~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_nod_aps(tjd, ipl, iflag, method)

   Calculate nodes and apsides of a planetary orbit.

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :param method: Method flag (SE_NODBIT_MEAN, SE_NODBIT_OSCU, etc.)
   :type method: int
   :returns: Tuple of four position tuples (ascending_node, descending_node,
             perihelion, aphelion)
   :rtype: tuple


swe_get_orbital_elements / get_orbital_elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_get_orbital_elements(tjd, ipl, iflag)

   Get Keplerian orbital elements for a planet.

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of orbital elements and flags
   :rtype: tuple


Planetary Phenomena
-------------------

swe_pheno / pheno
~~~~~~~~~~~~~~~~~

.. function:: swe_pheno(tjd, ipl, iflag)

   Calculate planetary phenomena (phase angle, elongation, etc.).

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of phenomena data
   :rtype: tuple


swe_pheno_ut / pheno_ut
~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_pheno_ut(tjd_ut, ipl, iflag)

   Calculate planetary phenomena for Universal Time.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of phenomena data
   :rtype: tuple


House Functions
---------------

swe_houses / houses
~~~~~~~~~~~~~~~~~~~

.. function:: swe_houses(tjd_ut, lat, lon, hsys)

   Calculate house cusps and angles (ASCMC).

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param lat: Geographic latitude in degrees (North positive)
   :type lat: float
   :param lon: Geographic longitude in degrees (East positive)
   :type lon: float
   :param hsys: House system as ASCII code (ord('P') for Placidus, etc.)
   :type hsys: int
   :returns: Tuple of (cusps, ascmc) where cusps is a list of 12 longitudes
             and ascmc contains angles [ASC, MC, ARMC, Vertex, Equasc, Co-asc Koch,
             Co-asc Munkasey, Polar Asc]
   :rtype: tuple[list[float], list[float]]

   **Supported House Systems:**

   - ``'P'`` - Placidus (most common, time-based)
   - ``'K'`` - Koch (birthplace system)
   - ``'O'`` - Porphyrius (space-based trisection)
   - ``'R'`` - Regiomontanus (medieval)
   - ``'C'`` - Campanus (prime vertical)
   - ``'A'`` or ``'E'`` - Equal (30 degree divisions from ASC)
   - ``'W'`` - Whole Sign
   - ``'X'`` - Meridian (equatorial)
   - ``'H'`` - Azimuthal/Horizontal
   - ``'T'`` - Polich-Page (Topocentric)
   - ``'B'`` - Alcabitus
   - ``'M'`` - Morinus
   - ``'U'`` - Krusinski-Pisa
   - ``'G'`` - Gauquelin sectors
   - ``'V'`` - Vehlow equal

   **Example:**

   >>> cusps, ascmc = swe_houses(2451545.0, 41.9, 12.5, ord('P'))
   >>> print(f"Ascendant: {ascmc[0]:.2f} degrees")
   >>> print(f"House 1 cusp: {cusps[0]:.2f} degrees")


swe_houses_ex / houses_ex
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_houses_ex(tjd_ut, iflag, lat, lon, hsys)

   Extended house calculation with sidereal support.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param iflag: Calculation flags (SEFLG_SIDEREAL for sidereal)
   :type iflag: int
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param hsys: House system code
   :type hsys: int
   :returns: Tuple of (cusps, ascmc)
   :rtype: tuple


swe_houses_ex2 / houses_ex2
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_houses_ex2(tjd_ut, iflag, lat, lon, hsys)

   Extended house calculation with velocities.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param iflag: Calculation flags
   :type iflag: int
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param hsys: House system code
   :type hsys: int
   :returns: Tuple of (cusps, ascmc, cusps_speed, ascmc_speed)
   :rtype: tuple


swe_house_pos / house_pos
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_house_pos(armc, geolat, eps, hsys, lon, lat)

   Calculate house position of a celestial point.

   :param armc: ARMC (sidereal time in degrees)
   :type armc: float
   :param geolat: Geographic latitude in degrees
   :type geolat: float
   :param eps: Obliquity of the ecliptic in degrees
   :type eps: float
   :param hsys: House system code
   :type hsys: int
   :param lon: Ecliptic longitude of the point
   :type lon: float
   :param lat: Ecliptic latitude of the point
   :type lat: float
   :returns: House position (1.0-12.999...)
   :rtype: float


swe_house_name / house_name
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_house_name(hsys)

   Get the name of a house system.

   :param hsys: House system code
   :type hsys: int
   :returns: Name of the house system
   :rtype: str


gauquelin_sector / swe_gauquelin_sector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: gauquelin_sector(tjd_ut, ipl, iflag, imeth, geolon, geolat, geoalt, atpress, attemp)

   Calculate Gauquelin sector position of a planet.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :param imeth: Method (0=sector, 1=sector/house position)
   :type imeth: int
   :param geolon: Geographic longitude
   :type geolon: float
   :param geolat: Geographic latitude
   :type geolat: float
   :param geoalt: Altitude in meters
   :type geoalt: float
   :param atpress: Atmospheric pressure in mbar
   :type atpress: float
   :param attemp: Temperature in Celsius
   :type attemp: float
   :returns: Gauquelin sector (1-36)
   :rtype: float


Ayanamsha (Sidereal) Functions
------------------------------

swe_set_sid_mode / set_sid_mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_set_sid_mode(sid_mode, t0=0.0, ayan_t0=0.0)

   Set the sidereal mode (ayanamsha system).

   :param sid_mode: Sidereal mode ID (SE_SIDM_* constants)
   :type sid_mode: int
   :param t0: Reference epoch for custom ayanamsha (Julian Day)
   :type t0: float
   :param ayan_t0: Ayanamsha value at t0 in degrees
   :type ayan_t0: float

   **Common Ayanamsha Systems:**

   - ``SE_SIDM_FAGAN_BRADLEY`` (0): Fagan-Bradley (Western sidereal)
   - ``SE_SIDM_LAHIRI`` (1): Lahiri (Indian standard)
   - ``SE_SIDM_RAMAN`` (3): B.V. Raman
   - ``SE_SIDM_KRISHNAMURTI`` (5): KP system
   - ``SE_SIDM_TRUE_CITRA`` (27): True position of Spica
   - ``SE_SIDM_USER`` (255): Custom user-defined

   **Example:**

   >>> swe_set_sid_mode(SE_SIDM_LAHIRI)
   >>> pos, _ = swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)


swe_get_ayanamsa_ut / get_ayanamsa_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_get_ayanamsa_ut(tjd_ut)

   Get the ayanamsha value for Universal Time.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :returns: Ayanamsha value in degrees
   :rtype: float

   **Example:**

   >>> ayan = swe_get_ayanamsa_ut(2451545.0)
   >>> print(f"Ayanamsha at J2000: {ayan:.4f} degrees")


swe_get_ayanamsa / get_ayanamsa
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_get_ayanamsa(tjd)

   Get the ayanamsha value for Ephemeris Time.

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :returns: Ayanamsha value in degrees
   :rtype: float


swe_get_ayanamsa_ex / get_ayanamsa_ex
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_get_ayanamsa_ex(tjd, iflag)

   Extended ayanamsha calculation with flags.

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (ayanamsha, return_flag)
   :rtype: tuple


swe_get_ayanamsa_name / get_ayanamsa_name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_get_ayanamsa_name(sid_mode)

   Get the name of an ayanamsha system.

   :param sid_mode: Sidereal mode ID
   :type sid_mode: int
   :returns: Name of the ayanamsha system
   :rtype: str


Fixed Star Functions
--------------------

swe_fixstar_ut / fixstar_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_fixstar_ut(star, tjd_ut, iflag)

   Calculate fixed star position for Universal Time.

   :param star: Star name or designation (e.g., "Regulus", "alLeo")
   :type star: str
   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (star_name, position_tuple, return_flag)
   :rtype: tuple
   :raises Error: If star not found in catalog

   **Example:**

   >>> name, pos, flag = swe_fixstar_ut("Regulus", 2451545.0, 0)
   >>> print(f"{name}: {pos[0]:.4f} longitude")


swe_fixstar / fixstar
~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_fixstar(star, tjd, iflag)

   Calculate fixed star position for Ephemeris Time.

   :param star: Star name or designation
   :type star: str
   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (star_name, position_tuple, return_flag)
   :rtype: tuple


swe_fixstar2 / fixstar2
~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_fixstar2(star, tjd, iflag)

   Enhanced fixed star lookup with catalog information.

   :param star: Star name, Bayer designation, or Hipparcos number
   :type star: str
   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (star_name, position_tuple, return_flag)
   :rtype: tuple


swe_fixstar_mag / fixstar_mag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_fixstar_mag(star)

   Get the visual magnitude of a fixed star.

   :param star: Star name or designation
   :type star: str
   :returns: Visual magnitude (apparent brightness)
   :rtype: float


Crossing Event Functions
------------------------

swe_solcross_ut / solcross_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_solcross_ut(x2cross, jd_start, iflag)

   Find when the Sun crosses a specific ecliptic longitude.

   :param x2cross: Target longitude in degrees
   :type x2cross: float
   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Julian Day of crossing
   :rtype: float

   **Example:**

   >>> # Find next vernal equinox (Sun at 0 Aries)
   >>> jd_equinox = swe_solcross_ut(0.0, 2451545.0, 0)


swe_mooncross_ut / mooncross_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_mooncross_ut(x2cross, jd_start, iflag)

   Find when the Moon crosses a specific ecliptic longitude.

   :param x2cross: Target longitude in degrees
   :type x2cross: float
   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Julian Day of crossing
   :rtype: float


swe_mooncross_node_ut / mooncross_node_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_mooncross_node_ut(jd_start, iflag)

   Find when the Moon crosses its ascending node.

   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (jd_crossing, longitude, latitude)
   :rtype: tuple


swe_cross_ut
~~~~~~~~~~~~

.. function:: swe_cross_ut(ipl, x2cross, jd_start, iflag)

   Find when a planet crosses a specific ecliptic longitude.

   :param ipl: Planet ID
   :type ipl: int
   :param x2cross: Target longitude in degrees
   :type x2cross: float
   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Julian Day of crossing
   :rtype: float


swe_helio_cross_ut / helio_cross_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_helio_cross_ut(ipl, x2cross, jd_start, iflag)

   Find when a planet crosses a heliocentric longitude.

   :param ipl: Planet ID
   :type ipl: int
   :param x2cross: Target heliocentric longitude in degrees
   :type x2cross: float
   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Julian Day of crossing
   :rtype: float


Eclipse Functions
-----------------

sol_eclipse_when_glob / swe_sol_eclipse_when_glob
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: sol_eclipse_when_glob(jd_start, iflag, ifltype)

   Find next global solar eclipse.

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :param ifltype: Eclipse type filter (SE_ECL_TOTAL, SE_ECL_ANNULAR, etc.)
   :type ifltype: int
   :returns: Tuple of (return_flag, time_array) with eclipse type and times
   :rtype: tuple

   **Time array indices:**

   - [0]: Time of maximum eclipse
   - [1]: Time of first contact (partial phase begins)
   - [2]: Time of second contact (total/annular phase begins)
   - [3]: Time of third contact (total/annular phase ends)
   - [4]: Time of fourth contact (partial phase ends)


sol_eclipse_when_loc / swe_sol_eclipse_when_loc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: sol_eclipse_when_loc(jd_start, iflag, lon, lat, alt)

   Find next solar eclipse visible at a specific location.

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param alt: Altitude in meters
   :type alt: float
   :returns: Tuple of (return_flag, time_array, attr_array)
   :rtype: tuple


sol_eclipse_where / swe_sol_eclipse_where
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: sol_eclipse_where(jd, iflag)

   Calculate geographic coordinates of eclipse maximum.

   :param jd: Julian Day of eclipse
   :type jd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (return_flag, geo_array, attr_array)
   :rtype: tuple


sol_eclipse_how / swe_sol_eclipse_how
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: sol_eclipse_how(jd, iflag, lon, lat, alt)

   Calculate eclipse circumstances at a specific location.

   :param jd: Julian Day
   :type jd: float
   :param iflag: Calculation flags
   :type iflag: int
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param alt: Altitude in meters
   :type alt: float
   :returns: Tuple of (return_flag, attr_array)
   :rtype: tuple


lun_eclipse_when / swe_lun_eclipse_when
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: lun_eclipse_when(jd_start, iflag, ifltype)

   Find next lunar eclipse.

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :param ifltype: Eclipse type filter (SE_ECL_TOTAL, SE_ECL_PARTIAL, SE_ECL_PENUMBRAL)
   :type ifltype: int
   :returns: Tuple of (return_flag, time_array)
   :rtype: tuple


Rise/Set/Transit Functions
--------------------------

rise_trans / swe_rise_trans
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: rise_trans(jd_ut, ipl, starname, rsmi, lon, lat, alt, atpress, attemp)

   Calculate rise, set, or transit times.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param ipl: Planet ID (or -1 for fixed star)
   :type ipl: int
   :param starname: Star name (if ipl=-1) or empty string
   :type starname: str
   :param rsmi: Event type flags (SE_CALC_RISE, SE_CALC_SET, SE_CALC_MTRANSIT)
   :type rsmi: int
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param alt: Altitude in meters
   :type alt: float
   :param atpress: Atmospheric pressure in mbar (0 for no refraction)
   :type atpress: float
   :param attemp: Temperature in Celsius
   :type attemp: float
   :returns: Tuple of (return_flag, jd_event)
   :rtype: tuple


rise_trans_true_hor / swe_rise_trans_true_hor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: rise_trans_true_hor(jd_ut, ipl, starname, rsmi, lon, lat, alt, atpress, attemp, horone)

   Calculate rise/set/transit with custom horizon altitude.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param starname: Star name or empty string
   :type starname: str
   :param rsmi: Event type flags
   :type rsmi: int
   :param lon: Geographic longitude
   :type lon: float
   :param lat: Geographic latitude
   :type lat: float
   :param alt: Altitude in meters
   :type alt: float
   :param atpress: Atmospheric pressure in mbar
   :type atpress: float
   :param attemp: Temperature in Celsius
   :type attemp: float
   :param horone: Custom horizon altitude in degrees
   :type horone: float
   :returns: Tuple of (return_flag, jd_event)
   :rtype: tuple


Heliacal Events
---------------

heliacal_ut / swe_heliacal_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: heliacal_ut(jd_start, geolon, geolat, geoalt, atm, obs, objname, event_type, helflag)

   Calculate heliacal rising/setting events.

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param geolon: Geographic longitude
   :type geolon: float
   :param geolat: Geographic latitude
   :type geolat: float
   :param geoalt: Altitude in meters
   :type geoalt: float
   :param atm: Atmospheric parameters tuple
   :type atm: tuple
   :param obs: Observer parameters tuple
   :type obs: tuple
   :param objname: Object name (planet name or star)
   :type objname: str
   :param event_type: Event type (SE_HELIACAL_RISING, SE_HELIACAL_SETTING, etc.)
   :type event_type: int
   :param helflag: Heliacal calculation flags
   :type helflag: int
   :returns: Tuple of result data
   :rtype: tuple


vis_limit_mag / swe_vis_limit_mag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: vis_limit_mag(jd_ut, geolon, geolat, geoalt, atm, obs, objname, helflag)

   Calculate visibility limit magnitude.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param geolon: Geographic longitude
   :type geolon: float
   :param geolat: Geographic latitude
   :type geolat: float
   :param geoalt: Altitude in meters
   :type geoalt: float
   :param atm: Atmospheric parameters tuple
   :type atm: tuple
   :param obs: Observer parameters tuple
   :type obs: tuple
   :param objname: Object name
   :type objname: str
   :param helflag: Heliacal calculation flags
   :type helflag: int
   :returns: Visibility limit magnitude
   :rtype: float


State Management Functions
--------------------------

swe_set_topo / set_topo
~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_set_topo(lon, lat, alt)

   Set observer's topocentric location.

   :param lon: Geographic longitude in degrees (East positive)
   :type lon: float
   :param lat: Geographic latitude in degrees (North positive)
   :type lat: float
   :param alt: Elevation above sea level in meters
   :type alt: float

   **Note:**

   Required for topocentric calculations (SEFLG_TOPOCTR),
   angles (Ascendant, MC), and Arabic parts.


swe_set_ephe_path / set_ephe_path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_set_ephe_path(path)

   Set the path to ephemeris files.

   :param path: Directory containing ephemeris files
   :type path: str


swe_set_jpl_file / set_jpl_file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_set_jpl_file(filename)

   Set the JPL ephemeris file to use.

   :param filename: Ephemeris file name (e.g., "de421.bsp", "de440.bsp")
   :type filename: str


swe_set_tid_acc / set_tid_acc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_set_tid_acc(acc)

   Set tidal acceleration value for Delta T calculations.

   :param acc: Tidal acceleration in arcsec/century^2
   :type acc: float


swe_get_tid_acc / get_tid_acc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_get_tid_acc()

   Get current tidal acceleration value.

   :returns: Tidal acceleration in arcsec/century^2
   :rtype: float


swe_set_delta_t_userdef / set_delta_t_userdef
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: swe_set_delta_t_userdef(dt)

   Set a user-defined Delta T value.

   :param dt: Delta T in days (or None to clear)
   :type dt: float or None


swe_close / close
~~~~~~~~~~~~~~~~~

.. function:: swe_close()

   Close ephemeris files and release resources.

   Call this to free memory in long-running applications or
   before switching ephemeris files.


Utility Functions
-----------------

Angle Normalization
~~~~~~~~~~~~~~~~~~~

.. function:: degnorm(x)

   Normalize angle to range [0, 360).

   :param x: Angle in degrees
   :type x: float
   :returns: Normalized angle
   :rtype: float

.. function:: radnorm(x)

   Normalize angle to range [0, 2*pi).

   :param x: Angle in radians
   :type x: float
   :returns: Normalized angle
   :rtype: float


Angular Differences
~~~~~~~~~~~~~~~~~~~

.. function:: difdeg2n(a, b)

   Calculate signed angular difference (a - b) in range [-180, 180].

   :param a: First angle in degrees
   :type a: float
   :param b: Second angle in degrees
   :type b: float
   :returns: Signed difference
   :rtype: float

.. function:: difdegn(a, b)

   Calculate positive angular difference (a - b) in range [0, 360).

   :param a: First angle in degrees
   :type a: float
   :param b: Second angle in degrees
   :type b: float
   :returns: Positive difference
   :rtype: float


Midpoint Calculations
~~~~~~~~~~~~~~~~~~~~~

.. function:: deg_midp(a, b)

   Calculate midpoint of two angles (shorter arc).

   :param a: First angle in degrees
   :type a: float
   :param b: Second angle in degrees
   :type b: float
   :returns: Midpoint angle
   :rtype: float


Coordinate Transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: cotrans(coord, obliquity)

   Transform between ecliptic and equatorial coordinates.

   :param coord: Tuple of (longitude, latitude, distance)
   :type coord: tuple[float, float, float]
   :param obliquity: Obliquity in degrees (positive: ecl->eq, negative: eq->ecl)
   :type obliquity: float
   :returns: Transformed coordinates
   :rtype: tuple[float, float, float]

.. function:: cotrans_sp(coord, speed, obliquity)

   Transform coordinates and velocities.

   :param coord: Tuple of (longitude, latitude, distance)
   :type coord: tuple[float, float, float]
   :param speed: Tuple of velocities
   :type speed: tuple[float, float, float]
   :param obliquity: Obliquity in degrees
   :type obliquity: float
   :returns: Tuple of (transformed_coord, transformed_speed)
   :rtype: tuple


Horizontal Coordinates
~~~~~~~~~~~~~~~~~~~~~~

.. function:: azalt(jd, calc_flag, lat, lon, altitude, pressure, temperature, coord)

   Convert to horizontal (azimuth/altitude) coordinates.

   :param jd: Julian Day in Universal Time
   :type jd: float
   :param calc_flag: SE_ECL2HOR (ecliptic) or SE_EQU2HOR (equatorial)
   :type calc_flag: int
   :param lat: Geographic latitude
   :type lat: float
   :param lon: Geographic longitude
   :type lon: float
   :param altitude: Observer altitude in meters
   :type altitude: float
   :param pressure: Atmospheric pressure in mbar (0 for no refraction)
   :type pressure: float
   :param temperature: Temperature in Celsius
   :type temperature: float
   :param coord: Input coordinates tuple
   :type coord: tuple
   :returns: Tuple of (azimuth, true_altitude, apparent_altitude)
   :rtype: tuple

.. function:: azalt_rev(jd, calc_flag, lat, lon, altitude, azalt)

   Convert from horizontal to celestial coordinates.

   :param jd: Julian Day in Universal Time
   :type jd: float
   :param calc_flag: SE_HOR2ECL or SE_HOR2EQU
   :type calc_flag: int
   :param lat: Geographic latitude
   :type lat: float
   :param lon: Geographic longitude
   :type lon: float
   :param altitude: Observer altitude
   :type altitude: float
   :param azalt: Tuple of (azimuth, altitude)
   :type azalt: tuple
   :returns: Celestial coordinates
   :rtype: tuple


Atmospheric Refraction
~~~~~~~~~~~~~~~~~~~~~~

.. function:: refrac(true_alt, atpress, attemp, calc_flag)

   Calculate atmospheric refraction.

   :param true_alt: True altitude in degrees
   :type true_alt: float
   :param atpress: Atmospheric pressure in mbar
   :type atpress: float
   :param attemp: Temperature in Celsius
   :type attemp: float
   :param calc_flag: SE_TRUE_TO_APP or SE_APP_TO_TRUE
   :type calc_flag: int
   :returns: Apparent altitude (or true altitude if calc_flag=SE_APP_TO_TRUE)
   :rtype: float


Degree Splitting
~~~~~~~~~~~~~~~~

.. function:: split_deg(deg, roundflag)

   Split a degree value into components.

   :param deg: Degree value
   :type deg: float
   :param roundflag: Rounding flags (SPLIT_DEG_ROUND_SEC, SPLIT_DEG_ZODIACAL, etc.)
   :type roundflag: int
   :returns: Tuple of (degrees, minutes, seconds, fraction, sign)
   :rtype: tuple


Arabic Parts
------------

calc_all_arabic_parts
~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_all_arabic_parts(positions)

   Calculate all standard Arabic parts.

   :param positions: Dictionary of celestial positions with keys:
                     'Asc', 'Sun', 'Moon', 'Mercury', 'Venus'
   :type positions: dict[str, float]
   :returns: Dictionary of Arabic parts:
             'Pars_Fortunae', 'Pars_Spiritus', 'Pars_Amoris', 'Pars_Fidei'
   :rtype: dict[str, float]

   **Example:**

   >>> positions = {
   ...     'Asc': 15.5, 'Sun': 120.0, 'Moon': 240.0,
   ...     'Mercury': 130.0, 'Venus': 110.0
   ... }
   >>> parts = calc_all_arabic_parts(positions)
   >>> print(f"Part of Fortune: {parts['Pars_Fortunae']:.2f}")


Constants
---------

Planet IDs
~~~~~~~~~~

.. data:: SE_SUN
   :value: 0

.. data:: SE_MOON
   :value: 1

.. data:: SE_MERCURY
   :value: 2

.. data:: SE_VENUS
   :value: 3

.. data:: SE_MARS
   :value: 4

.. data:: SE_JUPITER
   :value: 5

.. data:: SE_SATURN
   :value: 6

.. data:: SE_URANUS
   :value: 7

.. data:: SE_NEPTUNE
   :value: 8

.. data:: SE_PLUTO
   :value: 9

.. data:: SE_MEAN_NODE
   :value: 10

   Mean lunar ascending node (Dragon's Head)

.. data:: SE_TRUE_NODE
   :value: 11

   True (osculating) lunar ascending node

.. data:: SE_MEAN_APOG
   :value: 12

   Mean lunar apogee (Black Moon Lilith)

.. data:: SE_OSCU_APOG
   :value: 13

   Osculating (true) lunar apogee

.. data:: SE_EARTH
   :value: 14

.. data:: SE_CHIRON
   :value: 15

.. data:: SE_CERES
   :value: 17

.. data:: SE_PALLAS
   :value: 18

.. data:: SE_JUNO
   :value: 19

.. data:: SE_VESTA
   :value: 20


Calculation Flags
~~~~~~~~~~~~~~~~~

.. data:: SEFLG_SPEED
   :value: 256

   Calculate velocity

.. data:: SEFLG_HELCTR
   :value: 8

   Heliocentric position

.. data:: SEFLG_TOPOCTR
   :value: 32768

   Topocentric position (requires set_topo)

.. data:: SEFLG_SIDEREAL
   :value: 65536

   Sidereal zodiac positions

.. data:: SEFLG_EQUATORIAL
   :value: 2048

   Equatorial coordinates (RA/Dec)

.. data:: SEFLG_J2000
   :value: 32

   J2000.0 reference frame

.. data:: SEFLG_TRUEPOS
   :value: 16

   True geometric position (no light time)

.. data:: SEFLG_NOABERR
   :value: 1024

   No aberration

.. data:: SEFLG_NOGDEFL
   :value: 512

   No gravitational deflection


Sidereal Modes
~~~~~~~~~~~~~~

.. data:: SE_SIDM_FAGAN_BRADLEY
   :value: 0

.. data:: SE_SIDM_LAHIRI
   :value: 1

.. data:: SE_SIDM_RAMAN
   :value: 3

.. data:: SE_SIDM_KRISHNAMURTI
   :value: 5

.. data:: SE_SIDM_TRUE_CITRA
   :value: 27

.. data:: SE_SIDM_USER
   :value: 255


Calendar Flags
~~~~~~~~~~~~~~

.. data:: SE_GREG_CAL
   :value: 1

   Gregorian calendar

.. data:: SE_JUL_CAL
   :value: 0

   Julian calendar


Eclipse Flags
~~~~~~~~~~~~~

.. data:: SE_ECL_TOTAL
   :value: 4

.. data:: SE_ECL_ANNULAR
   :value: 8

.. data:: SE_ECL_PARTIAL
   :value: 16

.. data:: SE_ECL_PENUMBRAL
   :value: 64


Rise/Set Flags
~~~~~~~~~~~~~~

.. data:: SE_CALC_RISE
   :value: 1

.. data:: SE_CALC_SET
   :value: 2

.. data:: SE_CALC_MTRANSIT
   :value: 4

   Upper culmination (meridian transit)

.. data:: SE_CALC_ITRANSIT
   :value: 8

   Lower transit (anti-culmination)


SPK Kernel Functions
--------------------

These functions provide high-precision calculations for minor bodies using
SPK (SPICE kernel) files downloaded from NASA JPL Horizons.

download_spk
~~~~~~~~~~~~

.. function:: download_spk(body, start, end, path=None)

   Download an SPK kernel file from JPL Horizons for a minor body.

   :param body: Body identifier (asteroid number or name, e.g., "2060", "Chiron")
   :type body: str
   :param start: Start date for SPK coverage (e.g., "2000-01-01")
   :type start: str
   :param end: End date for SPK coverage (e.g., "2100-01-01")
   :type end: str
   :param path: Directory to save the SPK file (default: current directory)
   :type path: str, optional
   :returns: Path to the downloaded SPK file
   :rtype: str

   **Example:**

   >>> path = download_spk("2060", "2000-01-01", "2100-01-01", "./spk")


register_spk_body
~~~~~~~~~~~~~~~~~

.. function:: register_spk_body(ipl, spk_file, naif_id)

   Register an SPK file for a minor body.

   After registration, ``calc_ut()`` will use the SPK kernel for this body
   instead of Keplerian approximations.

   :param ipl: LibEphemeris body ID (e.g., SE_CHIRON)
   :type ipl: int
   :param spk_file: Path to the SPK file
   :type spk_file: str
   :param naif_id: NAIF SPICE ID for the body
   :type naif_id: int


unregister_spk_body
~~~~~~~~~~~~~~~~~~~

.. function:: unregister_spk_body(ipl)

   Remove SPK registration for a body, reverting to Keplerian calculations.

   :param ipl: LibEphemeris body ID
   :type ipl: int


download_and_register_spk
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: download_and_register_spk(body, ipl, start, end, path=None)

   Download an SPK kernel and register it in one step.

   :param body: Body identifier for JPL Horizons
   :type body: str
   :param ipl: LibEphemeris body ID
   :type ipl: int
   :param start: Start date for SPK coverage
   :type start: str
   :param end: End date for SPK coverage
   :type end: str
   :param path: Directory for SPK file (optional)
   :type path: str, optional
   :returns: Path to the downloaded SPK file
   :rtype: str


list_spk_bodies
~~~~~~~~~~~~~~~

.. function:: list_spk_bodies()

   List all bodies with registered SPK kernels.

   :returns: List of registered body IDs
   :rtype: list[int]


get_spk_body_info
~~~~~~~~~~~~~~~~~

.. function:: get_spk_body_info(ipl)

   Get information about an SPK registration.

   :param ipl: LibEphemeris body ID
   :type ipl: int
   :returns: Dictionary with 'spk_file' and 'naif_id', or None if not registered
   :rtype: dict or None


get_spk_coverage
~~~~~~~~~~~~~~~~

.. function:: get_spk_coverage(spk_path)

   Get the date range covered by an SPK file.

   :param spk_path: Path to SPK file
   :type spk_path: str
   :returns: Tuple of (start_jd, end_jd) Julian Day numbers
   :rtype: tuple[float, float]


SPK Auto-Download Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: set_auto_spk_download(enable)

   Enable or disable automatic SPK download for minor bodies.

   When enabled, calculations for supported minor bodies will automatically
   download and cache SPK kernels from JPL Horizons.

   :param enable: True to enable, False to disable
   :type enable: bool


.. function:: get_auto_spk_download()

   Check if automatic SPK download is enabled.

   :returns: Current auto-download setting
   :rtype: bool


.. function:: set_spk_cache_dir(path)

   Set the directory for caching downloaded SPK files.

   :param path: Directory path for SPK cache
   :type path: str


.. function:: get_spk_cache_dir()

   Get the current SPK cache directory.

   :returns: Path to SPK cache directory
   :rtype: str


.. function:: set_spk_date_padding(days)

   Set the date padding for auto-downloaded SPK files.

   When auto-downloading, the SPK will cover the requested date plus/minus
   this padding on each side.

   :param days: Number of days to pad the date range
   :type days: int


.. function:: get_spk_date_padding()

   Get the current SPK date padding setting.

   :returns: Padding in days
   :rtype: int


Minor Body IDs
~~~~~~~~~~~~~~

Centaurs
^^^^^^^^

.. data:: SE_CHIRON
   :value: 15

   Chiron (2060)

.. data:: SE_PHOLUS
   :value: 16

   Pholus (5145)

.. data:: SE_NESSUS

   Nessus (7066) - Centaur

.. data:: SE_ASBOLUS

   Asbolus (8405) - Centaur

.. data:: SE_CHARIKLO

   Chariklo (10199) - Largest known centaur, has ring system


Trans-Neptunian Objects (TNOs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. data:: SE_ORCUS

   Orcus (90482) - Plutino

.. data:: SE_IXION

   Ixion (28978) - Plutino

.. data:: SE_HAUMEA

   Haumea (136108) - Dwarf planet

.. data:: SE_QUAOAR

   Quaoar (50000) - Classical Kuiper belt object

.. data:: SE_MAKEMAKE

   Makemake (136472) - Dwarf planet

.. data:: SE_GONGGONG

   Gonggong (225088) - TNO, dwarf planet candidate

.. data:: SE_ERIS

   Eris (136199) - Largest known dwarf planet

.. data:: SE_SEDNA

   Sedna (90377) - Detached TNO


NAIF ID Constants
~~~~~~~~~~~~~~~~~

NAIF IDs are used for SPK kernel registration. For numbered asteroids,
the NAIF ID is ``asteroid_number + 2000000``.

.. data:: NAIF_ASTEROID_OFFSET
   :value: 2000000

   Add asteroid number to get NAIF ID

.. data:: NAIF_CHIRON
   :value: 2002060

.. data:: NAIF_PHOLUS
   :value: 2005145

.. data:: NAIF_NESSUS
   :value: 2007066

.. data:: NAIF_ASBOLUS
   :value: 2008405

.. data:: NAIF_CHARIKLO
   :value: 2010199

.. data:: NAIF_CERES
   :value: 2000001

.. data:: NAIF_PALLAS
   :value: 2000002

.. data:: NAIF_JUNO
   :value: 2000003

.. data:: NAIF_VESTA
   :value: 2000004

.. data:: NAIF_ERIS
   :value: 2136199

.. data:: NAIF_SEDNA
   :value: 2090377

.. data:: NAIF_GONGGONG
   :value: 2225088


Minor Body Resonance Detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: detect_mean_motion_resonance(elements, tolerance=0.02)

   Detect if a minor body is in a mean motion resonance with Neptune.

   :param elements: Orbital elements of the body
   :type elements: OrbitalElements
   :param tolerance: Fractional tolerance for resonance detection (default 0.02 = 2%)
   :type tolerance: float
   :returns: ResonanceResult with resonance info, or None if not in resonance
   :rtype: ResonanceResult or None

   **Example:**

   >>> from libephemeris.minor_bodies import detect_mean_motion_resonance, MINOR_BODY_ELEMENTS
   >>> from libephemeris.constants import SE_IXION
   >>> result = detect_mean_motion_resonance(MINOR_BODY_ELEMENTS[SE_IXION])
   >>> if result:
   ...     print(f"{result.resonance.name}: {result.resonance.p}:{result.resonance.q}")
