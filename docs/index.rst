libephemeris Documentation
==========================

libephemeris is a pure-Python astronomical ephemeris library compatible with
the Swiss Ephemeris API. It provides high-precision planetary positions,
house calculations, eclipse predictions, and more.

Features
--------

- **Swiss Ephemeris Compatible**: Drop-in replacement for pyswisseph
- **Pure Python**: No C extensions required
- **High Precision**: Uses NASA JPL DE421 ephemeris via Skyfield
- **Thread Safe**: EphemerisContext for concurrent calculations
- **19 House Systems**: Including Placidus, Koch, Whole Sign, and more
- **43 Ayanamshas**: Full sidereal zodiac support
- **Eclipses**: Solar and lunar eclipse calculations
- **Fixed Stars**: Regulus, Spica, and more with proper motion
- **Minor Bodies**: Asteroids, centaurs, and TNOs with SPK kernel support
- **High-Precision SPK**: Download SPK kernels from JPL Horizons for arcsecond-level precision

Quick Start
-----------

.. code-block:: python

   import libephemeris as ephem

   # Calculate planetary positions
   jd = ephem.julday(2024, 1, 1, 12.0)
   pos, flags = ephem.calc_ut(jd, ephem.SE_SUN, ephem.SEFLG_SPEED)
   print(f"Sun longitude: {pos[0]:.4f} degrees")

   # Calculate house cusps
   cusps, ascmc = ephem.houses(jd, 41.9, 12.5, ord('P'))
   print(f"Ascendant: {ascmc[0]:.2f} degrees")

   # Sidereal calculations
   ephem.set_sid_mode(ephem.SE_SIDM_LAHIRI)
   pos_sid, _ = ephem.calc_ut(jd, ephem.SE_MOON, ephem.SEFLG_SIDEREAL)
   print(f"Moon (sidereal): {pos_sid[0]:.4f} degrees")

Installation
------------

.. code-block:: bash

   pip install libephemeris


Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started:

   migration-guide

.. toctree::
   :maxdepth: 2
   :caption: Technical Documentation:

   api_reference
   PRECISION
   INTERPOLATED_APOGEE
   TRUE_LILITH_METHODS


Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
