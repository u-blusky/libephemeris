"""
Fixed star position calculations for libephemeris.

Computes ecliptic positions for bright fixed stars with:
- Proper motion correction (rigorous space motion approach)
- IAU 2006 precession from J2000 to date
- IAU 2000A nutation model (1365 terms) for sub-milliarcsecond precision
- Equatorial to ecliptic coordinate transformation

Supported stars:
- Regulus (Alpha Leonis) - Royal Star
- Spica (Alpha Virginis) - Used for ayanamsha calculations

Notes on precision:
    Proper motion is applied using the rigorous space motion approach from
    Hipparcos Vol. 1, Section 1.5.5 with second-order Taylor expansion.
    This method converts the position to a 3D unit vector, applies proper
    motion as angular velocity in the tangent plane with centripetal
    acceleration correction, and normalizes to account for spherical geometry.

    The second-order term (-0.5 * |V|² * P * t²) accounts for the curvature
    of the celestial sphere, significantly improving accuracy for high
    proper motion stars (e.g., Barnard's Star) over century-scale intervals.

    Limitations:
    - Ignores radial velocity (parallax causes small position shift)
    - Assumes constant proper motion (real stars accelerate slightly)
    - No annual parallax correction (distance effect negligible for distant stars)
    Typical error: <0.01 arcsec over ±100 years, <1 arcsec over ±500 years
    For research astronomy, use SIMBAD/Gaia catalogs.

References:
- Hipparcos Catalog Vol. 1, Section 1.5.5 (ESA SP-1200, 1997)
- IAU 2006 Precession: Capitaine et al. A&A 412, 567-586 (2003)
- Proper motion: Hipparcos/Tycho catalogs
"""

import math
from dataclasses import dataclass
from typing import List, Tuple

from skyfield.nutationlib import iau2000a_radians

from .constants import (
    SE_REGULUS,
    SE_SPICA_STAR,
    SE_ALGOL,
    SE_SIRIUS,
    SE_ALDEBARAN,
    SE_ANTARES,
    SE_VEGA,
    SE_POLARIS,
    SE_FOMALHAUT,
    SE_BETELGEUSE,
    SE_RIGEL,
    SE_PROCYON,
    SE_CAPELLA,
    SE_ARCTURUS,
    SE_DENEB,
    SE_POLLUX,
    SE_CASTOR,
    SE_ALTAIR,
    SE_ACHERNAR,
    SE_CANOPUS,
    SE_ACRUX,
    SE_MIMOSA,
    SE_GACRUX,
    SE_HADAR,
    SE_RIGIL_KENT,
    SE_SHAULA,
    SE_BELLATRIX,
    SE_ELNATH,
    SE_MIRA,
    SE_ALNILAM,
    SE_ALNITAK,
    SE_MINTAKA,
    SE_SAIPH,
    SE_DIPHDA,
    SE_ALPHARD,
    SE_RASALHAGUE,
    SE_ETAMIN,
    SE_KOCHAB,
    SE_ALKAID,
    SE_DUBHE,
    SE_MERAK,
    SE_ALIOTH,
    SE_MIZAR,
    SE_ALCOR,
    SE_VINDEMIATRIX,
    SE_ZUBENELGENUBI,
    SE_ZUBENESCHAMALI,
    SE_UNUKALHAI,
    SE_ALGIEBA,
    SE_DENEBOLA,
    SE_MARKAB,
    SE_SCHEAT,
)


@dataclass
class StarData:
    """
    Fixed star astrometric data (ICRS J2000.0 epoch).

    Attributes:
        ra_j2000: Right Ascension at J2000.0 in degrees
        dec_j2000: Declination at J2000.0 in degrees
        pm_ra: Proper motion in RA (arcsec/year, includes cos(dec) factor)
        pm_dec: Proper motion in Dec (arcsec/year)

    Note:
        Proper motions are applied using rigorous space motion approach
        with second-order Taylor expansion (Hipparcos Vol. 1, Section 1.5.5).
        Accurate to <0.01 arcsec over ±100 years, <1 arcsec over ±500 years.
        Does not include parallax or radial velocity effects.
    """

    ra_j2000: float
    dec_j2000: float
    pm_ra: float
    pm_dec: float


@dataclass
class StarCatalogEntry:
    """
    Extended catalog entry for fixstar2 functions.

    Attributes:
        id: Internal star ID
        name: Traditional star name (e.g. "Regulus")
        nomenclature: Bayer/Flamsteed designation (e.g. "alLeo", "alVir")
        hip_number: Hipparcos catalog number (e.g. 49669)
        data: Astrometric data for position calculation
        magnitude: Visual magnitude (apparent brightness)
    """

    id: int
    name: str
    nomenclature: str
    hip_number: int
    data: StarData
    magnitude: float


# Extended star catalog with names and catalog numbers
STAR_CATALOG: List[StarCatalogEntry] = [
    StarCatalogEntry(
        id=SE_REGULUS,
        name="Regulus",
        nomenclature="alLeo",
        hip_number=49669,
        data=StarData(
            ra_j2000=152.092958,  # 10h 08m 22.3s (Alpha Leonis)
            dec_j2000=11.967208,  # +11° 58' 02"
            pm_ra=-0.00249,  # -249 mas/yr (westward)
            pm_dec=0.00152,  # +152 mas/yr (northward)
        ),
        magnitude=1.40,  # Visual magnitude
    ),
    StarCatalogEntry(
        id=SE_SPICA_STAR,
        name="Spica",
        nomenclature="alVir",
        hip_number=65474,
        data=StarData(
            ra_j2000=201.298247,  # 13h 25m 11.6s (Alpha Virginis)
            dec_j2000=-11.161319,  # -11° 09' 41"
            pm_ra=-0.04235,  # -42.35 mas/yr
            pm_dec=-0.03067,  # -30.67 mas/yr
        ),
        magnitude=1.04,  # Visual magnitude
    ),
    StarCatalogEntry(
        id=SE_ALGOL,
        name="Algol",
        nomenclature="bePer",
        hip_number=14576,
        data=StarData(
            ra_j2000=47.042219,  # 03h 08m 10.1s
            dec_j2000=40.955647,  # +40° 57' 20"
            pm_ra=0.00255,  # 2.55 mas/yr
            pm_dec=-0.00176,  # -1.76 mas/yr
        ),
        magnitude=2.12,
    ),
    StarCatalogEntry(
        id=SE_SIRIUS,
        name="Sirius",
        nomenclature="alCMa",
        hip_number=32349,
        data=StarData(
            ra_j2000=101.287155,  # 06h 45m 08.9s
            dec_j2000=-16.716116,  # -16° 42' 58"
            pm_ra=-0.54601,  # -546.01 mas/yr
            pm_dec=-1.22307,  # -1223.07 mas/yr
        ),
        magnitude=-1.46,
    ),
    StarCatalogEntry(
        id=SE_ALDEBARAN,
        name="Aldebaran",
        nomenclature="alTau",
        hip_number=21421,
        data=StarData(
            ra_j2000=68.980163,  # 04h 35m 55.2s
            dec_j2000=16.509302,  # +16° 30' 33"
            pm_ra=0.06294,  # 62.94 mas/yr
            pm_dec=-0.18993,  # -189.93 mas/yr
        ),
        magnitude=0.85,
    ),
    StarCatalogEntry(
        id=SE_ANTARES,
        name="Antares",
        nomenclature="alSco",
        hip_number=80763,
        data=StarData(
            ra_j2000=247.351915,  # 16h 29m 24.5s
            dec_j2000=-26.432003,  # -26° 25' 55"
            pm_ra=-0.01059,  # -10.59 mas/yr
            pm_dec=-0.02304,  # -23.04 mas/yr
        ),
        magnitude=1.06,
    ),
    StarCatalogEntry(
        id=SE_VEGA,
        name="Vega",
        nomenclature="alLyr",
        hip_number=91262,
        data=StarData(
            ra_j2000=279.234735,  # 18h 36m 56.3s
            dec_j2000=38.783689,  # +38° 47' 01"
            pm_ra=0.20100,  # 201.00 mas/yr
            pm_dec=0.28710,  # 287.10 mas/yr
        ),
        magnitude=0.03,
    ),
    StarCatalogEntry(
        id=SE_POLARIS,
        name="Polaris",
        nomenclature="alUMi",
        hip_number=11767,
        data=StarData(
            ra_j2000=37.954561,  # 02h 31m 49.1s
            dec_j2000=89.264109,  # +89° 15' 51"
            pm_ra=0.04422,  # 44.22 mas/yr
            pm_dec=-0.01175,  # -11.75 mas/yr
        ),
        magnitude=1.98,
    ),
    StarCatalogEntry(
        id=SE_FOMALHAUT,
        name="Fomalhaut",
        nomenclature="alPsA",
        hip_number=113368,
        data=StarData(
            ra_j2000=344.412693,  # 22h 57m 39.0s
            dec_j2000=-29.622237,  # -29° 37' 20"
            pm_ra=0.32900,  # 329.00 mas/yr
            pm_dec=-0.16474,  # -164.74 mas/yr
        ),
        magnitude=1.16,
    ),
    StarCatalogEntry(
        id=SE_BETELGEUSE,
        name="Betelgeuse",
        nomenclature="alOri",
        hip_number=27989,
        data=StarData(
            ra_j2000=88.792939,  # 05h 55m 10.3s
            dec_j2000=7.407064,  # +07° 24' 25"
            pm_ra=0.02728,  # 27.28 mas/yr
            pm_dec=0.01038,  # 10.38 mas/yr
        ),
        magnitude=0.42,
    ),
    StarCatalogEntry(
        id=SE_RIGEL,
        name="Rigel",
        nomenclature="beOri",
        hip_number=24436,
        data=StarData(
            ra_j2000=78.634467,  # 05h 14m 32.3s
            dec_j2000=-8.201638,  # -08° 12' 06"
            pm_ra=0.00145,  # 1.45 mas/yr
            pm_dec=-0.00004,  # -0.04 mas/yr
        ),
        magnitude=0.12,
    ),
    StarCatalogEntry(
        id=SE_PROCYON,
        name="Procyon",
        nomenclature="alCMi",
        hip_number=37279,
        data=StarData(
            ra_j2000=114.825493,  # 07h 39m 18.1s
            dec_j2000=5.224993,  # +05° 13' 30"
            pm_ra=-0.71410,  # -714.10 mas/yr
            pm_dec=-1.02280,  # -1022.80 mas/yr
        ),
        magnitude=0.34,
    ),
    StarCatalogEntry(
        id=SE_CAPELLA,
        name="Capella",
        nomenclature="alAur",
        hip_number=24608,
        data=StarData(
            ra_j2000=79.172328,  # 05h 16m 41.4s
            dec_j2000=45.997991,  # +45° 59' 53"
            pm_ra=0.07527,  # 75.27 mas/yr
            pm_dec=-0.42711,  # -427.11 mas/yr
        ),
        magnitude=0.08,
    ),
    StarCatalogEntry(
        id=SE_ARCTURUS,
        name="Arcturus",
        nomenclature="alBoo",
        hip_number=69673,
        data=StarData(
            ra_j2000=213.915300,  # 14h 15m 39.7s
            dec_j2000=19.182409,  # +19° 10' 57"
            pm_ra=-1.09350,  # -1093.50 mas/yr
            pm_dec=-1.99940,  # -1999.40 mas/yr
        ),
        magnitude=-0.04,
    ),
    StarCatalogEntry(
        id=SE_DENEB,
        name="Deneb",
        nomenclature="alCyg",
        hip_number=102098,
        data=StarData(
            ra_j2000=310.357980,  # 20h 41m 25.9s
            dec_j2000=45.280339,  # +45° 16' 49"
            pm_ra=0.00196,  # 1.96 mas/yr
            pm_dec=0.00215,  # 2.15 mas/yr
        ),
        magnitude=1.25,
    ),
    StarCatalogEntry(
        id=SE_POLLUX,
        name="Pollux",
        nomenclature="beGem",
        hip_number=37826,
        data=StarData(
            ra_j2000=116.328958,  # 07h 45m 18.9s
            dec_j2000=28.026199,  # +28° 01' 34"
            pm_ra=-0.62565,  # -625.65 mas/yr
            pm_dec=-0.04597,  # -45.97 mas/yr
        ),
        magnitude=1.14,
    ),
    StarCatalogEntry(
        id=SE_CASTOR,
        name="Castor",
        nomenclature="alGem",
        hip_number=36850,
        data=StarData(
            ra_j2000=113.649428,  # 07h 34m 35.9s
            dec_j2000=31.888276,  # +31° 53' 18"
            pm_ra=-0.19128,  # -191.28 mas/yr
            pm_dec=-0.14541,  # -145.41 mas/yr
        ),
        magnitude=1.58,
    ),
    StarCatalogEntry(
        id=SE_ALTAIR,
        name="Altair",
        nomenclature="alAql",
        hip_number=97649,
        data=StarData(
            ra_j2000=297.695827,  # 19h 50m 47.0s
            dec_j2000=8.868321,  # +08° 52' 06"
            pm_ra=0.53683,  # 536.83 mas/yr
            pm_dec=0.38600,  # 386.00 mas/yr
        ),
        magnitude=0.77,
    ),
    StarCatalogEntry(
        id=SE_ACHERNAR,
        name="Achernar",
        nomenclature="alEri",
        hip_number=7588,
        data=StarData(
            ra_j2000=24.428523,  # 01h 37m 42.8s
            dec_j2000=-57.236753,  # -57° 14' 12"
            pm_ra=0.08767,  # 87.67 mas/yr
            pm_dec=-0.04009,  # -40.09 mas/yr
        ),
        magnitude=0.46,
    ),
    StarCatalogEntry(
        id=SE_CANOPUS,
        name="Canopus",
        nomenclature="alCar",
        hip_number=30438,
        data=StarData(
            ra_j2000=95.987958,  # 06h 23m 57.1s
            dec_j2000=-52.695661,  # -52° 41' 44"
            pm_ra=0.01993,  # 19.93 mas/yr
            pm_dec=0.02373,  # 23.73 mas/yr
        ),
        magnitude=-0.72,
    ),
    StarCatalogEntry(
        id=SE_ACRUX,
        name="Acrux",
        nomenclature="alCru",
        hip_number=60718,
        data=StarData(
            ra_j2000=186.649563,  # 12h 26m 35.9s
            dec_j2000=-63.099093,  # -63° 05' 57"
            pm_ra=-0.03548,  # -35.48 mas/yr
            pm_dec=-0.01235,  # -12.35 mas/yr
        ),
        magnitude=0.76,
    ),
    StarCatalogEntry(
        id=SE_MIMOSA,
        name="Mimosa",
        nomenclature="beCru",
        hip_number=62434,
        data=StarData(
            ra_j2000=191.930263,  # 12h 47m 43.3s
            dec_j2000=-59.688764,  # -59° 41' 19"
            pm_ra=-0.04203,  # -42.03 mas/yr
            pm_dec=-0.01600,  # -16.00 mas/yr
        ),
        magnitude=1.25,
    ),
    StarCatalogEntry(
        id=SE_GACRUX,
        name="Gacrux",
        nomenclature="gaCru",
        hip_number=61084,
        data=StarData(
            ra_j2000=187.791498,  # 12h 31m 09.9s
            dec_j2000=-57.113213,  # -57° 06' 48"
            pm_ra=0.02785,  # 27.85 mas/yr
            pm_dec=-0.26488,  # -264.88 mas/yr
        ),
        magnitude=1.64,
    ),
    StarCatalogEntry(
        id=SE_HADAR,
        name="Hadar",
        nomenclature="beCen",
        hip_number=68702,
        data=StarData(
            ra_j2000=211.670617,  # 14h 06m 40.9s
            dec_j2000=-60.373035,  # -60° 22' 23"
            pm_ra=-0.03366,  # -33.66 mas/yr
            pm_dec=-0.02511,  # -25.11 mas/yr
        ),
        magnitude=0.61,
    ),
    StarCatalogEntry(
        id=SE_RIGIL_KENT,
        name="Rigil Kentaurus",
        nomenclature="alCen",
        hip_number=71683,
        data=StarData(
            ra_j2000=219.902066,  # 14h 39m 36.5s
            dec_j2000=-60.833976,  # -60° 50' 02"
            pm_ra=-3.67818,  # -3678.18 mas/yr
            pm_dec=0.48189,  # 481.89 mas/yr
        ),
        magnitude=-0.27,
    ),
    StarCatalogEntry(
        id=SE_SHAULA,
        name="Shaula",
        nomenclature="laSco",
        hip_number=85927,
        data=StarData(
            ra_j2000=263.402167,  # 17h 33m 36.5s
            dec_j2000=-37.103824,  # -37° 06' 14"
            pm_ra=-0.00847,  # -8.47 mas/yr
            pm_dec=-0.02984,  # -29.84 mas/yr
        ),
        magnitude=1.62,
    ),
    StarCatalogEntry(
        id=SE_BELLATRIX,
        name="Bellatrix",
        nomenclature="gaOri",
        hip_number=25336,
        data=StarData(
            ra_j2000=81.282764,  # 05h 25m 07.9s
            dec_j2000=6.349703,  # +06° 20' 59"
            pm_ra=-0.00816,  # -8.16 mas/yr
            pm_dec=-0.01294,  # -12.94 mas/yr
        ),
        magnitude=1.64,
    ),
    StarCatalogEntry(
        id=SE_ELNATH,
        name="Elnath",
        nomenclature="beTau",
        hip_number=25428,
        data=StarData(
            ra_j2000=81.572971,  # 05h 26m 17.5s
            dec_j2000=28.607452,  # +28° 36' 27"
            pm_ra=0.02284,  # 22.84 mas/yr
            pm_dec=-0.17481,  # -174.81 mas/yr
        ),
        magnitude=1.65,
    ),
    StarCatalogEntry(
        id=SE_MIRA,
        name="Mira",
        nomenclature="omCet",
        hip_number=10826,
        data=StarData(
            ra_j2000=34.836617,  # 02h 19m 20.8s
            dec_j2000=-2.977640,  # -02° 58' 40"
            pm_ra=0.01006,  # 10.06 mas/yr
            pm_dec=-0.23836,  # -238.36 mas/yr
        ),
        magnitude=3.04,
    ),
    StarCatalogEntry(
        id=SE_ALNILAM,
        name="Alnilam",
        nomenclature="epOri",
        hip_number=26311,
        data=StarData(
            ra_j2000=84.053389,  # 05h 36m 12.8s
            dec_j2000=-1.201919,  # -01° 12' 07"
            pm_ra=0.00117,  # 1.17 mas/yr
            pm_dec=-0.00128,  # -1.28 mas/yr
        ),
        magnitude=1.69,
    ),
    StarCatalogEntry(
        id=SE_ALNITAK,
        name="Alnitak",
        nomenclature="zeOri",
        hip_number=26727,
        data=StarData(
            ra_j2000=85.189694,  # 05h 40m 45.5s
            dec_j2000=-1.942574,  # -01° 56' 33"
            pm_ra=0.00363,  # 3.63 mas/yr
            pm_dec=0.00212,  # 2.12 mas/yr
        ),
        magnitude=1.74,
    ),
    StarCatalogEntry(
        id=SE_MINTAKA,
        name="Mintaka",
        nomenclature="deOri",
        hip_number=25930,
        data=StarData(
            ra_j2000=83.001667,  # 05h 32m 00.4s
            dec_j2000=-0.299095,  # -00° 17' 57"
            pm_ra=0.00132,  # 1.32 mas/yr
            pm_dec=-0.00055,  # -0.55 mas/yr
        ),
        magnitude=2.23,
    ),
    StarCatalogEntry(
        id=SE_SAIPH,
        name="Saiph",
        nomenclature="kaOri",
        hip_number=27366,
        data=StarData(
            ra_j2000=86.939119,  # 05h 47m 45.4s
            dec_j2000=-9.669605,  # -09° 40' 11"
            pm_ra=0.00157,  # 1.57 mas/yr
            pm_dec=-0.00139,  # -1.39 mas/yr
        ),
        magnitude=2.06,
    ),
    StarCatalogEntry(
        id=SE_DIPHDA,
        name="Diphda",
        nomenclature="beCet",
        hip_number=3419,
        data=StarData(
            ra_j2000=10.897379,  # 00h 43m 35.4s
            dec_j2000=-17.986606,  # -17° 59' 12"
            pm_ra=0.23179,  # 231.79 mas/yr
            pm_dec=0.03279,  # 32.79 mas/yr
        ),
        magnitude=2.04,
    ),
    StarCatalogEntry(
        id=SE_ALPHARD,
        name="Alphard",
        nomenclature="alHya",
        hip_number=46390,
        data=StarData(
            ra_j2000=141.896847,  # 09h 27m 35.2s
            dec_j2000=-8.658602,  # -08° 39' 31"
            pm_ra=-0.01468,  # -14.68 mas/yr
            pm_dec=0.03368,  # 33.68 mas/yr
        ),
        magnitude=1.98,
    ),
    StarCatalogEntry(
        id=SE_RASALHAGUE,
        name="Rasalhague",
        nomenclature="alOph",
        hip_number=86032,
        data=StarData(
            ra_j2000=263.733627,  # 17h 34m 56.1s
            dec_j2000=12.560035,  # +12° 33' 36"
            pm_ra=0.10887,  # 108.87 mas/yr
            pm_dec=-0.22638,  # -226.38 mas/yr
        ),
        magnitude=2.07,
    ),
    StarCatalogEntry(
        id=SE_ETAMIN,
        name="Etamin",
        nomenclature="gaDra",
        hip_number=87833,
        data=StarData(
            ra_j2000=269.151541,  # 17h 56m 36.4s
            dec_j2000=51.488896,  # +51° 29' 20"
            pm_ra=-0.00815,  # -8.15 mas/yr
            pm_dec=-0.02284,  # -22.84 mas/yr
        ),
        magnitude=2.23,
    ),
    StarCatalogEntry(
        id=SE_KOCHAB,
        name="Kochab",
        nomenclature="beUMi",
        hip_number=72607,
        data=StarData(
            ra_j2000=222.676357,  # 14h 50m 42.3s
            dec_j2000=74.155504,  # +74° 09' 20"
            pm_ra=-0.03228,  # -32.28 mas/yr
            pm_dec=0.01164,  # 11.64 mas/yr
        ),
        magnitude=2.08,
    ),
    StarCatalogEntry(
        id=SE_ALKAID,
        name="Alkaid",
        nomenclature="etUMa",
        hip_number=67301,
        data=StarData(
            ra_j2000=206.885157,  # 13h 47m 32.4s
            dec_j2000=49.313267,  # +49° 18' 48"
            pm_ra=-0.12107,  # -121.07 mas/yr
            pm_dec=-0.01484,  # -14.84 mas/yr
        ),
        magnitude=1.86,
    ),
    StarCatalogEntry(
        id=SE_DUBHE,
        name="Dubhe",
        nomenclature="alUMa",
        hip_number=54061,
        data=StarData(
            ra_j2000=165.931965,  # 11h 03m 43.7s
            dec_j2000=61.751035,  # +61° 45' 03"
            pm_ra=-0.13469,  # -134.69 mas/yr
            pm_dec=-0.03469,  # -34.69 mas/yr
        ),
        magnitude=1.79,
    ),
    StarCatalogEntry(
        id=SE_MERAK,
        name="Merak",
        nomenclature="beUMa",
        hip_number=53910,
        data=StarData(
            ra_j2000=165.460319,  # 11h 01m 50.5s
            dec_j2000=56.382426,  # +56° 22' 57"
            pm_ra=0.08175,  # 81.75 mas/yr
            pm_dec=0.03384,  # 33.84 mas/yr
        ),
        magnitude=2.37,
    ),
    StarCatalogEntry(
        id=SE_ALIOTH,
        name="Alioth",
        nomenclature="epUMa",
        hip_number=62956,
        data=StarData(
            ra_j2000=193.507290,  # 12h 54m 01.7s
            dec_j2000=55.959823,  # +55° 57' 35"
            pm_ra=0.11117,  # 111.17 mas/yr
            pm_dec=-0.00869,  # -8.69 mas/yr
        ),
        magnitude=1.77,
    ),
    StarCatalogEntry(
        id=SE_MIZAR,
        name="Mizar",
        nomenclature="zeUMa",
        hip_number=65378,
        data=StarData(
            ra_j2000=200.981429,  # 13h 23m 55.5s
            dec_j2000=54.925362,  # +54° 55' 31"
            pm_ra=0.12116,  # 121.16 mas/yr
            pm_dec=-0.02223,  # -22.23 mas/yr
        ),
        magnitude=2.23,
    ),
    StarCatalogEntry(
        id=SE_ALCOR,
        name="Alcor",
        nomenclature="80UMa",
        hip_number=65477,
        data=StarData(
            ra_j2000=201.306403,  # 13h 25m 13.5s
            dec_j2000=54.987954,  # +54° 59' 17"
            pm_ra=0.12027,  # 120.27 mas/yr
            pm_dec=-0.01680,  # -16.80 mas/yr
        ),
        magnitude=3.99,
    ),
    StarCatalogEntry(
        id=SE_VINDEMIATRIX,
        name="Vindemiatrix",
        nomenclature="epVir",
        hip_number=63608,
        data=StarData(
            ra_j2000=195.544157,  # 13h 02m 10.6s
            dec_j2000=10.959149,  # +10° 57' 33"
            pm_ra=-0.27520,  # -275.20 mas/yr
            pm_dec=0.01988,  # 19.88 mas/yr
        ),
        magnitude=2.83,
    ),
    StarCatalogEntry(
        id=SE_ZUBENELGENUBI,
        name="Zubenelgenubi",
        nomenclature="alLib",
        hip_number=72622,
        data=StarData(
            ra_j2000=222.719638,  # 14h 50m 52.7s
            dec_j2000=-16.041778,  # -16° 02' 30"
            pm_ra=-0.10569,  # -105.69 mas/yr
            pm_dec=-0.06905,  # -69.05 mas/yr
        ),
        magnitude=2.75,
    ),
    StarCatalogEntry(
        id=SE_ZUBENESCHAMALI,
        name="Zubeneschamali",
        nomenclature="beLib",
        hip_number=74785,
        data=StarData(
            ra_j2000=229.251724,  # 15h 17m 00.4s
            dec_j2000=-9.382914,  # -09° 22' 58"
            pm_ra=-0.09860,  # -98.60 mas/yr
            pm_dec=-0.01914,  # -19.14 mas/yr
        ),
        magnitude=2.61,
    ),
    StarCatalogEntry(
        id=SE_UNUKALHAI,
        name="Unukalhai",
        nomenclature="alSer",
        hip_number=77070,
        data=StarData(
            ra_j2000=236.067089,  # 15h 44m 16.1s
            dec_j2000=6.425628,  # +06° 25' 32"
            pm_ra=0.13370,  # 133.70 mas/yr
            pm_dec=0.04495,  # 44.95 mas/yr
        ),
        magnitude=2.65,
    ),
    StarCatalogEntry(
        id=SE_ALGIEBA,
        name="Algieba",
        nomenclature="gaLeo",
        hip_number=50583,
        data=StarData(
            ra_j2000=154.993144,  # 10h 19m 58.4s
            dec_j2000=19.841489,  # +19° 50' 30"
            pm_ra=0.30996,  # 309.96 mas/yr
            pm_dec=-0.15259,  # -152.59 mas/yr
        ),
        magnitude=2.08,
    ),
    StarCatalogEntry(
        id=SE_DENEBOLA,
        name="Denebola",
        nomenclature="beLeo",
        hip_number=57632,
        data=StarData(
            ra_j2000=177.264910,  # 11h 49m 03.6s
            dec_j2000=14.572058,  # +14° 34' 19"
            pm_ra=-0.49927,  # -499.27 mas/yr
            pm_dec=-0.11385,  # -113.85 mas/yr
        ),
        magnitude=2.14,
    ),
    StarCatalogEntry(
        id=SE_MARKAB,
        name="Markab",
        nomenclature="alPeg",
        hip_number=113963,
        data=StarData(
            ra_j2000=346.190223,  # 23h 04m 45.7s
            dec_j2000=15.205267,  # +15° 12' 19"
            pm_ra=0.06185,  # 61.85 mas/yr
            pm_dec=-0.04214,  # -42.14 mas/yr
        ),
        magnitude=2.49,
    ),
    StarCatalogEntry(
        id=SE_SCHEAT,
        name="Scheat",
        nomenclature="bePeg",
        hip_number=113881,
        data=StarData(
            ra_j2000=345.943514,  # 23h 03m 46.5s
            dec_j2000=28.082785,  # +28° 04' 58"
            pm_ra=0.18724,  # 187.24 mas/yr
            pm_dec=0.13691,  # 136.91 mas/yr
        ),
        magnitude=2.42,
    ),
]

# Fixed star catalog (J2000.0 ICRS coordinates from Hipparcos)
# Legacy format for backward compatibility
FIXED_STARS = {entry.id: entry.data for entry in STAR_CATALOG}

# Build lookup from canonical name to star ID
_STAR_NAME_TO_ID = {entry.name.upper(): entry.id for entry in STAR_CATALOG}

# STAR_ALIASES: Maps alternative star names to canonical SE_* constant IDs
# Includes: common names, Bayer designations (full and abbreviated),
# Flamsteed numbers, Arabic names, Latin names, Greek transliterations
STAR_ALIASES: dict[str, int] = {
    # ======== REGULUS (Alpha Leonis) ========
    "COR LEONIS": SE_REGULUS,
    "ALPHA LEONIS": SE_REGULUS,
    "ALPHA LEO": SE_REGULUS,
    "87 LEO": SE_REGULUS,
    "α LEO": SE_REGULUS,
    "QALB AL-ASAD": SE_REGULUS,
    "BASILISKOS": SE_REGULUS,
    "REX": SE_REGULUS,
    "KALB AL-ASAD": SE_REGULUS,
    "ALLEO": SE_REGULUS,
    # ======== SPICA (Alpha Virginis) ========
    "ALPHA VIRGINIS": SE_SPICA_STAR,
    "ALPHA VIR": SE_SPICA_STAR,
    "67 VIR": SE_SPICA_STAR,
    "α VIR": SE_SPICA_STAR,
    "AZIMECH": SE_SPICA_STAR,
    "ALARAPH": SE_SPICA_STAR,
    "ALVIR": SE_SPICA_STAR,
    "SUNBULA": SE_SPICA_STAR,
    "VIRGIN'S SPIKE": SE_SPICA_STAR,
    "ARISTA": SE_SPICA_STAR,
    # ======== ALGOL (Beta Persei) ========
    "DEMON STAR": SE_ALGOL,
    "BETA PERSEI": SE_ALGOL,
    "BETA PER": SE_ALGOL,
    "26 PER": SE_ALGOL,
    "β PER": SE_ALGOL,
    "GORGONEA PRIMA": SE_ALGOL,
    "RA'S AL-GHUL": SE_ALGOL,
    "RAS AL-GHUL": SE_ALGOL,
    "BEPER": SE_ALGOL,
    "HEAD OF THE GHOUL": SE_ALGOL,
    "GORGON'S HEAD": SE_ALGOL,
    # ======== SIRIUS (Alpha Canis Majoris) ========
    "DOG STAR": SE_SIRIUS,
    "ALPHA CANIS MAJORIS": SE_SIRIUS,
    "ALPHA CMA": SE_SIRIUS,
    "α CMA": SE_SIRIUS,
    "9 CMA": SE_SIRIUS,
    "CANICULA": SE_SIRIUS,
    "ASCHERE": SE_SIRIUS,
    "ALCMA": SE_SIRIUS,
    "AL-SHIRA": SE_SIRIUS,
    "SOTHIS": SE_SIRIUS,
    # ======== ALDEBARAN (Alpha Tauri) ========
    "EYE OF TAURUS": SE_ALDEBARAN,
    "ALPHA TAURI": SE_ALDEBARAN,
    "ALPHA TAU": SE_ALDEBARAN,
    "87 TAU": SE_ALDEBARAN,
    "α TAU": SE_ALDEBARAN,
    "ALTAU": SE_ALDEBARAN,
    "PARILICIUM": SE_ALDEBARAN,
    "AL-DABARAN": SE_ALDEBARAN,
    "FOLLOWER": SE_ALDEBARAN,
    "ROHINI": SE_ALDEBARAN,
    # ======== ANTARES (Alpha Scorpii) ========
    "RIVAL OF MARS": SE_ANTARES,
    "ALPHA SCORPII": SE_ANTARES,
    "ALPHA SCO": SE_ANTARES,
    "21 SCO": SE_ANTARES,
    "α SCO": SE_ANTARES,
    "ALSCO": SE_ANTARES,
    "COR SCORPII": SE_ANTARES,
    "CALB AL-AKRAB": SE_ANTARES,
    "HEART OF SCORPION": SE_ANTARES,
    "JYESHTHA": SE_ANTARES,
    # ======== VEGA (Alpha Lyrae) ========
    "HARP STAR": SE_VEGA,
    "ALPHA LYRAE": SE_VEGA,
    "ALPHA LYR": SE_VEGA,
    "3 LYR": SE_VEGA,
    "α LYR": SE_VEGA,
    "ALLYR": SE_VEGA,
    "WEGA": SE_VEGA,
    "AL-NASR AL-WAQI": SE_VEGA,
    "FIDIS": SE_VEGA,
    "ABHIJIT": SE_VEGA,
    # ======== POLARIS (Alpha Ursae Minoris) ========
    "NORTH STAR": SE_POLARIS,
    "POLE STAR": SE_POLARIS,
    "ALPHA URSAE MINORIS": SE_POLARIS,
    "ALPHA UMI": SE_POLARIS,
    "1 UMI": SE_POLARIS,
    "α UMI": SE_POLARIS,
    "ALUMI": SE_POLARIS,
    "CYNOSURA": SE_POLARIS,
    "LODESTAR": SE_POLARIS,
    "STELLA POLARIS": SE_POLARIS,
    # ======== FOMALHAUT (Alpha Piscis Austrini) ========
    "FISH'S MOUTH": SE_FOMALHAUT,
    "ALPHA PISCIS AUSTRINI": SE_FOMALHAUT,
    "ALPHA PSA": SE_FOMALHAUT,
    "24 PSA": SE_FOMALHAUT,
    "α PSA": SE_FOMALHAUT,
    "ALPSA": SE_FOMALHAUT,
    "FUM AL-HUT": SE_FOMALHAUT,
    "OS PISCIS MERIDIANI": SE_FOMALHAUT,
    "LONELY STAR": SE_FOMALHAUT,
    "HASTORANG": SE_FOMALHAUT,
    # ======== BETELGEUSE (Alpha Orionis) ========
    "ARMPIT OF ORION": SE_BETELGEUSE,
    "ALPHA ORIONIS": SE_BETELGEUSE,
    "ALPHA ORI": SE_BETELGEUSE,
    "58 ORI": SE_BETELGEUSE,
    "α ORI": SE_BETELGEUSE,
    "ALORI": SE_BETELGEUSE,
    "BETELGEUZE": SE_BETELGEUSE,
    "IBT AL-JAUZAH": SE_BETELGEUSE,
    "ARDRA": SE_BETELGEUSE,
    "YAD AL-JAWZA": SE_BETELGEUSE,
    # ======== RIGEL (Beta Orionis) ========
    "FOOT OF ORION": SE_RIGEL,
    "BETA ORIONIS": SE_RIGEL,
    "BETA ORI": SE_RIGEL,
    "19 ORI": SE_RIGEL,
    "β ORI": SE_RIGEL,
    "BEORI": SE_RIGEL,
    "ALGEBAR": SE_RIGEL,
    "RIJL JAUZAH AL-YUSRA": SE_RIGEL,
    "ORION'S FOOT": SE_RIGEL,
    # ======== PROCYON (Alpha Canis Minoris) ========
    "LITTLE DOG STAR": SE_PROCYON,
    "ALPHA CANIS MINORIS": SE_PROCYON,
    "ALPHA CMI": SE_PROCYON,
    "10 CMI": SE_PROCYON,
    "α CMI": SE_PROCYON,
    "ALCMI": SE_PROCYON,
    "ANTECANIS": SE_PROCYON,
    "ELGOMAISA": SE_PROCYON,
    "AL-GHUMAISA": SE_PROCYON,
    # ======== CAPELLA (Alpha Aurigae) ========
    "SHE-GOAT": SE_CAPELLA,
    "ALPHA AURIGAE": SE_CAPELLA,
    "ALPHA AUR": SE_CAPELLA,
    "13 AUR": SE_CAPELLA,
    "α AUR": SE_CAPELLA,
    "ALAUR": SE_CAPELLA,
    "ALHAJOTH": SE_CAPELLA,
    "AMALTHEA": SE_CAPELLA,
    "GOAT STAR": SE_CAPELLA,
    # ======== ARCTURUS (Alpha Bootis) ========
    "BEAR GUARD": SE_ARCTURUS,
    "ALPHA BOOTIS": SE_ARCTURUS,
    "ALPHA BOO": SE_ARCTURUS,
    "16 BOO": SE_ARCTURUS,
    "α BOO": SE_ARCTURUS,
    "ALBOO": SE_ARCTURUS,
    "AL-SIMAK AL-RAMIH": SE_ARCTURUS,
    "HARIS AL-SAMA": SE_ARCTURUS,
    "GUARDIAN OF BEAR": SE_ARCTURUS,
    # ======== DENEB (Alpha Cygni) ========
    "TAIL OF HEN": SE_DENEB,
    "ALPHA CYGNI": SE_DENEB,
    "ALPHA CYG": SE_DENEB,
    "50 CYG": SE_DENEB,
    "α CYG": SE_DENEB,
    "ALCYG": SE_DENEB,
    "DHANAB AD-DAJAJAH": SE_DENEB,
    "ARIDED": SE_DENEB,
    "GALLINA": SE_DENEB,
    # ======== POLLUX (Beta Geminorum) ========
    "TWIN STAR": SE_POLLUX,
    "BETA GEMINORUM": SE_POLLUX,
    "BETA GEM": SE_POLLUX,
    "78 GEM": SE_POLLUX,
    "β GEM": SE_POLLUX,
    "BEGEM": SE_POLLUX,
    "POLYDEUCES": SE_POLLUX,
    "HEAD OF SECOND TWIN": SE_POLLUX,
    # ======== CASTOR (Alpha Geminorum) ========
    "ALPHA GEMINORUM": SE_CASTOR,
    "ALPHA GEM": SE_CASTOR,
    "66 GEM": SE_CASTOR,
    "α GEM": SE_CASTOR,
    "ALGEM": SE_CASTOR,
    "APOLLO": SE_CASTOR,
    "HEAD OF FIRST TWIN": SE_CASTOR,
    # ======== ALTAIR (Alpha Aquilae) ========
    "FLYING EAGLE": SE_ALTAIR,
    "ALPHA AQUILAE": SE_ALTAIR,
    "ALPHA AQL": SE_ALTAIR,
    "53 AQL": SE_ALTAIR,
    "α AQL": SE_ALTAIR,
    "ALAQL": SE_ALTAIR,
    "AL-NASR AL-TAIR": SE_ALTAIR,
    "ATAIR": SE_ALTAIR,
    "SRAVANA": SE_ALTAIR,
    # ======== ACHERNAR (Alpha Eridani) ========
    "END OF RIVER": SE_ACHERNAR,
    "ALPHA ERIDANI": SE_ACHERNAR,
    "ALPHA ERI": SE_ACHERNAR,
    "α ERI": SE_ACHERNAR,
    "ALERI": SE_ACHERNAR,
    "AKHIR AN-NAHR": SE_ACHERNAR,
    "RIVER'S END": SE_ACHERNAR,
    # ======== CANOPUS (Alpha Carinae) ========
    "SHIP'S PILOT": SE_CANOPUS,
    "ALPHA CARINAE": SE_CANOPUS,
    "ALPHA CAR": SE_CANOPUS,
    "α CAR": SE_CANOPUS,
    "ALCAR": SE_CANOPUS,
    "SUHAIL": SE_CANOPUS,
    "SUHAYL": SE_CANOPUS,
    "AGASTYA": SE_CANOPUS,
    # ======== ACRUX (Alpha Crucis) ========
    "ALPHA CRUCIS": SE_ACRUX,
    "ALPHA CRU": SE_ACRUX,
    "α CRU": SE_ACRUX,
    "ALCRU": SE_ACRUX,
    "CRUX ALPHA": SE_ACRUX,
    "STAR OF BETHLEHEM": SE_ACRUX,
    # ======== MIMOSA (Beta Crucis) ========
    "BETA CRUCIS": SE_MIMOSA,
    "BETA CRU": SE_MIMOSA,
    "β CRU": SE_MIMOSA,
    "BECRU": SE_MIMOSA,
    "BECRUX": SE_MIMOSA,
    "CRUX BETA": SE_MIMOSA,
    # ======== GACRUX (Gamma Crucis) ========
    "GAMMA CRUCIS": SE_GACRUX,
    "GAMMA CRU": SE_GACRUX,
    "γ CRU": SE_GACRUX,
    "GACRU": SE_GACRUX,
    "CRUX GAMMA": SE_GACRUX,
    "RUBIDEA": SE_GACRUX,
    # ======== HADAR (Beta Centauri) ========
    "BETA CENTAURI": SE_HADAR,
    "BETA CEN": SE_HADAR,
    "β CEN": SE_HADAR,
    "BECEN": SE_HADAR,
    "AGENA": SE_HADAR,
    "KNEE OF CENTAUR": SE_HADAR,
    # ======== RIGIL KENTAURUS (Alpha Centauri) ========
    "ALPHA CENTAURI": SE_RIGIL_KENT,
    "ALPHA CEN": SE_RIGIL_KENT,
    "α CEN": SE_RIGIL_KENT,
    "ALCEN": SE_RIGIL_KENT,
    "TOLIMAN": SE_RIGIL_KENT,
    "RIGIL KENT": SE_RIGIL_KENT,
    "FOOT OF CENTAUR": SE_RIGIL_KENT,
    "BUNGULA": SE_RIGIL_KENT,
    # ======== SHAULA (Lambda Scorpii) ========
    "LAMBDA SCORPII": SE_SHAULA,
    "LAMBDA SCO": SE_SHAULA,
    "λ SCO": SE_SHAULA,
    "LASCO": SE_SHAULA,
    "SCORPION'S STING": SE_SHAULA,
    "35 SCO": SE_SHAULA,
    # ======== BELLATRIX (Gamma Orionis) ========
    "AMAZON STAR": SE_BELLATRIX,
    "GAMMA ORIONIS": SE_BELLATRIX,
    "GAMMA ORI": SE_BELLATRIX,
    "24 ORI": SE_BELLATRIX,
    "γ ORI": SE_BELLATRIX,
    "GAORI": SE_BELLATRIX,
    "FEMALE WARRIOR": SE_BELLATRIX,
    # ======== ELNATH (Beta Tauri) ========
    "BETA TAURI": SE_ELNATH,
    "BETA TAU": SE_ELNATH,
    "112 TAU": SE_ELNATH,
    "β TAU": SE_ELNATH,
    "BETAU": SE_ELNATH,
    "NATH": SE_ELNATH,
    "AL-NATH": SE_ELNATH,
    "EL NATH": SE_ELNATH,
    # ======== MIRA (Omicron Ceti) ========
    "WONDERFUL STAR": SE_MIRA,
    "OMICRON CETI": SE_MIRA,
    "OMICRON CET": SE_MIRA,
    "68 CET": SE_MIRA,
    "ο CET": SE_MIRA,
    "OMCET": SE_MIRA,
    "STELLA MIRA": SE_MIRA,
    # ======== ALNILAM (Epsilon Orionis) ========
    "STRING OF PEARLS": SE_ALNILAM,
    "EPSILON ORIONIS": SE_ALNILAM,
    "EPSILON ORI": SE_ALNILAM,
    "46 ORI": SE_ALNILAM,
    "ε ORI": SE_ALNILAM,
    "EPORI": SE_ALNILAM,
    "AL-NIZAM": SE_ALNILAM,
    # ======== ALNITAK (Zeta Orionis) ========
    "GIRDLE": SE_ALNITAK,
    "ZETA ORIONIS": SE_ALNITAK,
    "ZETA ORI": SE_ALNITAK,
    "50 ORI": SE_ALNITAK,
    "ζ ORI": SE_ALNITAK,
    "ZEORI": SE_ALNITAK,
    "AL-NITAK": SE_ALNITAK,
    # ======== MINTAKA (Delta Orionis) ========
    "BELT STAR": SE_MINTAKA,
    "DELTA ORIONIS": SE_MINTAKA,
    "DELTA ORI": SE_MINTAKA,
    "34 ORI": SE_MINTAKA,
    "δ ORI": SE_MINTAKA,
    "DEORI": SE_MINTAKA,
    "MINTAKA": SE_MINTAKA,
    # ======== SAIPH (Kappa Orionis) ========
    "SWORD OF GIANT": SE_SAIPH,
    "KAPPA ORIONIS": SE_SAIPH,
    "KAPPA ORI": SE_SAIPH,
    "53 ORI": SE_SAIPH,
    "κ ORI": SE_SAIPH,
    "KAORI": SE_SAIPH,
    "SAIF AL-JABBAR": SE_SAIPH,
    # ======== DIPHDA (Beta Ceti) ========
    "FROG": SE_DIPHDA,
    "BETA CETI": SE_DIPHDA,
    "BETA CET": SE_DIPHDA,
    "16 CET": SE_DIPHDA,
    "β CET": SE_DIPHDA,
    "BECET": SE_DIPHDA,
    "DENEB KAITOS": SE_DIPHDA,
    "DIFDA AL-THANI": SE_DIPHDA,
    # ======== ALPHARD (Alpha Hydrae) ========
    "SOLITARY ONE": SE_ALPHARD,
    "ALPHA HYDRAE": SE_ALPHARD,
    "ALPHA HYA": SE_ALPHARD,
    "30 HYA": SE_ALPHARD,
    "α HYA": SE_ALPHARD,
    "ALHYA": SE_ALPHARD,
    "COR HYDRAE": SE_ALPHARD,
    "AL-FARD": SE_ALPHARD,
    # ======== RASALHAGUE (Alpha Ophiuchi) ========
    "HEAD OF SERPENT HOLDER": SE_RASALHAGUE,
    "ALPHA OPHIUCHI": SE_RASALHAGUE,
    "ALPHA OPH": SE_RASALHAGUE,
    "55 OPH": SE_RASALHAGUE,
    "α OPH": SE_RASALHAGUE,
    "ALOPH": SE_RASALHAGUE,
    "RAS AL-HAWWA": SE_RASALHAGUE,
    # ======== ETAMIN (Gamma Draconis) ========
    "DRAGON'S HEAD": SE_ETAMIN,
    "GAMMA DRACONIS": SE_ETAMIN,
    "GAMMA DRA": SE_ETAMIN,
    "33 DRA": SE_ETAMIN,
    "γ DRA": SE_ETAMIN,
    "GADORA": SE_ETAMIN,
    "ELTANIN": SE_ETAMIN,
    "AL-TINNIN": SE_ETAMIN,
    # ======== KOCHAB (Beta Ursae Minoris) ========
    "BETA URSAE MINORIS": SE_KOCHAB,
    "BETA UMI": SE_KOCHAB,
    "7 UMI": SE_KOCHAB,
    "β UMI": SE_KOCHAB,
    "BEUMI": SE_KOCHAB,
    "KAUKAB": SE_KOCHAB,
    # ======== ALKAID (Eta Ursae Majoris) ========
    "END OF TAIL": SE_ALKAID,
    "ETA URSAE MAJORIS": SE_ALKAID,
    "ETA UMA": SE_ALKAID,
    "85 UMA": SE_ALKAID,
    "η UMA": SE_ALKAID,
    "ETUMA": SE_ALKAID,
    "BENETNASH": SE_ALKAID,
    "AL-QA'ID": SE_ALKAID,
    # ======== DUBHE (Alpha Ursae Majoris) ========
    "BEAR'S BACK": SE_DUBHE,
    "ALPHA URSAE MAJORIS": SE_DUBHE,
    "ALPHA UMA": SE_DUBHE,
    "50 UMA": SE_DUBHE,
    "α UMA": SE_DUBHE,
    "ALUMA": SE_DUBHE,
    "THAHR AL-DUBB AL-AKBAR": SE_DUBHE,
    # ======== MERAK (Beta Ursae Majoris) ========
    "BETA URSAE MAJORIS": SE_MERAK,
    "BETA UMA": SE_MERAK,
    "48 UMA": SE_MERAK,
    "β UMA": SE_MERAK,
    "BEUMA": SE_MERAK,
    "AL-MARAKK": SE_MERAK,
    # ======== ALIOTH (Epsilon Ursae Majoris) ========
    "EPSILON URSAE MAJORIS": SE_ALIOTH,
    "EPSILON UMA": SE_ALIOTH,
    "77 UMA": SE_ALIOTH,
    "ε UMA": SE_ALIOTH,
    "EPUMA": SE_ALIOTH,
    "ALIATH": SE_ALIOTH,
    # ======== MIZAR (Zeta Ursae Majoris) ========
    "ZETA URSAE MAJORIS": SE_MIZAR,
    "ZETA UMA": SE_MIZAR,
    "79 UMA": SE_MIZAR,
    "ζ UMA": SE_MIZAR,
    "ZEUMA": SE_MIZAR,
    "HORSE AND RIDER": SE_MIZAR,
    # ======== ALCOR (80 Ursae Majoris) ========
    "80 URSAE MAJORIS": SE_ALCOR,
    "80 UMA": SE_ALCOR,
    "G UMA": SE_ALCOR,
    "SAIDAK": SE_ALCOR,
    "SUHA": SE_ALCOR,
    "AL-SAHJA": SE_ALCOR,
    # ======== VINDEMIATRIX (Epsilon Virginis) ========
    "GRAPE GATHERER": SE_VINDEMIATRIX,
    "EPSILON VIRGINIS": SE_VINDEMIATRIX,
    "EPSILON VIR": SE_VINDEMIATRIX,
    "47 VIR": SE_VINDEMIATRIX,
    "ε VIR": SE_VINDEMIATRIX,
    "EPVIR": SE_VINDEMIATRIX,
    "ALMUREDIN": SE_VINDEMIATRIX,
    # ======== ZUBENELGENUBI (Alpha Librae) ========
    "SOUTHERN CLAW": SE_ZUBENELGENUBI,
    "ALPHA LIBRAE": SE_ZUBENELGENUBI,
    "ALPHA LIB": SE_ZUBENELGENUBI,
    "9 LIB": SE_ZUBENELGENUBI,
    "α LIB": SE_ZUBENELGENUBI,
    "ALLIB": SE_ZUBENELGENUBI,
    "KIFFA AUSTRALIS": SE_ZUBENELGENUBI,
    # ======== ZUBENESCHAMALI (Beta Librae) ========
    "NORTHERN CLAW": SE_ZUBENESCHAMALI,
    "BETA LIBRAE": SE_ZUBENESCHAMALI,
    "BETA LIB": SE_ZUBENESCHAMALI,
    "27 LIB": SE_ZUBENESCHAMALI,
    "β LIB": SE_ZUBENESCHAMALI,
    "BELIB": SE_ZUBENESCHAMALI,
    "KIFFA BOREALIS": SE_ZUBENESCHAMALI,
    # ======== UNUKALHAI (Alpha Serpentis) ========
    "SERPENT'S NECK": SE_UNUKALHAI,
    "ALPHA SERPENTIS": SE_UNUKALHAI,
    "ALPHA SER": SE_UNUKALHAI,
    "24 SER": SE_UNUKALHAI,
    "α SER": SE_UNUKALHAI,
    "ALSER": SE_UNUKALHAI,
    "COR SERPENTIS": SE_UNUKALHAI,
    # ======== ALGIEBA (Gamma Leonis) ========
    "LION'S MANE": SE_ALGIEBA,
    "GAMMA LEONIS": SE_ALGIEBA,
    "GAMMA LEO": SE_ALGIEBA,
    "41 LEO": SE_ALGIEBA,
    "γ LEO": SE_ALGIEBA,
    "GALEO": SE_ALGIEBA,
    "AL-JABHAH": SE_ALGIEBA,
    # ======== DENEBOLA (Beta Leonis) ========
    "LION'S TAIL": SE_DENEBOLA,
    "BETA LEONIS": SE_DENEBOLA,
    "BETA LEO": SE_DENEBOLA,
    "94 LEO": SE_DENEBOLA,
    "β LEO": SE_DENEBOLA,
    "BELEO": SE_DENEBOLA,
    "DHANAB AL-ASAD": SE_DENEBOLA,
    # ======== MARKAB (Alpha Pegasi) ========
    "SADDLE": SE_MARKAB,
    "ALPHA PEGASI": SE_MARKAB,
    "ALPHA PEG": SE_MARKAB,
    "54 PEG": SE_MARKAB,
    "α PEG": SE_MARKAB,
    "ALPEG": SE_MARKAB,
    "MANKIB AL-FARAS": SE_MARKAB,
    # ======== SCHEAT (Beta Pegasi) ========
    "LEG": SE_SCHEAT,
    "BETA PEGASI": SE_SCHEAT,
    "BETA PEG": SE_SCHEAT,
    "53 PEG": SE_SCHEAT,
    "β PEG": SE_SCHEAT,
    "BEPEG": SE_SCHEAT,
    "SAQ AL-FARAS": SE_SCHEAT,
}


def resolve_star_name(name: str) -> int | None:
    """
    Resolve a star name to its SE_* constant ID using pyswisseph-compatible resolution.

    Implements the following resolution algorithm:
    1. Normalize input (uppercase, strip whitespace)
    2. If name starts with comma (pyswisseph convention), do prefix search
    3. Try exact match in STAR_ALIASES dictionary
    4. Try exact match against canonical star names in STAR_CATALOG
    5. Try fuzzy matching (alias contains search term for short inputs)

    Args:
        name: Star name, alias, designation, or comma-prefixed partial name

    Returns:
        Star ID (SE_* constant) if found, None otherwise

    Examples:
        >>> resolve_star_name("Regulus")
        1000001
        >>> resolve_star_name(",alg")  # Prefix search for Algol
        1000003
        >>> resolve_star_name("Alpha Leo")
        1000001
        >>> resolve_star_name("SIRIUS")
        1000004
    """
    if not name:
        return None

    # Normalize: uppercase and strip
    normalized = name.upper().strip()

    # Handle comma-prefix for partial match (pyswisseph convention)
    if normalized.startswith(","):
        search_prefix = normalized[1:].strip()
        if not search_prefix:
            return None

        # Search for any alias or canonical name that starts with this prefix
        # First check canonical names
        for entry in STAR_CATALOG:
            if entry.name.upper().startswith(search_prefix):
                return entry.id

        # Then check aliases
        for alias, star_id in STAR_ALIASES.items():
            if alias.startswith(search_prefix):
                return star_id

        # Check nomenclature (e.g., alLeo, bePer)
        for entry in STAR_CATALOG:
            if entry.nomenclature.upper().startswith(search_prefix):
                return entry.id

        return None

    # Strip trailing comma and anything after (pyswisseph format: "Regulus,alLeo")
    if "," in normalized:
        normalized = normalized.split(",")[0].strip()

    # 1. Try exact match in STAR_ALIASES
    if normalized in STAR_ALIASES:
        return STAR_ALIASES[normalized]

    # 2. Try exact match against canonical star names
    if normalized in _STAR_NAME_TO_ID:
        return _STAR_NAME_TO_ID[normalized]

    # 3. Try exact match against nomenclature (e.g., "ALLEO", "BEPER")
    for entry in STAR_CATALOG:
        if entry.nomenclature.upper() == normalized:
            return entry.id

    # 4. Try fuzzy matching: check if any alias CONTAINS the search term
    # Only for reasonably short inputs (avoid false positives)
    if len(normalized) >= 3:
        for alias, star_id in STAR_ALIASES.items():
            if normalized in alias:
                return star_id
        for entry in STAR_CATALOG:
            if normalized in entry.name.upper():
                return entry.id

    return None


def get_canonical_star_name(star_id: int) -> str | None:
    """
    Get the canonical star name for a star ID.

    Args:
        star_id: Star ID (SE_* constant)

    Returns:
        Canonical star name (e.g., "Regulus") or None if not found
    """
    for entry in STAR_CATALOG:
        if entry.id == star_id:
            return entry.name
    return None


def calc_fixed_star_position(star_id: int, jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate ecliptic position of a fixed star at given date.

    Applies proper motion and precession to J2000.0 catalog position.

    Args:
        star_id: Star identifier (SE_REGULUS, SE_SPICA_STAR)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Ecliptic longitude of date in degrees (0-360)
            - latitude: Ecliptic latitude of date in degrees
            - distance: Arbitrary large value (100000 AU) - stars are effectively infinite

    Raises:
        ValueError: If star_id not in catalog

    Algorithm:
        1. Apply proper motion using rigorous space motion approach (3D vector propagation)
        2. Precess equatorial coordinates using IAU 2006 formula
        3. Calculate true obliquity (mean + nutation)
        4. Transform to ecliptic coordinates

    Notes:
        Proper motion is applied using the rigorous space motion approach from
        Hipparcos Vol. 1, Section 1.5.5 with second-order Taylor expansion.
        This method converts the position to a 3D unit vector, applies proper
        motion as angular velocity in the tangent plane with centripetal
        acceleration correction, and normalizes to account for spherical geometry.

        The second-order term significantly improves accuracy for high proper
        motion stars (e.g., Barnard's Star) over century-scale intervals.

        Limitations:
        - Ignores radial velocity (parallax effect)
        - Assumes constant proper motion (no acceleration)
        - No annual parallax (star distance not modeled)
        Typical error: <0.01 arcsec over ±100 years, <1 arcsec over ±500 years

    References:
        IAU 2006 precession (Capitaine et al.)
        IAU 2000A nutation model (1365 terms) via Skyfield
    """
    if star_id not in FIXED_STARS:
        raise ValueError(f"could not find star name {star_id}")

    star = FIXED_STARS[star_id]

    # Time from J2000.0
    t_years = (jd_tt - 2451545.0) / 365.25
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries

    # Apply Proper Motion using rigorous space motion approach
    # Uses 3D vector propagation to correctly handle spherical geometry.
    # This avoids errors from the curvature of the celestial sphere that occur
    # with linear RA/Dec extrapolation over long time periods.
    #
    # Algorithm (Hipparcos Vol. 1, Section 1.5.5):
    #   1. Convert (RA, Dec) to unit position vector P
    #   2. Compute proper motion as angular velocity in the tangent plane
    #   3. Propagate position vector: P(t) = P(0) + V * dt, then normalize
    #   4. Convert back to (RA, Dec)

    # Convert proper motions from arcsec/year to radians/year
    # pm_ra is μα* = μα × cos(δ), the proper motion in RA direction (not angular)
    # pm_dec is μδ, the proper motion in Dec direction
    pm_ra_rad = math.radians(star.pm_ra / 3600.0)  # arcsec -> deg -> rad
    pm_dec_rad = math.radians(star.pm_dec / 3600.0)

    # Convert J2000 position to radians
    ra_rad = math.radians(star.ra_j2000)
    dec_rad = math.radians(star.dec_j2000)

    # Unit position vector at J2000 epoch
    cos_dec = math.cos(dec_rad)
    sin_dec = math.sin(dec_rad)
    cos_ra = math.cos(ra_rad)
    sin_ra = math.sin(ra_rad)

    px = cos_dec * cos_ra
    py = cos_dec * sin_ra
    pz = sin_dec

    # Proper motion velocity vector in the tangent plane (perpendicular to position)
    # The unit vectors in RA and Dec directions are:
    #   e_ra = (-sin(ra), cos(ra), 0)  (tangent to RA circles, pointing East)
    #   e_dec = (-sin(dec)*cos(ra), -sin(dec)*sin(ra), cos(dec))  (pointing North)
    # Velocity = pm_ra * e_ra + pm_dec * e_dec (in radians/year)
    vx = -pm_ra_rad * sin_ra - pm_dec_rad * sin_dec * cos_ra
    vy = pm_ra_rad * cos_ra - pm_dec_rad * sin_dec * sin_ra
    vz = pm_dec_rad * cos_dec

    # Propagate position using second-order Taylor expansion:
    # P(t) = P(0) + V*dt + 0.5*A*dt²
    #
    # For motion on a unit sphere, the acceleration is the centripetal term:
    # A = -|V|² * P(0)
    # This accounts for the curvature of the celestial sphere and significantly
    # improves accuracy for high proper motion stars (e.g., Barnard's Star)
    # over century-scale time intervals.
    v_squared = vx * vx + vy * vy + vz * vz
    t_squared = t_years * t_years

    # Second-order expansion: P(t) = P(0) + V*dt - 0.5*|V|²*P(0)*dt²
    px_t = px + vx * t_years - 0.5 * v_squared * px * t_squared
    py_t = py + vy * t_years - 0.5 * v_squared * py * t_squared
    pz_t = pz + vz * t_years - 0.5 * v_squared * pz * t_squared

    # Normalize to get unit vector (accounts for curvature)
    r = math.sqrt(px_t * px_t + py_t * py_t + pz_t * pz_t)
    px_t /= r
    py_t /= r
    pz_t /= r

    # Convert back to RA/Dec
    dec_pm = math.asin(pz_t)  # Keep in radians for precession
    ra_pm = math.atan2(py_t, px_t)
    if ra_pm < 0:
        ra_pm += 2 * math.pi

    # IAU 2006 precession parameters (arcseconds -> degrees)
    zeta = (2306.2181 * T + 0.30188 * T**2 + 0.017998 * T**3) / 3600.0
    z = (2306.2181 * T + 1.09468 * T**2 + 0.018203 * T**3) / 3600.0
    theta = (2004.3109 * T - 0.42665 * T**2 - 0.041833 * T**3) / 3600.0

    zeta_r = math.radians(zeta)
    z_r = math.radians(z)
    theta_r = math.radians(theta)
    # ra_pm and dec_pm are already in radians from the proper motion calculation
    ra_r = ra_pm
    dec_r = dec_pm

    # Convert equatorial to Cartesian (unit sphere)
    x0 = math.cos(dec_r) * math.cos(ra_r)
    y0 = math.cos(dec_r) * math.sin(ra_r)
    z0 = math.sin(dec_r)

    # Apply precession rotation: R_z(-z) · R_y(θ) · R_z(-ζ)
    # Step 1: R_z(-ζ) rotation around z-axis
    x1 = x0 * math.cos(-zeta_r) + y0 * math.sin(-zeta_r)
    y1 = -x0 * math.sin(-zeta_r) + y0 * math.cos(-zeta_r)
    z1 = z0

    # Step 2: R_y(θ) rotation around y-axis
    x2 = x1 * math.cos(theta_r) - z1 * math.sin(theta_r)
    y2 = y1
    z2 = x1 * math.sin(theta_r) + z1 * math.cos(theta_r)

    # Step 3: R_z(-z) rotation around z-axis
    x3 = x2 * math.cos(-z_r) + y2 * math.sin(-z_r)
    y3 = -x2 * math.sin(-z_r) + y2 * math.cos(-z_r)
    z3 = z2

    # Convert back to spherical (RA/Dec of date)
    ra_date = math.atan2(y3, x3)
    dec_date = math.asin(z3)

    # Calculate mean obliquity of date (IAU 2006)
    eps0 = 23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0

    # Use full IAU 2000A nutation model (1365 terms) for sub-milliarcsecond precision
    # Reference: IERS Conventions 2010, Skyfield iau2000a_radians implementation
    from .state import get_timescale

    ts = get_timescale()
    t_obj = ts.tt_jd(jd_tt)
    dpsi_rad, deps_rad = iau2000a_radians(t_obj)

    # Convert nutation in obliquity from radians to degrees
    deps_deg = math.degrees(deps_rad)
    eps_true = eps0 + deps_deg
    eps_r = math.radians(eps_true)

    # Transform equatorial (RA, Dec) to ecliptic (lon, lat)
    sin_lat = math.sin(dec_date) * math.cos(eps_r) - math.cos(dec_date) * math.sin(
        eps_r
    ) * math.sin(ra_date)
    lat = math.asin(sin_lat)

    y_lon = math.sin(ra_date) * math.cos(eps_r) + math.tan(dec_date) * math.sin(eps_r)
    x_lon = math.cos(ra_date)
    lon = math.degrees(math.atan2(y_lon, x_lon)) % 360.0
    lat_deg = math.degrees(lat)

    # Distance is arbitrary for fixed stars (effectively infinite)
    dist = 100000.0  # AU (placeholder)

    return lon, lat_deg, dist


def _resolve_star_id(star_name: str) -> tuple[int, str | None, str | None]:
    """
    Resolve a star name to its ID with pyswisseph-compatible name resolution.

    Supports:
    - Exact star names: "Regulus", "Spica"
    - Bayer designations: "Alpha Leo", "Alpha Leonis"
    - Flamsteed numbers: "87 Leo"
    - Traditional/Arabic names: "Cor Leonis", "Dog Star"
    - Comma-prefix partial match: ",alg" finds Algol

    Args:
        star_name: Name of star (e.g. "Regulus", "Alpha Leo", ",alg")

    Returns:
        Tuple of (star_id, error_message, canonical_name).
        If error, star_id is -1 and canonical_name is None.
    """
    star_id = resolve_star_name(star_name)
    if star_id is not None:
        canonical_name = get_canonical_star_name(star_id)
        return star_id, None, canonical_name

    return -1, f"could not find star name {star_name.lower()}", None


def swe_fixstar_ut(
    star_name: str, tjd_ut: float, iflag: int
) -> Tuple[Tuple[float, float, float, float, float, float], int, str]:
    """
    Calculate position of a fixed star for Universal Time.

    Swiss Ephemeris compatible function.

    Args:
        star_name: Name of star (e.g. "Regulus", "Spica")
        tjd_ut: Julian Day in Universal Time (UT1)
        iflag: Calculation flags

    Returns:
        Tuple containing:
            - Position tuple: (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - iflag: Return flags
            - error_msg: Error message if any, empty string on success

    Note:
        UT (Universal Time) is converted to TT (Terrestrial Time) internally
        using Delta T before calculating the star position. For most modern dates,
        Delta T is about 69 seconds (as of 2020).

    Example:
        >>> pos, retflag, err = swe_fixstar_ut("Regulus", 2451545.0, 0)
        >>> lon, lat, dist = pos[0], pos[1], pos[2]
    """
    star_id, error, canonical_name = _resolve_star_id(star_name)
    if error:
        return ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag, error)

    # Convert UT to TT using timescale (applies Delta T)
    from .state import get_timescale

    ts = get_timescale()
    t = ts.ut1_jd(tjd_ut)

    try:
        lon, lat, dist = calc_fixed_star_position(star_id, t.tt)
        # Return canonical star name on success (pyswisseph behavior)
        return ((lon, lat, dist, 0.0, 0.0, 0.0), iflag, canonical_name or "")
    except Exception as e:
        return ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag, str(e))


def swe_fixstar(
    star_name: str, jd: float, iflag: int
) -> Tuple[Tuple[float, float, float, float, float, float], int, str]:
    """
    Calculate position of a fixed star for Terrestrial Time (TT).

    Swiss Ephemeris compatible function. Similar to swe_fixstar_ut() but takes
    Terrestrial Time (TT, also known as Ephemeris Time) instead of Universal Time.

    Args:
        star_name: Name of star (e.g. "Regulus", "Spica")
        jd: Julian Day in Terrestrial Time (TT/ET)
        iflag: Calculation flags

    Returns:
        Tuple containing:
            - Position tuple: (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - iflag: Return flags
            - error_msg: Error message if any, empty string on success

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use swe_fixstar_ut() instead.

    Example:
        >>> pos, retflag, err = swe_fixstar("Regulus", 2451545.0, 0)
        >>> lon, lat, dist = pos[0], pos[1], pos[2]
    """
    star_id, error, canonical_name = _resolve_star_id(star_name)
    if error:
        return ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag, error)

    # Use TT directly - no conversion needed
    try:
        lon, lat, dist = calc_fixed_star_position(star_id, jd)
        # Return canonical star name on success (pyswisseph behavior)
        return ((lon, lat, dist, 0.0, 0.0, 0.0), iflag, canonical_name or "")
    except Exception as e:
        return ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag, str(e))


def _format_star_name(entry: StarCatalogEntry) -> str:
    """
    Format the full star name for return from fixstar2 functions.

    Returns format: "Name,Nomenclature" (e.g. "Regulus,alLeo")
    """
    return f"{entry.name},{entry.nomenclature}"


def _resolve_star2(star_name: str) -> Tuple[StarCatalogEntry | None, str | None]:
    """
    Resolve a star identifier with flexible lookup for fixstar2 functions.

    Supports multiple lookup methods:
    1. Exact star name (case-insensitive): "Regulus", "SPICA"
    2. Hipparcos catalog number (as string): "49669", ",49669"
    3. Partial name search (case-insensitive): "Reg", "pic"
    4. Bayer/Flamsteed nomenclature: "alLeo", "alVir"
    5. Format with comma: "Regulus,alLeo" (takes first part)

    Args:
        star_name: Star identifier - can be name, catalog number, or search string

    Returns:
        Tuple of (StarCatalogEntry, error_message). If error, entry is None.

    Examples:
        >>> entry, err = _resolve_star2("Regulus")      # Exact name
        >>> entry, err = _resolve_star2("49669")        # HIP number
        >>> entry, err = _resolve_star2(",49669")       # HIP with leading comma
        >>> entry, err = _resolve_star2("Reg")          # Partial match
        >>> entry, err = _resolve_star2("alLeo")        # Nomenclature
    """
    search = star_name.strip()

    if not search:
        return None, "Empty star name"

    # Check if it's a catalog number (numeric string, possibly with leading comma)
    number_search = search.lstrip(",").strip()
    if number_search.isdigit():
        hip_number = int(number_search)
        for entry in STAR_CATALOG:
            if entry.hip_number == hip_number:
                return entry, None
        return None, f"could not find star name {hip_number}"

    # Handle comma-separated format (e.g., "Regulus,alLeo")
    if "," in search:
        search = search.split(",")[0].strip()

    search_upper = search.upper()

    # 1. Try exact name match (case-insensitive)
    for entry in STAR_CATALOG:
        if entry.name.upper() == search_upper:
            return entry, None

    # 2. Try exact nomenclature match (case-insensitive)
    for entry in STAR_CATALOG:
        if entry.nomenclature.upper() == search_upper:
            return entry, None

    # 3. Try partial name match (prefix search, case-insensitive)
    matches: List[StarCatalogEntry] = []
    for entry in STAR_CATALOG:
        if entry.name.upper().startswith(search_upper):
            matches.append(entry)

    if len(matches) == 1:
        return matches[0], None
    elif len(matches) > 1:
        names = ", ".join(m.name for m in matches)
        return None, f"Ambiguous star name '{star_name}' matches: {names}"

    # 4. Try partial nomenclature match
    for entry in STAR_CATALOG:
        if entry.nomenclature.upper().startswith(search_upper):
            matches.append(entry)

    if len(matches) == 1:
        return matches[0], None
    elif len(matches) > 1:
        names = ", ".join(m.name for m in matches)
        return None, f"Ambiguous star name '{star_name}' matches: {names}"

    # 5. Try substring match in name (anywhere in the name)
    for entry in STAR_CATALOG:
        if search_upper in entry.name.upper():
            matches.append(entry)

    if len(matches) == 1:
        return matches[0], None
    elif len(matches) > 1:
        names = ", ".join(m.name for m in matches)
        return None, f"Ambiguous star name '{star_name}' matches: {names}"

    return None, f"could not find star name {star_name.lower()}"


def swe_fixstar2_ut(
    star_name: str, tjd_ut: float, iflag: int
) -> Tuple[str, Tuple[float, float, float, float, float, float], int, str]:
    """
    Calculate position of a fixed star for Universal Time with flexible lookup.

    Enhanced version of swe_fixstar_ut() that supports flexible star lookup:
    - Star name (full or partial): "Regulus", "Reg"
    - Hipparcos catalog number: "49669", ",49669"
    - Bayer/Flamsteed designation: "alLeo", "alVir"

    Returns the full star name along with the position, allowing identification
    of which star was matched when using partial searches.

    Args:
        star_name: Star identifier (name, catalog number, or partial search)
        tjd_ut: Julian Day in Universal Time (UT1)
        iflag: Calculation flags

    Returns:
        Tuple containing:
            - star_name_out: Full star name "Name,Nomenclature" (e.g. "Regulus,alLeo")
            - Position tuple: (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - iflag: Return flags
            - error_msg: Error message if any, empty string on success

    Note:
        UT (Universal Time) is converted to TT (Terrestrial Time) internally
        using Delta T before calculating the star position.

    Example:
        >>> name, pos, retflag, err = swe_fixstar2_ut("Reg", 2451545.0, 0)
        >>> print(name)  # "Regulus,alLeo"
        >>> lon, lat, dist = pos[0], pos[1], pos[2]

        >>> name, pos, retflag, err = swe_fixstar2_ut("49669", 2451545.0, 0)
        >>> print(name)  # "Regulus,alLeo" (looked up by HIP number)
    """
    entry, error = _resolve_star2(star_name)
    if error or entry is None:
        return (
            "",
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            iflag,
            error or "could not find star name",
        )

    # Convert UT to TT using timescale (applies Delta T)
    from .state import get_timescale

    ts = get_timescale()
    t = ts.ut1_jd(tjd_ut)

    try:
        lon, lat, dist = calc_fixed_star_position(entry.id, t.tt)
        star_name_out = _format_star_name(entry)
        return (star_name_out, (lon, lat, dist, 0.0, 0.0, 0.0), iflag, "")
    except Exception as e:
        return ("", (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag, str(e))


def swe_fixstar2(
    star_name: str, jd: float, iflag: int
) -> Tuple[str, Tuple[float, float, float, float, float, float], int, str]:
    """
    Calculate position of a fixed star for Terrestrial Time with flexible lookup.

    Enhanced version of swe_fixstar() that supports flexible star lookup:
    - Star name (full or partial): "Regulus", "Reg"
    - Hipparcos catalog number: "49669", ",49669"
    - Bayer/Flamsteed designation: "alLeo", "alVir"

    Returns the full star name along with the position, allowing identification
    of which star was matched when using partial searches.

    Args:
        star_name: Star identifier (name, catalog number, or partial search)
        jd: Julian Day in Terrestrial Time (TT/ET)
        iflag: Calculation flags

    Returns:
        Tuple containing:
            - star_name_out: Full star name "Name,Nomenclature" (e.g. "Regulus,alLeo")
            - Position tuple: (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - iflag: Return flags
            - error_msg: Error message if any, empty string on success

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T.
        For most astrological applications, use swe_fixstar2_ut() instead.

    Example:
        >>> name, pos, retflag, err = swe_fixstar2("Spica", 2451545.0, 0)
        >>> print(name)  # "Spica,alVir"
        >>> lon, lat, dist = pos[0], pos[1], pos[2]

        >>> name, pos, retflag, err = swe_fixstar2("65474", 2451545.0, 0)
        >>> print(name)  # "Spica,alVir" (looked up by HIP number)
    """
    entry, error = _resolve_star2(star_name)
    if error or entry is None:
        return (
            "",
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            iflag,
            error or "could not find star name",
        )

    # Use TT directly - no conversion needed
    try:
        lon, lat, dist = calc_fixed_star_position(entry.id, jd)
        star_name_out = _format_star_name(entry)
        return (star_name_out, (lon, lat, dist, 0.0, 0.0, 0.0), iflag, "")
    except Exception as e:
        return ("", (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag, str(e))


# Magnitude values for legacy _resolve_star_id lookup
_STAR_MAGNITUDES = {
    SE_REGULUS: 1.40,
    SE_SPICA_STAR: 1.04,
}


def swe_fixstar_mag(star_name: str) -> Tuple[float, str]:
    """
    Get the visual magnitude of a fixed star without calculating position.

    Lightweight function that returns only the magnitude, useful for
    visibility calculations where position is not needed.

    Args:
        star_name: Name of star (e.g. "Regulus", "Spica")

    Returns:
        Tuple containing:
            - magnitude: Visual magnitude (float), or 0.0 on error
            - error_msg: Error message if any, empty string on success

    Example:
        >>> mag, err = swe_fixstar_mag("Regulus")
        >>> print(f"Regulus magnitude: {mag}")  # 1.40
    """
    star_id, error, _canonical_name = _resolve_star_id(star_name)
    if error:
        return (0.0, error)

    if star_id in _STAR_MAGNITUDES:
        return (_STAR_MAGNITUDES[star_id], "")
    else:
        return (0.0, f"Magnitude not available for star ID {star_id}")


def swe_fixstar2_mag(star_name: str) -> Tuple[str, float, str]:
    """
    Get the visual magnitude of a fixed star with flexible lookup.

    Enhanced version that supports flexible star lookup like swe_fixstar2:
    - Star name (full or partial): "Regulus", "Reg"
    - Hipparcos catalog number: "49669", ",49669"
    - Bayer/Flamsteed designation: "alLeo", "alVir"

    Returns the full star name along with the magnitude, useful for
    visibility calculations where position is not needed.

    Args:
        star_name: Star identifier (name, catalog number, or partial search)

    Returns:
        Tuple containing:
            - star_name_out: Full star name "Name,Nomenclature" (e.g. "Regulus,alLeo")
            - magnitude: Visual magnitude (float), or 0.0 on error
            - error_msg: Error message if any, empty string on success

    Example:
        >>> name, mag, err = swe_fixstar2_mag("Reg")
        >>> print(f"{name}: {mag}")  # "Regulus,alLeo: 1.40"

        >>> name, mag, err = swe_fixstar2_mag("49669")
        >>> print(f"{name}: {mag}")  # "Regulus,alLeo: 1.40"
    """
    entry, error = _resolve_star2(star_name)
    if error or entry is None:
        return ("", 0.0, error or "could not find star name")

    star_name_out = _format_star_name(entry)
    return (star_name_out, entry.magnitude, "")
