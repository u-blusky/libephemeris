"""
Swiss Ephemeris API-compatible constants for libephemeris.

This module defines all constants used for planetary calculations, including:
- Planet/Body IDs: Numeric identifiers for celestial bodies
- Calculation Flags: Bitwise flags controlling observation parameters
- Sidereal Modes: Ayanamsha systems for sidereal astrology
- Calendar Systems: Julian vs Gregorian calendar selection
- Eclipse Types: Classification of solar and lunar eclipses

Constants are organized into logical groups for easy navigation.
Values match Swiss Ephemeris v2.x for API compatibility.
"""

from __future__ import annotations

# =============================================================================
# PLANET AND BODY IDENTIFIERS
# =============================================================================

# Special values
SE_ECL_NUT: int = -1  # Nutation and obliquity calculation

# Major planets (traditional + modern)
SE_SUN: int = 0
SE_MOON: int = 1
SE_MERCURY: int = 2
SE_VENUS: int = 3
SE_MARS: int = 4
SE_JUPITER: int = 5
SE_SATURN: int = 6
SE_URANUS: int = 7
SE_NEPTUNE: int = 8
SE_PLUTO: int = 9

# Planet aliases without SE_ prefix for pyswisseph compatibility
SUN: int = SE_SUN
MOON: int = SE_MOON
MERCURY: int = SE_MERCURY
VENUS: int = SE_VENUS
MARS: int = SE_MARS
JUPITER: int = SE_JUPITER
SATURN: int = SE_SATURN
URANUS: int = SE_URANUS
NEPTUNE: int = SE_NEPTUNE
PLUTO: int = SE_PLUTO

# Lunar nodes and apsides
SE_MEAN_NODE: int = 10  # Mean lunar node (Dragon's Head)
SE_TRUE_NODE: int = 11  # True (osculating) lunar node
SE_MEAN_APOG: int = 12  # Mean lunar apogee (Black Moon Lilith)
SE_OSCU_APOG: int = 13  # Osculating (true) lunar apogee

# Earth and centaurs
SE_EARTH: int = 14
SE_CHIRON: int = 15  # Centaur between Saturn and Uranus
SE_PHOLUS: int = 16  # Centaur beyond Saturn

# Main belt asteroids
SE_CERES: int = 17
SE_PALLAS: int = 18
SE_JUNO: int = 19
SE_VESTA: int = 20

# Interpolated lunar apsides
SE_INTP_APOG: int = 21  # Interpolated apogee
SE_INTP_PERG: int = 22  # Interpolated perigee

# Count and offsets
SE_NPLANETS: int = 23  # Total number of standard planet IDs
SE_AST_OFFSET: int = 10000  # Offset for asteroid catalog numbers
SE_VARUNA: int = SE_AST_OFFSET + 20000  # TNO Varuna
SE_FICT_OFFSET: int = 40  # Offset for fictitious bodies
SE_NFICT_ELEM: int = 15  # Number of fictitious elements
SE_COMET_OFFSET: int = 1000  # Offset for comet IDs

# Hamburg School Uranian planets (fictitious bodies)
SE_CUPIDO: int = SE_FICT_OFFSET + 0  # 40 - First Uranian planet (Cupido)
SE_HADES: int = SE_FICT_OFFSET + 1  # 41 - Second Uranian planet (Hades)
SE_ZEUS: int = SE_FICT_OFFSET + 2  # 42 - Third Uranian planet (Zeus)
SE_KRONOS: int = SE_FICT_OFFSET + 3  # 43 - Fourth Uranian planet (Kronos)
SE_APOLLON: int = SE_FICT_OFFSET + 4  # 44 - Fifth Uranian planet (Apollon)
SE_ADMETOS: int = SE_FICT_OFFSET + 5  # 45 - Sixth Uranian planet (Admetos)
SE_VULKANUS: int = SE_FICT_OFFSET + 6  # 46 - Seventh Uranian planet (Vulkanus)
SE_POSEIDON: int = SE_FICT_OFFSET + 7  # 47 - Eighth Uranian planet (Poseidon)

# Pyswisseph-compatible aliases for Uranian planets
CUPIDO: int = SE_CUPIDO
HADES: int = SE_HADES
ZEUS: int = SE_ZEUS
KRONOS: int = SE_KRONOS
APOLLON: int = SE_APOLLON
ADMETOS: int = SE_ADMETOS
VULKANUS: int = SE_VULKANUS
POSEIDON: int = SE_POSEIDON

# Transpluto (Isis) - hypothetical trans-Plutonian planet (section 2.7.2 of Swiss Ephemeris docs)
SE_ISIS: int = SE_FICT_OFFSET + 8  # 48 - Transpluto/Isis
SE_TRANSPLUTO: int = SE_ISIS  # Alias for SE_ISIS
ISIS: int = SE_ISIS  # Pyswisseph-compatible alias
TRANSPLUTO: int = SE_TRANSPLUTO  # Pyswisseph-compatible alias

# Vulcan - hypothetical intramercurial planet (section 2.7.5 of Swiss Ephemeris docs / seorbel.txt #16)
SE_VULCAN: int = SE_FICT_OFFSET + 15  # 55 - Intramercurial hypothetical planet
VULCAN: int = SE_VULCAN  # Pyswisseph-compatible alias

# Waldemath - Dr. Waldemath's hypothetical second moon of Earth (section 2.7.7 of Swiss Ephemeris docs / seorbel.txt #18)
# Note: This is different from Mean Lilith and True Lilith which are lunar apogee points
SE_WALDEMATH: int = SE_FICT_OFFSET + 18  # 58 - Waldemath's hypothetical Dark Moon
WALDEMATH: int = SE_WALDEMATH  # Pyswisseph-compatible alias

# Planet X Leverrier - Leverrier's calculated "Planet X" that led to Neptune's discovery (section 2.7.8 of Swiss Ephemeris docs / seorbel.txt #12)
# This is the orbital elements Leverrier computed in 1846 to predict the position of Neptune.
# Historical note: Leverrier called it "Planet X" before Neptune was discovered at that position.
SE_PLANET_X_LEVERRIER: int = (
    SE_FICT_OFFSET + 11
)  # 51 - Leverrier's Planet X (Neptune prediction)
PLANET_X_LEVERRIER: int = SE_PLANET_X_LEVERRIER  # Pyswisseph-compatible alias

# Planet X Adams - Adams' calculated "Planet X" (independently derived, similar to Leverrier's)
# John Couch Adams independently predicted Neptune's position around the same time as Leverrier.
# This uses Adams' orbital elements from seorbel.txt #13 (SE_FICT_OFFSET + 12 = 52).
SE_PLANET_X_ADAMS: int = (
    SE_FICT_OFFSET + 12
)  # 52 - Adams' Planet X (Neptune prediction)
PLANET_X_ADAMS: int = SE_PLANET_X_ADAMS  # Pyswisseph-compatible alias

# Planet X Lowell - Percival Lowell's predicted "Planet X" that led to Pluto's discovery
# Lowell predicted a trans-Neptunian planet based on perceived perturbations in Uranus's orbit.
# The search based on his predictions eventually led to Clyde Tombaugh's discovery of Pluto in 1930,
# though Pluto was too small to be Lowell's predicted Planet X.
# This uses Lowell's orbital elements from seorbel.txt #14 (SE_FICT_OFFSET + 13 = 53).
# Orbital elements (1915 prediction): a=43.0 AU, e=0.202, i=10°
SE_PLANET_X_LOWELL: int = (
    SE_FICT_OFFSET + 13
)  # 53 - Lowell's Planet X (Pluto prediction)
PLANET_X_LOWELL: int = SE_PLANET_X_LOWELL  # Pyswisseph-compatible alias

# Planet X Pickering - William H. Pickering's predicted "Planet O" (1919)
# Pickering proposed several trans-Neptunian planets (Planet O, P, Q, R, S, T, U).
# His most famous prediction, "Planet O" (1919), had the following orbital elements:
# - Semi-major axis: 51.9 AU
# - Eccentricity: 0.31
# - Inclination: 15°
# - Orbital period: ~373.5 years
# Like Lowell's Planet X, these predictions were based on supposed perturbations in outer
# planet orbits, which later proved to be observational errors.
# This uses seorbel.txt #15 (SE_FICT_OFFSET + 14 = 54).
SE_PLANET_X_PICKERING: int = (
    SE_FICT_OFFSET + 14
)  # 54 - Pickering's Planet O/X prediction
PLANET_X_PICKERING: int = SE_PLANET_X_PICKERING  # Pyswisseph-compatible alias

# White Moon (Selena) - Point opposite to Black Moon Lilith (lunar perigee = apogee + 180°)
# In Swiss Ephemeris, this is calculated as Mean Lilith + 180° (i.e., the mean lunar perigee)
# Note: Some systems use True Lilith + 180° instead; libephemeris supports both via calc_white_moon_position()
SE_WHITE_MOON: int = (
    SE_FICT_OFFSET + 16
)  # 56 - White Moon Selena (opposite to Black Moon Lilith)
WHITE_MOON: int = SE_WHITE_MOON  # Pyswisseph-compatible alias
SE_SELENA: int = SE_WHITE_MOON  # Alias - Selena is another name for White Moon
SELENA: int = SE_WHITE_MOON  # Pyswisseph-compatible alias

# Proserpina - hypothetical trans-Plutonian planet used by some astrologers
# This is a hypothetical body not in the standard Swiss Ephemeris seorbel.txt
# Orbital elements are based on traditional astrological sources
SE_PROSERPINA: int = SE_FICT_OFFSET + 17  # 57 - Hypothetical trans-Plutonian planet
PROSERPINA: int = SE_PROSERPINA  # Pyswisseph-compatible alias

SE_NALL_NAT_POINTS: int = SE_NPLANETS + SE_NFICT_ELEM + SE_AST_OFFSET + SE_COMET_OFFSET

# Trans-Neptunian Objects (TNOs) - Catalog number + offset
SE_ERIS: int = 136199 + SE_AST_OFFSET  # Largest known dwarf planet
SE_SEDNA: int = 90377 + SE_AST_OFFSET  # Detached TNO
SE_HAUMEA: int = 136108 + SE_AST_OFFSET  # Fast-rotating dwarf planet
SE_MAKEMAKE: int = 136472 + SE_AST_OFFSET  # Classical Kuiper belt object
SE_IXION: int = 28978 + SE_AST_OFFSET  # Plutino
SE_ORCUS: int = 90482 + SE_AST_OFFSET  # Plutino, "anti-Pluto"
SE_QUAOAR: int = 50000 + SE_AST_OFFSET  # Classical KBO
SE_NESSUS: int = 7066 + SE_AST_OFFSET  # Centaur, astrologically important
SE_ASBOLUS: int = 8405 + SE_AST_OFFSET  # Centaur, astrologically significant
SE_CHARIKLO: int = 10199 + SE_AST_OFFSET  # Centaur, largest known, has ring system
SE_GONGGONG: int = (
    225088 + SE_AST_OFFSET
)  # TNO, dwarf planet candidate (formerly 2007 OR10)
SE_APOPHIS: int = 99942 + SE_AST_OFFSET  # Near-Earth asteroid, close approach 2029
SE_HYGIEA: int = 10 + SE_AST_OFFSET  # Fourth largest asteroid, dwarf planet candidate
SE_INTERAMNIA: int = 704 + SE_AST_OFFSET  # Fifth largest asteroid, main belt
SE_DAVIDA: int = 511 + SE_AST_OFFSET  # Seventh largest asteroid, main belt
SE_EUROPA_AST: int = (
    52 + SE_AST_OFFSET
)  # 52 Europa (main belt asteroid, not Jupiter's moon)
SE_SYLVIA: int = (
    87 + SE_AST_OFFSET
)  # 87 Sylvia (triple asteroid system with moons Romulus and Remus)
SE_PSYCHE: int = (
    16 + SE_AST_OFFSET
)  # 16 Psyche (metallic M-type asteroid, NASA Psyche mission target)
SE_EROS: int = (
    433 + SE_AST_OFFSET
)  # 433 Eros (near-Earth asteroid, NEAR Shoemaker mission target)
SE_AMOR: int = (
    1221 + SE_AST_OFFSET
)  # 1221 Amor (prototype of the Amor near-Earth asteroid class)
SE_ICARUS: int = (
    1566 + SE_AST_OFFSET
)  # 1566 Icarus (Apollo asteroid, highly eccentric orbit reaches inside Mercury)
SE_TORO: int = 1685 + SE_AST_OFFSET  # 1685 Toro (Apollo asteroid, near-Earth asteroid)
SE_SAPPHO: int = (
    80 + SE_AST_OFFSET
)  # 80 Sappho (main belt asteroid, artistic expression and same-sex love)
SE_PANDORA_AST: int = (
    55 + SE_AST_OFFSET
)  # 55 Pandora (main belt asteroid, distinct from Saturn moon Pandora)
SE_LILITH_AST: int = (
    1181 + SE_AST_OFFSET
)  # 1181 Lilith (main belt asteroid, not to be confused with lunar apogee Lilith)
SE_HIDALGO: int = (
    944 + SE_AST_OFFSET
)  # 944 Hidalgo (Centaur-class asteroid with comet-like orbit, astrological research)
SE_TOUTATIS: int = (
    4179 + SE_AST_OFFSET
)  # 4179 Toutatis (Apollo PHA, radar and spacecraft target, tumbling rotation)
SE_ITOKAWA: int = (
    25143 + SE_AST_OFFSET
)  # 25143 Itokawa (Apollo PHA, Hayabusa sample return mission target)
SE_BENNU: int = (
    101955 + SE_AST_OFFSET
)  # 101955 Bennu (Apollo PHA, OSIRIS-REx sample return mission target)
SE_RYUGU: int = (
    162173 + SE_AST_OFFSET
)  # 162173 Ryugu (Apollo PHA, Hayabusa2 sample return mission target)

# =============================================================================
# NAIF IDS FOR SPK KERNELS
# =============================================================================
# NAIF IDs used in JPL SPK kernels for minor bodies.
# Convention: For numbered asteroids, NAIF ID = asteroid_number + 2000000
# These constants simplify registration of SPK bodies.

NAIF_ASTEROID_OFFSET: int = 2000000  # Add asteroid number to get NAIF ID

# Common minor body NAIF IDs (asteroid_number + 2000000)
NAIF_CERES: int = 2000001  # 1 Ceres
NAIF_PALLAS: int = 2000002  # 2 Pallas
NAIF_JUNO: int = 2000003  # 3 Juno
NAIF_VESTA: int = 2000004  # 4 Vesta
NAIF_CHIRON: int = 2002060  # 2060 Chiron
NAIF_PHOLUS: int = 2005145  # 5145 Pholus
NAIF_ERIS: int = 2136199  # 136199 Eris
NAIF_SEDNA: int = 2090377  # 90377 Sedna
NAIF_HAUMEA: int = 2136108  # 136108 Haumea
NAIF_MAKEMAKE: int = 2136472  # 136472 Makemake
NAIF_IXION: int = 2028978  # 28978 Ixion
NAIF_ORCUS: int = 2090482  # 90482 Orcus
NAIF_QUAOAR: int = 2050000  # 50000 Quaoar
NAIF_NESSUS: int = 2007066  # 7066 Nessus (Centaur)
NAIF_ASBOLUS: int = 2008405  # 8405 Asbolus (Centaur)
NAIF_CHARIKLO: int = 2010199  # 10199 Chariklo (Centaur, largest, has rings)
NAIF_GONGGONG: int = 2225088  # 225088 Gonggong (TNO, dwarf planet candidate)
NAIF_APOPHIS: int = 2099942  # 99942 Apophis (Near-Earth asteroid)
NAIF_HYGIEA: int = (
    2000010  # 10 Hygiea (fourth largest asteroid, dwarf planet candidate)
)
NAIF_EROS: int = 2000433  # 433 Eros (near-Earth asteroid, NEAR Shoemaker mission)

# =============================================================================
# SPK BODY NAME MAPPING
# =============================================================================
# Mapping from libephemeris body IDs (SE_*) to JPL Horizons target designations.
# This is used for automatic SPK downloads. The tuple contains:
#   (horizons_id: str, naif_id: int)
# where horizons_id is the identifier used in JPL Horizons queries (typically
# the asteroid catalog number), and naif_id is the NAIF SPICE ID.

SPK_BODY_NAME_MAP: dict[int, tuple[str, int]] = {
    SE_CHIRON: ("2060", NAIF_CHIRON),  # 2060 Chiron (centaur)
    SE_PHOLUS: ("5145", NAIF_PHOLUS),  # 5145 Pholus (centaur)
    SE_CERES: ("1", NAIF_CERES),  # 1 Ceres (dwarf planet)
    SE_PALLAS: ("2", NAIF_PALLAS),  # 2 Pallas (main belt asteroid)
    SE_JUNO: ("3", NAIF_JUNO),  # 3 Juno (main belt asteroid)
    SE_VESTA: ("4", NAIF_VESTA),  # 4 Vesta (main belt asteroid)
    SE_ERIS: ("136199", NAIF_ERIS),  # 136199 Eris (dwarf planet)
    SE_SEDNA: ("90377", NAIF_SEDNA),  # 90377 Sedna (detached TNO)
    SE_HAUMEA: ("136108", NAIF_HAUMEA),  # 136108 Haumea (dwarf planet)
    SE_MAKEMAKE: ("136472", NAIF_MAKEMAKE),  # 136472 Makemake (dwarf planet)
    SE_IXION: ("28978", NAIF_IXION),  # 28978 Ixion (plutino)
    SE_ORCUS: ("90482", NAIF_ORCUS),  # 90482 Orcus (plutino)
    SE_QUAOAR: ("50000", NAIF_QUAOAR),  # 50000 Quaoar (classical KBO)
    SE_VARUNA: ("20000", 2020000),  # 20000 Varuna (classical KBO)
    SE_NESSUS: ("7066", NAIF_NESSUS),  # 7066 Nessus (centaur)
    SE_ASBOLUS: ("8405", NAIF_ASBOLUS),  # 8405 Asbolus (centaur)
    SE_CHARIKLO: ("10199", NAIF_CHARIKLO),  # 10199 Chariklo (centaur, largest, rings)
    SE_GONGGONG: (
        "225088",
        NAIF_GONGGONG,
    ),  # 225088 Gonggong (TNO, dwarf planet candidate)
    SE_APOPHIS: ("99942", NAIF_APOPHIS),  # 99942 Apophis (Near-Earth asteroid)
    SE_HYGIEA: ("10", NAIF_HYGIEA),  # 10 Hygiea (fourth largest asteroid)
    SE_EROS: ("433", NAIF_EROS),  # 433 Eros (near-Earth asteroid, NEAR Shoemaker)
}


def get_horizons_id(ipl: int) -> str | None:
    """
    Get the JPL Horizons target identifier for a libephemeris body ID.

    Args:
        ipl: libephemeris body ID (e.g., SE_CHIRON, SE_ERIS)

    Returns:
        The Horizons target identifier as a string (e.g., "2060" for Chiron),
        or None if the body is not in the mapping.

    Example:
        >>> from libephemeris.constants import get_horizons_id, SE_CHIRON
        >>> get_horizons_id(SE_CHIRON)
        '2060'
    """
    if ipl in SPK_BODY_NAME_MAP:
        return SPK_BODY_NAME_MAP[ipl][0]
    return None


def get_naif_id_from_ipl(ipl: int) -> int | None:
    """
    Get the NAIF SPICE ID for a libephemeris body ID.

    Args:
        ipl: libephemeris body ID (e.g., SE_CHIRON, SE_ERIS)

    Returns:
        The NAIF ID (e.g., 2002060 for Chiron), or None if the body
        is not in the mapping.

    Example:
        >>> from libephemeris.constants import get_naif_id_from_ipl, SE_CHIRON
        >>> get_naif_id_from_ipl(SE_CHIRON)
        2002060
    """
    if ipl in SPK_BODY_NAME_MAP:
        return SPK_BODY_NAME_MAP[ipl][1]
    return None


def get_spk_body_info_from_map(ipl: int) -> tuple[str, int] | None:
    """
    Get both Horizons ID and NAIF ID for a libephemeris body.

    Args:
        ipl: libephemeris body ID (e.g., SE_CHIRON, SE_ERIS)

    Returns:
        A tuple of (horizons_id, naif_id), or None if the body
        is not in the mapping.

    Example:
        >>> from libephemeris.constants import get_spk_body_info_from_map, SE_ERIS
        >>> get_spk_body_info_from_map(SE_ERIS)
        ('136199', 2136199)
    """
    return SPK_BODY_NAME_MAP.get(ipl)


# =============================================================================
# VIRTUAL POINTS AND CALCULATED POSITIONS
# =============================================================================

# Fixed Stars (high offset to avoid ID collisions)
SE_FIXSTAR_OFFSET: int = 1000000
SE_REGULUS: int = SE_FIXSTAR_OFFSET + 1  # Alpha Leonis
SE_SPICA_STAR: int = SE_FIXSTAR_OFFSET + 2  # Alpha Virginis
SE_ALGOL: int = SE_FIXSTAR_OFFSET + 3  # Beta Persei
SE_SIRIUS: int = SE_FIXSTAR_OFFSET + 4  # Alpha Canis Majoris
SE_ALDEBARAN: int = SE_FIXSTAR_OFFSET + 5  # Alpha Tauri
SE_ANTARES: int = SE_FIXSTAR_OFFSET + 6  # Alpha Scorpii
SE_VEGA: int = SE_FIXSTAR_OFFSET + 7  # Alpha Lyrae
SE_POLARIS: int = SE_FIXSTAR_OFFSET + 8  # Alpha Ursae Minoris
SE_FOMALHAUT: int = SE_FIXSTAR_OFFSET + 9  # Alpha Piscis Austrini
SE_BETELGEUSE: int = SE_FIXSTAR_OFFSET + 10  # Alpha Orionis
SE_RIGEL: int = SE_FIXSTAR_OFFSET + 11  # Beta Orionis
SE_PROCYON: int = SE_FIXSTAR_OFFSET + 12  # Alpha Canis Minoris
SE_CAPELLA: int = SE_FIXSTAR_OFFSET + 13  # Alpha Aurigae
SE_ARCTURUS: int = SE_FIXSTAR_OFFSET + 14  # Alpha Bootis
SE_DENEB: int = SE_FIXSTAR_OFFSET + 15  # Alpha Cygni
SE_POLLUX: int = SE_FIXSTAR_OFFSET + 16  # Beta Geminorum
SE_CASTOR: int = SE_FIXSTAR_OFFSET + 17  # Alpha Geminorum
SE_ALTAIR: int = SE_FIXSTAR_OFFSET + 18  # Alpha Aquilae
SE_ACHERNAR: int = SE_FIXSTAR_OFFSET + 19  # Alpha Eridani
SE_CANOPUS: int = SE_FIXSTAR_OFFSET + 20  # Alpha Carinae
SE_ACRUX: int = SE_FIXSTAR_OFFSET + 21  # Alpha Crucis
SE_MIMOSA: int = SE_FIXSTAR_OFFSET + 22  # Beta Crucis
SE_GACRUX: int = SE_FIXSTAR_OFFSET + 23  # Gamma Crucis
SE_HADAR: int = SE_FIXSTAR_OFFSET + 24  # Beta Centauri
SE_RIGIL_KENT: int = SE_FIXSTAR_OFFSET + 25  # Alpha Centauri
SE_SHAULA: int = SE_FIXSTAR_OFFSET + 26  # Lambda Scorpii
SE_BELLATRIX: int = SE_FIXSTAR_OFFSET + 27  # Gamma Orionis
SE_ELNATH: int = SE_FIXSTAR_OFFSET + 28  # Beta Tauri
SE_MIRA: int = SE_FIXSTAR_OFFSET + 29  # Omicron Ceti
SE_ALNILAM: int = SE_FIXSTAR_OFFSET + 30  # Epsilon Orionis
SE_ALNITAK: int = SE_FIXSTAR_OFFSET + 31  # Zeta Orionis
SE_MINTAKA: int = SE_FIXSTAR_OFFSET + 32  # Delta Orionis
SE_SAIPH: int = SE_FIXSTAR_OFFSET + 33  # Kappa Orionis
SE_DIPHDA: int = SE_FIXSTAR_OFFSET + 34  # Beta Ceti
SE_ALPHARD: int = SE_FIXSTAR_OFFSET + 35  # Alpha Hydrae
SE_RASALHAGUE: int = SE_FIXSTAR_OFFSET + 36  # Alpha Ophiuchi
SE_ETAMIN: int = SE_FIXSTAR_OFFSET + 37  # Gamma Draconis
SE_KOCHAB: int = SE_FIXSTAR_OFFSET + 38  # Beta Ursae Minoris
SE_ALKAID: int = SE_FIXSTAR_OFFSET + 39  # Eta Ursae Majoris
SE_DUBHE: int = SE_FIXSTAR_OFFSET + 40  # Alpha Ursae Majoris
SE_MERAK: int = SE_FIXSTAR_OFFSET + 41  # Beta Ursae Majoris
SE_ALIOTH: int = SE_FIXSTAR_OFFSET + 42  # Epsilon Ursae Majoris
SE_MIZAR: int = SE_FIXSTAR_OFFSET + 43  # Zeta Ursae Majoris
SE_ALCOR: int = SE_FIXSTAR_OFFSET + 44  # 80 Ursae Majoris
SE_VINDEMIATRIX: int = SE_FIXSTAR_OFFSET + 45  # Epsilon Virginis
SE_ZUBENELGENUBI: int = SE_FIXSTAR_OFFSET + 46  # Alpha Librae
SE_ZUBENESCHAMALI: int = SE_FIXSTAR_OFFSET + 47  # Beta Librae
SE_UNUKALHAI: int = SE_FIXSTAR_OFFSET + 48  # Alpha Serpentis
SE_ALGIEBA: int = SE_FIXSTAR_OFFSET + 49  # Gamma Leonis
SE_DENEBOLA: int = SE_FIXSTAR_OFFSET + 50  # Beta Leonis
SE_MARKAB: int = SE_FIXSTAR_OFFSET + 51  # Alpha Pegasi
SE_SCHEAT: int = SE_FIXSTAR_OFFSET + 52  # Beta Pegasi
SE_ALCYONE: int = SE_FIXSTAR_OFFSET + 53  # Eta Tauri (Pleiades) - Behenian star
SE_ALGORAB: int = SE_FIXSTAR_OFFSET + 54  # Delta Corvi - Behenian star
SE_ALPHECCA: int = SE_FIXSTAR_OFFSET + 55  # Alpha Coronae Borealis - Behenian star
SE_DENEB_ALGEDI: int = SE_FIXSTAR_OFFSET + 56  # Delta Capricorni - Behenian star

# Pleiades cluster stars (visible members)
SE_ASTEROPE: int = SE_FIXSTAR_OFFSET + 57  # 21 Tauri - Pleiades member
SE_CELAENO: int = SE_FIXSTAR_OFFSET + 58  # 16 Tauri - Pleiades member
SE_ELECTRA: int = SE_FIXSTAR_OFFSET + 59  # 17 Tauri - Pleiades member
SE_MAIA: int = SE_FIXSTAR_OFFSET + 60  # 20 Tauri - Pleiades member
SE_MEROPE: int = SE_FIXSTAR_OFFSET + 61  # 23 Tauri - Pleiades member
SE_TAYGETA: int = SE_FIXSTAR_OFFSET + 62  # 19 Tauri - Pleiades member
SE_ATLAS: int = SE_FIXSTAR_OFFSET + 63  # 27 Tauri - Pleiades member
SE_PLEIONE: int = SE_FIXSTAR_OFFSET + 64  # 28 Tauri - Pleiades member

# Hyades cluster stars (visible members in Taurus)
SE_PRIMA_HYADUM: int = SE_FIXSTAR_OFFSET + 65  # Gamma Tauri - Hyades member
SE_SECUNDA_HYADUM: int = SE_FIXSTAR_OFFSET + 66  # Delta^1 Tauri - Hyades member
SE_THETA_TAURI: int = SE_FIXSTAR_OFFSET + 67  # Theta^2 Tauri - Hyades member
SE_AIN: int = SE_FIXSTAR_OFFSET + 68  # Epsilon Tauri - Hyades member

# Orion constellation - completing the major stars
SE_MEISSA: int = SE_FIXSTAR_OFFSET + 69  # Lambda Orionis - Orion's head

# Ursa Major (Big Dipper) - completing the asterism
SE_PHECDA: int = SE_FIXSTAR_OFFSET + 70  # Gamma Ursae Majoris - bowl star
SE_MEGREZ: int = SE_FIXSTAR_OFFSET + 71  # Delta Ursae Majoris - bowl-handle junction

# Crux (Southern Cross) - completing the constellation
SE_DELTA_CRUCIS: int = (
    SE_FIXSTAR_OFFSET + 72
)  # Delta Crucis - fourth star of Southern Cross

# Centaurus - completing the bright stars of the constellation
SE_MENKENT: int = SE_FIXSTAR_OFFSET + 73  # Theta Centauri
SE_MUHLIFAIN: int = SE_FIXSTAR_OFFSET + 74  # Gamma Centauri
SE_EPSILON_CENTAURI: int = SE_FIXSTAR_OFFSET + 75  # Epsilon Centauri
SE_ETA_CENTAURI: int = SE_FIXSTAR_OFFSET + 76  # Eta Centauri
SE_ZETA_CENTAURI: int = SE_FIXSTAR_OFFSET + 77  # Zeta Centauri

# Scorpius constellation - completing the bright stars
SE_SARGAS: int = SE_FIXSTAR_OFFSET + 78  # Theta Scorpii
SE_DSCHUBBA: int = SE_FIXSTAR_OFFSET + 79  # Delta Scorpii
SE_GRAFFIAS: int = SE_FIXSTAR_OFFSET + 80  # Beta Scorpii
SE_LESATH: int = SE_FIXSTAR_OFFSET + 81  # Upsilon Scorpii

# Leo constellation - completing the bright stars
SE_ZOSMA: int = SE_FIXSTAR_OFFSET + 82  # Delta Leonis

# ======== ZODIACAL CONSTELLATION BRIGHT STARS ========
# Stars from zodiacal constellations used in astrological interpretation

# Aries constellation - the Ram
SE_HAMAL: int = SE_FIXSTAR_OFFSET + 83  # Alpha Arietis - brightest in Aries
SE_SHERATAN: int = SE_FIXSTAR_OFFSET + 84  # Beta Arietis
SE_MESARTHIM: int = SE_FIXSTAR_OFFSET + 85  # Gamma Arietis

# Cancer constellation - the Crab
SE_ACUBENS: int = SE_FIXSTAR_OFFSET + 86  # Alpha Cancri
SE_TARF: int = SE_FIXSTAR_OFFSET + 87  # Beta Cancri - brightest in Cancer
SE_ASELLUS_BOREALIS: int = SE_FIXSTAR_OFFSET + 88  # Gamma Cancri - Northern Donkey
SE_ASELLUS_AUSTRALIS: int = SE_FIXSTAR_OFFSET + 89  # Delta Cancri - Southern Donkey

# Sagittarius constellation - the Archer
SE_KAUS_AUSTRALIS: int = SE_FIXSTAR_OFFSET + 90  # Epsilon Sagittarii - brightest
SE_NUNKI: int = SE_FIXSTAR_OFFSET + 91  # Sigma Sagittarii
SE_KAUS_MEDIA: int = SE_FIXSTAR_OFFSET + 92  # Delta Sagittarii
SE_KAUS_BOREALIS: int = SE_FIXSTAR_OFFSET + 93  # Lambda Sagittarii
SE_ASCELLA: int = SE_FIXSTAR_OFFSET + 94  # Zeta Sagittarii

# Capricornus constellation - the Sea Goat (complementing Deneb Algedi)
SE_ALGEDI: int = SE_FIXSTAR_OFFSET + 95  # Alpha Capricorni - the Goat
SE_DABIH: int = SE_FIXSTAR_OFFSET + 96  # Beta Capricorni
SE_NASHIRA: int = SE_FIXSTAR_OFFSET + 97  # Gamma Capricorni - the Fortunate One

# Aquarius constellation - the Water Bearer
SE_SADALSUUD: int = SE_FIXSTAR_OFFSET + 98  # Beta Aquarii - brightest in Aquarius
SE_SADALMELIK: int = SE_FIXSTAR_OFFSET + 99  # Alpha Aquarii
SE_SKAT: int = SE_FIXSTAR_OFFSET + 100  # Delta Aquarii

# Pisces constellation - the Fishes
SE_ETA_PISCIUM: int = SE_FIXSTAR_OFFSET + 101  # Eta Piscium - brightest in Pisces
SE_ALRESCHA: int = SE_FIXSTAR_OFFSET + 102  # Alpha Piscium - the Knot

# Astrological Angles (requires observer location)
SE_ANGLE_OFFSET: int = 2000000
SE_ASCENDANT: int = SE_ANGLE_OFFSET + 1  # Rising sign/degree
SE_MC: int = SE_ANGLE_OFFSET + 2  # Medium Coeli (Midheaven)
SE_DESCENDANT: int = SE_ANGLE_OFFSET + 3  # Setting point (Asc + 180°)
SE_IC: int = SE_ANGLE_OFFSET + 4  # Imum Coeli (MC + 180°)
SE_VERTEX: int = SE_ANGLE_OFFSET + 5  # Western intersection of prime vertical
SE_ANTIVERTEX: int = SE_ANGLE_OFFSET + 6  # Eastern intersection (Vertex + 180°)

# Arabic Parts (Lots) - Require pre-calculated planetary positions
SE_ARABIC_OFFSET: int = 3000000
SE_PARS_FORTUNAE: int = SE_ARABIC_OFFSET + 1  # Part of Fortune
SE_PARS_SPIRITUS: int = SE_ARABIC_OFFSET + 2  # Part of Spirit
SE_PARS_AMORIS: int = SE_ARABIC_OFFSET + 3  # Part of Eros/Love
SE_PARS_FIDEI: int = SE_ARABIC_OFFSET + 4  # Part of Faith

# =============================================================================
# CALCULATION FLAGS
# =============================================================================
# Ephemeris selection (currently only SWIEPH/JPL mode supported)
SEFLG_JPLEPH: int = 1  # Use JPL ephemeris
SEFLG_SWIEPH: int = 2  # Use Swiss Ephemeris (libephemeris uses Skyfield/JPL)
SEFLG_MOSEPH: int = 4  # Use Moshier ephemeris (not supported)

# Observer location and reference frame
SEFLG_HELCTR: int = 8  # Heliocentric position
SEFLG_TRUEPOS: int = 16  # True geometric position (no light time)
SEFLG_J2000: int = 32  # J2000.0 reference frame
SEFLG_NONUT: int = 64  # No nutation
SEFLG_SPEED3: int = 128  # High precision speed (3 calls)
SEFLG_SPEED: int = 256  # Calculate velocity
SEFLG_NOGDEFL: int = 512  # No gravitational deflection
SEFLG_NOABERR: int = 1024  # No aberration
SEFLG_ASTROMETRIC: int = SEFLG_NOABERR | SEFLG_NOGDEFL  # Astrometric position
SEFLG_EQUATORIAL: int = 2048  # Equatorial coordinates (RA/Dec)
SEFLG_XYZ: int = 4096  # Cartesian coordinates
SEFLG_RADIANS: int = 8192  # Return angles in radians
SEFLG_BARYCTR: int = 16384  # Barycentric position
SEFLG_TOPOCTR: int = 32768  # Topocentric position (requires swe_set_topo)
SEFLG_SIDEREAL: int = 65536  # Sidereal positions
SEFLG_ICRS: int = 131072  # ICRS reference frame

# =============================================================================
# PYSWISSEPH-COMPATIBLE FLAG ALIASES (FLG_* instead of SEFLG_*)
# =============================================================================
# pyswisseph uses FLG_* prefix while Swiss Ephemeris C library uses SEFLG_*
# These aliases provide full API compatibility with pyswisseph

FLG_JPLEPH: int = SEFLG_JPLEPH
FLG_SWIEPH: int = SEFLG_SWIEPH
FLG_MOSEPH: int = SEFLG_MOSEPH
FLG_HELCTR: int = SEFLG_HELCTR
FLG_TRUEPOS: int = SEFLG_TRUEPOS
FLG_J2000: int = SEFLG_J2000
FLG_NONUT: int = SEFLG_NONUT
FLG_SPEED3: int = SEFLG_SPEED3
FLG_SPEED: int = SEFLG_SPEED
FLG_NOGDEFL: int = SEFLG_NOGDEFL
FLG_NOABERR: int = SEFLG_NOABERR
FLG_ASTROMETRIC: int = SEFLG_ASTROMETRIC
FLG_EQUATORIAL: int = SEFLG_EQUATORIAL  # Equatorial coordinates (RA/Dec)
FLG_XYZ: int = SEFLG_XYZ
FLG_RADIANS: int = SEFLG_RADIANS
FLG_BARYCTR: int = SEFLG_BARYCTR
FLG_TOPOCTR: int = SEFLG_TOPOCTR
FLG_SIDEREAL: int = SEFLG_SIDEREAL
FLG_ICRS: int = SEFLG_ICRS

# Other aliases
AST_OFFSET: int = SE_AST_OFFSET

# =============================================================================
# SIDEREAL (AYANAMSHA) MODES
# =============================================================================

# Western sidereal traditions
SE_SIDM_FAGAN_BRADLEY: int = 0  # Fagan-Bradley (Synetic Vernal Point)
SE_SIDM_LAHIRI: int = 1  # Lahiri (Chitrapaksha, Indian standard)
SE_SIDM_DELUCE: int = 2  # De Luce
SE_SIDM_RAMAN: int = 3  # B.V. Raman
SE_SIDM_USHASHASHI: int = 4  # Ushashashi
SE_SIDM_KRISHNAMURTI: int = 5  # K.S. Krishnamurti (KP)
SE_SIDM_DJWHAL_KHUL: int = 6  # Djwhal Khul (Alice Bailey)
SE_SIDM_YUKTESHWAR: int = 7  # Yukteshwar
SE_SIDM_JN_BHASIN: int = 8  # J.N. Bhasin

# Babylonian traditions
SE_SIDM_BABYL_KUGLER1: int = 9  # Kugler variant 1
SE_SIDM_BABYL_KUGLER2: int = 10  # Kugler variant 2
SE_SIDM_BABYL_KUGLER3: int = 11  # Kugler variant 3
SE_SIDM_BABYL_HUBER: int = 12  # Huber
SE_SIDM_BABYL_ETPSC: int = 13  # ETPSC
SE_SIDM_BABYL_BRITTON: int = 38  # Britton

# Star-based ayanamshas
SE_SIDM_ALDEBARAN_15TAU: int = 14  # Aldebaran at 15° Taurus
SE_SIDM_TRUE_CITRA: int = 27  # True position of Spica (180° Citra)
SE_SIDM_TRUE_REVATI: int = 28  # True position of Revati
SE_SIDM_TRUE_PUSHYA: int = 29  # True position of Pushya (Cancri)
SE_SIDM_TRUE_MULA: int = 35  # True position of Mula (λ Scorpii)
SE_SIDM_TRUE_SHEORAN: int = 39  # True Sheoran

# Historical epochs
SE_SIDM_HIPPARCHOS: int = 15  # Hipparchos (128 BC)
SE_SIDM_SASSANIAN: int = 16  # Sassanian
SE_SIDM_J2000: int = 18  # J2000.0 (no ayanamsha)
SE_SIDM_J1900: int = 19  # J1900.0
SE_SIDM_B1950: int = 20  # B1950.0

# Indian traditions
SE_SIDM_SURYASIDDHANTA: int = 21  # Suryasiddhanta
SE_SIDM_SURYASIDDHANTA_MSUN: int = 22  # Suryasiddhanta (mean Sun)
SE_SIDM_ARYABHATA: int = 23  # Aryabhata
SE_SIDM_ARYABHATA_MSUN: int = 24  # Aryabhata (mean Sun)
SE_SIDM_ARYABHATA_522: int = 37  # Aryabhata 522
SE_SIDM_SS_REVATI: int = 25  # Suryasiddhanta Revati
SE_SIDM_SS_CITRA: int = 26  # Suryasiddhanta Citra

# Galactic alignment systems
SE_SIDM_GALCENT_0SAG: int = 17  # Galactic Center at 0° Sagittarius
SE_SIDM_GALCENT_RGILBRAND: int = 30  # Galactic Center (Gil Brand)
SE_SIDM_GALCENT_MULA_WILHELM: int = 36  # Galactic Center at Mula (Wilhelm)
SE_SIDM_GALCENT_COCHRANE: int = 40  # Galactic Center (Cochrane)
SE_SIDM_GALEQU_IAU1958: int = 31  # Galactic Equator (IAU 1958)
SE_SIDM_GALEQU_TRUE: int = 32  # Galactic Equator (True)
SE_SIDM_GALEQU_MULA: int = 33  # Galactic Equator at Mula
SE_SIDM_GALEQU_FIORENZA: int = 41  # Galactic Equator (Fiorenza)
SE_SIDM_GALALIGN_MARDYKS: int = 34  # Galactic Alignment (Mardyks)

# Other systems
SE_SIDM_VALENS_MOON: int = 42  # Vettius Valens (Moon-based)
SE_SIDM_USER: int = 255  # User-defined ayanamsha

# =============================================================================
# PYSWISSEPH-COMPATIBLE SIDEREAL MODE ALIASES (SIDM_* instead of SE_SIDM_*)
# =============================================================================
# pyswisseph uses SIDM_* prefix while Swiss Ephemeris C library uses SE_SIDM_*

# Western sidereal traditions
SIDM_FAGAN_BRADLEY: int = SE_SIDM_FAGAN_BRADLEY
SIDM_LAHIRI: int = SE_SIDM_LAHIRI
SIDM_DELUCE: int = SE_SIDM_DELUCE
SIDM_RAMAN: int = SE_SIDM_RAMAN
SIDM_USHASHASHI: int = SE_SIDM_USHASHASHI
SIDM_KRISHNAMURTI: int = SE_SIDM_KRISHNAMURTI
SIDM_DJWHAL_KHUL: int = SE_SIDM_DJWHAL_KHUL
SIDM_YUKTESHWAR: int = SE_SIDM_YUKTESHWAR
SIDM_JN_BHASIN: int = SE_SIDM_JN_BHASIN

# Babylonian traditions
SIDM_BABYL_KUGLER1: int = SE_SIDM_BABYL_KUGLER1
SIDM_BABYL_KUGLER2: int = SE_SIDM_BABYL_KUGLER2
SIDM_BABYL_KUGLER3: int = SE_SIDM_BABYL_KUGLER3
SIDM_BABYL_HUBER: int = SE_SIDM_BABYL_HUBER
SIDM_BABYL_ETPSC: int = SE_SIDM_BABYL_ETPSC
SIDM_BABYL_BRITTON: int = SE_SIDM_BABYL_BRITTON

# Star-based ayanamshas
SIDM_ALDEBARAN_15TAU: int = SE_SIDM_ALDEBARAN_15TAU
SIDM_TRUE_CITRA: int = SE_SIDM_TRUE_CITRA
SIDM_TRUE_REVATI: int = SE_SIDM_TRUE_REVATI
SIDM_TRUE_PUSHYA: int = SE_SIDM_TRUE_PUSHYA
SIDM_TRUE_MULA: int = SE_SIDM_TRUE_MULA
SIDM_TRUE_SHEORAN: int = SE_SIDM_TRUE_SHEORAN

# Historical epochs
SIDM_HIPPARCHOS: int = SE_SIDM_HIPPARCHOS
SIDM_SASSANIAN: int = SE_SIDM_SASSANIAN
SIDM_J2000: int = SE_SIDM_J2000
SIDM_J1900: int = SE_SIDM_J1900
SIDM_B1950: int = SE_SIDM_B1950

# Indian traditions
SIDM_SURYASIDDHANTA: int = SE_SIDM_SURYASIDDHANTA
SIDM_SURYASIDDHANTA_MSUN: int = SE_SIDM_SURYASIDDHANTA_MSUN
SIDM_ARYABHATA: int = SE_SIDM_ARYABHATA
SIDM_ARYABHATA_MSUN: int = SE_SIDM_ARYABHATA_MSUN
SIDM_ARYABHATA_522: int = SE_SIDM_ARYABHATA_522
SIDM_SS_REVATI: int = SE_SIDM_SS_REVATI
SIDM_SS_CITRA: int = SE_SIDM_SS_CITRA

# Galactic alignment systems
SIDM_GALCENT_0SAG: int = SE_SIDM_GALCENT_0SAG
SIDM_GALCENT_RGILBRAND: int = SE_SIDM_GALCENT_RGILBRAND
SIDM_GALCENT_MULA_WILHELM: int = SE_SIDM_GALCENT_MULA_WILHELM
SIDM_GALCENT_COCHRANE: int = SE_SIDM_GALCENT_COCHRANE
SIDM_GALEQU_IAU1958: int = SE_SIDM_GALEQU_IAU1958
SIDM_GALEQU_TRUE: int = SE_SIDM_GALEQU_TRUE
SIDM_GALEQU_MULA: int = SE_SIDM_GALEQU_MULA
SIDM_GALEQU_FIORENZA: int = SE_SIDM_GALEQU_FIORENZA
SIDM_GALALIGN_MARDYKS: int = SE_SIDM_GALALIGN_MARDYKS

# Other systems
SIDM_VALENS_MOON: int = SE_SIDM_VALENS_MOON
SIDM_USER: int = SE_SIDM_USER


# =============================================================================
# CALENDAR SYSTEMS
# =============================================================================

SE_JUL_CAL: int = 0  # Julian calendar
SE_GREG_CAL: int = 1  # Gregorian calendar

# =============================================================================
# ECLIPSE TYPES AND FLAGS
# =============================================================================
# Eclipse geometry types
SE_ECL_CENTRAL: int = 1  # Central eclipse (Moon's shadow axis crosses Earth)
SE_ECL_NONCENTRAL: int = 2  # Non-central eclipse
SE_ECL_TOTAL: int = 4  # Total eclipse (Sun/Earth completely covered)
SE_ECL_ANNULAR: int = 8  # Annular eclipse (ring of fire)
SE_ECL_PARTIAL: int = 16  # Partial eclipse
SE_ECL_ANNULAR_TOTAL: int = (
    32  # Hybrid eclipse (annular at some locations, total at others)
)
SE_ECL_PENUMBRAL: int = 64  # Penumbral lunar eclipse
SE_ECL_GRAZING: int = 65536  # Grazing occultation (target passes near lunar limb)

# Composite eclipse type masks
SE_ECL_ALLTYPES_SOLAR: int = (
    SE_ECL_CENTRAL
    | SE_ECL_NONCENTRAL
    | SE_ECL_TOTAL
    | SE_ECL_ANNULAR
    | SE_ECL_PARTIAL
    | SE_ECL_ANNULAR_TOTAL
)
SE_ECL_ALLTYPES_LUNAR: int = SE_ECL_TOTAL | SE_ECL_PARTIAL | SE_ECL_PENUMBRAL

# Eclipse visibility and contact flags
SE_ECL_VISIBLE: int = 128  # Eclipse visible at location
SE_ECL_MAX_VISIBLE: int = 256  # Maximum phase visible
SE_ECL_1ST_VISIBLE: int = 512  # First contact visible
SE_ECL_2ND_VISIBLE: int = 1024  # Second contact visible
SE_ECL_3RD_VISIBLE: int = 2048  # Third contact visible
SE_ECL_4TH_VISIBLE: int = 4096  # Fourth contact visible
SE_ECL_ONE_TRY: int = 32768  # Try only once (optimization flag)

# pyswisseph-compatible alias for grazing occultation
ECL_GRAZING: int = SE_ECL_GRAZING

# ============================================================================
# NODAL/APSIDAL CALCULATION METHOD FLAGS
# =============================================================================
# Method flags for nod_aps() and nod_aps_ut() functions

SE_NODBIT_MEAN: int = 1  # Mean orbital elements (averaged over perturbations)
SE_NODBIT_OSCU: int = 2  # Osculating elements (instantaneous/perturbed)
SE_NODBIT_OSCU_BAR: int = 4  # Barycentric osculating elements
SE_NODBIT_FOPOINT: int = 256  # Focal point (second focus of ellipse)

# pyswisseph-compatible aliases (without SE_ prefix)
NODBIT_MEAN: int = SE_NODBIT_MEAN
NODBIT_OSCU: int = SE_NODBIT_OSCU
NODBIT_OSCU_BAR: int = SE_NODBIT_OSCU_BAR
NODBIT_FOPOINT: int = SE_NODBIT_FOPOINT

# =============================================================================
# RISE/SET/TRANSIT CALCULATION FLAGS
# =============================================================================
# Event type flags for rise_trans() rsmi parameter

SE_CALC_RISE: int = 1  # Calculate rise time
SE_CALC_SET: int = 2  # Calculate set time
SE_CALC_MTRANSIT: int = 4  # Calculate meridian transit (upper culmination)
SE_CALC_ITRANSIT: int = 8  # Calculate lower transit (anti-culmination)

# Bitmask flags for additional rise/set options
SE_BIT_DISC_CENTER: int = 256  # Use center of disc instead of upper limb
SE_BIT_DISC_BOTTOM: int = 8192  # Use lower limb of disc
SE_BIT_NO_REFRACTION: int = 512  # No atmospheric refraction
SE_BIT_CIVIL_TWILIGHT: int = 1024  # Civil twilight (Sun at -6 degrees)
SE_BIT_NAUTIC_TWILIGHT: int = 2048  # Nautical twilight (Sun at -12 degrees)
SE_BIT_ASTRO_TWILIGHT: int = 4096  # Astronomical twilight (Sun at -18 degrees)
SE_BIT_FIXED_DISC_SIZE: int = 16384  # Use fixed disc size (ignore parallax)

# pyswisseph-compatible aliases (without SE_ prefix)
CALC_RISE: int = SE_CALC_RISE
CALC_SET: int = SE_CALC_SET
CALC_MTRANSIT: int = SE_CALC_MTRANSIT
CALC_ITRANSIT: int = SE_CALC_ITRANSIT
BIT_DISC_CENTER: int = SE_BIT_DISC_CENTER
BIT_DISC_BOTTOM: int = SE_BIT_DISC_BOTTOM
BIT_NO_REFRACTION: int = SE_BIT_NO_REFRACTION
BIT_CIVIL_TWILIGHT: int = SE_BIT_CIVIL_TWILIGHT
BIT_NAUTIC_TWILIGHT: int = SE_BIT_NAUTIC_TWILIGHT
BIT_ASTRO_TWILIGHT: int = SE_BIT_ASTRO_TWILIGHT
BIT_FIXED_DISC_SIZE: int = SE_BIT_FIXED_DISC_SIZE

# =============================================================================
# HELIACAL EVENT TYPES
# =============================================================================
# Event types for heliacal_ut() function

SE_HELIACAL_RISING: int = 1  # Heliacal rising (morning first)
SE_HELIACAL_SETTING: int = 2  # Heliacal setting (evening last)
SE_MORNING_FIRST: int = SE_HELIACAL_RISING  # Alias: first visibility at morning
SE_EVENING_LAST: int = SE_HELIACAL_SETTING  # Alias: last visibility at evening
SE_EVENING_FIRST: int = 3  # First visibility at evening (after superior conjunction)
SE_MORNING_LAST: int = 4  # Last visibility at morning (before superior conjunction)

# pyswisseph-compatible aliases (without SE_ prefix)
HELIACAL_RISING: int = SE_HELIACAL_RISING
HELIACAL_SETTING: int = SE_HELIACAL_SETTING
MORNING_FIRST: int = SE_MORNING_FIRST
EVENING_LAST: int = SE_EVENING_LAST
EVENING_FIRST: int = SE_EVENING_FIRST
MORNING_LAST: int = SE_MORNING_LAST

# =============================================================================
# HELIACAL VISIBILITY FLAGS
# =============================================================================
# Flags for vis_limit_mag() and heliacal calculations

SE_HELFLAG_OPTICAL_PARAMS: int = 1 << 9  # 512 - Use optical instrument parameters
SE_HELFLAG_NO_DETAILS: int = 1 << 10  # 1024 - Skip detailed calculations
SE_HELFLAG_VISLIM_DARK: int = 1 << 11  # 2048 - Assume Sun at nadir (dark sky)
SE_HELFLAG_VISLIM_NOMOON: int = 1 << 12  # 4096 - Exclude Moon's contribution

# pyswisseph-compatible aliases (without SE_ prefix)
HELFLAG_OPTICAL_PARAMS: int = SE_HELFLAG_OPTICAL_PARAMS
HELFLAG_NO_DETAILS: int = SE_HELFLAG_NO_DETAILS
HELFLAG_VISLIM_DARK: int = SE_HELFLAG_VISLIM_DARK
HELFLAG_VISLIM_NOMOON: int = SE_HELFLAG_VISLIM_NOMOON

# Visibility result codes for vis_limit_mag()
SE_HELFLAG_BELOW_HORIZON: int = -2  # Object is below horizon
SE_HELFLAG_PHOTOPIC: int = 0  # OK, photopic (daylight) vision
SE_HELFLAG_SCOTOPIC: int = 1  # OK, scotopic (night) vision
SE_HELFLAG_MIXED: int = 2  # OK, near limit photopic/scotopic vision

HELFLAG_BELOW_HORIZON: int = SE_HELFLAG_BELOW_HORIZON
HELFLAG_PHOTOPIC: int = SE_HELFLAG_PHOTOPIC
HELFLAG_SCOTOPIC: int = SE_HELFLAG_SCOTOPIC
HELFLAG_MIXED: int = SE_HELFLAG_MIXED

# =============================================================================
# SPLIT_DEG FLAGS
# =============================================================================
# Flags for split_deg() function to control output format and rounding

SPLIT_DEG_ROUND_SEC: int = 1  # Round to seconds
SPLIT_DEG_ROUND_MIN: int = 2  # Round to minutes
SPLIT_DEG_ROUND_DEG: int = 4  # Round to degrees
SPLIT_DEG_ZODIACAL: int = 8  # Return zodiac sign number (0-11)
SPLIT_DEG_NAKSHATRA: int = 1024  # Return nakshatra number (0-26)
SPLIT_DEG_KEEP_SIGN: int = 16  # Don't round to next zodiac sign/nakshatra
SPLIT_DEG_KEEP_DEG: int = 32  # Don't round to next degree

# =============================================================================
# TIDAL ACCELERATION CONSTANTS
# =============================================================================
# Tidal acceleration of the Moon in arcsec/century^2, for Delta T calculations.
# Different JPL ephemeris files use different tidal acceleration values.
# These affect the polynomial extrapolation of Delta T for historical dates.

SE_TIDAL_DE200: float = -23.8946  # DE200 (older ephemeris)
SE_TIDAL_DE403: float = -25.580  # DE403
SE_TIDAL_DE404: float = -25.580  # DE404
SE_TIDAL_DE405: float = -25.826  # DE405
SE_TIDAL_DE406: float = -25.826  # DE406
SE_TIDAL_DE421: float = -25.85  # DE421 (legacy default)
SE_TIDAL_DE422: float = -25.85  # DE422
SE_TIDAL_DE430: float = -25.82  # DE430
SE_TIDAL_DE431: float = -25.80  # DE431
SE_TIDAL_DE440: float = -25.936  # DE440 (current default)
SE_TIDAL_DE441: float = -25.936  # DE441 (latest, same as DE440)
SE_TIDAL_DEFAULT: float = SE_TIDAL_DE440  # Default value based on DE440
SE_TIDAL_AUTOMATIC: float = 0.0  # Let library choose based on ephemeris file

# pyswisseph-compatible aliases (without SE_ prefix)
TIDAL_DE200: float = SE_TIDAL_DE200
TIDAL_DE403: float = SE_TIDAL_DE403
TIDAL_DE404: float = SE_TIDAL_DE404
TIDAL_DE405: float = SE_TIDAL_DE405
TIDAL_DE406: float = SE_TIDAL_DE406
TIDAL_DE421: float = SE_TIDAL_DE421
TIDAL_DE422: float = SE_TIDAL_DE422
TIDAL_DE430: float = SE_TIDAL_DE430
TIDAL_DE431: float = SE_TIDAL_DE431
TIDAL_DE440: float = SE_TIDAL_DE440
TIDAL_DE441: float = SE_TIDAL_DE441
TIDAL_DEFAULT: float = SE_TIDAL_DEFAULT
TIDAL_AUTOMATIC: float = SE_TIDAL_AUTOMATIC

# =============================================================================
# STANDARD ASTRONOMICAL EPOCHS (JULIAN DAY)
# =============================================================================
# Reference epochs used for astrometric calculations and proper motion.
# All values are in Julian Day TT (Terrestrial Time).

# J2000.0 epoch: Jan 1, 2000, 12:00 TT - Standard reference epoch for ICRS
J2000: float = 2451545.0

# J1991.25 epoch: Apr 2, 1991, 13:30 TT - Hipparcos catalog reference epoch
# This is the mean epoch of observations for the Hipparcos mission.
# Proper motion values in Hipparcos are defined relative to this epoch.
# Reference: Hipparcos and Tycho Catalogues, ESA SP-1200, Vol. 1, Section 1.2.2
J1991_25: float = 2448349.0625

# J1900.0 epoch: Jan 0.5, 1900 TT - Historical reference epoch
J1900: float = 2415020.0

# B1950.0 epoch: Jan 0.923, 1950 - Besselian epoch for FK4 catalog
B1950: float = 2433282.4235

# Days per Julian year (exactly 365.25 days)
DAYS_PER_JULIAN_YEAR: float = 365.25

# Days per Julian century (exactly 36525 days)
DAYS_PER_JULIAN_CENTURY: float = 36525.0

# =============================================================================
# PLANETARY MOON IDENTIFIERS
# =============================================================================
# Body IDs for planetary satellites (moons) following Swiss Ephemeris 2.10+ convention
# These require satellite SPK files (jup365.bsp, sat441.bsp, etc.) to be registered
# using register_moon_spk() before calculation.

SE_MOON_OFFSET: int = 9000  # Base offset for planetary moon IDs

# Jupiter's Galilean Moons (discovered by Galileo Galilei in 1610)
SE_MOON_IO: int = SE_MOON_OFFSET + 1  # Jupiter I - innermost Galilean moon
SE_MOON_EUROPA: int = SE_MOON_OFFSET + 2  # Jupiter II - potential subsurface ocean
SE_MOON_GANYMEDE: int = SE_MOON_OFFSET + 3  # Jupiter III - largest moon in solar system
SE_MOON_CALLISTO: int = SE_MOON_OFFSET + 4  # Jupiter IV - heavily cratered

# Saturn's Major Moons
SE_MOON_MIMAS: int = SE_MOON_OFFSET + 11  # Saturn I - "Death Star" appearance
SE_MOON_ENCELADUS: int = SE_MOON_OFFSET + 12  # Saturn II - active geysers
SE_MOON_TETHYS: int = SE_MOON_OFFSET + 13  # Saturn III - icy body
SE_MOON_DIONE: int = SE_MOON_OFFSET + 14  # Saturn IV - trailing hemisphere features
SE_MOON_RHEA: int = SE_MOON_OFFSET + 15  # Saturn V - second largest Saturn moon
SE_MOON_TITAN: int = SE_MOON_OFFSET + 16  # Saturn VI - thick atmosphere, liquid lakes
SE_MOON_HYPERION: int = SE_MOON_OFFSET + 17  # Saturn VII - chaotic rotation
SE_MOON_IAPETUS: int = SE_MOON_OFFSET + 18  # Saturn VIII - two-toned surface

# Uranus' Major Moons
SE_MOON_MIRANDA: int = SE_MOON_OFFSET + 21  # Uranus V - extreme geological features
SE_MOON_ARIEL: int = SE_MOON_OFFSET + 22  # Uranus I - brightest Uranian moon
SE_MOON_UMBRIEL: int = SE_MOON_OFFSET + 23  # Uranus II - darkest Uranian moon
SE_MOON_TITANIA: int = SE_MOON_OFFSET + 24  # Uranus III - largest Uranian moon
SE_MOON_OBERON: int = SE_MOON_OFFSET + 25  # Uranus IV - outermost major Uranian moon

# Neptune's Major Moon
SE_MOON_TRITON: int = SE_MOON_OFFSET + 31  # Neptune I - retrograde orbit, captured KBO

# Mars' Moons
SE_MOON_PHOBOS: int = SE_MOON_OFFSET + 41  # Mars I - larger, closer moon
SE_MOON_DEIMOS: int = SE_MOON_OFFSET + 42  # Mars II - smaller, farther moon

# Pluto's Moon
SE_MOON_CHARON: int = SE_MOON_OFFSET + 51  # Pluto I - largest moon (binary system)

# Aliases without SE_ prefix for pyswisseph compatibility
MOON_IO: int = SE_MOON_IO
MOON_EUROPA: int = SE_MOON_EUROPA
MOON_GANYMEDE: int = SE_MOON_GANYMEDE
MOON_CALLISTO: int = SE_MOON_CALLISTO
MOON_MIMAS: int = SE_MOON_MIMAS
MOON_ENCELADUS: int = SE_MOON_ENCELADUS
MOON_TETHYS: int = SE_MOON_TETHYS
MOON_DIONE: int = SE_MOON_DIONE
MOON_RHEA: int = SE_MOON_RHEA
MOON_TITAN: int = SE_MOON_TITAN
MOON_HYPERION: int = SE_MOON_HYPERION
MOON_IAPETUS: int = SE_MOON_IAPETUS
MOON_MIRANDA: int = SE_MOON_MIRANDA
MOON_ARIEL: int = SE_MOON_ARIEL
MOON_UMBRIEL: int = SE_MOON_UMBRIEL
MOON_TITANIA: int = SE_MOON_TITANIA
MOON_OBERON: int = SE_MOON_OBERON
MOON_TRITON: int = SE_MOON_TRITON
MOON_PHOBOS: int = SE_MOON_PHOBOS
MOON_DEIMOS: int = SE_MOON_DEIMOS
MOON_CHARON: int = SE_MOON_CHARON

# =============================================================================
# NAIF IDS FOR PLANETARY MOONS
# =============================================================================
# Standard NAIF SPICE IDs for planetary satellites

# Jupiter system (5xx)
NAIF_IO: int = 501
NAIF_EUROPA: int = 502
NAIF_GANYMEDE: int = 503
NAIF_CALLISTO: int = 504

# Saturn system (6xx)
NAIF_MIMAS: int = 601
NAIF_ENCELADUS: int = 602
NAIF_TETHYS: int = 603
NAIF_DIONE: int = 604
NAIF_RHEA: int = 605
NAIF_TITAN: int = 606
NAIF_HYPERION: int = 607
NAIF_IAPETUS: int = 608

# Uranus system (7xx)
NAIF_MIRANDA: int = 705
NAIF_ARIEL: int = 701
NAIF_UMBRIEL: int = 702
NAIF_TITANIA: int = 703
NAIF_OBERON: int = 704

# Neptune system (8xx)
NAIF_TRITON: int = 801

# Mars system (4xx)
NAIF_PHOBOS: int = 401
NAIF_DEIMOS: int = 402

# Pluto system (9xx)
NAIF_CHARON: int = 901
