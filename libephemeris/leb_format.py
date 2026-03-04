"""
LEB (LibEphemeris Binary) file format definitions.

This module defines the binary file format constants, dataclasses, and
serialization helpers used by both the generator (write) and reader (read).

The .leb format stores precomputed Chebyshev polynomial coefficients for
fast ephemeris evaluation, providing 10-20x speedup over the Skyfield pipeline.

File layout:
    0x0000       File Header (64 bytes)
    0x0040       Section Directory (N x 24 bytes)
    variable     Section 0: Body Index
    variable     Section 1: Chebyshev Coefficient Data
    variable     Section 2: Nutation Chebyshev Data
    variable     Section 3: Delta-T Sparse Table
    variable     Section 4: Star Catalog
    variable     Section 5: Orbital Elements (reserved)
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from typing import Any, Optional

# =============================================================================
# FORMAT CONSTANTS
# =============================================================================

MAGIC = b"LEB1"
VERSION = 1

# Coordinate types
COORD_ICRS_BARY = 0  # ICRS barycentric (x, y, z) in AU
COORD_ECLIPTIC = 1  # Ecliptic of date (lon, lat, dist) in deg/deg/AU
COORD_HELIO_ECL = 2  # Heliocentric ecliptic (lon, lat, dist) in deg/deg/AU

# Section IDs
SECTION_BODY_INDEX = 0
SECTION_CHEBYSHEV = 1
SECTION_NUTATION = 2
SECTION_DELTA_T = 3
SECTION_STARS = 4
SECTION_ORBITAL_ELEMENTS = 5

# Struct formats (little-endian, no padding with '<' prefix)
HEADER_FMT = "<4sIIIdddI20s"
HEADER_SIZE = struct.calcsize(HEADER_FMT)  # 64

SECTION_DIR_FMT = "<IIQQ"
SECTION_DIR_SIZE = struct.calcsize(SECTION_DIR_FMT)  # 24

BODY_ENTRY_FMT = "<iIIdddIIQ"
BODY_ENTRY_SIZE = struct.calcsize(BODY_ENTRY_FMT)  # 52

STAR_ENTRY_FMT = "<iddddddd4s"
STAR_ENTRY_SIZE = struct.calcsize(STAR_ENTRY_FMT)  # 64

# Nutation section header format
NUTATION_HEADER_FMT = "<dddIIII"
NUTATION_HEADER_SIZE = struct.calcsize(NUTATION_HEADER_FMT)  # 40

# Delta-T entry format (jd, delta_t_days)
DELTA_T_ENTRY_FMT = "<dd"
DELTA_T_ENTRY_SIZE = struct.calcsize(DELTA_T_ENTRY_FMT)  # 16

# Delta-T section header
DELTA_T_HEADER_FMT = "<II"
DELTA_T_HEADER_SIZE = struct.calcsize(DELTA_T_HEADER_FMT)  # 8


# =============================================================================
# DATACLASSES
# =============================================================================


@dataclass
class FileHeader:
    """LEB file header (64 bytes)."""

    magic: bytes
    version: int
    section_count: int
    body_count: int
    jd_start: float
    jd_end: float
    generation_epoch: float
    flags: int


@dataclass
class SectionEntry:
    """Section directory entry (24 bytes)."""

    section_id: int
    offset: int  # Byte offset from file start
    size: int  # Section size in bytes


@dataclass
class BodyEntry:
    """Body index entry (52 bytes)."""

    body_id: int
    coord_type: int  # COORD_ICRS_BARY, COORD_ECLIPTIC, COORD_HELIO_ECL
    segment_count: int
    jd_start: float
    jd_end: float
    interval_days: float
    degree: int
    components: int  # 3 for xyz/lonlatdist
    data_offset: int  # Absolute byte offset to first coefficient


@dataclass
class StarEntry:
    """Fixed star catalog entry (64 bytes)."""

    star_id: int
    ra_j2000: float  # degrees
    dec_j2000: float  # degrees
    pm_ra: float  # deg/yr (includes cos(dec) factor)
    pm_dec: float  # deg/yr
    parallax: float  # arcsec
    rv: float  # km/s
    magnitude: float


@dataclass
class NutationHeader:
    """Nutation section header (40 bytes)."""

    jd_start: float
    jd_end: float
    interval_days: float
    degree: int
    components: int  # 2 (dpsi, deps)
    segment_count: int
    reserved: int  # Padding


# =============================================================================
# BODY PARAMETER TABLE
# =============================================================================
# Maps body_id -> (interval_days, degree, coord_type, components)
# This is the single source of truth for Chebyshev parameters.

BODY_PARAMS: dict[int, tuple[float, int, int, int]] = {
    # SE_SUN through SE_PLUTO, SE_EARTH: ICRS barycentric
    # Parameters: (interval_days, degree, coord_type, components)
    # Tuned for sub-arcsecond fitting error on all bodies.
    0: (32, 13, COORD_ICRS_BARY, 3),  # SE_SUN       — 0.00" (trivial)
    1: (4, 13, COORD_ICRS_BARY, 3),  # SE_MOON      — 0.00" (fast mover)
    2: (16, 15, COORD_ICRS_BARY, 3),  # SE_MERCURY   — 0.00" (eccentric)
    3: (16, 13, COORD_ICRS_BARY, 3),  # SE_VENUS     — was 32: 0.34" → ~0.01"
    4: (16, 13, COORD_ICRS_BARY, 3),  # SE_MARS      — was 32: 0.08" → ~0.003"
    5: (32, 13, COORD_ICRS_BARY, 3),  # SE_JUPITER   — was 64,11: 0.88" → ~0.01"
    6: (32, 13, COORD_ICRS_BARY, 3),  # SE_SATURN    — was 64,11: 0.34" → ~0.01"
    7: (64, 13, COORD_ICRS_BARY, 3),  # SE_URANUS    — 0.04" OK
    8: (64, 13, COORD_ICRS_BARY, 3),  # SE_NEPTUNE   — 0.12" OK
    9: (32, 13, COORD_ICRS_BARY, 3),  # SE_PLUTO     — was 64: 3.16" → ~0.1"
    14: (4, 13, COORD_ICRS_BARY, 3),  # SE_EARTH     — 0.00" (fast mover)
    # Lunar nodes/Lilith: ecliptic direct
    10: (8, 13, COORD_ECLIPTIC, 3),  # SE_MEAN_NODE  (lon, 0, 0)
    11: (8, 13, COORD_ECLIPTIC, 3),  # SE_TRUE_NODE  (lon, lat, dist)
    12: (8, 13, COORD_ECLIPTIC, 3),  # SE_MEAN_APOG  (lon, lat, 0)
    13: (8, 13, COORD_ECLIPTIC, 3),  # SE_OSCU_APOG  (lon, lat, dist)
    21: (8, 13, COORD_ECLIPTIC, 3),  # SE_INTP_APOG  (lon, lat, dist)
    22: (8, 13, COORD_ECLIPTIC, 3),  # SE_INTP_PERG  (lon, lat, dist)
    # Main asteroids: ICRS barycentric
    # Asteroids have eccentric/perturbed orbits — need shorter intervals.
    15: (8, 13, COORD_ICRS_BARY, 3),  # SE_CHIRON    — was 32: 177" → ~0.01"
    17: (8, 13, COORD_ICRS_BARY, 3),  # SE_CERES     — was 32: 317" → ~0.01"
    18: (8, 13, COORD_ICRS_BARY, 3),  # SE_PALLAS    — was 32: 254" → ~0.01"
    19: (8, 13, COORD_ICRS_BARY, 3),  # SE_JUNO      — was 32: 492" → ~0.01"
    20: (8, 13, COORD_ICRS_BARY, 3),  # SE_VESTA     — was 32: 378" → ~0.01"
    # Uranian hypotheticals: heliocentric ecliptic
    # Slow-moving outer bodies — increase degree for smoother fit.
    40: (32, 13, COORD_HELIO_ECL, 3),  # SE_CUPIDO   — was 64,11: bogus 766"
    41: (32, 13, COORD_HELIO_ECL, 3),  # SE_HADES    — was 64,11: bogus 562"
    42: (32, 13, COORD_HELIO_ECL, 3),  # SE_ZEUS     — was 64,11: bogus 448"
    43: (32, 13, COORD_HELIO_ECL, 3),  # SE_KRONOS   — was 64,11: bogus 390"
    44: (32, 13, COORD_HELIO_ECL, 3),  # SE_APOLLON  — was 64,11: bogus 347"
    45: (32, 13, COORD_HELIO_ECL, 3),  # SE_ADMETOS  — was 64,11: bogus 323"
    46: (32, 13, COORD_HELIO_ECL, 3),  # SE_VULKANUS — was 64,11: bogus 304"
    47: (32, 13, COORD_HELIO_ECL, 3),  # SE_POSEIDON — was 64,11: bogus 267"
    48: (32, 13, COORD_HELIO_ECL, 3),  # SE_ISIS     — was 64,11: bogus 136"
}


# =============================================================================
# SERIALIZATION HELPERS
# =============================================================================
# All use struct.pack_into / struct.unpack_from for zero-copy mmap compatibility.


def write_header(buf: bytearray, header: FileHeader) -> None:
    """Pack a FileHeader into a bytearray at offset 0."""
    struct.pack_into(
        HEADER_FMT,
        buf,
        0,
        header.magic,
        header.version,
        header.section_count,
        header.body_count,
        header.jd_start,
        header.jd_end,
        header.generation_epoch,
        header.flags,
        b"\x00" * 20,
    )


def read_header(data: Any, offset: int = 0) -> FileHeader:
    """Unpack a FileHeader from bytes."""
    fields = struct.unpack_from(HEADER_FMT, data, offset)
    return FileHeader(
        magic=fields[0],
        version=fields[1],
        section_count=fields[2],
        body_count=fields[3],
        jd_start=fields[4],
        jd_end=fields[5],
        generation_epoch=fields[6],
        flags=fields[7],
    )


def write_section_dir(buf: bytearray, offset: int, entry: SectionEntry) -> None:
    """Pack a SectionEntry into a bytearray at the given offset."""
    struct.pack_into(
        SECTION_DIR_FMT,
        buf,
        offset,
        entry.section_id,
        0,  # reserved padding
        entry.offset,
        entry.size,
    )


def read_section_dir(data: Any, offset: int) -> SectionEntry:
    """Unpack a SectionEntry from bytes at the given offset."""
    fields = struct.unpack_from(SECTION_DIR_FMT, data, offset)
    return SectionEntry(
        section_id=fields[0],
        # fields[1] is reserved padding
        offset=fields[2],
        size=fields[3],
    )


def write_body_entry(buf: bytearray, offset: int, entry: BodyEntry) -> None:
    """Pack a BodyEntry into a bytearray at the given offset."""
    struct.pack_into(
        BODY_ENTRY_FMT,
        buf,
        offset,
        entry.body_id,
        entry.coord_type,
        entry.segment_count,
        entry.jd_start,
        entry.jd_end,
        entry.interval_days,
        entry.degree,
        entry.components,
        entry.data_offset,
    )


def read_body_entry(data: Any, offset: int) -> BodyEntry:
    """Unpack a BodyEntry from bytes at the given offset."""
    fields = struct.unpack_from(BODY_ENTRY_FMT, data, offset)
    return BodyEntry(
        body_id=fields[0],
        coord_type=fields[1],
        segment_count=fields[2],
        jd_start=fields[3],
        jd_end=fields[4],
        interval_days=fields[5],
        degree=fields[6],
        components=fields[7],
        data_offset=fields[8],
    )


def write_star_entry(buf: bytearray, offset: int, entry: StarEntry) -> None:
    """Pack a StarEntry into a bytearray at the given offset."""
    struct.pack_into(
        STAR_ENTRY_FMT,
        buf,
        offset,
        entry.star_id,
        entry.ra_j2000,
        entry.dec_j2000,
        entry.pm_ra,
        entry.pm_dec,
        entry.parallax,
        entry.rv,
        entry.magnitude,
        b"\x00" * 4,
    )


def read_star_entry(data: Any, offset: int) -> StarEntry:
    """Unpack a StarEntry from bytes at the given offset."""
    fields = struct.unpack_from(STAR_ENTRY_FMT, data, offset)
    return StarEntry(
        star_id=fields[0],
        ra_j2000=fields[1],
        dec_j2000=fields[2],
        pm_ra=fields[3],
        pm_dec=fields[4],
        parallax=fields[5],
        rv=fields[6],
        magnitude=fields[7],
    )


def write_nutation_header(buf: bytearray, offset: int, header: NutationHeader) -> None:
    """Pack a NutationHeader into a bytearray at the given offset."""
    struct.pack_into(
        NUTATION_HEADER_FMT,
        buf,
        offset,
        header.jd_start,
        header.jd_end,
        header.interval_days,
        header.degree,
        header.components,
        header.segment_count,
        header.reserved,
    )


def read_nutation_header(data: Any, offset: int) -> NutationHeader:
    """Unpack a NutationHeader from bytes at the given offset."""
    fields = struct.unpack_from(NUTATION_HEADER_FMT, data, offset)
    return NutationHeader(
        jd_start=fields[0],
        jd_end=fields[1],
        interval_days=fields[2],
        degree=fields[3],
        components=fields[4],
        segment_count=fields[5],
        reserved=fields[6],
    )


def segment_byte_size(degree: int, components: int) -> int:
    """Calculate the byte size of a single Chebyshev segment.

    Args:
        degree: Polynomial degree.
        components: Number of components (3 for xyz or lonlatdist).

    Returns:
        Size in bytes of one segment's coefficient data.
    """
    return components * (degree + 1) * 8  # 8 bytes per f64
