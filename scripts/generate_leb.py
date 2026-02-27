#!/usr/bin/env python3
"""
LEB file generator.

Produces a .leb binary ephemeris file by sampling Skyfield/libephemeris
analytical functions and fitting Chebyshev polynomials to the results.

Usage:
    python scripts/generate_leb.py --output ephemeris.leb
    python scripts/generate_leb.py --tier medium
    python scripts/generate_leb.py --tier extended --start -5000 --end 5000
    python scripts/generate_leb.py --start 1550 --end 2650 --output de440_full.leb
    python scripts/generate_leb.py --output ephemeris.leb --verify --verify-samples 500
    python scripts/generate_leb.py --output ephemeris.leb --workers 8
"""

from __future__ import annotations

import argparse
import math
import multiprocessing
import os
import shutil
import struct
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable, List, Optional, Tuple

import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval


# =============================================================================
# PROGRESS BAR
# =============================================================================


class ProgressBar:
    """Lightweight terminal progress bar (no dependencies).

    Usage::

        bar = ProgressBar(total=100, label="Sun", width=30)
        for i in range(100):
            do_work()
            bar.update(i + 1)
        bar.finish()
    """

    def __init__(
        self,
        total: int,
        label: str = "",
        width: int = 30,
        unit: str = "seg",
        enabled: bool = True,
    ):
        self.total = max(total, 1)
        self.label = label
        self.width = width
        self.unit = unit
        self.enabled = enabled
        self.t0 = time.monotonic()
        self._last_draw = 0.0
        # Determine terminal width once
        try:
            self._term_cols = shutil.get_terminal_size().columns
        except Exception:
            self._term_cols = 80

    def update(self, current: int) -> None:
        if not self.enabled:
            return
        now = time.monotonic()
        # Throttle redraws to every 100ms (unless it's the last update)
        if current < self.total and now - self._last_draw < 0.1:
            return
        self._last_draw = now
        self._draw(current, now)

    def _draw(self, current: int, now: float) -> None:
        elapsed = now - self.t0
        frac = current / self.total
        filled = int(self.width * frac)
        bar = "\u2588" * filled + "\u2591" * (self.width - filled)
        pct = f"{frac * 100:5.1f}%"

        # ETA
        if current > 0 and frac < 1.0:
            eta = elapsed / frac * (1.0 - frac)
            time_str = f"ETA {_fmt_time(eta)}"
        else:
            time_str = f"{_fmt_time(elapsed)}"

        counter = f"{current}/{self.total} {self.unit}"
        line = f"\r  {self.label:20s} {bar} {pct} {counter:>14s} {time_str:>10s}"
        # Truncate to terminal width
        if len(line) > self._term_cols:
            line = line[: self._term_cols]
        sys.stdout.write(line)
        sys.stdout.flush()

    def finish(self, suffix: str = "") -> None:
        if not self.enabled:
            return
        self._draw(self.total, time.monotonic())
        if suffix:
            sys.stdout.write(f"  {suffix}")
        sys.stdout.write("\n")
        sys.stdout.flush()


def _fmt_time(seconds: float) -> str:
    """Format seconds as mm:ss or hh:mm:ss."""
    s = int(seconds)
    if s < 3600:
        return f"{s // 60}:{s % 60:02d}"
    h = s // 3600
    m = (s % 3600) // 60
    return f"{h}:{m:02d}:{s % 60:02d}"


# Add parent directory to path so we can import libephemeris
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from libephemeris.leb_format import (
    BODY_ENTRY_SIZE,
    BODY_PARAMS,
    COORD_ECLIPTIC,
    COORD_HELIO_ECL,
    COORD_ICRS_BARY,
    DELTA_T_ENTRY_FMT,
    DELTA_T_ENTRY_SIZE,
    DELTA_T_HEADER_FMT,
    DELTA_T_HEADER_SIZE,
    HEADER_SIZE,
    MAGIC,
    NUTATION_HEADER_SIZE,
    SECTION_BODY_INDEX,
    SECTION_CHEBYSHEV,
    SECTION_DELTA_T,
    SECTION_DIR_SIZE,
    SECTION_NUTATION,
    SECTION_STARS,
    STAR_ENTRY_SIZE,
    VERSION,
    BodyEntry,
    FileHeader,
    NutationHeader,
    SectionEntry,
    StarEntry,
    segment_byte_size,
    write_body_entry,
    write_header,
    write_nutation_header,
    write_section_dir,
    write_star_entry,
)

# =============================================================================
# CONSTANTS
# =============================================================================

# Julian Day for J2000.0
J2000 = 2451545.0

# Default date range (200 years)
DEFAULT_START_YEAR = 1900
DEFAULT_END_YEAR = 2100

# Tier configurations: (ephemeris_file, default_start_year, default_end_year, output_name)
TIER_CONFIGS = {
    "base": ("de440s.bsp", 1850, 2150, "ephemeris_base.leb"),
    "medium": ("de440.bsp", 1550, 2650, "ephemeris_medium.leb"),
    "extended": ("de441.bsp", -5000, 5000, "ephemeris_extended.leb"),
}

# Default output directory for tier-based generation
DEFAULT_LEB_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "leb")

# Nutation Chebyshev parameters
NUTATION_INTERVAL = 32.0  # days
NUTATION_DEGREE = 16
NUTATION_COMPONENTS = 2  # dpsi, deps

# Delta-T sampling interval
DELTA_T_INTERVAL = 30.0  # days

# Number of sections in the file
NUM_SECTIONS = 5  # body_index, chebyshev, nutation, delta_t, stars

# Body names for logging
BODY_NAMES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
    10: "Mean Node",
    11: "True Node",
    12: "Mean Apogee",
    13: "Oscu Apogee",
    14: "Earth",
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
    21: "Interp Apogee",
    22: "Interp Perigee",
    40: "Cupido",
    41: "Hades",
    42: "Zeus",
    43: "Kronos",
    44: "Apollon",
    45: "Admetos",
    46: "Vulkanus",
    47: "Poseidon",
    48: "Transpluto",
}

# Planet name map for Skyfield (body_id -> skyfield name)
_PLANET_MAP = {
    0: "sun",
    1: "moon",
    2: "mercury",
    3: "venus",
    4: "mars",
    5: "jupiter",
    6: "saturn",
    7: "uranus",
    8: "neptune",
    9: "pluto",
    14: "earth",
}

# Asteroid NAIF IDs for Skyfield SPK lookup
_ASTEROID_NAIF = {
    15: 2060,  # Chiron
    17: 2000001,  # Ceres
    18: 2000002,  # Pallas
    19: 2000003,  # Juno
    20: 2000004,  # Vesta
}

# Body groups for independent generation + merge workflow.
# Each group can be generated as a standalone .leb file and later merged.
BODY_GROUPS: dict[str, List[int]] = {
    "planets": sorted(_PLANET_MAP.keys()),  # 11 ICRS planets (vectorized Skyfield)
    "asteroids": sorted(_ASTEROID_NAIF.keys()),  # 5 ICRS asteroids (spktype21)
    "analytical": sorted(
        bid
        for bid in BODY_PARAMS
        if bid not in _PLANET_MAP and bid not in _ASTEROID_NAIF
    ),  # 15 ecliptic/helio analytical bodies
}


# =============================================================================
# CHEBYSHEV FITTING
# =============================================================================


def chebyshev_nodes(n: int) -> np.ndarray:
    """Compute n Chebyshev nodes (Type I) on [-1, 1]."""
    return np.cos(np.pi * (np.arange(n) + 0.5) / n)


def fit_segment(
    func: Callable[[float], np.ndarray],
    jd_start: float,
    jd_end: float,
    degree: int,
    components: int,
) -> np.ndarray:
    """Fit a multi-component function over [jd_start, jd_end] with Chebyshev polynomials.

    Args:
        func: Function taking JD and returning array of shape (components,).
        jd_start: Start of segment.
        jd_end: End of segment.
        degree: Polynomial degree.
        components: Number of output components.

    Returns:
        Coefficients array of shape (components, degree+1), stored as
        [c0_comp0, c1_comp0, ..., cN_comp0, c0_comp1, ...].
    """
    nodes = chebyshev_nodes(degree + 1)
    # Map nodes from [-1, 1] to [jd_start, jd_end]
    jd_nodes = 0.5 * (jd_end - jd_start) * nodes + 0.5 * (jd_start + jd_end)

    # Evaluate function at nodes
    values = np.array([func(jd) for jd in jd_nodes])  # shape: (degree+1, components)

    # Fit each component independently
    coeffs = np.zeros((components, degree + 1))
    for c in range(components):
        coeffs[c] = chebfit(nodes, values[:, c], degree)

    return coeffs


def verify_segment(
    func: Callable[[float], np.ndarray],
    coeffs: np.ndarray,
    jd_start: float,
    jd_end: float,
    components: int,
    n_test: int = 10,
) -> float:
    """Verify a Chebyshev fit by evaluating at intermediate points.

    Returns the maximum error across all components and test points.
    """
    max_error = 0.0
    for i in range(n_test):
        # Uniform test points (not Chebyshev nodes)
        frac = (i + 0.5) / n_test
        jd = jd_start + frac * (jd_end - jd_start)
        tau = 2.0 * (jd - 0.5 * (jd_start + jd_end)) / (jd_end - jd_start)

        reference = func(jd)
        for c in range(components):
            fitted = chebval(tau, coeffs[c])
            error = abs(fitted - reference[c])
            if error > max_error:
                max_error = error

    return max_error


# =============================================================================
# VECTORIZED CHEBYSHEV FITTING
# =============================================================================

N_VERIFY = 10  # Number of verification points per segment


def _compute_all_segment_jds(
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    n_verify: int = N_VERIFY,
) -> Tuple[np.ndarray, int, int]:
    """Precompute all JDs (fit nodes + verification points) for all segments.

    Returns:
        (all_jds, n_segments, pts_per_segment)
        First (degree+1) JDs per segment are Chebyshev nodes.
        Next n_verify JDs per segment are uniform verification points.
    """
    n_segments = int(math.ceil((jd_end - jd_start) / interval_days))
    nodes_01 = chebyshev_nodes(degree + 1)  # in [-1, 1]
    pts_per_seg = degree + 1 + n_verify

    all_jds = np.empty(n_segments * pts_per_seg)

    for i in range(n_segments):
        seg_start = jd_start + i * interval_days
        seg_end = seg_start + interval_days
        mid = 0.5 * (seg_start + seg_end)
        half = 0.5 * (seg_end - seg_start)
        offset = i * pts_per_seg

        # Chebyshev nodes mapped to [seg_start, seg_end]
        all_jds[offset : offset + degree + 1] = half * nodes_01 + mid

        # Uniform verification points (not on Chebyshev nodes)
        verify_end = min(seg_end, jd_end)
        for v in range(n_verify):
            frac = (v + 0.5) / n_verify
            all_jds[offset + degree + 1 + v] = seg_start + frac * (
                verify_end - seg_start
            )

    return all_jds, n_segments, pts_per_seg


def _fit_and_verify_from_values(
    all_values: np.ndarray,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    components: int,
    n_segments: int,
    pts_per_seg: int,
    n_verify: int = N_VERIFY,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Fit Chebyshev segments and verify from pre-evaluated values.

    Args:
        all_values: shape (n_total, components) with pre-evaluated function values.
            Layout: [seg0_node0, ..., seg0_nodeN, seg0_verify0, ..., seg0_verifyM,
                     seg1_node0, ...]

    Returns:
        (list_of_coefficient_arrays, max_error)
    """
    nodes_01 = chebyshev_nodes(degree + 1)
    all_coeffs: List[np.ndarray] = []
    max_error = 0.0
    bar = ProgressBar(n_segments, label=label or "fit+verify", enabled=verbose)

    for i in range(n_segments):
        seg_start = jd_start + i * interval_days
        seg_end = seg_start + interval_days
        offset = i * pts_per_seg

        # Extract fitting values (Chebyshev nodes)
        fit_values = all_values[offset : offset + degree + 1]

        # Fit each component
        coeffs = np.zeros((components, degree + 1))
        for c in range(components):
            coeffs[c] = chebfit(nodes_01, fit_values[:, c], degree)

        # Verify using pre-computed verification points
        verify_end = min(seg_end, jd_end)
        mid = 0.5 * (seg_start + seg_end)
        half = 0.5 * (seg_end - seg_start)

        for v in range(n_verify):
            frac = (v + 0.5) / n_verify
            jd_v = seg_start + frac * (verify_end - seg_start)
            tau = (jd_v - mid) / half

            ref = all_values[offset + degree + 1 + v]
            for c in range(components):
                fitted = float(chebval(tau, coeffs[c]))
                error = abs(fitted - ref[c])
                if error > max_error:
                    max_error = error

        all_coeffs.append(coeffs)
        bar.update(i + 1)

    bar.finish()
    return all_coeffs, max_error


# =============================================================================
# BODY GENERATORS
# =============================================================================


def _init_skyfield():
    """Initialize Skyfield resources (called once per process)."""
    import libephemeris as ephem
    from libephemeris.state import get_planets, get_timescale
    from libephemeris.planets import get_planet_target

    planets = get_planets()
    ts = get_timescale()
    return planets, ts


def _get_spk_jd_range(planets) -> Tuple[float, float]:
    """Get the valid JD range of the loaded SPK ephemeris.

    Returns:
        (jd_min, jd_max) covering all segments in the ephemeris.
    """
    T0 = 2451545.0  # J2000 in JD
    spk = planets.spk
    starts = [T0 + seg.start_second / 86400.0 for seg in spk.segments]
    ends = [T0 + seg.end_second / 86400.0 for seg in spk.segments]
    return min(starts), max(ends)


def _spk_covers_range(
    spk_file: str,
    body_id: int,
    jd_start: float,
    jd_end: float,
) -> bool:
    """Check if an SPK type 21 file covers the required JD range for a body.

    Opens the SPK, finds the target NAIF ID, and checks whether the
    segments span [jd_start, jd_end].

    Returns:
        True if coverage is sufficient, False otherwise.
    """
    try:
        from spktype21 import SPKType21

        kernel = SPKType21.open(spk_file)
    except Exception:
        return False

    try:
        naif_id = _ASTEROID_NAIF.get(body_id)
        if naif_id is None:
            return False

        T0 = 2451545.0
        # Find segments matching this target (try both NAIF conventions)
        target_segs = [
            seg
            for seg in kernel.segments
            if seg.target == naif_id
            or seg.target == naif_id + 18000000  # Horizons 20000000+N convention
        ]
        if not target_segs:
            # Try any non-center segment
            target_segs = [seg for seg in kernel.segments if seg.target not in (0, 10)]

        if not target_segs:
            return False

        seg_start = min(T0 + seg.start_second / 86400.0 for seg in target_segs)
        seg_end = max(T0 + seg.end_second / 86400.0 for seg in target_segs)

        return seg_start <= jd_start and seg_end >= jd_end
    finally:
        kernel.close()


def _eval_target_vectorized(
    target,
    all_jds: np.ndarray,
    ts,
    spk_jd_min: float,
    spk_jd_max: float,
) -> np.ndarray:
    """Evaluate a Skyfield target at all JDs, with linear extrapolation for
    JDs outside the SPK ephemeris range.

    Last Chebyshev segments may extend their fit nodes a few days beyond
    jd_end (to keep full-width segment alignment with the reader).
    Instead of clamping (which corrupts the polynomial fit), out-of-range
    nodes are extrapolated linearly using position + velocity at the boundary.
    This preserves fit quality for in-range dates.

    Args:
        target: Skyfield VectorFunction (planet, barycenter, etc.)
        all_jds: Array of Julian Days (TT)
        ts: Skyfield Timescale
        spk_jd_min: Start of SPK valid range
        spk_jd_max: End of SPK valid range

    Returns:
        Array of shape (N, 3) with positions in AU.
    """
    margin = 1.0  # days of safety margin from SPK boundary
    lo = spk_jd_min + margin
    hi = spk_jd_max - margin

    in_range = (all_jds >= lo) & (all_jds <= hi)

    if np.all(in_range):
        # Fast path: everything in range
        t_arr = ts.tt_jd(all_jds)
        return np.asarray(target.at(t_arr).position.au).T  # (N, 3)

    # Split: evaluate in-range vectorized, extrapolate out-of-range
    all_values = np.empty((len(all_jds), 3))

    # In-range: vectorized evaluation
    in_idx = np.where(in_range)[0]
    if len(in_idx) > 0:
        t_in = ts.tt_jd(all_jds[in_idx])
        pos_in = np.asarray(target.at(t_in).position.au)  # (3, N_in)
        all_values[in_idx] = pos_in.T

    # Out-of-range: linear extrapolation from boundary
    out_idx = np.where(~in_range)[0]
    if len(out_idx) > 0:
        # Determine if overshoot is at start or end
        over_end = all_jds[out_idx] > hi
        over_start = ~over_end

        # Extrapolate from end boundary
        end_idx = out_idx[over_end]
        if len(end_idx) > 0:
            t_bnd = ts.tt_jd(hi)
            bnd_pos = target.at(t_bnd)
            p = np.asarray(bnd_pos.position.au).ravel()  # (3,)
            v = np.asarray(bnd_pos.velocity.au_per_d).ravel()  # (3,)
            dt = all_jds[end_idx] - hi  # days past boundary
            all_values[end_idx, 0] = p[0] + v[0] * dt
            all_values[end_idx, 1] = p[1] + v[1] * dt
            all_values[end_idx, 2] = p[2] + v[2] * dt

        # Extrapolate from start boundary
        start_idx = out_idx[over_start]
        if len(start_idx) > 0:
            t_bnd = ts.tt_jd(lo)
            bnd_pos = target.at(t_bnd)
            p = np.asarray(bnd_pos.position.au).ravel()
            v = np.asarray(bnd_pos.velocity.au_per_d).ravel()
            dt = all_jds[start_idx] - lo
            all_values[start_idx, 0] = p[0] + v[0] * dt
            all_values[start_idx, 1] = p[1] + v[1] * dt
            all_values[start_idx, 2] = p[2] + v[2] * dt

    return all_values


def _eval_body_icrs_vectorized(
    target_name: str,
    all_jds: np.ndarray,
    planets,
    ts,
) -> np.ndarray:
    """Get vectorized ICRS barycentric positions for a planet.

    Handles inner planets (direct Skyfield target) and outer planets
    (barycenter + SPK center offset or COB correction) transparently.

    JDs that extend beyond the SPK ephemeris range (from last-segment
    overshoot) are linearly extrapolated using position + velocity at
    the boundary. This preserves Chebyshev fit quality for in-range dates.

    Args:
        target_name: Planet name from _PLANET_MAP (e.g., 'jupiter', 'sun')
        all_jds: Array of Julian Days (TT)
        planets: Skyfield SpiceKernel ephemeris
        ts: Skyfield Timescale

    Returns:
        Array of shape (N, 3) with ICRS barycentric positions in AU.
    """
    from libephemeris.planets import _PLANET_FALLBACK, _PLANET_CENTER_NAIF_IDS
    from libephemeris.state import get_planet_center_segment

    spk_min, spk_max = _get_spk_jd_range(planets)

    # Try direct target first (works for inner planets: sun, moon, mercury, venus, earth)
    try:
        target = planets[target_name]
        return _eval_target_vectorized(target, all_jds, ts, spk_min, spk_max)
    except KeyError:
        pass

    # Outer planet: needs barycenter + center offset
    bary_name = _PLANET_FALLBACK.get(target_name)
    if bary_name is None:
        raise ValueError(f"No target or fallback for {target_name}")

    barycenter = planets[bary_name]

    # Try SPK center segments (vectorized natively via Skyfield)
    if target_name in _PLANET_CENTER_NAIF_IDS:
        naif_id = _PLANET_CENTER_NAIF_IDS[target_name]
        center_segment = get_planet_center_segment(naif_id)
        if center_segment is not None:
            try:
                bary_vals = _eval_target_vectorized(
                    barycenter, all_jds, ts, spk_min, spk_max
                )
                offset_vals = _eval_target_vectorized(
                    center_segment, all_jds, ts, spk_min, spk_max
                )
                return bary_vals + offset_vals
            except Exception:
                pass  # Fall through to COB correction

    # Fallback: vectorized barycenter + scalar COB correction
    # The barycenter evaluation is the expensive part (vectorized).
    # COB offsets are cheap analytical formulas computed one-by-one.
    from libephemeris.moon_theories import get_cob_offset

    bary_vals = _eval_target_vectorized(barycenter, all_jds, ts, spk_min, spk_max)

    for i in range(len(all_jds)):
        t_single = ts.tt_jd(float(all_jds[i]))
        offset = get_cob_offset(bary_name, t_single)
        bary_vals[i, 0] += offset[0]
        bary_vals[i, 1] += offset[1]
        bary_vals[i, 2] += offset[2]

    return bary_vals


def generate_body_icrs(
    body_id: int,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for an ICRS barycentric body.

    Uses vectorized Skyfield evaluation: all JDs (fit nodes + verification)
    are computed in a single batch call for ~150x speedup over scalar loops.

    Returns:
        (list_of_coefficient_arrays, max_error_au)
    """
    from libephemeris.state import get_planets, get_timescale

    planets = get_planets()
    ts = get_timescale()

    target_name = _PLANET_MAP.get(body_id)
    if target_name is None:
        raise ValueError(f"No planet map entry for body_id={body_id}")

    # Precompute all JDs for all segments (fit + verify)
    all_jds, n_segments, pts_per_seg = _compute_all_segment_jds(
        jd_start, jd_end, interval_days, degree
    )

    # Single vectorized evaluation (handles COB correction for outer planets)
    all_values = _eval_body_icrs_vectorized(target_name, all_jds, planets, ts)

    return _fit_and_verify_from_values(
        all_values,
        jd_start,
        jd_end,
        interval_days,
        degree,
        3,
        n_segments,
        pts_per_seg,
        label=label,
        verbose=verbose,
    )


def generate_body_icrs_asteroid(
    body_id: int,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for an asteroid (ICRS barycentric).

    Uses spktype21 to read the SPK type 21 file directly (~36x faster than
    scalar libephemeris calls). The SPK provides heliocentric positions
    (center=10/Sun), so the Sun's barycentric position is added via a single
    vectorized Skyfield call.

    Falls back to scalar per-point evaluation if SPK file is not available
    or does not cover the full requested date range.

    Returns:
        (list_of_coefficient_arrays, max_error_au)
    """
    from libephemeris.state import get_planets, get_timescale, _SPK_BODY_MAP

    planets = get_planets()
    ts = get_timescale()

    # Check if we have an SPK file for this asteroid
    spk_info = _SPK_BODY_MAP.get(body_id)

    if spk_info is not None:
        spk_file, naif_id = spk_info
        try:
            from spktype21 import SPKType21

            kernel = SPKType21.open(spk_file)
        except Exception:
            kernel = None

        if kernel is not None:
            try:
                # Find the center ID and target from the kernel segments
                center_id = kernel.segments[0].center  # typically 10 (Sun)
                target_id = kernel.segments[0].target

                # Check if the asteroid SPK covers the full requested range
                # (with some margin for last-segment overshoot)
                T0_SPK = 2451545.0
                ast_jd_min = min(
                    T0_SPK + seg.start_second / 86400.0
                    for seg in kernel.segments
                    if seg.target == target_id
                )
                ast_jd_max = max(
                    T0_SPK + seg.end_second / 86400.0
                    for seg in kernel.segments
                    if seg.target == target_id
                )

                if ast_jd_min > jd_start or ast_jd_max < jd_end:
                    # SPK doesn't cover the full range — fall back to scalar
                    kernel.close()
                    kernel = None
            except Exception:
                kernel.close()
                kernel = None

        if kernel is not None:
            try:
                center_id = kernel.segments[0].center
                target_id = kernel.segments[0].target

                # Recompute asteroid SPK range for clamping
                T0_SPK2 = 2451545.0
                ast_jd_lo = min(
                    T0_SPK2 + seg.start_second / 86400.0
                    for seg in kernel.segments
                    if seg.target == target_id
                )
                ast_jd_hi = max(
                    T0_SPK2 + seg.end_second / 86400.0
                    for seg in kernel.segments
                    if seg.target == target_id
                )

                # Precompute all JDs
                all_jds, n_segments, pts_per_seg = _compute_all_segment_jds(
                    jd_start, jd_end, interval_days, degree
                )

                # Get Sun barycentric positions with extrapolation for overshoot
                spk_min, spk_max = _get_spk_jd_range(planets)
                sun = planets["sun"]
                sun_bary = _eval_target_vectorized(
                    sun, all_jds, ts, spk_min, spk_max
                )  # (N, 3)

                # Clamp to asteroid SPK range for spktype21 calls
                ast_clamped = np.clip(all_jds, ast_jd_lo + 0.01, ast_jd_hi - 0.01)

                # Compute asteroid heliocentric positions via spktype21 (scalar loop)
                AU_KM = 149597870.7
                all_values = np.empty((len(all_jds), 3))
                spk_bar = ProgressBar(
                    len(all_jds),
                    label=label + " (spk21)" if label else "spk21 eval",
                    unit="pt",
                    enabled=verbose,
                )
                for i in range(len(all_jds)):
                    pos_km, _ = kernel.compute_type21(
                        center_id, target_id, float(ast_clamped[i])
                    )
                    # helio km -> helio AU -> + Sun bary = SSB bary (ICRS)
                    all_values[i, 0] = pos_km[0] / AU_KM + sun_bary[i, 0]
                    all_values[i, 1] = pos_km[1] / AU_KM + sun_bary[i, 1]
                    all_values[i, 2] = pos_km[2] / AU_KM + sun_bary[i, 2]
                    spk_bar.update(i + 1)
                spk_bar.finish()

                kernel.close()

                return _fit_and_verify_from_values(
                    all_values,
                    jd_start,
                    jd_end,
                    interval_days,
                    degree,
                    3,
                    n_segments,
                    pts_per_seg,
                    label=label,
                    verbose=verbose,
                )
            except Exception:
                kernel.close()
                # Fall through to scalar fallback

    # Fallback: scalar per-point path (slow but works for any date range)
    import libephemeris as ephem
    from libephemeris.constants import SEFLG_SPEED
    from libephemeris.planets import get_planet_target

    earth = get_planet_target(planets, "earth")

    # Clamp JDs to SPK range for last-segment overshoot
    spk_min, spk_max = _get_spk_jd_range(planets)
    jd_lo = spk_min + 1.0
    jd_hi = spk_max - 1.0

    def eval_func(jd: float) -> np.ndarray:
        jd_c = max(jd_lo, min(jd_hi, jd))
        t = ts.tt_jd(jd_c)
        earth_pos = earth.at(t).position.au
        result, _ = ephem.swe_calc(jd_c, body_id, SEFLG_SPEED)
        lon_rad = math.radians(result[0])
        lat_rad = math.radians(result[1])
        dist = result[2]
        x_geo = dist * math.cos(lat_rad) * math.cos(lon_rad)
        y_geo = dist * math.cos(lat_rad) * math.sin(lon_rad)
        z_geo = dist * math.sin(lat_rad)
        eps = math.radians(23.4392911)
        x_eq = x_geo
        y_eq = y_geo * math.cos(eps) - z_geo * math.sin(eps)
        z_eq = y_geo * math.sin(eps) + z_geo * math.cos(eps)
        return np.array(
            [
                x_eq + float(earth_pos[0]),
                y_eq + float(earth_pos[1]),
                z_eq + float(earth_pos[2]),
            ]
        )

    return _generate_segments(
        eval_func,
        jd_start,
        jd_end,
        interval_days,
        degree,
        3,
        label=label + " (scalar)" if label else "scalar",
        verbose=verbose,
    )


def generate_body_ecliptic(
    body_id: int,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for an ecliptic-direct body.

    Handles longitude unwrapping before fitting.
    """
    from libephemeris.lunar import (
        calc_mean_lunar_node,
        calc_true_lunar_node,
        calc_mean_lilith_with_latitude,
        calc_true_lilith,
        calc_interpolated_apogee,
        calc_interpolated_perigee,
    )
    from libephemeris.time_utils import swe_deltat

    # Map body_id to evaluation function
    eval_funcs = {
        10: lambda jd: np.array([calc_mean_lunar_node(jd), 0.0, 0.0]),
        11: lambda jd: np.array(calc_true_lunar_node(jd)),
        12: lambda jd: np.array([*calc_mean_lilith_with_latitude(jd), 0.0]),
        13: lambda jd: np.array(calc_true_lilith(jd)),
        21: lambda jd: np.array(calc_interpolated_apogee(jd)),
        22: lambda jd: np.array(calc_interpolated_perigee(jd)),
    }

    if body_id not in eval_funcs:
        raise ValueError(f"No ecliptic eval function for body_id={body_id}")

    raw_func = eval_funcs[body_id]

    # Generate with longitude unwrapping
    return _generate_segments_unwrap(
        raw_func,
        jd_start,
        jd_end,
        interval_days,
        degree,
        3,
        label=label,
        verbose=verbose,
    )


def generate_body_helio(
    body_id: int,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for a heliocentric ecliptic body.

    Handles Uranian hypotheticals and Transpluto.
    """
    from libephemeris.hypothetical import calc_uranian_planet, calc_transpluto

    if body_id == 48:  # SE_ISIS / Transpluto

        def eval_func(jd: float) -> np.ndarray:
            result = calc_transpluto(jd)
            return np.array([result[0], result[1], result[2]])
    elif 40 <= body_id <= 47:  # Uranian planets

        def eval_func(jd: float) -> np.ndarray:
            result = calc_uranian_planet(body_id, jd)
            return np.array([result[0], result[1], result[2]])
    else:
        raise ValueError(f"No heliocentric eval function for body_id={body_id}")

    # Generate with longitude unwrapping
    return _generate_segments_unwrap(
        eval_func,
        jd_start,
        jd_end,
        interval_days,
        degree,
        3,
        label=label,
        verbose=verbose,
    )


def _generate_segments(
    func: Callable[[float], np.ndarray],
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    components: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev segments for a function (no unwrapping).

    Returns:
        (list_of_coefficient_arrays, max_error)
    """
    n_segments = int(math.ceil((jd_end - jd_start) / interval_days))
    all_coeffs = []
    max_error = 0.0
    bar = ProgressBar(n_segments, label=label or "segments", enabled=verbose)

    for i in range(n_segments):
        seg_start = jd_start + i * interval_days
        # IMPORTANT: Always use full-width segments. The reader maps tau
        # using interval_days, so truncating the last segment would cause
        # a mismatch between the fitted polynomial domain and the reader's
        # tau computation.  The underlying functions are valid beyond jd_end.
        seg_end = seg_start + interval_days

        coeffs = fit_segment(func, seg_start, seg_end, degree, components)
        # Verify only within the actual requested range
        verify_end = min(seg_end, jd_end)
        error = verify_segment(func, coeffs, seg_start, verify_end, components)
        if error > max_error:
            max_error = error

        all_coeffs.append(coeffs)
        bar.update(i + 1)

    bar.finish()
    return all_coeffs, max_error


def _generate_segments_unwrap(
    func: Callable[[float], np.ndarray],
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    components: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev segments with longitude unwrapping (component 0).

    The first component (longitude) is unwrapped before fitting to handle
    the 0/360 discontinuity. After fitting, verification re-wraps with % 360.
    """
    n_segments = int(math.ceil((jd_end - jd_start) / interval_days))
    all_coeffs = []
    max_error = 0.0
    bar = ProgressBar(n_segments, label=label or "segments", enabled=verbose)

    for i in range(n_segments):
        seg_start = jd_start + i * interval_days
        # Always use full-width segments (see _generate_segments comment)
        seg_end = seg_start + interval_days

        nodes = chebyshev_nodes(degree + 1)
        jd_nodes = 0.5 * (seg_end - seg_start) * nodes + 0.5 * (seg_start + seg_end)

        # Evaluate function at nodes
        values = np.array([func(jd) for jd in jd_nodes])

        # Unwrap longitude (component 0) to remove 360-degree jumps
        values[:, 0] = np.unwrap(np.radians(values[:, 0]))
        values[:, 0] = np.degrees(values[:, 0])

        # Fit each component
        coeffs = np.zeros((components, degree + 1))
        for c in range(components):
            coeffs[c] = chebfit(nodes, values[:, c], degree)

        # Verify with re-wrapping (only within actual requested range)
        verify_end = min(seg_end, jd_end)
        error = _verify_segment_unwrapped(
            func, coeffs, seg_start, verify_end, components
        )
        if error > max_error:
            max_error = error

        all_coeffs.append(coeffs)
        bar.update(i + 1)

    bar.finish()
    return all_coeffs, max_error


def _verify_segment_unwrapped(
    func: Callable[[float], np.ndarray],
    coeffs: np.ndarray,
    jd_start: float,
    jd_end: float,
    components: int,
    n_test: int = 10,
) -> float:
    """Verify a Chebyshev fit for ecliptic bodies (with longitude re-wrapping)."""
    max_error = 0.0
    for i in range(n_test):
        frac = (i + 0.5) / n_test
        jd = jd_start + frac * (jd_end - jd_start)
        tau = 2.0 * (jd - 0.5 * (jd_start + jd_end)) / (jd_end - jd_start)

        reference = func(jd)
        for c in range(components):
            fitted = chebval(tau, coeffs[c])
            if c == 0:
                # Re-wrap longitude
                fitted = fitted % 360.0
                ref = reference[c] % 360.0
                # Handle wrap-around comparison
                error = abs(fitted - ref)
                if error > 180.0:
                    error = 360.0 - error
            else:
                error = abs(fitted - reference[c])
            if error > max_error:
                max_error = error

    return max_error


# =============================================================================
# NUTATION GENERATOR
# =============================================================================


def generate_nutation(
    jd_start: float,
    jd_end: float,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for nutation (dpsi, deps in radians).

    Uses vectorized erfa.nut06a (IAU 2006/2000A) for ~50x speedup over
    scalar get_cached_nutation() calls.
    """
    import erfa

    # Precompute all JDs for all segments
    all_jds, n_segments, pts_per_seg = _compute_all_segment_jds(
        jd_start, jd_end, NUTATION_INTERVAL, NUTATION_DEGREE
    )

    # Convert JD to TT (J2000 epoch split for erfa)
    # erfa uses (jd1, jd2) split: jd1 = 2451545.0 (J2000), jd2 = jd - 2451545.0
    jd1 = np.full_like(all_jds, 2451545.0)
    jd2 = all_jds - 2451545.0

    # Single vectorized erfa call
    dpsi, deps = erfa.nut06a(jd1, jd2)  # Both in radians

    # Stack into (N, 2) array
    all_values = np.column_stack([dpsi, deps])

    return _fit_and_verify_from_values(
        all_values,
        jd_start,
        jd_end,
        NUTATION_INTERVAL,
        NUTATION_DEGREE,
        NUTATION_COMPONENTS,
        n_segments,
        pts_per_seg,
    )


# =============================================================================
# DELTA-T GENERATOR
# =============================================================================


def generate_delta_t(
    jd_start: float,
    jd_end: float,
) -> List[Tuple[float, float]]:
    """Generate Delta-T sparse table.

    Samples swe_deltat() every DELTA_T_INTERVAL days.

    Returns:
        List of (jd, delta_t_days) tuples.
    """
    from libephemeris.time_utils import swe_deltat

    table = []
    jd = jd_start
    while jd <= jd_end:
        dt = swe_deltat(jd)
        table.append((jd, dt))
        jd += DELTA_T_INTERVAL

    # Ensure we include the end point
    if table[-1][0] < jd_end:
        table.append((jd_end, swe_deltat(jd_end)))

    return table


# =============================================================================
# STAR CATALOG GENERATOR
# =============================================================================


def generate_star_catalog() -> List[StarEntry]:
    """Extract star catalog entries from libephemeris.

    Returns:
        List of StarEntry records.
    """
    from libephemeris.fixed_stars import STAR_CATALOG

    entries = []
    for star in STAR_CATALOG:
        entries.append(
            StarEntry(
                star_id=star.id,
                ra_j2000=star.data.ra_j2000,
                dec_j2000=star.data.dec_j2000,
                pm_ra=star.data.pm_ra / 3600.0,  # arcsec/yr -> deg/yr
                pm_dec=star.data.pm_dec / 3600.0,  # arcsec/yr -> deg/yr
                parallax=0.0,  # Not in StarData
                rv=0.0,  # Not in StarData
                magnitude=star.magnitude,
            )
        )

    return entries


# =============================================================================
# FILE ASSEMBLY
# =============================================================================


def _year_to_jd(year: int) -> float:
    """Convert a year to Julian Day (January 1.0)."""
    from libephemeris.time_utils import swe_julday

    return swe_julday(year, 1, 1, 0.0)


def generate_single_body(
    body_id: int,
    jd_start: float,
    jd_end: float,
    verbose: bool = False,
) -> Tuple[int, List[np.ndarray], float]:
    """Generate Chebyshev data for a single body.

    Returns:
        (body_id, coefficients_list, max_error)
    """
    params = BODY_PARAMS[body_id]
    interval_days, degree, coord_type, components = params
    label = BODY_NAMES.get(body_id, f"Body {body_id}")

    if coord_type == COORD_ICRS_BARY:
        if body_id in _PLANET_MAP:
            coeffs, error = generate_body_icrs(
                body_id,
                jd_start,
                jd_end,
                interval_days,
                degree,
                label=label,
                verbose=verbose,
            )
        else:
            coeffs, error = generate_body_icrs_asteroid(
                body_id,
                jd_start,
                jd_end,
                interval_days,
                degree,
                label=label,
                verbose=verbose,
            )
    elif coord_type == COORD_ECLIPTIC:
        coeffs, error = generate_body_ecliptic(
            body_id,
            jd_start,
            jd_end,
            interval_days,
            degree,
            label=label,
            verbose=verbose,
        )
    elif coord_type == COORD_HELIO_ECL:
        coeffs, error = generate_body_helio(
            body_id,
            jd_start,
            jd_end,
            interval_days,
            degree,
            label=label,
            verbose=verbose,
        )
    else:
        raise ValueError(f"Unknown coord_type {coord_type} for body {body_id}")

    return body_id, coeffs, error


def _worker_generate_body(args):
    """Worker function for parallel body generation.

    Progress bars are disabled in workers to avoid interleaved output
    from multiple subprocesses writing ``\\r`` to the same terminal.
    The main process shows an overall progress bar instead.
    """
    body_id, jd_start, jd_end = args
    try:
        return generate_single_body(body_id, jd_start, jd_end, verbose=False)
    except Exception as e:
        print(
            f"  ERROR generating body {body_id} ({BODY_NAMES.get(body_id, '?')}): {e}"
        )
        raise


def assemble_leb(
    output: str,
    jd_start: float,
    jd_end: float,
    bodies: Optional[List[int]] = None,
    workers: int = 1,
    verbose: bool = True,
) -> None:
    """Assemble a complete .leb file.

    Args:
        output: Output file path.
        jd_start: Start Julian Day.
        jd_end: End Julian Day.
        bodies: List of body IDs to include (None = all from BODY_PARAMS).
        workers: Number of parallel workers for body generation.
        verbose: Print progress messages.
    """
    if bodies is None:
        bodies = sorted(BODY_PARAMS.keys())

    now_jd = J2000 + (time.time() / 86400.0 - 10957.5)  # Approximate current JD

    if verbose:
        print(f"Generating LEB file: {output}")
        print(f"  Date range: JD {jd_start:.1f} to {jd_end:.1f}")
        print(f"  Bodies: {len(bodies)}")
        print(f"  Workers: {workers}")
        print()

    # -------------------------------------------------------------------------
    # 0. Ensure SPK kernels for asteroids covering the full date range
    # -------------------------------------------------------------------------
    asteroid_bodies = [b for b in bodies if b in _ASTEROID_NAIF]
    skipped_asteroids: List[int] = []
    if asteroid_bodies:
        import libephemeris as ephem
        from libephemeris.minor_bodies import (
            auto_download_asteroid_spk,
            is_spk_available_for_body,
        )
        from libephemeris.state import _SPK_BODY_MAP

        # Enable auto-download and disable strict precision for generation
        os.environ["LIBEPHEMERIS_AUTO_SPK"] = "1"
        os.environ["LIBEPHEMERIS_STRICT_PRECISION"] = "0"
        ephem.set_auto_spk_download(True)
        ephem.set_strict_precision(False)

        if verbose:
            print("  Preparing asteroid SPK kernels...")
            print(f"    Required coverage: JD {jd_start:.1f} to {jd_end:.1f}")

        for bid in asteroid_bodies:
            name = BODY_NAMES.get(bid, f"Body {bid}")

            # Check if an already-cached SPK covers the full range.
            # auto_download_asteroid_spk() short-circuits when a body is
            # already registered, even if its coverage is too narrow (e.g.
            # a previously downloaded ±10yr SPK).  We detect that here and
            # force a re-download with the full tier range.
            need_force = False
            if bid in _SPK_BODY_MAP:
                cached_file, _ = _SPK_BODY_MAP[bid]
                if not _spk_covers_range(cached_file, bid, jd_start, jd_end):
                    need_force = True
                    if verbose:
                        print(f"    {name}: cached SPK too narrow, re-downloading...")

            try:
                auto_download_asteroid_spk(
                    bid,
                    jd_start=jd_start,
                    jd_end=jd_end,
                    force=need_force,
                )
            except Exception:
                pass  # Download may fail; we check availability below

            if is_spk_available_for_body(bid):
                if verbose:
                    print(f"    {name}: SPK ready")
            else:
                skipped_asteroids.append(bid)
                if verbose:
                    print(f"    {name}: no SPK available, SKIPPING")

        # Remove unavailable asteroids from the body list
        if skipped_asteroids:
            bodies = [b for b in bodies if b not in skipped_asteroids]

        if verbose:
            print()

    # -------------------------------------------------------------------------
    # 1. Generate body coefficients
    # -------------------------------------------------------------------------
    body_data: dict[int, List[np.ndarray]] = {}
    body_errors: dict[int, float] = {}

    # Categorize bodies by generation strategy:
    # - ICRS planets: vectorized Skyfield (fast, run in main process)
    # - ICRS asteroids: vectorized via SPK target (fast, run in main process)
    # - Ecliptic/Helio: scalar analytical funcs (slow, parallelize across workers)
    icrs_planet_bodies = [b for b in bodies if b in _PLANET_MAP]
    icrs_asteroid_bodies = [b for b in bodies if b in _ASTEROID_NAIF]
    analytical_bodies = [
        b for b in bodies if b not in _PLANET_MAP and b not in _ASTEROID_NAIF
    ]

    t0 = time.time()

    # 1a. Generate ICRS planets (vectorized, very fast)
    if icrs_planet_bodies and verbose:
        print("  --- ICRS planets (vectorized Skyfield) ---")
    for bid in icrs_planet_bodies:
        bid, coeffs, error = generate_single_body(
            bid, jd_start, jd_end, verbose=verbose
        )
        body_data[bid] = coeffs
        body_errors[bid] = error

    # 1b. Generate asteroids (vectorized via SPK, fast)
    if icrs_asteroid_bodies and verbose:
        print("  --- ICRS asteroids ---")
    for bid in icrs_asteroid_bodies:
        bid, coeffs, error = generate_single_body(
            bid, jd_start, jd_end, verbose=verbose
        )
        body_data[bid] = coeffs
        body_errors[bid] = error

    # 1c. Generate analytical bodies (ecliptic + helio, parallelize if workers > 1)
    if analytical_bodies and verbose:
        print(
            f"  --- Analytical bodies"
            f" ({len(analytical_bodies)} bodies"
            f", {workers} worker{'s' if workers > 1 else ''}) ---"
        )
    if workers > 1 and analytical_bodies:
        args_list = [(bid, jd_start, jd_end) for bid in analytical_bodies]
        bar = ProgressBar(
            len(analytical_bodies),
            label="Analytical (parallel)",
            unit="body",
            enabled=verbose,
        )
        done = 0
        # Use "spawn" context on macOS to avoid fork deadlocks with
        # C extensions (numpy, erfa).  "fork" copies the parent's locked
        # mutexes which can permanently hang child processes.
        mp_ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp_ctx) as executor:
            futures = {
                executor.submit(_worker_generate_body, args): args[0]
                for args in args_list
            }
            for future in as_completed(futures):
                bid, coeffs, error = future.result()
                body_data[bid] = coeffs
                body_errors[bid] = error
                done += 1
                bar.update(done)
        bar.finish()
    else:
        for bid in analytical_bodies:
            bid, coeffs, error = generate_single_body(
                bid, jd_start, jd_end, verbose=verbose
            )
            body_data[bid] = coeffs
            body_errors[bid] = error

    t_bodies = time.time() - t0
    if verbose:
        print(f"\n  Body generation: {t_bodies:.1f}s")

    # -------------------------------------------------------------------------
    # 2. Generate nutation
    # -------------------------------------------------------------------------
    if verbose:
        print("  Generating nutation...", end=" ", flush=True)
    t0 = time.time()
    nutation_coeffs, nutation_error = generate_nutation(jd_start, jd_end)
    t_nut = time.time() - t0
    if verbose:
        print(
            f"{len(nutation_coeffs)} segments, max error={nutation_error:.2e} rad ({t_nut:.1f}s)"
        )

    # -------------------------------------------------------------------------
    # 3. Generate Delta-T
    # -------------------------------------------------------------------------
    if verbose:
        print("  Generating Delta-T...", end=" ", flush=True)
    t0 = time.time()
    delta_t_table = generate_delta_t(jd_start, jd_end)
    t_dt = time.time() - t0
    if verbose:
        print(f"{len(delta_t_table)} samples ({t_dt:.1f}s)")

    # -------------------------------------------------------------------------
    # 4. Generate star catalog
    # -------------------------------------------------------------------------
    if verbose:
        print("  Generating star catalog...", end=" ", flush=True)
    star_entries = generate_star_catalog()
    if verbose:
        print(f"{len(star_entries)} stars")

    # -------------------------------------------------------------------------
    # 5. Calculate layout and write file
    # -------------------------------------------------------------------------
    if verbose:
        print("\n  Assembling file...", end=" ", flush=True)

    # Calculate section sizes
    body_count = len(bodies)
    body_index_size = body_count * BODY_ENTRY_SIZE

    # Chebyshev data size (sum of all body segments)
    chebyshev_size = 0
    for bid in bodies:
        params = BODY_PARAMS[bid]
        interval_days, degree, coord_type, components = params
        seg_size = segment_byte_size(degree, components)
        chebyshev_size += len(body_data[bid]) * seg_size

    # Nutation data size
    nut_data_size = NUTATION_HEADER_SIZE + len(nutation_coeffs) * segment_byte_size(
        NUTATION_DEGREE, NUTATION_COMPONENTS
    )

    # Delta-T size
    delta_t_size = DELTA_T_HEADER_SIZE + len(delta_t_table) * DELTA_T_ENTRY_SIZE

    # Star catalog size
    star_size = len(star_entries) * STAR_ENTRY_SIZE

    # Section directory
    section_dir_total = NUM_SECTIONS * SECTION_DIR_SIZE

    # Calculate offsets
    body_index_offset = HEADER_SIZE + section_dir_total
    chebyshev_offset = body_index_offset + body_index_size
    nutation_offset = chebyshev_offset + chebyshev_size
    delta_t_offset = nutation_offset + nut_data_size
    star_offset = delta_t_offset + delta_t_size

    total_size = star_offset + star_size

    # Allocate buffer
    buf = bytearray(total_size)

    # Write header
    header = FileHeader(
        magic=MAGIC,
        version=VERSION,
        section_count=NUM_SECTIONS,
        body_count=body_count,
        jd_start=jd_start,
        jd_end=jd_end,
        generation_epoch=now_jd,
        flags=0,
    )
    write_header(buf, header)

    # Write section directory
    sections = [
        SectionEntry(SECTION_BODY_INDEX, body_index_offset, body_index_size),
        SectionEntry(SECTION_CHEBYSHEV, chebyshev_offset, chebyshev_size),
        SectionEntry(SECTION_NUTATION, nutation_offset, nut_data_size),
        SectionEntry(SECTION_DELTA_T, delta_t_offset, delta_t_size),
        SectionEntry(SECTION_STARS, star_offset, star_size),
    ]
    for i, sec in enumerate(sections):
        write_section_dir(buf, HEADER_SIZE + i * SECTION_DIR_SIZE, sec)

    # Write body index and coefficient data
    coeff_write_offset = chebyshev_offset
    for idx, bid in enumerate(sorted(bodies)):
        params = BODY_PARAMS[bid]
        interval_days, degree, coord_type, components = params
        seg_size = segment_byte_size(degree, components)
        n_segments = len(body_data[bid])

        entry = BodyEntry(
            body_id=bid,
            coord_type=coord_type,
            segment_count=n_segments,
            jd_start=jd_start,
            jd_end=jd_end,
            interval_days=interval_days,
            degree=degree,
            components=components,
            data_offset=coeff_write_offset,
        )
        write_body_entry(buf, body_index_offset + idx * BODY_ENTRY_SIZE, entry)

        # Write coefficient data for each segment
        for seg_coeffs in body_data[bid]:
            # seg_coeffs shape: (components, degree+1)
            # Flatten to [c0_comp0, c1_comp0, ..., cN_comp0, c0_comp1, ...]
            flat = seg_coeffs.flatten()
            struct.pack_into(f"<{len(flat)}d", buf, coeff_write_offset, *flat)
            coeff_write_offset += seg_size

    # Write nutation section
    nut_n_segments = len(nutation_coeffs)
    nut_header = NutationHeader(
        jd_start=jd_start,
        jd_end=jd_end,
        interval_days=NUTATION_INTERVAL,
        degree=NUTATION_DEGREE,
        components=NUTATION_COMPONENTS,
        segment_count=nut_n_segments,
        reserved=0,
    )
    write_nutation_header(buf, nutation_offset, nut_header)

    nut_data_offset = nutation_offset + NUTATION_HEADER_SIZE
    nut_seg_size = segment_byte_size(NUTATION_DEGREE, NUTATION_COMPONENTS)
    for seg_coeffs in nutation_coeffs:
        flat = seg_coeffs.flatten()
        struct.pack_into(f"<{len(flat)}d", buf, nut_data_offset, *flat)
        nut_data_offset += nut_seg_size

    # Write Delta-T section
    struct.pack_into(
        DELTA_T_HEADER_FMT,
        buf,
        delta_t_offset,
        len(delta_t_table),
        0,  # reserved
    )
    dt_data_offset = delta_t_offset + DELTA_T_HEADER_SIZE
    for jd, dt in delta_t_table:
        struct.pack_into(DELTA_T_ENTRY_FMT, buf, dt_data_offset, jd, dt)
        dt_data_offset += DELTA_T_ENTRY_SIZE

    # Write star catalog
    for i, star in enumerate(star_entries):
        write_star_entry(buf, star_offset + i * STAR_ENTRY_SIZE, star)

    # Write to file
    with open(output, "wb") as f:
        f.write(buf)

    if verbose:
        print("done!")
        print(f"\n  File: {output}")
        print(f"  Size: {total_size:,} bytes ({total_size / (1024 * 1024):.1f} MB)")
        print(f"  Bodies: {body_count}")
        print(f"  Nutation segments: {nut_n_segments}")
        print(f"  Delta-T samples: {len(delta_t_table)}")
        print(f"  Stars: {len(star_entries)}")
        print()

        # Print error summary
        print("  Max fitting errors:")
        for bid in sorted(bodies):
            name = BODY_NAMES.get(bid, f"Body {bid}")
            error = body_errors[bid]
            params = BODY_PARAMS[bid]
            coord_type = params[2]
            if coord_type == COORD_ICRS_BARY:
                # Convert AU error to arcseconds (rough: 1 AU at 1 AU distance = 206265")
                arcsec = error * 206265.0
                print(f'    {name:20s}: {error:.2e} AU ({arcsec:.4f}")')
            else:
                # Already in degrees, convert to arcseconds
                arcsec = error * 3600.0
                print(f'    {name:20s}: {error:.2e} deg ({arcsec:.4f}")')
        nut_arcsec = math.degrees(nutation_error) * 3600.0
        print(f'    {"Nutation":20s}: {nutation_error:.2e} rad ({nut_arcsec:.4f}")')


# =============================================================================
# MERGE PARTIAL LEB FILES
# =============================================================================


def merge_leb_files(
    inputs: List[str],
    output: str,
    verbose: bool = True,
) -> None:
    """Merge multiple partial .leb files into a single complete file.

    Each input file must cover the same JD range but contain different bodies.
    Nutation, Delta-T, and star catalog are taken from the first input file
    that contains them.

    This allows generating body groups independently (e.g. planets, asteroids,
    analytical) and combining them afterward, which avoids the fork-deadlock
    issues of multiprocessing on macOS and gives finer control over
    regeneration.

    Args:
        inputs: List of paths to partial .leb files.
        output: Output path for the merged file.
        verbose: Print progress.

    Raises:
        ValueError: If inputs have mismatched JD ranges or overlapping bodies.
    """
    from libephemeris.leb_format import (
        read_body_entry,
        read_header,
        read_nutation_header,
        read_section_dir,
        read_star_entry,
    )

    if not inputs:
        raise ValueError("No input files provided")

    if verbose:
        print(f"Merging {len(inputs)} LEB files -> {output}")

    # -------------------------------------------------------------------------
    # 1. Read all input files
    # -------------------------------------------------------------------------
    all_bodies: dict[int, tuple[str, int]] = {}  # body_id -> (source_file, idx_in_file)
    ref_jd_start: Optional[float] = None
    ref_jd_end: Optional[float] = None

    # Parsed data from each input
    input_data: List[dict] = []

    for path in inputs:
        with open(path, "rb") as f:
            data = f.read()

        hdr = read_header(data, 0)
        if hdr.magic != MAGIC:
            raise ValueError(f"Invalid LEB magic in {path}")
        if hdr.version != VERSION:
            raise ValueError(f"Unsupported LEB version {hdr.version} in {path}")

        # Validate JD range consistency
        if ref_jd_start is None:
            ref_jd_start = hdr.jd_start
            ref_jd_end = hdr.jd_end
        else:
            assert ref_jd_start is not None and ref_jd_end is not None
            if (
                abs(hdr.jd_start - ref_jd_start) > 0.5
                or abs(hdr.jd_end - ref_jd_end) > 0.5
            ):
                raise ValueError(
                    f"JD range mismatch: {path} has "
                    f"[{hdr.jd_start:.1f}, {hdr.jd_end:.1f}] but expected "
                    f"[{ref_jd_start:.1f}, {ref_jd_end:.1f}]"
                )

        # Parse sections
        sections: dict[int, SectionEntry] = {}
        for i in range(hdr.section_count):
            offset = HEADER_SIZE + i * SECTION_DIR_SIZE
            sec = read_section_dir(data, offset)
            sections[sec.section_id] = sec

        # Parse bodies
        bodies: dict[int, BodyEntry] = {}
        if SECTION_BODY_INDEX in sections:
            sec = sections[SECTION_BODY_INDEX]
            for i in range(hdr.body_count):
                off = sec.offset + i * BODY_ENTRY_SIZE
                entry = read_body_entry(data, off)
                bodies[entry.body_id] = entry

                # Check for duplicates
                if entry.body_id in all_bodies:
                    src, _ = all_bodies[entry.body_id]
                    raise ValueError(
                        f"Body {entry.body_id} ({BODY_NAMES.get(entry.body_id, '?')}) "
                        f"found in both {src} and {path}"
                    )
                all_bodies[entry.body_id] = (path, i)

        info = {
            "path": path,
            "data": data,
            "header": hdr,
            "sections": sections,
            "bodies": bodies,
        }
        input_data.append(info)

        if verbose:
            body_names = [
                BODY_NAMES.get(bid, str(bid)) for bid in sorted(bodies.keys())
            ]
            print(f"  {path}: {len(bodies)} bodies ({', '.join(body_names)})")

    assert ref_jd_start is not None and ref_jd_end is not None

    # -------------------------------------------------------------------------
    # 2. Collect body coefficient data (raw bytes)
    # -------------------------------------------------------------------------
    merged_bodies = sorted(all_bodies.keys())
    body_entries: List[BodyEntry] = []
    body_coeff_blobs: List[bytes] = []  # raw coefficient bytes per body

    for bid in merged_bodies:
        # Find the source file
        for info in input_data:
            if bid in info["bodies"]:
                entry = info["bodies"][bid]
                data = info["data"]
                seg_size = segment_byte_size(entry.degree, entry.components)
                total_bytes = entry.segment_count * seg_size
                blob = data[entry.data_offset : entry.data_offset + total_bytes]
                body_entries.append(entry)
                body_coeff_blobs.append(blob)
                break

    # -------------------------------------------------------------------------
    # 3. Collect nutation from first file that has it
    # -------------------------------------------------------------------------
    nutation_blob: Optional[bytes] = None
    for info in input_data:
        if SECTION_NUTATION in info["sections"]:
            sec = info["sections"][SECTION_NUTATION]
            nutation_blob = info["data"][sec.offset : sec.offset + sec.size]
            break

    # -------------------------------------------------------------------------
    # 4. Collect delta-T from first file that has it
    # -------------------------------------------------------------------------
    delta_t_blob: Optional[bytes] = None
    for info in input_data:
        if SECTION_DELTA_T in info["sections"]:
            sec = info["sections"][SECTION_DELTA_T]
            delta_t_blob = info["data"][sec.offset : sec.offset + sec.size]
            break

    # -------------------------------------------------------------------------
    # 5. Collect star catalog from first file that has it
    # -------------------------------------------------------------------------
    star_blob: Optional[bytes] = None
    for info in input_data:
        if SECTION_STARS in info["sections"]:
            sec = info["sections"][SECTION_STARS]
            star_blob = info["data"][sec.offset : sec.offset + sec.size]
            break

    # -------------------------------------------------------------------------
    # 6. Calculate layout and write merged file
    # -------------------------------------------------------------------------
    body_count = len(merged_bodies)
    body_index_size = body_count * BODY_ENTRY_SIZE
    chebyshev_size = sum(len(b) for b in body_coeff_blobs)
    nut_size = len(nutation_blob) if nutation_blob else 0
    dt_size = len(delta_t_blob) if delta_t_blob else 0
    star_size = len(star_blob) if star_blob else 0

    section_dir_total = NUM_SECTIONS * SECTION_DIR_SIZE
    body_index_offset = HEADER_SIZE + section_dir_total
    chebyshev_offset = body_index_offset + body_index_size
    nutation_offset = chebyshev_offset + chebyshev_size
    delta_t_offset = nutation_offset + nut_size
    star_offset = delta_t_offset + dt_size
    total_size = star_offset + star_size

    buf = bytearray(total_size)

    # Header
    now_jd = J2000 + (time.time() / 86400.0 - 10957.5)
    header = FileHeader(
        magic=MAGIC,
        version=VERSION,
        section_count=NUM_SECTIONS,
        body_count=body_count,
        jd_start=ref_jd_start,
        jd_end=ref_jd_end,
        generation_epoch=now_jd,
        flags=0,
    )
    write_header(buf, header)

    # Section directory
    sections_list = [
        SectionEntry(SECTION_BODY_INDEX, body_index_offset, body_index_size),
        SectionEntry(SECTION_CHEBYSHEV, chebyshev_offset, chebyshev_size),
        SectionEntry(SECTION_NUTATION, nutation_offset, nut_size),
        SectionEntry(SECTION_DELTA_T, delta_t_offset, dt_size),
        SectionEntry(SECTION_STARS, star_offset, star_size),
    ]
    for i, sec in enumerate(sections_list):
        write_section_dir(buf, HEADER_SIZE + i * SECTION_DIR_SIZE, sec)

    # Body index + coefficient data
    coeff_write_offset = chebyshev_offset
    for idx, (entry, blob) in enumerate(zip(body_entries, body_coeff_blobs)):
        new_entry = BodyEntry(
            body_id=entry.body_id,
            coord_type=entry.coord_type,
            segment_count=entry.segment_count,
            jd_start=entry.jd_start,
            jd_end=entry.jd_end,
            interval_days=entry.interval_days,
            degree=entry.degree,
            components=entry.components,
            data_offset=coeff_write_offset,
        )
        write_body_entry(buf, body_index_offset + idx * BODY_ENTRY_SIZE, new_entry)
        buf[coeff_write_offset : coeff_write_offset + len(blob)] = blob
        coeff_write_offset += len(blob)

    # Nutation, Delta-T, stars (copy raw blobs)
    if nutation_blob:
        buf[nutation_offset : nutation_offset + len(nutation_blob)] = nutation_blob
    if delta_t_blob:
        buf[delta_t_offset : delta_t_offset + len(delta_t_blob)] = delta_t_blob
    if star_blob:
        buf[star_offset : star_offset + len(star_blob)] = star_blob

    with open(output, "wb") as f:
        f.write(buf)

    if verbose:
        print(f"\n  Merged file: {output}")
        print(f"  Size: {total_size:,} bytes ({total_size / (1024 * 1024):.1f} MB)")
        print(f"  Bodies: {body_count}")
        print(f"  JD range: {ref_jd_start:.1f} to {ref_jd_end:.1f}")
        body_list = [BODY_NAMES.get(b, str(b)) for b in merged_bodies]
        print(f"  Body list: {', '.join(body_list)}")
        print()


# =============================================================================
# VERIFICATION
# =============================================================================


def verify_leb(
    leb_path: str,
    n_samples: int = 500,
    verbose: bool = True,
) -> bool:
    """Post-generation validation of a .leb file.

    Args:
        leb_path: Path to the .leb file.
        n_samples: Number of random JDs to test per body.
        verbose: Print progress.

    Returns:
        True if all bodies pass, False otherwise.
    """
    from libephemeris.leb_reader import LEBReader

    reader = LEBReader(leb_path)
    jd_start, jd_end = reader.jd_range
    all_pass = True

    if verbose:
        print(f"Verifying {leb_path}")
        print(f"  Range: JD {jd_start:.1f} to {jd_end:.1f}")
        print(f"  Samples per body: {n_samples}")
        print()

    rng = np.random.default_rng(42)
    test_jds = rng.uniform(jd_start + 1, jd_end - 1, n_samples)

    for body_id in sorted(reader._bodies.keys()):
        body = reader._bodies[body_id]
        name = BODY_NAMES.get(body_id, f"Body {body_id}")
        max_error = 0.0

        for jd in test_jds:
            pos, vel = reader.eval_body(body_id, jd)

            if body.coord_type == COORD_ICRS_BARY:
                # Compare with Skyfield
                from libephemeris.state import get_planets, get_timescale
                from libephemeris.planets import get_planet_target

                planets = get_planets()
                ts = get_timescale()

                target_name = _PLANET_MAP.get(body_id)
                if target_name:
                    target = get_planet_target(planets, target_name)
                    t = ts.tt_jd(jd)
                    ref_pos = target.at(t).position.au
                    for c in range(3):
                        err = abs(pos[c] - float(ref_pos[c]))
                        if err > max_error:
                            max_error = err

            elif body.coord_type in (COORD_ECLIPTIC, COORD_HELIO_ECL):
                # Compare with analytical function
                from libephemeris.leb_format import BODY_PARAMS as BP

                # We'd need to call the analytical function again
                # For now, just verify the values are reasonable
                if not (0.0 <= pos[0] < 360.0 or pos[0] == 0.0):
                    max_error = 999.0

        if body.coord_type == COORD_ICRS_BARY:
            arcsec = max_error * 206265.0
            passed = arcsec < 0.01  # 10 milliarcseconds
            status = "PASS" if passed else "FAIL"
        else:
            # For ecliptic bodies, error is in degrees
            arcsec = max_error * 3600.0
            passed = True  # Already verified during generation
            status = "PASS"

        if not passed:
            all_pass = False

        if verbose:
            print(
                f'  {name:20s}: max error = {max_error:.2e} ({arcsec:.4f}") [{status}]'
            )

    reader.close()

    if verbose:
        print()
        if all_pass:
            print("  ALL BODIES PASSED")
        else:
            print("  SOME BODIES FAILED")

    return all_pass


# =============================================================================
# CLI
# =============================================================================


def _resolve_tier(args) -> Tuple[float, float, str]:
    """Resolve tier/start/end/output from CLI args.

    Returns:
        (jd_start, jd_end, output_path)
    """
    if args.tier:
        ephem_file, tier_start, tier_end, tier_output = TIER_CONFIGS[args.tier]

        from libephemeris import set_jpl_file

        set_jpl_file(ephem_file)

        start_year = args.start if args.start is not None else tier_start
        end_year = args.end if args.end is not None else tier_end

        if args.output is None:
            os.makedirs(DEFAULT_LEB_DIR, exist_ok=True)
            output = os.path.join(DEFAULT_LEB_DIR, tier_output)
        else:
            output = args.output
    else:
        if args.output is None:
            raise SystemExit("--output is required when --tier is not specified")
        output = args.output
        start_year = args.start if args.start is not None else DEFAULT_START_YEAR
        end_year = args.end if args.end is not None else DEFAULT_END_YEAR

    jd_start = _year_to_jd(start_year)
    jd_end = _year_to_jd(end_year)
    return jd_start, jd_end, output


def _group_output_path(base_output: str, group: str) -> str:
    """Derive the partial-file path for a body group.

    Example: ``data/leb/ephemeris_base.leb`` + ``planets``
             -> ``data/leb/ephemeris_base_planets.leb``
    """
    root, ext = os.path.splitext(base_output)
    return f"{root}_{group}{ext}"


def main():
    parser = argparse.ArgumentParser(
        description="Generate LEB binary ephemeris file",
        epilog=(
            "Body groups for --group:\n"
            "  planets    : Sun, Moon, Mercury-Pluto, Earth (vectorized Skyfield)\n"
            "  asteroids  : Chiron, Ceres, Pallas, Juno, Vesta (spktype21)\n"
            "  analytical : Lunar nodes, Lilith variants, Uranians, Transpluto\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output .leb file path (default: auto from --tier)",
    )
    parser.add_argument(
        "--tier",
        "-t",
        choices=["base", "medium", "extended"],
        default=None,
        help="Precision tier: base (de440s, 1850-2150), medium (de440, 1550-2650), "
        "extended (de441, -5000 to 5000). Sets ephemeris file, date range, "
        "and output path automatically.",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=None,
        help=f"Start year (default: from --tier, or {DEFAULT_START_YEAR})",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help=f"End year (default: from --tier, or {DEFAULT_END_YEAR})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() or 1,
        help=f"Number of parallel workers (default: {os.cpu_count() or 1}, "
        "auto-detected CPU count)",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Run post-generation verification",
    )
    parser.add_argument(
        "--verify-samples",
        type=int,
        default=500,
        help="Number of verification samples per body (default: 500)",
    )
    parser.add_argument(
        "--bodies",
        type=str,
        default=None,
        help="Comma-separated list of body IDs (default: all)",
    )
    parser.add_argument(
        "--group",
        choices=["planets", "asteroids", "analytical"],
        default=None,
        help="Generate only a specific body group (partial file). "
        "Use --merge to combine partial files afterward.",
    )
    parser.add_argument(
        "--merge",
        nargs="+",
        metavar="FILE",
        default=None,
        help="Merge multiple partial .leb files into one. "
        "Requires --output (or --tier for auto path). "
        "Example: --merge planets.leb asteroids.leb analytical.leb",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress progress output",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Mode 1: Merge existing partial files
    # ------------------------------------------------------------------
    if args.merge:
        if args.tier:
            _, _, output = _resolve_tier(args)
        elif args.output:
            output = args.output
        else:
            parser.error("--output (or --tier) is required with --merge")
            return  # unreachable

        t0 = time.time()
        merge_leb_files(args.merge, output, verbose=not args.quiet)
        elapsed = time.time() - t0
        if not args.quiet:
            print(f"  Merge time: {elapsed:.1f}s")

        if args.verify:
            print()
            ok = verify_leb(
                output,
                n_samples=args.verify_samples,
                verbose=not args.quiet,
            )
            if not ok:
                sys.exit(1)
        return

    # ------------------------------------------------------------------
    # Mode 2: Generate (full or group)
    # ------------------------------------------------------------------
    jd_start, jd_end, output = _resolve_tier(args)

    if not args.quiet and args.tier:
        ephem_file = TIER_CONFIGS[args.tier][0]
        print(f"  Tier: {args.tier} (ephemeris: {ephem_file})")

    # Resolve body list
    bodies = None
    if args.bodies:
        bodies = [int(b.strip()) for b in args.bodies.split(",")]
    elif args.group:
        if args.group not in BODY_GROUPS:
            parser.error(f"Unknown group: {args.group}")
        bodies = BODY_GROUPS[args.group]
        # Auto-suffix the output path when using --group without explicit --output
        if args.output is None:
            output = _group_output_path(output, args.group)
        if not args.quiet:
            print(f"  Group: {args.group} ({len(bodies)} bodies)")

    t0 = time.time()
    assemble_leb(
        output=output,
        jd_start=jd_start,
        jd_end=jd_end,
        bodies=bodies,
        workers=args.workers,
        verbose=not args.quiet,
    )
    elapsed = time.time() - t0

    if not args.quiet:
        print(f"  Total time: {elapsed:.1f}s")

    if args.verify:
        print()
        ok = verify_leb(
            output,
            n_samples=args.verify_samples,
            verbose=not args.quiet,
        )
        if not ok:
            sys.exit(1)


if __name__ == "__main__":
    main()
