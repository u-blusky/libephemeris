"""
Keplerian precision benchmark: measures accuracy of the Keplerian fallback
against SPK reference positions for the 5 major asteroids.

This test serves as a baseline for evaluating improvements to the Keplerian
propagation (secular perturbations, multi-epoch elements, etc.).

Results are printed as a table of errors in arcseconds at various time
offsets from the orbital element epoch.

Usage:
    pytest tests/test_keplerian_precision_benchmark.py -v -s
"""

from __future__ import annotations

import math
from typing import Tuple

import pytest

from libephemeris.constants import (
    SE_CERES,
    SE_CHIRON,
    SE_JUNO,
    SE_PALLAS,
    SE_VESTA,
)
from libephemeris.minor_bodies import (
    MINOR_BODY_ELEMENTS,
    _get_closest_epoch_elements,
    calc_minor_body_position,
)


# The 5 major asteroids with SPK support
BENCHMARK_BODIES = [
    (SE_CERES, "Ceres"),
    (SE_PALLAS, "Pallas"),
    (SE_JUNO, "Juno"),
    (SE_VESTA, "Vesta"),
    (SE_CHIRON, "Chiron"),
]

# Time offsets from epoch to test (in years)
TIME_OFFSETS_YEARS = [0.0, 0.083, 0.25, 0.5, 1, 2, 5, 10, 25, 50, 100]

# Element epoch
EPOCH_JD = 2461000.5  # 2025-Sep-19 TDB


def _angular_separation_arcsec(
    lon1: float, lat1: float, lon2: float, lat2: float
) -> float:
    """Compute angular separation in arcseconds between two ecliptic positions.

    Uses the Vincenty formula for better accuracy at small separations.
    """
    lon1_r = math.radians(lon1)
    lat1_r = math.radians(lat1)
    lon2_r = math.radians(lon2)
    lat2_r = math.radians(lat2)

    dlon = lon2_r - lon1_r
    cos_lat1 = math.cos(lat1_r)
    cos_lat2 = math.cos(lat2_r)
    sin_lat1 = math.sin(lat1_r)
    sin_lat2 = math.sin(lat2_r)

    # Vincenty formula
    num = math.sqrt(
        (cos_lat2 * math.sin(dlon)) ** 2
        + (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * math.cos(dlon)) ** 2
    )
    den = sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * math.cos(dlon)

    sep_rad = math.atan2(num, den)
    return math.degrees(sep_rad) * 3600.0


def _get_spk_heliocentric_j2000(
    body_id: int, jd_tt: float
) -> Tuple[float, float, float]:
    """Get heliocentric ecliptic J2000 (lon, lat, dist) via spktype21.

    Returns (lon_deg, lat_deg, dist_au) in ecliptic J2000 frame,
    matching the frame used by calc_minor_body_position().

    The SPK type 21 file provides heliocentric positions in the J2000
    equatorial frame (ICRS). We rotate to ecliptic J2000 using the
    standard obliquity.
    """
    from libephemeris.state import _SPK_BODY_MAP

    if body_id not in _SPK_BODY_MAP:
        raise RuntimeError(f"No SPK registered for body {body_id}")

    spk_file, naif_id = _SPK_BODY_MAP[body_id]

    from spktype21 import SPKType21

    kernel = SPKType21.open(spk_file)
    try:
        center_id = kernel.segments[0].center
        target_id = kernel.segments[0].target

        AU_KM = 149597870.7
        pos_km, _ = kernel.compute_type21(center_id, target_id, jd_tt)

        # SPK gives heliocentric ICRS (J2000 equatorial) in km
        x_eq = pos_km[0] / AU_KM
        y_eq = pos_km[1] / AU_KM
        z_eq = pos_km[2] / AU_KM

        # Rotate from equatorial J2000 to ecliptic J2000
        # Using the standard J2000 obliquity: 23.4392911 degrees
        eps = math.radians(23.4392911)
        cos_eps = math.cos(eps)
        sin_eps = math.sin(eps)

        x_ecl = x_eq
        y_ecl = y_eq * cos_eps + z_eq * sin_eps
        z_ecl = -y_eq * sin_eps + z_eq * cos_eps

        # Convert to spherical
        r = math.sqrt(x_ecl**2 + y_ecl**2 + z_ecl**2)
        lon = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
        lat = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0

        return lon, lat, r
    finally:
        kernel.close()


def _get_keplerian_heliocentric(
    body_id: int, jd_tt: float
) -> Tuple[float, float, float]:
    """Get heliocentric ecliptic (lon, lat, dist) via pure Keplerian."""
    elements = _get_closest_epoch_elements(body_id, jd_tt)
    x, y, z = calc_minor_body_position(elements, jd_tt, body_id=body_id)

    r = math.sqrt(x**2 + y**2 + z**2)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(z / r)) if r > 0 else 0.0

    return lon, lat, r


def _ensure_spk_available(body_id: int) -> bool:
    """Ensure SPK is downloaded and registered for a body. Returns True on success."""
    try:
        from libephemeris.minor_bodies import auto_download_asteroid_spk

        import libephemeris

        libephemeris.set_auto_spk_download(True)

        # Download SPK for a wide range
        auto_download_asteroid_spk(
            body_id,
            jd_start=EPOCH_JD - 100 * 365.25,
            jd_end=EPOCH_JD + 100 * 365.25,
            force=False,
        )
        from libephemeris.state import _SPK_BODY_MAP

        return body_id in _SPK_BODY_MAP
    except Exception:
        return False


@pytest.mark.slow
class TestKeplerianPrecisionBenchmark:
    """Benchmark Keplerian fallback accuracy against SPK reference."""

    @pytest.fixture(autouse=True)
    def setup_spk(self):
        """Ensure SPK kernels are available for all benchmark bodies."""
        self.available_bodies = []
        for body_id, name in BENCHMARK_BODIES:
            if _ensure_spk_available(body_id):
                self.available_bodies.append((body_id, name))

        if not self.available_bodies:
            pytest.skip("No SPK kernels available for benchmark")

    def test_keplerian_precision_table(self):
        """Print a comprehensive precision table.

        Measures Keplerian error vs SPK at multiple time offsets for each body.
        This is a diagnostic test — it always passes but prints the results.
        """
        print("\n" + "=" * 90)
        print("KEPLERIAN PRECISION BENCHMARK")
        print(f"Epoch: JD {EPOCH_JD:.1f} (2025-Sep-19)")
        print("Reference: JPL SPK type 21 (sub-arcsecond accuracy)")
        print("=" * 90)

        # Header
        header = f"{'Body':>10s}"
        for yr in TIME_OFFSETS_YEARS:
            if yr == 0:
                header += f"  {'epoch':>8s}"
            elif yr < 1:
                header += f"  {f'{yr * 12:.0f}mo':>8s}"
            else:
                header += f"  {f'{yr:.0f}yr':>8s}"
        print(header)
        print("-" * len(header))

        # Results dict for assertions
        results: dict[int, dict[float, float]] = {}

        for body_id, name in self.available_bodies:
            results[body_id] = {}
            row = f"{name:>10s}"

            for yr_offset in TIME_OFFSETS_YEARS:
                dt_days = yr_offset * 365.25
                # Test both forward and backward
                errors = []
                for sign in [1, -1]:
                    jd = EPOCH_JD + sign * dt_days
                    try:
                        spk_lon, spk_lat, spk_dist = _get_spk_heliocentric_j2000(
                            body_id, jd
                        )
                        kep_lon, kep_lat, kep_dist = _get_keplerian_heliocentric(
                            body_id, jd
                        )
                        err = _angular_separation_arcsec(
                            spk_lon, spk_lat, kep_lon, kep_lat
                        )
                        errors.append(err)
                    except Exception:
                        pass

                if errors:
                    max_err = max(errors)
                    results[body_id][yr_offset] = max_err

                    # Format based on magnitude
                    if max_err < 1.0:
                        row += f'  {max_err:7.3f}"'
                    elif max_err < 60.0:
                        row += f'  {max_err:7.1f}"'
                    elif max_err < 3600.0:
                        row += f"  {max_err / 60.0:6.1f}'"
                    else:
                        row += f"  {max_err / 3600.0:5.1f}deg"
                else:
                    row += "      N/A"

            print(row)

        print("-" * len(header))
        print("Units: \" = arcseconds, ' = arcminutes, deg = degrees")
        print()

        # Summary statistics
        print("SUMMARY (max error across all bodies at each offset):")
        for yr_offset in TIME_OFFSETS_YEARS:
            all_errors = [
                results[bid].get(yr_offset, 0)
                for bid in results
                if yr_offset in results[bid]
            ]
            if all_errors:
                max_all = max(all_errors)
                if yr_offset == 0:
                    label = "at epoch"
                elif yr_offset < 1:
                    label = f"at {yr_offset * 12:.0f} months"
                else:
                    label = f"at {yr_offset:.0f} years"

                if max_all < 1.0:
                    print(f'  {label:>15s}: {max_all:.3f}"')
                elif max_all < 60.0:
                    print(f'  {label:>15s}: {max_all:.1f}"')
                elif max_all < 3600.0:
                    print(f"  {label:>15s}: {max_all / 60.0:.1f}'")
                else:
                    print(f"  {label:>15s}: {max_all / 3600.0:.1f} deg")

        print()

    def test_keplerian_near_epoch_accuracy(self):
        """Verify Keplerian is accurate within 1" at epoch."""
        for body_id, name in self.available_bodies:
            jd = EPOCH_JD
            try:
                spk_lon, spk_lat, _ = _get_spk_heliocentric_j2000(body_id, jd)
                kep_lon, kep_lat, _ = _get_keplerian_heliocentric(body_id, jd)
                err = _angular_separation_arcsec(spk_lon, spk_lat, kep_lon, kep_lat)
                assert err < 1.0, (
                    f'{name} at epoch: Keplerian error {err:.3f}" exceeds 1"'
                )
            except RuntimeError:
                pass  # SPK unavailable

    def test_keplerian_1month_accuracy(self):
        """Verify Keplerian is accurate within 60" at ±1 month."""
        for body_id, name in self.available_bodies:
            for sign in [1, -1]:
                jd = EPOCH_JD + sign * 30.0
                try:
                    spk_lon, spk_lat, _ = _get_spk_heliocentric_j2000(body_id, jd)
                    kep_lon, kep_lat, _ = _get_keplerian_heliocentric(body_id, jd)
                    err = _angular_separation_arcsec(spk_lon, spk_lat, kep_lon, kep_lat)
                    assert err < 60.0, (
                        f'{name} at ±1mo: Keplerian error {err:.1f}" exceeds 60"'
                    )
                except RuntimeError:
                    pass

    def test_keplerian_1year_accuracy(self):
        """Verify Keplerian is within 5' at ±1 year."""
        for body_id, name in self.available_bodies:
            for sign in [1, -1]:
                jd = EPOCH_JD + sign * 365.25
                try:
                    spk_lon, spk_lat, _ = _get_spk_heliocentric_j2000(body_id, jd)
                    kep_lon, kep_lat, _ = _get_keplerian_heliocentric(body_id, jd)
                    err = _angular_separation_arcsec(spk_lon, spk_lat, kep_lon, kep_lat)
                    assert err < 300.0, (
                        f"{name} at ±1yr: Keplerian error {err:.1f}\" exceeds 5'"
                    )
                except RuntimeError:
                    pass

    def test_distance_error(self):
        """Measure distance (AU) error at various offsets."""
        print("\n" + "=" * 70)
        print("DISTANCE ERROR (AU)")
        print("=" * 70)

        header = f"{'Body':>10s}"
        for yr in [0, 1, 5, 10, 50, 100]:
            if yr == 0:
                header += f"  {'epoch':>10s}"
            else:
                header += f"  {f'{yr}yr':>10s}"
        print(header)
        print("-" * len(header))

        for body_id, name in self.available_bodies:
            row = f"{name:>10s}"
            for yr_offset in [0, 1, 5, 10, 50, 100]:
                dt_days = yr_offset * 365.25
                errors = []
                for sign in [1, -1]:
                    jd = EPOCH_JD + sign * dt_days
                    try:
                        _, _, spk_dist = _get_spk_heliocentric_j2000(body_id, jd)
                        _, _, kep_dist = _get_keplerian_heliocentric(body_id, jd)
                        errors.append(abs(spk_dist - kep_dist))
                    except Exception:
                        pass
                if errors:
                    max_err = max(errors)
                    row += f"  {max_err:10.6f}"
                else:
                    row += "        N/A"
            print(row)

        print()
