"""Deep precision audit: fixed star calculations vs pyswisseph."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path("swisseph/ephe")

# --- Star list -----------------------------------------------------------
STARS = [
    # Bright navigational stars
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Betelgeuse",
    "Aldebaran",
    "Spica",
    "Antares",
    "Pollux",
    "Fomalhaut",
    "Deneb",
    "Regulus",
    # Astrologically important
    "Algol",
    "Alcyone",
    "Vindemiatrix",
    "Scheat",
]

# --- Test dates -----------------------------------------------------------
DATES = {
    "J2000": 2451545.0,
    "Current": 2460400.0,
    "Y1900": 2415020.5,
}

# --- Run audit ------------------------------------------------------------
print("=" * 90)
print("DEEP AUDIT: Fixed Star Precision — libephemeris vs pyswisseph")
print("=" * 90)
print()

# Column headers
dlon_h = 'dLon(")'
dlat_h = 'dLat(")'
hdr = (
    f"{'Star':<14} {'Date':<8} "
    f"{dlon_h:>10} {dlat_h:>10} {'dDist':>14} "
    f"{'Lib Lon':>12} {'Swe Lon':>12} {'Flag':>5}"
)
print(hdr)
print("-" * len(hdr))

concerns: list[str] = []
all_lon_diffs: list[float] = []
all_lat_diffs: list[float] = []
errors: list[str] = []

for star in STARS:
    for date_label, jd in DATES.items():
        # --- libephemeris ---
        try:
            lib_result = ephem.swe_fixstar_ut(star, jd, 0)
            lib_pos = lib_result[0]  # (lon, lat, dist, ...)
            lib_err = lib_result[2]  # error string
            if "could not find" in lib_err.lower():
                errors.append(f"  LIB  {star:<14} {date_label}: {lib_err}")
                continue
        except Exception as e:
            errors.append(f"  LIB  {star:<14} {date_label}: {e}")
            continue

        # --- pyswisseph ---
        try:
            swe_result = swe.fixstar_ut(star, jd, 0)
            swe_pos = swe_result[0]  # (lon, lat, dist, ...)
        except Exception as e:
            errors.append(f"  SWE  {star:<14} {date_label}: {e}")
            continue

        # --- Differences ---
        dlon = (lib_pos[0] - swe_pos[0]) * 3600.0  # degrees -> arcsec
        dlat = (lib_pos[1] - swe_pos[1]) * 3600.0
        ddist = lib_pos[2] - swe_pos[2]

        all_lon_diffs.append(abs(dlon))
        all_lat_diffs.append(abs(dlat))

        flag = ""
        if abs(dlon) > 1.0:
            flag = " ***"
            concerns.append(
                f'{star:<14} {date_label:<8}  dLon={dlon:+.4f}"  dLat={dlat:+.4f}"'
            )
        elif abs(dlon) > 0.5:
            flag = "  * "

        print(
            f"{star:<14} {date_label:<8} "
            f"{dlon:>+10.4f} {dlat:>+10.4f} {ddist:>+14.2f} "
            f"{lib_pos[0]:>12.6f} {swe_pos[0]:>12.6f} {flag:>5}"
        )

print()

# --- Summary --------------------------------------------------------------
print("=" * 90)
print("SUMMARY")
print("=" * 90)
print(f"Stars tested:       {len(STARS)}")
print(f"Dates per star:     {len(DATES)}")
print(f"Total comparisons:  {len(all_lon_diffs)}")
print()

if all_lon_diffs:
    print(
        f'Longitude  |  max: {max(all_lon_diffs):.4f}"   mean: {sum(all_lon_diffs) / len(all_lon_diffs):.4f}"   median: {sorted(all_lon_diffs)[len(all_lon_diffs) // 2]:.4f}"'
    )
if all_lat_diffs:
    print(
        f'Latitude   |  max: {max(all_lat_diffs):.4f}"   mean: {sum(all_lat_diffs) / len(all_lat_diffs):.4f}"   median: {sorted(all_lat_diffs)[len(all_lat_diffs) // 2]:.4f}"'
    )

print()
if errors:
    print(f"ERRORS ({len(errors)}):")
    for e in errors:
        print(e)
    print()

if concerns:
    print(f'CONCERNING RESULTS (>1" longitude difference): {len(concerns)}')
    for c in concerns:
        print(f"  {c}")
else:
    print('No concerning results (all longitude differences <= 1").')

# Count by threshold
over_1 = sum(1 for d in all_lon_diffs if d > 1.0)
over_05 = sum(1 for d in all_lon_diffs if d > 0.5)
over_01 = sum(1 for d in all_lon_diffs if d > 0.1)
under_01 = sum(1 for d in all_lon_diffs if d <= 0.1)

print()
print("Distribution of |dLon|:")
print(f'  > 1.0":   {over_1:>3} / {len(all_lon_diffs)}')
print(f'  > 0.5":   {over_05:>3} / {len(all_lon_diffs)}')
print(f'  > 0.1":   {over_01:>3} / {len(all_lon_diffs)}')
print(f'  <= 0.1":  {under_01:>3} / {len(all_lon_diffs)}')
print()
print('Legend:  *** = >1"   * = >0.5"')
