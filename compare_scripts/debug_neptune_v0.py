"""Debug Neptune magnitude calculation with direct V0 inspection."""

from __future__ import annotations

import os

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import importlib
import libephemeris
import libephemeris.planets as planets_mod

importlib.reload(planets_mod)

import swisseph as swe
import math

swe.set_ephe_path("swisseph/ephe")

# Check _PLANET_MAG_PARAMS
print("_PLANET_MAG_PARAMS:", planets_mod._PLANET_MAG_PARAMS)
print()

# Check if SE_NEPTUNE is in it
from libephemeris.constants import SE_NEPTUNE

print(f"SE_NEPTUNE = {SE_NEPTUNE}")
print(
    f"SE_NEPTUNE in _PLANET_MAG_PARAMS: {SE_NEPTUNE in planets_mod._PLANET_MAG_PARAMS}"
)
print()

# Check what _calc_planet_magnitude does for Neptune
# Test with tjd for year 2000 and 1970
_J2000 = 2451545.0

for year_label, tjd in [("2000", _J2000), ("1979", 2444000.5), ("1970", 2440587.5)]:
    result = planets_mod._calc_planet_magnitude(
        ipl=SE_NEPTUNE,
        helio_dist=30.0,
        geo_dist=30.0,
        phase_angle=1.0,
        tjd=tjd,
    )
    expected_year = 2000.0 + (tjd - _J2000) / 365.25
    print(
        f"Year ~{year_label} (tjd={tjd}, computed_year={expected_year:.1f}): magnitude={result:.4f}"
    )
