#!/usr/bin/env python3
"""
Tier diagnostic: base (de440s.bsp) -- SPK range 1850-2150.

Calculates all planets, lunar points, and minor bodies showing
ecliptic/equatorial coordinates, velocities, and data source
(DE440s / Analytical / SPK / Keplerian).

Usage:
    python scripts/tier_diagnostic_base.py
    python scripts/tier_diagnostic_base.py 1985-06-15
    python scripts/tier_diagnostic_base.py 2000-01-01 1950-03-15
    python scripts/tier_diagnostic_base.py --jd 2451545.0
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts._tier_diagnostic import run_diagnostic  # noqa: E402

if __name__ == "__main__":
    run_diagnostic("base")
