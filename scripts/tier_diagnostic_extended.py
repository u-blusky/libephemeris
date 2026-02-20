#!/usr/bin/env python3
"""
Tier diagnostic: extended (de441.bsp) -- SPK range 1550-2650.

Calculates all planets, lunar points, and minor bodies showing
ecliptic/equatorial coordinates, velocities, and data source
(DE441 / Analytical / SPK / Keplerian).

Usage:
    python scripts/tier_diagnostic_extended.py
    python scripts/tier_diagnostic_extended.py 1985-06-15
    python scripts/tier_diagnostic_extended.py 2000-01-01 1550-01-01 0800-01-01
    python scripts/tier_diagnostic_extended.py --jd 2451545.0
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts._tier_diagnostic import run_diagnostic  # noqa: E402

if __name__ == "__main__":
    run_diagnostic("extended")
