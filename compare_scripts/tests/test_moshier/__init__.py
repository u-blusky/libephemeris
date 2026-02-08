"""
Test suite for Moshier semi-analytical ephemeris comparison with JPL/Skyfield.

This package validates:
1. Moshier vs JPL precision across overlapping date range (1550-2650)
2. Moshier range handling for dates outside JPL coverage
3. Error handling for unsupported bodies in Moshier mode
"""
