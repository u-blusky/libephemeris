#!/usr/bin/env python3
"""
Sidereal Calculations Example for libephemeris.

This script demonstrates how to work with sidereal (Vedic) astrology:
1. Setting different ayanamsha modes
2. Calculating sidereal planetary positions
3. Comparing tropical vs sidereal positions
4. Working with popular Indian ayanamshas

Requirements:
    pip install libephemeris

Usage:
    python examples/sidereal_calculations.py
"""

from __future__ import annotations

import libephemeris as eph
from libephemeris.constants import (
    # Planet IDs
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
    SE_TRUE_NODE,
    # Calculation flags
    SEFLG_SWIEPH,
    SEFLG_SPEED,
    SEFLG_SIDEREAL,
    # Sidereal modes (ayanamshas)
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_LAHIRI,
    SE_SIDM_RAMAN,
    SE_SIDM_USHASHASHI,
    SE_SIDM_KRISHNAMURTI,
    SE_SIDM_DJWHAL_KHUL,
    SE_SIDM_YUKTESHWAR,
    SE_SIDM_JN_BHASIN,
    SE_SIDM_BABYL_KUGLER1,
    SE_SIDM_TRUE_CITRA,
    SE_SIDM_TRUE_REVATI,
    SE_SIDM_TRUE_PUSHYA,
    SE_SIDM_GALCENT_0SAG,
)


# Zodiac signs - same for both tropical and sidereal
SIGNS = [
    "Aries",
    "Taurus",
    "Gemini",
    "Cancer",
    "Leo",
    "Virgo",
    "Libra",
    "Scorpio",
    "Sagittarius",
    "Capricorn",
    "Aquarius",
    "Pisces",
]

# Vedic names for signs (nakshatras use different divisions)
RASHIS = [
    "Mesha",
    "Vrishabha",
    "Mithuna",
    "Karka",
    "Simha",
    "Kanya",
    "Tula",
    "Vrischika",
    "Dhanu",
    "Makara",
    "Kumbha",
    "Meena",
]


def format_zodiac(longitude: float, use_vedic_names: bool = False) -> str:
    """Convert ecliptic longitude to zodiac notation."""
    signs = RASHIS if use_vedic_names else SIGNS
    sign_num = int(longitude / 30) % 12
    pos_in_sign = longitude % 30
    degrees = int(pos_in_sign)
    minutes = int((pos_in_sign - degrees) * 60)
    return f"{degrees:02d}° {minutes:02d}' {signs[sign_num]}"


def example_ayanamsha_overview() -> None:
    """Example 1: Overview of available ayanamshas."""
    print("=" * 60)
    print("Example 1: Ayanamsha Overview")
    print("=" * 60)

    print("\nThe ayanamsha is the angular difference between the tropical")
    print("and sidereal zodiacs, caused by the precession of equinoxes.")
    print("\nPopular ayanamsha systems:\n")

    jd = eph.swe_julday(2024, 1, 1, 0.0)

    # List of common ayanamshas with descriptions
    ayanamshas = [
        (SE_SIDM_LAHIRI, "Lahiri", "Official Indian government standard"),
        (SE_SIDM_RAMAN, "Raman", "B.V. Raman's ayanamsha"),
        (SE_SIDM_KRISHNAMURTI, "Krishnamurti", "KP astrology system"),
        (SE_SIDM_FAGAN_BRADLEY, "Fagan-Bradley", "Western sidereal standard"),
        (SE_SIDM_TRUE_CITRA, "True Citra", "Spica at 0° Libra"),
        (SE_SIDM_YUKTESHWAR, "Yukteshwar", "Sri Yukteshwar's calculation"),
        (SE_SIDM_JN_BHASIN, "JN Bhasin", "J.N. Bhasin's ayanamsha"),
        (SE_SIDM_USHASHASHI, "Ushashashi", "Based on Surya Siddhanta"),
        (SE_SIDM_DJWHAL_KHUL, "Djwhal Khul", "Theosophical tradition"),
        (SE_SIDM_TRUE_REVATI, "True Revati", "Zeta Piscium at 29°50' Pisces"),
        (SE_SIDM_TRUE_PUSHYA, "True Pushya", "Delta Cancri at 16° Cancer"),
        (SE_SIDM_BABYL_KUGLER1, "Babylonian (Kugler 1)", "Ancient Babylonian"),
        (SE_SIDM_GALCENT_0SAG, "Galactic Center", "Galactic center at 0° Sagittarius"),
    ]

    print(f"{'Ayanamsha':<25} {'Value':>12}  Description")
    print("-" * 65)

    for mode, name, desc in ayanamshas:
        eph.swe_set_sid_mode(mode)
        aya = eph.swe_get_ayanamsa_ut(jd)
        print(f"{name:<25} {aya:>11.6f}°  {desc}")


def example_tropical_vs_sidereal() -> None:
    """Example 2: Compare tropical and sidereal positions."""
    print("\n" + "=" * 60)
    print("Example 2: Tropical vs Sidereal Positions")
    print("=" * 60)

    jd = eph.swe_julday(2024, 3, 21, 12.0)  # Spring equinox 2024

    print(f"\nComparison at JD {jd} (Spring Equinox 2024):\n")

    # Set Lahiri ayanamsha
    eph.swe_set_sid_mode(SE_SIDM_LAHIRI)
    ayanamsha = eph.swe_get_ayanamsa_ut(jd)
    print(f"Lahiri Ayanamsha: {ayanamsha:.6f}°\n")

    planets = [
        (SE_SUN, "Sun"),
        (SE_MOON, "Moon"),
        (SE_MARS, "Mars"),
        (SE_JUPITER, "Jupiter"),
        (SE_SATURN, "Saturn"),
        (SE_TRUE_NODE, "Rahu"),
    ]

    print(f"{'Planet':<12} {'Tropical':>22} {'Sidereal (Lahiri)':>22}")
    print("-" * 58)

    for planet_id, name in planets:
        # Tropical position
        pos_trop, _ = eph.swe_calc_ut(jd, planet_id, SEFLG_SWIEPH | SEFLG_SPEED)

        # Sidereal position
        pos_sid, _ = eph.swe_calc_ut(
            jd, planet_id, SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_SIDEREAL
        )

        trop_str = format_zodiac(pos_trop[0])
        sid_str = format_zodiac(pos_sid[0], use_vedic_names=True)

        print(f"{name:<12} {trop_str:>22} {sid_str:>22}")

    print("\nNote: Ketu (South Node) is 180° opposite from Rahu")


def example_vedic_chart() -> None:
    """Example 3: Calculate a Vedic (Sidereal) natal chart."""
    print("\n" + "=" * 60)
    print("Example 3: Vedic Natal Chart Calculation")
    print("=" * 60)

    # Example birth data: Mumbai, India
    year, month, day, hour = 1990, 5, 15, 10.5
    latitude = 19.0760  # Mumbai latitude
    longitude = 72.8777  # Mumbai longitude

    jd = eph.swe_julday(year, month, day, hour)

    print(
        f"\nBirth data: {year}-{month:02d}-{day:02d} at {int(hour):02d}:{int((hour % 1) * 60):02d}"
    )
    print(f"Location: Mumbai ({latitude}°N, {longitude}°E)")

    # Set Lahiri ayanamsha (most common in India)
    eph.swe_set_sid_mode(SE_SIDM_LAHIRI)
    ayanamsha = eph.swe_get_ayanamsa_ut(jd)
    print(f"Ayanamsha (Lahiri): {ayanamsha:.6f}°\n")

    # Calculate sidereal planets
    print("Grahas (Planets):")
    print("-" * 40)

    planets = [
        (SE_SUN, "Surya (Sun)"),
        (SE_MOON, "Chandra (Moon)"),
        (SE_MERCURY, "Budha (Mercury)"),
        (SE_VENUS, "Shukra (Venus)"),
        (SE_MARS, "Mangal (Mars)"),
        (SE_JUPITER, "Guru (Jupiter)"),
        (SE_SATURN, "Shani (Saturn)"),
        (SE_TRUE_NODE, "Rahu"),
    ]

    flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_SIDEREAL

    for planet_id, name in planets:
        pos, _ = eph.swe_calc_ut(jd, planet_id, flags)
        formatted = format_zodiac(pos[0], use_vedic_names=True)
        retro = " (Vakri)" if pos[3] < 0 else ""
        print(f"  {name:20} {formatted}{retro}")

    # Calculate Ketu (opposite of Rahu)
    rahu_pos, _ = eph.swe_calc_ut(jd, SE_TRUE_NODE, flags)
    ketu_lon = (rahu_pos[0] + 180) % 360
    print(f"  {'Ketu':20} {format_zodiac(ketu_lon, use_vedic_names=True)}")

    # Calculate houses (Bhava) - using Whole Sign houses (common in Vedic)
    print("\nBhava Sphutas (House Cusps) - Whole Sign:")
    print("-" * 40)

    # Get ascendant in sidereal
    cusps, ascmc = eph.swe_houses_ex(jd, latitude, longitude, ord("W"), SEFLG_SIDEREAL)

    for i in range(12):
        print(f"  Bhava {i + 1:2}: {format_zodiac(cusps[i], use_vedic_names=True)}")

    print(f"\n  Lagna (Ascendant): {format_zodiac(ascmc[0], use_vedic_names=True)}")


def example_ayanamsha_change_over_time() -> None:
    """Example 4: Show ayanamsha change over centuries."""
    print("\n" + "=" * 60)
    print("Example 4: Ayanamsha Change Over Time")
    print("=" * 60)

    print("\nThe ayanamsha increases by approximately 50.3 arcseconds per year")
    print("due to the precession of the Earth's axis.\n")

    eph.swe_set_sid_mode(SE_SIDM_LAHIRI)

    years = [1900, 1950, 2000, 2024, 2050, 2100]

    print(f"{'Year':>6}  {'Ayanamsha':>12}  {'Change from 2000':>18}")
    print("-" * 40)

    jd_2000 = eph.swe_julday(2000, 1, 1, 0.0)
    aya_2000 = eph.swe_get_ayanamsa_ut(jd_2000)

    for year in years:
        jd = eph.swe_julday(year, 1, 1, 0.0)
        aya = eph.swe_get_ayanamsa_ut(jd)
        change = aya - aya_2000
        sign = "+" if change >= 0 else ""
        print(f"{year:>6}  {aya:>11.6f}°  {sign}{change:>17.6f}°")

    # Rate of precession
    jd1 = eph.swe_julday(2024, 1, 1, 0.0)
    jd2 = eph.swe_julday(2025, 1, 1, 0.0)
    aya1 = eph.swe_get_ayanamsa_ut(jd1)
    aya2 = eph.swe_get_ayanamsa_ut(jd2)
    annual_precession = aya2 - aya1

    print(f"\nAnnual precession rate: {annual_precession * 3600:.2f} arcseconds")
    print(f"                        ({annual_precession:.6f}°)")


def example_custom_ayanamsha() -> None:
    """Example 5: Set a custom ayanamsha."""
    print("\n" + "=" * 60)
    print("Example 5: Custom Ayanamsha")
    print("=" * 60)

    from libephemeris.constants import SE_SIDM_USER

    jd = eph.swe_julday(2024, 1, 1, 0.0)

    print("\nYou can define a custom ayanamsha by specifying:")
    print("- Reference epoch (Julian Day)")
    print("- Ayanamsha value at that epoch")
    print("- (Optionally) precession rate\n")

    # Set a custom ayanamsha: 23.5° at J2000
    j2000 = 2451545.0  # J2000.0 epoch
    custom_aya = 23.5  # Custom value at J2000

    # SE_SIDM_USER (255) enables user-defined ayanamsha
    eph.swe_set_sid_mode(SE_SIDM_USER, j2000, custom_aya)

    aya = eph.swe_get_ayanamsa_ut(jd)
    print("Custom ayanamsha (23.5° at J2000):")
    print(f"  Value at 2024-01-01: {aya:.6f}°")

    # Compare with standard Lahiri
    eph.swe_set_sid_mode(SE_SIDM_LAHIRI)
    lahiri = eph.swe_get_ayanamsa_ut(jd)
    print("\nLahiri for comparison:")
    print(f"  Value at 2024-01-01: {lahiri:.6f}°")


def example_nakshatra_calculation() -> None:
    """Example 6: Calculate Nakshatra position."""
    print("\n" + "=" * 60)
    print("Example 6: Nakshatra Calculation")
    print("=" * 60)

    # Nakshatras (27 lunar mansions of 13°20' each)
    NAKSHATRAS = [
        "Ashwini",
        "Bharani",
        "Krittika",
        "Rohini",
        "Mrigashira",
        "Ardra",
        "Punarvasu",
        "Pushya",
        "Ashlesha",
        "Magha",
        "Purva Phalguni",
        "Uttara Phalguni",
        "Hasta",
        "Chitra",
        "Swati",
        "Vishakha",
        "Anuradha",
        "Jyeshtha",
        "Mula",
        "Purva Ashadha",
        "Uttara Ashadha",
        "Shravana",
        "Dhanishta",
        "Shatabhisha",
        "Purva Bhadrapada",
        "Uttara Bhadrapada",
        "Revati",
    ]

    jd = eph.swe_julday(2024, 1, 1, 12.0)

    # Set Lahiri ayanamsha
    eph.swe_set_sid_mode(SE_SIDM_LAHIRI)

    print(f"\nNakshatra positions at JD {jd}:\n")

    flags = SEFLG_SWIEPH | SEFLG_SIDEREAL

    planets = [
        (SE_SUN, "Sun"),
        (SE_MOON, "Moon"),
        (SE_MARS, "Mars"),
        (SE_MERCURY, "Mercury"),
        (SE_JUPITER, "Jupiter"),
        (SE_VENUS, "Venus"),
        (SE_SATURN, "Saturn"),
    ]

    print(f"{'Planet':<10} {'Nakshatra':<20} {'Pada':>6} {'Position':>12}")
    print("-" * 52)

    for planet_id, name in planets:
        pos, _ = eph.swe_calc_ut(jd, planet_id, flags)
        longitude = pos[0]

        # Calculate nakshatra (each is 13°20' = 13.3333°)
        nakshatra_size = 360 / 27  # 13.3333...
        nakshatra_num = int(longitude / nakshatra_size)
        nakshatra_name = NAKSHATRAS[nakshatra_num]

        # Calculate pada (each nakshatra has 4 padas of 3°20' each)
        pos_in_nakshatra = longitude - (nakshatra_num * nakshatra_size)
        pada = int(pos_in_nakshatra / (nakshatra_size / 4)) + 1

        print(f"{name:<10} {nakshatra_name:<20} {pada:>6} {longitude:>11.4f}°")


def main() -> None:
    """Run all sidereal calculation examples."""
    print("\n" + "#" * 60)
    print("# libephemeris Sidereal Calculations Examples")
    print("#" * 60)

    example_ayanamsha_overview()
    example_tropical_vs_sidereal()
    example_vedic_chart()
    example_ayanamsha_change_over_time()
    example_custom_ayanamsha()
    example_nakshatra_calculation()

    print("\n" + "#" * 60)
    print("# Sidereal examples completed!")
    print("#" * 60 + "\n")


if __name__ == "__main__":
    main()
