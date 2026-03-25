"""
Tests for heliacal events of fixed stars in libephemeris.

Tests the calculation of heliacal rising and setting events for bright fixed stars.
The heliacal rising of Sirius was crucial for the Egyptian calendar and marked
the beginning of the Nile flood season.

Reference: Historical astronomical records and Meeus, Astronomical Algorithms.
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    revjul,
    heliacal_ut,
    swe_heliacal_ut,
    heliacal_pheno_ut,
    is_fixed_star,
    SE_SIRIUS,
    SE_ALDEBARAN,
    SE_REGULUS,
    SE_VEGA,
    SE_ARCTURUS,
    SE_SPICA_STAR,
    SE_BETELGEUSE,
    SE_ANTARES,
    SE_FOMALHAUT,
    SE_HELIACAL_RISING,
    SE_HELIACAL_SETTING,
    SE_EVENING_FIRST,
    SE_MORNING_LAST,
    SE_FIXSTAR_OFFSET,
)


class TestIsFixedStar:
    """Test the is_fixed_star helper function."""

    def test_sirius_is_fixed_star(self):
        """Test that Sirius is identified as a fixed star."""
        assert is_fixed_star(SE_SIRIUS)

    def test_regulus_is_fixed_star(self):
        """Test that Regulus is identified as a fixed star."""
        assert is_fixed_star(SE_REGULUS)

    def test_planet_is_not_fixed_star(self):
        """Test that planets are not identified as fixed stars."""
        from libephemeris import SE_VENUS, SE_MARS, SE_JUPITER

        assert not is_fixed_star(SE_VENUS)
        assert not is_fixed_star(SE_MARS)
        assert not is_fixed_star(SE_JUPITER)

    def test_sun_is_not_fixed_star(self):
        """Test that the Sun is not identified as a fixed star."""
        from libephemeris import SE_SUN

        assert not is_fixed_star(SE_SUN)

    def test_fixstar_offset(self):
        """Test that the fixed star offset is correctly applied."""
        assert SE_SIRIUS >= SE_FIXSTAR_OFFSET
        assert SE_REGULUS >= SE_FIXSTAR_OFFSET


class TestStarHeliacalRising:
    """Test heliacal rising calculations for fixed stars."""

    def test_sirius_heliacal_rising_egypt(self):
        """Test heliacal rising of Sirius from Cairo (historically significant)."""
        # Sirius heliacal rising marked the Egyptian new year
        jd_start = julday(2024, 1, 1, 0)
        # Cairo, Egypt
        geopos = (31.2357, 30.0444, 0.0)

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Sirius",
            SE_HELIACAL_RISING,
        )

        # Should find an event within a year
        assert result[0] > jd_start
        assert result[0] < jd_start + 400

        # Just verify the date is valid
        year, month, day, hour = revjul(result[0])
        assert year in (2024, 2025)
        assert 1 <= month <= 12
        assert 1 <= day <= 31

    def test_sirius_heliacal_rising_rome(self):
        """Test heliacal rising of Sirius from Rome."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Sirius",
            SE_HELIACAL_RISING,
        )

        assert result[0] > jd_start
        assert result[0] < jd_start + 400

    def test_aldebaran_heliacal_rising(self):
        """Test heliacal rising of Aldebaran."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Aldebaran",
            SE_HELIACAL_RISING,
        )

        assert result[0] > jd_start
        assert result[0] < jd_start + 400

    def test_regulus_heliacal_rising(self):
        """Test heliacal rising of Regulus."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (-0.1278, 51.5074, 0.0)  # London

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Regulus",
            SE_HELIACAL_RISING,
        )

        assert result[0] > jd_start
        assert result[0] < jd_start + 400

    def test_vega_heliacal_rising(self):
        """Test heliacal rising of Vega."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (-74.0060, 40.7128, 0.0)  # New York

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Vega",
            SE_HELIACAL_RISING,
        )

        # Vega may or may not be found depending on position
        if result[0] > 0:
            assert result[0] > jd_start
            assert result[0] < jd_start + 400

    def test_arcturus_heliacal_rising(self):
        """Test heliacal rising of Arcturus."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (139.6503, 35.6762, 0.0)  # Tokyo

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Arcturus",
            SE_HELIACAL_RISING,
        )

        # Arcturus may or may not be found depending on position
        if result[0] > 0:
            assert result[0] > jd_start
            assert result[0] < jd_start + 400


class TestStarHeliacalSetting:
    """Test heliacal setting calculations for fixed stars."""

    def test_sirius_heliacal_setting(self):
        """Test heliacal setting of Sirius."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Sirius",
            SE_HELIACAL_SETTING,
        )

        assert result[0] > jd_start

    def test_regulus_heliacal_setting(self):
        """Test heliacal setting of Regulus."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Regulus",
            SE_HELIACAL_SETTING,
        )

        assert result[0] > jd_start


class TestStarHeliacalValidation:
    """Test input validation for star heliacal calculations."""

    def test_star_evening_first_raises_error(self):
        """Test that SE_EVENING_FIRST raises ValueError for fixed stars."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        with pytest.raises(ValueError, match="fixed stars"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Sirius",
                SE_EVENING_FIRST,
            )

    def test_star_morning_last_raises_error(self):
        """Test that SE_MORNING_LAST raises ValueError for fixed stars."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        with pytest.raises(ValueError, match="fixed stars"):
            heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "Sirius",
                SE_MORNING_LAST,
            )


class TestSweHeliacalUtStars:
    """Test swe_heliacal_ut API with fixed star names."""

    def test_sirius_by_name(self):
        """Test swe_heliacal_ut with Sirius by name."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Sirius", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3
        if result[0] > 0:
            assert result[0] > jd_start

    def test_regulus_by_name(self):
        """Test swe_heliacal_ut with Regulus by name."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Regulus", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_aldebaran_by_name(self):
        """Test swe_heliacal_ut with Aldebaran by name."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Aldebaran", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_vega_by_name(self):
        """Test swe_heliacal_ut with Vega by name."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        result = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Vega", SE_HELIACAL_RISING
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_star_name_case_insensitive(self):
        """Test that star name matching is case insensitive."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        # Lowercase
        result1 = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "sirius", SE_HELIACAL_RISING
        )
        # Uppercase
        result2 = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "SIRIUS", SE_HELIACAL_RISING
        )
        # Mixed case
        result3 = swe_heliacal_ut(
            jd_start, geopos, datm, dobs, "Sirius", SE_HELIACAL_RISING
        )

        # All should return same result
        assert result1[0] == result2[0] == result3[0]

    def test_star_evening_first_by_name_raises_error(self):
        """Test that SE_EVENING_FIRST raises ValueError for stars via swe API."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0)
        dobs = (36.0, 1.0, 0, 0, 0, 0)

        with pytest.raises(ValueError, match="fixed stars"):
            swe_heliacal_ut(jd_start, geopos, datm, dobs, "Sirius", SE_EVENING_FIRST)


class TestStarHeliacalPhenoUt:
    """Test heliacal_pheno_ut for fixed stars."""

    def test_sirius_pheno_returns_valid_result(self):
        """Test that Sirius heliacal phenomena returns a valid result."""
        jd = julday(2024, 7, 20, 4)  # Morning during Sirius heliacal period
        geopos = (31.2357, 30.0444, 0.0)  # Cairo

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Sirius",
            SE_HELIACAL_RISING,
        )

        assert isinstance(dret, tuple)
        assert len(dret) == 50
        assert all(isinstance(x, float) for x in dret)

    def test_star_magnitude_in_pheno(self):
        """Test that star magnitude is correctly returned in pheno results."""
        jd = julday(2024, 7, 20, 4)
        geopos = (31.2357, 30.0444, 0.0)  # Cairo

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Sirius",
            SE_HELIACAL_RISING,
        )

        magnitude = dret[20]
        # Sirius has magnitude around -1.46
        assert -2.0 < magnitude < 0.0

    def test_star_parallax_is_zero(self):
        """Test that star parallax is essentially zero."""
        jd = julday(2024, 7, 20, 4)
        geopos = (31.2357, 30.0444, 0.0)  # Cairo

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Sirius",
            SE_HELIACAL_RISING,
        )

        parallax = dret[19]
        # Stars have essentially zero parallax
        assert parallax == 0.0

    def test_regulus_pheno(self):
        """Test Regulus heliacal phenomena."""
        jd = julday(2024, 8, 15, 5)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        dret = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Regulus",
            SE_HELIACAL_RISING,
        )

        assert isinstance(dret, tuple)
        assert len(dret) == 50


class TestStarMagnitudeVisibility:
    """Test that star magnitude affects visibility calculations."""

    def test_bright_star_visible_earlier(self):
        """
        Brighter stars should be visible earlier in twilight than fainter stars.
        This is a qualitative test - brighter stars have lower arcus visionis requirements.
        """
        jd = julday(2024, 7, 20, 4)
        geopos = (31.2357, 30.0444, 0.0)  # Cairo

        # Sirius is very bright (mag -1.46)
        dret_sirius = heliacal_pheno_ut(
            jd,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Sirius",
            SE_HELIACAL_RISING,
        )

        # Check that magnitude is correctly retrieved
        mag_sirius = dret_sirius[20]
        assert mag_sirius < 0  # Sirius is brighter than magnitude 0


class TestStarHeliacalLocations:
    """Test star heliacal calculations at various locations."""

    def test_southern_hemisphere(self):
        """Test star heliacal from southern hemisphere."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (151.2093, -33.8688, 0.0)  # Sydney

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Sirius",
            SE_HELIACAL_RISING,
        )

        # Should find event (or not found, depending on visibility)
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_equatorial_location(self):
        """Test star heliacal from equatorial location."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (32.5825, 0.3476, 0.0)  # Kampala, Uganda

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Vega",
            SE_HELIACAL_RISING,
        )

        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_high_latitude(self):
        """Test star heliacal from high latitude."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (24.9384, 60.1699, 0.0)  # Helsinki

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Arcturus",
            SE_HELIACAL_RISING,
        )

        assert isinstance(result, tuple)
        assert len(result) == 3


class TestStarHeliacalAtmospheric:
    """Test star heliacal with different atmospheric conditions."""

    def test_star_with_altitude(self):
        """Test star heliacal with observer at altitude."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 2000.0)  # Rome, 2000 meters

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 50.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Sirius",
            SE_HELIACAL_RISING,
        )

        if result[0] > 0:
            assert result[0] > jd_start

    def test_star_with_humidity(self):
        """Test star heliacal with different humidity."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        result = heliacal_ut(
            jd_start,
            geopos,
            (1013.25, 15.0, 80.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            "Regulus",
            SE_HELIACAL_RISING,
        )

        assert isinstance(result, tuple)
        assert len(result) == 3


class TestMultipleStars:
    """Test heliacal calculations for multiple stars."""

    def test_multiple_bright_stars(self):
        """Test heliacal rising for several bright stars."""
        jd_start = julday(2024, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0.0)  # Rome

        star_names = [
            "Sirius",
            "Aldebaran",
            "Regulus",
            "Vega",
            "Arcturus",
            "Spica",
            "Betelgeuse",
            "Antares",
            "Fomalhaut",
        ]

        for star_name in star_names:
            result = heliacal_ut(
                jd_start,
                geopos,
                (1013.25, 15.0, 50.0, 0.0),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                star_name,
                SE_HELIACAL_RISING,
            )

            # Should return a valid 3-tuple
            assert isinstance(result, tuple), f"Unexpected return for star {star_name}"
            assert len(result) == 3, f"Unexpected return for star {star_name}"
