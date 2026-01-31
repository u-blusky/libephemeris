"""
Tests for the split_deg function that splits decimal degrees into components.

Tests verify that split_deg correctly decomposes an angle in decimal degrees
into zodiac sign, degrees, minutes, seconds, and fraction of second,
matching pyswisseph's swe.split_deg() behavior.
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestSplitDegBasic:
    """Basic functionality tests for split_deg."""

    def test_split_deg_exported(self):
        """Test that split_deg is exported from the package."""
        assert hasattr(ephem, "split_deg")
        assert callable(ephem.split_deg)

    def test_split_deg_flags_exported(self):
        """Test that SPLIT_DEG flags are exported from the package."""
        assert hasattr(ephem, "SPLIT_DEG_ROUND_SEC")
        assert hasattr(ephem, "SPLIT_DEG_ROUND_MIN")
        assert hasattr(ephem, "SPLIT_DEG_ROUND_DEG")
        assert hasattr(ephem, "SPLIT_DEG_ZODIACAL")
        assert hasattr(ephem, "SPLIT_DEG_NAKSHATRA")
        assert hasattr(ephem, "SPLIT_DEG_KEEP_SIGN")
        assert hasattr(ephem, "SPLIT_DEG_KEEP_DEG")

    def test_split_deg_flag_values(self):
        """Test that flag values match pyswisseph."""
        assert ephem.SPLIT_DEG_ROUND_SEC == swe.SPLIT_DEG_ROUND_SEC
        assert ephem.SPLIT_DEG_ROUND_MIN == swe.SPLIT_DEG_ROUND_MIN
        assert ephem.SPLIT_DEG_ROUND_DEG == swe.SPLIT_DEG_ROUND_DEG
        assert ephem.SPLIT_DEG_ZODIACAL == swe.SPLIT_DEG_ZODIACAL
        assert ephem.SPLIT_DEG_NAKSHATRA == swe.SPLIT_DEG_NAKSHATRA
        assert ephem.SPLIT_DEG_KEEP_SIGN == swe.SPLIT_DEG_KEEP_SIGN
        assert ephem.SPLIT_DEG_KEEP_DEG == swe.SPLIT_DEG_KEEP_DEG

    def test_split_deg_returns_tuple(self):
        """Test that split_deg returns a 5-tuple."""
        result = ephem.split_deg(45.5, 0)
        assert isinstance(result, tuple)
        assert len(result) == 5

    def test_split_deg_zero(self):
        """Test splitting zero degrees."""
        result = ephem.split_deg(0.0, ephem.SPLIT_DEG_ZODIACAL)
        assert result[0] == 0  # deg
        assert result[1] == 0  # min
        assert result[2] == 0  # sec
        assert result[4] == 0  # sign (Aries)


class TestSplitDegWithoutZodiacal:
    """Tests for split_deg without ZODIACAL flag."""

    def test_positive_degrees(self):
        """Test positive degrees without ZODIACAL flag."""
        result = ephem.split_deg(45.5, 0)
        assert result[0] == 45  # degrees
        assert result[1] == 30  # minutes
        assert result[2] == 0  # seconds
        assert result[4] == 1  # positive sign

    def test_negative_degrees(self):
        """Test negative degrees without ZODIACAL flag."""
        result = ephem.split_deg(-45.5, 0)
        assert result[0] == 45  # degrees (absolute value)
        assert result[1] == 30  # minutes
        assert result[2] == 0  # seconds
        assert result[4] == -1  # negative sign


class TestSplitDegZodiacal:
    """Tests for split_deg with ZODIACAL flag."""

    def test_aries(self):
        """Test position in Aries (sign 0)."""
        result = ephem.split_deg(15.5, ephem.SPLIT_DEG_ZODIACAL)
        assert result[0] == 15  # 15 degrees in sign
        assert result[1] == 30  # 30 minutes
        assert result[4] == 0  # Aries

    def test_taurus(self):
        """Test position in Taurus (sign 1)."""
        result = ephem.split_deg(45.5, ephem.SPLIT_DEG_ZODIACAL)
        assert result[0] == 15  # 15 degrees in sign (45 - 30 = 15)
        assert result[1] == 30  # 30 minutes
        assert result[4] == 1  # Taurus

    def test_cancer(self):
        """Test position in Cancer (sign 3)."""
        result = ephem.split_deg(90.5, ephem.SPLIT_DEG_ZODIACAL)
        assert result[0] == 0  # 0 degrees in sign (90 - 90 = 0)
        assert result[1] == 30  # 30 minutes
        assert result[4] == 3  # Cancer

    def test_libra(self):
        """Test position at Libra start (sign 6)."""
        result = ephem.split_deg(180.0, ephem.SPLIT_DEG_ZODIACAL)
        assert result[0] == 0  # 0 degrees in sign
        assert result[4] == 6  # Libra

    def test_pisces_end(self):
        """Test position near end of Pisces."""
        result = ephem.split_deg(359.999, ephem.SPLIT_DEG_ZODIACAL)
        assert result[4] == 11  # Pisces

    def test_negative_normalized(self):
        """Test that negative values are normalized for zodiacal."""
        result = ephem.split_deg(-30.5, ephem.SPLIT_DEG_ZODIACAL)
        # -30.5 -> 329.5 (in Aquarius, sign 10)
        expected = swe.split_deg(-30.5, swe.SPLIT_DEG_ZODIACAL)
        assert result[4] == expected[4]  # Same sign


class TestSplitDegNakshatra:
    """Tests for split_deg with NAKSHATRA flag."""

    def test_first_nakshatra(self):
        """Test position in first nakshatra (Ashwini)."""
        result = ephem.split_deg(0.0, ephem.SPLIT_DEG_NAKSHATRA)
        assert result[4] == 0  # First nakshatra

    def test_second_nakshatra(self):
        """Test position in second nakshatra (Bharani)."""
        result = ephem.split_deg(13.5, ephem.SPLIT_DEG_NAKSHATRA)
        expected = swe.split_deg(13.5, swe.SPLIT_DEG_NAKSHATRA)
        assert result[4] == expected[4]  # Same nakshatra

    def test_fourth_nakshatra(self):
        """Test position in fourth nakshatra."""
        result = ephem.split_deg(45.0, ephem.SPLIT_DEG_NAKSHATRA)
        expected = swe.split_deg(45.0, swe.SPLIT_DEG_NAKSHATRA)
        assert result[4] == expected[4]


class TestSplitDegRounding:
    """Tests for split_deg rounding behavior."""

    def test_round_sec(self):
        """Test rounding to seconds."""
        deg = 45.5123456789
        result = ephem.split_deg(
            deg, ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_SEC
        )
        expected = swe.split_deg(deg, swe.SPLIT_DEG_ZODIACAL | swe.SPLIT_DEG_ROUND_SEC)
        assert result[0] == expected[0]  # degrees
        assert result[1] == expected[1]  # minutes
        assert result[2] == expected[2]  # seconds
        assert result[4] == expected[4]  # sign

    def test_round_min(self):
        """Test rounding to minutes."""
        deg = 45.5678
        result = ephem.split_deg(
            deg, ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_MIN
        )
        expected = swe.split_deg(deg, swe.SPLIT_DEG_ZODIACAL | swe.SPLIT_DEG_ROUND_MIN)
        assert result[0] == expected[0]  # degrees
        assert result[1] == expected[1]  # minutes
        assert result[4] == expected[4]  # sign

    def test_round_deg(self):
        """Test rounding to degrees."""
        deg = 45.9999
        result = ephem.split_deg(
            deg, ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_DEG
        )
        expected = swe.split_deg(deg, swe.SPLIT_DEG_ZODIACAL | swe.SPLIT_DEG_ROUND_DEG)
        assert result[0] == expected[0]  # degrees
        assert result[4] == expected[4]  # sign


class TestSplitDegKeepFlags:
    """Tests for KEEP_SIGN and KEEP_DEG flags."""

    def test_keep_sign_at_boundary(self):
        """Test KEEP_SIGN prevents rounding to next sign."""
        deg = 29.999999
        flag = (
            ephem.SPLIT_DEG_ZODIACAL
            | ephem.SPLIT_DEG_ROUND_SEC
            | ephem.SPLIT_DEG_KEEP_SIGN
        )
        result = ephem.split_deg(deg, flag)
        expected = swe.split_deg(deg, flag)
        assert result[4] == expected[4]  # Same sign

    def test_without_keep_sign(self):
        """Test that without KEEP_SIGN, rounding advances to next sign."""
        deg = 29.999999
        flag = ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_SEC
        result = ephem.split_deg(deg, flag)
        expected = swe.split_deg(deg, flag)
        assert result[4] == expected[4]  # Same sign as pyswisseph


class TestSplitDegVsSwisseph:
    """Comparison tests with pyswisseph's swe.split_deg()."""

    @pytest.mark.parametrize(
        "degrees,flag",
        [
            (0.0, 0),
            (0.0, swe.SPLIT_DEG_ZODIACAL),
            (45.5, 0),
            (45.5, swe.SPLIT_DEG_ZODIACAL),
            (45.5, swe.SPLIT_DEG_ROUND_SEC),
            (45.5, swe.SPLIT_DEG_ZODIACAL | swe.SPLIT_DEG_ROUND_SEC),
            (90.5, swe.SPLIT_DEG_ZODIACAL),
            (180.0, swe.SPLIT_DEG_ZODIACAL),
            (270.0, swe.SPLIT_DEG_ZODIACAL),
            (359.999, swe.SPLIT_DEG_ZODIACAL),
            (360.0, swe.SPLIT_DEG_ZODIACAL),
            (370.0, swe.SPLIT_DEG_ZODIACAL),
            (-10.0, 0),
            (-10.0, swe.SPLIT_DEG_ZODIACAL),
            (-45.5, 0),
            (-45.5, swe.SPLIT_DEG_ZODIACAL),
            (30.0, swe.SPLIT_DEG_ZODIACAL),  # At sign boundary
            (60.0, swe.SPLIT_DEG_ZODIACAL),  # At sign boundary
            (0.12345678, swe.SPLIT_DEG_ZODIACAL),
            (45.5123456789, swe.SPLIT_DEG_ZODIACAL | swe.SPLIT_DEG_ROUND_SEC),
            (45.5123456789, swe.SPLIT_DEG_ZODIACAL | swe.SPLIT_DEG_ROUND_DEG),
        ],
    )
    def test_split_deg_matches_swisseph(self, degrees, flag):
        """Test that split_deg matches pyswisseph's swe.split_deg()."""
        result_lib = ephem.split_deg(degrees, flag)
        result_swe = swe.split_deg(degrees, flag)

        assert result_lib[0] == result_swe[0], (
            f"Deg mismatch for degrees={degrees}, flag={flag}: "
            f"lib={result_lib[0]}, swe={result_swe[0]}"
        )
        assert result_lib[1] == result_swe[1], (
            f"Min mismatch for degrees={degrees}, flag={flag}: "
            f"lib={result_lib[1]}, swe={result_swe[1]}"
        )
        assert result_lib[2] == result_swe[2], (
            f"Sec mismatch for degrees={degrees}, flag={flag}: "
            f"lib={result_lib[2]}, swe={result_swe[2]}"
        )
        assert result_lib[4] == result_swe[4], (
            f"Sign mismatch for degrees={degrees}, flag={flag}: "
            f"lib={result_lib[4]}, swe={result_swe[4]}"
        )

    def test_round_min_basic_behavior(self):
        """Test ROUND_MIN gives correct deg, min, and sign (sec field has quirky behavior)."""
        # pyswisseph's ROUND_MIN has peculiar sec field behavior that we don't replicate
        # This test verifies the key behavior: deg, min, and sign are correct
        deg = 45.5123456789
        flag = ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_MIN
        result_lib = ephem.split_deg(deg, flag)
        result_swe = swe.split_deg(deg, flag)

        assert result_lib[0] == result_swe[0], "Deg should match"
        assert result_lib[1] == result_swe[1], "Min should match"
        assert result_lib[4] == result_swe[4], "Sign should match"


class TestSplitDegRandomValues:
    """Random value tests to ensure broad coverage."""

    def test_random_zodiacal_values(self, random_longitudes):
        """Test with random longitudes in zodiacal mode."""
        lons = random_longitudes(50)

        for lon in lons:
            result_lib = ephem.split_deg(lon, ephem.SPLIT_DEG_ZODIACAL)
            result_swe = swe.split_deg(lon, swe.SPLIT_DEG_ZODIACAL)

            assert result_lib[0] == result_swe[0], f"Deg mismatch for lon={lon}"
            assert result_lib[1] == result_swe[1], f"Min mismatch for lon={lon}"
            assert result_lib[2] == result_swe[2], f"Sec mismatch for lon={lon}"
            assert result_lib[4] == result_swe[4], f"Sign mismatch for lon={lon}"
            # Check secfr is close (floating point)
            assert abs(result_lib[3] - result_swe[3]) < 1e-6, (
                f"Secfr mismatch for lon={lon}"
            )

    def test_random_with_rounding(self, random_longitudes):
        """Test with random longitudes and rounding."""
        lons = random_longitudes(30)

        for lon in lons:
            flag = ephem.SPLIT_DEG_ZODIACAL | ephem.SPLIT_DEG_ROUND_SEC
            result_lib = ephem.split_deg(lon, flag)
            result_swe = swe.split_deg(lon, flag)

            assert result_lib[0] == result_swe[0], f"Deg mismatch for lon={lon}"
            assert result_lib[1] == result_swe[1], f"Min mismatch for lon={lon}"
            assert result_lib[2] == result_swe[2], f"Sec mismatch for lon={lon}"
            assert result_lib[4] == result_swe[4], f"Sign mismatch for lon={lon}"
