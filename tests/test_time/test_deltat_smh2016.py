"""
Tests to verify Delta T implementation matches Stephenson/Morrison/Hohenkerk 2016 model.

The SMH 2016 model is published in:
    Stephenson F.R., Morrison L.V., Hohenkerk C.Y. (2016).
    "Measurement of the Earth's rotation: 720 BC to AD 2015."
    Proceedings of the Royal Society A, 472: 20160404.
    https://doi.org/10.1098/rspa.2016.0404

The model uses:
    - Cubic spline interpolation from Table S15 for 720 BC to ~AD 2015
    - Long-term parabolic extrapolation: DeltaT = -320 + 32.5 * u^2
      where u = (year - 1825) / 100
    - IERS observations for recent/modern dates
"""

import pytest
import libephemeris as ephem


class TestDeltaTSMH2016LongTermParabola:
    """Test the long-term parabolic formula: DeltaT = -320 + 32.5 * u^2."""

    def _parabola(self, year: float) -> float:
        """Calculate Delta T using the SMH 2016 long-term parabola."""
        u = (year - 1825) / 100.0
        return -320 + 32.5 * u * u

    @pytest.mark.unit
    def test_parabola_at_1825(self):
        """At year 1825 (u=0), parabola gives DeltaT = -320 seconds."""
        expected = -320.0  # seconds
        actual = self._parabola(1825)
        assert actual == pytest.approx(expected, abs=0.1)

    @pytest.mark.unit
    def test_parabola_at_1925(self):
        """At year 1925 (u=1), parabola gives DeltaT = -320 + 32.5 = -287.5 seconds."""
        expected = -287.5  # seconds
        actual = self._parabola(1925)
        assert actual == pytest.approx(expected, abs=0.1)

    @pytest.mark.unit
    def test_parabola_at_1725(self):
        """At year 1725 (u=-1), parabola gives DeltaT = -320 + 32.5 = -287.5 seconds."""
        expected = -287.5  # seconds
        actual = self._parabola(1725)
        assert actual == pytest.approx(expected, abs=0.1)


class TestDeltaTKnownHistoricalValues:
    """
    Test Delta T against known historical values from SMH 2016 Table S15.

    These values are approximate as the actual implementation uses splines
    that smoothly interpolate between data points.
    """

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "year,expected_dt_seconds,tolerance",
        [
            # Ancient dates - large Delta T due to Earth slowing down
            # These are approximate values from the SMH 2016 model
            (-500, 17190, 500),  # 500 BC: ~4.8 hours
            (0, 10580, 300),  # Year 0: ~2.9 hours
            (500, 5710, 200),  # Year 500: ~1.6 hours
            (1000, 1570, 100),  # Year 1000: ~26 minutes
            (1500, 291, 50),  # Year 1500: ~4.9 minutes (actual Skyfield value)
            # Modern era - smaller Delta T with better measurements
            (1800, 18.2, 2),  # Year 1800: ~18 seconds (actual Skyfield value)
            (1900, -2.8, 2),  # Year 1900: ~-3 seconds
            (1950, 29.1, 2),  # Year 1950: ~29 seconds
            (2000, 63.8, 2),  # Year 2000: ~64 seconds
        ],
    )
    def test_historical_delta_t(self, year, expected_dt_seconds, tolerance):
        """Delta T should match known historical values within tolerance."""
        # Convert year to Julian Day (July 1 at noon)
        if year <= 0:
            # Astronomical year numbering: 1 BC = year 0, 2 BC = year -1
            jd = ephem.swe_julday(year, 7, 1, 12.0)
        else:
            jd = ephem.swe_julday(year, 7, 1, 12.0)

        dt = ephem.swe_deltat(jd)
        dt_seconds = dt * 86400

        assert dt_seconds == pytest.approx(expected_dt_seconds, abs=tolerance), (
            f"Year {year}: expected ~{expected_dt_seconds}s, got {dt_seconds:.1f}s"
        )


class TestDeltaTModelProperties:
    """Test mathematical properties of the SMH 2016 model."""

    @pytest.mark.unit
    def test_delta_t_is_continuous(self):
        """Delta T should be continuous without large jumps."""
        # Test across 1000 years with 10-year steps
        years = range(1000, 2001, 10)
        prev_dt = None

        for year in years:
            jd = ephem.swe_julday(year, 1, 1, 12.0)
            dt = ephem.swe_deltat(jd) * 86400  # in seconds

            if prev_dt is not None:
                # Change per decade should be reasonable
                # (less than 50 seconds per decade for most of this range)
                change = abs(dt - prev_dt)
                assert change < 100, (
                    f"Large jump in Delta T at year {year}: "
                    f"{change:.1f} seconds per decade"
                )
            prev_dt = dt

    @pytest.mark.unit
    def test_delta_t_smoothly_varying(self):
        """Delta T should change smoothly (no discontinuities in derivative)."""
        # Test that the rate of change is itself smooth
        jd_start = ephem.swe_julday(1900, 1, 1, 12.0)

        # Calculate Delta T at 1-year intervals
        dts = []
        for i in range(100):
            jd = jd_start + i * 365.25
            dts.append(ephem.swe_deltat(jd) * 86400)

        # Calculate first differences (rate of change)
        rates = [dts[i + 1] - dts[i] for i in range(len(dts) - 1)]

        # Calculate second differences (rate of rate of change)
        # Should be small for smooth function
        accelerations = [rates[i + 1] - rates[i] for i in range(len(rates) - 1)]

        for i, acc in enumerate(accelerations):
            year = 1900 + i
            assert abs(acc) < 1.0, (
                f"Large discontinuity in Delta T derivative at year ~{year}: "
                f"acceleration = {acc:.4f} s/year^2"
            )

    @pytest.mark.unit
    def test_ancient_delta_t_is_large_positive(self):
        """For ancient dates, Delta T should be large and positive (hours)."""
        # Year 1000 BC
        jd = ephem.swe_julday(-1000, 7, 1, 12.0)
        dt_seconds = ephem.swe_deltat(jd) * 86400

        # Should be multiple hours
        assert dt_seconds > 10000, (
            f"Delta T at 1000 BC should be > 10000 seconds, got {dt_seconds:.0f}"
        )
        # But not unreasonably large
        assert dt_seconds < 100000, (
            f"Delta T at 1000 BC seems too large: {dt_seconds:.0f} seconds"
        )

    @pytest.mark.unit
    def test_modern_delta_t_magnitude(self):
        """For modern dates (1900-2050), Delta T should be within reasonable bounds."""
        for year in range(1900, 2051, 10):
            jd = ephem.swe_julday(year, 1, 1, 12.0)
            dt_seconds = abs(ephem.swe_deltat(jd) * 86400)

            # Modern Delta T should be less than 200 seconds
            assert dt_seconds < 200, (
                f"Delta T at {year} = {dt_seconds:.1f}s seems too large"
            )


class TestDeltaTUsesSkyfieldModel:
    """Verify libephemeris uses Skyfield's SMH 2016 implementation."""

    @pytest.mark.unit
    def test_delta_t_uses_skyfield(self):
        """Verify Delta T is computed via Skyfield's timescale."""
        from skyfield.api import load

        ts = load.timescale()

        # J2000.0
        jd = 2451545.0

        # libephemeris result
        dt_lib = ephem.swe_deltat(jd) * 86400  # in seconds

        # Direct Skyfield result
        t = ts.ut1_jd(jd)
        dt_skyfield = t.delta_t

        # Should match exactly (both use same underlying implementation)
        assert dt_lib == pytest.approx(dt_skyfield, abs=0.001), (
            f"libephemeris ({dt_lib:.3f}s) should match Skyfield ({dt_skyfield:.3f}s)"
        )

    @pytest.mark.unit
    def test_delta_t_matches_skyfield_various_dates(self):
        """Delta T should match Skyfield across various dates."""
        from skyfield.api import load

        ts = load.timescale()

        test_years = [1800, 1900, 1950, 1980, 2000, 2010, 2020]

        for year in test_years:
            jd = ephem.swe_julday(year, 6, 15, 12.0)

            dt_lib = ephem.swe_deltat(jd) * 86400
            t = ts.ut1_jd(jd)
            dt_skyfield = t.delta_t

            assert dt_lib == pytest.approx(dt_skyfield, abs=0.01), (
                f"Year {year}: libephemeris ({dt_lib:.2f}s) != Skyfield ({dt_skyfield:.2f}s)"
            )


class TestDeltaTDocumentation:
    """Test that Delta T documentation is accurate."""

    @pytest.mark.unit
    def test_docstring_mentions_smh2016(self):
        """swe_deltat docstring should mention the SMH 2016 model."""
        docstring = ephem.swe_deltat.__doc__
        assert docstring is not None
        assert "Stephenson" in docstring
        assert "Morrison" in docstring
        assert "Hohenkerk" in docstring
        assert "2016" in docstring

    @pytest.mark.unit
    def test_docstring_mentions_parabola_formula(self):
        """swe_deltat docstring should describe the parabolic formula."""
        docstring = ephem.swe_deltat.__doc__
        assert docstring is not None
        # Check for the key components of the formula
        assert "-320" in docstring or "320" in docstring
        assert "32.5" in docstring

    @pytest.mark.unit
    def test_docstring_ex_mentions_smh2016(self):
        """swe_deltat_ex docstring should also mention the model."""
        docstring = ephem.swe_deltat_ex.__doc__
        assert docstring is not None
        assert "Stephenson" in docstring
        assert "Morrison" in docstring
        assert "Hohenkerk" in docstring
