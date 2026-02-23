# Lunar Apsides: The Scientific Approach of LibEphemeris

One of the most complex challenges in computational astronomy is determining the position of the lunar apsides (the perigee and apogee, also known in astrology as True/Black Moon Lilith).

Unlike the planets, the Moon's orbit is heavily deformed by the Sun's gravity. This means the instantaneous (osculating) perigee is highly volatile, swinging wildly by up to ±30° over the course of a single month.

To provide a useful value for daily calculations, ephemeris software must compute an "Interpolated" or "Natural" perigee—a mathematically smoothed curve that removes this short-period volatility while remaining physically accurate to the Moon's actual orbital dynamics.

Here, `libephemeris` deliberately departs from the legacy methodology used by Swiss Ephemeris (`pyswisseph`) to prioritize modern scientific and physical rigor.

## The Legacy Approach (Swiss Ephemeris)

In the 1990s, the authors of Swiss Ephemeris faced the challenge of computing this smoothed curve without the computing power available today.

Their solution was based on **ELP2000-82B**, an analytical lunar theory developed by Chapront-Touzé and Chapront in the 1980s. ELP2000 uses thousands of trigonometric terms to approximate the Moon's position. To create a "smooth" perigee, the developers manually isolated specific terms in the theory and mathematically "switched off" any term that contained the mean lunar anomaly (the Moon's monthly motion).

While this creates a mathematically smooth curve, **it is an artificial mathematical construct**. Because it relies on a truncated analytical theory from 1988, this curve can deviate from the *actual physical geometry* of the Earth-Moon system by up to **5 degrees**.

## The Modern Scientific Approach (LibEphemeris)

Instead of relying on truncated 1980s analytical theories, `libephemeris` leverages modern computing power and the **NASA JPL DE440/441** numerical integrations (the same models used to navigate spacecraft).

Our approach to the Interpolated Perigee is physically grounded:

1. **Physical Reality:** We calculate the exact moments in time when the Moon physically reaches its closest point to Earth (the actual perigee passage) using the JPL DE440/441 state vectors. At these exact moments, the volatile osculating perigee perfectly aligns with physical reality.
2. **Mathematical Smoothing:** We use modern mathematical interpolation (splines and harmonic series fits calibrated over thousands of years) to draw a smooth, continuous curve directly through these true physical passages.

### The Result

The result is a Natural Perigee that accurately reflects the physical reality of the solar system as defined by modern space agencies, rather than an artifact of a truncated mathematical theory.

**Because of this commitment to scientific accuracy, the `libephemeris` Interpolated Perigee (`SE_INTP_PERG`) will differ from the Swiss Ephemeris output by up to 5 degrees.**

This is not a bug or a lack of precision. It is a deliberate enhancement. We consider the JPL numerical integrations to be the absolute ground truth of the solar system, and we prioritize physical reality over backwards compatibility with legacy analytical approximations.