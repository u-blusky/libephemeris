# Testing Guide

## Test Results Overview

When running `pytest`, you may see results like:

```
6336 passed
   4 xpassed
   3 xfailed
  57 skipped
```

### What do these mean?

| Status | Meaning |
|--------|---------|
| **passed** | Test completed successfully |
| **xpassed** | Test was expected to fail, but passed (unexpected success) |
| **xfailed** | Test was expected to fail and did fail (expected behavior) |
| **skipped** | Test was intentionally skipped due to missing dependencies |

> **Note:** `xfailed` tests are NOT failures! They document known limitations or future goals.

---

## Expected Failures (xfail)

These tests are marked as "expected to fail" for documented reasons:

| Test | Reason |
|------|--------|
| `test_deltat_future_years_match_swe` | Delta T for future dates (>2030) is a prediction - different algorithms give divergent results |
| `test_interpolated_apogee_target_precision` | Target precision <0.1° requires implementing analytical Moshier method |
| `test_interpolated_perigee_target_precision` | Target precision <0.1° requires implementing analytical Moshier method |
| Benchmark tests | Performance comparisons between pure Python and C (pyswisseph) - Python is inherently slower |

---

## Running All Tests (Without Skips)

To run the complete test suite without skips, you need to satisfy several dependencies.

### 1. Install pyswisseph

Most skipped tests require the C-based Swiss Ephemeris library for comparison:

```bash
pip install pyswisseph
```

### 2. Install Optional Dependencies

```bash
pip install requests  # Required for orbital elements update tests
```

### 3. Enable Network Tests

Some tests require network access to download data from JPL Horizons. Enable them with environment variables:

```bash
# Enable SPK auto-download tests
export LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1

# Enable SPK download tests
export LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1
```

### 4. Swiss Ephemeris Data Files

Some tests require ephemeris data files in `swisseph/ephe/`:

- `seorbel.txt` - orbital elements
- Asteroid/TNO data files

---

## Full Test Command

Run all tests with network tests enabled:

```bash
LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1 \
LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 \
pytest
```

Or set the variables permanently in your shell:

```bash
export LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1
export LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1
pytest
```

---

## Categories of Skipped Tests

### pyswisseph-dependent (majority of skips)

Tests that compare libephemeris results against pyswisseph:

- `test_lunar/test_true_lilith_latitude.py`
- `test_lunar/test_true_lilith_precision.py`
- `test_minor_bodies/test_tno_validation.py`
- `test_sol_eclipse_where_how.py`
- And many more...

### Network-dependent

Tests that download data from external sources:

| Test File | Environment Variable |
|-----------|---------------------|
| `test_spk_auto.py` | `LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1` |
| `test_spk.py` | `LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1` |
| `test_asteroid_by_number.py` | Always skipped (requires JPL SBDB) |

### Data-dependent (runtime skips)

These tests skip dynamically based on available data:

- Eclipse tests: skip if no eclipse found in search period
- Star tests: skip if star not in catalog
- TNO tests: skip if SwissEph data not available
- Historical tests: skip if date outside ephemeris range

---

## Considerations

1. **Network tests are slow** - they download data from JPL Horizons
2. **Some skips are unavoidable** - dynamic skips depend on runtime data availability
3. **pyswisseph versions vary** - some versions don't include eclipse functions
4. **Benchmark xfails are expected** - pure Python will always be slower than C

---

## Quick Reference

| Goal | Command |
|------|---------|
| Run all tests | `pytest` |
| Run with verbose output | `pytest -v` |
| Run with network tests | `LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1 LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 pytest` |
| Run specific test file | `pytest tests/test_file.py` |
| Run tests matching pattern | `pytest -k "pattern"` |
| Show skip reasons | `pytest -rs` |
| Show xfail reasons | `pytest -rx` |
