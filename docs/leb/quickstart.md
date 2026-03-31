# LEB Generation Quickstart

> Quick reference for generating and regenerating LEB1 and LEB2 files.
> For full technical details, see [guide.md](guide.md).

---

## Prerequisites

```bash
uv pip install -e ".[dev]"
```

For asteroids, JPL Horizons SPK files are needed (downloaded automatically
during generation if `LIBEPHEMERIS_AUTO_SPK=1`).

---

## 1. Generate LEB1 Files

### Via poe (recommended)

```bash
# Group workflow (recommended — avoids macOS deadlocks, allows partial regen)
poe leb:generate:base:groups       # 3 groups + merge → ephemeris_base.leb
poe leb:generate:medium:groups     # same for medium tier
poe leb:generate:extended:groups   # same for extended tier

# Or all at once (slower, no partial regen)
poe leb:generate:base              # → data/leb/ephemeris_base.leb
poe leb:generate:medium
poe leb:generate:extended
```

### Step-by-step group workflow

```bash
# 1. Generate each group separately
poe leb:generate:base:planets      # → ephemeris_base_planets.leb    (~1s)
poe leb:generate:base:asteroids    # → ephemeris_base_asteroids.leb  (~15-60s)
poe leb:generate:base:analytical   # → ephemeris_base_analytical.leb (~2-3min)

# 2. Merge + verify
poe leb:generate:base:merge        # → ephemeris_base.leb (with --verify)
```

### Regenerate a single group

If only asteroids changed (e.g. SPK update):

```bash
poe leb:generate:base:asteroids    # regenerate asteroids only
poe leb:generate:base:merge        # re-merge all 3 partial files
```

### Direct CLI

```bash
# Full tier
python scripts/generate_leb.py --tier base --verify

# Custom range
python scripts/generate_leb.py --output custom.leb --start 1900 --end 2100 --verify

# Specific bodies only
python scripts/generate_leb.py --output test.leb --start 2000 --end 2030 --bodies 0,1,2,3,4

# Single group + manual merge
python scripts/generate_leb.py --tier base --group planets
python scripts/generate_leb.py --tier base --group asteroids
python scripts/generate_leb.py --tier base --group analytical
python scripts/generate_leb.py --tier base --merge \
  data/leb/ephemeris_base_planets.leb \
  data/leb/ephemeris_base_asteroids.leb \
  data/leb/ephemeris_base_analytical.leb \
  --verify
```

### Output files

| Tier | File | Size | Coverage |
|------|------|------|----------|
| Base | `data/leb/ephemeris_base.leb` | ~102 MB | 1850-2150 |
| Medium | `data/leb/ephemeris_medium.leb` | ~377 MB | 1550-2650 |
| Extended | `data/leb/ephemeris_extended.leb` | ~2.8 GB | -5000 / +5000 |

### LEB1 body groups

| Group | Bodies | Method | Time (base) |
|-------|--------|--------|-------------|
| `planets` | Sun-Pluto, Earth (11) | Vectorized Skyfield | ~1s |
| `asteroids` | Chiron, Ceres, Pallas, Juno, Vesta (5) | spktype21 scalar | ~15-60s |
| `analytical` | Nodes, Lilith, Uranians, Apogees (15) | Scalar Python | ~2-3min |

---

## 2. Convert LEB1 → LEB2

> LEB2 requires an existing LEB1 file as input. Generate one first
> (section 1) if you don't have one.

### Via poe

```bash
# All 4 groups for a tier
poe leb2:convert:base              # → data/leb2/base_{core,asteroids,apogee,uranians}.leb

# Single group
poe leb2:convert:base:core         # → data/leb2/base_core.leb
poe leb2:convert:base:asteroids    # → data/leb2/base_asteroids.leb
poe leb2:convert:base:apogee       # → data/leb2/base_apogee.leb
poe leb2:convert:base:uranians     # → data/leb2/base_uranians.leb

# Medium / Extended
poe leb2:convert:medium
poe leb2:convert:extended
```

### Direct CLI

```bash
# Single group
python scripts/generate_leb2.py convert data/leb/ephemeris_base.leb \
  -o data/leb2/base_core.leb --group core

# All groups at once
python scripts/generate_leb2.py convert-all data/leb/ephemeris_base.leb \
  -o data/leb2/ --tier-name base
```

### Generate LEB2 from scratch (no pre-existing LEB1)

Creates a temporary LEB1 and converts it automatically:

```bash
python scripts/generate_leb2.py generate --tier base --group core \
  -o data/leb2/base_core.leb
```

### Output files (base tier)

| Group | File | Size | Bodies |
|-------|------|------|--------|
| core | `data/leb2/base_core.leb` | 7.7 MB | Sun-Pluto, Earth, Nodes, Mean Apogee (14) |
| asteroids | `data/leb2/base_asteroids.leb` | 7.7 MB | Chiron, Ceres, Pallas, Juno, Vesta (5) |
| apogee | `data/leb2/base_apogee.leb` | 10.3 MB | OscuApog, IntpApog, IntpPerig (3) |
| uranians | `data/leb2/base_uranians.leb` | 2.0 MB | Cupido-Transpluto (9) |
| **Total** | | **28.0 MB** | **31 bodies** |

### LEB2 groups vs LEB1 groups

LEB2 groups differ from LEB1 groups:

| LEB1 | LEB2 | Difference |
|------|------|------------|
| `planets` (11) | `core` (14) | Core also includes Mean/True Node and Mean Apogee |
| `asteroids` (5) | `asteroids` (5) | Identical |
| `analytical` (15) | `apogee` (3) + `uranians` (9) | Analytical split into 2 groups |

---

## 3. Verification and Testing

### Verify LEB1

```bash
# Verification built into generation (add --verify)
python scripts/generate_leb.py --tier base --verify --verify-samples 1000

# Quick manual smoke test
python3 -c "
from libephemeris.leb_reader import LEBReader
reader = LEBReader('data/leb/ephemeris_base.leb')
pos, vel = reader.eval_body(0, 2451545.0)  # Sun at J2000
print(f'Sun ICRS: x={pos[0]:.10f} y={pos[1]:.10f} z={pos[2]:.10f} AU')
reader.close()
"
```

### Verify LEB2

```bash
# Against LEB1 reference
poe leb2:verify:base

# Direct CLI
python scripts/generate_leb2.py verify data/leb2/base_core.leb \
  --reference data/leb/ephemeris_base.leb --samples 500
```

### Test suites

```bash
# LEB1
poe test:leb                       # Unit tests (no @slow)
poe test:leb:precision             # Full precision suite
poe test:leb:precision:quick       # Medium tier only

# LEB2
poe test:leb2                      # Compression round-trip + reader tests
poe test:leb2:precision:base       # End-to-end precision, base (~15s)
poe test:leb2:precision:all        # All tiers (~45s)

# With LEB active (any format)
poe test:unit:leb                  # Unit tests in LEB mode
poe test:unit:leb:fast             # Same, parallel
```

### Individual test files

```bash
pytest tests/test_leb/test_leb_format.py -v
pytest tests/test_leb/test_leb_reader.py -v
pytest tests/test_leb/test_fast_calc.py -v
pytest tests/test_leb/test_leb_compression.py -v
pytest tests/test_leb/test_leb2_reader.py -v
```

---

## 4. Release

```bash
# 1. Generate the LEB files
poe leb:generate:medium:groups

# 2. Dry run
poe release:leb:dry-run 0.22.0

# 3. Upload + update hashes
poe release:leb:medium 0.22.0

# 4. Commit
git add libephemeris/download.py
git commit -m "update LEB medium tier hash after regeneration"
```

---

## 5. Runtime Usage

```python
from libephemeris import set_leb_file, set_calc_mode

# LEB1 (single file)
set_leb_file("data/leb/ephemeris_base.leb")

# LEB2 (auto-discovers companion files)
set_leb_file("data/leb2/base_core.leb")

# Or via environment variable
# export LIBEPHEMERIS_LEB=/path/to/file.leb

# Force LEB mode (error if no file available)
set_calc_mode("leb")
```

Auto-discovery: files in `~/.libephemeris/leb/` are found automatically.

```bash
# End-user download (no poe needed)
libephemeris download:leb:base       # ~53 MB
libephemeris download:leb:medium     # ~175 MB
```

---

## Command Summary

| Action | Command |
|--------|---------|
| **LEB1 base (recommended)** | `poe leb:generate:base:groups` |
| **LEB1 medium** | `poe leb:generate:medium:groups` |
| **LEB1 extended** | `poe leb:generate:extended:groups` |
| **LEB2 base (all groups)** | `poe leb2:convert:base` |
| **LEB2 base core only** | `poe leb2:convert:base:core` |
| **LEB2 medium** | `poe leb2:convert:medium` |
| **Verify LEB2** | `poe leb2:verify:base` |
| **Test LEB1** | `poe test:leb` |
| **Test LEB2** | `poe test:leb2` |
| **Test LEB2 precision** | `poe test:leb2:precision:base` |
| **Release** | `poe release:leb 0.22.0` |
