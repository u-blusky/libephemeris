# LEB vs Skyfield — Guida ai Test di Comparazione

Questa guida spiega come eseguire e interpretare i test che confrontano i risultati del calcolo LEB (LibEphemeris Binary) con quelli di Skyfield (calcolo diretto da effemeridi NASA JPL).

## Prerequisiti

### File LEB

I test richiedono file `.leb` pre-generati. Ci sono due posizioni possibili:

| Posizione | Path | Note |
|-----------|------|------|
| **Locale** (default) | `data/leb/ephemeris_{tier}.leb` | Usata automaticamente dai test |
| **Esterna** | `/Volumes/Data/libephemeris/leb/ephemeris_{tier}.leb` | Va indicata via env var |

Per usare un file LEB esterno (ad esempio su un disco separato):

```bash
export LIBEPHEMERIS_LEB=/Volumes/Data/libephemeris/leb/ephemeris_medium.leb
```

Se il file LEB non viene trovato, i test vengono **skippati** (non falliscono).

### Tier disponibili

| Tier | Effemeridi | Range | File |
|------|-----------|-------|------|
| `base` | de440s.bsp | 1850–2150 | `ephemeris_base.leb` (~94 MB) |
| `medium` | de440.bsp | 1550–2650 | `ephemeris_medium.leb` (~315 MB) |
| `extended` | de441.bsp | -13200–+17191 | `ephemeris_extended.leb` |

### Dipendenze

```bash
uv pip install -e ".[dev]"
```

I test degli asteroidi richiedono accesso alla rete per scaricare automaticamente i file SPK21 da JPL Horizons (gestito da `set_auto_spk_download(True)` nel setup dei test).

---

## Comandi rapidi

### Medium tier (default)

```bash
# Tutti i test medium (completo, ~90 secondi con -n 4)
poe test:leb:compare

# Solo i test veloci (no @pytest.mark.slow)
poe test:leb:compare:quick

# Con parallelismo (raccomandato)
LIBEPHEMERIS_LEB=/path/to/ephemeris_medium.leb \
  pytest tests/test_leb/compare/ -m "leb_compare" -v --tb=short -n 4
```

### Base tier

```bash
# Tutti i test base (~90 secondi con -n 4)
poe test:leb:compare:base

# Solo veloci
poe test:leb:compare:base:quick

# Con parallelismo
pytest tests/test_leb/compare/base/ -m "leb_compare_base" -v --tb=short -n 4
```

### Extended tier

```bash
poe test:leb:compare:extended
poe test:leb:compare:extended:quick
```

### Tutti i tier insieme

```bash
poe test:leb:compare:all
```

### Test singolo

```bash
# Un file specifico
pytest tests/test_leb/compare/test_compare_leb_planets.py -v -n 4

# Un test specifico
pytest tests/test_leb/compare/test_compare_leb_planets.py::TestPlanetLongitude::test_longitude[6-Saturn] -v

# Per keyword
pytest tests/test_leb/compare/ -m "leb_compare" -k "asteroid" -v -n 4
```

---

## Struttura dei test

### Directory

```
tests/test_leb/compare/
├── conftest.py                          # Infrastruttura condivisa (tolleranze, helper, fixture)
├── test_compare_leb_planets.py          # Lon, lat, dist, speed dei pianeti ICRS
├── test_compare_leb_asteroids.py        # Posizione, speed, distanza asteroidi
├── test_compare_leb_hypothetical.py     # Corpi Uraniani (Cupido, Hades, Zeus, ...)
├── test_compare_leb_velocities.py       # Speed lon/lat/dist per tutti i 30 corpi
├── test_compare_leb_distances.py        # Distanza geocentrica e eliocentrica
├── test_compare_leb_crossings.py        # swe_cross_ut, swe_solcross_ut, ...
├── test_compare_leb_eclipses_solar.py   # Eclissi solari
├── test_compare_leb_eclipses_lunar.py   # Eclissi lunari
├── test_compare_leb_nutation.py         # Nutazione
├── test_compare_leb_deltat.py           # Delta-T
├── test_compare_leb_ayanamsha.py        # Ayanamsha (27 modi siderali)
├── test_compare_leb_sidereal.py         # Posizioni siderali
├── test_compare_leb_observations.py     # Coordinate equatoriali, J2000
├── test_compare_leb_houses.py           # Case (Placidus, Koch, ...)
├── test_compare_leb_rise_transit.py     # Alba, tramonto, transito
├── test_compare_leb_stations.py         # Stazioni (retrogradazione)
├── test_compare_leb_elongation.py       # Elongazione
├── test_compare_leb_gauquelin.py        # Settori Gauquelin
├── test_compare_leb_lunar.py            # Funzioni lunari specifiche
├── base/                                # Test base tier (de440s)
│   ├── conftest.py                      # Tolleranze e fixture base tier
│   ├── test_base_planets.py
│   ├── test_base_asteroids.py
│   ├── test_base_velocities.py
│   ├── test_base_distances.py
│   ├── test_base_hypothetical.py
│   ├── test_base_sidereal.py
│   ├── test_base_flags.py
│   └── test_base_lunar.py
├── extended/                            # Test extended tier (de441)
└── crosstier/                           # Test di consistenza cross-tier
```

### Marker pytest

| Marker | Significato |
|--------|-------------|
| `leb_compare` | Test medium tier |
| `leb_compare_base` | Test base tier |
| `leb_compare_extended` | Test extended tier |
| `leb_compare_crosstier` | Test cross-tier |
| `slow` | Test intensivi (100-200 date per corpo) |

---

## Come funziona il confronto

### CompareHelper

Il cuore dell'infrastruttura è la classe `CompareHelper` in `conftest.py`. Esegue la stessa funzione in due modalità:

```python
# Modalità Skyfield (riferimento): calcolo diretto da effemeridi NASA
ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

# Modalità LEB: calcolo via Chebyshev precomputed
leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
```

Internamente, `CompareHelper`:
1. **Salva** lo stato globale di libephemeris (mode, tier, LEB file)
2. **skyfield()**: forza `set_calc_mode("skyfield")` e il precision tier corretto
3. **leb()**: forza il file LEB specificato e `set_calc_mode("auto")`
4. **teardown()**: ripristina lo stato originale

### Pattern tipico di un test

```python
@pytest.mark.leb_compare
@pytest.mark.slow
@pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
def test_longitude(self, compare, test_dates_200, body_id, body_name):
    max_err = 0.0
    worst_jd = 0.0

    for jd in test_dates_200:
        ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
        leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

        err = lon_error_arcsec(ref[0], leb[0])
        if err > max_err:
            max_err = err
            worst_jd = jd

    assert max_err < TOLS.POSITION_ARCSEC, (
        f'{body_name}: max lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
    )
```

Il test:
1. Itera su N date uniformemente distribuite nel range del tier
2. Per ogni data, calcola con Skyfield e con LEB
3. Misura l'errore massimo (worst case)
4. Confronta con la tolleranza del tier

### Date di test

Le date sono generate uniformemente nel range del tier con un margine di 30 giorni ai bordi:

| Fixture | N date | Uso |
|---------|--------|-----|
| `test_dates_200` | 200 | Posizione pianeti, asteroidi |
| `test_dates_100` | 100 | Velocità, distanze |
| `test_dates_50` | 50 | Equatoriali, siderali |
| `test_dates_20` | 20 | Test rapidi |

Per il base tier: `base_dates_300`, `base_dates_150`, `base_dates_100`, `base_dates_50`.

### Filtro date asteroidi

Gli asteroidi hanno dati SPK validi solo per ~1900-2100 CE. Le date vengono filtrate automaticamente:

```python
dates = filter_asteroid_dates(test_dates_200, body_id)
```

Per corpi non-asteroidi, restituisce tutte le date invariate.

---

## Tolleranze

### Struttura

Le tolleranze sono definite nella dataclass `TierTolerances` e configurate per tier in `TIER_DEFAULTS`:

```python
TOLS = TierTolerances.for_tier("medium")  # Carica i default del medium tier
```

### Override via variabili d'ambiente

```bash
# Override per un tier specifico
LEB_TOL_BASE_POSITION_ARCSEC=10.0 pytest ...

# Override globale (fallback per tutti i tier)
LEB_TOL_POSITION_ARCSEC=10.0 pytest ...
```

Ordine di priorità (dal più alto al più basso):
1. Override esplicito via kwargs
2. `LEB_TOL_{TIER}_{FIELD}` (env var per-tier)
3. `LEB_TOL_{FIELD}` (env var globale)
4. `TIER_DEFAULTS[tier]` (default nel codice)
5. Default della dataclass

### Tolleranze correnti

#### Posizione

| Campo | Base | Medium | Unità | Note |
|-------|------|--------|-------|------|
| `POSITION_ARCSEC` | 5.0 | 5.0 | arcsec | Saturn/Uranus (limite architetturale) |
| `ASTEROID_ARCSEC` | 0.5 | 0.5 | arcsec | |
| `ECLIPTIC_ARCSEC` | 0.05 | 0.05 | arcsec | Nodi, Lilith |
| `HYPOTHETICAL_ARCSEC` | 0.001 | 0.001 | arcsec | Uraniani (errore ~0) |
| `EQUATORIAL_ARCSEC` | 1.0 | 0.5 | arcsec | |
| `J2000_ARCSEC` | 1.0 | 0.5 | arcsec | |
| `SIDEREAL_ARCSEC` | 5.0 | 5.0 | arcsec | |
| `DISTANCE_AU` | 3e-5 | 3e-5 | AU | |

#### Velocità

| Campo | Base | Medium | Unità | Note |
|-------|------|--------|-------|------|
| `SPEED_LON_DEG_DAY` | 0.045 | 0.045 | deg/day | OscuApogee domina |
| `SPEED_LAT_DEG_DAY` | 0.005 | 0.004 | deg/day | |
| `SPEED_DIST_AU_DAY` | 1e-4 | 3e-5 | AU/day | |
| `ASTEROID_SPEED_LON_DEG_DAY` | — | 0.001 | deg/day | |
| `ASTEROID_SPEED_LAT_DEG_DAY` | 0.75 | 0.40 | deg/day | Limite architetturale |
| `ASTEROID_SPEED_DIST_AU_DAY` | — | 1e-6 | AU/day | |

#### Timing (funzioni indirette)

| Campo | Valore | Unità | Funzione testata |
|-------|--------|-------|-----------------|
| `CROSSING_SUN_SEC` | 1.0 | sec | `swe_solcross_ut` |
| `CROSSING_MOON_SEC` | 5.0 | sec | `swe_mooncross_ut` |
| `CROSSING_PLANET_SEC` | 30.0 | sec | `swe_cross_ut` |
| `ECLIPSE_TIMING_SEC` | 1.0 | sec | `swe_sol_eclipse_*` |
| `STATION_TIMING_SEC` | 1.0 | sec | Stazioni retrograde |
| `RISE_TRANSIT_SEC` | 1.0 | sec | `swe_rise_trans` |

---

## Corpi testati (31 totali)

### Pipeline A — ICRS (11 corpi)

| ID | Nome | Note |
|----|------|------|
| 0 | Sun | |
| 1 | Moon | |
| 2 | Mercury | |
| 3 | Venus | |
| 4 | Mars | |
| 5 | Jupiter | |
| 6 | Saturn | |
| 7 | Uranus | |
| 8 | Neptune | |
| 9 | Pluto | |
| 14 | Earth | |

### Pipeline A — Asteroidi (5 corpi)

| ID | Nome | Note |
|----|------|------|
| 15 | Chiron | SPK valido solo 1900–2100 |
| 17 | Ceres | SPK valido solo 1900–2100 |
| 18 | Pallas | SPK valido solo 1900–2100, orbita inclinata 34.8° |
| 19 | Juno | SPK valido solo 1900–2100 |
| 20 | Vesta | SPK valido solo 1900–2100 |

### Pipeline B — Eclittici (6 corpi)

| ID | Nome | Note |
|----|------|------|
| 10 | MeanNode | Errore ~0 (formula analitica) |
| 11 | TrueNode | |
| 12 | MeanApogee | Errore ~0 (formula analitica) |
| 13 | OscuApogee | Errore maggiore in velocità |
| 21 | InterpApogee | |
| 22 | InterpPerigee | |

### Pipeline B — Ipotetici/Uraniani (9 corpi)

| ID | Nome |
|----|------|
| 40 | Cupido |
| 41 | Hades |
| 42 | Zeus |
| 43 | Kronos |
| 44 | Apollon |
| 45 | Admetos |
| 46 | Vulkanus |
| 47 | Poseidon |
| 48 | Transpluto |

---

## Rigenerazione dei file LEB

Se modifichi i parametri Chebyshev o la pipeline di calcolo, devi rigenerare i file LEB.

### Generazione per gruppi (raccomandato)

Su macOS usare sempre la generazione per gruppi (evita deadlock del multiprocessing):

```bash
# Base tier
poe leb:generate:base:groups

# Medium tier
poe leb:generate:medium:groups

# Extended tier
poe leb:generate:extended:groups
```

Ogni comando esegue in sequenza:
1. `planets` — Sun-Pluto, Earth (11 corpi)
2. `asteroids` — Chiron, Ceres, Pallas, Juno, Vesta (5 corpi)
3. `analytical` — Nodi, Lilith, Uraniani (15 corpi)
4. `merge` — Unisce i 3 file parziali + verifica

### Generazione singola (Linux)

```bash
poe leb:generate:base
poe leb:generate:medium
poe leb:generate:extended
```

### Copia su disco esterno

Dopo la generazione, il file è in `data/leb/`. Per usarlo dai test con env var:

```bash
cp data/leb/ephemeris_medium.leb /Volumes/Data/libephemeris/leb/
```

---

## Test xfailed e skipped

### xfail (fallimenti attesi)

| Test | Motivo |
|------|--------|
| Jupiter/Saturn crossing geocentrico | Bug pre-esistente nel solver `crossing.py` (non LEB) |
| Mars 180° crossing geocentrico | `RuntimeError: Maximum iterations reached` (non LEB) |
| Saturn 180°/270° crossing eliocentrico | `RuntimeError: Heliocentric crossing search diverged` (non LEB) |

### skip

I test vengono skippati se il file LEB del tier corrispondente non è trovato.

---

## Interpretare i fallimenti

### Errore di posizione (arcsec)

```
AssertionError: Saturn: max lon error = 6.2345" at JD 2451544.5
assert 6.2345 < 5.0
```

L'errore è in **arcosecondi**. Per convertire:
- In gradi: dividere per 3600 (6.23" = 0.0017°)
- In minuti d'arco: dividere per 60 (6.23" = 0.10')

### Errore di velocità (deg/day)

```
AssertionError: Pallas: max lat speed error = 0.450000 deg/day at JD 2451544.5
assert 0.45 < 0.40
```

### Errore di distanza (AU)

```
AssertionError: Pluto: max dist error = 4.50e-05 AU at JD 2451544.5
assert 4.5e-05 < 3e-05
```

1 AU = ~150 milioni di km. Un errore di 3e-5 AU = ~4500 km.

### Cosa fare se un test fallisce

1. **Controllare la data** — errori grandi a date estreme (inizio/fine del tier) suggeriscono contaminazione SPK o segmenti Chebyshev di bordo
2. **Controllare il corpo** — Saturn, Uranus, asteroidi hanno limiti architetturali noti
3. **Allargare la tolleranza?** — solo se l'errore è vicino al limite e il margine di sicurezza (2x) è confermato su molte date
4. **Rigenerare il LEB?** — se hai cambiato parametri in `leb_format.py` o `fast_calc.py`
