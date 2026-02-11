# PLAN: Rimozione Moshier + Supporto DE441

## Obiettivo

Semplificare libephemeris eliminando le effemeridi Moshier (VSOP87, ELP2000-82B, Pluto analytical) e impostando tutto su JPL DE440/DE441 a prescindere dai flag. Aggiungere supporto per DE441 via variabile d'ambiente per coprire il range temporale esteso.

## Decisioni di design

| Aspetto | Decisione |
|---------|-----------|
| `SEFLG_MOSEPH` | Accettato e **ignorato silenziosamente** (compatibilità API pyswisseph) |
| Costanti `SEFLG_MOSEPH` / `FLG_MOSEPH` | **Mantenute** in `constants.py` e negli export |
| Fallback DE421 | **Rimosso** — chi vuole un range diverso usa DE441 esplicitamente |
| Fallback kepleriano (minor bodies) | **Non toccato** — resta com'è |
| `strict_precision` | **Non toccato** — resta com'è |
| Auto-SPK download | **Non toccato** — resta com'è |
| Pianeti uraniani/ipotietici | **Non toccati** — funzionano solo via kepleriano |
| Range default | DE440 1550–2650 CE — **invariato** |
| Range esteso | DE441 -13200 – +17191 CE via env var `LIBEPHEMERIS_EPHEMERIS` |

## Impatto

| Aspetto | Prima | Dopo |
|---------|-------|------|
| Range default | 1550–2650 CE (DE440) | 1550–2650 CE (DE440) — invariato |
| Range esteso | Moshier -3000/+3000 CE | **DE441 -13200/+17191 CE** (env var) |
| Config effemeridi | Solo `set_ephemeris_file()` | `set_ephemeris_file()` + `LIBEPHEMERIS_EPHEMERIS` env var |
| `SEFLG_MOSEPH` | Attiva Moshier | Accettato e ignorato silenziosamente |
| Codice rimosso | — | ~8 file moshier/ + ~300 righe planets.py + ~50 file test |
| Kepleriani | Funzionanti | **Invariati** |

## Dipendenza critica: `spk.py` → `moshier/`

`spk.py` (righe 855-856) importa dal pacchetto `moshier/` funzioni usate anche nel path JPL:

- `moshier.precession.nutation_angles`
- `moshier.precession.precess_from_j2000`
- `moshier.utils.apply_aberration_to_position`

Queste vanno estratte in un nuovo modulo `libephemeris/astrometry.py` **prima** di eliminare il pacchetto `moshier/`.

---

## Fasi di implementazione

### Fase 1 — `state.py`: Env var `LIBEPHEMERIS_EPHEMERIS` + rimozione fallback DE421

**File**: `libephemeris/state.py`

- Aggiungere `_EPHEMERIS_ENV_VAR = "LIBEPHEMERIS_EPHEMERIS"`
- Modificare `get_planets()` per leggere l'env var con priorità: `set_ephemeris_file()` > env var > default `de440.bsp`
- Rimuovere il blocco fallback DE421 (righe 153-175)

### Fase 2 — Nuovo `astrometry.py` + aggiornamento `spk.py`

**File nuovi**: `libephemeris/astrometry.py`
**File modificati**: `libephemeris/spk.py`

Estrarre da `moshier/precession.py` e `moshier/utils.py`:
- `nutation_angles`, `precess_from_j2000` (da `precession.py`)
- `apply_aberration_to_position` (da `utils.py`)
- Costanti e helper necessari: `ARCSEC_TO_RAD`, `J2000`, `JD_PER_CENTURY`, `jd_to_julian_centuries`, `normalize_angle`, `cartesian_to_spherical`, ecc.

Aggiornare import in `spk.py` (righe 855-856):
```python
# Prima:
from .moshier.precession import nutation_angles, precess_from_j2000
from .moshier.utils import apply_aberration_to_position

# Dopo:
from .astrometry import nutation_angles, precess_from_j2000, apply_aberration_to_position
```

### Fase 3 — Rimuovere routing Moshier dal core

**File modificati**: `libephemeris/planets.py`, `libephemeris/houses.py`, `libephemeris/time_utils.py`

- **`planets.py`**:
  - Eliminare `_calc_body_moshier()` (~300 righe, 733-1046)
  - Eliminare `_MOSHIER_SUPPORTED_BODIES` (righe 167-189)
  - Rimuovere gate `if iflag & SEFLG_MOSEPH:` in `swe_calc_ut()` (righe 1115-1118) e `swe_calc()` (righe 1166-1169)
  - Rimuovere import di `SEFLG_MOSEPH` (riga 102)
  - Mascherare il flag prima del calcolo (strip bit 4)
  - Aggiornare docstrings

- **`houses.py`**:
  - Rimuovere `SEFLG_MOSEPH` da import (riga 81)
  - Rimuovere dalla logica di propagazione flag (righe 762-765, 1548-1549)

- **`time_utils.py`**:
  - Rimuovere `SEFLG_MOSEPH` da import (riga 15)
  - Rimuovere riferimenti da `swe_deltat_ex()` (righe 213, 239-240, 266, 269)

### Fase 4 — Pulizia `exceptions.py`

**File modificato**: `libephemeris/exceptions.py`

- Rimuovere costanti: `MOSHIER_JD_START`, `MOSHIER_JD_END`, `MOSHIER_START_YEAR`, `MOSHIER_END_YEAR` (righe 966-975)
- Rimuovere funzione: `validate_jd_range_moshier()` (righe 978-1047)

### Fase 5 — Eliminare pacchetto `libephemeris/moshier/`

**File eliminati** (8 file):
- `libephemeris/moshier/__init__.py`
- `libephemeris/moshier/vsop87.py`
- `libephemeris/moshier/elp82b.py`
- `libephemeris/moshier/elp82b_data.py`
- `libephemeris/moshier/pluto.py`
- `libephemeris/moshier/pluto_data.py`
- `libephemeris/moshier/precession.py`
- `libephemeris/moshier/utils.py`

### Fase 6 — Eliminare test Moshier

**File eliminati in `tests/`** (~10 file):
- `tests/test_moshier_package.py`
- `tests/test_calc_body_moshier.py`
- `tests/test_moshier_compare_planets.py`
- `tests/test_moshier_routing.py`
- `tests/test_moshier_vsop87_precision.py`
- `tests/test_moshier_precession.py`
- `tests/test_moshier_pluto.py`
- `tests/test_moshier_utils.py`
- `tests/test_seflg_moseph_documentation.py`
- `tests/test_houses/test_houses_moshier.py`
- `tests/test_lunar/test_moshier_apogee.py`

**File eliminati in `compare_scripts/tests/test_moshier/`** (~39 file):
- Tutta la directory `compare_scripts/tests/test_moshier/`

### Fase 7 — Aggiornare `__init__.py`

**File modificato**: `libephemeris/__init__.py`

- Rimuovere export di `validate_jd_range_moshier` (se presente)
- **Mantenere** `SEFLG_MOSEPH` e `FLG_MOSEPH` nella lista export e in `__all__`

### Fase 8 — Aggiornare `README.md`

**File modificato**: `README.md`

- Rimuovere bullet Moshier analytical ephemeris (riga 114)
- Aggiornare struttura progetto — rimuovere `moshier/` (riga 587)
- Aggiungere `LIBEPHEMERIS_EPHEMERIS` nella sezione ephemeris file selection
- Aggiornare sezione configurazione con tabella env var completa

### Fase 9 — Aggiornare `AGENTS.md`

**File modificato**: `AGENTS.md`

- Rimuovere sezione su Moshier (righe 182-186)
- Aggiornare struttura progetto — rimuovere `moshier/`
- Aggiungere nota su env var `LIBEPHEMERIS_EPHEMERIS` e supporto DE441

### Fase 10 — Aggiornare `CHANGELOG.md`

**File modificato**: `CHANGELOG.md`

Nuova entry versione con:
- BREAKING: `SEFLG_MOSEPH` accettato ma ignorato (usa sempre JPL)
- Aggiunta env var `LIBEPHEMERIS_EPHEMERIS` per selezione effemeridi (es. `de441.bsp`)
- Rimozione fallback DE421
- Rimozione pacchetto `moshier/`

### Fase 11 — Aggiornare `docs/`

| File | Modifica |
|------|----------|
| `docs/PRECISION_TUNING.md` | Rimuovere Moshier da tabella (riga 27), aggiornare sezione ephemeris file selection, aggiungere env var |
| `docs/PRECISION.md` | Aggiornare tabella range (riga 258), rimuovere riferimenti Moshier |
| `docs/migration-guide.md` | Aggiornare sezione ephemeris configuration, aggiungere env var |
| `docs/api_reference.rst` | Rimuovere `SEFLG_MOSEPH` dal parametro `ephe_flag` (riga 263) |
| `docs/PLANET_CENTERS_SPK.md` | Rimuovere riferimenti a DE421 |

### Fase 12 — Verifica

- `poe test` — tutti i test devono passare
- `poe typecheck` — nessun errore mypy
- Verifica manuale che `SEFLG_MOSEPH` sia accettato e ignorato
