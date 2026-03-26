# LEB2: Piano di Implementazione

**Status:** Implementato in v1.0.0a2. Fase 7 (packaging PyPI) pianificata per v1.0.0a3.
**Data:** 2026-03-26
**Basato su:** `proposals/leb2-compressed-format.md` (proposta originale, superata dal benchmark)

---

## 1. Contesto Generale

### Cos'e' LibEphemeris
LibEphemeris e' una reimplementazione Python 1:1 API-compatibile della Swiss Ephemeris.
Usa Skyfield + JPL DE440/DE441 come backend, con un fast-path opzionale basato su file
`.leb` contenenti coefficienti Chebyshev precomputed (~14x speedup vs Skyfield).

### Il Problema
I file `.leb` non entrano su PyPI (limite 100 MB per file):

| Tier | Copertura | File corrente | Post-ottimizzazione params |
|------|-----------|---------------|---------------------------|
| Base | 1850-2150 | 107 MB | ~97 MB (stima) |
| Medium | 1550-2650 | 377 MB | ~340 MB (stima) |
| Extended | -5000/+5000 | 2.8 GB | ~2.5 GB (stima) |

Nemmeno il base tier entra in PyPI. L'utente deve scaricare DE440 (~114 MB) o file
`.leb` separatamente prima di poter usare la libreria.

### L'Obiettivo
Spedire un file `.leb` compresso direttamente nel wheel PyPI, cosi':
```
pip install libephemeris
python -c "import libephemeris as swe; print(swe.calc_ut(2451545.0, 0))"
# Funziona immediatamente, senza download extra
```

---

## 2. Lavoro Gia' Fatto

### 2.1 Ottimizzazione parametri Chebyshev (completato)
**Commit:** `0fb2497` su `feat/leb-param-optimization`
**File modificato:** `libephemeris/leb_format.py` (BODY_PARAMS)

| Corpo | Vecchio (int/deg) | Nuovo (int/deg) | Risparmio |
|-------|-------------------|-----------------|-----------|
| Uraniani (40-48) | 32d / deg 13 | 256d / deg 7 | 9.5 MB |
| Pluto | 32d / deg 13 | 64d / deg 11 | 0.6 MB |
| **Totale** | | | **~10.1 MB** |

I parametri di tutti gli altri corpi sono stati testati e documentati in
`proposals/leb-optimization-findings.md`: Mercury, Venus, Mars, Jupiter,
Saturn, Uranus, Neptune, asteroidi, nodi lunari — nessuno puo' essere
ulteriormente ottimizzato senza superare la soglia di 0.001" di errore.

### 2.2 Benchmark compressione (completato, non committato)
Testato nell'attuale sessione. Risultati chiave:

**Compressione lossless** (byte-shuffle + zstd-19): solo **1.24x** sul file
corrente. La stima di 3:1 nella proposta originale era errata: i coefficienti
Chebyshev in coordinate ICRS baricentricitiche hanno mantisse troppo "casuali"
per comprimere bene anche dopo byte-shuffle.

**Compressione lossy** (mantissa truncation + coeff-major reorder + shuffle + zstd):
risultati eccellenti. La chiave e' che i coefficienti di ordine alto (c6-c13) hanno
magnitudini infinitesimali e possono essere arrotondati aggressivamente senza perdere
precisione:

| Corpo | Coefficienti "utili" | Ratio lossy |
|-------|---------------------|-------------|
| Sun | c0-c3 (21-3 bit) | 13.5x |
| Moon | c0-c5 (28-4 bit) | 8.7x |
| Earth | c0-c4 (28-3 bit) | 12.0x |
| MeanNode | c0-c1 (37-26 bit) | 134x |
| Chiron | c0-c12 (32-13 bit) | 4.2x |
| OscuApog | c0-c15 (misti) | 3.3x |

### 2.3 Analisi composizione per gruppi (completato)
Benchmark sui 4 gruppi di body:

| Gruppo | Body | Raw MB | Lossy MB | Ratio |
|--------|------|--------|----------|-------|
| core | 14 (pianeti+nodi+Lilith) | 44.5 | 4.6 | 9.6x |
| apogee | 3 (OscuApog, IntpApog/Perig) | 31.6 | 8.4 | 3.8x |
| asteroids | 5 (Chiron, Ceres, Pallas, Juno, Vesta) | 23.0 | 5.7 | 4.0x |
| uranians | 9 (Cupido-Transpluto) | 0.7 | 0.1 | 8.6x |

Scenari packaging con overhead:

| Scenario | Size | vs LEB1 |
|----------|------|---------|
| **core_base (PyPI)** | **6.6 MB** | **15.5x** |
| core+asteroids_base | 12.3 MB | 8.3x |
| full_base (31 body) | 20.8 MB | 4.9x |

---

## 3. Decisioni Prese

1. **Pacchetto PyPI**: solo `core` (14 body, ~6.6 MB base) nel wheel
2. **Precisione**: 0.001" (identica al target LEB1)
3. **Compressione**: lossy error-bounded (mantissa truncation + coeff-major reorder + byte-shuffle + zstd-19)
4. **Architettura modulare**: file `.leb` separati per gruppo, stessi 3 tier
5. **Compatibilita'**: LEB1 resta supportato, `open_leb()` auto-detecta via magic bytes

### Struttura file risultante

```
# Nella distribuzione PyPI (wheel):
libephemeris/data/core_base.leb          ~6.6 MB

# Scaricabili separatamente (download command):
~/.libephemeris/leb/core_base.leb        ~6.6 MB  (oppure gia' nel wheel)
~/.libephemeris/leb/core_medium.leb      ~23 MB   (stima)
~/.libephemeris/leb/core_extended.leb    ~180 MB  (stima)
~/.libephemeris/leb/asteroids_base.leb   ~5.7 MB
~/.libephemeris/leb/asteroids_medium.leb ~20 MB   (stima)
~/.libephemeris/leb/asteroids_extended.leb ~80 MB (stima)
~/.libephemeris/leb/apogee_base.leb      ~8.4 MB
~/.libephemeris/leb/uranians_base.leb    ~0.1 MB
... (ogni gruppo x ogni tier)
```

Quando l'utente chiede un body non presente nel file caricato, il sistema:
1. Cerca file `.leb` aggiuntivi nella stessa directory
2. Se non trovato, fallback a Skyfield (come oggi)

---

## 4. Tecnica di Compressione (dettaglio)

### 4.1 Mantissa Truncation (lossy, error-bounded)

Per ogni ordine di coefficiente `k` di ogni body, calcoliamo il numero minimo di
bit di mantissa necessari affinche' l'errore di arrotondamento resti sotto il
target di precisione (0.001" = 5e-9 AU):

```
bits_needed(k) = ceil(-log2(target / max|c_k|))
```

Se `max|c_k| < target`, il coefficiente viene azzerato (0 bit necessari).

I bit in eccesso nella mantissa IEEE 754 vengono azzerati tramite bitmask:
```python
mask = 0xFFFFFFFFFFFFFFFF << (52 - keep_bits)
uint64_repr &= mask
```

Questo e' lossy ma l'errore e' rigorosamente bounded: ogni coefficiente perde
al massimo `max|c_k| * 2^(-bits)` di precisione, e la somma degli errori su
tutti i coefficienti resta sotto il target.

### 4.2 Coefficient-Major Reorder

I coefficienti vengono riordinati da segment-major:
```
seg0: [c0_x, c1_x, ..., c13_x, c0_y, ..., c13_y, c0_z, ..., c13_z]
seg1: [c0_x, c1_x, ..., c13_x, ...]
```
A coefficient-major:
```
[c0_x_seg0, c0_x_seg1, ...., c0_y_seg0, c0_y_seg1, ..., c1_x_seg0, ...]
```

Questo raggruppa valori simili (tutti i c0 insieme, tutti i c1, ecc.),
creando sequenze con esponenti ripetitivi che comprimono meglio.

### 4.3 Byte Shuffle + zstd

Dopo truncation e reorder, i byte vengono trasposti (byte lane 0 di tutti
i float insieme, poi lane 1, ecc.) e compressi con zstd level 19. I byte
azzerati dalla truncation comprimono a quasi nulla.

### 4.4 Decompressione (runtime)

Reversibile: zstd decompress -> unshuffle -> inverse reorder.
I float risultanti hanno mantisse troncate ma restano float64 validi.
Il reader poi usa lo stesso algoritmo Clenshaw di LEB1.

**Nota**: NON serve "de-truncare" — i float troncati sono gia' validi per
la valutazione. L'unica differenza e' che gli ultimi bit di mantissa sono zero.

---

## 5. Piano di Implementazione

### Fase 1: Primitives (compressione/decompressione)
**Nuovo file:** `libephemeris/leb_compression.py`

Funzioni:
- `compute_mantissa_bits(coeffs_3d, target_precision) -> list[int]`
  Calcola i bit necessari per ordine di coefficiente.
- `truncate_mantissa(data, keep_bits) -> ndarray`
  Azzera bit di mantissa in eccesso.
- `reorder_coeff_major(data, segment_count, degree, components) -> bytes`
  Riordina da segment-major a coefficient-major.
- `reorder_segment_major(data, segment_count, degree, components) -> bytes`
  Inverso del precedente (per decompressione).
- `shuffle_bytes(data, element_size=8) -> bytes`
- `unshuffle_bytes(data, element_size=8) -> bytes`
- `compress_body(raw_coeffs, segment_count, degree, components, target_prec) -> tuple[bytes, list[int]]`
  Pipeline completa: truncation -> reorder -> shuffle -> zstd. Ritorna (blob, bits_per_order).
- `decompress_body(compressed, uncompressed_size, segment_count, degree, components) -> bytes`
  Pipeline inversa: zstd -> unshuffle -> inverse reorder.

**Dipendenza:** `zstandard>=0.22.0` in pyproject.toml

**Test:** `tests/test_leb/test_leb_compression.py`
- Round-trip per ogni primitiva
- Round-trip compress/decompress su dati realistici
- Verifica che max error e' sotto il target
- Edge cases

### Fase 2: Formato LEB2 (costanti + serializzazione)
**File da editare:** `libephemeris/leb_format.py`

Aggiungere:
- `LEB2_MAGIC = b"LEB2"`, `LEB2_VERSION = 1`
- `SECTION_COMPRESSED_CHEBYSHEV = 6`
- `COMPRESSION_ZSTD_TRUNC_SHUFFLE = 1`
- Dataclass `CompressedBodyEntry` (68 bytes: 52 LEB1 + compressed_size + uncompressed_size)
- `COMPRESSED_BODY_ENTRY_FMT`, `COMPRESSED_BODY_ENTRY_SIZE`
- Helpers `write_compressed_body_entry()` / `read_compressed_body_entry()`
- Campo `bits_per_order` in CompressedBodyEntry (o in header aggiuntivo)

### Fase 3: LEB2 Reader
**Nuovo file:** `libephemeris/leb2_reader.py`

Classe `LEB2Reader` — stessa interfaccia di `LEBReader`:
- Constructor: open, mmap, parse header/sections/body index
- `_decompress_body(body_id)`: legge blob compresso, chiama `decompress_body()`,
  cache in `self._cache: dict[int, bytes]`
- `eval_body(body_id, jd)`: identico a LEBReader dopo decompressione
- Riusa `_clenshaw`, `_clenshaw_with_derivative` da `leb_reader.py`
- Nutation/delta-T/stars: letti da mmap (non compressi)

**Test:** `tests/test_leb/test_leb2_reader.py`
- Eval body: risultati entro 0.001" dal riferimento Skyfield
- Nutation/delta-T/stars identici a LEB1
- Context manager, error cases

### Fase 4: LEB2 Writer (generazione)
**File da editare:** `scripts/generate_leb.py`

Aggiungere flag `--format {leb1,leb2}` e `--group {core,asteroids,apogee,uranians,all}`:
- Stessa pipeline di fitting Chebyshev
- Per LEB2: dopo il fitting, applica `compress_body()` per-body
- Genera file con magic `LEB2`, sezione `SECTION_COMPRESSED_CHEBYSHEV`
- Supporto output per gruppo singolo

### Fase 5: Integrazione nel runtime
**File da editare:**

`libephemeris/leb_reader.py`:
- Protocol `LEBReaderLike` (typing.Protocol)
- Factory `open_leb(path) -> LEBReader | LEB2Reader` (auto-detect via magic)

`libephemeris/state.py` (linee ~334-336):
```python
# da: from .leb_reader import LEBReader; _LEB_READER = LEBReader(path)
# a:  from .leb_reader import open_leb; _LEB_READER = open_leb(path)
```

`libephemeris/context.py` (linee ~212-214):
```python
# da: from .leb_reader import LEBReader; self._leb_reader = LEBReader(...)
# a:  from .leb_reader import open_leb; self._leb_reader = open_leb(...)
```

`libephemeris/fast_calc.py` (linea ~56):
```python
# da: from .leb_reader import LEBReader
# a:  from .leb_reader import LEBReaderLike as LEBReader
```

### Fase 6: Multi-file discovery
**File da editare:** `libephemeris/state.py`

Estendere `get_leb_reader()` per supportare piu' file `.leb` nella stessa directory:
- Il core file e' sempre caricato (dal wheel o da download)
- File aggiuntivi (asteroids, apogee, uranians) vengono scoperti automaticamente
- Un `CompositeLEBReader` delega `eval_body(body_id, jd)` al reader corretto

### Fase 7: Packaging
**File da editare:** `pyproject.toml`

- `zstandard>=0.22.0` nelle dependencies
- Includere `core_base.leb` nel package data (o generarlo nel build step)
- Aggiornare download command per i file aggiuntivi

---

## 6. File Coinvolti (riepilogo)

| File | Azione | Fase |
|------|--------|------|
| `pyproject.toml` | Edit: add zstandard dep + package data | 1, 7 |
| `libephemeris/leb_compression.py` | **Nuovo** | 1 |
| `tests/test_leb/test_leb_compression.py` | **Nuovo** | 1 |
| `libephemeris/leb_format.py` | Edit: LEB2 constants + CompressedBodyEntry | 2 |
| `libephemeris/leb2_reader.py` | **Nuovo** | 3 |
| `tests/test_leb/test_leb2_reader.py` | **Nuovo** | 3 |
| `scripts/generate_leb.py` | Edit: --format leb2 + --group | 4 |
| `libephemeris/leb_reader.py` | Edit: Protocol + open_leb() factory | 5 |
| `libephemeris/state.py` | Edit: open_leb() + multi-file discovery | 5, 6 |
| `libephemeris/context.py` | Edit: open_leb() | 5 |
| `libephemeris/fast_calc.py` | Edit: TYPE_CHECKING import | 5 |

---

## 7. Verifica

1. **Unit test compressione**: round-trip, error bounds
2. **Test reader LEB2**: eval_body entro 0.001" dal riferimento
3. **Test integrazione**: `set_leb_file("core_base.leb")` -> `calc_ut()` corretto
4. **Test multi-file**: body nel core -> OK, body negli asteroidi -> caricamento automatico
5. **Regressione LEB1**: tutti i test LEB esistenti passano con file LEB1
6. **Size check**: `core_base.leb` < 10 MB

---

## 8. Rischi e Mitigazioni

| Rischio | Mitigazione |
|---------|-------------|
| Truncation introduce errore > 0.001" | Verifica end-to-end su campione denso di date. Il calcolo dei bit e' conservativo (usa max|c_k| globale). |
| zstandard come nuova dep | ~200 KB wheel, ben mantenuto, nessun rischio reale |
| Complessita' multi-file | CompositeLEBReader e' semplice (dict di reader per body_id). Fallback a Skyfield se body non trovato. |
| File LEB2 non rigenerabile senza Skyfield | Il writer riusa la stessa pipeline di fitting. L'unica differenza e' il passo di serializzazione. |
