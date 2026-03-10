# Piano Integrato: API Alignment, Cusp Speeds, e Miglioramenti di Precisione

## Panoramica

Tre aree di intervento prioritarie su libephemeris, ordinate per urgenza:

1. **Allineamento signature API** — 6 funzioni con signature incompatibili con pyswisseph
2. **Velocita cuspidi case** — 3 sistemi (Whole Sign, Porphyry, Koch) con derivate numeriche errate
3. **Miglioramenti di precisione** — Eliminare xfail per minor bodies, pianeti ipotetici, nod_aps

---

## FASE 1: Allineamento Signature API

Obiettivo: compatibilita 1:1 con pyswisseph. Breaking changes accettate.

### 1.1 `swe_deltat_ex` (Semplice)

| | Attuale | Target |
|---|---------|--------|
| **File** | `time_utils.py:201` | |
| **Signature** | `(tjd, ephe_flag) -> (float, str)` | `(tjd, ephe_flag) -> float` |
| **Callers** | 11 in tests/, 16 in compare_scripts/ | |

**Modifiche:**
- Cambiare il return da `(delta_t, serr)` a solo `delta_t`
- Aggiornare tutti i 27 callers che fanno unpacking `dt, serr = swe_deltat_ex(...)`
- Il messaggio `serr` non e usato da nessun caller — non serve preservarlo

**Verifica:** `pytest tests/test_utc_leap_seconds.py tests/test_time/ -v -k deltat`

### 1.2 `swe_get_ayanamsa_ex_ut` e `swe_get_ayanamsa_ex` (Complesso)

| | Attuale | Target |
|---|---------|--------|
| **File** | `planets.py:3351` e `planets.py:3312` | |
| **Signature** | `(tjd, sid_mode, flags) -> (ayan, eps, nut)` | `(tjd, flags) -> (flags, ayan)` |
| **Callers** | ~9 in tests/, ~20 in compare_scripts/ | |

**Modifiche:**
- Rimuovere parametro `sid_mode` — deve venire dallo stato globale via `swe_set_sid_mode()`
- Cambiare return da `(ayan, eps, nut)` a `(retflag, ayan)`
- Verificare che `swe_set_sid_mode()` sia gia implementato e funzionante
- Aggiornare tutti i callers

**Verifica:** `pytest tests/test_sidereal/ -v` e `pytest compare_scripts/tests/test_sidereal/ -v`

### 1.3 `swe_get_orbital_elements` (Medio)

| | Attuale | Target |
|---|---------|--------|
| **File** | `planets.py:3832` | |
| **Signature** | `(tjd_et, ipl, iflag) -> (tuple_17, int)` | `(tjd_et, ipl, iflag) -> tuple_50` |
| **Callers** | ~21 in tests/, ~4 in compare_scripts/ | |

**Modifiche:**
- Rimuovere il wrapping `(elements, retflag)` — ritornare una flat tuple
- Espandere da 17 a 50 elementi (pad con 0.0 per gli indici non calcolati)
- Aggiornare i 25 callers che fanno `result[0][0]` -> `result[0]`

**Verifica:** `pytest tests/test_planets/test_orbital_elements.py -v`

### 1.4 `swe_heliacal_ut` (Medio)

| | Attuale | Target |
|---|---------|--------|
| **File** | `heliacal.py:1325` | |
| **Signature** | `(jd, geopos, datm, dobs, name, type, flag) -> (tuple_50, int)` | stesso -> `(jd1, jd2, jd3)` |
| **Callers** | ~25 in tests/, ~5 in compare_scripts/ | |

**Modifiche:**
- Cambiare return da `(dret_50, retflag)` a `(jd1, jd2, jd3)`
- Aggiornare callers che accedono `result[0][0]` -> `result[0]`
- Allineare anche il duplicato in `eclipse.py:8056`

**Verifica:** `pytest tests/test_heliacal.py tests/test_heliacal_stars.py -v`

### 1.5 `swe_orbit_max_min_true_distance` — Solo test stale

| | Attuale | Target |
|---|---------|--------|
| **File** | `tests/test_planets/test_orbit_max_min_true_distance.py` | |
| **Impl** | Gia corretta (3-tuple) | |
| **Test** | 17 call sites aspettano 2-tuple | |

**Modifiche:**
- Aggiornare tutti i 17 call sites nei test unitari
- `assert len(result) == 2` -> `assert len(result) == 3`
- Aggiornare unpacking: `min_dist, max_dist = ...` -> `max_dist, min_dist, true_dist = ...`

**Verifica:** `pytest tests/test_planets/test_orbit_max_min_true_distance.py -v`

### 1.6 `swe_houses_armc ascmc[3]` — Nessun fix necessario

L'analisi ha confermato che pyswisseph e libephemeris concordano: entrambi mettono il Vertex
all'indice 3. I test di comparazione confermano la corrispondenza entro 0.001 deg.

---

## FASE 2: Velocita Cuspidi Case

Obiettivo: ridurre l'errore delle derivate numeriche per Koch, Porphyry e Whole Sign.

Codice di differenziazione centrale:
- `houses.py:1608-1641` (`swe_houses_ex2`)
- `houses.py:1417-1452` (`swe_houses_armc_ex2`)

Override per sistema specifico: dentro il blocco `if flags & SEFLG_SPEED:`, prima del fallback generico.

### 2.1 Whole Sign (W) — Override analitico

**Problema:** `math.floor(asc / 30.0) * 30.0` e una step function. Derivata 0 ovunque tranne
ai punti di discontinuita. Central difference ritorna ~0 o spikes ~280-550 deg/day.

**Soluzione:**
- `_whole_sign_cusp_speeds()` -> `(0.0,) * 12` (le cuspidi Whole Sign sono fisse ai 0 dei segni)
- Le velocita ASCMC vanno comunque calcolate con central difference

### 2.2 Porphyry (O) — Formula algebrica

**Problema:** Cuspidi sono combo lineari di ASC/MC, ma central difference introduce errore O(h^2).

**Soluzione:** Velocita analitiche:
```
cusps[11] speed = (2/3)*mc_speed + (1/3)*asc_speed
cusps[12] speed = (1/3)*mc_speed + (2/3)*asc_speed
cusps[2]  speed = (2/3)*asc_speed + (1/3)*mc_speed
cusps[3]  speed = (1/3)*asc_speed + (2/3)*mc_speed
cusps[1]  speed = asc_speed
cusps[10] speed = mc_speed
cusps[4]  speed = mc_speed  (IC)
cusps[7]  speed = asc_speed (DESC)
cusps[5-6,8-9] = stesse velocita dei cuspidi opposti
```

### 2.3 Koch (K) — Timestep ridotto

**Problema:** 5 funzioni trigonometriche annidate amplificano errore O(h^2) a 1 minuto.

**Soluzione:** Ridurre `dt` a `1.0/86400.0` (1 secondo) solo per Koch. Il costo computazionale
e trascurabile (3 chiamate extra a `_houses_koch` che e pura trigonometria). Riduce errore ~3600x.

### 2.4 Aggiornamento test di comparazione

- Aggiungere K, O, W ai sistemi testati in `test_deep_validation_4.py:753` e `test_all_systems.py:691`
- Rimuovere commenti di esclusione (`test_all_systems.py:701-703`)
- Definire tolleranze appropriate per sistema

---

## FASE 3: Miglioramenti di Precisione

### 3.1 Minor Bodies — SPK auto-download

**Stato attuale:** SPK auto-download gia ON di default (`state.py:1406`). Strict precision
anche True di default.

**Azione:**
- Verificare che i big 6 siano in `SPK_BODY_NAME_MAP`
- Abilitare auto-download nel conftest dei compare_scripts
- Rimuovere xfail dai test che passano con SPK
- Marker `@pytest.mark.network` per test che richiedono rete

**Verifica:** `pytest compare_scripts/tests/test_compare_minor_bodies.py -v`

### 3.2 Pianeti Ipotetici — Allineamento coefficienti

**Stato attuale:** Tre sistemi di coefficienti paralleli in `hypothetical.py`. Valori differiscono.

**Azione:**
1. Generare posizioni di riferimento con pyswisseph
2. Identificare coefficienti causa dello scostamento
3. Aggiornare usando fonti indipendenti (Witte 1928, Neely 1988, Hamburg School)
4. Unificare i tre sistemi paralleli

**Vincolo:** Coefficienti NON da Swiss Ephemeris `fictdata.txt`. Solo fonti indipendenti.

### 3.3 nod_aps Geocentrico

**Stato attuale:** `_calc_nod_aps()` usa elementi eliocentrici. Differenza fino a ~250 per Mercurio/Venere.

**Azione:**
1. Implementare nodi osculanti geocentrici (posizione+velocita geocentrica -> piano orbitale -> intersezione eclittica)
2. Mantenere metodo eliocentrico come opzione
3. Default geocentrico

### 3.4 Fallback Kepleriani — Miglioramenti secondari (futuro)

- Perturbazioni short-period per Marte/Urano/Nettuno
- Teoria di Brouwer & Clemence completa (secondo ordine)
- Yarkovsky piu accurato

---

## Ordine di Esecuzione

```
Fase 1 — API Alignment
  1.1 swe_deltat_ex                    [30 min]
  1.5 test stale orbit_max_min         [20 min]
  1.2 swe_get_ayanamsa_ex_ut/ex       [1h]
  1.3 swe_get_orbital_elements         [45 min]
  1.4 swe_heliacal_ut                  [1h]

Fase 2 — Cusp Speeds
  2.1 Whole Sign analytical override   [30 min]
  2.2 Porphyry analytical override     [30 min]
  2.3 Koch timestep reduction          [1h]
  2.4 Update comparison tests          [30 min]

Fase 3 — Precision
  3.1 SPK auto-download in tests       [1h]
  3.2 Coefficienti ipotetici           [2 giorni]
  3.3 nod_aps geocentrico              [3 giorni]
```

---

## Rischi e Decisioni Aperte

| # | Rischio/Decisione | Raccomandazione |
|---|-------------------|-----------------|
| 1 | `swe_heliacal_ut` duplicato in `eclipse.py:8056` | Allineare entrambi |
| 2 | Koch: derivata analitica vs timestep ridotto | Timestep ridotto — piu semplice, robusto |
| 3 | Whole Sign: pyswisseph ritorna 0 o velocita ASC? | Verificare empiricamente |
| 4 | Coefficienti ipotetici: fonte indipendente? | Witte 1928, Neely 1988 |
| 5 | nod_aps: cambiare default o aggiungere flag? | Default geocentrico + flag eliocentrico |
| 6 | Test SPK richiedono rete — CI offline? | Marker `@pytest.mark.network` |
