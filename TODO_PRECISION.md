# TODO: Precisione Scientifica di LibEphemeris

> Obiettivo: rendere LibEphemeris **più** scientificamente precisa di SwissEphemeris in ogni aspetto.
>
> Nuova dipendenza obbligatoria: **pyerfa** (binding Python ufficiale della libreria ERFA dell'IAU).

---

## TODO 1 — CRITICO: Aberrazione mancante nel path ayanamsa star-based

**Impatto: fino a ~20.5″ di errore**

**File:** `libephemeris/planets.py`, funzione `_get_star_position_ecliptic()` (linee 2322–2535)

### Problema

Questa funzione calcola la posizione eclittica delle stelle di riferimento (Spica, Revati, ecc.)
per le ayanamsa star-based (True Citra, True Revati, True Pushya, ecc.), ma **non applica
l'aberrazione annuale**. L'aberrazione annuale sposta le posizioni stellari fino a ~20.5″
(la costante di aberrazione). Per confronto, `fixed_stars.py` applica correttamente
l'aberrazione usando `astrometric.apparent()` di Skyfield.

### Fix

**Opzione A:** Usare `erfa.ab()` (aberrazione annuale) dopo la trasformazione di precessione.

**Opzione B (raccomandata):** Riscrivere la funzione per usare il pipeline Skyfield:

```python
earth.at(t).observe(star).apparent().frame_latlon(ecliptic_frame)
```

Questo include automaticamente aberrazione, precessione, e nutazione.
La funzione `_calc_star_based_ayanamsha()` (linee 2923–3053) già implementa questo
approccio ma non viene usata — unificare i due path.

---

## TODO 2 — ALTO: Unificare la nutazione su IAU 2006/2000A (via pyerfa)

**Impatto: ~1 mas di consistenza, sub-milliarcsecond precision ovunque**

**File:** `libephemeris/cache.py` (linee 34–62)

### Problema

L'hot path della nutazione (usato per obliquità, case, ayanamsa true) usa
`iau2000b_radians` di Skyfield (77 termini, ~1 mas), mentre Skyfield internamente
(per `t.gast`, `ecliptic_frame`) usa `iau2000a_radians` (1365 termini, ~0.1 mas).
Questo crea un'inconsistenza di ~1 mas tra GAST e obliquità nel calcolo delle case.

### Fix

In `get_cached_nutation()` (linea 61), sostituire:

```python
dpsi_rad, deps_rad = iau2000b_radians(t)
```

con:

```python
import erfa
dpsi, deps = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
# dpsi e deps sono già in radianti
dpsi_rad, deps_rad = dpsi, deps
```

Questo usa il modello IAU 2006/2000A di pyerfa (il più preciso disponibile, ~0.01–0.05 mas).

### Propagazione automatica

- `get_cached_obliquity()` (cache.py:79)
- `get_true_obliquity()` (cache.py:119) → case
- `_get_true_ayanamsa()` (planets.py:2918)
- `_calc_ayanamsa()` star-based modes (planets.py:2695)
- `azalt()` / `azalt_rev()` (utils.py)

---

## TODO 3 — ALTO: Unificare l'obliquità su IAU 2006

**Impatto: 0.042″ costante + crescente con T**

**File:** `libephemeris/cache.py` (linee 99–104)

### Problema

Tre formule diverse per l'obliquità media:

| Posizione | Formula | Costante |
|-----------|---------|----------|
| `cache.py:103` (case, ayanamsa) | Laskar 1986 | 84381.448″ |
| `planets.py:1180` (SE_ECL_NUT) | IAU 2006 | 84381.406″ |
| `astrometry.py:47`, `lunar.py:1570` | IAU 2006 | 84381.406″ |

Il path critico (case, ayanamsa) usa Laskar 1986, con un offset costante di **0.042″**
dall'IAU 2006 e termini mancanti (T⁴, T⁵).

### Fix

In `get_cached_obliquity()`, sostituire:

```python
eps0_arcsec = 84381.448 - 46.8150 * T - 0.00059 * T**2 + 0.001813 * T**3
```

con una delle due opzioni:

**Opzione A — formula diretta IAU 2006:**

```python
eps0_arcsec = (84381.406
               - 46.836769 * T
               - 0.0001831 * T**2
               + 0.00200340 * T**3
               - 0.000000576 * T**4
               - 0.0000000434 * T**5)
```

**Opzione B — pyerfa (raccomandata):**

```python
eps0_rad = erfa.obl06(2451545.0, jd_tt - 2451545.0)
eps0 = math.degrees(eps0_rad)
```

Applicare lo stesso fix in `_calc_ayanamsa()` (linea 2688) dove l'obliquità è calcolata
per i modi siderali star-based.

---

## TODO 4 — ALTO: Correggere i coefficienti di precessione e aggiungere termini superiori

**Impatto: ~0.08 mas/cy immediato + crescente per date distanti da J2000**

**File:** `libephemeris/planets.py`, `_calc_ayanamsa()` (linee 2600–2603, 2884–2892)

### Problema 4a — Coefficienti leggermente errati

```python
# Attuale (linea 2600–2601):
PREC_RATE = 5028.796273       # off by +0.000078″/cy dal valore IAU 2006
PREC_RATE_QUAD = 1.105608     # off by +0.000173″/cy²

# Corretto (IAU 2006, Capitaine et al. 2003, A&A 412):
PREC_RATE = 5028.796195
PREC_RATE_QUAD = 1.1054348
```

### Problema 4b — Termini mancanti T³, T⁴, T⁵

```python
# Attuale (linea 2888):
ayanamsa = aya_j2000 + (precession * T + PREC_RATE_QUAD * T * T) / 3600.0

# Completo IAU 2006:
ayanamsa = aya_j2000 + (5028.796195 * T
                        + 1.1054348 * T**2
                        + 0.00007964 * T**3
                        - 0.000023857 * T**4
                        - 0.0000000383 * T**5) / 3600.0
```

### Problema 4c — `SE_SIDM_USER` troppo approssimativo

Linea 2867: usa `5027.8″/cy` senza termine quadratico (~1″/secolo di errore).
Aggiungere almeno il termine quadratico e aggiornare il rate al valore IAU 2006.

---

## TODO 5 — ALTO: Sostituire Lieske 1977 con IAU 2006 nella precessione stellare

**Impatto: ~0.1″/secolo + frame bias ~23 mas mancante**

**File:** `libephemeris/planets.py`, `_get_star_position_ecliptic()` (linee 2464–2466)

### Problema

I coefficienti zeta/z/theta sono di Lieske 1977 (IAU 1976), **non** di IAU 2006 come
dichiarato nel commento. Inoltre manca il frame bias GCRS→J2000 (~23 mas).

### Fix

**Opzione A — coefficienti IAU 2006 con frame bias incluso:**

```python
# IAU 2006 precession angles (Capitaine et al. 2003, A&A 412, Table 1)
# I termini costanti (±2.650545″) includono il frame bias GCRS→J2000
zeta  = (2.650545 + 2306.083227 * T + 1.0967790 * T**2
         + 0.01860606 * T**3 - 0.000013 * T**4
         - 0.0000005 * T**5) / 3600.0
z     = (-2.650545 + 2306.077181 * T + 1.0927348 * T**2
         + 0.01826837 * T**3 - 0.000028 * T**4
         - 0.0000003 * T**5) / 3600.0
theta = (2004.191903 * T - 0.4294934 * T**2
         - 0.04182264 * T**3 - 0.000007089 * T**4
         - 0.0000001274 * T**5) / 3600.0
```

**Opzione B (raccomandata) — pyerfa:**

```python
# Matrice di precessione IAU 2006 completa (include frame bias)
rbp = erfa.pmat06(2451545.0, (jd_tt - 2451545.0))
# Applicare: pos_precessed = rbp @ pos_j2000
```

### Cleanup associato

Rimuovere il dead code delle formule scalari A/B/C (linee 2471–2487) che calcolano la
precessione ma vengono sovrascritte dalla matrice di rotazione successiva.

---

## TODO 6 — MEDIO: Frame bias GCRS→J2000

**Impatto: ~23 mas sistematico**

**Risolto da TODO 5** se si usano i coefficienti IAU 2006 con i termini costanti
(±2.650545″ in zeta/z), oppure se si usa `erfa.pmat06()`.

Se si preferisce un fix separato, applicare la matrice di frame bias `erfa.bi00()`
prima della precessione in `_get_star_position_ecliptic()`.

---

## TODO 7 — MEDIO: Velocità analitiche invece di finite-difference

**Impatto: Luna ~0.001 deg/giorno; pianeti minore ma sistematico. Performance ~3×.**

**File:** `libephemeris/planets.py` (linee 1934–1991)

### Problema

Tutte le velocità planetarie sono calcolate con finite-difference numerico
(central difference O(h²)). Skyfield/jplephem forniscono velocità analitiche ICRS dalla
differenziazione dei polinomi Chebyshev del DE440, che sono già usate per le correzioni
COB ma **non** per le velocità eclittiche angolari.

### Fix per pianeti standard (Sun–Pluto)

Ottenere posizione e velocità eclittiche da Skyfield e calcolare la velocità angolare
con le derivate esatte:

```python
r_ecl, v_ecl = pos.frame_xyz_and_velocity(ecliptic_frame)
x, y, z = r_ecl.au
vx, vy, vz = v_ecl.au_per_d

xy_sq = x * x + y * y
r_sq = xy_sq + z * z
xy = math.sqrt(xy_sq)
r = math.sqrt(r_sq)

speed_lon = math.degrees((x * vy - y * vx) / xy_sq)         # dλ/dt
speed_lat = math.degrees((z * (x*vx + y*vy) / xy - xy * vz) / r_sq)  # dβ/dt
speed_dist = (x * vx + y * vy + z * vz) / r                  # dr/dt
```

Questo approccio è già implementato per:
- SPK Type 21 (`spk.py:1006–1023`)
- Stelle fisse (`fixed_stars.py:3491–3530`)

Va esteso a tutti i pianeti standard — eliminare le 2 chiamate ricorsive `_calc_body()`
(performance ~3×).

### Sub-fix 7b — Forward-difference → central-difference

Sostituire forward-difference con central-difference per:

| Posizione | File:linea | dt attuale |
|-----------|------------|------------|
| Corpi uraniani | `hypothetical.py:2158` | 1/86400 (fwd) |
| Transpluto/Vulcan/Keplerian | `hypothetical.py:2473, 2710, 3345, 3493` | 1.0 giorno (fwd) |
| SPK Type 2/3 fallback | `spk.py:1248–1276` | 1/86400 (fwd) |
| Correzione rate ayanamsha | `planets.py:1987` | 1/86400 (fwd) |
| Lune planetarie | `planetary_moons.py:515` | 1/86400 (fwd) |

---

## TODO 8 — MEDIO: True Node — applicare correzioni perturbative

**Impatto: da ~8.9 arcsec a <1 arcsec vs SWE**

**File:** `libephemeris/lunar.py`, `calc_true_lunar_node()` (linee 1839–1870)

### Problema

Il calcolo usa solo `h = r × v` (meccanica orbitale pura). La funzione
`_calc_elp2000_node_perturbations()` (linee 334–1247, **900+ righe, 170+ termini**)
esiste nel file ma **non viene mai chiamata** dal calcolo del True Node.

### Fix proposto — approccio ibrido

1. Calcolare il nodo geometrico con `h = r × v` (come ora)
2. Calibrare un offset sistematico confrontando con SWE su un campione di date
3. Applicare la funzione `_calc_elp2000_node_perturbations()` come correzione residua
   alle perturbazioni a breve periodo non catturate dal vettore momento angolare
   istantaneo

Approccio alternativo: confrontare direttamente il nodo geometrico con i risultati SWE
su un set di date campione su tutto il range DE440 e derivare una funzione di correzione
empirica.

---

## TODO 9 — MEDIO: True Lilith — includere perturbazioni solari

**Impatto: da ~54 arcsec a <10 arcsec vs SWE**

**File:** `libephemeris/lunar.py`, `calc_true_lilith()` (linee 2135–2188)

### Problema

Il vettore eccentricità `e = (v × h)/μ − r/|r|` usa un modello a due corpi
(Terra-Luna), ignorando la perturbazione solare che causa oscillazioni fino a ~30°
nel vettore eccentricità. SWE stessa considera l'apogeo osculante "somewhat artificial"
ma produce comunque risultati più stabili.

### Fix proposto

1. Applicare `_calc_elp2000_apogee_perturbations()` (linee 1250–1421, ~50 termini) al
   risultato del vettore eccentricità come correzione perturbativa
2. Oppure includere la forza di marea solare nel calcolo del vettore eccentricità,
   trasformando il problema da due corpi a tre corpi ristretto

---

## TODO 10 — BASSO: Cleanup e consistenza

### 10a — Dead code in `_get_star_position_ecliptic()`

Rimuovere il calcolo scalare A/B/C di precessione (linee 2471–2487) — mai usato,
la matrice di rotazione (linee 2489–2510) lo sovrascrive.

### 10b — Funzione morta `_calc_star_based_ayanamsha()`

Rimuovere o integrare `_calc_star_based_ayanamsha()` (linee 2923–3053).
Se TODO 1 Opzione B viene implementato, questa funzione diventa la base del nuovo
path e il vecchio `_get_star_position_ecliptic()` può essere eliminato.

### 10c — Commento stale in `planets.py:34`

```python
# Attuale (errato):
# - Ecliptic frame uses J2000.0 for performance (true date would add ~0.01" precision but 2x slower)

# Corretto:
# - Ecliptic frame uses true ecliptic of date (Skyfield's ecliptic_frame with IAU 2006 precession + IAU 2000A nutation)
```

Il codice a linea 999 usa effettivamente `ecliptic_frame` che è la true ecliptic of date.

### 10d — Obliquità in `_calc_ayanamsa()`

Unificare anche l'obliquità a linea 2688 con IAU 2006 (contestualmente a TODO 3).

---

## Dipendenze

### Nuova dipendenza obbligatoria

**pyerfa** — da aggiungere come `required` in `pyproject.toml`:

```toml
dependencies = [
    "pyerfa>=2.0",
    # ... existing deps
]
```

Funzioni pyerfa usate:

| Funzione | TODO | Scopo |
|----------|------|-------|
| `erfa.nut06a()` | 2 | Nutazione IAU 2006/2000A |
| `erfa.obl06()` | 3 | Obliquità media IAU 2006 |
| `erfa.pmat06()` | 5 | Matrice precessione IAU 2006 |
| `erfa.ab()` | 1 | Aberrazione annuale (se Opzione A) |
| `erfa.bi00()` | 6 | Frame bias GCRS→J2000 (se fix separato) |

**Peso:** ~2 MB, puro C (wrapper ERFA/SOFA dell'IAU), nessuna dipendenza transitiva
significativa. È già usato opzionalmente in `erfa_nutation.py`.

---

## Grafo delle dipendenze tra TODO

```
TODO 1 (aberrazione)     ─── indipendente, fix più grande
TODO 2 (nutazione)       ─── prerequisito per TODO 3
TODO 3 (obliquità)       ─── dopo TODO 2
TODO 4 (precessione)     ─── indipendente
TODO 5 (Lieske→IAU2006)  ─── indipendente, risolve anche TODO 6
TODO 6 (frame bias)      ─── risolto da TODO 5
TODO 7 (velocità)        ─── indipendente, il più complesso
TODO 8 (true node)       ─── indipendente
TODO 9 (true lilith)     ─── indipendente
TODO 10 (cleanup)        ─── dopo tutto il resto
```

## Stima dell'impatto complessivo

| Metrica | Prima | Dopo |
|---------|-------|------|
| Ayanamsa star-based | fino a 20.5″ errore | < 0.1″ |
| Nutazione | ~1 mas (IAU 2000B, 77 termini) | ~0.01–0.05 mas (IAU 2006/2000A) |
| Obliquità | 0.042″ offset + inconsistenza | Consistente IAU 2006 ovunque |
| Precessione ayanamsa | ~0.08 mas/cy errore sistematico | Esatta IAU 2006 fino a T⁵ |
| Precessione stellare | Lieske 1977 (~0.1″/cy) | IAU 2006 + frame bias |
| Velocità Luna | ~0.001 deg/giorno (finite-diff) | Analitica (sub-arcsec/giorno) |
| True Node vs SWE | ~8.9 arcsec | Target < 1 arcsec |
| True Lilith vs SWE | ~54 arcsec | Target < 10 arcsec |
