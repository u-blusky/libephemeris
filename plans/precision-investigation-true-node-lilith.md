# Investigazione Precisione: Nodo Lunare Vero e Lilith Vera

## Contesto

Dopo 8 sessioni di lavoro e 15 commit su `dev/precision-v3`, libephemeris raggiunge
precisione sub-arcsecondo per tutti i corpi principali (Sole-Plutone), case, ayanamsha,
stelle fisse. Restano due scarti significativi:

| Corpo | Scarto vs pyswisseph | Tolleranza test |
|-------|---------------------|-----------------|
| **SE_TRUE_NODE** (Nodo Vero) | fino a ~510" (~0.14°) | 0.15° |
| **SE_OSCU_APOG** (Lilith Vera) | fino a ~235" (~0.065°) | 0.1° |

Questi sono gli scarti più grandi per corpi astrologicamente importanti.

**La domanda chiave**: lo scarto significa che noi siamo meno precisi, o che siamo
più precisi e SE ha un modello inferiore?

### Come calcoliamo noi

Entrambi usano **meccanica orbitale pura dai vettori di stato JPL DE440**:

**Nodo Vero** (`lunar.py:1696-1815`):
1. Posizione `r` e velocità `v` della Luna nel frame eclittico vero del giorno (Skyfield `ecliptic_frame`)
2. Momento angolare: `h = r × v`
3. Longitudine nodo ascendente: `Ω = atan2(h_x, -h_y)`

**Lilith Vera** (`lunar.py:1943-2062`):
1. Stessi `r` e `v` nel frame eclittico
2. Vettore eccentricità: `e = (v × h)/μ - r/|r|` (punta al perigeo)
3. Apogeo = direzione opposta: `lon = atan2(-e_y, -e_x)`
4. μ = G(M_Terra + M_Luna) con costanti IAU 2015

### Come calcola SE

Swiss Ephemeris usa un modello semi-analitico interno per la Luna (non JPL DE440
direttamente), e da quello estrae gli elementi osculanti. Il modello è basato sulla
teoria lunare ELP/MPP02, con correzioni proprietarie.

### Perché non sappiamo chi ha ragione

Entrambi gli approcci sono legittimi ma diversi:
- **Noi**: dati JPL DE440 (numerici, massima precisione) + estrazione geometrica diretta
- **SE**: modello semi-analitico (serie trigonometriche) + estrazione degli elementi

Lo scarto potrebbe venire da:
1. Differenza nel modello lunare di base (JPL vs ELP)
2. Differenza nella definizione dell'eclittica di riferimento
3. Bug nel nostro algoritmo di estrazione
4. Bug nel loro algoritmo di estrazione
5. Differenza nel frame di riferimento (eclittica vera vs media, precessione)

---

## Metodologia: Triangolazione con 3 Fonti Indipendenti

### Le 3 fonti

| # | Fonte | Tipo | Accesso |
|---|-------|------|---------|
| 1 | **libephemeris** | JPL DE440 via Skyfield + geometria vettoriale | Locale |
| 2 | **pyswisseph** | Modello semi-analitico SE | Locale |
| 3 | **JPL Horizons** | Integrazione numerica JPL, fonte di verità | API REST |

JPL Horizons è la fonte di verità perché:
- Usa gli stessi dati JPL DE440/DE441
- Restituisce direttamente gli **elementi osculanti** (Ω, ω, e, a, i, M)
- È il servizio ufficiale NASA/JPL
- Non ha bisogno di nessun algoritmo di estrazione — gli elementi sono calcolati
  internamente dal loro integratore

### Accesso a JPL Horizons

API REST: `https://ssd.jpl.nasa.gov/api/horizons.api`

Parametri per elementi osculanti lunari:
```
format=json
COMMAND='301'           # Luna (NAIF ID)
CENTER='500@399'        # Geocentrico (centro Terra)
EPHEM_TYPE=ELEMENTS     # Elementi osculanti
REF_PLANE=ECLIPTIC      # Piano eclittico
REF_SYSTEM=J2000        # Frame J2000 (da convertire a eclittica del giorno per confronto)
TLIST=<JD1>,<JD2>,...   # Date specifiche in JD TT
```

**Campi restituiti** che ci interessano:
- `OM` = Longitude of Ascending Node (Ω) → confronto con Nodo Vero
- `W` = Argument of Perifocus (ω) → apogeo = Ω + ω + 180°
- `EC` = Eccentricity (e) → verifica consistenza
- `A` = Semi-major axis (a) → verifica distanza apogeo
- `IN` = Inclination (i) → verifica piano orbitale

**Attenzione al frame**: Horizons restituisce Ω nel frame **J2000 eclittico** (non
eclittica del giorno). Bisogna:
1. O chiedere a Horizons il frame eclittica del giorno (`REF_SYSTEM=ECLIPDATE`)
2. O convertire i nostri risultati da eclittica del giorno a J2000
3. O applicare la precessione al valore Horizons

L'opzione più pulita è usare `REF_SYSTEM=ECLIPDATE` se Horizons lo supporta,
altrimenti fare il confronto in J2000 calcolando i nostri valori con `SEFLG_J2000`.

---

## Piano di Esecuzione

### Passo 0: Preparazione date di test

Scegliere ~30 date distribuite su 100 anni (1950-2050), incluse:
- 10 date "normali" distribuite uniformemente
- 5 date vicine a eclissi (nodo vicino al Sole → massima importanza astrologica)
- 5 date di massima oscillazione del nodo vero (quando d(Ω)/dt cambia segno)
- 5 date di perigeo/apogeo lunare (massima eccentricità → Lilith più sensibile)
- 5 date estreme (1950, 1970, 2030, 2050) per testare stabilità nel tempo

### Passo 1: Query JPL Horizons

Script: `scripts/investigate_true_node_horizons.py`

```python
# Pseudocodice
import requests
import json

dates_jd = [...]  # 30 JD in TT

# Query per elementi osculanti
params = {
    "format": "json",
    "COMMAND": "'301'",
    "CENTER": "'500@399'",
    "EPHEM_TYPE": "ELEMENTS",
    "REF_PLANE": "ECLIPTIC",
    "REF_SYSTEM": "ECLIPDATE",  # eclittica del giorno
    "TLIST": ",".join(str(jd) for jd in dates_jd),
    "CSV_FORMAT": "YES",
    "CAL_FORMAT": "JD",
}

response = requests.get("https://ssd.jpl.nasa.gov/api/horizons.api", params=params)
data = response.json()

# Estrarre: Ω (nodo), ω (arg perigeo), e, a, i per ogni data
# Salvare in CSV per riproducibilità
```

Output atteso per ogni data:
```
JD_TT, Omega_Horizons, omega_Horizons, e_Horizons, a_Horizons, i_Horizons
```

### Passo 2: Calcolo con libephemeris

Per ogni data:
```python
import libephemeris as ephem
from libephemeris.lunar import calc_true_lunar_node, calc_true_lilith

jd_tt = ...
# Nodo Vero
node_lon, node_lat, node_dist = calc_true_lunar_node(jd_tt)

# Lilith Vera
lilith_lon, lilith_lat, lilith_dist = calc_true_lilith(jd_tt)
```

**Attenzione**: i nostri valori sono in **eclittica del giorno** (true ecliptic of date).
Se Horizons restituisce in J2000, bisogna convertire. Verificare con una data nota.

### Passo 3: Calcolo con pyswisseph

```python
import swisseph as swe
swe.set_ephe_path("swisseph/ephe/")

# Nodo Vero (SE_TRUE_NODE = 11)
node_se, _ = swe.calc_ut(jd_ut, 11, swe.FLG_SPEED)
# node_se[0] = longitudine eclittica del giorno

# Lilith Vera (SE_OSCU_APOG = 13)
lilith_se, _ = swe.calc_ut(jd_ut, 13, swe.FLG_SPEED)
# lilith_se[0] = longitudine eclittica del giorno
```

**Attenzione**: pyswisseph `calc_ut` vuole JD in UT, non TT.
Convertire: `jd_ut = jd_tt - deltat`.

### Passo 4: Confronto triangolare

Per ogni data, calcolare:
```
err_noi    = |Ω_libephemeris - Ω_Horizons|
err_se     = |Ω_pyswisseph   - Ω_Horizons|

err_noi_lilith = |Lilith_libephemeris - Lilith_Horizons|
err_se_lilith  = |Lilith_pyswisseph   - Lilith_Horizons|
```

Dove `Lilith_Horizons = (Ω + ω + 180°) mod 360°`.

Output: tabella con tutte le differenze + statistiche (media, mediana, max, RMS).

### Passo 5: Verifica indipendente con vettori grezzi

Come ulteriore controllo, ricalcolare Ω e l'apogeo **direttamente dai vettori
di stato JPL** usando Skyfield raw, senza passare per le nostre funzioni
`calc_true_lunar_node` / `calc_true_lilith`:

```python
from skyfield.api import load
from skyfield.framelib import ecliptic_frame
import numpy as np

ts = load.timescale()
planets = load('de440.bsp')
earth = planets['earth']
moon = planets['moon']

t = ts.tt_jd(jd_tt)
pos = (moon - earth).at(t)
r = pos.frame_xyz(ecliptic_frame).au
v = pos.frame_xyz_and_velocity(ecliptic_frame)[1].au_per_d

# h = r × v
h = np.cross(r, v)

# Nodo: Ω = atan2(h_x, -h_y)
omega_raw = np.degrees(np.arctan2(h[0], -h[1])) % 360

# Eccentricità: e = (v × h)/μ - r/|r|
mu = 8.997011e-10  # GM_earth_moon in AU³/day²
r_mag = np.linalg.norm(r)
vxh = np.cross(v, h)
e_vec = vxh / mu - r / r_mag
apogee_lon_raw = np.degrees(np.arctan2(-e_vec[1], -e_vec[0])) % 360
```

Se `omega_raw` differisce dalla nostra funzione → bug nel nostro codice.
Se concorda → il nostro algoritmo è corretto, lo scarto con SE è dovuto al
loro modello lunare.

---

## Scenari Possibili e Azioni

### Scenario A: Noi ≈ Horizons, SE lontano
**Significato**: Noi siamo più precisi. SE usa un modello lunare meno accurato.
**Azione**:
- Documentare come vantaggio competitivo
- Stringere le tolleranze dei test (da 0.15° a ~0.01° o meno)
- Aggiungere test di confronto vs Horizons come reference

### Scenario B: SE ≈ Horizons, noi lontani
**Significato**: Abbiamo un bug nell'estrazione degli elementi osculanti.
**Azione**:
- Identificare la causa (frame? costanti? formula?)
- Fixare il nostro algoritmo
- Verificare che il fix riduca lo scarto sotto 1"

### Scenario C: Entrambi lontani da Horizons
**Significato**: Problema di frame di riferimento nel confronto (eclittica
del giorno vs J2000 vs ICRS) o Horizons usa una definizione diversa.
**Azione**:
- Verificare attentamente il frame di riferimento
- Ripetere con frame esplicito J2000 per tutti e tre
- Se il frame era il problema, il confronto diretto diventa valido

### Scenario D: Noi ≈ SE ≈ Horizons (lo scarto sparisce)
**Significato**: Lo scarto di 0.15° nei test era dovuto a date specifiche
o a un bug nel test di comparazione, non a un errore sistematico.
**Azione**:
- Analizzare perché i test mostrano 0.15° e il confronto diretto no
- Potenziale problema di conversione TT/UT nei test

---

## Rischi e Cose da Tenere d'Occhio

### 1. Frame di riferimento (CRITICO)
Il confronto è valido SOLO se tutti e tre usano lo stesso frame.
- **libephemeris**: eclittica vera del giorno (Skyfield `ecliptic_frame` = IAU 2006 precessione + IAU 2000A nutazione)
- **pyswisseph**: eclittica vera del giorno (con il suo modello di precessione/nutazione)
- **Horizons**: dipende da `REF_SYSTEM` — verificare se `ECLIPDATE` restituisce
  la stessa eclittica

**Mitigazione**: fare il confronto anche in J2000 (`SEFLG_J2000` per noi e SE,
`REF_SYSTEM=J2000` per Horizons). In J2000 il frame è fisso e non c'è ambiguità
di modello di precessione/nutazione.

### 2. Conversione TT/UT
- Horizons vuole JD in TT (per ELEMENTS)
- pyswisseph `calc_ut` vuole JD in UT
- Le nostre funzioni `calc_true_lunar_node` vogliono JD in TT
- Errore di 69s su Delta-T → errore di ~34" sulla Luna → potrebbe spiegare parte dello scarto

**Mitigazione**: usare lo stesso JD TT per tutti. Per pyswisseph, convertire
esplicitamente: `jd_ut = jd_tt - swe.deltat(jd_tt)`.

### 3. Definizione di "eclittica"
Horizons potrebbe usare una definizione leggermente diversa del piano eclittico
(es. eclittica media J2000 vs eclittica inerziale). La differenza è di ~0.5".

**Mitigazione**: confrontare prima le posizioni della Luna stessa (lon, lat) tra
le 3 fonti. Se concordano a <1", il frame è allineato.

### 4. μ (parametro gravitazionale)
Il nostro μ = G(M_Terra + M_Luna) usa costanti IAU 2015. SE potrebbe usare valori
leggermente diversi. Questo influenza il vettore eccentricità e quindi Lilith.

**Mitigazione**: testare la sensibilità variando μ di ±0.1% e misurando l'effetto
sulla longitudine di Lilith.

---

## Deliverable

1. **Script** `scripts/investigate_true_node_horizons.py` — eseguibile standalone
2. **Tabella risultati** — CSV con tutti i confronti
3. **Report** — chi ha ragione, con quanta confidenza, e cosa fare

## Stima Tempo

| Passo | Tempo |
|-------|-------|
| Preparazione date + script query Horizons | 1h |
| Calcolo con libephemeris + pyswisseph | 30min |
| Confronto + analisi | 1h |
| Verifica vettori grezzi | 30min |
| Fix se necessario (Scenario B) | 2-4h |
| **Totale** | **3-7h** |

---

## Bonus: Bug Stelle Fisse (da fixare comunque)

Indipendentemente dall'investigazione sopra, `fixed_stars.py` ha un bug di
correttezza: i flag `SEFLG_SIDEREAL`, `SEFLG_J2000`, `SEFLG_NONUT`, `SEFLG_TRUEPOS`,
`SEFLG_XYZ`, `SEFLG_RADIANS`, `SEFLG_BARYCTR`, `SEFLG_TOPOCTR`, `SEFLG_ICRS` sono
**tutti silenziosamente ignorati**. I soli flag attivi sono `SEFLG_SPEED`,
`SEFLG_NOABERR`, `SEFLG_NOGDEFL`, `SEFLG_EQUATORIAL`.

Il più grave è `SEFLG_SIDEREAL`: un utente che chiama `swe_fixstar_ut(star, jd,
SEFLG_SIDEREAL)` riceve coordinate **tropicali** senza nessun avviso.

Fix necessario: applicare la sottrazione dell'ayanamsha dopo il calcolo della
posizione, come viene fatto in `planets.py` per i pianeti.
