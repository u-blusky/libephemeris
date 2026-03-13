# LibEphemeris — Manuale Pratico di Efemeridi

## Piano di scrittura del manuale

Questo documento delinea la struttura completa del manuale. Ogni sezione e un task da implementare come capitolo markdown in `docs/manual/`. Il manuale e pensato come libro tascabile (~100-150 pagine), rivolto a chi parte da zero e vuole capire sia i concetti astronomici/astrologici sia come usare la libreria per calcoli reali.

**Formato finale**: Serie di file markdown in `docs/manual/`, uno per capitolo.
**Tono**: Didattico, con esempi concreti e riferimenti al codice.
**Lingua**: Italiano.

---

## Indice del manuale

| # | Capitolo | Pagine stimate |
|---|----------|----------------|
| 0 | Prefazione e come usare questo manuale | 3 |
| 1 | Il cielo visto dalla Terra | 8 |
| 2 | Misurare il tempo in astronomia | 10 |
| 3 | Dove si trova un pianeta: coordinate celesti | 12 |
| 4 | Le efemeridi: cosa sono e come funzionano | 8 |
| 5 | Calcolare la posizione dei pianeti | 10 |
| 6 | La Luna: nodi, apogeo, perigeo e Lilith | 10 |
| 7 | Le case astrologiche | 10 |
| 8 | Le stelle fisse | 6 |
| 9 | Eclissi solari e lunari | 10 |
| 10 | Alba, tramonto e visibilita | 8 |
| 11 | Lo zodiaco siderale e le ayanamsha | 8 |
| 12 | Corpi minori: asteroidi, centauri, TNO | 8 |
| 13 | Pianeti ipotetici e parti arabe | 6 |
| 14 | Precisione, performance e modalita LEB | 8 |
| 15 | Ricettario: calcoli pratici dalla A alla Z | 15 |
| | **Totale stimato** | **~140** |

---

## Dettaglio capitoli

### Capitolo 0 — Prefazione e come usare questo manuale (~3 pagine)

**Obiettivo**: Orientare il lettore, spiegare a chi serve il manuale e come e organizzato.

**Contenuti**:
- A chi e rivolto: sviluppatori, astrologi, curiosi di astronomia
- Cosa imparerai: i concetti dietro le efemeridi, non solo "come chiamare una funzione"
- Come e organizzato: ogni capitolo ha teoria, esempio di vita reale, esempio di codice
- Prerequisiti: Python base, nessuna conoscenza astronomica richiesta
- Installazione rapida: `pip install libephemeris`, primo calcolo in 3 righe
- Convenzione del manuale: icone/simboli per teoria, codice, vita reale

**Esempio di codice introduttivo**:
```python
import libephemeris as ephem
jd = ephem.julday(2024, 4, 8, 12.0)  # Eclissi solare 8 aprile 2024
pos, flag = ephem.calc_ut(jd, ephem.SE_SUN)
print(f"Sole a {pos[0]:.4f} gradi in Ariete")
```

**Riferimenti codebase**: `__init__.py` (API pubblica), `swe_julday`, `swe_calc_ut`

---

### Capitolo 1 — Il cielo visto dalla Terra (~8 pagine)

**Obiettivo**: Costruire l'intuizione geometrica di base. Il lettore deve "vedere" mentalmente la sfera celeste.

**Sezioni**:

1. **La sfera celeste** (~1.5 pag)
   - Immagina di essere dentro una sfera enorme con le stelle dipinte sopra
   - I poli celesti: prolungamento dell'asse terrestre
   - L'equatore celeste: prolungamento dell'equatore terrestre
   - Esempio vita reale: la Stella Polare indica il polo nord celeste
   - Perche le stelle "girano" durante la notte (rotazione terrestre)

2. **L'eclittica: il cammino del Sole** (~1.5 pag)
   - Il Sole non sta fermo sulla sfera celeste: si muove di ~1 grado al giorno
   - Il percorso annuale del Sole = eclittica
   - Perche e inclinata di 23.4 gradi rispetto all'equatore (obliquita)
   - Le stagioni nascono da questa inclinazione
   - Esempio vita reale: d'estate il Sole sale piu alto nel cielo
   - Codice: `ephem.calc_ut(jd, ephem.SE_ECL_NUT)` per ottenere l'obliquita

3. **Lo zodiaco: 12 settori da 30 gradi** (~1.5 pag)
   - L'eclittica divisa in 12 segni da 30 gradi ciascuno
   - Punto zero = equinozio di primavera (punto vernale, punto gamma)
   - Ariete 0-30, Toro 30-60, Gemelli 60-90... Pesci 330-360
   - Differenza tra segni (settori geometrici) e costellazioni (gruppi di stelle)
   - La precessione: il punto zero si sposta di ~50"/anno (anticipazione cap. 11)
   - Esempio vita reale: "Sole in Ariete" = il Sole e tra 0 e 30 gradi di longitudine eclittica

4. **L'orizzonte locale** (~1.5 pag)
   - Il cielo che vedi dipende da dove sei sulla Terra
   - Zenit (sopra la testa), nadir (sotto i piedi), orizzonte
   - Azimut: direzione sull'orizzonte (N, E, S, O)
   - Altezza: angolo sopra l'orizzonte (0 = orizzonte, 90 = zenit)
   - Convenzione astronomica vs navigazionale dell'azimut (S=0 vs N=0)
   - Codice: `ephem.azalt(jd, flag, geopos, atpress, attemp, pos)` per altezza/azimut

5. **I quattro angoli: Ascendente, MC, Discendente, IC** (~2 pag)
   - L'Ascendente: il grado dell'eclittica che sorge a est in quel momento
   - Il Medio Cielo (MC): il grado dell'eclittica al punto piu alto
   - Discendente e IC: i punti opposti
   - Perche cambiano ogni 4 minuti (~1 grado)
   - Esempio vita reale: il tuo segno ascendente dipende dall'ora esatta di nascita
   - Codice: `cusps, angles = ephem.houses(jd, lat, lon, b'P')`

**Riferimenti codebase**: `utils.py:azalt()`, `houses.py:swe_houses()`, `planets.py:swe_calc_ut()`

---

### Capitolo 2 — Misurare il tempo in astronomia (~10 pagine)

**Obiettivo**: Capire perche il tempo e complicato e come la libreria lo gestisce.

**Sezioni**:

1. **Il Giorno Giuliano (JD)** (~2 pag)
   - Il problema: date con calendari diversi (giuliano, gregoriano, islamico...)
   - Soluzione: un contatore unico che parte dal 1 gennaio 4713 a.C.
   - JD 2460408.5 = mezzogiorno del 8 aprile 2024
   - Perche .5: il giorno giuliano inizia a mezzogiorno, non a mezzanotte
   - Come convertire: `jd = ephem.julday(anno, mese, giorno, ore_decimali)`
   - Come riconvertire: `anno, mese, giorno, ore = ephem.revjul(jd)`
   - Il calendario giuliano vs gregoriano: `SE_GREG_CAL` vs `SE_JUL_CAL`
   - Esempio vita reale: calcolare quanti giorni mancano a un evento

2. **UT, TT e Delta-T** (~2.5 pag)
   - UT (Universal Time): basato sulla rotazione terrestre — il "tempo del meridiano"
   - La Terra rallenta: un giorno non e sempre uguale
   - TT (Terrestrial Time): tempo perfettamente uniforme, basato su orologi atomici
   - Delta-T = TT - UT: oggi circa 69 secondi, nel 1900 era ~3 secondi
   - Perche serve: i pianeti si muovono con il tempo fisico (TT), ma noi viviamo con UT
   - Funzioni `_ut` vs senza suffisso: `calc_ut(jd_ut, ...)` vs `calc(jd_tt, ...)`
   - Codice: `delta_t = ephem.deltat(jd)` — restituisce la differenza in giorni
   - Esempio vita reale: per una carta natale del 1800, Delta-T fa spostare la Luna di ~30"

3. **TAI, UTC e i secondi intercalari** (~1.5 pag)
   - UTC: il tempo dei nostri orologi, con i "leap seconds" aggiunti ogni tanto
   - TAI: tempo atomico puro, senza aggiustamenti — TAI = UTC + N secondi
   - TT = TAI + 32.184 secondi (costante fissa per sempre)
   - Codice: `ephem.utc_to_tai_jd(anno, mese, giorno, ore, min, sec)`
   - Esempio vita reale: perche il tuo orologio ogni tanto ha un secondo in piu

4. **Tempo siderale** (~1.5 pag)
   - Il giorno siderale: quanto ci mette la Terra a fare un giro rispetto alle stelle
   - 23h 56m 4s — circa 4 minuti piu corto del giorno solare
   - ARMC (Ascensione Retta del Medio Cielo): l'orologio siderale locale
   - Perche serve: le case astrologiche dipendono dal tempo siderale
   - Codice: `st = ephem.sidtime(jd)` — tempo siderale in ore
   - Esempio vita reale: perche le costellazioni "si alzano prima" ogni notte

5. **Tempo locale, fusi orari e LMT** (~1.5 pag)
   - LMT (Local Mean Time): il tempo solare medio per la tua longitudine esatta
   - La differenza tra fuso orario e tempo locale reale
   - Conversione: `ephem.lmt_to_lat(jd_lmt, longitudine)` (LMT -> tempo universale approssimato)
   - Equazione del tempo: perche il mezzogiorno solare non e sempre alle 12:00
   - Codice: `eot = ephem.time_equ(jd)` — restituisce la differenza in ore
   - Esempio vita reale: a Milano (lon 9.19) il mezzogiorno solare e ~37 minuti dopo quello di Greenwich

6. **IERS e Delta-T osservato** (~1 pag)
   - Il servizio IERS pubblica i valori misurati di Delta-T
   - Differenza tra Delta-T calcolato (polinomi storici) e Delta-T osservato (dati IERS)
   - Codice: `ephem.set_iers_delta_t_enabled(True)`, `ephem.download_delta_t_data()`
   - Quando serve: per massima precisione su date recenti/attuali

**Riferimenti codebase**: `time_utils.py` (tutto), `iers_data.py`, costanti `SE_GREG_CAL`/`SE_JUL_CAL`

---

### Capitolo 3 — Dove si trova un pianeta: coordinate celesti (~12 pagine)

**Obiettivo**: Capire tutti i sistemi di coordinate usati in astronomia e astrologia, e come convertirli.

**Sezioni**:

1. **Coordinate eclittiche: longitudine e latitudine** (~2 pag)
   - Il sistema piu usato in astrologia
   - Longitudine eclittica: 0-360 gradi lungo l'eclittica (il "segno zodiacale")
   - Latitudine eclittica: distanza dall'eclittica (-90 a +90)
   - Il Sole ha sempre latitudine ~0 (per definizione, e lui che traccia l'eclittica)
   - La Luna puo arrivare a +-5.1 gradi di latitudine
   - Codice: `pos, flag = ephem.calc_ut(jd, ephem.SE_MOON)` — pos[0]=lon, pos[1]=lat, pos[2]=dist
   - Esempio vita reale: "Luna a 15 gradi Cancro" = longitudine eclittica 105.x gradi

2. **Coordinate equatoriali: ascensione retta e declinazione** (~2 pag)
   - Il sistema usato dai telescopi
   - Ascensione Retta (RA): 0-360 gradi (o 0-24 ore) lungo l'equatore celeste
   - Declinazione (Dec): distanza dall'equatore celeste (-90 a +90)
   - Perche i telescopi preferiscono questo sistema (si allinea con la rotazione terrestre)
   - Codice: `pos, flag = ephem.calc_ut(jd, ephem.SE_MOON, ephem.SEFLG_EQUATORIAL)`
   - Conversione: `ephem.cotrans((lon, lat, dist), -obliquity)` eclittica -> equatoriale

3. **Coordinate orizzontali: altezza e azimut** (~1.5 pag)
   - Il sistema "pratico": dove guardo nel cielo?
   - Altezza: sopra l'orizzonte (> 0) o sotto (< 0)
   - Azimut: in che direzione sull'orizzonte
   - Servono latitudine, longitudine e altitudine dell'osservatore
   - Codice: `ephem.azalt(jd, ephem.SE_EQU2HOR, geopos, atpress, attemp, (ra, dec, dist))`
   - Esempio vita reale: "Giove e visibile stasera?" — calcola altezza, se > 0 e sopra l'orizzonte

4. **Sistemi di riferimento: ICRS, J2000, equinozio del giorno** (~2 pag)
   - Equinozio del giorno ("of date"): le coordinate tengono conto della precessione attuale
   - J2000: le coordinate sono fissate al 1 gennaio 2000 — il "default"
   - ICRS: il sistema moderno basato su quasar lontanissimi, quasi identico a J2000
   - La differenza tra ICRS e J2000: il "frame bias" di ~0.02 secondi d'arco
   - Quando serve quale: astrologia usa "of date", cataloghi astronomici usano J2000/ICRS
   - Flag: `SEFLG_J2000`, `SEFLG_ICRS`, nessun flag = of date
   - Codice: `ephem.calc_ut(jd, body, ephem.SEFLG_EQUATORIAL | ephem.SEFLG_ICRS)`

5. **Precessione, nutazione e aberrazione** (~2 pag)
   - **Precessione**: l'asse terrestre disegna un cerchio in ~26000 anni (come una trottola)
   - Effetto: il punto vernale si sposta di ~50"/anno
   - **Nutazione**: l'asse terrestre "oscilla" leggermente (periodo principale ~18.6 anni)
   - Ampiezza: ~9" in longitudine, ~17" in obliquita
   - Flag `SEFLG_NONUT`: calcola senza nutazione (coordinate "medie" invece di "vere")
   - **Aberrazione**: la luce impiega tempo ad arrivare, e noi ci muoviamo
   - Effetto: fino a ~20" di spostamento apparente
   - Flag `SEFLG_NOABERR`, `SEFLG_TRUEPOS`: per posizioni geometriche pure
   - Codice: confronto tra `calc_ut(jd, body, 0)` e `calc_ut(jd, body, SEFLG_NONUT)`

6. **Coordinate cartesiane (XYZ)** (~1 pag)
   - A volte servono X, Y, Z in unita astronomiche invece di angoli
   - Flag `SEFLG_XYZ`: restituisce (x, y, z, vx, vy, vz)
   - Utile per: calcoli di distanza 3D, angoli di fase, problemi geometrici
   - Codice: `ephem.calc_ut(jd, body, ephem.SEFLG_XYZ | ephem.SEFLG_EQUATORIAL)`

7. **Convertire tra sistemi** (~1.5 pag)
   - `cotrans()`: conversione eclittica <-> equatoriale
   - `cotrans_sp()`: conversione con velocita
   - `azalt()`: eclittica/equatoriale -> orizzontale
   - `azalt_rev()`: orizzontale -> eclittica/equatoriale
   - Tabella riassuntiva: quale flag per quale output
   - Esempio completo: data una posizione eclittica, trovare altezza e azimut

**Riferimenti codebase**: `planets.py:_calc_body()`, `utils.py:cotrans()`, `utils.py:azalt()`, costanti `SEFLG_*`

---

### Capitolo 4 — Le efemeridi: cosa sono e come funzionano (~8 pagine)

**Obiettivo**: Capire cosa c'e "dentro" una efemeride, da dove vengono i dati e come la libreria li usa.

**Sezioni**:

1. **Cos'e un'efemeride** (~1.5 pag)
   - Definizione: una tabella che dice "a questa data, questo corpo celeste e qui"
   - Dalle tavolette babilonesi alle effemeridi digitali: 3000 anni di storia in 1 pagina
   - Efemeride cartacea vs digitale: stessa idea, precisione diversa
   - Esempio vita reale: l'Almanacco Nautico usato dai marinai per navigare

2. **JPL DE440: la nostra sorgente dati** (~2 pag)
   - NASA Jet Propulsion Laboratory: chi sono e cosa fanno
   - DE440: integrazione numerica del Sistema Solare per 800 anni (1550-2650)
   - DE441: versione estesa per 30000 anni (-13200 a +17191)
   - DE440s: versione compatta per 300 anni (1849-2150)
   - Come funziona internamente: polinomi di Chebyshev per ogni corpo, ogni 8-32 giorni
   - Precisione: sub-millimetrica per i pianeti interni, pochi km per Plutone
   - Skyfield: la libreria Python che legge i file JPL e fa i conti
   - Codice: `ephem.set_precision_tier('medium')` — sceglie DE440 (default)

3. **I tre livelli di precisione** (~1.5 pag)
   - `base`: DE440s — leggero (17 MB), 1849-2150, ideale per astrologia moderna
   - `medium`: DE440 — bilanciato (114 MB), 1550-2650, default
   - `extended`: DE441 — massimo (3.1 GB), -13200 a +17191, per ricerca storica
   - Come scegliere: `ephem.set_precision_tier('extended')`
   - Codice: `ephem.download_for_tier('medium')` — scarica i file necessari

4. **Geocentrico, eliocentrico e baricentrico** (~1.5 pag)
   - **Geocentrico** (default): visto dalla Terra — quello che serve in astrologia
   - **Eliocentrico**: visto dal Sole — flag `SEFLG_HELCTR`
   - **Baricentrico**: visto dal centro di massa del Sistema Solare — flag `SEFLG_BARYCTR`
   - La differenza tra Sole e baricentro: il Sole "oscilla" di ~2 raggi solari
   - **Planetocentrico**: visto da un altro pianeta — `calc_pctr(jd, body, center, flag)`
   - Esempio vita reale: Marte visto da Giove (per curiosita)
   - Codice confronto: stessa data, stessa Luna, tre punti di vista diversi

5. **Topocentrico: la posizione conta** (~1.5 pag)
   - Geocentrico = dal centro della Terra; topocentrico = dalla tua posizione esatta
   - La differenza: fino a ~1 grado per la Luna (parallasse lunare)
   - Quando serve: eclissi, occultazioni, transiti — la geometria dipende da dove sei
   - Setup: `ephem.set_topo(lat, lon, altitudine_metri)`
   - Flag: `SEFLG_TOPOCTR`
   - Esempio vita reale: un'eclissi solare e totale a Milano ma parziale a Roma

**Riferimenti codebase**: `planets.py:swe_calc_ut()`, `planets.py:swe_calc_pctr()`, `state.py:set_topo()`, `state.py:PrecisionTier`

---

### Capitolo 5 — Calcolare la posizione dei pianeti (~10 pagine)

**Obiettivo**: Padroneggiare `calc_ut` e le funzioni correlate. Capire cosa restituisce e tutti i flag.

**Sezioni**:

1. **La funzione principale: `calc_ut`** (~2 pag)
   - Firma: `calc_ut(jd_ut, body, flag) -> (tuple_6, flag_restituito)`
   - La tupla: `(longitudine, latitudine, distanza, vel_lon, vel_lat, vel_dist)`
   - Unita: gradi, gradi, UA, gradi/giorno, gradi/giorno, UA/giorno
   - Il body: `SE_SUN=0`, `SE_MOON=1`, `SE_MERCURY=2`, ... `SE_PLUTO=9`
   - `SE_MEAN_NODE=10`, `SE_TRUE_NODE=11`, `SE_MEAN_APOG=12`, `SE_OSCU_APOG=13`
   - `SE_CHIRON=15`, `SE_PHOLUS=16`, `SE_CERES=17`...
   - `calc_ut` vs `calc`: la prima vuole JD in UT, la seconda in TT
   - Codice: calcolare dove sono tutti i pianeti adesso

2. **I flag di calcolo: la tabella completa** (~2.5 pag)
   - Tabella con tutti i `SEFLG_*`, raggruppati per categoria:
     - **Coordinate**: `EQUATORIAL`, `XYZ`, `RADIANS`
     - **Centro**: `HELCTR`, `BARYCTR`, `TOPOCTR`
     - **Frame**: `J2000`, `ICRS`, `SIDEREAL`, `NONUT`
     - **Correzioni**: `NOABERR`, `NOGDEFL`, `TRUEPOS`, `SPEED`
   - Combinazioni comuni: `SEFLG_SPEED` (quasi sempre), `SEFLG_EQUATORIAL | SEFLG_SPEED`
   - Esempio: stessa Luna con 5 flag diversi — confronto dei risultati

3. **Velocita e moto retrogrado** (~1.5 pag)
   - La velocita in longitudine: positiva = diretto, negativa = retrogrado
   - Cos'e il moto retrogrado: un effetto ottico, come sorpassare un'auto
   - Le stazioni: quando la velocita e zero (il pianeta "si ferma")
   - Codice: `vel = pos[3]` — controllare se `vel < 0`
   - Funzioni dedicate: `is_retrograde()`, `get_station_info()`, `swe_find_station_ut()`
   - Esempio vita reale: "Mercurio retrogrado" — calcolare quando inizia e finisce

4. **Fenomeni planetari: magnitudine, fase, elongazione** (~2 pag)
   - `pheno_ut(jd, body, flag)` restituisce: angolo di fase, fase (0-1), elongazione, diametro apparente, magnitudine
   - **Angolo di fase**: angolo Sole-Pianeta-Terra — serve per la luminosita
   - **Elongazione**: distanza angolare dal Sole — il pianeta e visibile?
   - **Magnitudine**: quanto e luminoso? (numeri piu bassi = piu luminoso)
   - Funzioni helper: `is_morning_star()`, `is_evening_star()`, `get_elongation_type()`
   - Esempio vita reale: "Venere e visibile stasera?" — elongazione > 10 gradi e altezza > 0

5. **Elementi orbitali** (~1 pag)
   - `get_orbital_elements(jd, body, flag)`: restituisce i 6 elementi kepleriani
   - Semi-asse maggiore, eccentricita, inclinazione, nodo, argomento del perielio, anomalia media
   - Distanza massima e minima: `orbit_max_min_true_distance()`
   - Esempio: confrontare l'orbita di Marte (e=0.093) con quella della Terra (e=0.017)

6. **Nodi e apsidi** (~1 pag)
   - `nod_aps_ut(jd, body, flag, method)`: nodo ascendente, discendente, perielio, afelio
   - Differenza tra nodo medio e nodo vero (oscillazioni a breve termine)
   - Differenza tra apsidi medie e osculanti
   - Esempio: il perielio della Terra (inizio gennaio) — calcolare la data esatta

**Riferimenti codebase**: `planets.py` (tutto), `crossing.py:is_retrograde()`, costanti body `SE_*`

---

### Capitolo 6 — La Luna: nodi, apogeo, perigeo e Lilith (~10 pagine)

**Obiettivo**: Capire i punti lunari speciali, fondamentali in astrologia, e la complessita del moto lunare.

**Sezioni**:

1. **Perche la Luna e speciale** (~1.5 pag)
   - L'orbita lunare e la piu complessa del Sistema Solare (dopo Mercurio)
   - Perturbata dal Sole: l'orbita cambia forma, inclinazione e orientamento
   - Velocita: ~12-15 gradi al giorno (30x piu veloce dei pianeti)
   - Importanza in astrologia: la Luna si muove cosi tanto che l'ora di nascita conta molto
   - La parallasse lunare: ~1 grado — l'unico corpo per cui topocentrico vs geocentrico fa una differenza enorme

2. **I nodi lunari: dove la Luna attraversa l'eclittica** (~2 pag)
   - L'orbita lunare e inclinata di ~5.1 gradi rispetto all'eclittica
   - Nodo ascendente: dove la Luna sale sopra l'eclittica (da sud a nord)
   - Nodo discendente: dove la Luna scende sotto (sempre opposto, +180 gradi)
   - **Nodo medio**: posizione calcolata con un polinomio regolare (~18.6 anni per un giro)
   - **Nodo vero**: posizione reale con tutte le oscillazioni (~1.5 gradi di ampiezza)
   - Perche i nodi contano: le eclissi avvengono solo vicino ai nodi
   - Codice: `ephem.calc_ut(jd, ephem.SE_MEAN_NODE)` vs `ephem.calc_ut(jd, ephem.SE_TRUE_NODE)`
   - Codice: `ephem.calc_mean_lunar_node(jd_tt)`, `ephem.calc_true_lunar_node(jd_tt)`
   - Esempio vita reale: il nodo nord nel tema natale e dove "devi andare"

3. **Apogeo e perigeo: la Luna vicina e lontana** (~2 pag)
   - L'orbita lunare e un'ellisse: a volte la Luna e piu vicina (perigeo), a volte piu lontana (apogeo)
   - Distanza: ~356000 km (perigeo) a ~407000 km (apogeo)
   - Effetto visibile: la "superluna" (perigeo + luna piena) appare ~14% piu grande
   - **Apogeo medio (Lilith Nera media)**: vedi sezione successiva
   - **Apogeo osculante**: posizione istantanea dell'apogeo tenendo conto delle perturbazioni
   - **Apogeo interpolato**: lisciato per eliminare le oscillazioni troppo rapide
   - Codice: `ephem.calc_ut(jd, ephem.SE_MEAN_APOG)` — apogeo medio
   - Codice: `ephem.calc_ut(jd, ephem.SE_OSCU_APOG)` — apogeo osculante
   - Codice: `ephem.calc_interpolated_apogee(jd_tt)` — apogeo interpolato
   - Esempio vita reale: prevedere la prossima "superluna"

4. **Lilith: la Luna Nera** (~2 pag)
   - **Lilith media** (la piu usata in astrologia): = apogeo medio dell'orbita lunare
   - L'apogeo medio fa un giro completo dello zodiaco in ~8.85 anni
   - **Lilith vera/osculante**: l'apogeo istantaneo, con oscillazioni di diversi gradi
   - **Lilith interpolata**: via di mezzo — lisciata ma piu precisa della media
   - Differenza tipica: media vs vera puo essere 20-30 gradi
   - Codice: `ephem.calc_mean_lilith(jd_tt)` — solo longitudine
   - Codice: `ephem.calc_true_lilith(jd_tt)` — longitudine, velocita, distanza
   - **Luna Bianca (Selena)**: il punto opposto a Lilith Nera (apogeo -> perigeo)
   - Codice: `ephem.calc_white_moon_position(jd_tt)` — opposto alla Lilith media
   - Perche "Lilith": il mito della prima moglie di Adamo, associata al lato ombra

5. **Perigeo interpolato e calibrazione** (~1 pag)
   - Il perigeo ha le stesse varianti: medio, osculante, interpolato
   - `ephem.calc_interpolated_perigee(jd_tt)` — per i transiti di precisione
   - La calibrazione: come abbiamo migliorato la precisione (cenno al workflow in AGENTS.md)
   - Esempio: distanza della Luna al perigeo vs distanza all'apogeo

6. **Attraversamenti del nodo** (~1.5 pag)
   - Trovare il momento esatto in cui la Luna attraversa il suo nodo
   - Codice: `ephem.mooncross_node_ut(jd_start, flag)` — restituisce (jd, xlon, xlat)
   - La questione TT/UT: perche il nostro risultato e piu preciso (il bug Delta-T trovato in pyswisseph)
   - Esempio vita reale: i nodi lunari e le eclissi — se Luna piena vicino al nodo = eclissi lunare

**Riferimenti codebase**: `lunar.py` (tutto), `crossing.py:swe_mooncross_node_ut()`, costanti `SE_MEAN_NODE`, `SE_TRUE_NODE`, `SE_MEAN_APOG`, `SE_OSCU_APOG`

---

### Capitolo 7 — Le case astrologiche (~10 pagine)

**Obiettivo**: Capire cosa sono le case, come funzionano i diversi sistemi, e come calcolarle.

**Sezioni**:

1. **Cosa sono le case** (~1.5 pag)
   - Le case dividono il cielo in 12 settori, come fette di una torta
   - Ogni casa "governa" un ambito della vita (1a = se stessi, 7a = relazioni, 10a = carriera...)
   - La differenza fondamentale tra segni (divisione dell'eclittica) e case (divisione dello spazio locale)
   - Le case dipendono da: data, ora, e luogo — i segni dipendono solo dalla data
   - Cuspide: il punto di inizio di ogni casa
   - Esempio vita reale: due persone nate lo stesso giorno ma a ore diverse hanno case diverse

2. **L'Ascendente e il Medio Cielo** (~1.5 pag)
   - L'Ascendente (ASC): dove l'eclittica interseca l'orizzonte a est — la cuspide della 1a casa
   - Il Medio Cielo (MC): dove l'eclittica interseca il meridiano superiore — la cuspide della 10a casa
   - Il Discendente (DSC): opposto all'ASC — cuspide della 7a
   - L'Imum Coeli (IC): opposto al MC — cuspide della 4a
   - Il Vertice (Vertex): dove il primo verticale interseca l'eclittica a ovest
   - Codice: `cusps, ascmc = ephem.houses(jd, lat, lon, b'P')`
   - `ascmc`: (ASC, MC, ARMC, Vertex, equatorial ASC, co-ASC Koch, co-ASC Munkasey, polar ASC)

3. **I sistemi di case: quale scegliere?** (~3 pag)
   - Tabella completa dei 20+ sistemi supportati con lettera identificativa:
     - `P` Placidus — il piu diffuso in Occidente, divide il tempo
     - `K` Koch — simile a Placidus, usa il semi-arco del MC
     - `R` Regiomontanus — divide l'equatore celeste
     - `C` Campanus — divide il primo verticale
     - `E` Uguale dall'ASC — tutte le case di 30 gradi
     - `W` Segno intero — ogni segno = una casa
     - `O` Porfirio — divide proporzionalmente i quadranti
     - `B` Alcabizio — divide i semi-archi
     - `M` Morino — divide l'equatore dal MC
     - `T` Polich/Page (topocentriche) — proiezione topografica
     - `G` Gauquelin (36 settori) — per ricerca statistica
     - `V` Vehlow — uguale ma spostato di 15 gradi
     - `X` Meridiano — basato sul meridiano
     - `F` Carter — pari-passu
     - `S` Sripati — media tra Porfirio e uguale
     - `L`/`Q` Pullen — sinusoidi/ratio
   - Quando usare quale: tradizione, tipo di astrologia, latitudine
   - Il problema delle latitudini estreme: Placidus e Koch non funzionano sopra ~66 gradi
   - Codice: `ephem.house_name(ord('P'))` — restituisce "Placidus"

4. **Posizione di un pianeta nelle case** (~2 pag)
   - `house_pos(armc, lat, obliquity, hsys, (lon, lat_ecl))` — restituisce posizione 1.0-12.999
   - Il numero intero = casa, la parte decimale = posizione nella casa
   - Esempio: 7.5 = meta della 7a casa
   - Differenza tra metodi: Placidus usa il semi-arco del corpo, Koch usa quello del MC
   - Il problema della latitudine eclittica: non tutti i sistemi la considerano
   - Settori Gauquelin: 36 settori per analisi statistica delle "zone Gauquelin"
   - Codice: `ephem.gauquelin_sector(jd, body, flag, method, geopos, atpress, attemp)`
   - Esempio vita reale: in che casa cade il tuo Sole?

5. **Latitudini estreme e circolo polare** (~1 pag)
   - Sopra ~66.5 gradi: il Sole non sorge o non tramonta per giorni
   - Placidus/Koch: impossibile calcolare — la libreria lancia `PolarCircleError`
   - Fallback: `houses_with_fallback()` — prova Placidus, se fallisce usa Porfirio
   - Codice: `ephem.get_extreme_latitude_info(lat)` — dice quali sistemi funzionano
   - Esempio: calcolare le case a Tromso (69.6 N) — Placidus fallisce, Porfirio funziona

**Riferimenti codebase**: `houses.py` (tutto), `houses.py:swe_house_pos()`, `houses.py:swe_houses()`, costanti nelle funzioni `_houses_*`

---

### Capitolo 8 — Le stelle fisse (~6 pagine)

**Obiettivo**: Capire come le stelle fisse sono catalogate, come si muovono, e come usarle nei calcoli.

**Sezioni**:

1. **Le stelle non sono cosi fisse** (~1.5 pag)
   - Le stelle si muovono, ma molto lentamente: il "moto proprio"
   - Esempio: Sirio si sposta di ~1.3"/anno — in 2000 anni, quasi mezzo grado
   - Due componenti: moto proprio in RA e in Dec (come latitudine/longitudine sulla sfera)
   - La precessione sposta TUTTE le stelle di ~50"/anno (effetto molto piu grande del moto proprio)
   - Codice: `ephem.propagate_proper_motion(ra, dec, pm_ra, pm_dec, parallax, rv, jd_from, jd_to)`

2. **Il catalogo: Hipparcos e oltre** (~1.5 pag)
   - Il satellite Hipparcos (ESA, 1989-1993): ha misurato 118218 stelle con precisione sub-mas
   - La libreria include un catalogo completo con nomi tradizionali, designazioni Bayer (alpha Centauri) e Flamsteed (61 Cygni)
   - Cercare una stella: per nome tradizionale, per numero HIP, per designazione
   - Codice: `pos, name, flag = ephem.fixstar2_ut(jd, "Sirius")` — cerca per nome
   - Codice: `pos, name, flag = ephem.fixstar2_ut(jd, ",HIP70890")` — cerca per HIP
   - Il ritorno: posizione eclittica (lon, lat, dist, vel_lon, vel_lat, vel_dist) + nome formattato
   - Magnitudine: `nome, mag, err = ephem.fixstar2_mag("Regulus")`

3. **Stelle in astrologia** (~1.5 pag)
   - Le stelle fisse "classiche": le 15 stelle di Agrippa, le stelle regali (Aldebaran, Regulus, Antares, Fomalhaut)
   - Congiunzione con un pianeta: quando un pianeta e entro ~1 grado di una stella fissa
   - Orbe tipico per le stelle fisse: 1-2 gradi (molto piu stretto dei pianeti)
   - Esempio vita reale: Regulus e a ~0 gradi Vergine — un pianeta a 29 Leone "congiunge" Regulus
   - Codice: calcolare se oggi un pianeta qualsiasi e congiunto a una stella reale

4. **Ricerca fuzzy e risoluzione nomi** (~1.5 pag)
   - La libreria accetta nomi in molte forme: "Sirius", "sirius", "Alpha CMa", "alpha Canis Majoris"
   - Ricerca fonetica: trova la stella anche con errori di spelling
   - Codice: `hip = ephem.get_hip_from_star_name("Aldebaran")` — restituisce numero HIP
   - La funzione `resolve_star_name()` per ricerche avanzate
   - Esempio: cercare "Betelgeuze" (con la z) — la libreria trova comunque Betelgeuse

**Riferimenti codebase**: `fixed_stars.py` (tutto), catalogo interno con dati Hipparcos/ESA

---

### Capitolo 9 — Eclissi solari e lunari (~10 pagine)

**Obiettivo**: Capire la geometria delle eclissi, i tipi, e come calcolarle con la libreria.

**Sezioni**:

1. **Come avviene un'eclissi** (~1.5 pag)
   - Eclissi solare: la Luna passa davanti al Sole (Luna nuova vicino a un nodo)
   - Eclissi lunare: la Terra blocca la luce del Sole sulla Luna (Luna piena vicino a un nodo)
   - Perche non succede ogni mese: l'orbita lunare e inclinata di 5.1 gradi
   - La "stagione delle eclissi": finestra di ~35 giorni quando i nodi sono allineati con il Sole
   - Ogni anno: 2-5 eclissi solari, 0-3 eclissi lunari
   - Esempio vita reale: l'eclissi totale di Sole dell'8 aprile 2024 in USA

2. **Tipi di eclissi** (~1.5 pag)
   - **Solare totale** (`SE_ECL_TOTAL`): la Luna copre tutto il Sole — la corona diventa visibile
   - **Solare anulare** (`SE_ECL_ANNULAR`): la Luna e troppo lontana, resta un anello di fuoco
   - **Solare parziale** (`SE_ECL_PARTIAL`): la Luna copre solo una parte del Sole
   - **Solare ibrida** (`SE_ECL_ANNULAR_TOTAL`): anulare in certi punti, totale in altri
   - **Lunare totale** (`SE_ECL_TOTAL`): la Luna entra completamente nell'ombra — diventa rossa
   - **Lunare parziale** (`SE_ECL_PARTIAL`): solo una parte della Luna nell'ombra
   - **Lunare penombrale** (`SE_ECL_PENUMBRAL`): la Luna nella penombra — quasi invisibile
   - Magnitudine: quanta parte del diametro e coperta (0 = niente, 1 = totale)
   - Oscuramento: quanta superficie e coperta — puo superare 1.0 per eclissi totali

3. **Trovare la prossima eclissi** (~2 pag)
   - Eclissi solare globale: `sol_eclipse_when_glob(jd_start, flag, ecl_type)` — cerca nel mondo
   - Eclissi solare locale: `sol_eclipse_when_loc(jd_start, flag, geopos)` — visibile da qui?
   - Eclissi lunare: `lun_eclipse_when(jd_start, flag, ecl_type)` — cerca globalmente
   - Eclissi lunare locale: `lun_eclipse_when_loc(jd_start, flag, geopos)` — visibile da qui?
   - Il risultato: tupla con tempi dei vari contatti (inizio, massimo, fine)
   - Codice completo: trovare le prossime 5 eclissi solari visibili da Roma
   - Esempio vita reale: "quando sara la prossima eclissi visibile da casa mia?"

4. **Dettagli di un'eclissi** (~2 pag)
   - `sol_eclipse_how(jd, flag, geopos)` — tipo e magnitudine in un dato momento e luogo
   - `sol_eclipse_where(jd, flag)` — dove nel mondo e visibile al massimo
   - `lun_eclipse_how(jd, flag, geopos)` — dettagli eclissi lunare
   - I 4 contatti solari: C1 (inizio), C2 (totalita inizio), C3 (totalita fine), C4 (fine)
   - I 6 contatti lunari: P1, U1, U2 (penombra, ombra inizio, totalita inizio), U3, U4, P4
   - Larghezza del percorso: `calc_eclipse_path_width(jd)` — in km
   - Linea centrale: `calc_eclipse_central_line(jd)` — coordinate geografiche
   - Limiti nord/sud: `calc_eclipse_northern_limit(jd)`, `calc_eclipse_southern_limit(jd)`
   - Elementi besseliani: la parametrizzazione matematica classica dell'eclissi

5. **Saros e Inex: i cicli delle eclissi** (~1.5 pag)
   - Il ciclo di Saros: 6585.32 giorni (~18 anni, 11 giorni, 8 ore)
   - Ogni Saros un'eclissi simile si ripete, spostata di ~120 gradi in longitudine
   - Serie Saros: una "famiglia" di eclissi correlate (es. Serie 145: include l'eclissi USA 2024)
   - Il ciclo Inex: 10571.95 giorni (~29 anni) — collega eclissi di serie Saros diverse
   - Codice: `ephem.get_saros_number(jd, ecl_type)` — restituisce (serie, membro)
   - Codice: `ephem.get_inex_number(jd, ecl_type)` — restituisce numero Inex
   - Esempio vita reale: l'eclissi del 2024 e la sorella del 2006 (stessa serie Saros, membro successivo)

6. **Occultazioni stellari e planetarie** (~1.5 pag)
   - La Luna puo passare davanti a una stella o un pianeta: occultazione
   - `lun_occult_when_glob(jd, body, flag)` — trova la prossima occultazione
   - `lun_occult_when_loc(jd, body, flag, geopos)` — visibile da qui?
   - Occultazioni planetarie: `planet_occult_when_glob()`, `planet_occult_when_loc()`
   - Esempio vita reale: la Luna che occulta Giove — sparisce per un'ora dietro il disco lunare

**Riferimenti codebase**: `eclipse.py` (tutto), costanti `SE_ECL_*`, `SAROS_CYCLE_DAYS`, `INEX_CYCLE_DAYS`

---

### Capitolo 10 — Alba, tramonto e visibilita (~8 pagine)

**Obiettivo**: Calcolare quando un corpo celeste sorge, transita e tramonta, e quando e visibile a occhio nudo.

**Sezioni**:

1. **Alba e tramonto** (~2 pag)
   - La definizione astronomica: il centro del disco tocca l'orizzonte
   - La rifrazione atmosferica: il Sole appare ~34' piu alto di quanto sia realmente
   - Il Sole "sorge" quando geometricamente e ancora sotto l'orizzonte
   - `rise_trans(jd, body, flag, rsmi, geopos, atpress, attemp)` — trova il prossimo evento
   - I flag `rsmi`: `SE_CALC_RISE=1`, `SE_CALC_SET=2`, `SE_CALC_MTRANSIT=4`, `SE_CALC_ITRANSIT=8`
   - `rise_trans_true_hor()` — con orizzonte reale (montagne, mare)
   - Codice completo: alba e tramonto del Sole oggi a Milano
   - Esempio vita reale: a che ora sorge il Sole il giorno del tuo compleanno?

2. **Transiti al meridiano** (~1 pag)
   - Transito superiore (culminazione): il corpo raggiunge il punto piu alto nel cielo
   - Transito inferiore: il corpo raggiunge il punto piu basso (sotto l'orizzonte, di solito)
   - Il transito del Sole al meridiano = mezzogiorno solare (non sempre alle 12:00!)
   - Codice: `rise_trans(jd, body, flag, SE_CALC_MTRANSIT, geopos, ...)`
   - Esempio: calcolare il mezzogiorno solare esatto per una data e un luogo

3. **I crepuscoli** (~1.5 pag)
   - Crepuscolo civile: Sole tra 0 e -6 gradi — si legge ancora il giornale
   - Crepuscolo nautico: Sole tra -6 e -12 gradi — si vede l'orizzonte marino
   - Crepuscolo astronomico: Sole tra -12 e -18 gradi — iniziano a comparire le stelle deboli
   - Notte astronomica: Sole sotto -18 gradi — cielo completamente buio
   - Codice: usare `rise_trans` con `disc_center` e un angolo orizzontale personalizzato
   - Costanti: `TWILIGHT_CIVIL_START`, `TWILIGHT_NAUTICAL_END`, `TWILIGHT_ASTRONOMICAL_END`
   - Esempio vita reale: a che ora inizia la notte astronomica stasera? (per osservare le stelle)

4. **Rifrazione atmosferica** (~1 pag)
   - L'atmosfera curva la luce: gli oggetti vicini all'orizzonte appaiono piu alti
   - A 0 gradi di altezza: ~34 arcminuti di rifrazione
   - A 10 gradi: ~5 arcminuti. A 45 gradi: ~1 arcminuto. Allo zenit: zero
   - Dipende da pressione e temperatura
   - `refrac(alt, atpress, attemp, SE_TRUE_TO_APP)` — da altezza vera ad apparente
   - `refrac(alt, atpress, attemp, SE_APP_TO_TRUE)` — da apparente a vera
   - `refrac_extended(alt, alt_geo, atpress, attemp, lapse_rate, flag)` — modello avanzato

5. **Visibilita eliacale** (~2.5 pag)
   - **Levata eliaca**: la prima volta che una stella/pianeta e visibile all'alba dopo la congiunzione col Sole
   - Importanza storica: la levata eliaca di Sirio segnava l'inizio dell'anno egizio
   - Modello di Schaefer: luminosita del cielo, estinzione atmosferica, acuita visiva dell'osservatore
   - `heliacal_ut(jd, geopos, atmo, observer, object_name, event_type, flag)` — trova la data
   - Tipi di evento: levata eliaca mattutina, tramonto eliaco serale, levata acronica, tramonto cosmico
   - `heliacal_pheno_ut()` — dettagli fenomenologici (magnitudine, altezza, azimut, ...)
   - `vis_limit_mag()` — magnitudine limite di visibilita per date condizioni
   - Estinzione atmosferica: `calc_airmass()`, `calc_extinction_magnitude()`
   - Modello di brillanza del cielo al crepuscolo: `calc_twilight_sky_brightness()`
   - Esempio vita reale: quando sara visibile Venere come "stella del mattino" quest'anno?

**Riferimenti codebase**: `eclipse.py:rise_trans()`, `heliacal.py`, `extinction.py`, `schaefer.py`, `utils.py:refrac()`

---

### Capitolo 11 — Lo zodiaco siderale e le ayanamsha (~8 pagine)

**Obiettivo**: Capire la differenza tra zodiaco tropicale e siderale, e le diverse scuole di ayanamsha.

**Sezioni**:

1. **Due zodiaci: tropicale e siderale** (~2 pag)
   - Zodiaco tropicale (Occidente): il punto zero e l'equinozio di primavera
   - Zodiaco siderale (India, astrologia vedica): il punto zero e legato alle stelle fisse
   - 2000 anni fa coincidevano — oggi sono sfasati di ~24 gradi
   - La precessione degli equinozi: il punto vernale si sposta di ~50.3"/anno rispetto alle stelle
   - Effetto pratico: se in tropicale sei "Ariete", in siderale potresti essere "Pesci"
   - Esempio vita reale: confronto tra un tema natale tropicale e lo stesso in siderale

2. **Cos'e l'ayanamsha** (~1.5 pag)
   - Ayanamsha = lo sfasamento attuale tra punto vernale e punto zero siderale
   - Oggi circa 24 gradi — cresce di ~50"/anno
   - Per convertire: longitudine siderale = longitudine tropicale - ayanamsha
   - Codice: `ephem.get_ayanamsa_ut(jd)` — restituisce l'ayanamsha in gradi
   - Codice: `ephem.set_sid_mode(ephem.SE_SIDM_LAHIRI)` — imposta il sistema siderale
   - Flag `SEFLG_SIDEREAL`: per avere le posizioni direttamente in siderale
   - Esempio: il Sole oggi in tropicale e in siderale Lahiri — la differenza e l'ayanamsha

3. **Le scuole di ayanamsha** (~2.5 pag)
   - Il problema: nessuno e d'accordo su dove mettere il punto zero siderale
   - **Lahiri** (`SE_SIDM_LAHIRI`): la piu usata, standard del governo indiano — Spica a 0 Bilancia
   - **Fagan-Bradley** (`SE_SIDM_FAGAN_BRADLEY`): usata in astrologia siderale occidentale
   - **Raman** (`SE_SIDM_RAMAN`): popolare in India del sud
   - **Krishnamurti** (`SE_SIDM_KRISHNAMURTI`): per il sistema KP
   - **Varianti Lahiri**: `SE_SIDM_LAHIRI_1940` (commissione 1940), `SE_SIDM_LAHIRI_VP285`, `SE_SIDM_LAHIRI_ICRC`
   - **Basate su stelle vere**: `SE_SIDM_TRUE_CITRA` (Spica a 180), `SE_SIDM_TRUE_REVATI`, `SE_SIDM_TRUE_PUSHYA`
   - **Babilonesi**: Kugler 1/2/3, Huber, ETPSC, Britton — per ricerca storica
   - **Galattiche**: centro galattico a 0 Sagittario, equatore galattico, allineamenti galattici
   - **User-defined**: `SE_SIDM_USER` con `set_sid_mode(255, t0, ayan_t0)` — definisci il tuo
   - Tabella completa dei 47 modi supportati con valore J2000

4. **Ayanamsha basate su stelle vere** (~1 pag)
   - Le "true star" ayanamsha usano la posizione reale (con moto proprio) di una stella di riferimento
   - `SE_SIDM_TRUE_CITRA`: posizione vera di Spica fissata a 180 gradi
   - Vantaggio: seguono il moto proprio reale della stella
   - Svantaggio: la stella potrebbe avere un moto proprio impreciso
   - Codice: `ephem.get_ayanamsa_name(sid_mode)` — restituisce il nome leggibile

5. **Calcoli pratici in siderale** (~1 pag)
   - Due approcci: usare `SEFLG_SIDEREAL` oppure sottrarre manualmente l'ayanamsha
   - `get_ayanamsa_ex(jd_tt, flag)` — versione estesa con flag
   - `get_ayanamsa_ex_ut(jd_ut, flag)` — versione UT
   - Esempio completo: carta natale vedica con Lahiri — tutti i pianeti in siderale
   - Esempio: calcolare il Nakshatra (divisione in 27 settori di 13.33 gradi) del pianeta

**Riferimenti codebase**: `planets.py:swe_set_sid_mode()`, `planets.py:_calc_ayanamsa()`, costanti `SE_SIDM_*`

---

### Capitolo 12 — Corpi minori: asteroidi, centauri, TNO (~8 pagine)

**Obiettivo**: Capire come la libreria gestisce migliaia di corpi oltre ai pianeti principali.

**Sezioni**:

1. **La gerarchia dei corpi celesti** (~1.5 pag)
   - Pianeti (Mercurio-Nettuno): sempre disponibili con precisione massima (JPL DE440)
   - Plutone e Chirone: inclusi nelle efemeridi JPL come i pianeti
   - Asteroidi maggiori (Cerere, Pallade, Giunone, Vesta...): disponibili via kernel SPK dedicati
   - Centauri (Pholus, Nessus, Chariklo...): orbite tra Giove e Nettuno
   - Trans-Nettuniani (Eris, Haumea, Makemake, Sedna...): oltre Nettuno
   - NEA (Near-Earth Asteroids): Apophis, Bennu, Ryugu... — orbite vicine alla Terra
   - ID dei corpi: `SE_CHIRON=15`, `SE_PHOLUS=16`, `SE_CERES=17`, ... `SE_AST_OFFSET + numero`

2. **La catena di calcolo** (~2 pag)
   - La libreria prova 4 metodi in ordine, dal piu preciso al meno preciso:
   - **1. Kernel SPK**: file binario JPL con posizioni precise — sub-arcsecondo
   - **2. Download SPK automatico**: se il kernel non c'e, lo scarica da JPL Horizons
   - **3. ASSIST/REBOUND**: integrazione numerica n-corpi — preciso ma serve il pacchetto opzionale
   - **4. Fallback kepleriano**: leggi di Keplero + perturbazioni — funziona sempre ma meno preciso
   - Codice: `ephem.set_auto_spk_download(True)` — abilita il download automatico
   - Codice: `ephem.calc_ut(jd, ephem.SE_CERES)` — la libreria sceglie il metodo migliore
   - Esempio: calcolare la posizione di Cerere — confronto tra i metodi

3. **Kernel SPK: il gold standard** (~1.5 pag)
   - Cosa sono: file binari NASA con polinomi di Chebyshev per traiettorie precise
   - `ephem.download_and_register_spk(body_id, jd_start, jd_end)` — scarica e registra
   - `ephem.register_spk_body(body_id, spk_path, naif_id)` — registra un SPK esistente
   - `ephem.list_spk_bodies()` — quali corpi hanno un SPK caricato
   - Asteroidi "maggiori" con SPK preconfigurato: `ephem.list_major_asteroids()`
   - `ephem.ensure_major_asteroid_spk(body_id)` — scarica se necessario
   - Esempio: scaricare l'SPK di Apophis e calcolare la sua posizione nel 2029

4. **Cercare un asteroide per nome o numero** (~1 pag)
   - `ephem.get_asteroid_number("Vesta")` — restituisce 4
   - `ephem.calc_asteroid_by_number(jd, 4)` — calcola Vesta per numero
   - `ephem.fetch_orbital_elements_from_sbdb(body_id)` — scarica gli elementi orbitali da JPL SBDB
   - La cache: `ephem.clear_asteroid_elements_cache()`, `ephem.clear_asteroid_name_cache()`
   - Esempio: "dov'e Eros adesso?" — cercare per nome e calcolare la posizione

5. **Il fallback kepleriano: come funziona** (~2 pag)
   - 37 corpi con elementi orbitali pre-caricati: 13 cintura principale, 6 centauri, 9 TNO, 9 NEA
   - Elementi osculanti: i 6 numeri che descrivono un'ellisse (a, e, i, omega, Omega, M)
   - Equazione di Keplero: dato il tempo, trova la posizione sull'ellisse
   - Perturbazioni secolari: Giove e Saturno "spingono" l'orbita, facendola ruotare
   - Perturbazioni di Laplace-Lagrange: eccentricita e inclinazione oscillano
   - Multi-epoca: per 6 corpi, elementi a intervalli di 50 anni (1650-2450) — si sceglie l'epoca piu vicina
   - Modello di librazione: per i plutini (Ixion, Orcus) in risonanza 2:3 con Nettuno
   - Precisione tipica: ~7" a 1 mese, ~2' a 1 anno, ~3.6 gradi a 50 anni (dal singolo epoch)
   - Esempio: confronto SPK vs Keplerian per Cerere a +10 anni — quanto sbaglia?

**Riferimenti codebase**: `minor_bodies.py` (tutto), `spk.py`, `spk_auto.py`, `rebound_integration.py`

---

### Capitolo 13 — Pianeti ipotetici e parti arabe (~6 pagine)

**Obiettivo**: Capire i corpi celesti non fisici usati in varie tradizioni astrologiche.

**Sezioni**:

1. **I pianeti uraniani (Scuola di Amburgo)** (~2 pag)
   - La Scuola di Amburgo (Alfred Witte, 1920s): ha postulato 8 "pianeti" trans-nettuniani
   - Non esistono fisicamente — sono punti matematici con orbite ipotetiche
   - Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus, Poseidon
   - Ogni pianeta ha elementi orbitali definiti (semi-asse, eccentricita, inclinazione, etc.)
   - Usati nell'astrologia uraniana e nella tecnica dei "punti medi"
   - Codice: `ephem.calc_cupido(jd_tt)`, `ephem.calc_hades(jd_tt)`, etc.
   - Codice generico: `ephem.calc_uranian_planet(jd_tt, planet_id)`
   - Anche accessibili via `calc_ut(jd, SE_CUPIDO)`, `calc_ut(jd, SE_HADES)`, etc.
   - Esempio vita reale: calcolare la posizione di Kronos e confrontarla con il MC

2. **Altri corpi ipotetici** (~1.5 pag)
   - **Transpluto (Isis)**: ipotetico pianeta oltre Plutone (prima della scoperta di Eris)
   - **Vulcano**: ipotetico pianeta tra Mercurio e il Sole (cercato nel 1800, mai trovato)
   - **Luna di Waldemath**: ipotetico secondo satellite della Terra (1898, mai confermato)
   - **Pianeta X di Pickering**: previsione del 1919 per un pianeta oltre Nettuno
   - **Proserpina**: altro ipotetico trans-plutoniano
   - **Luna Bianca (Selena)**: il punto opposto alla Lilith Nera — non un corpo fisico
   - Codice: `ephem.calc_transpluto(jd_tt)`, `ephem.calc_vulcan(jd_tt)`, etc.
   - Codice: `ephem.calc_white_moon_position(jd_tt)` — Selena

3. **Corpi ipotetici personalizzati** (~1 pag)
   - File di orbite fittizie: formato testo con elementi orbitali
   - `ephem.parse_orbital_elements("mio_file.txt")` — carica orbite custom
   - `ephem.get_orbital_body_by_name(elements, "NomePianeta")` — cerca per nome
   - `ephem.calc_orbital_position(jd_tt, elements)` — calcola la posizione
   - Il file bundled: `ephem.load_bundled_fictitious_orbits()` — orbite predefinite
   - Polinomi T: `TPolynomial` — per elementi orbitali che variano con il tempo
   - Esempio: definire un corpo ipotetico custom e calcolarne la posizione

4. **Le parti arabe (lotti)** (~1.5 pag)
   - Antichissima tradizione (Persia, poi Medioevo): punti calcolati da combinazioni di pianeti
   - Formula base: Parte = ASC + Pianeta_A - Pianeta_B
   - **Parte di Fortuna**: ASC + Luna - Sole (di giorno) / ASC + Sole - Luna (di notte)
   - **Parte di Spirito**: l'inverso della Parte di Fortuna
   - **Parte dell'Amore**: ASC + Venere - Sole
   - **Parte della Fede**: ASC + Mercurio - Luna
   - Giorno vs notte: alcune parti invertono la formula di notte
   - `ephem.calc_all_arabic_parts(jd, lat, lon)` — calcola tutte le parti
   - `ephem.is_day_chart(jd, lat, lon)` — il Sole e sopra l'orizzonte?
   - Esempio vita reale: la Parte di Fortuna nel tuo tema natale — calcolala

**Riferimenti codebase**: `hypothetical.py` (tutto), `arabic_parts.py`, costanti `SE_CUPIDO`...`SE_POSEIDON`

---

### Capitolo 14 — Precisione, performance e modalita LEB (~8 pagine)

**Obiettivo**: Capire i compromessi tra precisione e velocita, e come ottimizzare per il proprio uso.

**Sezioni**:

1. **Quanto e precisa questa libreria?** (~2 pag)
   - Pianeti principali: sub-arcsecondo rispetto a JPL DE440 (la sorgente e la stessa)
   - Luna: sub-arcsecondo — la piu complessa per via delle perturbazioni
   - Stelle fisse: ~0.5" — limitato dalla precisione del catalogo Hipparcos
   - Case: ~0.02" per tutti i sistemi — praticamente identico a qualsiasi altra implementazione
   - Ayanamsha: ~0.07" — differenze tra implementazioni diverse
   - Confronto con altri software: le differenze sono dovute a scelte di modello, non a errori
   - La nostra filosofia: non copiamo un altro software, usiamo le fonti astronomiche migliori
   - Dove siamo PIU precisi: magnitudini planetarie (Mallama & Hilton 2018), diametri apparenti (IAU 2015)
   - Delta-T: dove la precisione cala — lontano dal presente, Delta-T e meno certo

2. **Le scale temporali di errore** (~1.5 pag)
   - Le efemeridi JPL sono valide per periodi limitati (DE440: 1550-2650)
   - Fuori range: la libreria usa DE441 (-13200 a +17191) con precisione gradualmente inferiore
   - Delta-T: prima del 1600, l'incertezza su Delta-T domina l'errore sulla Luna
   - Stelle fisse: il moto proprio e preciso per ~pochi secoli dal catalogo J2000
   - Corpi minori (Keplerian): si degradano rapidamente (arcminuti in mesi)
   - Regola pratica: per astrologia moderna (1900-2100), la precisione e sub-arcsecondo per tutto

3. **Modalita LEB: efemeridi binarie precompilate** (~2 pag)
   - Il problema: Skyfield e preciso ma lento — ~1ms per calcolo
   - LEB (LibEphemeris Binary): approssimazioni con polinomi di Chebyshev precompilati
   - Speedup: ~14x piu veloce di Skyfield (~70 microsec per calcolo)
   - Precisione LEB: sub-arcsecondo rispetto a Skyfield — nessuna differenza pratica
   - Come funziona: il file .leb contiene coefficienti per ogni corpo, ogni intervallo temporale
   - Setup: `ephem.set_leb_file("/path/to/ephemeris.leb")`
   - Oppure: variabile ambiente `LIBEPHEMERIS_LEB=/path/to/ephemeris.leb`
   - Generazione: `poe leb:generate:medium:groups` — crea il file .leb per il tier medio
   - Tre tier: base (~piccolo), medium (~medio), extended (~grande)
   - Fallback automatico: se il corpo/flag non e supportato da LEB, usa Skyfield
   - Codice: `ephem.set_calc_mode("auto")` — usa LEB se disponibile, altrimenti Skyfield

4. **Modalita di calcolo e configurazione** (~1.5 pag)
   - `set_calc_mode("auto")`: default — LEB se disponibile, altrimenti Skyfield
   - `set_calc_mode("skyfield")`: forza sempre Skyfield (utile per debug/confronti)
   - `set_calc_mode("leb")`: forza LEB — errore se non disponibile
   - `set_precision_tier("medium")`: sceglie quale set di efemeridi JPL usare
   - `set_strict_precision(True)`: modalita di precisione massima
   - `EphemerisContext`: contesto thread-safe per calcoli paralleli
   - Codice: `with EphemerisContext(sid_mode=SE_SIDM_LAHIRI) as ctx: ...`

5. **Download e gestione dati** (~1 pag)
   - `ephem.download_for_tier("medium")` — scarica i file JPL necessari
   - `ephem.download_leb_for_tier("medium")` — scarica il file LEB precompilato
   - `ephem.download_iers_finals()` — scarica dati IERS per Delta-T preciso
   - `ephem.discover_local_spks()` — trova SPK gia scaricati nella cache
   - `ephem.ensure_all_ephemerides()` — scarica tutto il necessario
   - Cache directory: `ephem.set_spk_cache_dir("/path")`, `ephem.set_iers_cache_dir("/path")`

**Riferimenti codebase**: `fast_calc.py`, `leb_reader.py`, `leb_format.py`, `state.py`, `context.py`, `download.py`

---

### Capitolo 15 — Ricettario: calcoli pratici dalla A alla Z (~15 pagine)

**Obiettivo**: Ricette complete copia-incolla per i calcoli piu comuni. Ogni ricetta: problema, codice, spiegazione del risultato.

**Ricette**:

1. **Tema natale completo** (~2 pag)
   - Input: data, ora, luogo di nascita
   - Calcolare: JD, tutti i pianeti, case, angoli, nodi lunari, Lilith, Chirone
   - Formattare: segno, gradi, minuti per ogni posizione
   - Codice completo funzionante (30-40 righe)
   - Uso di `split_deg()` per formattazione zodiacale

2. **Transiti del giorno** (~1.5 pag)
   - Input: data e posizioni natali
   - Trovare: quali pianeti di transito sono in aspetto (congiunzione, opposizione, quadratura, trigono, sestile) con i pianeti natali
   - Calcolo dell'orbe (distanza dall'aspetto esatto)
   - Codice: ciclo su tutti i pianeti, `difdeg2n()` per differenze angolari

3. **Attraversamenti (ingressi nei segni)** (~1 pag)
   - Quando il Sole entra in Ariete (equinozio di primavera)?
   - `solcross_ut(0.0, jd_start)` — trova il prossimo attraversamento di 0 gradi
   - Quando la Luna entra in un segno qualsiasi?
   - `mooncross_ut(target_deg, jd_start)` — trova l'attraversamento
   - Esempio: tutti gli ingressi del Sole nei 12 segni per l'anno corrente

4. **Retrogradazione di Mercurio** (~1 pag)
   - Trovare inizio (stazione retrograda) e fine (stazione diretta)
   - `swe_find_station_ut(SE_MERCURY, jd_start)` — prossima stazione
   - `is_retrograde(SE_MERCURY, jd)` — e retrogrado in questo momento?
   - Codice: tutte le retrogradazioni di Mercurio dell'anno

5. **Luna Nuova e Luna Piena** (~1 pag)
   - Luna Nuova: quando `longitudine_Luna == longitudine_Sole` (congiunzione)
   - Luna Piena: quando la differenza e 180 gradi (opposizione)
   - Usare `mooncross_ut()` con la longitudine del Sole come target
   - Codice: tutte le Lune Nuove e Piene dell'anno

6. **Prossima eclissi visibile dalla mia citta** (~1.5 pag)
   - `sol_eclipse_when_loc(jd, flag, (lon, lat, alt))` in un ciclo
   - Filtro: solo eclissi con magnitudine > 0.5
   - Output: data, tipo, magnitudine, orari dei contatti
   - Codice completo per Roma, Milano, o qualsiasi citta

7. **Carta del cielo siderale (vedica)** (~1 pag)
   - Come la ricetta 1, ma con `set_sid_mode(SE_SIDM_LAHIRI)` e `SEFLG_SIDEREAL`
   - Calcolo dei Nakshatra (27 settori)
   - Calcolo del Dasha (cenno — serve solo la Luna siderale)

8. **Visibilita dei pianeti stasera** (~1 pag)
   - Per ogni pianeta: calcolare altezza al tramonto del Sole e magnitudine
   - Filtrare: altezza > 10 gradi e magnitudine < 6 (visibile a occhio nudo)
   - Elongazione dal Sole: > 15 gradi per non essere perso nel bagliore
   - Codice con output tipo: "Giove: visibile, alt 45, mag -2.3, a sud-ovest"

9. **Distanza tra due pianeti (aspetti)** (~1 pag)
   - Calcolare la distanza angolare tra due pianeti qualsiasi
   - Determinare l'aspetto piu vicino (0, 60, 90, 120, 180)
   - Calcolare l'orbe (quanto manca all'aspetto esatto)
   - Se l'aspetto e applicante (si avvicina) o separante (si allontana)

10. **Rivoluzione solare (Solar Return)** (~1 pag)
    - Trovare il momento esatto in cui il Sole torna alla stessa longitudine di nascita
    - `solcross_ut(sun_natal_lon, jd_birthday)` — trova il ritorno solare
    - Calcolare la carta per quel momento esatto
    - Codice: rivoluzione solare per l'anno corrente

11. **Ore planetarie** (~1 pag)
    - Dividere il giorno (alba-tramonto) e la notte (tramonto-alba successiva) in 12 parti
    - Ogni "ora" e governata da un pianeta nella sequenza caldea
    - Codice: calcolare alba e tramonto, dividere, assegnare i pianeti

12. **Efemeride mensile stampabile** (~1 pag)
    - Generare una tabella con la posizione di tutti i pianeti per ogni giorno del mese
    - Formato: data, Sole, Luna, Mercurio... Plutone, Nodo, Lilith
    - Output CSV o tabella formattata
    - Codice completo con formattazione zodiacale

---

## Note di implementazione

### Convenzioni per tutto il manuale

- Ogni capitolo inizia con un paragrafo "Cosa imparerai" e finisce con un "Riepilogo"
- Gli esempi di codice sono sempre completi e funzionanti (copia-incolla)
- `import libephemeris as ephem` e l'unico import usato ovunque
- Le date degli esempi usano eventi astronomici reali e riconoscibili (eclissi 2024, etc.)
- I concetti tecnici hanno sempre un'analogia con la vita quotidiana
- I numeri hanno sempre le unita: gradi, arcminuti ('), arcsecondi ("), UA, km

### Struttura dei file

```
docs/manual/
  00-prefazione.md
  01-cielo-dalla-terra.md
  02-misurare-il-tempo.md
  03-coordinate-celesti.md
  04-efemeridi.md
  05-posizione-pianeti.md
  06-luna.md
  07-case-astrologiche.md
  08-stelle-fisse.md
  09-eclissi.md
  10-alba-tramonto-visibilita.md
  11-zodiaco-siderale.md
  12-corpi-minori.md
  13-pianeti-ipotetici-parti-arabe.md
  14-precisione-performance.md
  15-ricettario.md
```

### Priorita di scrittura consigliata

1. **Capitoli 1-3** (fondamenta): senza questi nulla ha senso
2. **Capitolo 5** (calc_ut): la funzione piu usata
3. **Capitolo 7** (case): fondamentale per astrologia
4. **Capitolo 15** (ricettario): il piu utile per uso pratico immediato
5. **Capitoli 2, 4, 6** (tempo, efemeridi, Luna): completano le basi
6. **Capitoli 9, 10, 11** (eclissi, visibilita, siderale): argomenti avanzati
7. **Capitoli 8, 12, 13** (stelle, asteroidi, ipotetici): specialistici
8. **Capitolo 14** (performance): per utenti avanzati
9. **Capitolo 0** (prefazione): scrivere per ultimo, quando sai cosa c'e dentro
