# LibEphemeris — Piano di Verifica Esaustiva

## Contesto Generale del Progetto

LibEphemeris e' una reimplementazione clean-room in puro Python della Swiss Ephemeris,
la libreria standard de facto per calcoli astronomici ed astrologici. Usa esclusivamente
dati NASA JPL (DE440/DE441) tramite Skyfield, anziche' i file di effemeridi proprietari
della Swiss Ephemeris. L'obiettivo e' la compatibilita' 1:1 con l'API PySwissEph
mantenendo una precisione sub-arcsecond rispetto ai dati JPL.

### Architettura a 3 Backend

La libreria supporta 3 backend di calcolo con fallback automatico:

1. **LEB** (LibEphemeris Binary) — Polinomi Chebyshev precomputed dai dati JPL.
   ~5 microsecondi per valutazione. File `.leb` (LEB1) o `.leb` compresso (LEB2).
   Il file `base_core.leb` (8.7 MB, 14 corpi) e' incluso nel wheel PyPI.

2. **Horizons** — REST API di NASA JPL Horizons. Zero file locali necessari.
   ~300ms per il primo calcolo, poi cached. Richiede internet.

3. **Skyfield** — Calcolo diretto dai file DE440/DE441 via Skyfield.
   ~120 microsecondi per valutazione. Il gold standard per precisione.

Modalita' di calcolo (`set_calc_mode()` o `LIBEPHEMERIS_MODE`):
- `"auto"` (default): LEB -> Horizons (se no DE440) -> Skyfield
- `"leb"`: solo LEB, errore se non configurato
- `"horizons"`: solo Horizons API
- `"skyfield"`: solo Skyfield/DE440

### Corpi Celesti Supportati

| ID | Corpo | Note |
|----|-------|------|
| 0 | Sole | Geocentrico (apparente), eliocentrico = (0,0,0) |
| 1 | Luna | Piu' alta precisione richiesta (moto rapido ~13 deg/day) |
| 2-9 | Mercurio-Plutone | Include correzione baricentro->centro fisico per 5-9 |
| 10 | Nodo Medio | Moto retrogrado liscio ~-0.053 deg/day |
| 11 | Nodo Vero | Oscillazione ±1.5 deg attorno al medio, orbita osculante |
| 12 | Apogeo Medio | Black Moon Lilith, ~40 deg/anno diretto |
| 13 | Apogeo Osculante | True Lilith, oscillazione rapida |
| 14 | Terra | Geocentrico = (0,0,0), utile per eliocentrico |
| 15 | Chirone | Centauro, orbita caotica, richiede SPK per alta precisione |
| 17-20 | Cerere-Vesta | Asteroidi cintura principale |
| 21-22 | Apogeo/Perigeo Interpolati | Apsidi lunari da passaggi fisici distanza |
| 40-47 | Uraniani | Cupido-Poseidon, Keplerian helio -> geocentrico |
| 48 | Transpluto/Isis | Ipotetico, Keplerian |
| 10000+ | Asteroidi via offset | SE_AST_OFFSET + numero |

### Tier di Effemeridi

| Tier | File JPL | Range date | Dimensione |
|------|----------|-----------|------------|
| base | de440s.bsp | 1849-2150 | ~31 MB |
| medium | de440.bsp | 1550-2650 | ~128 MB (default) |
| extended | de441.bsp | -13200 a +17191 | ~3.1 GB |

### Flag di Calcolo

I flag controllano il tipo di calcolo e il frame di riferimento. Sono bitmask
combinabili con `|`.

| Flag | Valore | Effetto |
|------|--------|---------|
| SEFLG_SPEED | 256 | Calcola velocita' (quasi sempre necessario) |
| SEFLG_HELCTR | 8 | Osservatore al Sole (eliocentrico) |
| SEFLG_BARYCTR | 16384 | Osservatore al baricentro SSB |
| SEFLG_TOPOCTR | 32768 | Osservatore sulla superficie terrestre |
| SEFLG_EQUATORIAL | 2048 | Output in RA/Dec anziche' lon/lat eclittica |
| SEFLG_J2000 | 32 | Frame J2000 (no precessione a data) |
| SEFLG_NONUT | 64 | Eclittica/equatore medio (no nutazione) |
| SEFLG_SIDEREAL | 65536 | Zodiaco siderale (richiede set_sid_mode) |
| SEFLG_TRUEPOS | 16 | Posizione geometrica (no light-time/aberrazione) |
| SEFLG_NOABERR | 1024 | No aberrazione (astrometrica) |
| SEFLG_NOGDEFL | 512 | No deflessione gravitazionale |
| SEFLG_XYZ | 4096 | Output cartesiano (x,y,z in AU) |
| SEFLG_RADIANS | 8192 | Angoli in radianti |
| SEFLG_ICRS | 131072 | Frame ICRS |

### Sistemi di Case (24)

P=Placidus, K=Koch, O=Porphyrius, R=Regiomontanus, C=Campanus,
E/A=Uguale, W=Segno Intero, X=Meridiano, M=Morinus, H=Orizzontale,
T=Topocentric, B=Alcabitus, G=Gauquelin, I/i=Sunshine, L=Pullen SD,
N=Pullen SR, Q=Carter, Y=APC, F=Carter Poli-Equatorial, U=Krusinski,
D=Sripati, J=vedic

### Ayanamsha (43 + custom)

Modi 0-42 predefiniti (Fagan-Bradley, Lahiri, Raman, etc.)
Modo 255 = custom con t0 e ayan_t0 utente.
Modi stellari (27-29) usano posizioni reali di stelle.

### Precisione Attesa

| Categoria | Tipica | Max | Note |
|-----------|--------|-----|------|
| Pianeti (Sole-Plutone) | 0.04-0.26" | 0.75" | Sub-arcsecond |
| Luna | 0.70" | 3.32" | Modelli lunari diversi |
| Cuspidi case | < 0.01" | 0.02" | 24 sistemi testati |
| Stelle fisse | < 0.1" | 0.51" | 116 stelle Hipparcos |
| Eclissi solari | — | < 6s | Timing |
| Eclissi lunari | — | < 8s | Timing |
| LEB vs Skyfield | < 0.005" | < 0.01" | Compressione Chebyshev |
| Horizons vs Skyfield | < 0.001" | 0.03" | Eliocentrico offset |

### File Chiave dell'Implementazione

| File | Righe | Scopo |
|------|-------|-------|
| `planets.py` | ~5500 | Core: `_calc_body()`, flag dispatch, velocity, planet centers |
| `houses.py` | ~5000 | 24 sistemi di case, cusps, angles (ASC, MC, Vertex) |
| `eclipse.py` | ~14000 | Eclissi solari/lunari, occultazioni, Besselian elements |
| `lunar.py` | ~1900 | Nodi, apsidi, True Node osculante, perturbazioni ELP2000 |
| `heliacal.py` | ~2700 | Visibilita' eliacale, magnitudine limite |
| `fast_calc.py` | ~1100 | LEB fast path pipeline (Clenshaw, correzioni, frame) |
| `horizons_backend.py` | ~720 | HTTP client, pipeline geocentrica, cache |
| `state.py` | ~2000 | Stato globale, setter/getter, lock, init lazy |
| `context.py` | ~600 | EphemerisContext thread-safe |
| `time_utils.py` | ~800 | Julian day, Delta T, sidereal time, UTC |
| `fixed_stars.py` | ~4800 | Catalogo Hipparcos, moto proprio |
| `crossing.py` | ~1200 | Equinozi, solstizi, crossings, stazioni |
| `constants.py` | ~1430 | 785 costanti + 3 funzioni |

---

## Piano di Verifica

Ogni sezione descrive verifiche standalone. Ogni check e' un'asserzione Python
individuale che puo' essere eseguita senza dipendenze da test suite esistenti.
Il riferimento e' pyswisseph 2.10.

### Notazione

- `jd` = Julian Day (es. 2451545.0 = J2000.0 = 1 Gennaio 2000 12:00 TT)
- `body` = ID corpo celeste (SE_SUN=0, SE_MOON=1, etc.)
- `flags` = bitmask flag di calcolo
- `"<"` = strettamente minore (tolleranza)
- Quando si dice "per N date casuali" si intende con seed fisso per riproducibilita'

---

## 1. Posizioni Planetarie — Accuratezza Cross-Backend

Questa e' la sezione piu' critica. Verifica che tutti e 3 i backend producano
posizioni consistenti tra loro e con la reference implementation.

### 1.1 Skyfield vs pyswisseph — tutti i 22 corpi

Corpi: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22
Date: 500 JD casuali (seed=42) in [2415020.5, 2488069.5] (1900-2100)
Flag: SEFLG_SPEED

Per ogni (corpo, data):
- [ ] Longitudine: |lib - ref| < 1.0 arcsec (pianeti), < 4.0 arcsec (Luna)
- [ ] Latitudine: |lib - ref| < 1.0 arcsec
- [ ] Distanza: |lib - ref| < 1e-5 AU
- [ ] Velocita' lon: |lib - ref| < 0.01 deg/day
- [ ] Nessun crash o eccezione inattesa

Stima: 22 corpi x 500 date x 5 check = **55.000 check**

### 1.2 Skyfield vs pyswisseph — flag varianti

Per ogni flag in [SEFLG_SPEED, SPEED|SIDEREAL, SPEED|EQUATORIAL, SPEED|J2000,
SPEED|NOABERR, SPEED|HELCTR, SPEED|BARYCTR, SPEED|TRUEPOS, SPEED|NONUT,
SPEED|NOGDEFL]:
Per 10 corpi (0,1,2,4,5,6,8,9,14,15):
Per 100 date casuali:
- [ ] Longitudine: |lib - ref| < 2.0 arcsec
- [ ] Latitudine: |lib - ref| < 2.0 arcsec

Stima: 10 flag x 10 corpi x 100 date x 2 check = **20.000 check**

### 1.3 LEB vs Skyfield — 14 corpi core

Corpi: 0,1,2,3,4,5,6,7,8,9,10,11,12,14 (quelli nel base_core.leb)
Date: 500 JD in [2396760, 2506330] (range base tier)
Metodo: `set_calc_mode("leb")` vs `set_calc_mode("skyfield")`

Per ogni (corpo, data):
- [ ] |lon_leb - lon_sky| < 0.005 arcsec
- [ ] |lat_leb - lat_sky| < 0.005 arcsec
- [ ] |dist_leb - dist_sky| < 1e-7 AU
- [ ] |speed_leb - speed_sky| < 0.001 deg/day

Stima: 14 x 500 x 4 = **28.000 check**

### 1.4 LEB vs Skyfield — con flag siderale/equatoriale/J2000

Per 14 corpi x 200 date x 4 flag (sidereal, equatorial, J2000, noaberr):
- [ ] |lon_diff| < 0.005 arcsec

Stima: 14 x 200 x 4 = **11.200 check**

### 1.5 Horizons vs Skyfield — 13 corpi geocentrici

Corpi: 0,1,2,3,4,5,6,7,8,9,14,15,17
Date: 200 JD in [2430000, 2470000]
6 flag: default, sidereal, equatorial, J2000, noaberr, truepos

Per ogni (corpo, data, flag):
- [ ] |lon_hz - lon_sky| < 0.001 arcsec (geocentrico)
- [ ] |lat_hz - lat_sky| < 0.001 arcsec

Stima: 13 x 200 x 6 x 2 = **31.200 check**

### 1.6 Horizons vs Skyfield — eliocentrici

Per 10 corpi x 100 date:
- [ ] |lon_hz - lon_sky| < 0.03 arcsec (offset sistematico noto)

Stima: 10 x 100 = **1.000 check**

### 1.7 LEB2 vs LEB1 — integrita' compressione

Per 14 corpi x 500 date:
- [ ] |pos_leb2 - pos_leb1| < 0.001 arcsec

Stima: 14 x 500 = **7.000 check**

### 1.8 Triple cross-validation

Per 11 corpi (comuni a tutti i backend) x 100 date:
- [ ] |Skyfield - LEB| < 0.01 arcsec
- [ ] |Skyfield - Horizons| < 0.01 arcsec
- [ ] |LEB - Horizons| < 0.01 arcsec
(per lon e lat = 6 check per combinazione)

Stima: 11 x 100 x 6 = **6.600 check**

---

## 2. Flag Combinazioni — Robustezza

Verifica che OGNI combinazione di flag con OGNI corpo non crashi e produca
output valido.

### 2.1 Singoli flag — no crash, output valido

13 flag x 22 corpi x 3 backend (skyfield, leb, horizons) x 50 date

Per ogni (flag, corpo, backend, data):
- [ ] Nessun crash
- [ ] Tupla a 6 elementi
- [ ] Tutti i valori math.isfinite()
- [ ] lon/RA nel range atteso
- [ ] dist >= 0 (escluso XYZ)
- [ ] speed finita

Stima: 13 x 22 x 3 x 50 x 6 = **257.400 check**

### 2.2 Coppie di flag compatibili

78 coppie x 5 corpi x 20 date:
- [ ] No crash
- [ ] Output valido (6 valori finiti)

Stima: 78 x 5 x 20 x 2 = **15.600 check**

### 2.3 SEFLG_TOPOCTR con diverse localita'

5 localita' x 10 corpi x 20 date:
- [ ] set_topo() + calc_ut() con SEFLG_TOPOCTR: no crash
- [ ] Risultato differisce dal geocentrico (almeno per la Luna)

Stima: 5 x 10 x 20 x 2 = **2.000 check**

### 2.4 SEFLG_XYZ coerenza con sferico

10 corpi x 50 date:
- [ ] Conversione manuale lon/lat/dist -> xyz matches SEFLG_XYZ output (< 1e-8)

Stima: 10 x 50 x 3 = **1.500 check**

### 2.5 SEFLG_RADIANS coerenza con gradi

10 corpi x 50 date:
- [ ] math.radians(lon_deg) == lon_rad (< 1e-10)

Stima: 10 x 50 x 2 = **1.000 check**

### 2.6 Sun heliocentric invariante

50 date:
- [ ] calc_ut(jd, SE_SUN, SEFLG_SPEED|SEFLG_HELCTR) == (0,0,0,0,0,0)

Stima: 50 x 6 = **300 check**

---

## 3. Sistemi di Case (24 sistemi)

Verifica tutti i 24 sistemi di case per correttezza geometrica e
compatibilita' con la reference implementation.

### 3.1 Validita' output per ogni sistema

24 sistemi x 6 localita' x 20 date:
Localita': (0,0), (12.5,41.9), (139.7,35.7), (-73.9,40.7), (0,60), (0,-60)

Per ogni (sistema, localita', data):
- [ ] len(cusps) >= 12
- [ ] Tutti i cusps in [0, 360)
- [ ] ASC (ascmc[0]) in [0, 360)
- [ ] MC (ascmc[1]) in [0, 360)
- [ ] ASC != MC

Stima: 24 x 6 x 20 x 5 = **14.400 check**

### 3.2 Confronto vs pyswisseph

24 sistemi x 3 localita' x 20 date:
- [ ] |cusp_i_lib - cusp_i_ref| < 0.1 deg per tutti i cusps (i=1..12)
- [ ] |ASC_lib - ASC_ref| < 0.01 deg
- [ ] |MC_lib - MC_ref| < 0.01 deg

Stima: 24 x 3 x 20 x 14 = **20.160 check**

### 3.3 houses_ex con flag siderale

5 sistemi x 3 ayanamsha x 20 date:
- [ ] Cuspidi siderali != tropicali
- [ ] Differenza ~ ayanamsha

Stima: 5 x 3 x 20 x 2 = **600 check**

### 3.4 Latitudini polari

Per Placidus e Koch a lat 70, 80, 85, 89:
- [ ] No crash (fallback a Porphyry)
Per Equal e Whole Sign a lat 90:
- [ ] Funziona normalmente

Stima: 2 x 4 x 10 + 2 x 10 = **100 check**

### 3.5 houses_armc

5 sistemi x 10 ARMC x 3 eps:
- [ ] Output valido
- [ ] Coerente con houses() per stesso momento

Stima: 5 x 10 x 3 x 5 = **750 check**

### 3.6 house_pos

10 lon planetarie x 5 sistemi x 5 date:
- [ ] Risultato in [1.0, 13.0)

Stima: 10 x 5 x 5 = **250 check**

### 3.7 Gauquelin sectors

5 corpi x 10 date:
- [ ] Settore in [1, 36]

Stima: 5 x 10 = **50 check**

---

## 4. Ayanamsha (43 modi + custom)

### 4.1 Tutti i modi producono valori distinti

A J2000:
- [ ] get_ayanamsa_ut(jd) per mode 0-42: 43 valori tutti diversi
- [ ] Tutti in range [-5, 30] deg

Stima: **43 + 43 = 86 check**

### 4.2 Coerenza UT/TT

10 date x 5 modi:
- [ ] |get_ayanamsa_ut(jd_ut) - get_ayanamsa(jd_tt)| ragionevole

Stima: **50 check**

### 4.3 get_ayanamsa_ex/get_ayanamsa_ex_ut

5 modi x 5 date:
- [ ] Ritorna (ayanamsa, retflag)
- [ ] ayanamsa matches get_ayanamsa_ut

Stima: **50 check**

### 4.4 get_ayanamsa_name

- [ ] 43 nomi non vuoti

Stima: **43 check**

### 4.5 Custom ayanamsha (mode 255)

- [ ] set_sid_mode(255, t0=J2000, ayan_t0=23.5): get_ayanamsa_ut(J2000) == 23.5

Stima: **5 check**

### 4.6 Star-based (modes 27, 28, 29)

3 modi x 5 date:
- [ ] No crash, risultato finito

Stima: **15 check**

### 4.7 Posizioni siderali per tutti i 43 modi

43 modi x 3 corpi x 10 date:
- [ ] lon_siderale in [0, 360)
- [ ] lon_siderale != lon_tropicale

Stima: 43 x 3 x 10 x 2 = **2.580 check**

---

## 5. Eclissi Solari

### 5.1 Ricerca globale

26 anni (2000-2025):
- [ ] sol_eclipse_when_glob() trova almeno 1 eclissi per anno
- [ ] JD eclissi nell'anno
- [ ] Tipo valido (TOTAL|ANNULAR|PARTIAL|HYBRID)

Stima: 26 x 3 = **78 check**

### 5.2 Ricerca locale

5 eclissi note x 5 localita':
- [ ] sol_eclipse_when_loc() trova eclissi
- [ ] Timing entro 60 secondi dal dato noto

Stima: 25 x 2 = **50 check**

### 5.3 Geometria

5 eclissi:
- [ ] sol_eclipse_where(): coordinate centrali valide
- [ ] sol_eclipse_how(): magnitudine > 0
- [ ] sol_eclipse_how_details(): dettagli finiti
- [ ] Obscuration in [0, 1]

Stima: 5 x 4 = **20 check**

### 5.4 Besselian elements

3 eclissi:
- [ ] x, y, d, l1, l2, mu tutti finiti
- [ ] Derivate temporali finite

Stima: 3 x 12 = **36 check**

### 5.5 Saros/Inex

5 eclissi:
- [ ] get_saros_number() > 0
- [ ] get_inex_number() finito

Stima: 5 x 2 = **10 check**

---

## 6. Eclissi Lunari

### 6.1 Ricerca globale

26 anni:
- [ ] lun_eclipse_when() trova almeno 1 per anno
- [ ] Tipo valido

Stima: 26 x 2 = **52 check**

### 6.2 Dettagli

5 eclissi note:
- [ ] lun_eclipse_how(): gamma, magnitudine
- [ ] Umbral/penumbral magnitude finiti

Stima: 5 x 4 = **20 check**

### 6.3 Occultazioni

3 occultazioni:
- [ ] lun_occult_when_glob() trova evento
- [ ] planet_occult_when_glob() funziona

Stima: 3 x 2 = **6 check**

---

## 7. Alba/Tramonto/Transito

### 7.1 Sole

5 localita' x 20 date:
- [ ] rise_trans(CALC_RISE): JD alba valido
- [ ] rise_trans(CALC_SET): JD tramonto valido
- [ ] Alba < Tramonto
- [ ] Durata giorno ragionevole (6-18 ore)
- [ ] rise_trans(CALC_MTRANSIT): mezzogiorno ~12h locale

Stima: 5 x 20 x 5 = **500 check**

### 7.2 Luna

3 localita' x 20 date:
- [ ] Rise e set calcolati
- [ ] Progressione ~ +50 min/giorno

Stima: 3 x 20 x 3 = **180 check**

### 7.3 Pianeti

5 pianeti x 5 date:
- [ ] Rise, set, transit: no crash

Stima: 5 x 5 x 3 = **75 check**

### 7.4 Crepuscoli

3 tipi (civil, nautic, astro) x 5 date:
- [ ] Angoli corretti (-6, -12, -18 deg)

Stima: 3 x 5 x 2 = **30 check**

### 7.5 Latitudini polari

5 date estive/invernali a 70N:
- [ ] Notte bianca estate (no set)
- [ ] Notte polare inverno (no rise)

Stima: **10 check**

---

## 8. Visibilita' Eliacale

### 8.1 heliacal_ut

5 pianeti x 4 tipi di evento x 3 date:
- [ ] JD evento valido

Stima: 5 x 4 x 3 = **60 check**

### 8.2 vis_limit_mag

5 corpi x 3 condizioni atmosferiche:
- [ ] Magnitudine limite finita

Stima: 5 x 3 = **15 check**

---

## 9. Stelle Fisse

### 9.1 Posizioni

20 stelle note (Sirius, Vega, Regulus, Spica, Aldebaran, Antares, Fomalhaut,
Pollux, Arcturus, Capella, Rigel, Procyon, Betelgeuse, Altair, Deneb,
Canopus, Achernar, Hadar, Acrux, Mimosa):
Per 5 date:
- [ ] fixstar_ut(): lon in [0, 360), lat in [-90, 90]
- [ ] dist > 0

Stima: 20 x 5 x 3 = **300 check**

### 9.2 Magnitudini

20 stelle:
- [ ] fixstar_mag(): magnitudine finita
- [ ] Sirius ~ -1.46 (entro 0.2)

Stima: 20 x 2 = **40 check**

### 9.3 Moto proprio

Per Barnard's Star e Sirius (alto moto proprio):
- [ ] Posizione J2000 vs J2100 differisce in modo misurabile

Stima: 2 x 2 = **4 check**

---

## 10. Trasformazioni di Coordinate

### 10.1 cotrans round-trip

3 obliquita' x 200 punti:
- [ ] ecl -> eq -> ecl: |diff| < 1e-10 deg

Stima: 3 x 200 x 2 = **1.200 check**

### 10.2 cotrans_sp

20 punti con velocita':
- [ ] Round-trip preserva posizione e velocita'

Stima: 20 x 4 = **80 check**

### 10.3 azalt

5 localita' x 10 date x 5 corpi:
- [ ] az in [0, 360), alt in [-90, 90]

Stima: 5 x 10 x 5 x 2 = **500 check**

### 10.4 azalt_rev round-trip

20 punti:
- [ ] azalt -> azalt_rev -> azalt: consistente

Stima: 20 x 2 = **40 check**

### 10.5 refrac

- [ ] All'orizzonte: refrazione ~ 34 arcmin
- [ ] round-trip: refrac(refrac(alt, TRUE_TO_APP), APP_TO_TRUE) ~ alt

Stima: **20 check**

---

## 11. Funzioni Tempo

### 11.1 julday/revjul

2000 date casuali (-5000 a +5000):
- [ ] Round-trip esatto

Stima: **2.000 check**

### 11.2 deltat

500 JD in [1900, 2100]:
- [ ] > 0, < 0.01 giorni, smooth

Stima: 500 x 3 = **1.500 check**

### 11.3 sidtime

100 JD:
- [ ] in [0, 24) ore
- [ ] Incremento ~3m56s/giorno

Stima: 100 x 2 = **200 check**

### 11.4 utc_to_jd / inversi

20 date:
- [ ] Round-trip esatto

Stima: 20 x 2 = **40 check**

### 11.5 day_of_week

10 date note:
- [ ] Giorno corretto

Stima: **10 check**

### 11.6 TAI functions

10 date:
- [ ] Round-trip utc->tai->utc, tt->tai->tt

Stima: 10 x 4 = **40 check**

---

## 12. Nodi e Apsidi Lunari

### 12.1 Nodo Medio

200 date:
- [ ] in [0, 360), retrogrado, ~-0.053 deg/day

Stima: 200 x 3 = **600 check**

### 12.2 Nodo Vero

200 date:
- [ ] in [0, 360), oscillazione ±1.5 deg, distanza ~0.0025 AU

Stima: 200 x 3 = **600 check**

### 12.3 Lilith Media

200 date:
- [ ] in [0, 360), diretto, ~40 deg/anno

Stima: 200 x 3 = **600 check**

### 12.4 Apogeo osculante, interpolati

100 date ciascuno per OSCU_APOG, INTP_APOG, INTP_PERG:
- [ ] Posizione valida
- [ ] INTP_APOG e INTP_PERG opposti (~180 deg)

Stima: 300 x 3 = **900 check**

### 12.5 nod_aps_ut

5 corpi x 5 date x 3 flag:
- [ ] 4 tuple (ascending, descending, perihelion, aphelion) valide

Stima: 5 x 5 x 3 x 4 = **300 check**

---

## 13. Corpi Uraniani e Ipotetici

### 13.1 Uraniani geocentrici (40-47)

8 corpi x 20 date:
- [ ] Posizione valida, dist > 10 AU

Stima: 8 x 20 x 2 = **320 check**

### 13.2 Uraniani eliocentrici

8 corpi x 10 date:
- [ ] Matches geocentrico entro ~1 deg

Stima: 8 x 10 = **80 check**

### 13.3 Transpluto

20 date:
- [ ] Geocentrico e eliocentrico funzionano

Stima: 20 x 2 = **40 check**

### 13.4 Siderali

8 corpi x 3 ayanamsha x 5 date:
- [ ] lon siderale != tropicale

Stima: 8 x 3 x 5 = **120 check**

---

## 14. Asteroidi

### 14.1 Principali (5)

5 corpi x 50 date:
- [ ] Posizione valida, confronto ref < 2 arcsec

Stima: 5 x 50 x 2 = **500 check**

### 14.2 Via SE_AST_OFFSET

5 mapping (10001=17, 10002=18, 10003=19, 10004=20, 12060=15):
- [ ] Risultato identico

Stima: 5 x 5 = **25 check**

---

## 15. Lune Planetarie

Per 6 lune (Io, Europa, Ganymede, Callisto, Titan, Triton):
- [ ] is_planetary_moon() == True
- [ ] calc_ut() con SPK registrato: posizione valida

Stima: 6 x 5 = **30 check**

---

## 16. Elementi Orbitali

5 corpi x 5 date:
- [ ] get_orbital_elements_ut(): semi-asse, eccentricita', inclinazione validi
- [ ] orbit_max_min_true_distance(): max > min

Stima: 5 x 5 x 4 + 5 = **105 check**

---

## 17. Fenomeni e Elongazione

10 corpi x 10 date:
- [ ] pheno_ut(): fase, elongazione, magnitudine valide

Stima: 10 x 10 x 3 = **300 check**

Mercury e Venus x 50 date:
- [ ] get_elongation_from_sun() in [0, 180]
- [ ] Mercury max ~28 deg, Venus max ~47 deg

Stima: 2 x 50 x 2 = **200 check**

---

## 18. Crossings e Stazioni

### 18.1 Equinozi/solstizi

6 anni x 4 punti (0, 90, 180, 270 deg):
- [ ] solcross_ut() trova crossing corretto

Stima: 6 x 4 x 2 = **48 check**

### 18.2 Moon crossings

10 crossing targets:
- [ ] mooncross_ut() trova crossing

Stima: **10 check**

### 18.3 Stazioni planetarie

Mars e Mercury x 5 date:
- [ ] find_station_ut(): stazione trovata
- [ ] Velocita' ~ 0 alla stazione

Stima: 2 x 5 x 2 = **20 check**

---

## 19. LEB-Specifico

### 19.1 open_leb factory

- [ ] LEB1 -> LEBReader
- [ ] LEB2 -> LEB2Reader/CompositeLEBReader
- [ ] File inesistente -> FileNotFoundError
- [ ] File non-LEB -> ValueError

Stima: **4 check**

### 19.2 LEBReader API

- [ ] .path, .jd_range, .has_body(), .eval_body(), .eval_nutation(), .delta_t(), .close()
- [ ] Context manager

Stima: **20 check**

### 19.3 LEB2 compression

- [ ] shuffle/unshuffle round-trip
- [ ] compress/decompress round-trip
- [ ] Compressione effettiva (compressed < raw)

Stima: **10 check**

### 19.4 CompositeLEBReader

- [ ] from_file_with_companions
- [ ] has_body dispatch corretto

Stima: **10 check**

### 19.5 fast_calc pipeline

14 corpi x 50 date:
- [ ] fast_calc_ut matches calc_ut skyfield < 0.005 arcsec

Stima: 14 x 50 = **700 check**

### 19.6 Bundled LEB2

- [ ] File esiste nel pacchetto
- [ ] _discover_leb_file() lo trova

Stima: **2 check**

---

## 20. Horizons-Specifico

### 20.1 Client HTTP

- [ ] fetch_state_vector: 6 componenti
- [ ] Cache hit
- [ ] clear_cache, shutdown idempotente

Stima: **10 check**

### 20.2 Pipeline geocentrica

5 corpi x 10 date:
- [ ] Risultato valido con correzioni (light-time, aberrazione, deflessione)

Stima: 5 x 10 = **50 check**

### 20.3 Corpi analitici

2 corpi (Mean Node, Mean Apogee) x 20 date:
- [ ] Nessuna HTTP, risultato matches Skyfield

Stima: 2 x 20 x 2 = **80 check**

### 20.4 Fallback

4 corpi non supportati (11, 13, 21, 22):
- [ ] KeyError -> Skyfield fallback

Stima: 4 x 2 = **8 check**

---

## 21. State Management

### 21.1 Setter/getter round-trip

- [ ] calc_mode: 4 valori + None + invalid
- [ ] topo: set/get
- [ ] sid_mode: set/get, full=True
- [ ] precision_tier: 3 tiers

Stima: **30 check**

### 21.2 close() reset

- [ ] Tutto resettato
- [ ] Re-init automatico dopo close

Stima: **10 check**

### 21.3 Environment variables

- [ ] LIBEPHEMERIS_MODE, LIBEPHEMERIS_PRECISION, LIBEPHEMERIS_LEB

Stima: **6 check**

---

## 22. EphemerisContext

### 22.1 Funzionalita' base

- [ ] calc_ut, houses matches globale
- [ ] set_topo isolato
- [ ] set_sid_mode isolato

Stima: **20 check**

### 22.2 Concorrenza

2 context x 50 iterazioni:
- [ ] Risultati consistenti per context

Stima: **100 check**

---

## 23. Edge Cases

### 23.1 Date limite

5 corpi x 5 date boundary:
- [ ] No crash, gestione corretta

Stima: **25 check**

### 23.2 Input invalidi

- [ ] corpo 999, NaN, Inf, lat 91, mode -1: errore corretto

Stima: **10 check**

### 23.3 Corpi speciali

- [ ] ECL_NUT, Earth geocentrico, Sun helio

Stima: **10 check**

### 23.4 Velocita' vs differenze finite

11 corpi x 50 date:
- [ ] |speed - numerical| < 0.05 deg/day

Stima: 11 x 50 = **550 check**

### 23.5 Wrap-around 360

20 date vicino a 360/0:
- [ ] lon in [0, 360)

Stima: **20 check**

---

## 24. Utility Functions

- [ ] degnorm, radnorm, difdeg2n, split_deg, get_planet_name, swe_version
- [ ] Valori noti per verifica

Stima: **50 check**

---

## 25. Parti Arabe

5 date x 3 localita':
- [ ] calc_all_arabic_parts: Pars Fortunae, Pars Spiritus presenti
- [ ] Valori in [0, 360)

Stima: 5 x 3 x 3 = **45 check**

---

## 26. Golden Regression

100 combinazioni fissate (corpo, data, flag):
- [ ] Risultato identico a referenza salvata

Stima: **100 check**

---

## Riepilogo

| Sezione | Check stimati |
|---------|--------------|
| 1. Posizioni cross-backend | 160.000 |
| 2. Flag combinazioni | 277.800 |
| 3. Case (24 sistemi) | 36.310 |
| 4. Ayanamsha (43 modi) | 2.869 |
| 5. Eclissi solari | 194 |
| 6. Eclissi lunari | 78 |
| 7. Alba/tramonto | 795 |
| 8. Visibilita' eliacale | 75 |
| 9. Stelle fisse | 344 |
| 10. Coordinate | 1.840 |
| 11. Tempo | 3.790 |
| 12. Nodi/apsidi lunari | 3.000 |
| 13. Uraniani/ipotetici | 560 |
| 14. Asteroidi | 525 |
| 15. Lune planetarie | 30 |
| 16. Elementi orbitali | 105 |
| 17. Fenomeni/elongazione | 500 |
| 18. Crossings/stazioni | 78 |
| 19. LEB specifico | 746 |
| 20. Horizons specifico | 148 |
| 21. State management | 46 |
| 22. EphemerisContext | 120 |
| 23. Edge cases | 615 |
| 24. Utility | 50 |
| 25. Parti arabe | 45 |
| 26. Golden regression | 100 |
| **TOTALE** | **~490.763** |
