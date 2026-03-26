# Exhaustive Verification Plan — LibEphemeris v1.0.0a3+

Ogni check e' una verifica standalone eseguibile come script Python.
Nessun riferimento a unit test, poe, o pytest. Ogni sezione e' autonoma.

**Contesto:** LibEphemeris e' una libreria di effemeridi astronomiche con 3 backend
(Skyfield/DE440, LEB Chebyshev precomputed, Horizons API), 40+ corpi celesti,
24 sistemi di case, 43 ayanamsha, eclissi, visibilita' eliacale, stelle fisse,
coordinate, e altro. Questo piano verifica la correttezza di OGNI funzione pubblica.

**Riferimento:** pyswisseph 2.10 (Swiss Ephemeris) come implementazione di riferimento.

---

## 1. Posizioni Planetarie — Cross-Backend

### 1.1 Skyfield vs pyswisseph — tutti i corpi principali

Per ogni corpo in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22]:
Per ogni data in 500 JD casuali [2415020, 2488069] (1900-2100):
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` lon diff < 1" vs pyswisseph
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` lat diff < 1" vs pyswisseph
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` dist diff < 1e-5 AU vs pyswisseph
- [ ] `swe_calc_ut(jd, body, SEFLG_SPEED)` speed_lon diff < 0.01 deg/day vs pyswisseph

### 1.2 LEB vs Skyfield — corpi core (14 body)

Per ogni corpo in [0,1,2,3,4,5,6,7,8,9,10,11,12,14]:
Per ogni data in 500 JD in [2396760, 2506330] (base tier):
- [ ] `set_calc_mode("leb")` vs `set_calc_mode("skyfield")`: lon diff < 0.005"
- [ ] lat diff < 0.005"
- [ ] dist diff < 1e-7 AU
- [ ] speed_lon diff < 0.001 deg/day

### 1.3 Horizons vs Skyfield — 13 corpi geocentrici

Per ogni corpo in [0,1,2,3,4,5,6,7,8,9,14,15,17]:
Per ogni data in 200 JD in [2430000, 2470000]:
- [ ] `set_calc_mode("horizons")` vs `set_calc_mode("skyfield")`: lon diff < 0.001"
- [ ] lat diff < 0.001"
- [ ] dist diff < 1e-7 AU

### 1.4 LEB2 vs LEB1 — compressione intatta

Per ogni corpo in LEB2 core (14 body):
Per ogni data in 500 JD:
- [ ] Posizione da LEB2Reader.eval_body() matches LEBReader.eval_body() < 0.001"

### 1.5 Triple cross-validation (Skyfield x LEB x Horizons)

Per 11 corpi x 100 date:
- [ ] |Skyfield - LEB| < 0.005" (lon e lat)
- [ ] |Skyfield - Horizons| < 0.005" (lon e lat)
- [ ] |LEB - Horizons| < 0.005" (lon e lat)

---

## 2. Flag Combinazioni — OGNI flag con OGNI backend

### 2.1 Singoli flag (13 flag x 22 corpi x 3 backend x 50 date)

Flag: SEFLG_SPEED, SEFLG_SPEED|SEFLG_SIDEREAL, SEFLG_SPEED|SEFLG_EQUATORIAL,
SEFLG_SPEED|SEFLG_J2000, SEFLG_SPEED|SEFLG_NOABERR, SEFLG_SPEED|SEFLG_HELCTR,
SEFLG_SPEED|SEFLG_BARYCTR, SEFLG_SPEED|SEFLG_TRUEPOS, SEFLG_SPEED|SEFLG_NONUT,
SEFLG_SPEED|SEFLG_NOGDEFL, SEFLG_SPEED|SEFLG_XYZ, SEFLG_SPEED|SEFLG_RADIANS,
SEFLG_SPEED|SEFLG_TOPOCTR

Per ogni (flag, corpo, backend, data):
- [ ] Nessun crash
- [ ] Ritorna tupla a 6 elementi
- [ ] Tutti i valori finiti (math.isfinite)
- [ ] lon/RA in range valido
- [ ] dist >= 0 (tranne XYZ)

### 2.2 Flag combinati (coppie compatibili)

78 coppie di flag x 5 corpi x 20 date:
- [ ] SEFLG_SPEED|SEFLG_SIDEREAL|SEFLG_EQUATORIAL — nessun crash
- [ ] SEFLG_SPEED|SEFLG_J2000|SEFLG_NOABERR — nessun crash
- [ ] SEFLG_SPEED|SEFLG_HELCTR|SEFLG_J2000 — nessun crash
- [ ] SEFLG_SPEED|SEFLG_SIDEREAL|SEFLG_J2000 — nessun crash
- [ ] ... tutte le 78 coppie

### 2.3 SEFLG_TOPOCTR con set_topo

Per 5 localita' (Roma, Tokyo, NYC, Polo Nord, Equatore):
Per 10 date:
- [ ] set_topo(lon, lat, alt) + calc_ut(jd, SE_MOON, SEFLG_SPEED|SEFLG_TOPOCTR)
- [ ] Risultato diverso dal geocentrico (Moon parallax ~1 deg)
- [ ] Risultato consistente se ripetuto

### 2.4 SEFLG_SIDEREAL con tutte le 43 ayanamsha

Per ogni mode 0-42:
Per 10 date:
Per 3 corpi (Sun, Moon, Mars):
- [ ] set_sid_mode(mode); calc_ut(jd, body, SEFLG_SPEED|SEFLG_SIDEREAL)
- [ ] Risultato != tropicale
- [ ] lon siderale in [0, 360)
- [ ] Ayanamsha = get_ayanamsa_ut(jd) e' finito e ragionevole (0-30 deg per mode moderni)

### 2.5 Sun heliocentric = (0,0,0)

- [ ] calc_ut(jd, SE_SUN, SEFLG_SPEED|SEFLG_HELCTR) == (0,0,0,0,0,0) per qualsiasi jd

### 2.6 SEFLG_XYZ vs conversione manuale

Per 10 corpi x 20 date:
- [ ] pos_ecl = calc_ut(jd, body, SEFLG_SPEED)
- [ ] pos_xyz = calc_ut(jd, body, SEFLG_SPEED|SEFLG_XYZ)
- [ ] x = dist*cos(lat)*cos(lon), y = dist*cos(lat)*sin(lon), z = dist*sin(lat)
- [ ] |x - pos_xyz[0]| < 1e-8

### 2.7 SEFLG_RADIANS vs gradi

Per 10 corpi x 20 date:
- [ ] pos_deg = calc_ut(jd, body, SEFLG_SPEED)
- [ ] pos_rad = calc_ut(jd, body, SEFLG_SPEED|SEFLG_RADIANS)
- [ ] |pos_rad[0] - math.radians(pos_deg[0])| < 1e-10

---

## 3. Sistemi di Case (24 sistemi)

### 3.1 Per ogni sistema: P K O R C E W X M H T B G I i A L N Q Y F U D J

Per ogni sistema:
Per 20 date x 6 localita' (equatore, Roma, Tokyo, NYC, 60N, 60S):
- [ ] houses(jd, lat, lon, system) ritorna (cusps, ascmc)
- [ ] len(cusps) >= 12
- [ ] Tutti i cusps in [0, 360)
- [ ] ASC (ascmc[0]) in [0, 360)
- [ ] MC (ascmc[1]) in [0, 360)
- [ ] ASC != MC (tranne equatore)
- [ ] Cusps ordinati crescenti (con wrap-around a 360)

### 3.2 Confronto houses vs pyswisseph

Per ogni sistema x 10 date x 3 localita':
- [ ] cusps diff < 0.1 deg vs pyswisseph per tutti i cusps
- [ ] ASC diff < 0.01 deg vs pyswisseph
- [ ] MC diff < 0.01 deg vs pyswisseph

### 3.3 houses_ex con flag siderale

Per 5 sistemi (P, K, W, E, B):
Per 3 ayanamsha (Lahiri, Fagan-Bradley, Raman):
Per 10 date:
- [ ] houses_ex(jd, lat, lon, system, SEFLG_SIDEREAL) != houses(jd, lat, lon, system)
- [ ] Differenza = ayanamsha (entro 0.01 deg)

### 3.4 Latitudini polari

- [ ] Placidus a 70N: fallback a Porphyry (no crash)
- [ ] Koch a 70N: fallback a Porphyry (no crash)
- [ ] Equal a 90N: funziona normalmente
- [ ] Whole Sign a 90N: funziona normalmente

### 3.5 houses_armc

Per 5 sistemi x 10 ARMC values x 3 eps values:
- [ ] houses_armc(armc, lat, eps, system) ritorna risultato valido
- [ ] Consistente con houses() per lo stesso momento

### 3.6 house_pos

Per 10 longitudini planetarie x 5 sistemi:
- [ ] house_pos(armc, lat, eps, system, planet_lon) in [1.0, 13.0)
- [ ] Pianeta al cusp 1 -> house_pos ~ 1.0

### 3.7 Gauquelin sectors

Per 5 corpi x 5 date:
- [ ] gauquelin_sector(jd, body, lon, lat, SEFLG_SPEED) in [1, 36]

---

## 4. Ayanamsha (43 modi + custom)

### 4.1 Tutti i modi producono valori diversi

Per JD = 2451545.0 (J2000):
- [ ] get_ayanamsa_ut(jd) per mode 0-42: tutti diversi tra loro
- [ ] Tutti in range [-5, 30] gradi (ragionevole per epoca moderna)

### 4.2 get_ayanamsa_ut vs get_ayanamsa

Per 10 date:
- [ ] |get_ayanamsa_ut(jd_ut) - get_ayanamsa(jd_tt)| < delta_t * derivata (coerenza UT/TT)

### 4.3 get_ayanamsa_ex e get_ayanamsa_ex_ut

Per 5 modi:
- [ ] Ritorna (ayanamsa, retflag) tuple
- [ ] ayanamsa matches get_ayanamsa_ut()

### 4.4 get_ayanamsa_name

Per ogni mode 0-42:
- [ ] get_ayanamsa_name(mode) ritorna stringa non vuota

### 4.5 Custom ayanamsha (mode 255)

- [ ] set_sid_mode(255, t0=2451545.0, ayan_t0=23.5)
- [ ] get_ayanamsa_ut(2451545.0) == 23.5
- [ ] Varia linearmente nel tempo

### 4.6 Star-based ayanamsha (modes 27, 28, 29)

Per ognuno:
- [ ] Calcolo non crasha (richiede star positions)
- [ ] Risultato finito e ragionevole

---

## 5. Eclissi

### 5.1 Eclissi solari — ricerca globale

Per ogni anno 2000-2025 (26 anni):
- [ ] sol_eclipse_when_glob(jd_start) trova almeno 1 eclissi nell'anno
- [ ] JD eclissi e' nell'anno
- [ ] Tipo eclissi e' valido (SE_ECL_TOTAL|ANNULAR|PARTIAL)
- [ ] Coordinate centrali nell'emisfero corretto

### 5.2 Eclissi solari — ricerca locale

Per 3 eclissi note (es. 2024-04-08 USA, 2023-10-14 Americas, 2020-06-21 Africa):
- [ ] sol_eclipse_when_loc(jd, lon, lat) trova l'eclissi
- [ ] Timing entro 1 minuto dal dato noto

### 5.3 Eclissi solari — geometria

Per 5 eclissi note:
- [ ] sol_eclipse_where(jd) ritorna coordinate centrali valide
- [ ] sol_eclipse_how(jd, lon, lat) ritorna magnitudine > 0
- [ ] sol_eclipse_how_details(jd, lon, lat) ritorna dettagli
- [ ] Magnitude, obscuration finiti

### 5.4 Eclissi lunari — ricerca globale

Per ogni anno 2000-2025:
- [ ] lun_eclipse_when(jd_start) trova almeno 1 eclissi
- [ ] Tipo valido (TOTAL|PARTIAL|PENUMBRAL)

### 5.5 Eclissi lunari — dettagli

Per 5 eclissi note:
- [ ] lun_eclipse_how(jd) ritorna gamma, magnitudine
- [ ] Umbral/penumbral magnitude finiti

### 5.6 Occultazioni

Per 3 occultazioni note:
- [ ] lun_occult_when_glob(jd, body) trova occultazione
- [ ] planet_occult_when_glob(jd, body1, body2) funziona

### 5.7 Elementi Besseliani

Per 3 eclissi solari:
- [ ] BesselianElements calcolati (x, y, d, l1, l2, mu)
- [ ] Tutti i valori finiti
- [ ] Derivate temporali (dx_dt, dy_dt) calcolate

### 5.8 Saros e Inex

Per 5 eclissi:
- [ ] get_saros_number(jd) ritorna intero > 0
- [ ] get_inex_number(jd) ritorna intero

---

## 6. Alba/Tramonto/Transito

### 6.1 rise_trans per Sole

Per 5 localita' x 10 date:
- [ ] rise_trans(jd, SE_SUN, lon, lat, SE_CALC_RISE) ritorna JD alba
- [ ] rise_trans(jd, SE_SUN, lon, lat, SE_CALC_SET) ritorna JD tramonto
- [ ] Alba < Tramonto (stesso giorno)
- [ ] Differenza alba-tramonto ragionevole (6-18 ore)
- [ ] rise_trans(jd, SE_SUN, lon, lat, SE_CALC_MTRANSIT) ritorna mezzogiorno ~12h locale

### 6.2 rise_trans per Luna

Per 3 localita' x 10 date:
- [ ] rise_trans(jd, SE_MOON, ..., SE_CALC_RISE) ritorna JD
- [ ] Luna sorge ~50 min dopo ogni giorno

### 6.3 rise_trans per pianeti

Per Mercury, Venus, Mars, Jupiter, Saturn:
Per 5 date:
- [ ] Rise, set, transit calcolati senza crash
- [ ] Transit altitude ragionevole

### 6.4 Crepuscolo

- [ ] BIT_CIVIL_TWILIGHT: sole a -6 deg
- [ ] BIT_NAUTIC_TWILIGHT: sole a -12 deg
- [ ] BIT_ASTRO_TWILIGHT: sole a -18 deg

### 6.5 Latitudini polari

- [ ] Sole a 70N in estate: notte bianca (no set)
- [ ] Sole a 70N in inverno: notte polare (no rise)

---

## 7. Visibilita' Eliacale

### 7.1 heliacal_ut per pianeti

Per Mercury, Venus, Mars, Jupiter, Saturn:
- [ ] heliacal_ut(jd, lon, lat, atmo, observer, body, SE_HELIACAL_RISING) ritorna JD
- [ ] JD e' nel futuro rispetto a jd_start
- [ ] Differenza < 365 giorni

### 7.2 vis_limit_mag

Per 5 corpi x 3 condizioni atmosferiche:
- [ ] vis_limit_mag(jd, geopos, atmo, observer, body) ritorna risultato
- [ ] Magnitudine limite finita

### 7.3 Tipi di evento eliacale

Per Venus:
- [ ] SE_HELIACAL_RISING (levata eliacale mattutina)
- [ ] SE_HELIACAL_SETTING (tramonto eliacale serale)
- [ ] SE_EVENING_FIRST
- [ ] SE_MORNING_LAST

---

## 8. Stelle Fisse

### 8.1 fixstar_ut per stelle note

Per Regulus, Spica, Aldebaran, Sirius, Vega, Antares, Fomalhaut, Pollux:
- [ ] fixstar_ut(jd, name, SEFLG_SPEED) ritorna posizione valida
- [ ] lon in [0, 360), lat in [-90, 90]
- [ ] dist > 0

### 8.2 fixstar_mag

Per 20 stelle:
- [ ] fixstar_mag(name) ritorna magnitudine finita
- [ ] Sirius ~ -1.46, Vega ~ 0.03 (entro 0.1)

### 8.3 Moto proprio

Per Sirius e Barnard's Star (alto moto proprio):
- [ ] Posizione a J2000 vs J2100: differenza misurabile
- [ ] Barnard's Star: moto > 1" per anno

### 8.4 fixstar2 (ricerca per numero)

Per 10 stelle Hipparcos:
- [ ] fixstar2(jd, hip_number, flags) ritorna posizione
- [ ] Matches fixstar_ut con stesso nome

---

## 9. Trasformazioni di Coordinate

### 9.1 cotrans round-trip

Per 3 obliquita' x 100 punti casuali (lon, lat):
- [ ] ecl -> eq -> ecl: |lon_finale - lon_iniziale| < 1e-10 deg
- [ ] ecl -> eq -> ecl: |lat_finale - lat_iniziale| < 1e-10 deg

### 9.2 cotrans_sp con velocita'

Per 10 punti + velocita':
- [ ] Round-trip preserva posizione e velocita'

### 9.3 azalt

Per 5 localita' x 10 date x 5 corpi:
- [ ] azalt(jd, SE_ECL2HOR, lon, lat, alt, ecl_lon, ecl_lat) ritorna (az, alt_true, alt_app)
- [ ] Azimut in [0, 360)
- [ ] Altitudine in [-90, 90]

### 9.4 azalt_rev

Per 10 punti (az, alt):
- [ ] azalt_rev(jd, SE_HOR2ECL, lon, lat, alt, az, true_alt) ritorna (ecl_lon, ecl_lat)
- [ ] Round-trip: azalt -> azalt_rev -> azalt: consistente

### 9.5 refrac

- [ ] refrac(alt, SE_TRUE_TO_APP) > alt per alt > 0
- [ ] refrac(alt, SE_APP_TO_TRUE) < alt per alt > 0
- [ ] All'orizzonte (alt=0): refrazione ~ 34'

---

## 10. Funzioni Tempo

### 10.1 julday/revjul round-trip

Per 2000 date casuali (anno -5000 a +5000):
- [ ] revjul(julday(y, m, d, h)) == (y, m, d, h) entro 1e-6

### 10.2 deltat

Per 500 JD in [1900, 2100]:
- [ ] deltat(jd) > 0
- [ ] deltat(jd) < 0.01 giorni (~14 min)
- [ ] Variazione smooth (no discontinuita')

### 10.3 deltat_ex

Per 10 JD:
- [ ] deltat_ex(jd, SEFLG_JPLEPH) ritorna (delta_t, retflag)
- [ ] delta_t matches deltat(jd)

### 10.4 sidtime

Per 50 JD:
- [ ] sidtime(jd) in [0, 24) ore
- [ ] Incremento ~3m56s per giorno solare

### 10.5 utc_to_jd e inversi

Per 10 date UTC:
- [ ] utc_to_jd(y, m, d, h, min, sec) ritorna (jd_et, jd_ut)
- [ ] jdet_to_utc(jd_et) ritorna data originale
- [ ] jdut1_to_utc(jd_ut) ritorna data originale

### 10.6 day_of_week

- [ ] day_of_week(julday(2024, 1, 1, 0)) == 1 (Monday)
- [ ] day_of_week(julday(2024, 3, 26, 0)) == 2 (Tuesday)

### 10.7 time_equ

Per 10 date:
- [ ] time_equ(jd) in [-0.02, 0.02] giorni (equazione del tempo < 17 min)

### 10.8 TAI functions

- [ ] utc_to_tai_jd(jd) > jd (TAI ahead of UTC)
- [ ] tai_jd_to_utc(utc_to_tai_jd(jd)) == jd (round-trip)
- [ ] tt_to_tai_jd(jd) < jd (TT = TAI + 32.184s)
- [ ] tai_to_tt_jd(tt_to_tai_jd(jd)) == jd (round-trip)

---

## 11. Nodi e Apsidi Lunari

### 11.1 Nodo Medio

Per 200 date:
- [ ] calc_mean_lunar_node(jd) in [0, 360)
- [ ] Moto retrogrado: nodo(jd+1) < nodo(jd) (in media)
- [ ] Velocita' ~ -0.053 deg/day

### 11.2 Nodo Vero

Per 200 date:
- [ ] calc_true_lunar_node(jd) in [0, 360)
- [ ] Oscilla attorno al nodo medio (entro ~1.5 deg)
- [ ] Distanza (3o elemento) ragionevole (~0.0025 AU)

### 11.3 Lilith Media (Black Moon)

Per 200 date:
- [ ] calc_mean_lilith(jd) in [0, 360)
- [ ] Moto diretto: ~40 deg/anno

### 11.4 Lilith Vera (Osculating Apogee)

Per 100 date:
- [ ] calc_ut(jd, SE_OSCU_APOG, SEFLG_SPEED): lon in [0, 360)

### 11.5 Apogee/Perigee interpolati

Per 50 date:
- [ ] calc_ut(jd, SE_INTP_APOG, SEFLG_SPEED): lon in [0, 360)
- [ ] calc_ut(jd, SE_INTP_PERG, SEFLG_SPEED): lon in [0, 360)
- [ ] INTP_APOG e INTP_PERG differiscono di ~180 deg (opposti)

### 11.6 nod_aps_ut

Per Sun, Moon, Mars, Jupiter:
- [ ] nod_aps_ut(jd, body, SEFLG_SPEED, SE_NODBIT_MEAN) ritorna 4 tuple
- [ ] Ascending node, descending node, perihelion, aphelion
- [ ] Tutti i valori finiti

---

## 12. Corpi Uraniani / Ipotetici

### 12.1 Uraniani geocentrici (40-47)

Per Cupido(40), Hades(41), Zeus(42), Kronos(43), Apollon(44), Admetos(45), Vulkanus(46), Poseidon(47):
Per 20 date:
- [ ] calc_ut(jd, body, SEFLG_SPEED): lon in [0, 360)
- [ ] dist > 10 AU (sono lontani)

### 12.2 Uraniani eliocentrici

Per tutti i body 40-47:
Per 10 date:
- [ ] calc_ut(jd, body, SEFLG_SPEED|SEFLG_HELCTR): lon in [0, 360)
- [ ] Matches geocentric entro ~1 deg (parallasse trascurabile a 40+ AU)

### 12.3 Transpluto (48)

Per 20 date:
- [ ] calc_ut(jd, 48, SEFLG_SPEED): geocentrico funziona
- [ ] calc_ut(jd, 48, SEFLG_SPEED|SEFLG_HELCTR): eliocentrico funziona

### 12.4 Uraniani siderali

Per body 40-47:
Per 3 ayanamsha:
- [ ] calc_ut(jd, body, SEFLG_SPEED|SEFLG_SIDEREAL): lon diverso da tropicale

---

## 13. Asteroidi e Corpi Minori

### 13.1 Asteroidi principali

Per Chiron(15), Ceres(17), Pallas(18), Juno(19), Vesta(20):
Per 50 date:
- [ ] calc_ut(jd, body, SEFLG_SPEED): posizione valida
- [ ] Confronto vs pyswisseph: lon diff < 2"

### 13.2 Via SE_AST_OFFSET

- [ ] calc_ut(jd, 10001, flags) == calc_ut(jd, 17, flags) (Ceres)
- [ ] calc_ut(jd, 10002, flags) == calc_ut(jd, 18, flags) (Pallas)
- [ ] calc_ut(jd, 10003, flags) == calc_ut(jd, 19, flags) (Juno)
- [ ] calc_ut(jd, 10004, flags) == calc_ut(jd, 20, flags) (Vesta)
- [ ] calc_ut(jd, 12060, flags) == calc_ut(jd, 15, flags) (Chiron)

### 13.3 TNO noti

Per SE_ERIS, SE_SEDNA, SE_MAKEMAKE, SE_HAUMEA (se SPK disponibili):
- [ ] calc_ut(jd, body, SEFLG_SPEED): posizione valida o errore esplicito

---

## 14. Lune Planetarie

Per SE_MOON_IO, SE_MOON_EUROPA, SE_MOON_GANYMEDE, SE_MOON_CALLISTO:
Per SE_MOON_TITAN:
Per SE_MOON_TRITON:
- [ ] calc_ut(jd, body, SEFLG_SPEED): posizione valida (se SPK registrato) o errore
- [ ] is_planetary_moon(body) == True

---

## 15. Elementi Orbitali

### 15.1 get_orbital_elements_ut

Per Sun, Mercury, Venus, Mars, Jupiter, Saturn:
Per 5 date:
- [ ] Ritorna tupla con semi-asse maggiore, eccentricita', inclinazione, etc.
- [ ] Semi-asse maggiore ragionevole (Marte ~1.52 AU, Giove ~5.2 AU)
- [ ] Eccentricita' in [0, 1)

### 15.2 orbit_max_min_true_distance

Per 5 corpi:
- [ ] max_dist > min_dist
- [ ] Rapporto max/min ragionevole (Terra: ~1.034)

### 15.3 nod_aps_ut

Per 5 corpi x 3 flag (NODBIT_MEAN, NODBIT_OSCU, NODBIT_FOPOINT):
- [ ] Ritorna (ascending_node, descending_node, perihelion, aphelion)
- [ ] Tutti valori finiti

---

## 16. Fenomeni e Elongazione

### 16.1 pheno_ut

Per 10 corpi x 5 date:
- [ ] pheno_ut(jd, body, flags) ritorna tupla con phase angle, elongation, magnitude, etc.
- [ ] Phase angle in [0, 180]
- [ ] Elongation in [0, 180]

### 16.2 Elongazione dal Sole

Per Mercury e Venus:
Per 50 date:
- [ ] get_elongation_from_sun(jd, body) in [0, 180]
- [ ] Mercury: max ~28 deg
- [ ] Venus: max ~47 deg

### 16.3 Stella del mattino/sera

Per Venus:
Per 10 date:
- [ ] is_morning_star(jd, SE_VENUS) o is_evening_star(jd, SE_VENUS) == True
- [ ] Non entrambi True

---

## 17. Crossings e Stazioni

### 17.1 Equinozi e solstizi

Per anno 2020-2025:
- [ ] solcross_ut(0, jd_start) trova equinozio di primavera (lon=0)
- [ ] solcross_ut(90, jd_start) trova solstizio d'estate (lon=90)
- [ ] solcross_ut(180, jd_start) trova equinozio d'autunno
- [ ] solcross_ut(270, jd_start) trova solstizio d'inverno
- [ ] Date entro 1 giorno dai valori noti

### 17.2 Moon crossings

Per 5 date:
- [ ] mooncross_ut(lon_target, jd) trova il crossing
- [ ] La Luna e' effettivamente a lon_target al JD ritornato

### 17.3 Stazioni e retrogradazioni

Per Marte e Mercurio:
- [ ] find_station_ut(jd, body) trova la prossima stazione
- [ ] next_retrograde_ut(jd, body) trova inizio retrogradazione
- [ ] Velocita' ~ 0 alla stazione

---

## 18. LEB-Specifico

### 18.1 open_leb() factory

- [ ] open_leb("data/leb/ephemeris_base.leb") ritorna LEBReader
- [ ] open_leb("data/leb2/base_core.leb") ritorna LEB2Reader o CompositeLEBReader
- [ ] open_leb("nonexistent.leb") solleva FileNotFoundError
- [ ] open_leb("README.md") solleva ValueError (magic bytes)

### 18.2 LEBReader

- [ ] .path ritorna il path corretto
- [ ] .jd_range ritorna (start, end) valido
- [ ] .has_body(0) == True (Sun)
- [ ] .has_body(999) == False
- [ ] .eval_body(0, 2451545.0) ritorna (pos, vel) tuple
- [ ] .eval_nutation(2451545.0) ritorna (dpsi, deps) finiti
- [ ] .delta_t(2451545.0) > 0
- [ ] .close() non solleva eccezione
- [ ] Context manager funziona

### 18.3 LEB2Reader

Stessi test del LEBReader +
- [ ] Lazy decompression: prima chiamata eval_body piu' lenta
- [ ] Seconda chiamata stessa body: dalla cache

### 18.4 CompositeLEBReader

- [ ] from_file_with_companions("base_core.leb") carica companion files se presenti
- [ ] has_body dispatch al file corretto
- [ ] eval_body cross-file: core body vs asteroid body

### 18.5 LEB fast_calc pipeline

Per 14 corpi core x 50 date:
- [ ] fast_calc_ut(reader, jd, body, flags) matches calc_ut skyfield < 0.005"
- [ ] Nutation correction applicata (lon nuda + dpsi)
- [ ] Sidereal correction applicata quando richiesto

### 18.6 Bundled LEB2

- [ ] `libephemeris/data/leb2/base_core.leb` esiste nel pacchetto
- [ ] `_discover_leb_file()` lo trova come fallback
- [ ] Dopo `pip install`, LEB fast path funziona senza download

---

## 19. Horizons-Specifico

### 19.1 HorizonsClient

- [ ] fetch_state_vector("399", jd, "@0", "TDB") ritorna 6 componenti
- [ ] Cache hit: seconda chiamata stessi parametri ritorna stesso risultato
- [ ] clear_cache() svuota la cache
- [ ] shutdown() e' idempotente

### 19.2 Pipeline geocentrica

Per 5 corpi x 10 date:
- [ ] horizons_calc_ut(jd, body, SEFLG_SPEED): risultato valido
- [ ] Light-time correction applicata
- [ ] Aberration correction applicata
- [ ] Frame rotation corretta (ICRS -> ecliptic of date)

### 19.3 Corpi analitici (no HTTP)

- [ ] horizons_calc_ut(jd, 10, flags): Mean Node via formula analitica
- [ ] horizons_calc_ut(jd, 12, flags): Mean Apogee via formula analitica
- [ ] Nessuna HTTP call per questi corpi

### 19.4 Fallback a Skyfield

- [ ] horizons_calc_ut(jd, 11, flags): True Node -> KeyError -> Skyfield
- [ ] horizons_calc_ut(jd, 13, flags): Oscu Apogee -> KeyError -> Skyfield
- [ ] SEFLG_TOPOCTR -> KeyError -> Skyfield

### 19.5 Auto mode

- [ ] set_calc_mode("auto"): con LEB -> usa LEB
- [ ] set_calc_mode("auto"): senza LEB, senza DE440 -> usa Horizons
- [ ] set_calc_mode("auto"): senza LEB, con DE440 -> usa Skyfield

---

## 20. State Management

### 20.1 Setter/getter round-trip

- [ ] set_calc_mode("skyfield"); get_calc_mode() == "skyfield"
- [ ] set_calc_mode("leb"); get_calc_mode() == "leb"
- [ ] set_calc_mode("horizons"); get_calc_mode() == "horizons"
- [ ] set_calc_mode("auto"); get_calc_mode() == "auto"
- [ ] set_calc_mode(None); get_calc_mode() == "auto" (default)
- [ ] set_calc_mode("invalid") solleva ValueError
- [ ] set_topo(12.5, 41.9, 0); get_topo() != None
- [ ] set_sid_mode(1); get_sid_mode() == 1
- [ ] set_sid_mode(0, 2451545.0, 23.5); get_sid_mode(full=True) == (0, 2451545.0, 23.5)

### 20.2 close() reset

- [ ] close() resetta topo, sidereal mode, LEB reader, Horizons client
- [ ] Dopo close(), calc_ut funziona ancora (re-init automatico)

### 20.3 Environment variables

- [ ] LIBEPHEMERIS_MODE=skyfield: get_calc_mode() == "skyfield"
- [ ] LIBEPHEMERIS_PRECISION=base: get_precision_tier() == "base"
- [ ] LIBEPHEMERIS_LEB=/path: LEB reader usa quel file

### 20.4 Precision tier

- [ ] set_precision_tier("base"); get_precision_tier() == "base"
- [ ] set_precision_tier("medium"); get_precision_tier() == "medium"
- [ ] set_precision_tier("extended"); get_precision_tier() == "extended"
- [ ] Tier base: range ~1850-2150
- [ ] Tier medium: range ~1550-2650

---

## 21. EphemerisContext (Thread-Safe)

### 21.1 Funzionalita' base

- [ ] ctx = EphemerisContext()
- [ ] ctx.calc_ut(jd, SE_SUN, SEFLG_SPEED) matches swe_calc_ut()
- [ ] ctx.houses(jd, lat, lon, b"P") matches swe_houses()

### 21.2 Isolamento stato

- [ ] ctx1.set_topo(Roma); ctx2.set_topo(Tokyo)
- [ ] ctx1.calc_ut(jd, SE_MOON, SEFLG_TOPOCTR) != ctx2.calc_ut()
- [ ] Stato globale non influenzato

### 21.3 Sidereal isolation

- [ ] ctx1.set_sid_mode(1); ctx2.set_sid_mode(0)
- [ ] ctx1 risultati diversi da ctx2

### 21.4 LEB per-context

- [ ] ctx1.set_leb_file("path1.leb"); ctx2 usa LEB diverso

---

## 22. Edge Cases

### 22.1 Date limite

- [ ] calc_ut(2396758.5, SE_SUN, 0): inizio base tier -> no crash
- [ ] calc_ut(2506331.5, SE_SUN, 0): fine base tier -> no crash
- [ ] calc_ut(2451545.0, SE_SUN, 0): J2000 -> risultato noto (~280.37 deg)
- [ ] calc_ut(1000000.0, SE_SUN, 0): fuori range -> EphemerisRangeError o fallback

### 22.2 Input invalidi

- [ ] calc_ut(jd, 999, 0): corpo sconosciuto -> errore
- [ ] calc_ut(float('nan'), SE_SUN, 0): NaN -> errore o gestione
- [ ] calc_ut(float('inf'), SE_SUN, 0): Inf -> errore o gestione
- [ ] houses(jd, 91, 0, b"P"): lat fuori range -> CoordinateError
- [ ] set_sid_mode(-1): modo invalido -> gestione

### 22.3 Corpi speciali

- [ ] calc_ut(jd, SE_ECL_NUT, 0): ritorna nutation/obliquity
- [ ] calc_ut(jd, SE_EARTH, SEFLG_SPEED): geocentrico Earth = (0,0,0)
- [ ] calc_ut(jd, SE_SUN, SEFLG_HELCTR): helio Sun = (0,0,0)

### 22.4 Velocita' vs differenze finite

Per 11 corpi x 50 date:
- [ ] speed = calc_ut(jd, body, SEFLG_SPEED)[0][3]
- [ ] numerical_speed = (calc_ut(jd+0.001)[0][0] - calc_ut(jd-0.001)[0][0]) / 0.002
- [ ] |speed - numerical_speed| < 0.05 deg/day

### 22.5 Wrap-around 360 deg

Per 10 date dove il Sole e' vicino a 360/0:
- [ ] lon in [0, 360) (mai negativo, mai >= 360)

---

## 23. Utility Functions

### 23.1 degnorm / radnorm

- [ ] degnorm(361) == 1
- [ ] degnorm(-1) == 359
- [ ] radnorm(7) ~ 7 - 2*pi

### 23.2 difdeg2n

- [ ] difdeg2n(350, 10) == -20 (shortest arc)
- [ ] difdeg2n(10, 350) == 20

### 23.3 split_deg

- [ ] split_deg(123.456, SPLIT_DEG_ZODIACAL): ritorna (deg, sign, min, sec, ...)
- [ ] split_deg(90.5, 0): ritorna componenti

### 23.4 get_planet_name

Per ogni corpo 0-22:
- [ ] get_planet_name(body) ritorna stringa non vuota

### 23.5 swe_version

- [ ] swe_version() ritorna stringa versione

---

## 24. Parti Arabe

### 24.1 calc_all_arabic_parts

Per 5 date x 3 localita':
- [ ] calc_all_arabic_parts(jd, lon, lat) ritorna dizionario
- [ ] Pars Fortunae presente
- [ ] Pars Spiritus presente
- [ ] Tutti i valori in [0, 360)

---

## 25. Golden Regression

### 25.1 Stabilita' posizioni

Per 100 combinazioni (corpo, data, flag) fissate:
- [ ] Risultato identico bit-per-bit a referenza salvata
- [ ] Nessuna regressione rispetto all'ultima release

---

## Riepilogo Stime

| Sezione | Check stimati |
|---------|--------------|
| 1. Posizioni cross-backend | ~50,000 |
| 2. Flag combinazioni | ~100,000 |
| 3. Case (24 sistemi) | ~30,000 |
| 4. Ayanamsha (43 modi) | ~5,000 |
| 5. Eclissi | ~500 |
| 6. Alba/tramonto | ~1,000 |
| 7. Visibilita' eliacale | ~200 |
| 8. Stelle fisse | ~500 |
| 9. Coordinate | ~3,000 |
| 10. Tempo | ~5,000 |
| 11. Nodi/apsidi lunari | ~2,000 |
| 12. Uraniani | ~1,000 |
| 13. Asteroidi | ~1,000 |
| 14. Lune planetarie | ~100 |
| 15. Elementi orbitali | ~500 |
| 16. Fenomeni | ~500 |
| 17. Crossings | ~200 |
| 18. LEB specifico | ~2,000 |
| 19. Horizons specifico | ~500 |
| 20. State management | ~200 |
| 21. EphemerisContext | ~200 |
| 22. Edge cases | ~500 |
| 23. Utility | ~100 |
| 24. Parti arabe | ~100 |
| 25. Golden regression | ~100 |
| **TOTALE** | **~204,200** |
