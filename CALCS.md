# Specifica Tecnica: Libreria Sostitutiva Swiss Ephemeris

## Obiettivo

Creare una libreria Python **pure Python** (nessuna estensione C) che sostituisca 1:1 le API di Swiss Ephemeris per il calcolo delle posizioni dei corpi celesti. La libreria userà **Skyfield** come motore astronomico principale e dati **NASA JPL** come fonte ephemeris.

---

## Architettura generale

### Backend di calcolo

- **Skyfield**: Libreria astronomica pure Python di Brandon Rhodes
- **JPL Ephemeris DE440**: File binari SPK della NASA per posizioni planetarie
- **Formule analitiche**: Per punti derivati (nodi lunari, Lilith, ecc.)
- **Catalogo Hipparcos**: Per stelle fisse

### Dati richiesti

| Tipo | Fonte | Formato | Dimensione |
|------|-------|---------|------------|
| Pianeti principali | NASA JPL | DE440.bsp | ~115 MB |
| Asteroidi/TNO | NASA JPL Horizons | SPK individuali | ~10-50 KB ciascuno |
| Stelle fisse | ESA Hipparcos | CSV/Pickle | ~2 MB |

---

## Specifiche per ogni corpo celeste

### 1. SOLE (Sun)

**Metodo**: Skyfield con DE440

**Procedura**:

1. Caricare ephemeris DE440 contenente posizione baricentrica del Sole
2. Calcolare posizione geocentrica: vettore Terra→Sole
3. Applicare correzione per tempo-luce (aberrazione)
4. Convertire da coordinate cartesiane ICRS a eclittiche

**Output**:

- Longitudine eclittica: conversione da RA/Dec usando obliquità eclittica
- Latitudine eclittica: sempre ~0° per il Sole (per definizione sull'eclittica)
- Distanza: modulo del vettore in AU
- Velocità: differenza finita tra posizioni a ±1 minuto, convertita in °/giorno
- Declinazione: coordinata equatoriale diretta da Skyfield

**Precisione attesa**: < 0.001" (sub-milliarcsecondo)

---

### 2. LUNA (Moon)

**Metodo**: Skyfield con DE440

**Procedura**:

1. DE440 contiene ephemeris lunare ad alta precisione (basata su Lunar Laser Ranging)
2. Calcolare posizione geocentrica diretta
3. La Luna è sufficientemente vicina che il tempo-luce è trascurabile (~1.3 sec)
4. Convertire a coordinate eclittiche

**Velocità**: La Luna si muove ~13°/giorno, calcolare con differenza finita su intervallo breve (1-2 minuti)

**Note speciali**:

- La latitudine eclittica della Luna varia significativamente (±5.15°)
- Importante per calcoli di eclissi e nodi

---

### 3. MERCURIO, VENERE, MARTE, GIOVE, SATURNO, URANO, NETTUNO, PLUTONE

**Metodo**: Skyfield con DE440

**Procedura identica per tutti**:

1. Identificare il target nel file DE440 (es. "mercury", "venus", ecc.)
2. Per i pianeti esterni, utilizzare il **centro planetario** (NAIF ID x99) non il baricentro del sistema
3. Calcolare posizione apparente geocentrica con correzione tempo-luce
4. Convertire a eclittiche

**Mappatura target DE440**:

| Pianeta | Target DE440 | NAIF ID |
|---------|--------------|---------|
| Mercurio | `mercury` | 199 |
| Venere | `venus` | 299 |
| Marte | `mars` | 499 |
| Giove | `jupiter` | 599 |
| Saturno | `saturn` | 699 |
| Urano | `uranus` | 799 |
| Nettuno | `neptune` | 899 |
| Plutone | `pluto` | 999 |

**Nota importante sui baricentri**:
- I baricentri del sistema (NAIF ID x, es. 5 per Jupiter system) includono la massa delle lune
- Per i giganti gassosi, il baricentro può essere migliaia di km dal centro del pianeta
- Il baricentro di Giove-sistema è ~7000 km dal centro di Giove a causa dell'influenza di Io, Europa, Ganimede e Callisto
- Per posizioni planetarie accurate, usare sempre i centri planetari (NAIF ID x99)

**Retrogradazione**: Determinata dal segno della velocità (negativo = retrogrado)

---

### 4. TERRA (Earth)

**Metodo**: Posizione sempre 0 in coordinate geocentriche

Per calcoli eliocentrici:

1. Posizione Terra dal baricentro solare disponibile in DE440
2. Invertire vettore Sole→Terra

---

### 5. NODO NORD LUNARE MEDIO (Mean North Node)

**Metodo**: Formula analitica (Meeus, Astronomical Algorithms, Cap. 47)

**Formula**:

```
T = (JD - 2451545.0) / 36525.0  [secoli giuliani da J2000.0]

Ω = 125.0445479°
    - 1934.1362891° × T
    + 0.0020754° × T²
    + T³ / 467441.0
    - T⁴ / 60616000.0

Risultato = Ω mod 360°
```

**Caratteristiche**:

- Si muove retrogradamente di ~19.35°/anno
- Compie un ciclo completo in ~18.6 anni
- Velocità media: -0.0529539°/giorno

**Velocità**: Derivata della formula o differenza finita

---

### 6. NODO NORD LUNARE VERO (True North Node)

**Metodo**: Calcolo geometrico dall'orbita lunare istantanea

**Procedura**:

1. Ottenere posizione e velocità della Luna in coordinate eclittiche
2. Calcolare il piano orbitale istantaneo della Luna
3. Trovare l'intersezione di questo piano con il piano dell'eclittica
4. Il nodo ascendente è dove la Luna passa da latitudine negativa a positiva

**Calcolo dettagliato**:

1. Vettore posizione Luna: **r** = (x, y, z) in eclittiche
2. Vettore velocità Luna: **v** = (vx, vy, vz)
3. Momento angolare: **h** = **r** × **v** (prodotto vettoriale)
4. Vettore nodale: **n** = **k** × **h** dove **k** = (0, 0, 1) è il polo eclittico
5. Longitudine nodo: Ω = atan2(ny, nx)

**Oscillazione**: Il nodo vero oscilla ±1.5° attorno al nodo medio con periodo ~27 giorni

---

### 7. NODO SUD LUNARE (Mean e True)

**Metodo**: Derivato dal nodo nord

**Formula**:

```
South Node = (North Node + 180°) mod 360°
```

**Velocità e declinazione**: Stesso valore assoluto del nodo nord, segno invertito per declinazione

---

### 8. LILITH MEDIA (Mean Black Moon / Mean Apogee)

**Metodo**: Formula analitica per l'apogeo medio dell'orbita lunare

**Formula** (Swiss Ephemeris documentation):

```
T = (JD - 2451545.0) / 36525.0

Lilith = 83.3532430°
         + 4069.0137111° × T
         - 0.0103238° × T²
         - T³ / 80053.0
         + T⁴ / 18999000.0

Risultato = Lilith mod 360°
```

**Caratteristiche**:

- Velocità media: ~40.7°/anno (0.1114°/giorno)
- Ciclo completo: ~8.85 anni
- Si muove in senso diretto (prograde)

---

### 9. LILITH VERA (True Black Moon / Osculating Apogee)

**Metodo**: Calcolo dell'apogeo osculante dell'orbita lunare

**Procedura**:

1. Calcolare gli elementi orbitali kepleriani istantanei della Luna
2. L'apogeo è il punto più lontano dalla Terra sull'ellisse osculante
3. Longitudine apogeo = longitudine del perigeo + 180°

**Calcolo elementi orbitali**:

1. Dalla posizione **r** e velocità **v** della Luna
2. Calcolare: semi-asse maggiore `a`, eccentricità `e`, argomento del perigeo `ω`
3. Longitudine del perigeo: `ϖ = Ω + ω`
4. Longitudine apogeo: `ϖ + 180°`

**Oscillazione**: Varia significativamente rispetto alla media (±10-15°)

---

### 10. CHIRONE (Chiron)

**Metodo**: File SPK da NASA JPL Horizons

**Procedura per generare SPK**:

1. Accedere a https://ssd.jpl.nasa.gov/horizons/app.html
2. Target: "Chiron" o "2060"
3. Ephemeris Type: "Vector Table"
4. Output: SPK/BSP file
5. Intervallo: -3000 a +3000 CE (o range desiderato)

**Calcolo posizione**:

1. Caricare file SPK in Skyfield
2. Calcolare posizione geocentrica con tempo-luce
3. Convertire a eclittiche

**Limitazioni temporali**: Chirone ha orbita caotica, precisione degrada prima del 650 CE e dopo 4650 CE per incontri ravvicinati con Saturno

---

### 11. PHOLUS

**Metodo**: Identico a Chirone (file SPK)

**ID JPL**: 5145 Pholus

**Limitazioni**: Orbita caotica, limitazioni temporali simili a Chirone

---

### 12. CERERE, PALLADE, GIUNONE, VESTA

**Metodo**: File SPK individuali o file combinato asteroidi principali

**Procedura**:

1. NASA fornisce SPK per i 4 asteroidi principali
2. Alternativamente, DE440 potrebbe includere Cerere (verificare)
3. Stessa procedura dei pianeti

**ID JPL**:

| Asteroide | Numero |
|-----------|--------|
| Cerere | 1 |
| Pallade | 2 |
| Giunone | 3 |
| Vesta | 4 |

---

### 13. TNO: ERIS

**Metodo**: File SPK da JPL Horizons

**ID JPL**: 136199 Eris

**Generazione SPK**:

1. JPL Horizons → Target: "136199"
2. Intervallo temporale: 1800-2200 CE (sufficiente per astrologia moderna)
3. Formato: SPK/BSP

**Caratteristiche orbitali**:

- Semi-asse maggiore: 67.8 AU
- Periodo orbitale: ~559 anni
- Eccentricità: 0.44
- Inclinazione: 44°

---

### 14. TNO: SEDNA

**Metodo**: File SPK da JPL Horizons

**ID JPL**: 90377 Sedna

**Caratteristiche**:

- Orbita estremamente eccentrica (e = 0.855)
- Perielio: 76 AU, Afelio: ~937 AU
- Periodo: ~11,400 anni

**Nota**: Per orbite così eccentriche, gli elementi kepleriani sono imprecisi. File SPK essenziale.

---

### 15. TNO: HAUMEA

**ID JPL**: 136108 Haumea

**Metodo**: File SPK

**Caratteristiche**:

- Pianeta nano con forma ellissoidale
- Periodo: ~284 anni

---

### 16. TNO: MAKEMAKE

**ID JPL**: 136472 Makemake

**Metodo**: File SPK

---

### 17. TNO: IXION

**ID JPL**: 28978 Ixion

**Metodo**: File SPK

**Tipo**: Plutino (risonanza 2:3 con Nettuno)

---

### 18. TNO: ORCUS

**ID JPL**: 90482 Orcus

**Metodo**: File SPK

**Tipo**: Plutino, "anti-Plutone"

---

### 19. TNO: QUAOAR

**ID JPL**: 50000 Quaoar

**Metodo**: File SPK

**Tipo**: Cubewano (oggetto classico della fascia di Kuiper)

---

### 20. STELLE FISSE: REGULUS, SPICA

**Metodo**: Catalogo Hipparcos + correzione per moto proprio e precessione

**Procedura**:

1. Ottenere coordinate J2000 dal catalogo Hipparcos
   - Regulus: HIP 49669 (α Leonis)
   - Spica: HIP 65474 (α Virginis)
2. Applicare moto proprio per epoca di osservazione
3. Applicare precessione da J2000 a data di osservazione
4. Convertire RA/Dec equatoriali a longitudine/latitudine eclittiche

**Dati Hipparcos richiesti**:

- RA, Dec (J2000)
- Moto proprio in RA e Dec (mas/anno)
- Parallasse (per distanza)

**Formula moto proprio**:

```
RA(t) = RA(J2000) + μ_RA × (t - 2000) / cos(Dec)
Dec(t) = Dec(J2000) + μ_Dec × (t - 2000)
```

**Conversione a eclittiche**:

```
ε = obliquità eclittica alla data

sin(β) = sin(Dec) × cos(ε) - cos(Dec) × sin(ε) × sin(RA)
tan(λ) = (sin(RA) × cos(ε) + tan(Dec) × sin(ε)) / cos(RA)
```

Dove β = latitudine eclittica, λ = longitudine eclittica

---

## Calcoli trasversali

### Calcolo della velocità (speed)

**Metodo universale**: Differenza finita centrata

**Procedura**:

1. Calcolare posizione a tempo t - Δt
2. Calcolare posizione a tempo t + Δt
3. Velocità = (pos₂ - pos₁) / (2 × Δt)

**Intervallo raccomandato**: Δt = 1 minuto (1/1440 giorni)

**Gestione wrap-around 360°**:

```
diff = pos₂ - pos₁
if diff > 180: diff -= 360
if diff < -180: diff += 360
velocity = diff / (2 × Δt)
```

---

### Conversione coordinate ICRS → Eclittiche

**Skyfield fornisce** coordinate ICRS (quasi identiche a J2000 equatoriali).

**Conversione a eclittiche**:

1. Ottenere obliquità dell'eclittica per la data: ε
2. Applicare rotazione attorno all'asse X di angolo ε

**Formula obliquità** (IAU 2006):

```
T = (JD - 2451545.0) / 36525.0

ε = 84381.406"
    - 46.836769" × T
    - 0.0001831" × T²
    + 0.00200340" × T³
    - 0.000000576" × T⁴
    - 0.0000000434" × T⁵
```

Convertire da arcsec a gradi dividendo per 3600.

---

### Calcolo declinazione

**Fonte**: Coordinate equatoriali dirette da Skyfield (metodo `radec()`)

La declinazione è già calcolata da Skyfield come parte delle coordinate equatoriali apparenti.

---

### Funzione difdeg2n (differenza angolare)

**Scopo**: Calcolare la differenza angolare più breve tra due posizioni

**Algoritmo**:

```
input: p1, p2 (gradi, 0-360)
diff = (p1 - p2) mod 360
if diff > 180:
    diff = diff - 360
return diff  (range: -180 a +180)
```

---

## Gestione zodiaco siderale

**Nota**: Questo riguarda solo l'offset delle posizioni, non il calcolo.

**Procedura**:

1. Calcolare posizione tropicale (come descritto sopra)
2. Calcolare ayanamsa per la data
3. Longitudine siderale = longitudine tropicale - ayanamsa

**Formula Ayanamsa Lahiri** (esempio):

```
T = (JD - 2451545.0) / 36525.0
ayanamsa = 23.85° + 50.27" × T / 3600
```

Ogni ayanamsa ha la propria formula.

---

## Struttura dei dati ephemeris

### File da includere/scaricare

| File | Contenuto | Fonte | Priorità |
|------|-----------|-------|----------|
| `de440.bsp` | Pianeti principali, Luna, Sole | NASA NAIF | Obbligatorio |
| `chiron.bsp` | Chirone | JPL Horizons | Obbligatorio |
| `pholus.bsp` | Pholus | JPL Horizons | Obbligatorio |
| `asteroids_main.bsp` | Cerere, Pallade, Giunone, Vesta | JPL Horizons | Obbligatorio |
| `eris.bsp` | Eris | JPL Horizons | Obbligatorio |
| `sedna.bsp` | Sedna | JPL Horizons | Obbligatorio |
| `haumea.bsp` | Haumea | JPL Horizons | Obbligatorio |
| `makemake.bsp` | Makemake | JPL Horizons | Obbligatorio |
| `ixion.bsp` | Ixion | JPL Horizons | Obbligatorio |
| `orcus.bsp` | Orcus | JPL Horizons | Obbligatorio |
| `quaoar.bsp` | Quaoar | JPL Horizons | Obbligatorio |
| `hip_main.csv` | Catalogo Hipparcos | ESA | Obbligatorio |

---

## API da implementare (1:1 con Swiss Ephemeris)

### Funzione principale: calc_ut

**Equivalente SWE**: `swe.calc_ut(jd, body_id, flags)`

**Input**:

- Julian Day (UT)
- ID corpo celeste
- Flags (coordinate eclittiche/equatoriali, velocità, ecc.)

**Output** (array 6 elementi):

- [0] Longitudine (o RA)
- [1] Latitudine (o Dec)
- [2] Distanza (AU)
- [3] Velocità in longitudine (°/giorno)
- [4] Velocità in latitudine
- [5] Velocità in distanza

### Altre funzioni da implementare

| Funzione SWE | Descrizione | Metodo |
|--------------|-------------|--------|
| `difdeg2n` | Differenza angolare | Formula matematica |
| `set_ephe_path` | Imposta percorso dati | Gestione path |
| `get_planet_name` | Nome da ID | Lookup table |
| `julday` | Data → Julian Day | Formula standard |
| `revjul` | Julian Day → Data | Formula inversa |

---

## Riepilogo fonti dati

### NASA NAIF/JPL

- **DE440**: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
- **Horizons**: https://ssd.jpl.nasa.gov/horizons/

### ESA

- **Hipparcos**: https://www.cosmos.esa.int/web/hipparcos/catalogues

### Skyfield (libreria Python)

- **PyPI**: `pip install skyfield`
- **Documentazione**: https://rhodesmill.org/skyfield/

---

## Limitazioni rispetto a Swiss Ephemeris

| Funzionalità | Swiss Ephemeris | Nuova libreria |
|--------------|-----------------|----------------|
| Sistemi di case | Tutti | Deve essere implementato separatamente |
| Solar/Lunar returns | `solcross_ut` | Richiede algoritmo di ricerca |
| Coordinate topocentriche | `set_topo` | Supportato da Skyfield |
| Coordinate eliocentriche | Flag | Supportato da Skyfield |
| Zodiaco siderale | Built-in | Solo offset da implementare |

---

## Riepilogo corpi celesti

### Tabella completa

| # | Corpo | Metodo | Fonte dati | Note |
|---|-------|--------|------------|------|
| 0 | Sun | Skyfield | DE440 | Posizione geocentrica |
| 1 | Moon | Skyfield | DE440 | Alta precisione LLR |
| 2 | Mercury | Skyfield | DE440 | Centro planetario (199) |
| 3 | Venus | Skyfield | DE440 | Centro planetario (299) |
| 4 | Mars | Skyfield | DE440 | Centro planetario (499) |
| 5 | Jupiter | Skyfield | DE440 | Centro planetario (599) |
| 6 | Saturn | Skyfield | DE440 | Centro planetario (699) |
| 7 | Uranus | Skyfield | DE440 | Centro planetario (799) |
| 8 | Neptune | Skyfield | DE440 | Centro planetario (899) |
| 9 | Pluto | Skyfield | DE440 | Centro planetario (999) |
| 10 | Mean Node | Formula | Meeus Cap.47 | Analitico |
| 11 | True Node | Geometria | Orbita Luna | Osculante |
| 12 | Mean Lilith | Formula | SWE docs | Apogeo medio |
| 13 | True Lilith | Geometria | Orbita Luna | Apogeo osculante |
| 14 | Earth | Skyfield | DE440 | Per eliocentriche |
| 15 | Chiron | Skyfield | SPK Horizons | Orbita caotica |
| 16 | Pholus | Skyfield | SPK Horizons | |
| 17 | Ceres | Skyfield | SPK Horizons | Pianeta nano |
| 18 | Pallas | Skyfield | SPK Horizons | |
| 19 | Juno | Skyfield | SPK Horizons | |
| 20 | Vesta | Skyfield | SPK Horizons | |
| - | Eris | Skyfield | SPK Horizons | TNO, pianeta nano |
| - | Sedna | Skyfield | SPK Horizons | TNO, orbita estrema |
| - | Haumea | Skyfield | SPK Horizons | TNO, pianeta nano |
| - | Makemake | Skyfield | SPK Horizons | TNO, pianeta nano |
| - | Ixion | Skyfield | SPK Horizons | TNO, plutino |
| - | Orcus | Skyfield | SPK Horizons | TNO, plutino |
| - | Quaoar | Skyfield | SPK Horizons | TNO, cubewano |
| - | Regulus | Hipparcos | ESA | Stella fissa |
| - | Spica | Hipparcos | ESA | Stella fissa |

---

## Dipendenze Python

```toml
[project]
dependencies = [
    "skyfield>=1.54",
    "jplephem>=2.22",
    "numpy>=1.21",
]
```

Tutte le dipendenze sono **pure Python** (nessuna compilazione C richiesta).
