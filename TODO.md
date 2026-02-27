# TODO — LibEphemeris Roadmap

> Ultimo aggiornamento: Febbraio 2026
>
> Questo file elenca tutti i task rimanenti in ordine di priorita.
> Ogni task include contesto, file coinvolti, e passi di implementazione.

---

## Indice

1. [LEB: Copertura asteroidi per-body (Opzione 2)](#1-leb-copertura-asteroidi-per-body-opzione-2)
2. [LEB: Velocita di generazione](#2-leb-velocità-di-generazione)
3. [Architettura: Disaccoppiamento modalita SPK+Kepler vs LEB](#3-architettura-disaccoppiamento-modalità-spkkepler-vs-leb)
4. [Precisione: Miglioramento fallback Kepleriano](#4-precisione-miglioramento-fallback-kepleriano)

---

## 1. LEB: Copertura asteroidi per-body (Opzione 2)

### Contesto

I file SPK degli asteroidi (scaricati da JPL Horizons) coprono un range limitato
(circa 1600-2500 CE), mentre le effemeridi planetarie (DE440/DE441) coprono range
molto piu ampi (DE441: -13200 a +17191). Oggi il generatore LEB prova a scaricare
SPK per l'intero range del tier: se il tier e "extended" (-5000 a +5000), chiede a
Horizons 10.000 anni di copertura per Ceres — Horizons rifiuta, e l'asteroide viene
silenziosamente escluso dal file LEB o cade nel fallback Kepleriano (errori di gradi).

### Soluzione scelta (Opzione 2)

Ogni body nel file LEB avra il proprio range di date indipendente. I pianeti
copriranno l'intero range del tier; gli asteroidi copriranno solo il range
effettivamente supportato da SPK (circa 1600-2500). A runtime, se la data
richiesta e fuori dal range di un asteroide nel LEB, la libreria fara fallback
al calcolo al volo (SPK se disponibile, altrimenti Kepler).

Il formato binario LEB **gia supporta** range per-body — ogni body entry ha i
campi `jd_start` e `jd_end` indipendenti (`leb_format.py:BodyEntry`).

### Task

#### 1.1 Scoprire il range massimo SPK da Horizons

**File:** `scripts/generate_leb.py`, `libephemeris/minor_bodies.py`

Oggi `auto_download_asteroid_spk()` richiede il range esatto passato dal chiamante.
Se Horizons non puo fornire quel range, il download fallisce silenziosamente.

Implementare una funzione che:
1. Prova a scaricare l'SPK per il range richiesto
2. Se fallisce, prova con range via via piu stretti (binary search o range noti)
3. Dopo il download, ispeziona il file SPK scaricato per determinare il range
   effettivo (gia possibile con `spktype21.SPKType21.open()` + lettura `seg.start_second`
   e `seg.end_second`)
4. Ritorna il range effettivo

Alternativa piu semplice: hardcodare i range noti per i 5 asteroidi principali
(Ceres, Pallas, Juno, Vesta, Chiron) basandosi sui dati empirici dei file gia
scaricati (~1600-2500). Meno elegante ma affidabile.

**File coinvolti:**
- `libephemeris/minor_bodies.py` — `auto_download_asteroid_spk()` (linea 2284)
- `libephemeris/spk.py` — `download_and_register_spk()` (linea 1073)

#### 1.2 Generare asteroidi con range indipendente

**File:** `scripts/generate_leb.py`

Modificare `generate_body_icrs_asteroid()` (linea 733) per accettare un range
diverso da quello del tier. La funzione oggi riceve `jd_start`/`jd_end` dal tier;
deve invece ricevere il range effettivo dell'SPK.

Modificare `assemble_leb()` (linea 1341) per:
1. Per ogni asteroide, determinare il range SPK disponibile (da task 1.1)
2. Passare quel range a `generate_body_icrs_asteroid()` invece del range del tier
3. Scrivere nel body index entry il range effettivo dell'asteroide

**Attenzione:** la nutation e delta-T devono coprire il range piu ampio (quello
del tier), non quello degli asteroidi.

#### 1.3 Aggiornare il reader per range per-body

**File:** `libephemeris/leb_reader.py`

Verificare che `eval_body()` (linea ~200) gestisca correttamente il caso in cui
il JD richiesto e dentro il range globale del file ma fuori dal range di un
body specifico. Oggi lancia `EphemerisRangeError` — deve invece ritornare un
segnale che permetta al chiamante di fare fallback a Skyfield.

**File coinvolti:**
- `libephemeris/leb_reader.py` — `eval_body()` 
- `libephemeris/fast_calc.py` — pipeline A/B/C, gestione fallback

#### 1.4 Fallback trasparente per body fuori range

**File:** `libephemeris/fast_calc.py`, `libephemeris/planets.py`

Quando `eval_body()` segnala "fuori range per questo body", il codice deve
fare fallback al path Skyfield/SPK/Kepler **solo per quel body**, senza
disabilitare il LEB per gli altri body nella stessa chiamata.

Oggi il dispatch in `swe_calc_ut()` (`planets.py:772-782`) e binario: o usa
LEB o usa Skyfield. Deve diventare per-body: LEB se disponibile e in range,
Skyfield altrimenti.

#### 1.5 Abortire se SPK manca (strict mode per generatore)

**File:** `scripts/generate_leb.py`

Oggi il generatore fa `set_strict_precision(False)` (linea 1369) e puo
silenziosamente usare il fallback Keplerian per generare i coefficienti
degli asteroidi. Questo produce dati imprecisi nel LEB senza alcun avviso.

Cambiare il comportamento: se un asteroide non ha SPK con copertura
sufficiente, il generatore deve **abortire con errore chiaro** invece di
generare dati Kepleriani. L'asteroide deve essere escluso dal file LEB
oppure generato con il range ridotto (task 1.2), mai con dati Kepleriani.

Stampare un warning chiaro nell'output che indica per ogni asteroide:
- `(spk21)` — dati SPK, precisione sub-arcsecond
- `(EXCLUDED)` — escluso per mancanza SPK

#### 1.6 Migliorare verify_leb() per asteroidi e eclittici

**File:** `scripts/generate_leb.py`

La funzione `verify_leb()` (linea 1933) oggi:
- Per i pianeti ICRS: confronta con Skyfield live — **funziona**
- Per gli asteroidi ICRS: non li confronta (non sono in `_PLANET_MAP`) — **buco**
- Per i body eclittici: controlla solo `0 <= lon < 360` — **insufficiente**

Migliorare:
1. Asteroidi: confrontare con `spktype21` o con `swe_calc()` con strict on
2. Eclittici: richiamare la funzione analitica di riferimento e confrontare
3. Report chiaro con errore massimo in arcsec per ogni body

#### 1.7 Comandi poe per medium e extended con range adattivo

**File:** `pyproject.toml`

I comandi poe per group generation (`leb:generate:medium:groups`,
`leb:generate:extended:groups`) gia esistono. Dopo l'implementazione del
range per-body, funzioneranno automaticamente perche la logica e nel
generatore. Verificare e aggiornare la documentazione.

#### 1.8 Aggiornare documentazione

**File:** `docs/LEB_GUIDE.md`, `docs/LEB_PLAN.md`

Documentare:
- Il concetto di range per-body nel formato LEB
- I limiti di copertura SPK per gli asteroidi
- Il comportamento di fallback quando un body e fuori range
- Tabella dei range effettivi per ogni tier e body

---

## 2. LEB: Velocita di generazione

### Contesto

La generazione del file LEB ha tre livelli di velocita molto diversi:

| Gruppo | Tempo (base 300yr) | Motivo |
|--------|-------------------|--------|
| Pianeti (11 body) | ~15s | Vettorizzato: un singolo `target.at(t_array)` Skyfield |
| Asteroidi (5 body) | ~3-5 min | Scalare: `spktype21.compute_type21()` non ha API batch |
| Analytical eclittici (6 body) | ~5-8 min | Scalare: 328K chiamate Skyfield punto per punto |
| Analytical Uraniani (9 body) | ~10s | Scalare ma pura aritmetica, veloce |

Il bottleneck principale sono i 6 body eclittici (true node, true lilith,
osculating apogee, interp apogee, interp perigee, mean apogee). Ognuno chiama
`(moon - earth).at(t)` di Skyfield ~328.000 volte (una per ogni punto di fitting +
verifica). Skyfield supporta nativamente array di tempi, ma le nostre funzioni
(`calc_true_lunar_node()`, `calc_true_lilith()`, ecc.) accettano solo scalari.

Il pattern per la vettorizzazione **esiste gia** nel generatore:
`generate_nutation()` usa `_compute_all_segment_jds()` per pre-calcolare tutti i
JD, `erfa.nut06a()` per valutarli in batch, e `_fit_and_verify_from_values()` per
fittare senza ulteriori chiamate alla funzione.

### Task

#### 2.1 Vettorizzare calc_true_lunar_node (body 11)

**File:** `scripts/generate_leb.py` (nuova funzione batch, privata del generatore)

Scrivere `_eval_true_node_batch(jd_array)` nel generatore che:
1. Chiama `(moon - earth).at(t_array)` una sola volta con tutti i JD
2. Ottiene `frame_xyz(ecliptic_frame)` e `frame_xyz_and_velocity(ecliptic_frame)` — vettorizzati
3. Calcola il cross product h = r x v con numpy (element-wise su array `(3, N)`)
4. Calcola `atan2(h_x, -h_y)` con `np.arctan2` — vettorizzato
5. Ritorna array `(N, 3)` di `[lon, lat, dist]`

**Non toccare** `libephemeris/lunar.py` — la funzione pubblica `calc_true_lunar_node()`
resta scalare per compatibilita API. Il wrapper batch e privato del generatore.

**Riferimento:** `lunar.py:2007-2043` per la logica scalare da replicare.

**Speedup atteso:** da ~2 min a ~2s per questo body (100-150x).

#### 2.2 Vettorizzare calc_true_lilith (body 13)

**File:** `scripts/generate_leb.py` (nuova funzione batch)

Stessa struttura del task 2.1. `calc_true_lilith()` (`lunar.py:2188-2249`) fa:
1. Skyfield: `(moon - earth).at(t)` — vettorizzabile
2. Cross product h = r x v — numpy element-wise
3. Eccentricity vector e = (v x h)/mu - r/|r| — numpy element-wise
4. Apogee = -e, coordinate sferiche — `np.arctan2`, `np.arcsin`

Tutto il post-processing e aritmetica vettoriale, zero branching.

**Speedup atteso:** ~100x.

#### 2.3 Vettorizzare calc_interpolated_apogee e perigee (body 21, 22)

**File:** `scripts/generate_leb.py`

`calc_interpolated_apogee()` (`lunar.py:2658`) calcola:
- Longitudine: pura math polinomiale (trivialmente vettorizzabile)
- Latitudine e distanza: chiama `calc_true_lilith()` (vettorizzabile con 2.2)

`calc_interpolated_perigee()` (`lunar.py:3000`):
- Longitudine: pura math polinomiale
- Latitudine e distanza: chiama `calc_osculating_perigee()` → che chiama Skyfield

Scrivere batch wrappers che riusano `_eval_true_lilith_batch()` dal task 2.2
per la parte lat/dist, e calcolano la longitudine con numpy.

#### 2.4 Vettorizzare body polinomiali (body 10, 12)

**File:** `scripts/generate_leb.py`

`calc_mean_lunar_node()` e `calc_mean_lilith()` sono puri polinomi in `T`:
```python
Omega = 125.04 - 1934.14*T + 0.0021*T**2 + T**3/467441 - T**4/60616000
```

Trivialmente vettorizzabili: passare un array di T a numpy, ottenere un array
di risultati. Speedup enorme ma il tempo assoluto e gia basso (~secondi).

#### 2.5 Vettorizzare osculating apogee (body 13) 

**File:** `scripts/generate_leb.py`

`calc_oscu_apog()` in `lunar.py` chiama Skyfield + calcolo vettore
eccentricita. Stessa struttura di `calc_true_lilith()` — riusare il batch
wrapper del task 2.2 con piccole variazioni.

#### 2.6 Nuovo generate_body_ecliptic_vectorized()

**File:** `scripts/generate_leb.py`

Creare una nuova funzione che segue lo schema gia usato per i pianeti:
1. `_compute_all_segment_jds()` — pre-calcola tutti i JD (fit + verify)
2. Batch eval con i wrapper dei task 2.1-2.5 — una singola chiamata per body
3. Unwrap longitudine sull'array (numpy: `np.unwrap(np.radians(lons))`)
4. `_fit_and_verify_from_values()` — fitta senza ulteriori chiamate

**Attenzione all'unwrapping:** oggi `_generate_segments_unwrap()` fa l'unwrap
per-segmento. Con il path vettorizzato, l'unwrap va fatto sull'array completo
di tutti i JD del segmento corrente, non globalmente. Riusare la logica di
`_fit_and_verify_from_values()` ma aggiungendo l'unwrap per ogni blocco di
`degree + 1` fit values.

**Riferimento:** `generate_body_icrs()` (linea 685) per il pattern completo.

#### 2.7 Aggiornare assemble_leb() per usare il path vettorizzato

**File:** `scripts/generate_leb.py`

Nel blocco "1c. Generate analytical bodies" (linea ~1454), sostituire il loop
sequenziale che chiama `generate_single_body()` con il nuovo
`generate_body_ecliptic_vectorized()` per i 6 body eclittici.
I 9 Uraniani restano sequenziali (gia veloci, pura aritmetica).

#### 2.8 Benchmark e validazione

Dopo le ottimizzazioni, verificare:
1. I tempi di generazione per ogni gruppo
2. Gli errori di fitting sono identici (o migliori) rispetto al path scalare
3. I test esistenti passano (`pytest tests/test_leb/ -v -m "not slow"`)

**Tempi attesi dopo ottimizzazione:**

| Gruppo | Prima | Dopo | Speedup |
|--------|-------|------|---------|
| Pianeti | ~15s | ~15s | — |
| Asteroidi | ~3-5 min | ~3-5 min | — (limite spktype21) |
| Eclittici (6 body) | ~5-8 min | **~5-10s** | **50-100x** |
| Uraniani (9 body) | ~10s | ~10s | — |
| **Totale (base tier)** | **~10-15 min** | **~3-5 min** | **~3-5x** |

Il bottleneck residuo dopo queste ottimizzazioni saranno gli asteroidi,
limitati dalla natura scalare di `spktype21.compute_type21()`.

---

## 3. Architettura: Disaccoppiamento modalita SPK+Kepler vs LEB

### Contesto

Oggi la libreria ha due modalita di funzionamento:

1. **Modalita Skyfield** (default): ogni chiamata a `swe_calc_ut()` interroga
   Skyfield in tempo reale (SPK + frame rotation + nutation). Per gli asteroidi
   senza SPK, fallback a Keplerian. Precisa ma lenta (~120us/chiamata).

2. **Modalita LEB**: se un file `.leb` e attivato (via `set_leb_file()` o
   `LIBEPHEMERIS_LEB`), le chiamate vengono servite dai coefficienti Chebyshev
   precompilati (~2us/chiamata, ~14x speedup). Fallback a Skyfield per body
   non nel LEB.

L'attivazione avviene in modo esplicito tramite:
```python
from libephemeris import set_leb_file
set_leb_file("/path/to/ephemeris.leb")
```
oppure:
```bash
export LIBEPHEMERIS_LEB=/path/to/ephemeris.leb
```

### Obiettivo

Rendere completamente trasparente l'uso della libreria senza LEB. Un utente
deve poter usare `libephemeris` installando solo il pacchetto e senza mai
toccare file `.leb`, con la certezza che:
- Tutti i pianeti funzionano (via Skyfield/DE440)
- Gli asteroidi funzionano (via auto-download SPK + fallback Kepler)
- Non ci sono errori o warning relativi a LEB

### Task

#### 3.1 Verificare che la modalita senza LEB sia completa

**File:** `libephemeris/planets.py`, `libephemeris/state.py`

Verificare che TUTTE le funzioni pubbliche funzionino correttamente quando
nessun file LEB e caricato. In particolare:
- `swe_calc_ut()` / `swe_calc()` — dispatch a Skyfield, non a LEB
- `swe_houses()` — calcolo case senza LEB
- `swe_fixstar_ut()` — stelle fisse senza LEB
- Tutte le ayanamsa (Lahiri, True Citra, ecc.)
- `swe_deltat()` — delta-T senza LEB (da tabella IERS)

Creare test specifici che eseguono le funzioni principali **senza** LEB
attivato e verificano che i risultati siano corretti.

#### 3.2 Documentare le due modalita

**File:** `docs/LEB_GUIDE.md`, `README.md`

Creare una sezione chiara che spieghi:

**Modalita 1: Skyfield (default, zero configurazione)**
```bash
pip install libephemeris
```
```python
import libephemeris as eph
result, flags = eph.swe_calc_ut(2451545.0, 0, 0)  # Funziona subito
```
- Pro: nessun file da gestire, sempre aggiornato
- Contro: piu lento (~120us/chiamata)

**Modalita 2: LEB (massima velocita)**
```bash
pip install libephemeris
# Generare o scaricare il file .leb
export LIBEPHEMERIS_LEB=/path/to/ephemeris_base.leb
```
```python
import libephemeris as eph
result, flags = eph.swe_calc_ut(2451545.0, 0, 0)  # ~14x piu veloce
```
- Pro: ~14x speedup
- Contro: richiede file .leb precompilato (~50MB per base tier)

#### 3.3 Variabile d'ambiente per modalita esplicita

**File:** `libephemeris/state.py`

Considerare l'aggiunta di `LIBEPHEMERIS_MODE` come variabile d'ambiente:

```bash
export LIBEPHEMERIS_MODE=skyfield  # Forza modalita Skyfield (ignora LEB anche se presente)
export LIBEPHEMERIS_MODE=leb       # Forza modalita LEB (errore se file non trovato)
export LIBEPHEMERIS_MODE=auto      # Default: usa LEB se disponibile, altrimenti Skyfield
```

Questo e opzionale — il comportamento attuale (auto-detection basata su
`LIBEPHEMERIS_LEB`) e gia ragionevole. Valutare se la complessita aggiuntiva
vale la pena.

#### 3.4 Gestione graceful di LEB mancante

**File:** `libephemeris/state.py`

Oggi `get_leb_reader()` (linea 191) logga un warning se il file LEB non
esiste e ritorna `None`. Verificare che questo sia il comportamento in tutti
i casi edge:
- File LEB specificato ma non esistente → warning + fallback Skyfield
- File LEB corrotto → warning + fallback Skyfield
- File LEB con range troppo stretto per la data richiesta → fallback Skyfield per quel body
- Nessun file LEB specificato → modalita Skyfield pura, zero warning

#### 3.5 Distribuzione file LEB pre-generati

Valutare se distribuire file `.leb` pre-generati:
- Come asset di GitHub release
- Come pacchetto PyPI separato (`libephemeris-data`)
- Come download on-demand (simile a come Skyfield scarica i file DE)

Questo e a bassa priorita ma migliorerebbe l'esperienza utente per chi
vuole la modalita LEB senza dover generare il file localmente.

---

## 4. Precisione: Miglioramento fallback Kepleriano

### Contesto

Il fallback Kepleriano attuale (`minor_bodies.py`) usa:
- Elementi orbitali osculanti a un'epoca fissa (JD 2461000.5 = Sep 2025)
- Perturbazioni secolari di primo ordine da Giove, Saturno, Urano, Nettuno
- Solo `omega`, `Omega`, `n` vengono perturbati; `a`, `e`, `i` sono costanti
- Correzione di librazione risonante per plutini (Ixion, Orcus)

Errori tipici:
- Asteroidi fascia principale: 10-30" su mesi, gradi su decenni
- TNO: 1-3' su mesi, gradi su anni
- Causa principale: `a`, `e`, `i` costanti + assenza perturbazioni a breve periodo

Esiste gia nel codice un hook per REBOUND/ASSIST (`rebound_integration.py`),
che e il metodo piu preciso disponibile (sub-arcsecond). Ma richiede
dipendenze opzionali (`pip install rebound assist`) e file dati (~1GB).

La cascata di fallback attuale in `planets.py:1643-1728` e:
```
1. SPK kernel (registrato)            → sub-arcsecond
2. Auto SPK download (da Horizons)    → sub-arcsecond
3. Strict precision check             → errore se attivo
4. REBOUND/ASSIST (se installato)     → sub-arcsecond
5. Keplerian con perturbazioni        → gradi su decenni (ULTIMO RESORT)
```

### Task

#### 4.1 Perturbazioni secolari di secondo ordine

**Priorita:** Alta — miglioramento significativo con complessita bassa.

**File:** `libephemeris/minor_bodies.py`

La funzione `apply_secular_perturbations()` (linea 908) oggi perturba solo
`omega`, `Omega`, `n`. Aggiungere la variazione secolare anche per:
- `a` (semi-asse maggiore) — drift secolare dovuto alle risonanze
- `e` (eccentricita) — variazione forzata da Giove
- `i` (inclinazione) — variazione forzata da Giove

Le formule sono nella teoria di Laplace-Lagrange (Murray & Dermott,
"Solar System Dynamics", Cap. 7). Per il caso piu semplice (singolo
pianeta perturbatore):

```
de/dt = -(n_j * m_j * alpha * b_{3/2}^{(2)}(alpha)) / (4 * a_j) * sin(omega - omega_j)
di/dt = -(n_j * m_j * alpha * b_{3/2}^{(1)}(alpha)) / (4 * a_j) * sin(Omega - Omega_j)
```

dove `b_{s}^{(j)}` sono i coefficienti di Laplace (gia implementati come
`_calc_laplace_coefficients()` a linea 670).

Modificare `calc_secular_perturbation_rates()` (linea 713) per calcolare
anche `d_a`, `d_e`, `d_i` e ritornare 6 rate invece di 3.

Modificare `apply_secular_perturbations()` per applicare:
```python
a_pert = elements.a + d_a * dt
e_pert = elements.e + d_e * dt
i_pert = elements.i + d_i * dt
```

**Guadagno atteso:** errori da gradi a decine di minuti d'arco su secoli.

**Test:** confrontare con le posizioni SPK su intervalli di 10, 50, 100 anni
e verificare che l'errore migliori rispetto al Keplerian attuale.

#### 4.2 Perturbazioni a breve periodo

**Priorita:** Media — miglioramento notevole ma implementazione complessa.

**File:** `libephemeris/minor_bodies.py`

Le perturbazioni a breve periodo sono oscillazioni nell'orbita con periodo
comparabile al periodo orbitale del pianeta perturbatore. Per Giove (periodo
~12 anni), queste oscillazioni valgono diversi minuti d'arco e si ripetono
ciclicamente.

Le formule classiche (Brouwer & Clemence) danno le perturbazioni come serie
trigonometriche nell'anomalia media dell'asteroide e del pianeta:

```
delta_lon = sum_j sum_k C_{j,k} * cos(j*M + k*M_J + phi_{j,k})
```

Per i 5 asteroidi principali, i termini dominanti sono ~10-20 per Giove +
~5-10 per Saturno.

Implementazione:
1. Aggiungere le frequenze `M_J(t)`, `M_S(t)` dei pianeti perturbatori
   (anomalie medie di Giove e Saturno come funzione del tempo — formule
   semplici, lineari in T)
2. Per ogni asteroide, pre-computare i coefficienti `C_{j,k}` e fasi 
   `phi_{j,k}` (da letteratura o da fit sui dati SPK)
3. Sommare le correzioni alla longitudine, latitudine, distanza

**Guadagno atteso:** errori da minuti d'arco a decine di secondi d'arco su decenni.

**Alternativa:** invece di implementare la teoria completa, fare un fit
empirico dei residui (SPK - Keplerian) con serie di Fourier sui primi N
armonici. Questo cattura automaticamente le perturbazioni dominanti senza
dover derivare le formule analitiche.

#### 4.3 Elementi orbitali multi-epoca

**Priorita:** Media — semplice e molto efficace.

**File:** `libephemeris/minor_bodies.py`

Oggi tutti i 34 asteroidi in `MINOR_BODY_ELEMENTS` (linea 973) hanno elementi
osculanti a un'unica epoca (JD 2461000.5 = Sep 2025). La propagazione parte
sempre da li, accumulando errore nel tempo.

Alternativa: scaricare da JPL SBDB gli elementi osculanti a **piu epoche**
(es. ogni 50 anni dal 1600 al 2500) e scegliere l'epoca piu vicina alla data
richiesta. Questo riduce drasticamente il `dt` di propagazione e quindi l'errore.

```python
MINOR_BODY_ELEMENTS_MULTI = {
    SE_CERES: [
        OrbitalElements(epoch=2305447.5, ...),  # ~1600
        OrbitalElements(epoch=2323711.5, ...),  # ~1650
        OrbitalElements(epoch=2341975.5, ...),  # ~1700
        # ...
        OrbitalElements(epoch=2506331.5, ...),  # ~2150
    ],
}
```

La funzione `calc_minor_body_position()` sceglierebbe l'epoca piu vicina
e propagherebbe da li. Il `dt` massimo passerebbe da ~500 anni a ~25 anni,
riducendo l'errore di ~20x.

**Sorgente dati:** JPL SBDB API (`https://ssd.jpl.nasa.gov/api/sbdb.api`)
permette di richiedere elementi a epoche specifiche.

**Guadagno atteso:** errori da gradi a minuti d'arco su secoli, senza
alcuna teoria delle perturbazioni aggiuntiva.

#### 4.4 Integrazione REBOUND/ASSIST come fallback intermedio

**Priorita:** Alta — il codice esiste gia, va solo attivato e documentato.

**File:** `libephemeris/rebound_integration.py`, `libephemeris/planets.py`

Il modulo `rebound_integration.py` (740 linee) e completo e funzionale:
- `propagate_orbit_assist()` — integrazione ephemeris-quality con ASSIST
- `propagate_orbit_rebound()` — integrazione 2-body con REBOUND
- `propagate_trajectory()` — multi-point con fallback automatico
- `compare_with_keplerian()` — confronto di precisione

Il hook in `planets.py:1677-1710` chiama `check_assist_available()` e se
REBOUND/ASSIST e installato, lo usa prima del Keplerian.

Task:
1. Verificare che il path ASSIST funzioni end-to-end per i 5 asteroidi
   principali su date dentro e fuori il range SPK
2. Aggiungere test automatici
3. Documentare come installare: `pip install libephemeris[nbody]`
4. Verificare la performance (ASSIST e lento per singole valutazioni —
   potrebbe essere piu lento del path Skyfield diretto)
5. Considerare un meccanismo di caching: integrare una volta, salvare i
   risultati, e riusarli per chiamate future (mini-LEB runtime)

**Dipendenze:** `rebound>=4.0.0`, `assist>=1.1.0` (opzionali, gia in
`pyproject.toml` come extra `nbody`)

**File dati richiesti da ASSIST:**
- `linux_p1550p2650.440` — effemeridi planetarie (~100MB)
- `sb441-n16.bsp` — 16 asteroidi massicci (~900MB)

#### 4.5 Caching dei risultati REBOUND per asteroidi

**Priorita:** Bassa — ottimizzazione di performance.

**File:** nuovo modulo o estensione di `rebound_integration.py`

L'integrazione REBOUND/ASSIST e lenta per singole valutazioni (~ms).
Per un uso efficiente come fallback runtime, implementare un meccanismo di
caching:
1. Al primo accesso a un asteroide senza SPK, integrare l'orbita su una
   griglia di punti (es. ogni 10 giorni per 100 anni)
2. Salvare i risultati in un file cache locale
3. Per richieste successive, interpolare dalla griglia (Chebyshev o spline)

Questo e in pratica un "mini-LEB" generato on-demand per un singolo asteroide.

#### 4.6 Validazione sistematica della precisione Kepleriana

**Priorita:** Alta — prerequisito per valutare i miglioramenti.

**File:** nuovo test in `tests/`

Creare un test di benchmark che confronta il Keplerian con le posizioni SPK
per i 5 asteroidi principali a intervalli di:
- 1 mese, 1 anno, 5 anni, 10 anni, 50 anni, 100 anni dall'epoca

Registrare gli errori in longitudine, latitudine, distanza. Questo serve come
baseline per misurare i miglioramenti dei task 4.1-4.3.

Usare le posizioni SPK come "verita" (scaricate via `auto_download_asteroid_spk()`).

---

## Note

### Ordine di implementazione raccomandato

1. **Task 1.5** — Abortire se SPK manca (prerequisito per tutto)
2. **Task 1.6** — Migliorare verify_leb() (prerequisito per validazione)
3. **Task 1.1-1.4** — Range per-body per asteroidi nel LEB
4. **Task 2.1-2.7** — Vettorizzazione generatore (speedup eclittici)
5. **Task 3.1-3.4** — Disaccoppiamento modalita
6. **Task 4.6** — Benchmark Keplerian (baseline)
7. **Task 4.1** — Perturbazioni secolari secondo ordine
8. **Task 4.3** — Elementi multi-epoca
9. **Task 4.4** — Attivazione REBOUND/ASSIST
10. **Task 4.2** — Perturbazioni a breve periodo (piu complesso)
11. **Task 4.5** — Caching REBOUND

### File principali coinvolti

| File | Task |
|------|------|
| `scripts/generate_leb.py` | 1.1-1.7, 2.1-2.8 |
| `libephemeris/leb_reader.py` | 1.3 |
| `libephemeris/fast_calc.py` | 1.4 |
| `libephemeris/planets.py` | 1.4, 3.1, 4.4 |
| `libephemeris/state.py` | 3.3, 3.4 |
| `libephemeris/minor_bodies.py` | 1.1, 4.1, 4.2, 4.3, 4.6 |
| `libephemeris/rebound_integration.py` | 4.4, 4.5 |
| `libephemeris/spk.py` | 1.1 |
| `libephemeris/lunar.py` | 2.1-2.5 (riferimento, non modificato) |
| `docs/LEB_GUIDE.md` | 1.8, 3.2 |
| `pyproject.toml` | 1.7 |
| `tests/` | 2.8, 3.1, 4.6 |
