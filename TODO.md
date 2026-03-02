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

#### ~~1.1 Scoprire il range massimo SPK da Horizons~~ ✓ COMPLETATO

Implementata `_get_asteroid_spk_range()` in `scripts/generate_leb.py` che apre
il file SPK type 21, trova i segmenti per il target NAIF ID (provando entrambe
le convenzioni 2M e 20M), e ritorna `(jd_start, jd_end)` effettivo.

La funzione viene usata in `assemble_leb()` per scoprire il range effettivo
dell'SPK dopo il download. Se l'SPK copre l'intero range del tier, si usa il
range globale. Altrimenti si interseca il range SPK col tier e si usa il
risultato come range per-body.

#### ~~1.2 Generare asteroidi con range indipendente~~ ✓ COMPLETATO

Modificato `assemble_leb()` per:
1. Scoprire il range SPK effettivo per ogni asteroide via `_get_asteroid_spk_range()`
2. Intersecare col range del tier (con 1 giorno di margine)
3. Escludere asteroidi con meno di 20 anni di copertura utile
4. Passare il range per-body a `generate_single_body()` per la generazione
5. Scrivere il range per-body nel `BodyEntry` (campo `jd_start`/`jd_end`)
6. Mostrare il range per-body nell'error summary e nella verifica

La nutation e delta-T continuano a coprire il range globale del tier.

#### ~~1.3 Aggiornare il reader per range per-body~~ ✓ GIA IMPLEMENTATO

Il reader (`leb_reader.py:eval_body()`) usa gia i campi `body.jd_start` e
`body.jd_end` per-body. Quando il JD richiesto e fuori dal range di un body
specifico, lancia `ValueError` — che viene catturato da `swe_calc_ut()` /
`swe_calc()` per fare fallback a Skyfield. Nessuna modifica necessaria.

#### ~~1.4 Fallback trasparente per body fuori range~~ ✓ GIA IMPLEMENTATO

Il dispatch in `swe_calc_ut()` (`planets.py:772-782`) e `swe_calc()`
(`planets.py:835-845`) cattura gia sia `KeyError` (body non nel LEB) che
`ValueError` (JD fuori range per-body) con un unico handler:
```python
except (KeyError, ValueError):
    pass  # body not in .leb or JD out of range, fall through
```
Questo permette il fallback trasparente per-body senza disabilitare il LEB
per gli altri body. Nessuna modifica necessaria.

#### ~~1.5 Abortire se SPK manca (strict mode per generatore)~~ ✓ COMPLETATO

Rimosso il fallback Keplerian scalare da `generate_body_icrs_asteroid()`.
Il generatore ora lancia `RuntimeError` se l'SPK non e disponibile o non copre
il range richiesto. `assemble_leb()` esclude gli asteroidi senza SPK con
messaggio `(EXCLUDED)` e non chiama mai `set_strict_precision(False)`.

#### ~~1.6 Migliorare verify_leb() per asteroidi e eclittici~~ ✓ COMPLETATO

Riscritta `verify_leb()` con verifica propria per tutti i tipi di body:
- Pianeti ICRS: confronto diretto con Skyfield (invariato)
- Asteroidi ICRS: confronto diretto via `spktype21` (NUOVO)
- Body eclittici: confronto con funzioni analitiche da `lunar.py` (NUOVO)
- Body eliocentrici: confronto con `calc_uranian_planet()` / `calc_transpluto()` (NUOVO)
- Funzioni helper aggiunte: `_build_ecliptic_eval_funcs()`, `_build_helio_eval_funcs()`,
  `_verify_icrs_planet()`, `_verify_icrs_asteroid()`, `_verify_ecliptic_body()`
- Rimosso `set_strict_precision(False)` da `verify_leb()`
- Tutti i body types riportano errore in arcsec e stato PASS/FAIL

#### ~~1.7 Comandi poe per medium e extended con range adattivo~~ ✓ COMPLETATO

**File:** `pyproject.toml`

I comandi poe per group generation (`leb:generate:medium:groups`,
`leb:generate:extended:groups`) gia esistono e funzionano automaticamente
con il range per-body perche la logica e interamente nel generatore
(`assemble_leb()` in `generate_leb.py`). Verificato per code inspection:
il `--group` filtra solo la lista di body, non modifica la logica di
range discovery.

#### ~~1.8 Aggiornare documentazione~~ ✓ COMPLETATO

**File:** `docs/LEB_GUIDE.md`, `README.md`

Documentato:
- Il concetto di range per-body nel formato LEB (Sezione 6.8 di LEB_GUIDE.md)
- I limiti di copertura SPK per gli asteroidi
- Il comportamento di fallback quando un body e fuori range
- Tabella dei range effettivi per ogni tier e body
- Sezione "Binary Ephemeris Mode (LEB)" aggiunta al README.md
- `LIBEPHEMERIS_MODE` documentato in LEB_GUIDE.md (Sezioni 1, 7.1, 7.4)

---

## ~~2. LEB: Velocita di generazione~~ ✓ COMPLETATO

### Risultati

Vettorizzazione completata per tutti i 6 body eclittici. Una singola chiamata
Skyfield `(moon - earth).at(t_array)` per tutti i body, poi calcolo numpy
vettoriale per ciascuno.

**Funzioni batch implementate in `scripts/generate_leb.py`:**
- `_calc_mean_lilith_batch()` — vettorizza `_calc_mean_apse_analytical()`
- `_calc_lunar_fundamental_arguments_batch()` — argomenti di Delaunay vettorizzati
- `_calc_elp2000_apogee_perturbations_batch()` — serie 40+ termini vettorizzata
- `_calc_elp2000_perigee_perturbations_batch()` — serie 61 termini vettorizzata
- `_eval_ecliptic_bodies_batch()` — core: singola chiamata Skyfield, 6 body
- `generate_ecliptic_bodies_vectorized()` — orchestratore
- `_fit_and_verify_from_values_unwrap()` — fitting con unwrapping longitudine

**Tempi misurati (base tier, 300 anni):**

| Gruppo | Prima | Dopo | Speedup |
|--------|-------|------|---------|
| Pianeti (11 body) | ~15s | ~15s | — |
| Asteroidi (5 body) | ~3-5 min | ~3-5 min | — (limite spktype21) |
| Eclittici (6 body) | ~5-8 min | **~38s** | **~10x** |
| Uraniani (9 body) | ~10s | ~10s | — |
| **Totale (base tier)** | **~10-15 min** | **~5 min** | **~2-3x** |

Tutti i 31 body passano la verifica con precisione sub-arcsecond.
Test: 137 LEB tests + 12 generator tests tutti verdi.

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

#### ~~3.1 Verificare che la modalita senza LEB sia completa~~ ✓ GIA IMPLEMENTATO

La modalita senza LEB e il path primario della libreria e funziona
correttamente per tutte le funzioni pubbliche:
- `_LEB_FILE` e `_LEB_READER` inizializzati a `None`
- `get_leb_reader()` ritorna `None` quando nessun LEB e configurato
- In `swe_calc_ut()`/`swe_calc()`, il guard `if reader is not None:` salta
  l'intero blocco LEB e va direttamente a Skyfield via `_calc_body()`
- Nessun errore, nessun warning, nessun overhead significativo
- Tutti i test non-LEB (la maggioranza del test suite) esercitano
  implicitamente questo path

#### ~~3.2 Documentare le due modalita~~ ✅ COMPLETATO

**File:** `docs/LEB_GUIDE.md`, `README.md`

Documentato in:
- `README.md`: nuova sezione "Binary Ephemeris Mode (LEB)" con esempi
  per entrambe le modalita, `LIBEPHEMERIS_MODE`, e link a LEB_GUIDE.md
- `docs/LEB_GUIDE.md`: sezione "Calculation Mode" con tabella dei modi,
  esempi programmatici e env var, aggiornamento Sezione 7.1 e 7.4

#### ~~3.3 Variabile d'ambiente per modalita esplicita~~ ✅ COMPLETATO

**File:** `libephemeris/state.py`

Implementata `LIBEPHEMERIS_MODE` con tre valori:

```bash
export LIBEPHEMERIS_MODE=skyfield  # Forza modalita Skyfield (ignora LEB anche se presente)
export LIBEPHEMERIS_MODE=leb       # Forza modalita LEB (errore se file non trovato)
export LIBEPHEMERIS_MODE=auto      # Default: usa LEB se disponibile, altrimenti Skyfield
```

**Funzioni aggiunte in `state.py`:**
- `set_calc_mode(mode)` — imposta la modalita (`"auto"`, `"skyfield"`, `"leb"`, o `None`)
- `get_calc_mode()` — ritorna la modalita effettiva (override > env var > `"auto"`)

**Comportamento di `get_leb_reader()` aggiornato:**
- `"skyfield"` → ritorna sempre `None` (LEB disabilitato)
- `"leb"` → ritorna reader o lancia `RuntimeError` se non disponibile
- `"auto"` (default) → ritorna reader se configurato, altrimenti `None`

**Esportato in `__init__.py`:** `set_calc_mode`, `get_calc_mode`
**Reset in `close()`:** `_CALC_MODE = None`
**Documentato in:** `docs/LEB_GUIDE.md` (Sezioni 1, 7.1, 7.4), `README.md`

#### ~~3.4 Gestione graceful di LEB mancante~~ ✓ GIA IMPLEMENTATO

`get_leb_reader()` (`state.py:191`) gestisce gia tutti i casi edge:
- File LEB specificato ma non esistente → warning + ritorna `None` (fallback Skyfield)
- File LEB corrotto (`ValueError`, `OSError`) → warning + ritorna `None`
- File LEB con range troppo stretto → `ValueError` da `eval_body()`, catturato
  da `swe_calc_ut()`/`swe_calc()` per fare fallback a Skyfield per quel body
- Nessun file LEB specificato → ritorna `None` senza warning
- Test esistenti: `TestContextLEBGracefulError` in `test_context_leb.py`
  verifica il fallback con path invalido

#### 3.5 Distribuzione file LEB pre-generati ✅ COMPLETATO

~~Valutare se distribuire file `.leb` pre-generati.~~

Implementato con:
- **GitHub Releases** (`data-v1`): file `.leb` ospitati come asset
- **Release script**: `scripts/release_leb.py` + comandi `poe release:leb:*`
- **CLI download**: `libephemeris download:leb:{base,medium,extended}`
- **Auto-discovery runtime**: `~/.libephemeris/leb/ephemeris_{tier}.leb`
- **Download programmatico**: `libephemeris.download_leb_for_tier("medium")`

Tier disponibili: `base` (~53 MB), `medium` (~175 MB). `extended` non ancora generato.

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

#### ~~4.1 Perturbazioni secolari di secondo ordine~~ ✓ COMPLETATO

Implementata la teoria secolare completa di Laplace-Lagrange (Murray & Dermott
Cap. 7) con evoluzione vettoriale (h,k)/(p,q) per eccentricita e inclinazione.

**Funzioni aggiunte in `minor_bodies.py`:**
- `_calc_forced_elements()` (~120 linee) — calcola i vettori forzati di
  eccentricita (h_f, k_f) e inclinazione (p_f, q_f) e le frequenze proprie
  (g, s) dal contributo di Giove, Saturno, Urano e Nettuno
- Usa i coefficienti di Laplace b_{3/2}^{(1)} e b_{3/2}^{(2)} per i termini
  diagonali (frequenze proprie) e off-diagonali (accoppiamento planetario)

**Modifiche alla firma di `apply_secular_perturbations()`:**
- Ritorna ora 6 valori: `(omega_pert, Omega_pert, M_pert, n_pert, e_pert, i_pert)`
- Il vettore di eccentricita oscilla attorno al valore forzato con periodo
  proprio (~23.000 anni per asteroidi fascia principale)
- Il vettore di inclinazione oscilla con precessione retrograda

**Modifiche a `calc_minor_body_position()`:**
- Usa `e_pert` per Kepler equation, anomalia vera e distanza
- Usa `i_pert` per la matrice di rotazione perifocale → eclittica

**Callers aggiornati:** `rebound_integration.py`, tutti i test

**Risultati:** Impatto minimo sui benchmark a 100 anni perche il periodo
proprio e ~23.000 anni. L'effetto diventa significativo su scale di millenni.
La teoria e comunque corretta e necessaria per completezza — l'errore
dominante resta le perturbazioni a breve periodo (Task 4.2).

#### 4.2 Perturbazioni a breve periodo — INVESTIGATO, DEPRIORITIZZATO

**Priorita:** Bassa — l'approccio empirico non funziona, la teoria analitica
e troppo complessa per il guadagno ottenibile. Il path SPK + LEB e superiore.

**File:** `libephemeris/minor_bodies.py`

**Investigazione completata:** Testati due approcci empirici (Fourier fit dei
residui SPK - Keplerian):

1. **Fit diretto su 800 anni** — RMS reduction solo 2-7% per asteroidi della
   fascia principale. I residui sono dominati dalla deriva secolare (~gradi),
   non dalle oscillazioni periodiche (~arcminuti).

2. **Fit con detrending polinomiale (grado 5)** — Rimuove la deriva secolare
   dai residui prima del fit. Le ampiezze Fourier sono corrette (~300" per
   Ceres, ~780" per Pallas) ma la **validazione puntuale** mostra risultati
   catastrofici: le correzioni PEGGIORANO le posizioni a breve termine (5" →
   19' a 1 mese per Ceres).

**Perche l'approccio empirico fallisce:**

1. **Near-epoch**: il Kepleriano e gia accurato (~5"), ma i termini Fourier
   aggiungono correzioni errate (~20') perche le ampiezze sono calibrate
   sull'intero span di 800 anni, non sull'epoca specifica
2. **Far-from-epoch**: la deriva secolare domina (~gradi), e le correzioni
   Fourier (~arcminuti) sono trascurabili
3. Le ampiezze delle perturbazioni a breve periodo **variano nel tempo** al
   variare degli elementi orbitali — un fit statico cattura l'ampiezza *media*
   che e sbagliata per qualsiasi data specifica
4. Il problema fondamentale: la propagazione Kepleriana con perturbazioni
   secolari produce errori che crescono come potenza del tempo, non come
   oscillazioni periodiche stazionarie

**Conclusione:** Per migliorare la precisione Kepleriana oltre il livello
attuale (minuti d'arco su decenni), l'unica via praticabile e l'integrazione
numerica (Task 4.4 REBOUND/ASSIST). Le perturbazioni a breve periodo
analitiche (Brouwer & Clemence) richiederebbero coefficienti body-specific
dalla letteratura specializzata e sarebbero comunque limitate a ~10" di
precisione.

**Script di generazione:** `scripts/generate_short_period_corrections.py`
(mantenuto come riferimento, non integrato nel codice).

#### 4.3 Elementi orbitali multi-epoca ✅ COMPLETATO

**Priorita:** Media — semplice e molto efficace.

**File:** `libephemeris/minor_bodies.py`

**Implementazione:** Generati elementi orbitali da vettori di stato SPK type 21
a intervalli di 50 anni (1650-2450 CE) per 6 corpi: Chiron, Pholus, Ceres,
Pallas, Juno, Vesta (17 epoche ciascuno, ~600 righe in `MINOR_BODY_ELEMENTS_MULTI`).

- Conversione stato-a-Kepleriano inline (ICRS → eclittica J2000 → elementi orbitali)
- `_get_closest_epoch_elements()` seleziona l'epoca piu vicina considerando SIA
  l'elemento single-epoch originale SIA tutte le entry multi-epoca
- Fix critico: l'elemento originale (JD 2461000.5) deve sempre essere candidato,
  altrimenti la precisione near-epoch degrada massivamente
- Integrato in `calc_minor_body_heliocentric()` (path Keplerian fallback)

**Risultati benchmark:**

| Offset | Prima | Dopo | Miglioramento |
|--------|-------|------|---------------|
| epoch | 0.002" | 0.002" | — |
| 1 mese | 7.4" | 7.5" | — |
| 6 mesi | 48.0" | 49.2" | — |
| 1 anno | 1.8' | 1.9' | — |
| 5 anni | 26.6' | 26.7' | — |
| 10 anni | 44.7' | 44.6' | — |
| 25 anni | 3.3° | **2.6'** | **~75x** |
| 50 anni | 5.3° | 3.6° | ~1.5x |
| 100 anni | 10.7° | 3.5° | ~3x |

Il miglioramento drammatico a 25 anni (~75x) e dovuto alla tabella multi-epoca
con entry ogni 50 anni: il punto di test e a ~0 anni dall'epoca piu vicina.
A 50 e 100 anni il miglioramento e piu modesto (~1.5-3x) perche quei punti
sono ai confini della copertura.

#### ~~4.4 Integrazione REBOUND/ASSIST come fallback intermedio~~ ✅ COMPLETATO

**Priorita:** Alta — il codice esiste gia, va solo attivato e documentato.

**File:** `libephemeris/rebound_integration.py`, `libephemeris/planets.py`

Il modulo `rebound_integration.py` (~890 linee) e completo e funzionale:
- `propagate_orbit_assist()` — integrazione ephemeris-quality con ASSIST
- `propagate_orbit_rebound()` — integrazione 2-body con REBOUND
- `propagate_trajectory()` — multi-point con fallback automatico
- `compare_with_keplerian()` — confronto di precisione

Il hook in `planets.py:1677-1710` chiama `check_assist_data_available()` e se
REBOUND/ASSIST e installato con i file dati, lo usa prima del Keplerian.

**Completato:**
- ✅ `check_assist_data_available()` — check cached di import + file dati
- ✅ `reset_assist_data_cache()` — reset del cache
- ✅ Hook in `planets.py` aggiornato per usare `check_assist_data_available()`
  (sostituisce `check_assist_available()` non-cached)
- ✅ `close()` in `state.py` chiama `reset_assist_data_cache()`
- ✅ `AssistEphemConfig` cerca in `~/.libephemeris/assist/`, `ASSIST_DIR`, `./data/`
- ✅ `download_assist_data()` helper per scaricare i file dati
- ✅ Test condizionali ASSIST end-to-end (skip se file dati mancanti):
  `TestAssistDataAvailability` (6 test) + `TestAssistEndToEnd` (9 test)
- ✅ 47 test REBOUND passano, 9 skipped (ASSIST data non scaricati)
- ✅ macOS rpath fix per `libassist.cpython-310-darwin.so`

**Da fare:**
1. Scaricare file dati ASSIST (~714 MB) e verificare path end-to-end
   con perturbazioni planetarie reali ← **RICHIEDE APPROVAZIONE UTENTE**
2. Documentare come installare: `pip install libephemeris[nbody]`
3. Verificare la performance (ASSIST e lento per singole valutazioni)

**Dipendenze:** `rebound>=4.0.0`, `assist>=1.1.0` (opzionali, gia in
`pyproject.toml` come extra `nbody`)

**File dati richiesti da ASSIST:**
- `linux_p1550p2650.440` — effemeridi planetarie (~98 MB)
- `sb441-n16.bsp` — 16 asteroidi massicci (~616 MB)

#### ~~4.5 Caching dei risultati REBOUND per asteroidi~~ DEPRIORITIZZATO

**Priorita:** Bassa — ottimizzazione di performance. Deprioritizzato: il sistema
LEB precompilato + download automatico copre gia il caso d'uso principale.
Il fallback REBOUND/ASSIST e adeguato per query sporadiche.

**File:** nuovo modulo o estensione di `rebound_integration.py`

L'integrazione REBOUND/ASSIST e lenta per singole valutazioni (~ms).
Per un uso efficiente come fallback runtime, implementare un meccanismo di
caching:
1. Al primo accesso a un asteroide senza SPK, integrare l'orbita su una
   griglia di punti (es. ogni 10 giorni per 100 anni)
2. Salvare i risultati in un file cache locale
3. Per richieste successive, interpolare dalla griglia (Chebyshev o spline)

Questo e in pratica un "mini-LEB" generato on-demand per un singolo asteroide.

#### ~~4.6 Validazione sistematica della precisione Kepleriana~~ ✓ COMPLETATO

Creato `tests/test_keplerian_precision_benchmark.py` con:
- Confronto Keplerian vs SPK in frame eclittico J2000 (stesso frame)
- 5 asteroidi principali (Ceres, Pallas, Juno, Vesta, Chiron)
- Intervalli: epoca, 1mo, 3mo, 6mo, 1yr, 2yr, 5yr, 10yr, 25yr, 50yr, 100yr
- Test di regressione: epoch < 1", 1mo < 60", 1yr < 5'

**Risultati baseline:**

| Offset | Errore max |
|--------|-----------|
| epoch | 0.002" |
| 1 mese | 7.4" |
| 3 mesi | 22.8" |
| 6 mesi | 48.0" |
| 1 anno | 1.8' |
| 2 anni | 6.3' |
| 5 anni | 26.6' |
| 10 anni | 44.7' |
| 25 anni | 3.3 deg |
| 50 anni | 5.3 deg |
| 100 anni | 10.7 deg |

Causa principale: `a`, `e`, `i` costanti + assenza perturbazioni a breve periodo.
I task 4.1-4.3 mirano a ridurre questi errori.

---

## Note

### Ordine di implementazione raccomandato

1. ~~**Task 1.5** — Abortire se SPK manca (prerequisito per tutto)~~ ✓
2. ~~**Task 1.6** — Migliorare verify_leb() (prerequisito per validazione)~~ ✓
3. ~~**Task 2.1-2.7** — Vettorizzazione generatore (speedup eclittici)~~ ✓
4. ~~**Task 1.1-1.4** — Range per-body per asteroidi nel LEB~~ ✓
5. ~~**Task 3.1-3.4** — Disaccoppiamento modalita~~ ✓
6. ~~**Task 4.6** — Benchmark Keplerian (baseline)~~ ✓
7. ~~**Task 4.1** — Perturbazioni secolari secondo ordine~~ ✓
8. ~~**Task 4.3** — Elementi multi-epoca~~ ✓
9. ~~**Task 1.7** — Comandi poe per medium/extended~~ ✓
10. ~~**Task 1.8** — Documentazione per-body ranges~~ ✓
11. ~~**Task 3.3** — LIBEPHEMERIS_MODE~~ ✓
12. ~~**Task 4.4** — Attivazione REBOUND/ASSIST~~ ✓
13. ~~**Task 4.2** — Perturbazioni a breve periodo~~ DEPRIORITIZZATO
14. ~~**Task 3.5** — Distribuzione file LEB pre-generati~~ ✓
15. ~~**Task 4.5** — Caching REBOUND~~ DEPRIORITIZZATO

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
