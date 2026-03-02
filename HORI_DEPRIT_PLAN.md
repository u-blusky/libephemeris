# Hori-Deprit Perturbation Theory — Implementation Plan

## Obiettivo

Implementare la teoria di perturbazione di Hori-Deprit (trasformazione di Lie) per
il fallback Kepleriano degli asteroidi in `libephemeris/minor_bodies.py`.

Precisione target: da ~15-30" (stato attuale) a ~1-5" a 1 anno di propagazione.

## Stato attuale del pipeline perturbativo

Il pipeline in `minor_bodies.py` applica attualmente:

| Componente | Funzione | Cosa fa |
|------------|----------|---------|
| Perturbazioni secolari lineari | `calc_secular_perturbation_rates()` L897 | dω/dt, dΩ/dt, dn via Laplace-Lagrange 1° ordine |
| Eccentricità/inclinazione forzate | `_calc_forced_elements()` L1150 | Vettori (h,k)/(p,q) con oscillazione libera + forzata |
| Elementi planetari evoluti | `_get_planet_elements_at_time()` L176 | Tassi lineari Simon et al. (1994) |
| Coefficienti di Laplace | `_calc_laplace_coefficients()` L851 | Integrazione numerica b_s^(j)(α) |
| Correzione short-period (M1) | `_calc_short_period_correction()` L3597 | Termini sinodici 1° ordine da Giove/Saturno |
| Amplificazione near-resonance | in `calc_secular_perturbation_rates()` | Fattore ≤3x per ω/Ω near resonance |
| Orchestrazione | `apply_secular_perturbations()` L1379 | Combina tutto, restituisce elementi perturbati |
| Posizione finale | `calc_minor_body_position()` L3748 | Kepler + perturbazioni → (x,y,z) eclittica |

### Cosa manca (gap coperto da Hori-Deprit)

1. **Termini a lungo periodo**: oscillazioni con periodo ∝ 1/e (decenni-secoli), non
   catturate dalla teoria secolare lineare
2. **Accoppiamento secolare-periodico**: l'ampiezza dei termini short-period varia col
   tempo man mano che ω e Ω precèdono — M1 usa ampiezze istantanee ma non il loro
   drift sistematico
3. **Termini di 2° ordine nella funzione disturbante**: le correzioni O(μ²) alle
   frequenze proprie e alle ampiezze forzate
4. **Correzioni short-period a tutti e 6 gli elementi orbitali**: M1 corregge solo la
   longitudine media, non (a, e, i, ω, Ω) separatamente

## Cos'è la trasformazione di Hori-Deprit

La teoria di Deprit (1969) usa una **trasformazione canonica generata da una serie di Lie**
per separare il moto in componenti a breve periodo, lungo periodo e secolare.

Data l'Hamiltoniana H = H₀ + εH₁ + ε²H₂ + ..., la trasformazione produce:

- **Elementi medi** (ā, ē, ī, ω̄, Ω̄, M̄) che evolvono solo secolarmente e a lungo periodo
- **Funzioni generatrici** W₁, W₂ che danno le correzioni periodiche:
  - `elemento_osculante = elemento_medio + ∂W₁/∂coniugato + ½ ∂²W₂/∂... + ...`

Il vantaggio rispetto a von Zeipel (Brouwer): la trasformazione è **ricorsiva** — l'ordine
N si calcola meccanicamente dall'ordine N-1, senza riderivare tutto.

## Architettura dell'implementazione

### Strategia: sostituzione incrementale

Non riscriviamo tutto. Manteniamo il pipeline esistente e **sostituiamo/estendiamo** i
pezzi rilevanti:

```
PRIMA (attuale):
  elementi osculanti
  → apply_secular_perturbations()     [secolare 1° ordine, vettori (h,k)/(p,q)]
  → calc_minor_body_position()        [Kepler + _calc_short_period_correction()]
  → (x, y, z)

DOPO (Hori-Deprit):
  elementi osculanti
  → hori_deprit_mean_elements()       [NUOVO: osculanti → medi via W₁]
  → propagate_mean_elements()         [MODIFICATO: evoluzione secolare + lungo periodo]
  → hori_deprit_osculating()          [NUOVO: medi → osculanti via W₁ inversa]
  → calc_minor_body_position()        [Kepler standard, NO _calc_short_period_correction()]
  → (x, y, z)
```

La chiave: i termini short-period sono ora **dentro la trasformazione** invece che
applicati ad hoc dopo.

### File da modificare

- `libephemeris/minor_bodies.py` — unico file (tutte le funzioni perturbative sono qui)

### File da NON modificare

- `constants.py` — nessun nuovo dato necessario
- `tests/test_keplerian_precision_benchmark.py` — già espanso a 37 corpi
- API pubblica — `calc_minor_body_position()` mantiene la stessa firma

## Piano dettagliato per fase

### Fase 1: Espansione della funzione disturbante (fondamenta)

**Cosa**: Calcolo sistematico dei coefficienti della funzione disturbante R.

**Funzioni da creare**:

```python
def _disturbing_function_coefficients(
    alpha: float, e: float, e_p: float, i: float, i_p: float
) -> DisturbingFunctionCoeffs:
    """Espansione della funzione disturbante R fino a 2° ordine in e, i.

    Usa i coefficienti di Laplace b_{s}^{(j)}(α) già implementati in
    _calc_laplace_coefficients() per calcolare i coefficienti della
    serie di Fourier:

    R = Σ C_{j,k,l,m,n} · cos(j·λ + k·λ' + l·ϖ + m·ϖ' + n·Ω + ...)

    Dove includiamo termini fino a:
    - Ordine 0 in e: termine secolare diretto
    - Ordine 1 in e: termini e·cos(ϖ - ϖ'), e·cos(λ - λ'), etc.
    - Ordine 2 in e: termini e²·cos(2ϖ - 2ϖ'), e·e'·cos(...), etc.
    - Ordine 2 in sin(i/2): termini analoghi per l'inclinazione

    Coefficienti dal Cap. 6 di Murray & Dermott (1999), Tabella 6.1-6.4.
    """
```

**Dipendenze**: `_calc_laplace_coefficients()` (L851, già implementata)

**Stima**: ~150 righe

### Fase 2: Funzione generatrice W₁ (1° ordine)

**Cosa**: La funzione generatrice di Lie al 1° ordine che separa i termini short-period
dalla parte secolare.

**Teoria**:
```
W₁ = -Σ R_{j≠0}(short-period) / (j·n)
```

Dove R_{j≠0} sono i termini della funzione disturbante che dipendono dall'anomalia
media M (i.e. hanno j ≠ 0 nella serie di Fourier). Il divisore j·n è la frequenza
del termine — è qui che i piccoli divisori diventano problematici near resonance.

**Funzioni da creare**:

```python
def _lie_generator_w1(
    elements: OrbitalElements,
    planet_elements: list[PlanetElements],
) -> LieGenerator:
    """Calcola W₁ — la funzione generatrice di Lie al 1° ordine.

    W₁ contiene i termini short-period della funzione disturbante
    divisi per le rispettive frequenze.

    Termini inclusi:
    - Termine diretto: -2μ'α·b_{1/2}^{(1)}·sin(λ-λ') / (n-n')
    - Termine indiretto: +μ'·(a/a')²·sin(λ-λ') / (n-n')
    - Termini in eccentricità: e·sin(2λ-λ'-ϖ), e·sin(λ'-ϖ), etc.

    Nota: per corpi near-resonance, i termini con piccolo divisore
    (|j·n + k·n'| < threshold) vengono esclusi da W₁ e trattati
    come termini a lungo periodo nella Fase 3.
    """
```

```python
def _apply_lie_transform(
    elements: OrbitalElements,
    W1: LieGenerator,
    direction: int,  # +1 = osculating→mean, -1 = mean→osculating
) -> OrbitalElements:
    """Applica la trasformazione di Lie generata da W₁.

    Per ogni elemento orbitale x:
        x_mean = x_osc - {x, W₁}  (direction = +1)
        x_osc = x_mean + {x, W₁}  (direction = -1)

    Dove {x, W₁} è la parentesi di Poisson:
        {x, W₁} = ∂x/∂q · ∂W₁/∂p - ∂x/∂p · ∂W₁/∂q

    Per gli elementi di Delaunay (L, G, H, l, g, h):
        δl = ∂W₁/∂L,  δL = -∂W₁/∂l
        δg = ∂W₁/∂G,  δG = -∂W₁/∂g
        δh = ∂W₁/∂H,  δH = -∂W₁/∂h

    In pratica, per gli elementi Kepleriani standard:
        δa = (2/na) · ∂W₁/∂M
        δe = -(1-e²)^½/(na²e) · ∂W₁/∂M + (1-e²)^½/(na²e) · ∂W₁/∂ω
        δi = -1/(na²(1-e²)^½·sin(i)) · ∂W₁/∂Ω + cos(i)/... · ∂W₁/∂ω
        δM = -2/(na) · ∂W₁/∂a - (1-e²)^½·(...)
        δω = (1-e²)^½/(na²e) · ∂W₁/∂e - cos(i)/... · ∂W₁/∂i
        δΩ = 1/(na²(1-e²)^½·sin(i)) · ∂W₁/∂i
    """
```

**Stima**: ~250 righe (W₁ + trasformazione diretta/inversa)

### Fase 3: Hamiltoniana media e propagazione secolare + lungo periodo

**Cosa**: L'Hamiltoniana mediata H* contiene solo termini secolari e a lungo periodo.
Sostituisce `apply_secular_perturbations()` per i tassi di precessione.

**Teoria**:
```
H* = H₀* + εH₁* + ε²H₂*
H₁* = <R>_M  (media della funzione disturbante sull'anomalia media)
H₂* = ½<{R_sp, W₁}>_M  (termine di 2° ordine dal commutatore)
```

Dove:
- H₁* dà i tassi secolari di 1° ordine (= ciò che abbiamo già in `calc_secular_perturbation_rates()`)
- H₂* aggiunge le correzioni di 2° ordine ai tassi + i **termini a lungo periodo**

**Funzioni da creare/modificare**:

```python
def _calc_mean_hamiltonian_corrections(
    elements: OrbitalElements,
    planet_elements: list[PlanetElements],
    W1: LieGenerator,
) -> SecularCorrections:
    """Calcola le correzioni di 2° ordine dall'Hamiltoniana media H₂*.

    H₂* = ½ · <{R_{short-period}, W₁}>_M

    Questo produce:
    1. Correzioni ai tassi secolari di ω, Ω (migliorano i tassi di 1° ordine)
    2. Termini a LUNGO PERIODO che oscillano con frequenza ~ g - s (decenni)
    3. Correzione secolare al semi-asse maggiore medio ā

    I termini a lungo periodo hanno la forma:
        e_p · sin(ϖ - ϖ_Jupiter) con periodo ~ 1/|g - g_Jupiter|
        sin(Ω - Ω_Jupiter) · sin(i/2) · sin(i_Jupiter/2) con periodo ~ 1/|s - s_Jupiter|

    Questi sono i termini che mancano nella nostra implementazione attuale
    e che causano l'errore residuo di ~15-30" a 1 anno.
    """
```

```python
def _propagate_mean_elements(
    mean_elements: OrbitalElements,
    jd_tt: float,
    secular_corrections: SecularCorrections,
) -> OrbitalElements:
    """Propaga gli elementi medi dal loro epoch al tempo target.

    Combina:
    1. Tassi secolari di 1° ordine (calc_secular_perturbation_rates esistente)
    2. Correzioni di 2° ordine ai tassi (da H₂*)
    3. Oscillazioni a lungo periodo (da H₂*, con sin/cos di frequenze lente)
    4. Evoluzione (h,k)/(p,q) esistente (_calc_forced_elements)

    Restituisce elementi medi al tempo target.
    """
```

**Stima**: ~200 righe

### Fase 4: Integrazione nel pipeline

**Cosa**: Collegare tutto in `calc_minor_body_position()` e `apply_secular_perturbations()`.

**Modifiche**:

1. **`apply_secular_perturbations()`** (L1379): aggiungere parametro `use_hori_deprit=True`
   - Se True: usa il nuovo pipeline (mean→propagate→osculating)
   - Se False: pipeline attuale (backward compatible)

2. **`calc_minor_body_position()`** (L3748): quando Hori-Deprit è attivo, NON chiamare
   `_calc_short_period_correction()` (i termini SP sono già nella trasformazione di Lie)

3. **`_calc_short_period_correction()`** (L3597): mantenerla per backward compatibility
   ma bypassarla nel nuovo pipeline

**Stima**: ~100 righe di modifiche

### Fase 5: Tabelle di coefficienti precalcolati

**Cosa**: Per performance, precalcolare i coefficienti della funzione disturbante
che dipendono solo da α (rapporto semi-assi) e cache-arli.

**Motivazione**: `_calc_laplace_coefficients()` con 200-500 passi è lenta. In Hori-Deprit
servono ~10 coefficienti diversi per ogni coppia asteroide-pianeta, e la funzione
generatrice W₁ va valutata due volte (forward e inverse) per ogni posizione.

**Approccio**:

```python
# Cache per rapporto semi-assi quantizzato (risoluzione 0.001 AU)
_LAPLACE_CACHE: dict[tuple[float, float, int], float] = {}

def _calc_laplace_cached(alpha: float, s: float, j: int) -> float:
    """Versione cached di _calc_laplace_coefficients()."""
    key = (round(alpha, 4), s, j)
    if key not in _LAPLACE_CACHE:
        _LAPLACE_CACHE[key] = _calc_laplace_coefficients(alpha, s, j)
    return _LAPLACE_CACHE[key]
```

**Stima**: ~50 righe

### Fase 6: Validazione e benchmark

**Cosa**: Verificare la precisione rispetto alle posizioni SPK di riferimento.

**Passi**:
1. Eseguire `poe test:keplerian:benchmark` (il benchmark a 37 corpi)
2. Confrontare gli errori prima e dopo Hori-Deprit
3. Target: errore < 5" a 1 anno per i main belt, < 30" per TNO/centauri

**Test specifici**:
- Ceres a ±1 anno: target < 3" (attualmente ~20")
- Vesta a ±1 anno: target < 3" (attualmente ~15")
- Chiron a ±1 anno: target < 10" (attualmente ~30")
- Pallas a ±1 anno: target < 5" (attualmente ~25")
- Juno a ±1 anno: target < 5" (attualmente ~20")

**Stima**: ~100 righe di test aggiuntivi

## Strutture dati

```python
@dataclass
class DisturbingFunctionCoeffs:
    """Coefficienti dell'espansione della funzione disturbante per una coppia
    asteroide-pianeta."""
    alpha: float              # rapporto semi-assi
    # Coefficienti di Laplace necessari
    b_half_0: float           # b_{1/2}^{(0)}(α)
    b_half_1: float           # b_{1/2}^{(1)}(α)
    b_half_2: float           # b_{1/2}^{(2)}(α)
    b_3half_1: float          # b_{3/2}^{(1)}(α) — già usato
    b_3half_2: float          # b_{3/2}^{(2)}(α) — già usato
    # Derivate per la trasformazione inversa
    db_half_1_dalpha: float   # d/dα [b_{1/2}^{(1)}(α)]

@dataclass
class LieGenerator:
    """Funzione generatrice W₁ parametrizzata per valutazione rapida."""
    # Per ogni pianeta perturbatore:
    planet_terms: list[PlanetW1Terms]

@dataclass
class PlanetW1Terms:
    """Termini di W₁ per un singolo pianeta perturbatore."""
    # Amplitudini dei termini di Fourier in W₁
    # Forma: W₁ = Σ A_k · sin(k·M + j·M_planet + l·ω + ...)
    amplitudes: list[float]
    frequencies: list[tuple[int, int, int, int]]  # (k_M, j_M', l_ω, m_Ω)
    planet_n: float   # mean motion del pianeta (rad/day)
    planet_lambda: float  # longitudine media del pianeta al tempo t

@dataclass
class SecularCorrections:
    """Correzioni di 2° ordine ai tassi secolari + termini a lungo periodo."""
    # Correzioni ai tassi di 1° ordine
    delta_d_omega: float  # correzione a dω/dt (rad/day)
    delta_d_Omega: float  # correzione a dΩ/dt (rad/day)
    delta_d_n: float      # correzione a dn (rad/day)
    # Termini a lungo periodo (lista di sinusoidi)
    lp_terms_omega: list[tuple[float, float, float]]  # (ampiezza, frequenza, fase)
    lp_terms_Omega: list[tuple[float, float, float]]
    lp_terms_e: list[tuple[float, float, float]]
    lp_terms_i: list[tuple[float, float, float]]
```

## Termini della funzione disturbante da includere

### 1° ordine in massa (ε = μ_planet)

| Termine | Tipo | Frequenza | Ampiezza tipica (Ceres) |
|---------|------|-----------|-------------------------|
| b_{1/2}^{(1)} · cos(λ-λ_J) | Short-period | n - n_J | ~49" |
| b_{1/2}^{(1)} · cos(λ-λ_S) | Short-period | n - n_S | ~7" |
| e · b_{1/2}^{(2)} · cos(2λ-λ_J-ϖ) | Short-period | 2n-n_J | ~12" |
| e · b_{3/2}^{(2)} · cos(ϖ-ϖ_J) | Secolare | g - g_J | ~5"/yr (lungo periodo) |
| sin(i/2)·sin(i_J/2) · cos(Ω-Ω_J) | Secolare | s - s_J | ~2"/yr (lungo periodo) |

### 2° ordine (ε² = μ²)

| Termine | Tipo | Effetto |
|---------|------|---------|
| {R_sp, W₁}_M averaged | Secolare | Correzione ~10% ai tassi di ω, Ω |
| {R_sp, W₁}_M long-period | Lungo periodo | Oscillazione ~1-3" con P ~ decades |
| e² · b_{1/2}^{(3)} terms | Short-period | Correzione ~2" alla longitudine |

## Stima complessiva

| Fase | Righe di codice | Complessità | Dipendenze |
|------|----------------|-------------|------------|
| 1. Funzione disturbante | ~150 | Media | `_calc_laplace_coefficients()` |
| 2. W₁ + trasformazione Lie | ~250 | Alta | Fase 1 |
| 3. Hamiltoniana media + propagazione | ~200 | Alta | Fase 1, 2 |
| 4. Integrazione pipeline | ~100 | Bassa | Fase 2, 3 |
| 5. Cache coefficienti | ~50 | Bassa | Fase 1 |
| 6. Validazione | ~100 | Media | Tutto |
| **Totale** | **~850** | | |

## Riferimenti

1. **Deprit, A.** (1969). "Canonical transformations depending on a small parameter."
   *Celestial Mechanics*, 1, 12-30. — Paper fondamentale sulla trasformazione di Lie.

2. **Hori, G.** (1966). "Theory of general perturbations with unspecified canonical
   variables." *Publications of the Astronomical Society of Japan*, 18, 287-296.

3. **Murray, C.D. & Dermott, S.F.** (1999). *Solar System Dynamics*. Cambridge
   University Press. Cap. 6 (funzione disturbante), Cap. 7 (teoria secolare).

4. **Brouwer, D.** (1959). "Solution of the problem of artificial satellite theory
   without drag." *Astronomical Journal*, 64, 378-397. — Teoria originale (von Zeipel).

5. **Ferraz-Mello, S.** (2007). *Canonical Perturbation Theories*. Springer.
   Cap. 3-4. — Trattamento moderno e completo di Hori-Deprit.

6. **Ellis, K.M. & Murray, C.D.** (2000). "The disturbing function in solar system
   dynamics." *Icarus*, 147, 129-144. — Espansione sistematica dei coefficienti.

7. **Simon, J.L. et al.** (1994). "Numerical expressions for precession formulae and
   mean elements for the Moon and the planets." *A&A*, 282, 663-683. — Elementi
   planetari usati per le perturbazioni forzate (già implementato come L4).
