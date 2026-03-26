# LEB Optimization: Findings & Known Issues

## Risultato ottimizzazione parametri

### Modifiche applicate (base tier)
| Corpo | Vecchio (int/deg) | Nuovo (int/deg) | Risparmio |
|-------|-------------------|-----------------|-----------|
| Pluto | 32/13 | 64/11 | 0.63 MB |
| 9 Uraniani (40-48) | 32/13 | 256/7 | 9.5 MB |
| **Totale** | | | **~10.1 MB** |

File base: 107 MB -> 97 MB (-9.3%)

### Parametri che NON funzionano (testati e falliti)
| Corpo | Proposto | Fitting Error | Motivo |
|-------|----------|---------------|--------|
| Mercury | 32/13 | 1.4" | Orbita eccentrica (e=0.206), alta frequenza |
| Venus | 32/11 | 15.9" | Necessita intervalli corti per perturbazioni |
| Mars | 32/11 | 2.8" | Stesse ragioni di Venus |
| Jupiter | 64/11 | 0.26" | Perturbazioni da Saturno significative |
| Saturn | 64/11 | 0.044" | Perturbazioni reciproche con Giove |
| Uranus | 128/9 | 0.026" | Grado 9 insufficiente per perturbazioni |
| Neptune | 128/9 | 0.007" | Borderline, ma supera il worst-case target |
| Asteroidi (16d) | 16/13 | 15-321" | Orbite irregolari, perturbazioni forti |
| Mean Node | 128/7 | 0.82" | Precessione aggiunge complessita |
| True Node | 32/11 | 24.6" | Perturbazioni lunari ad alta frequenza |

### Insight chiave
1. I margini di precisione riportati dall'accuracy sweep (es. "3800x margine per Jupiter") sono fuorvianti perche lo sweep non testa ai momenti di minima distanza geocentrica (worst case).
2. I fitting errors dei parametri ORIGINALI gia superano 0.001" per alcuni corpi (Jupiter 0.0015", Uranus 0.0025") — il target documentato non e raggiunto al worst case nemmeno ora.
3. Le traiettorie ICRS baricentriche NON sono cosi lisce come ipotizzato — le perturbazioni planetarie richiedono intervalli corti e gradi alti.
4. I corpi eclittici (nodi lunari, apogei) hanno oscillazioni complesse da perturbazioni lunari che richiedono parametri conservativi.

## Bug preesistenti scoperti

### 1. Uraniani geocentrici non supportati (bodies 40-47)
**File**: `libephemeris/planets.py:1854-1879`
**Problema**: `swe_calc(jd, 40, SEFLG_SPEED)` senza `SEFLG_HELCTR` cade nel fallthrough `raise UnknownBodyError`. Il blocco `if SE_CUPIDO <= ipl <= SE_POSEIDON:` ha solo il path `if is_helio:` e non gestisce il caso geocentrico (a differenza di Transpluto/body 48 che ha la conversione geocentrica).
**Impatto**: 16 test failures nella test suite LEB + tutti gli Uraniani SKIP nell'accuracy sweep.
**Severity**: Medium — gli Uraniani sono usati quasi esclusivamente in modo eliocentrico in astrologia, ma la mancanza di supporto geocentrico e un'incompatibilita con pyswisseph.

### 2. Sun eliocentrico (SEFLG_HELCTR) restituisce dati errati
**File**: `tests/test_leb/compare/base/test_base_flags.py:85`
**Problema**: `swe_calc(jd, SE_SUN, SEFLG_SPEED | SEFLG_HELCTR)` produce errore di 646800" (~180 gradi).
**Impatto**: 1 test failure nella suite flags.

### 3. TrueNode distance failure
**File**: `tests/test_leb/compare/base/test_base_lunar.py`
**Problema**: La distanza del True Node supera la tolleranza DISTANCE_AU.
**Impatto**: 1 test failure.
