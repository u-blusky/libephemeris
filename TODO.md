# LEB Precision Improvement — TODO

> Piano di riferimento: `docs/leb/precision-improvement-plan.md`
> Branch: `leb-precision-improvement`

## Contesto

Il lavoro di precision improvement del formato LEB (LibEphemeris Binary) è suddiviso in 5 fasi (0-4), eseguite nell'ordine 0 → 1 → 2 → 4 → 3.

Il **base tier** (de440s, 1850-2150) è completo: tutte le 5 fasi sono state eseguite, il file LEB rigenerato (94.3 MB), e i test passano 404/404 con tolleranze strette.

### Cosa è stato fatto (base tier)

- **Phase 0**: Rimossi xfail dai test asteroidi, abilitato SPK auto-download nel CompareHelper
- **Phase 1**: Parametri Chebyshev aggressivi (Moon/Earth 4d, Venus/Mars 16d, Jupiter/Saturn 32d, Asteroidi 8d, Ipotetici 32d, Nutation 16d)
- **Phase 2**: Velocità analitica via derivata Chebyshev trasformata attraverso le matrici di rotazione (sostituisce central difference)
- **Phase 4**: File LEB base rigenerato con tutti i 30 corpi verificati
- **Phase 3**: Tolleranze strette al minimo, documentazione aggiornata (guide.md + design.md)
- **Bug fix**: Tau bug in verify_segment del generatore, arcsec reporting con distanza geocentrica

### Errori misurati (base tier, worst case)

- Posizione pianeti: da 0.0002" (Sole) a 4.85" (Saturno lat) — sub-arcsecond per tutti tranne Saturno lat
- Velocità pianeti: < 0.014 deg/day (Saturno lon è il peggiore)
- Asteroidi posizione: ~0.44"
- Asteroidi velocità lat: 0.19-0.71 deg/day (limite architetturale della pipeline ICRS→eclittica)
- Corpi eclittici: < 0.028"
- Ipotetici: ~0 (1e-14)

### Limiti architetturali noti (non risolvibili senza cambiare formato)

1. **Saturno latitudine (~5")** — La pipeline ICRS→eclittica amplifica gli errori per 1/distanza_geocentrica
2. **Velocità latitudine asteroidi (0.2-0.7 deg/day)** — Stesso meccanismo, peggiore per corpi vicini
3. **Plutone velocità distanza (~9e-5 AU/day)** — Combinazione di distanza estrema e orbita eccentrica

---

## TODO

### Medium Tier (de440, 1550-2650)

- [ ] Rigenerare il file LEB medium: `poe leb:generate:medium:groups`
- [ ] Copiare il file generato in `/Volumes/Data/libephemeris/leb/ephemeris_medium.leb`
- [ ] Runnare i test medium: `pytest tests/test_leb/compare/ -m leb_compare -v`
- [ ] Stringere le tolleranze `TIER_DEFAULTS["medium"]` in `tests/test_leb/compare/conftest.py` al minimo basandosi sugli errori osservati
- [ ] Gestire eventuali fallimenti asteroidi (copertura SPK più stretta per 1550-2650)
- [ ] Verificare che tutti i test passino con le tolleranze strette

### Extended Tier (de441, -5000 a +5000) — Opzionale

- [ ] Rigenerare il file LEB extended: `poe leb:generate:extended:groups`
- [ ] Runnare i test extended e stringere le tolleranze
- [ ] Gestire range SPK asteroidi (coprono solo ~1600-2500)

### Release

- [ ] Upload dei file LEB aggiornati su GitHub Releases: `poe release:leb <version>`
- [ ] Aggiornare gli hash in `libephemeris/download.py`
- [ ] Commit degli hash aggiornati

### Cleanup

- [ ] Verificare che `poe typecheck` passi (mypy)
- [ ] Merge del branch `leb-precision-improvement` in main
