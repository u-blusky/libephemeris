# Test Fixes TODO

Questo file contiene i TODO per correggere i test falliti. libephemeris deve essere un drop-in replacement di pyswisseph, quindi quando c'è discrepanza tra l'implementazione e pyswisseph, VA CORRETTA L'IMPLEMENTAZIONE (non i test, a meno che i test stessi siano sbagliati). pyswisseph è usato solo per test incrociati di precisione, non come dipendenza runtime.

---

## Categoria 1: Valori Letterali Errati nei Test delle Costanti (CORREGGERE I TEST)

- [ ] Nel file tests/test_constants.py alle linee 137-162 ci sono 5 test che confrontano le costanti SEFLG_* di libephemeris con i valori di pyswisseph swe.FLG_*, ma hanno anche un terzo valore letterale hardcodato che è SBAGLIATO. L'implementazione di libephemeris è CORRETTA (i valori matchano pyswisseph), ma i valori letterali nei test sono invertiti. Per esempio test_flg_baryctr alla linea 138 dice "assert constants.SEFLG_BARYCTR == swe.FLG_BARYCTR == 4" ma swe.FLG_BARYCTR vale 16384 (non 4) e anche constants.SEFLG_BARYCTR vale 16384 quindi l'implementazione è corretta. La soluzione è rimuovere il terzo operando (il valore letterale) da tutte queste asserzioni lasciando solo il confronto tra libephemeris e pyswisseph. I test da correggere sono: test_flg_baryctr (linea 138, rimuovere "== 4"), test_flg_j2000 (linea 150, rimuovere "== 16"), test_flg_truepos (linea 154, rimuovere "== 16384"), test_flg_noaberr (linea 158, rimuovere "== 512"), test_flg_nogdefl (linea 162, rimuovere "== 1024").

---

## Categoria 2: Alias Pianeti Senza Prefisso SE_ (CORREGGERE L'IMPLEMENTAZIONE)

- [ ] Nel file tests/test_constants.py alla linea 262 c'è un test chiamato test_planet_names_without_se_prefix che verifica "assert constants.SUN == constants.SE_SUN" ma l'alias SUN non esiste nel modulo libephemeris/constants.py. Se si vuole compatibilità con codice che usa nomi senza prefisso SE_, bisogna aggiungere questi alias nel file libephemeris/constants.py. Aggiungere dopo le definizioni dei pianeti (dopo linea 35 circa): SUN = SE_SUN, MOON = SE_MOON, MERCURY = SE_MERCURY, VENUS = SE_VENUS, MARS = SE_MARS, JUPITER = SE_JUPITER, SATURN = SE_SATURN, URANUS = SE_URANUS, NEPTUNE = SE_NEPTUNE, PLUTO = SE_PLUTO, e tutti gli altri che il test verifica. Controllare il test completo per vedere tutti gli alias richiesti.

---

## Categoria 3: Semantica Return Value delle Funzioni fixstar (CORREGGERE I TEST)

- [ ] Nei file tests/test_fixed_stars/test_fixstar.py, tests/test_fixed_stars/test_fixstar2.py e tests/test_fixed_stars/test_nutation_iau2000a.py ci sono circa 13 test che verificano che il terzo valore di ritorno delle funzioni swe_fixstar e swe_fixstar_ut sia una stringa vuota "" su successo. Ma secondo la documentazione di pyswisseph, swe.fixstar ritorna "(xx, stnam, retflags)" dove stnam è il NOME della stella trovata, NON un messaggio di errore. L'implementazione attuale di libephemeris in fixed_stars.py alle linee 1736 e 1789 ritorna il nome canonico della stella su successo, che è il comportamento CORRETTO. I test sono sbagliati perché aspettano stringa vuota. La soluzione è modificare tutti questi test per verificare che il valore ritornato NON contenga messaggi di errore come "could not find" invece di aspettarsi stringa vuota. Cambiare da 'assert err == "", f"Unexpected error: {err}"' a 'assert "could not find" not in err.lower(), f"Unexpected error: {err}"'. I test da correggere sono in test_fixstar.py: test_fixstar_basic_regulus (linea 21), test_fixstar_basic_spica (linea 30), test_fixstar_case_insensitive (linea 50), test_fixstar_vs_swisseph_regulus (linea 78), test_fixstar_vs_swisseph_spica (linea 97), test_fixstar_multiple_dates (linea 158), test_fixstar_with_comma_in_name (linea 171), test_both_functions_work (linea 214). In test_fixstar2.py: test_fixstar2_vs_fixstar_consistency (linea 166), test_fixstar2_ut_vs_fixstar_ut_consistency (linea 232). In test_nutation_iau2000a.py: tutti i test che verificano err == "".

---

## Categoria 4: _resolve_star_id Ritorna 3 Valori Non 2 (CORREGGERE L'IMPLEMENTAZIONE)

- [ ] La funzione _resolve_star_id definita in libephemeris/fixed_stars.py alla linea 1659 ritorna una tupla di 3 elementi (star_id, error_message, canonical_name) ma nel file libephemeris/eclipse.py ci sono 5 chiamate che aspettano solo 2 valori causando "ValueError: too many values to unpack (expected 2)". Le chiamate errate sono alle linee 3434, 3824, 3951, 4203, 4333 dove il codice dice "star_id, _ = _resolve_star_id(star_name)". La soluzione è correggere tutte queste chiamate in eclipse.py cambiandole in "star_id, err, _ = _resolve_star_id(star_name)" e poi aggiungere la gestione dell'errore se err non è None sollevando l'eccezione appropriata. Questo fix risolverà tutti i test in test_lun_occult.py, test_lun_occult_loc.py e test_lun_occult_where.py che falliscono con "ValueError: could not find star name regulus".

- [ ] Nel file tests/test_error_messages.py alla linea 100 c'è un test che chiama "_resolve_star_id('NonExistentStar')" e aspetta solo 2 valori di ritorno con "star_id, error = _resolve_star_id('NonExistentStar')" ma la funzione ne ritorna 3. Questo è un bug nel TEST non nell'implementazione. La soluzione è cambiare in "star_id, error, _ = _resolve_star_id('NonExistentStar')".

---

## Categoria 5: sol_eclipse_where Deve Ritornare geopos[10] e attr[20] (CORREGGERE L'IMPLEMENTAZIONE)

- [ ] IMPORTANTE: Secondo la documentazione di pyswisseph, la funzione sol_eclipse_where deve ritornare "(retflags, geopos, attr)" dove geopos è una tupla di 10 elementi e attr è una tupla di 20 elementi. L'implementazione attuale in libephemeris/eclipse.py alla linea 1827 ritorna geopos come tupla di soli 3 elementi "(central_lon, central_lat, 0.0)" che è SBAGLIATO. I test in tests/test_eclipse/test_sol_eclipse.py aspettano 2 elementi (anche sbagliato). La soluzione corretta è modificare l'implementazione in eclipse.py per ritornare geopos come tupla di 10 float contenente: [0] longitude central line, [1] latitude central line, [2] lon northern limit umbra, [3] lat northern limit umbra, [4] lon southern limit umbra, [5] lat southern limit umbra, [6] lon northern limit penumbra, [7] lat northern limit penumbra, [8] lon southern limit penumbra, [9] lat southern limit penumbra. Stesso discorso per attr che deve avere 20 elementi. Poi aggiornare i test per aspettarsi le dimensioni corrette.

---

## Categoria 6: sol_eclipse_how Deve Ritornare attr[20] (CORREGGERE L'IMPLEMENTAZIONE)

- [ ] IMPORTANTE: Secondo la documentazione di pyswisseph, la funzione sol_eclipse_how deve ritornare "(retflags, attr)" dove attr è una tupla di 20 elementi contenente: [0] fraction of solar diameter covered, [1] ratio lunar/solar diameter, [2] fraction of disc covered (obscuration), [3] diameter core shadow km, [4] azimuth of sun, [5] true altitude of sun, [6] apparent altitude of sun, [7] elongation of moon, [8] magnitude acc. to NASA, [9] saros series number, [10] saros series member number, e altri campi riservati fino a [19]. L'implementazione attuale in libephemeris/eclipse.py alle linee 2122-2132 ritorna attr con solo 8 elementi che è SBAGLIATO. I test aspettano 11 elementi (anche sbagliato). La soluzione è modificare l'implementazione per ritornare esattamente 20 elementi come pyswisseph, poi aggiornare i test per verificare len(attr) == 20.

---

## Categoria 7: Funzioni swe_* Devono Essere Alias o Avere Stessa Signature (CORREGGERE L'IMPLEMENTAZIONE O I TEST)

- [ ] Nel file tests/test_eclipse/test_sol_eclipse.py ci sono test che verificano l'identità delle funzioni con l'operatore "is" (es. linea 361: "assert swe_sol_eclipse_when_loc is sol_eclipse_when_loc"). In pyswisseph NON esistono funzioni con prefisso "swe_", esistono solo "sol_eclipse_when_loc", "sol_eclipse_where", etc. Il prefisso "swe_" è una convenzione di libephemeris per distinguere le funzioni compatibili pyswisseph. Se le funzioni swe_* hanno signature diversa dalle funzioni base, i test con "is" falliranno sempre perché sono oggetti diversi. La soluzione è: (A) se swe_sol_eclipse_when_loc e sol_eclipse_when_loc devono essere la stessa funzione, definire "swe_sol_eclipse_when_loc = sol_eclipse_when_loc" come alias; (B) se hanno signature diversa intenzionalmente, modificare i test per verificare equivalenza funzionale invece di identità oggetto, o rimuovere questi test.

---

## Categoria 8: deltat Ritorna Giorni - Test Legacy Aspetta Secondi (CORREGGERE I TEST)

- [ ] Nel file tests/test_time_legacy.py alla linea 83 il test test_deltat_exists verifica "assert 60 < dt_py < 70" aspettandosi Delta T in secondi, ma sia pyswisseph che libephemeris ritornano Delta T in GIORNI. Il valore ritornato da swe.deltat(2451545.0) è circa 0.00073876 giorni che moltiplicato per 86400 dà circa 63.8 secondi. L'implementazione è CORRETTA, il test è sbagliato. La soluzione è correggere il test per convertire in secondi prima del confronto: cambiare "dt_py = ephem.swe_deltat(jd)" e "assert 60 < dt_py < 70" in "dt_py = ephem.swe_deltat(jd) * 86400" e "assert 60 < dt_py < 70".

---

## Categoria 9: compare_with_swisseph Fixture Ritorna 3 Valori (CORREGGERE I TEST)

- [ ] Nel file tests/test_planets_legacy.py alle linee 17 e 25 ci sono test che chiamano la fixture compare_with_swisseph e aspettano 2 valori di ritorno con "passed, diffs = compare_with_swisseph(...)" ma la fixture definita in tests/conftest.py alle linee 335-366 ritorna 3 valori "(passed, diffs, raw_results)". Questo è un bug nei TEST, la fixture è corretta. La soluzione è correggere i test in test_planets_legacy.py cambiando "passed, diffs = compare_with_swisseph(...)" in "passed, diffs, _ = compare_with_swisseph(...)". I test da correggere sono test_sun_position_j2000 (linea 17) e test_all_planets_geocentric (linea 25).

---

## Categoria 10: Messaggi di Errore Non Corrispondono alle Regex (DECIDERE E ALLINEARE)

- [ ] Nel file tests/test_rise_trans.py alla linea 389 il test test_invalid_planet_raises_error usa "pytest.raises(ValueError, match='Invalid planet ID')" ma l'implementazione in libephemeris/eclipse.py alla linea 4551 solleva "ValueError(f'illegal planet number {planet}.')". Bisogna decidere quale messaggio è quello canonico. Se si vuole compatibilità con pyswisseph, verificare quale messaggio usa pyswisseph per pianeti invalidi e usare quello. In ogni caso allineare test e implementazione. Stesso problema in tests/test_rise_trans_true_hor.py linea 322 e tests/test_heliacal.py linea 540. La soluzione più semplice è cambiare le regex nei test da "Invalid planet ID" a "illegal planet number" per matchare l'implementazione attuale.

---

## Categoria 11: Validazione Latitudine ±90° (DECIDERE COMPORTAMENTO)

- [ ] Nel file tests/test_edge_cases/test_boundaries.py alle linee 78 e 89 i test test_latitude_over_90 e test_latitude_under_minus_90 passano latitudini 91.0 e -91.0 alla funzione swe_houses. L'implementazione solleva "Error('swe_houses: within polar circle, switched to Porphyry')". Bisogna decidere il comportamento corretto: (A) se latitudini > 90 o < -90 sono invalide, l'implementazione dovrebbe sollevare un errore specifico di validazione prima di arrivare al calcolo; (B) se l'errore attuale è corretto, modificare i test per aspettarsi questa eccezione con "pytest.raises(Error, match='within polar circle')". Verificare cosa fa pyswisseph con latitudini fuori range.

---

## Categoria 12: Tolleranze Numeriche Troppo Strette (CORREGGERE I TEST)

- [ ] Nel file tests/test_time/test_julday.py alla linea 81 il test test_fractional_hours usa una tolleranza di 1e-10 che è troppo stretta: "assert jd - jd_midnight == pytest.approx(expected_fraction, abs=1e-10)". La differenza reale è circa 1.2e-07 che supera la tolleranza ma è comunque precisione eccellente per Julian Day. La soluzione è aumentare la tolleranza a 1e-6 o 1e-7.

- [ ] Nel file tests/test_time_legacy.py alla linea 72 il test test_revjul_roundtrip usa una tolleranza di 1e-10 per confrontare le ore: "assert abs(hour - hour_orig) < 1e-10". La differenza reale è circa 1.6e-09 che supera di poco la tolleranza. La soluzione è aumentare la tolleranza a 1e-8.

---

## Categoria 13: Discrepanze Delta T per Date Future (CORREGGERE I TEST O DOCUMENTARE)

- [ ] Nel file tests/test_time/test_deltat.py alle linee 108 e 123 i test test_deltat_various_years_match_swe[2040] e test_deltat_100_dates_match_swe falliscono perché libephemeris usa Skyfield per calcolare Delta T mentre pyswisseph usa un algoritmo diverso, e per date future (dopo ~2030) le predizioni divergono significativamente. Questa è una limitazione nota: Delta T per il futuro è una PREDIZIONE e diversi algoritmi danno risultati diversi. La soluzione è aumentare la tolleranza per date future (permettere 5+ secondi di differenza per anni > 2030) oppure usare pytest.mark.xfail per questi casi con una nota che spiega la divergenza algoritmica per date future.

---

## Categoria 14: True Node Distance Non Esattamente Zero (CORREGGERE I TEST)

- [ ] Nel file tests/test_lunar/test_dynamic_obliquity.py alla linea 95 il test test_true_node_returns_valid_longitude verifica "assert dist == 0.0" aspettandosi che la distanza del True Node sia esattamente zero. Ma il True Node è un punto matematico sull'orbita lunare e può avere una piccola distanza non-zero nell'implementazione. La soluzione è cambiare l'asserzione per verificare che la distanza sia molto piccola invece di esattamente zero: "assert dist < 0.01" o "assert dist == pytest.approx(0.0, abs=0.01)".

---

## Categoria 15: Velocità Stelle Fisse (CORREGGERE I TEST O DOCUMENTARE)

- [ ] Nel file tests/test_precision/test_precision_docs.py alla linea 495 e nel file tests/test_migration_guide.py alla linea 212 i test verificano che le velocità delle stelle fisse siano esattamente zero con "assert pos[3] == 0.0". Ma le stelle fisse hanno un moto proprio dovuto alla precessione degli equinozi, quindi la velocità NON è zero (è circa 3.8e-05 gradi/giorno). pyswisseph ritorna anche velocità non-zero per le stelle fisse. La soluzione è modificare i test per aspettarsi velocità molto piccole ma non zero: "assert abs(pos[3]) < 0.001" oppure documentare che questa è una limitazione nota.

- [ ] Nel file tests/test_fixed_star_velocity.py alla linea 217 il test test_spica_velocity_vs_pyswisseph fallisce perché la velocità calcolata da libephemeris (0.00003803) differisce significativamente da quella di pyswisseph (0.00013941) con un errore del 72%. Questo indica un problema REALE nell'algoritmo di calcolo della velocità delle stelle fisse in libephemeris. La soluzione richiede investigazione dell'implementazione in fixed_stars.py per capire perché le velocità divergono così tanto da pyswisseph e correggere l'algoritmo.

---

## Categoria 16: lun_occult_when_glob Signature e Return Value (VERIFICARE E CORREGGERE)

- [ ] Secondo la documentazione di pyswisseph, lun_occult_when_glob ha signature "lun_occult_when_glob(tjdut, body, flags, ecltype, backwards)" e ritorna "(retflags, tret)" dove tret è una tupla di 10 elementi. Verificare che l'implementazione in libephemeris/eclipse.py corrisponda a questa signature e a questi return values. I test in test_lun_occult.py falliscono principalmente per il problema di _resolve_star_id (Categoria 4), ma dopo quel fix verificare che la struttura dei dati ritornati sia corretta.

---

## Categoria 17: lun_occult_when_loc Signature e Return Value (VERIFICARE E CORREGGERE)

- [ ] Secondo la documentazione di pyswisseph, lun_occult_when_loc ha signature "lun_occult_when_loc(tjdut, body, geopos, flags, backwards)" e ritorna "(retflags, tret, attr)" dove tret è tupla di 10 elementi e attr è tupla di 20 elementi. Verificare che l'implementazione in libephemeris/eclipse.py corrisponda. I test in test_lun_occult_loc.py aspettano una struttura diversa, potrebbero essere i test sbagliati o l'implementazione.

---

## Categoria 18: lun_occult_where Signature e Return Value (VERIFICARE E CORREGGERE)

- [ ] Secondo la documentazione di pyswisseph, lun_occult_where ha signature "lun_occult_where(tjdut, body, flags)" e ritorna "(retflags, geopos, attr)" dove geopos è tupla di 10 elementi e attr è tupla di 20 elementi. Verificare l'implementazione in libephemeris/eclipse.py. In particolare alla linea 4203 c'è il bug di _resolve_star_id che va corretto prima di poter testare correttamente.
