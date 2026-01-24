# libephemeris TODO

Questa libreria è un drop-in replacement puro Python per pyswisseph (Swiss Ephemeris). Usa Skyfield e le effemeridi NASA JPL DE440/DE441 invece del codice C di Swiss Ephemeris. L'obiettivo è compatibilità API 1:1 con pyswisseph e precisione sub-arcsecond.

## Missing pyswisseph Functions

- [ ] Implementare la funzione `calc_pctr(jd, planet, center, flags)` che calcola la posizione di un pianeta rispetto a un altro corpo come centro (planet-centric), invece che rispetto alla Terra (geocentric) o al Sole (heliocentric). Questa funzione è usata per calcolare ad esempio la posizione della Luna vista da Marte. pyswisseph la espone come `swe.calc_pctr()`. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare le funzioni `nod_aps(jd, planet, flags, method)` e `nod_aps_ut(jd_ut, planet, flags, method)` che calcolano i nodi (ascending/descending) e le apsidi (perihelion/aphelion) orbitali di qualsiasi pianeta. Il parametro method specifica se usare nodi medi o osculanti. pyswisseph ritorna una tupla con (ascending_node, descending_node, perihelion, aphelion). Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `get_orbital_elements(jd, planet, flags)` che ritorna gli elementi orbitali kepleriani di un corpo: semiasse maggiore (a), eccentricità (e), inclinazione (i), longitudine del nodo ascendente (Ω), argomento del perielio (ω), anomalia media (M), e altri parametri derivati. pyswisseph ritorna questi come una tupla di 17 elementi. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `orbit_max_min_true_distance(jd, planet, flags)` che calcola la distanza massima e minima dalla Terra durante l'orbita del pianeta, utile per determinare quando un pianeta è al perigeo o apogeo. pyswisseph ritorna (min_distance, max_distance). Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `deltat_ex(jd, ephe_flag)` che è una versione estesa di deltat() che permette di specificare esplicitamente quale fonte usare per il calcolo di Delta-T (differenza tra TT e UT). Il parametro ephe_flag può essere SEFLG_SWIEPH, SEFLG_JPLEPH, o SEFLG_MOSEPH. Attualmente deltat() usa sempre un valore di default.

- [ ] Implementare la funzione `date_conversion(year, month, day, hour, calendar)` che converte date tra calendario Giuliano e Gregoriano. Il parametro calendar è 'j' per Giuliano o 'g' per Gregoriano. Ritorna (year, month, day, hour) nel calendario richiesto. La transizione standard è il 15 ottobre 1582.

- [ ] Implementare la funzione `day_of_week(jd)` che ritorna il giorno della settimana (0=Lunedì, 1=Martedì, ..., 6=Domenica) per qualsiasi Julian Day. Questa è una semplice operazione matematica: (floor(jd + 1.5) % 7) ma pyswisseph la espone come funzione per comodità.

- [ ] Implementare la funzione `utc_to_jd(year, month, day, hour, min, sec, calendar)` che converte una data/ora UTC in Julian Day, gestendo correttamente i leap seconds. A differenza di julday() che assume UT1, questa funzione tiene conto della differenza tra UTC e UT1. Ritorna una tupla (jd_et, jd_ut) con entrambe le versioni.

- [ ] Implementare la funzione `jdet_to_utc(jd_et, calendar)` che converte un Julian Day in Ephemeris Time (TT/ET) a data/ora UTC, applicando Delta-T e leap seconds. Ritorna (year, month, day, hour, min, sec). Il parametro calendar specifica se usare Giuliano o Gregoriano.

- [ ] Implementare la funzione `jdut1_to_utc(jd_ut1, calendar)` che converte un Julian Day in UT1 a data/ora UTC. La differenza tra UT1 e UTC è piccola (< 0.9 secondi) ma importante per precisione astronomica. Ritorna (year, month, day, hour, min, sec).

- [ ] Implementare la funzione `utc_time_zone(year, month, day, hour, min, sec, timezone_offset)` che applica un offset di fuso orario a una data/ora UTC. L'offset è in ore (es. +1 per CET, -5 per EST). Ritorna la data/ora locale come tupla.

- [ ] Implementare la funzione `time_equ(jd)` che calcola l'Equazione del Tempo, cioè la differenza tra tempo solare apparente e tempo solare medio. Ritorna il valore in giorni (moltiplicare per 1440 per minuti). L'equazione del tempo varia da circa -14 a +16 minuti durante l'anno.

- [ ] Implementare la funzione `lat_to_lmt(jd_lat, longitude)` che converte Local Apparent Time (tempo solare vero) a Local Mean Time (tempo solare medio) per una data longitudine. La differenza è l'equazione del tempo. Ritorna il Julian Day in LMT.

- [ ] Implementare la funzione `lmt_to_lat(jd_lmt, longitude)` che converte Local Mean Time a Local Apparent Time per una data longitudine. È l'operazione inversa di lat_to_lmt(). Ritorna il Julian Day in LAT.

- [ ] Implementare la funzione `sidtime(jd, longitude, obliquity, nutation)` che calcola il tempo siderale locale per una data posizione geografica. Il tempo siderale è l'angolo orario del punto vernale e serve per calcolare quali stelle sono visibili. Ritorna il valore in ore (0-24).

- [ ] Implementare la funzione `sidtime0(jd, obliquity, nutation)` che calcola il tempo siderale a Greenwich (longitudine 0°) a mezzanotte UT. È la base per calcolare il tempo siderale locale. Ritorna il valore in ore. pyswisseph usa l'algoritmo IAU 2006.

- [ ] Implementare la funzione `houses_ex2(jd, lat, lon, hsys, flags)` che è una versione estesa di houses_ex() che ritorna anche le velocità (derivate) delle cuspidi delle case. pyswisseph ritorna (cusps, ascmc, cusps_speed, ascmc_speed). Questo è utile per calcolare quando una cuspide cambia segno. Vedere CALCS.md per i dettagli sui calcoli delle case.

- [ ] Implementare la funzione `houses_armc(armc, lat, obliquity, hsys)` che calcola le case a partire dall'ARMC (Ascensione Retta del Medium Coeli) invece che dal Julian Day. Questo è utile quando si ha già l'ARMC precalcolato o per sistemi di case che dipendono solo dall'ARMC. Vedere CALCS.md per i dettagli sui calcoli delle case.

- [ ] Implementare la funzione `houses_armc_ex2(armc, lat, obliquity, hsys, flags)` che combina houses_armc() con il calcolo delle velocità come houses_ex2(). Ritorna (cusps, ascmc, cusps_speed, ascmc_speed). Vedere CALCS.md per i dettagli sui calcoli delle case.

- [ ] Implementare la funzione `house_pos(armc, lat, obliquity, hsys, lon, lat_body)` che determina in quale casa si trova un corpo celeste data la sua longitudine e latitudine eclittica. Ritorna un valore decimale dove la parte intera è il numero della casa (1-12) e la parte decimale indica la posizione all'interno della casa (0.0 = inizio cuspide, 0.999 = fine casa). Vedere CALCS.md per i dettagli sui calcoli delle case.

- [ ] Implementare la funzione `gauquelin_sector(jd, planet, lat, lon, altitude, pressure, temperature, flags, method)` che calcola il settore Gauquelin (1-36) in cui si trova un pianeta. I settori Gauquelin sono una divisione dell'eclittica usata nella ricerca statistica astrologica di Michel Gauquelin. Vedere CALCS.md per i dettagli sui calcoli delle case.

- [ ] Implementare la funzione `get_ayanamsa_ex(jd, sid_mode, flags)` che ritorna l'ayanamsa (precessione siderale) con componenti aggiuntive. A differenza di get_ayanamsa() che ritorna solo il valore totale, questa ritorna una tupla con (ayanamsa, eps_true, nut_long) includendo obliquità vera e nutazione in longitudine. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `get_ayanamsa_ex_ut(jd_ut, sid_mode, flags)` che è la versione UT di get_ayanamsa_ex(). Converte internamente da UT a TT prima di calcolare. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `fixstar(star_name, jd, flags)` che calcola la posizione di una stella fissa in Terrestrial Time (TT). Attualmente è implementata solo fixstar_ut() che usa Universal Time. La differenza è l'applicazione di Delta-T. Vedere CALCS.md per i dettagli sui calcoli di posizione delle stelle fisse.

- [ ] Implementare le funzioni `fixstar2(star_name, jd, flags)` e `fixstar2_ut(star_name, jd_ut, flags)` che sono versioni migliorate di fixstar() con lookup più flessibile. Accettano il nome della stella, un numero di catalogo, o una stringa di ricerca parziale. Ritornano anche il nome completo della stella trovata. Vedere CALCS.md per i dettagli sui calcoli di posizione delle stelle fisse.

- [ ] Implementare le funzioni `fixstar_mag(star_name)` e `fixstar2_mag(star_name)` che ritornano solo la magnitudine visuale di una stella fissa senza calcolare la posizione. Utile quando serve solo la luminosità per calcoli di visibilità.

- [ ] Implementare la funzione `sol_eclipse_when_glob(jd_start, flags, eclipse_type)` che cerca la prossima eclissi solare globale dopo una data specificata. Il parametro eclipse_type filtra per tipo (totale, anulare, parziale). Ritorna una tupla con i tempi delle varie fasi dell'eclissi e informazioni sul tipo. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `sol_eclipse_when_loc(jd_start, lat, lon, altitude, flags)` che cerca la prossima eclissi solare visibile da una specifica località. Ritorna i tempi delle fasi dell'eclissi come viste da quella posizione, inclusi primo contatto, massimo, ultimo contatto. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `sol_eclipse_where(jd, flags)` che calcola la posizione geografica (lat, lon) dove un'eclissi solare è visibile come totale/anulare a un dato momento. Ritorna anche la larghezza del percorso di totalità. Usata per tracciare il percorso dell'eclissi sulla Terra. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `sol_eclipse_how(jd, lat, lon, altitude, flags)` che calcola le circostanze di un'eclissi solare per una specifica località a un dato momento. Ritorna magnitudine, oscuramento, angolo del sole, ecc. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `lun_eclipse_when(jd_start, flags, eclipse_type)` che cerca la prossima eclissi lunare dopo una data specificata. Ritorna i tempi delle fasi (penombra, ombra parziale, ombra totale) e il tipo di eclissi. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `lun_eclipse_when_loc(jd_start, lat, lon, altitude, flags)` che cerca la prossima eclissi lunare visibile da una specifica località, tenendo conto che la Luna deve essere sopra l'orizzonte durante l'eclissi. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `lun_eclipse_how(jd, lat, lon, altitude, flags)` che calcola le circostanze di un'eclissi lunare per una specifica località. Ritorna magnitudine, posizione della Luna, altezza sopra l'orizzonte, ecc. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `lun_occult_when_glob(jd_start, planet, star_name, flags)` che cerca la prossima occultazione lunare di un pianeta o stella fissa. Un'occultazione avviene quando la Luna passa davanti a un altro corpo celeste. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `lun_occult_when_loc(jd_start, planet, star_name, lat, lon, altitude, flags)` che cerca la prossima occultazione lunare visibile da una specifica località. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `lun_occult_where(jd, planet, star_name, flags)` che calcola dove sulla Terra un'occultazione lunare è visibile a un dato momento. Ritorna le coordinate geografiche del centro e dei limiti del percorso. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `rise_trans(jd, planet, lat, lon, altitude, pressure, temperature, flags, rsmi)` che calcola i tempi di levata, tramonto e transito (culminazione) di un corpo celeste per una specifica località. Il parametro rsmi specifica quale evento cercare (rise=1, set=2, transit=4, lower transit=8). Gestisce correttamente i corpi circumpolari che non tramontano mai. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `rise_trans_true_hor(jd, planet, lat, lon, altitude, pressure, temperature, horizon_altitude, flags, rsmi)` che è come rise_trans() ma permette di specificare un'altezza dell'orizzonte diversa da 0°, utile per località con montagne o edifici che occludono l'orizzonte reale. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `solcross(target_lon, jd_start, flags)` che è la versione TT di solcross_ut(). Cerca quando il Sole attraversa una specifica longitudine eclittica, ma usa Terrestrial Time invece di Universal Time per l'input e l'output. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `mooncross(target_lon, jd_start, flags)` che è la versione TT di mooncross_ut(). Cerca quando la Luna attraversa una specifica longitudine eclittica usando Terrestrial Time. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare le funzioni `mooncross_node(jd_start, flags)` e `mooncross_node_ut(jd_start, flags)` che cercano quando la Luna attraversa il proprio nodo (ascending o descending). Questo è importante per calcolare le eclissi, che avvengono solo quando Sole e Luna sono vicini ai nodi lunari. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare le funzioni `helio_cross(planet, target_lon, jd_start, flags)` e `helio_cross_ut(planet, target_lon, jd_start, flags)` che cercano quando un pianeta attraversa una specifica longitudine eclittica eliocentrica (vista dal Sole). Utile per calcoli di astrologia eliocentrica. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `cotrans(coord, obliquity)` che trasforma coordinate tra sistema eclittico e equatoriale. Il parametro coord è una tupla (lon, lat, dist) e obliquity è l'obliquità dell'eclittica. Ritorna le coordinate trasformate. La direzione della trasformazione dipende dal segno dell'obliquità. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `cotrans_sp(coord, speed, obliquity)` che è come cotrans() ma trasforma anche le velocità (derivate) delle coordinate. Il parametro speed è una tupla (lon_speed, lat_speed, dist_speed). Ritorna (coord_trasf, speed_trasf). Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `azalt(jd, calc_flag, lat, lon, altitude, pressure, temperature, coord)` che converte coordinate equatoriali (RA, Dec) in coordinate orizzontali (Azimut, Altezza). Tiene conto della rifrazione atmosferica. Ritorna (azimut, altezza_vera, altezza_apparente). Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `azalt_rev(jd, calc_flag, lat, lon, altitude, azimut, altezza)` che è l'inverso di azalt(): converte coordinate orizzontali in equatoriali. Utile per determinare quale stella si trova in una specifica direzione del cielo. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `refrac(altitude, pressure, temperature, calc_flag)` che calcola la rifrazione atmosferica per una data altezza sopra l'orizzonte. La rifrazione fa apparire gli oggetti celesti più alti di quanto siano realmente, specialmente vicino all'orizzonte (circa 34 arcminuti a 0°). Ritorna la correzione in gradi.

- [ ] Implementare la funzione `refrac_extended(altitude, altitude_geo, pressure, temperature, lapse_rate, calc_flag)` che è una versione estesa di refrac() che usa un modello atmosferico più sofisticato con lapse rate (variazione di temperatura con l'altitudine) personalizzabile.

- [ ] Implementare le funzioni `pheno(jd, planet, flags)` e `pheno_ut(jd_ut, planet, flags)` che calcolano i fenomeni planetari: fase (0-1 dove 0=nuova, 0.5=piena), angolo di fase, elongazione dal Sole, diametro apparente, magnitudine visuale. Ritorna una tupla con tutti questi valori. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `heliacal_ut(jd_start, lat, lon, altitude, pressure, temperature, humidity, body, event_type, flags)` che calcola quando un corpo celeste ha la sua levata o tramonto eliacale (prima/ultima visibilità all'alba/tramonto). Questo era fondamentale per i calendari antichi. Ritorna il Julian Day dell'evento. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `heliacal_pheno_ut(jd, lat, lon, altitude, pressure, temperature, humidity, body, flags)` che calcola le circostanze dettagliate di un evento eliacale: altezza del corpo, altezza del Sole, azimut, magnitudine limite, visibilità. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare la funzione `vis_limit_mag(jd, lat, lon, altitude, pressure, temperature, humidity, body, flags)` che calcola la magnitudine limite per la visibilità di un corpo celeste date le condizioni atmosferiche. Utile per determinare se una stella o pianeta debole sarà visibile.

- [ ] Implementare la funzione `degnorm(angle)` che normalizza un angolo in gradi all'intervallo 0-360°. Equivalente a `angle % 360` ma gestisce correttamente i numeri negativi. pyswisseph la espone per comodità e consistenza.

- [ ] Implementare la funzione `radnorm(angle)` che normalizza un angolo in radianti all'intervallo 0-2π. Equivalente a `angle % (2*pi)` ma gestisce correttamente i numeri negativi.

- [ ] Implementare la funzione `deg_midp(a, b)` che calcola il punto medio tra due angoli in gradi, gestendo correttamente il wraparound a 360°. Per esempio, il punto medio tra 350° e 10° è 0° (o 360°), non 180°.

- [ ] Implementare la funzione `rad_midp(a, b)` che calcola il punto medio tra due angoli in radianti, gestendo correttamente il wraparound a 2π.

- [ ] Implementare la funzione `difdegn(a, b)` che calcola la differenza tra due angoli normalizzata all'intervallo 0-360° (sempre positiva). A differenza di difdeg2n() che ritorna -180..+180, questa ritorna sempre un valore positivo.

- [ ] Implementare la funzione `difrad2n(a, b)` che calcola la differenza tra due angoli in radianti normalizzata all'intervallo -π..+π. È l'equivalente in radianti di difdeg2n().

- [ ] Implementare la funzione `difcs2n(a, b)` che calcola la differenza tra due angoli in centesimi di secondo d'arco (centiseconds, 1/100 di arcsec) normalizzata all'intervallo -180°..+180° (in centiseconds).

- [ ] Implementare la funzione `difcsn(a, b)` che calcola la differenza tra due angoli in centiseconds normalizzata all'intervallo 0-360° (in centiseconds, sempre positiva).

- [ ] Implementare la funzione `csnorm(cs)` che normalizza un valore in centiseconds all'intervallo 0-360° (0 - 360*3600*100 centiseconds).

- [ ] Implementare la funzione `csroundsec(cs)` che arrotonda un valore in centiseconds al secondo d'arco più vicino. Divide per 100 e arrotonda.

- [ ] Implementare la funzione `d2l(value)` che converte un double a long integer con arrotondamento. pyswisseph la usa internamente ma la espone anche pubblicamente.

- [ ] Implementare la funzione `split_deg(degrees, roundflag)` che divide un angolo in gradi nelle sue componenti: segno zodiacale (0-11), gradi (0-29), minuti (0-59), secondi (0-59), frazioni di secondo. Il roundflag controlla l'arrotondamento. Ritorna una tupla con tutte le componenti.

- [ ] Implementare la funzione `cs2degstr(cs)` che converte un valore in centiseconds a una stringa formattata in gradi, minuti, secondi (es. "123°45'67.89\"").

- [ ] Implementare la funzione `cs2lonlatstr(cs, plus_char, minus_char)` che converte un valore in centiseconds a una stringa di longitudine o latitudine con carattere direzionale (es. "45°30'00\" N" o "122°15'30\" W").

- [ ] Implementare la funzione `cs2timestr(cs)` che converte un valore in centiseconds a una stringa di tempo formattata come ore, minuti, secondi (es. "12:34:56").

- [ ] Implementare la funzione `set_jpl_file(filename)` che specifica quale file di effemeridi JPL usare (es. "de441.bsp"). libephemeris usa Skyfield che scarica automaticamente i file, ma questa funzione dovrebbe permettere di specificare un file locale.

- [ ] Implementare le funzioni `set_tid_acc(value)` e `get_tid_acc()` per impostare e leggere l'accelerazione mareale usata nei calcoli di Delta-T. Il valore di default è basato su DE431, ma può essere personalizzato per studi storici.

- [ ] Implementare la funzione `set_delta_t_userdef(dt)` che permette di specificare un valore di Delta-T definito dall'utente invece di usare i valori tabulati. Utile per test o per date molto antiche/future dove Delta-T è incerto.

- [ ] Implementare la funzione `set_lapse_rate(lapse_rate)` che imposta il tasso di variazione della temperatura atmosferica con l'altitudine (lapse rate) per i calcoli di rifrazione. Il default è circa 0.0065 K/m.

- [ ] Implementare la funzione `close()` che chiude tutti i file di effemeridi aperti e libera le risorse. In libephemeris con Skyfield questo potrebbe significare rilasciare la cache dei kernel SPK.

- [ ] Implementare la funzione `get_library_path()` che ritorna il percorso dove sono memorizzati i file di effemeridi (per pyswisseph i file .se1, per libephemeris i file .bsp di Skyfield).

- [ ] Implementare la funzione `get_current_file_data(ifno)` che ritorna informazioni sul file di effemeridi attualmente in uso: nome del file, intervallo di date coperto, tipo di effemeridi. Il parametro ifno seleziona quale file (0=pianeti, 1=Luna, 2=asteroidi).

- [ ] Implementare la funzione `get_planet_name(planet_id)` che ritorna il nome di un pianeta dato il suo ID numerico. Per esempio, get_planet_name(0) ritorna "Sun", get_planet_name(1) ritorna "Moon", ecc. Utile per messaggi di errore e debug.

- [ ] Implementare la classe eccezione `Error` che corrisponde a `swe.Error` di pyswisseph. Deve essere sollevata in tutti i casi dove pyswisseph solleva swe.Error (effemeridi non trovate, pianeta non supportato, date fuori range, ecc.) per garantire compatibilità del codice client.

## Precision Improvements

- [ ] Correggere il fallimento del calcolo delle case per latitudini polari (>66.5°) in houses.py linea 678. Attualmente alcune case non vengono calcolate correttamente quando l'eclittica non interseca l'orizzonte in modo normale. Bisogna implementare la gestione speciale per case circumpolari. Vedere CALCS.md per i dettagli sui calcoli delle case.

- [ ] Correggere il secondo caso di fallimento per latitudini polari in houses.py linea 1003. Stesso problema del precedente ma in un diverso sistema di case. Vedere CALCS.md per i dettagli sui calcoli delle case.

- [ ] Correggere l'approssimazione nel sistema di case Alcabitus in houses.py linea 1918 dove viene usata un'approssimazione perché non si ha accesso al vero Ascendente calcolato. Bisogna refactorizzare per passare l'Ascendente già calcolato. Vedere CALCS.md per i dettagli sui calcoli delle case.

- [ ] Correggere l'uso dei baricentri invece dei centri planetari per i giganti gassosi in planets.py linea 61. Attualmente Skyfield ritorna la posizione del baricentro Jupiter-system invece del centro di Giove, causando errori di alcune migliaia di km. Bisogna usare gli ID corretti per i centri planetari (599 per Giove invece di 5). Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Correggere l'approssimazione lineare del moto proprio delle stelle in planets.py linea 757. Attualmente si usa `pos = pos_j2000 + proper_motion * years` ma per precisione sub-arcsecond su lunghi periodi serve un'espansione polinomiale che tenga conto della curvatura del moto. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Correggere la stessa approssimazione lineare del moto proprio in planets.py linea 769. È lo stesso problema del precedente ma in un punto diverso del codice. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Migliorare la tolleranza di convergenza Newton-Raphson in crossing.py linea 12. Attualmente la tolleranza è 1 arcsecondo ma per compatibilità con pyswisseph serve precisione sub-arcsecond (0.1" o meglio). Ridurre la tolleranza e aumentare le iterazioni massime se necessario. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Stringere la tolleranza di convergenza da 1 arcsec in crossing.py linea 57 per il calcolo degli attraversamenti solari (solcross). pyswisseph raggiunge precisione di circa 0.001 arcsec. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Aumentare il numero massimo di iterazioni da 20 in crossing.py linea 97. Per pianeti lenti come Plutone o per attraversamenti vicini a stazioni retrograde, 20 iterazioni potrebbero non bastare per convergere. Usare un valore adattivo basato sulla velocità del pianeta. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Stringere la tolleranza di convergenza in crossing.py linea 113 che attualmente è 1 arcsecondo. Stesso problema dei precedenti. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Correggere la tolleranza di 1 arcsec per gli attraversamenti lunari in crossing.py linea 155. La Luna si muove velocemente (~13°/giorno) quindi la convergenza dovrebbe essere veloce, ma la tolleranza finale deve comunque essere sub-arcsecond. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Rivedere il numero di iterazioni (30) per la Luna in crossing.py linea 191. Potrebbe essere eccessivo dato che la Luna converge rapidamente, o potrebbe servire ancora di più per casi edge vicino ai nodi. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Ottimizzare la gestione delle stazioni retrograde in crossing.py linea 248. Quando un pianeta è vicino a una stazione (velocità ~0), l'algoritmo Newton-Raphson può avere problemi di convergenza. Bisogna detectare questa situazione e usare un metodo alternativo (bisection o Brent). Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Implementare iterazioni adattive per pianeti lenti in crossing.py linea 298. Plutone si muove di soli 0.01°/giorno quindi serve più precisione e più iterazioni rispetto a Mercurio. Il numero di iterazioni dovrebbe dipendere dalla velocità media del pianeta. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Correggere l'approssimazione lineare del moto proprio in fixed_stars.py linea 14. Per le stelle con alto moto proprio (come Barnard's Star) l'errore diventa significativo su scale di secoli. Usare espansione al secondo ordine. Vedere CALCS.md per i dettagli sui calcoli di posizione delle stelle fisse.

- [ ] Correggere la stessa approssimazione in fixed_stars.py linea 99. Vedere CALCS.md per i dettagli sui calcoli di posizione delle stelle fisse.

- [ ] Correggere il moto proprio che ignora la curvatura in fixed_stars.py linea 118. Per precisione massima bisogna considerare che il moto proprio cambia direzione nel tempo a causa della curvatura della sfera celeste. Vedere CALCS.md per i dettagli sui calcoli di posizione delle stelle fisse.

- [ ] Sostituire l'approssimazione a 2 termini della nutazione con il modello IAU completo in fixed_stars.py linea 162. Attualmente si usano solo i termini principali (18.6 anni e 6 mesi) ma il modello IAU 2000A ha 1365 termini per precisione sub-milliarcsecond. Vedere CALCS.md per i dettagli sui calcoli di posizione delle stelle fisse.

- [ ] Rivedere il range di validità del polinomio di Meeus per il nodo lunare medio in lunar.py linea 49. Il polinomio è ottimizzato per date vicine a J2000; per date molto antiche o future l'errore cresce. Documentare i limiti o usare una formula più robusta. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Sostituire l'obliquità fissa J2000 (23.4392911°) con l'obliquità dinamica in lunar.py linea 114. L'obliquità dell'eclittica cambia nel tempo a causa della precessione; usare il valore fisso introduce un errore che cresce con la distanza da J2000. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Aggiungere le perturbazioni planetarie alla formula semplificata in lunar.py linea 157. La posizione della Luna è perturbata principalmente da Giove e dal Sole; ignorare queste perturbazioni può causare errori di diversi arcminuti. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Usare un parametro gravitazionale più preciso in lunar.py linea 234. Il valore di GM (costante gravitazionale × massa) della Terra è noto con precisione molto alta da misure satellitari; usare il valore più recente IAU. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Sostituire l'obliquità fissa J2000 in lunar.py linea 250. Stesso problema della linea 114. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Aggiornare gli elementi orbitali osculanti dalla epoca 2023.0 in minor_bodies.py linea 81. Gli elementi kepleriani degli asteroidi cambiano nel tempo a causa delle perturbazioni; bisogna o usare elementi più recenti o implementare l'integrazione delle perturbazioni. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Gestire le orbite paraboliche e iperboliche (e >= 1) in minor_bodies.py linea 259. L'equazione di Keplero standard funziona solo per orbite ellittiche (e < 1); per comete con orbite paraboliche o iperboliche serve una formula diversa. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

- [ ] Aggiungere le perturbazioni alla propagazione kepleriana 2-body in minor_bodies.py linea 298. La propagazione kepleriana pura ignora l'influenza gravitazionale di Giove e altri pianeti, causando errori crescenti nel tempo. Implementare almeno le perturbazioni principali o usare elementi medi invece di osculanti. Vedere CALCS.md per i dettagli sui calcoli di posizione planetaria.

## Test Coverage

- [ ] Aggiungere test per tutti i 43 modi ayanamsha definiti in pyswisseph, confrontando i valori calcolati da libephemeris con quelli di pyswisseph per diverse date (J2000, 1900, 2100). Usare tolleranze appropriate: 0.001° per ayanamsha basati su formule, 1.0° per quelli basati su stelle fisse che hanno definizioni ambigue.

- [ ] Aggiungere test per tutti i 13+ sistemi di case a latitudini polari (>66.5° N/S) dove alcuni sistemi falliscono o producono risultati degeneri. Verificare che libephemeris gestisca questi casi come pyswisseph (ritornando case "piatte" o sollevando errori appropriati).

- [ ] Aggiungere test completi per la rilevazione dei periodi retrogradi di Mercurio, Venere, Marte, Giove e Saturno. Verificare che le date di inizio/fine retrogradazione calcolate da libephemeris corrispondano a pyswisseph entro 1 minuto.

- [ ] Aggiungere test di confronto di precisione sub-arcsecond contro pyswisseph per Sole, Luna, Mercurio, Venere, Marte, Giove, Saturno. Calcolare posizioni per 1000+ date casuali e verificare che la differenza massima sia < 1 arcsecondo per i pianeti interni e < 5 arcsecondo per i pianeti esterni.

- [ ] Aggiungere test di thread-safety per l'uso concorrente di EphemerisContext. Creare 10+ thread che calcolano simultaneamente posizioni planetarie con diverse configurazioni (diversi ayanamsha, diverse località) e verificare che i risultati siano corretti e che non ci siano race conditions.

- [ ] Aggiungere test per i casi edge delle date: transizione Giuliano/Gregoriano (4-15 ottobre 1582), anno 0 (che non esiste nel calendario storico), anni negativi (BCE), date molto antiche (-3000) e molto future (+3000). Verificare che julday() e revjul() gestiscano tutti questi casi come pyswisseph.

- [ ] Aggiungere test per le funzioni di crossing (solcross_ut, mooncross_ut) verificando che i solstizi e equinozi calcolati corrispondano alle date note (es. equinozio di primavera 2024 = 20 marzo 03:06 UTC). La precisione deve essere < 1 secondo.

- [ ] Aggiungere test di integrazione usando dati di nascita di persone famose per calcolare carte natali complete (posizioni planetarie, case, aspetti) e confrontare con risultati noti da software astrologico di riferimento.

- [ ] Aggiungere test di performance/benchmark che misurano il tempo di calcolo di 10000 posizioni planetarie e lo confrontano con pyswisseph. libephemeris dovrebbe essere al massimo 10x più lento di pyswisseph (che è in C).

## Documentation

- [ ] Scrivere una reference API completa che documenti tutte le funzioni implementate con signature, parametri, valori di ritorno, eccezioni possibili, e esempi di utilizzo. Seguire il formato di docstring Google/NumPy per integrazione con Sphinx.

- [ ] Scrivere una guida di migrazione per utenti che passano da pyswisseph a libephemeris, documentando: differenze API (se presenti), differenze di precisione note, funzioni non ancora implementate, e come usare EphemerisContext per thread-safety.

- [ ] Documentare tutte le limitazioni di precisione note, con numeri specifici: "posizioni planetarie: ±0.5 arcsec vs pyswisseph", "case polari: non supportate per lat > 85°", ecc. Questo aiuta gli utenti a decidere se libephemeris è adatto al loro caso d'uso.

- [ ] Creare un cookbook con esempi pratici: calcolo carta natale completa, ricerca prossimo transito, calcolo sinastria, ricerca eclissi, calcolo effemeridi mensili. Ogni esempio deve essere eseguibile e ben commentato.

- [ ] Creare e mantenere un CHANGELOG.md che documenti tutte le modifiche per ogni versione: nuove funzioni, bug fix, breaking changes, miglioramenti di precisione. Seguire il formato Keep a Changelog.

## Quality and Infrastructure

- [ ] Correggere tutti gli errori LSP/mypy nei file context.py, state.py, houses.py, planets.py. Gli errori principali sono: tipi incompatibili con SpiceKernel, operatori non supportati per tipi Skyfield, variabili potenzialmente non definite. Aggiungere type hints corretti e cast dove necessario.

- [ ] Aggiungere il file py.typed nella root del package per indicare supporto PEP 561 per type hints. Questo permette a mypy e altri type checker di usare i type hints di libephemeris quando è usato come dipendenza.

- [ ] Configurare GitHub Actions per CI/CD automatico: eseguire pytest su ogni push/PR, testare su Python 3.9/3.10/3.11/3.12, eseguire mypy per type checking, eseguire ruff/black per linting/formatting.

- [ ] Aggiungere code coverage con pytest-cov e pubblicare i risultati su Codecov o Coveralls. Obiettivo: >80% coverage per il codice core, >60% per utilities.

- [ ] Creare una suite di benchmark che confronta le performance di libephemeris vs pyswisseph per le operazioni più comuni: calc_ut, houses, get_ayanamsa. Pubblicare i risultati nel README.

- [ ] Configurare pre-commit hooks per eseguire automaticamente prima di ogni commit: ruff (linting), black (formatting), mypy (type checking), pytest (test veloci).

- [ ] Preparare il package per la distribuzione su PyPI: pyproject.toml completo con metadata, README.md come descrizione, LICENSE, classifiers appropriati. Testare con `pip install -e .` e `python -m build`.

- [ ] Profilare l'uso di memoria per calcoli lunghi (es. 10000 posizioni in loop) per verificare che non ci siano memory leak. Skyfield carica file di effemeridi grandi; verificare che la cache sia gestita correttamente.

- [ ] Audit completo della thread-safety: identificare tutte le variabili globali/modulo-level (in state.py, constants.py), verificare che EphemerisContext le isoli correttamente, aggiungere lock dove necessario.

- [ ] Assicurarsi che tutti i messaggi di errore corrispondano al formato di pyswisseph per permettere al codice client esistente di fare pattern matching sugli errori. Es. "illegal planet number" invece di "invalid planet".
