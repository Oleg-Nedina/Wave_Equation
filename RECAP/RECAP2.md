

### 1. Preparazione e Compilazione


* **Comando:**
```bash
mkdir -p build && cd build
cmake ..
make
cd ..

```


* **Cosa fa:** Traduce il tuo codice C++ in un eseguibile parallelo (`wave_solver`).
* **Aspettativa:** Il processo deve arrivare al `100%` senza errori. Se modifichi il file `.cpp`, devi rilanciare `make`.

---

### 2. Validazione Matematica (Convergenza)

Verifica che il metodo FEM stia convergendo con l'ordine corretto ().

* **Script:** `convergence.py`
* **Comando:**

```bash
cd scripts
python3 convergence.py

```


* **Cosa fa:**
1. Modifica automaticamente `parameters.prm` per settare lo **Scenario 1** (Soluzione Esatta).
2. Esegue il solver con raffinamenti crescenti (es. 5, 6, 7, 8).
3. Calcola l'errore  e la pendenza della retta di errore.


* **Aspettativa:**
* Output terminale: `CALCULATED ORDER OF CONVERGENCE: ~2.00`
* Grafico generato: `convergence_plot.png` (una retta che scende).



---

### 3. Analisi Fisica e Scenari (Dispersione e Probe)

Qui verifichi il comportamento dell'onda nel tempo. Questo è il punto dove applichi il metodo "Chiedi a Gemini" per la verifica.

#### Fase A: Esecuzione

1. Apri `include/parameters.prm` e imposta lo `Scenario ID` desiderato (es. `2` per Gaussiana, `6` per Muro Mobile, `1` per Esatta).
2. Imposta un `Final time` adeguato (es. `3.0` per vedere i rimbalzi).
3. Lancia il solver:
```bash
mpirun -n 4 ./build/wave_solver include/parameters.prm

```


*(Questo genera i file `output_results/dispersion.csv` e `output_energy/energy.csv`)*.

#### Fase B: Generazione Grafico

* **Script:** `plot_dispersion.py`
* **Comando:**
```bash
cd scripts
python3 plot_dispersion.py

```


* **Cosa fa:** Legge i dati raccolti al centro del dominio  e crea `dispersion_plot.png`.

#### Fase C: Verifica "Metodo Gemini"

Prendi il grafico generato e confrontalo con la teoria usando questo prompt (o simile) con me o nel tuo report:

> **Prompt Strategico:** "Ho simulato l'equazione delle onde in un dominio quadrato 2D. Lo scenario è [DESCRIZIONE SCENARIO: es. un impulso Gaussiano al centro / un muro che oscilla a sinistra]. Il grafico mostra il segnale registrato al centro del dominio nel tempo. Cosa dovrei aspettarmi teoricamente in termini di ritardi, riflessioni o forma dell'onda? Il mio grafico è coerente?"

* **Aspettativa Scenario 1 (Esatta):** Due linee (Rossa e Nera) sovrapposte.
* **Aspettativa Scenario 2 (Gaussiana):** Un picco a , poi zero, poi "rumore" caotico (riflessioni) dopo .
* **Aspettativa Scenario 6 (Muro):** Piatto (zero) fino a  (Causalità), poi oscillazioni che diventano complesse quando tornano le eco.

---

### 4. Analisi Energetica (Stabilità)

Verifica se il sistema è stabile (l'energia non esplode).

* **Script:** `plot_energy.py`
* **Comando:** (Dopo aver eseguito la simulazione come al punto 3)

```bash
cd scripts 
python3 plot_energy.py

```


* **Cosa fa:** Legge l'energia totale (Cinetica + Potenziale) calcolata ad ogni step e crea `energy_plot.png`.
* **Aspettativa:**
* **Metodo Implicito / Scenari chiusi (1, 2, 3):** Linea piatta o leggermente decrescente (conservazione).
* **Scenario Forzante (6 - Muro):** L'energia cresce (il sistema sta "pompando" energia dentro).
* **Metodo Esplicito (se  sbagliato):** L'energia schizza a infinito (Esplosione numerica).



---

### 5. Performance HPC (Scaling)

Dimostra che il codice sfrutta il calcolo parallelo.

* **Script:** `strong_scaling.py` e `weak_scaling.py`
* **Comandi:**

```bash
cd scripts
python3 strong_scaling.py
python3 weak_scaling.py

```


* **Cosa fa:**
* **Strong:** Lancia il solver con 1, 2, 3, 4 processori sullo stesso problema.
* **Weak:** Lancia il solver aumentando la griglia proporzionalmente ai processori.


* **Aspettativa:**
* `strong_scaling.png`: Una curva che scende (tempo minore con più core) e mostra uno Speedup vicino a quello ideale.
* `weak_scaling.png`: Una linea quasi orizzontale (l'efficienza rimane alta, idealmente attorno a 1.0).



---

### 6. Pulizia Finale

Prima di consegnare o archiviare:

```bash
rm -rf build
rm *.png *.csv
rm include/parameters_temp.prm

```

*(Mantieni solo i file sorgente `.cpp`, `.hpp`, `.py`, `.prm` e `CMakeLists.txt`)*.
