# NMPDE Project 2: Parallel Wave Equation Solver

**Course:** Numerical Methods for Partial Differential Equations (A.Y. 2025/2026)  
**Topic:** Scalable Finite Element Solver for the 2D Wave Equation with Time-Dependent Boundary Conditions.

## 📌 Project Overview

Questo progetto implementa un solutore agli Elementi Finiti (FEM) per l'equazione delle onde scalare in un dominio bidimensionale. Il codice è scritto in C++ utilizzando la libreria **deal.II** ed è parallelizzato tramite **MPI** per operare su architetture a memoria distribuita.

Il core del progetto è il confronto tra due strategie di integrazione temporale della famiglia **Newmark-$\beta$**:
1.  **Metodo Implicito (A-Stable):** Robusto, incondizionatamente stabile, risolto tramite solver lineari paralleli (CG + AMG).
2.  **Metodo Esplicito (Mass Lumping):** Ottimizzato per HPC, condizionatamente stabile (CFL), risolto tramite inversione diagonale della matrice di massa.

---

## 📐 Mathematical Formulation

Il problema forte è definito come segue:

$$
\begin{cases}
\frac{\partial^2 u}{\partial t^2} - \Delta u = f & \text{in } \Omega \times (0, T] \\
u = g(t) & \text{su } \partial\Omega \\
u(\mathbf{x},0) = u_0, \quad \dot{u}(\mathbf{x},0) = v_0 & \text{in } \Omega
\end{cases}
$$

### Semi-Discrete System & Lifting
Discretizzando nello spazio (Galerkin FEM), otteniamo il sistema matriciale:
$$\mathbf{M}\ddot{\mathbf{U}} + \mathbf{K}\mathbf{U} = \mathbf{F}(t)$$

Poiché le condizioni al bordo $g(t)$ sono non omogenee, applichiamo una tecnica di **Lifting Algebrico**. Il vettore delle incognite viene partizionato in nodi interni $\mathcal{I}$ e nodi di bordo $\Gamma$. I termini noti vengono spostati al membro di destra:
$$\mathbf{M}_{\mathcal{I}\mathcal{I}} \ddot{\mathbf{U}}_{\mathcal{I}} + \mathbf{K}_{\mathcal{I}\mathcal{I}} \mathbf{U}_{\mathcal{I}} = \mathbf{F}_{\mathcal{I}} - (\mathbf{K}_{\mathcal{I}\Gamma} \mathbf{g}(t) + \mathbf{M}_{\mathcal{I}\Gamma} \ddot{\mathbf{g}}(t))$$

In `deal.II`, questo è gestito tramite `AffineConstraints` e `MatrixTools::apply_boundary_values`.

---

## 📂 Repository Structure

```text
├── CMakeLists.txt         # Configurazione di build (linka deal.II, Trilinos, MPI)
├── main.cpp               # Entry point: Inizializza MPI e lancia il solver
├── wave_equation.hpp      # Header: Definizione classe, variabili MPI e metodi
├── wave_equation.cpp      # Source: Implementazione logica (Assemblaggio, Newmark)
├── parameters.prm         # File di input: Parametri modificabili a runtime
└── README.md              # Questo File
```

🚀 Build & Run Instructions
Prerequisiti
deal.II (versione 9.3 o superiore) configurato con MPI e Trilinos.

Compilatore C++17 compatibile.

CMake.

Compilazione
Bash

mkdir build
cd build
cmake ..
make -j4  # Compila usando 4 core
Esecuzione
Il codice deve essere lanciato tramite mpirun. Il file dei parametri è opzionale (default: parameters.prm).

# Esempio su 4 processori
mpirun -n 4 ./wave_solver parameters.prm

🛠 Implementation Guide (Developer Notes)Questa sezione dettaglia le scelte architetturali basate sulla Master Guide del progetto.
1. Parallelismo e Gestione della MemoriaUtilizziamo parallel::distributed::Triangulation. La distinzione critica nella gestione dei vettori è:Locally Owned DoFs: I gradi di libertà di cui il processore è responsabile per il calcolo. Usati per inizializzare le matrici e il vettore RHS.Locally Relevant DoFs (Ghosts): Include i DoF "owned" più quelli dei processori vicini necessari per visualizzazione e calcolo gradienti. Usati per i vettori soluzione ($U, V, A$).
2. Time Stepping: Newmark-betaLa classe WaveEquation gestisce due modalità tramite un parametro method_type:Caso A: Implicito (beta = 0.25, gamma = 0.5)Workflow: Assembla matrice $\mathbf{A}_{sys} = \mathbf{M} + \beta \Delta t^2 \mathbf{K}$.Solver: Usa TrilinosWrappers::SolverCG precondizionato con TrilinosWrappers::PreconditionAMG. Questo è l'unico modo per garantire scalabilità su mesh fini.Lifting: Ad ogni step, constraints viene aggiornato con il tempo corrente $t_{n+1}$ e applicato al sistema lineare.Caso B: Esplicito (beta = 0.0, gamma = 0.5)Mass Lumping: Per evitare di invertire la matrice di massa, usiamo una regola di quadratura QGaussLobatto durante l'assemblaggio di $\mathbf{M}$. Questo rende la matrice diagonale (o quasi, per elementi Q1).Solver: L'inversione è banale: $a_i = RHS_i / M_{ii}$. Nessun solver iterativo è necessario.CFL: Il codice calcola automaticamente $\Delta t_{critical} \le h_{min}/c$ per garantire la stabilità numerica (rispetto dei Coni di Luce).
3. Inizializzazione ConsistenteLo schema di Newmark richiede l'accelerazione iniziale $\ddot{U}_0$, che non è data. Il codice implementa un solutore statico al tempo $t=0$:$$\mathbf{M}\ddot{U}_0 = \mathbf{F}(0) - \mathbf{K}U_0$$


🗺 Roadmap di Sviluppo
Se stai collaborando a questo codice, segui questo ordine di implementazione per minimizzare i bug MPI.

Phase 0: Infrastructure

Setup di make_grid (HyperCube) e setup_system con corretta distribuzione dei DoF.

Verifica output .pvtu vuoto su Paraview.

Phase 1: Static Assembly

Implementare assemble_matrices() per M e K (sono costanti).

Implementare il calcolo di  ddot{U}_0 
 .

Phase 2: Implicit Solver (Reference)

Implementare il loop temporale per il caso Implicito.

Gestire AffineConstraints per g(t).

Verifica: Onda stazionaria.

Phase 3: Explicit Solver (HPC)

Implementare il Mass Lumping (quadratura Lobatto).

Sostituire il solver CG con l'inversione vettoriale.

Benchmark: Confrontare il tempo di esecuzione per step rispetto all'implicito.

📚 References
[1] A. Quarteroni, Numerical Models for Differential Problems, Springer. (Ch. 13: Time Discretization).

[2] S. Salsa, Partial Differential Equations in Action, Springer. (Ch. 5: Waves).

[3] deal.II Step-23 (Wave Equation) & Step-40 (MPI Parallelism).

[4] Project Specification Document (projects.pdf).


