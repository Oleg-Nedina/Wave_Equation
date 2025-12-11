### Running the Solver
The executable `wave_solver` accepts command-line arguments to select the dimension and the parameter file.

**Syntax:**
```bash
mpirun -n <N_PROCS> ./wave_solver [DIMENSION] [PARAM_FILE]
```

./wave_solver ../include/parameters.prm


### Examples:
**Run in 2D (Default):**
```bash
mpirun -n 4 ./wave_solver
```

**Run in 3D:**
```bash
mpirun -n 4 ./wave_solver 3
```

**Run in 3D with specific parameters:**
```bash
mpirun -n 8 ./wave_solver 3 parameters.prm
```

---

## 📊 Visualization (ParaView)
1. Open ParaView.  
2. Navigate to the `build/output_results` folder.  
3. Open the group file `solution_..pvtu` (ParaView should group them automatically).  
4. Click **Apply**.

### For 2D visualization:
- Add filter: **Warp By Scalar**  
- Scalars: `displacement`  
- Normal: `0 0 1`  
- Disable **2D mode** to freely rotate the camera  

---

## 🔢 Mathematical Formulation

The full mathematical formulation (strong form, weak form, and the semi-discrete Newmark scheme)  
is documented in detail inside the repository.

➡️ Please refer to the files:

- `./repo` — mathematical formulation and comments  
- `wave_equation.cpp` — implementation details of the scheme  

For readability reasons, the full equations are not rendered in this Markdown file.

---

## 📝 To-Do List (Future Development)

- [ ] **Explicit Time Stepping:** Implement the explicit Newmark scheme (β = 0) with mass lumping to avoid solving linear systems.  
- [ ] **Parameter Parsing:** Connect the `parameters.prm` file to the C++ solver (currently values are hard-coded).  
- [ ] **Time-Dependent Boundary Conditions:** Implement non-homogeneous Dirichlet BCs \( g(t) 
eq 0 \).  
- [ ] **CFL Condition Check:** Automate time-step selection based on mesh size \( h_{\min} \).  

---

## 📁 Repository Structure
```
├── CMakeLists.txt         # Build configuration
├── cmake-common.cmake     # Shared CMake settings
├── main.cpp               # Entry point (MPI init & Arg parsing)
├── wave_equation.hpp      # Class declaration & Documentation
├── wave_equation.cpp      # Implementation (Logic & Math)
└── parameters.prm         # Runtime parameters configuration
```

