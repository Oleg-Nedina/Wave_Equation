# Simulation Scenarios and Test Cases

The `test_cases.hpp` file serves as the "physics engine" of the solver. It contains a collection of classes that inherit from `dealii::Function<dim>`, allowing the solver to evaluate initial conditions, source terms, and exact solutions at any point in space and time.

## 1. Numerical Validation (Scenario 1)

**Classes:** `ExactSolution_Case1`, `InitialValuesV_Case1`, `RightHandSide_Case1`.

This scenario is designed for **convergence verification**. It uses the *Method of Manufactured Solutions* to define a known analytical result:

* **Utility:** It allows for the calculation of the  error between the numerical approximation and the exact solution.
* **Technical Detail:** The source term  is derived by substituting  into the wave equation, resulting in .

---

## 2. Propagation and Superposition (Scenarios 2 & 3)

**Classes:** `InitialValuesU_Case2` (Gaussian Pulse), `InitialValuesU_Case3` (Double Pulse).

These scenarios test how the solver handles traveling waves and reflecting Dirichlet boundaries.

* **Scenario 2 (Gaussian Pulse):** A single pulse centered at . It is used to observe radial symmetry and wave speed accuracy.
* **Scenario 3 (Double Pulse):** Two pulses starting at opposite corners. It demonstrates the **principle of linear superposition** as the waves pass through each other.

---

## 3. Complexity and Energy Concentration (Scenarios 4 & 5)

**Classes:** `InitialValuesU_Raindrops`, `InitialValuesU_Ring`.

* **Raindrop Symphony:** Simulates multiple random drops with varying amplitudes. This creates complex constructive and destructive interference patterns.
* **Ring Implosion:** A Gaussian ring that collapses toward the center. This is a stress test for the solver’s stability, as energy concentrates at a single point (the center).

---

## 4. Dynamic Boundary Conditions (Scenario 6)

**Classes:** `BoundaryFunction_MovingWall`, `InitialValuesV_MovingWall`.

This is the most advanced scenario, featuring a **Pumping Wall**.

* **Mechanism:** The left wall () moves according to the function: .
* **Initial Consistency:** To prevent numerical shocks, `InitialValuesV_MovingWall` ensures that the initial velocity of the domain matches the wall's velocity at  ().

---

## Scenario Summary Table

| ID | Name | Primary Goal | Source Term () |
| --- | --- | --- | --- |
| **1** | Manufactured | Order of convergence verification |  |
| **2** | Gaussian Pulse | Propagation and reflection | 0 |
| **3** | Collision | Wave superposition | 0 |
| **4** | Raindrops | Complex interference patterns | 0 |
| **5** | Ring | Energy concentration (Implosion) | 0 |
| **6** | Moving Wall | Time-dependent Dirichlet BCs | 0 |

---
