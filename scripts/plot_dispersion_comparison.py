import os
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

PARAM_FILE = "../include/parameters.prm"
DISPERSION_FILE = "../output_results/dispersion.csv"

def set_method_and_scenario(method):
    """Update file .prm to use Scenario 1, selected method and sparse grid"""
    with open(PARAM_FILE, 'r') as f:
        content = f.read()
    
    content = re.sub(r"(set Scenario ID\s*=\s*)(\d+)", r"\g<1>1", content)
    content = re.sub(r"(set Method\s*=\s*)([A-Z]+)", f"\\g<1>{method}", content)
    content = re.sub(r"(set Refinement level\s*=\s*)(\d+)", r"\g<1>3", content)
    
    with open(PARAM_FILE, 'w') as f:
        f.write(content)

print("=== Generating Numerical Dispersion Graph ===")

# Execution EXPLICIT
print("\n1/2: Execution simulation EXPLICIT (FE rank 1) over sparse grid...")
set_method_and_scenario("EXPLICIT")
subprocess.run(["mpirun", "-n", "1", "build/wave_solver", "include/parameters.prm"], cwd="..", check=True)

if os.path.exists(DISPERSION_FILE):
    df_fe1 = pd.read_csv(DISPERSION_FILE)
else:
    print(f"Error: file {DISPERSION_FILE} not found!")
    exit(1)

# Execution IMPLICIT
print("\n2/2: Execution simulation IMPLICIT (FE rank 2) over sparse grid...")
set_method_and_scenario("IMPLICIT")
subprocess.run(["mpirun", "-n", "1", "build/wave_solver", "include/parameters.prm"], cwd="..", check=True)

if os.path.exists(DISPERSION_FILE):
    df_fe2 = pd.read_csv(DISPERSION_FILE)
else:
    print(f"Error: file {DISPERSION_FILE} not found!")
    exit(1)

# 3. Graph creation
print("\nCreating graph...")
plt.figure(figsize=(10, 6))

# Plot Exact Solution
plt.plot(df_fe1["time"], df_fe1["u_exact"], 'k--', label="Exact Solution", linewidth=2)

# Plot Solution fe(1)
plt.plot(df_fe1["time"], df_fe1["u_numerical"], 'r-', label="Numeric FE(1) (Explicit)", alpha=0.8)

# Plot Solution fe(2)
plt.plot(df_fe2["time"], df_fe2["u_numerical"], 'b-', label="Numeric FE(2) (Implicit)", alpha=0.8)

plt.xlabel("Time (t)")
plt.ylabel("u(0.5, 0.5)")
plt.title("Comparison Numerical Dispersion over Sparse Grid (Refinement 3)")
plt.legend(loc="upper right")
plt.grid(True)

plt.savefig("dispersion_comparison_sparse.png", dpi=300)
print("Done. Graph saved as 'dispersion_comparison_sparse.png'")