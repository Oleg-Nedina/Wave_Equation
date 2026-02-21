import os
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

PARAM_FILE = "../include/parameters.prm"
ENERGY_FILE = "../output_energy/energy.csv"

def set_energy_scenario():
    """Update .prm file to use Scenario 5 and IMPLICIT method"""
    with open(PARAM_FILE, 'r') as f:
        content = f.read()
    
    content = re.sub(r"(set Scenario ID\s*=\s*)(\d+)", r"\g<1>6", content)
    content = re.sub(r"(set Method\s*=\s*)([A-Z]+)", r"\g<1>IMPLICIT", content)
    content = re.sub(r"(set Refinement level\s*=\s*)(\d+)", r"\g<1>6", content)
    
    with open(PARAM_FILE, 'w') as f:
        f.write(content)

print("=== Generation Grafh Energy Conservation ===")

print("\n1/2: Parameter Configuration (Scenario 5, IMPLICIT) and execution...")
set_energy_scenario()
subprocess.run(["mpirun", "-n", "1", "./build/wave_solver", "include/parameters.prm"], cwd="..", check=True)

print("\n2/2: Reading resultsand creating graph...")
if not os.path.exists(ENERGY_FILE):
    print(f"Error: file {ENERGY_FILE} not found. Failed simulation?")
    exit(1)

df = pd.read_csv(ENERGY_FILE)

t = df.iloc[:, 0].to_numpy()
E = df.iloc[:, 1].to_numpy()

# Compute normalized energy with respect to non-zero start value
idx0 = np.argmax(E > 0) 
E0 = E[idx0] if E[idx0] > 0 else 1.0
En = E / E0

# Create Graph
fig, ax1 = plt.subplots(figsize=(10, 6))

# Curve 1: Absolute Energy (Left Axis - Red)
ax1.plot(t, E, "r-", label="E(t) [Absolute Energy]", linewidth=2)
ax1.set_xlabel("Time (t)")
ax1.set_ylabel("Energy E(t)", color="r")
ax1.tick_params(axis="y", labelcolor="r")
ax1.grid(True)

# Curve 2: Normalized Energy (Right Axis - Blue)
ax2 = ax1.twinx()
ax2.plot(t, En, "b--", label="E(t)/E0 [Normalized Energy]", alpha=0.8, linewidth=2)
ax2.set_ylabel("E(t)/E0", color="b")
ax2.tick_params(axis="y", labelcolor="b")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right")

plt.title("Energy Conservation through Time (Implicit Newmark)")
plt.tight_layout()

# Saving
plt.savefig("energy_plot_conservation.png", dpi=300)
print("Done. Graph saved as 'energy_plot_conservation.png'")