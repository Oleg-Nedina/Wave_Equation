import os
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

PARAM_FILE = "../include/parameters.prm"
ENERGY_FILE = "../output_energy/energy.csv"

def set_forcing_scenario():
    """Update .prm file to use Scenario 6 (Pumping Wall)"""
    with open(PARAM_FILE, 'r') as f:
        content = f.read()
    
    content = re.sub(r"(set Scenario ID\s*=\s*)(\d+)", r"\g<1>6", content)
    content = re.sub(r"(set Method\s*=\s*)([A-Z]+)", r"\g<1>IMPLICIT", content)
    content = re.sub(r"(set Refinement level\s*=\s*)(\d+)", r"\g<1>6", content)
    content = re.sub(r"(set Final time\s*=\s*)(\d+\.?\d*)", r"\g<1>1.0", content)
    
    with open(PARAM_FILE, 'w') as f:
        f.write(content)

print("=== Generation Graph: Overlapped Energy and Work (Scenario 6) ===")

print("\n1/2: Parameter configuration and execution...")
set_forcing_scenario()
subprocess.run(["mpirun", "-n", "1", "./build/wave_solver", "include/parameters.prm"], cwd="..", check=True)

print("\n2/2: Reading data and generation graph...")
if not os.path.exists(ENERGY_FILE):
    print(f"Error: file {ENERGY_FILE} not found. Failed simulation?")
    exit(1)

df = pd.read_csv(ENERGY_FILE)
t = df.iloc[:, 0].to_numpy()
E = df.iloc[:, 1].to_numpy()

W_theorical = (t / 2.0) + (np.sin(10.0 * t) / 20.0)

W_theorical = W_theorical * (np.max(E) / np.max(W_theorical))

# Generation Graph
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plotting both curves
ax1.plot(t, E, "r-", label="Total Computed Energy $E(t)$", linewidth=3, alpha=0.8)
ax1.plot(t, W_theorical, "b--", label=r"Theoretical Work $W(t) \propto \int \dot{g}^2 dt$", linewidth=2)

ax1.set_xlabel("Time (t)", fontsize=11)
ax1.set_ylabel("EEnergy / Work", fontsize=11)
ax1.tick_params(axis="both", labelsize=10)
ax1.grid(True)

definition = r"BC (x=0): $g(t,y) = 0.5 \sin(5t) \sin(\pi y)$" 
ax1.legend(loc="upper left", title=definition, fontsize=10, title_fontsize=11)

plt.title("Overlapping of Energy and Work (Scenario 6)", fontsize=12)
plt.tight_layout()

# Salvataggio
plt.savefig("energy_work_overlap.png", dpi=300)
print("Done. Graph saved as 'energy_work_overlap.png'")