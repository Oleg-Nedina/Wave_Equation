import os
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

PARAM_FILE = "../include/parameters.prm"
ENERGY_FILE = "../output_energy/energy.csv"

def set_dissipative_scenario():
    """Set Gamma > 0.5 to obtain artificial numerical dissipation"""
    with open(PARAM_FILE, 'r') as f:
        content = f.read()
    
    content = re.sub(r"(set Scenario ID\s*=\s*)(\d+)", r"\g<1>5", content)
    content = re.sub(r"(set Method\s*=\s*)([A-Z]+)", r"\g<1>IMPLICIT", content)
    content = re.sub(r"(set Refinement level\s*=\s*)(\d+)", r"\g<1>6", content)
    content = re.sub(r"(set Final time\s*=\s*)(\d+\.?\d*)", r"\g<1>3.0", content)
    
    content = re.sub(r"(set Newmark Gamma\s*=\s*)(\d+\.?\d*)", r"\g<1>0.6", content)
    content = re.sub(r"(set Newmark Beta\s*=\s*)(\d+\.?\d*)", r"\g<1>0.3025", content)
    
    with open(PARAM_FILE, 'w') as f:
        f.write(content)

def restore_conservative_scenario():
    """Reset default parameters"""
    with open(PARAM_FILE, 'r') as f:
        content = f.read()
    
    content = re.sub(r"(set Newmark Gamma\s*=\s*)(\d+\.?\d*)", r"\g<1>0.5", content)
    content = re.sub(r"(set Newmark Beta\s*=\s*)(\d+\.?\d*)", r"\g<1>0.25", content)
    
    with open(PARAM_FILE, 'w') as f:
        f.write(content)

print("=== Graph Generation: Numerical Dissipation (Newmark Gamma > 0.5) ===")

print("\n1/2: Parameter Configuration (Gamma=0.6) and execution...")
set_dissipative_scenario()
subprocess.run(["mpirun", "-n", "1", "./build/wave_solver", "include/parameters.prm"], cwd="..", check=True)

print("\n2/2: Reading data and generating graph...")
if not os.path.exists(ENERGY_FILE):
    print(f"Error: file {ENERGY_FILE} not found. Failed execution?")
    exit(1)

df = pd.read_csv(ENERGY_FILE)
t = df.iloc[:, 0].to_numpy()
E = df.iloc[:, 1].to_numpy()

idx0 = np.argmax(E > 0) 
E0 = E[idx0] if E[idx0] > 0 else 1.0
En = E / E0

# Graph Creation
fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.plot(t, E, "r-", label="Absolute Energy $E(t)$", linewidth=2)
ax1.set_xlabel("Time (t)")
ax1.set_ylabel("Energy $E(t)$", color="r")
ax1.tick_params(axis="y", labelcolor="r")
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(t, En, "b--", label="Normalized Energy $E(t)/E_0$", alpha=0.8, linewidth=2)
ax2.set_ylabel("Normalized Energy", color="b")
ax2.tick_params(axis="y", labelcolor="b")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

definition = r"Numerical Dissipation. ($\gamma=0.6, \beta=0.3$)" 
ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", title=definition)

plt.title("Algorithmic Dissipation of Energy in Time")
plt.tight_layout()

# Salvataggio e ripristino
plt.savefig("energy_dissipation_plot.png", dpi=300)
restore_conservative_scenario()
print("Done. Graph saved as 'energy_dissipation_plot.png'")
print("Default parameters Gamma=0.5 and Beta=0.25 has been restored.")