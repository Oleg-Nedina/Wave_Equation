import os
import re
import subprocess
import matplotlib.pyplot as plt
import numpy as np

# Setup
refinements = [5, 6, 7, 8]  # We will test these grid sizes
errors = []
hs = []

prm_file = "../include/parameters.prm"

def update_prm(ref_level):
    """Updates parameters.prm with the new refinement level and auto-dt."""
    with open(prm_file, 'r') as f:
        lines = f.readlines()
    
    with open(prm_file, 'w') as f:
        for line in lines:
            if "set Refinement level" in line:
                f.write(f"    set Refinement level = {ref_level}\n")
            elif "set Scenario ID" in line:
                f.write("    set Scenario ID      = 1\n") # Force Scenario 1 for error check
            elif "set Time step" in line:
                f.write("    set Time step        = 0.0\n") # Force Auto-CFL
            elif "set Final time" in line:
                 f.write("    set Final time       = 0.1\n") # Short run is enough
            else:
                f.write(line)

print("Starting Convergence Test...")

for ref in refinements:
    print(f"--- Running Refinement {ref} ---")
    
    # 1. Update Parameter File
    update_prm(ref)
    
    # 2. Run Solver (using 5 processors)
    # Note: Capture stdout to find the error
    result = subprocess.run(["mpirun", "-n", "5", "../build/wave_solver", prm_file], capture_output=True, text=True)
    
    # 3. Parse L2 Error from output
    # Looking for: "L2 Error: 0.00123..."
    output = result.stdout
    match = re.search(r"L2 Error:\s+([0-9\.eE+-]+)", output)
    
    if match:
        err = float(match.group(1))
        errors.append(err)
        # h = 1 / 2^ref
        h = 1.0 / (2**ref)
        hs.append(h)
        print(f"  -> h={h:.4f}, Error={err:.4e}")
    else:
        print("  -> ERROR: Could not find L2 Error in output!")
        print(output) # Print output for debugging
        break

# 4. Plotting
if len(errors) > 1:
    hs = np.array(hs)
    errors = np.array(errors)

    # Calculate Slope (Order of Convergence)
    slope, intercept = np.polyfit(np.log(hs), np.log(errors), 1)
    
    plt.figure(figsize=(8,6))
    plt.loglog(hs, errors, 'o-', label=f'Simulation (Slope={slope:.2f})', linewidth=2)
    
    # Reference Line (Order 2)
    plt.loglog(hs, errors[0] * (hs/hs[0])**2, 'k--', label='Theory Order 2')
    
    plt.xlabel('Grid Size (h)')
    plt.ylabel('L2 Error')
    plt.title(f'Convergence Analysis (Slope = {slope:.2f})')
    plt.grid(True, which="both", ls="-")
    plt.legend()
    plt.savefig("convergence_plot.png")
    print(f"Done! Plot saved to convergence_plot.png. Calculated Order: {slope:.2f}")
else:
    print("Not enough data to plot.")
