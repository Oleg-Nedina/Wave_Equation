import subprocess
import re
import os
import json

try:
    import matplotlib.pyplot as plt
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

EXECUTABLE = "build/wave_solver"
PARAM_FILE = "include/parameters.prm"
TEMP_PARAM_FILE = "include/parameters_temp.prm"
DATA_FILE = "strong_data.json"
PROCS_LIST = [1, 2, 3, 4] 
FIXED_REFINEMENT = 8  

def update_parameter_file(ref):
    with open(PARAM_FILE, 'r') as f:
        content = f.read()
    content = re.sub(r"(set Refinement level\s*=\s*)(\d+)", f"\\g<1>{ref}", content)
    content = re.sub(r"(set Final time\s*=\s*)(\d+\.?\d*)", "set Final time = 0.5", content)
    content = re.sub(r"(set Output frequency\s*=\s*)(\d+)", "set Output frequency = 1000000", content)
    with open(TEMP_PARAM_FILE, 'w') as f:
        f.write(content)

def parse_time(output):
    match = re.search(r"\| Total wallclock time elapsed since start\s*\|\s*([\d\.]+)s", output)
    return float(match.group(1)) if match else None

def run_tests():
    update_parameter_file(FIXED_REFINEMENT)
    results = {}
    for p in PROCS_LIST:
        print(f"Running {p} procs...", end=" ", flush=True)
        cmd = ["mpirun", "-n", str(p), EXECUTABLE, TEMP_PARAM_FILE]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            elapsed = parse_time(result.stdout)
            print(f"Done: {elapsed}s")
            results[p] = elapsed
        else:
            print("FAILED")
    with open(DATA_FILE, 'w') as f:
        json.dump(results, f)

def plot_results():
    with open(DATA_FILE, 'r') as f:
        results = json.load(f)
    procs = sorted([int(p) for p in results.keys()])
    times = [results[str(p)] for p in procs]
    speedups = [times[0] / t for t in times]
    
    plt.figure(figsize=(8, 6))
    plt.plot(procs, speedups, 'o-', color='blue', label='Measured Speedup')
    plt.plot(procs, procs, 'k--', label='Ideal Speedup')
    plt.xlabel('Processors')
    plt.ylabel('Speedup')
    plt.title(f'Strong Scaling (Ref {FIXED_REFINEMENT})')
    plt.grid(True)
    plt.legend()
    plt.xticks(PROCS_LIST)
    plt.savefig("strong_scaling.png")
    print("Grafico salvato: strong_scaling.png")

if __name__ == "__main__":
    if not HAS_PLOT:
        run_tests()
    else:
        plot_results()
