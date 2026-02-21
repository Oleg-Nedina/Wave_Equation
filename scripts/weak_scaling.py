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
DATA_FILE = "weak_data.json"

SCALING_CASES = [(1, 7), (4, 8), (16, 9), (64, 10)]

def update_parameter_file(ref):
    with open(PARAM_FILE, 'r') as f:
        content = f.read()
    
    times = {7: 0.5, 8: 0.25, 9: 0.125, 10: 0.0625}
    final_time = times.get(ref, 0.5)

    content = re.sub(r"(set Refinement level\s*=\s*)(\d+)", f"\\g<1>{ref}", content)
    content = re.sub(r"(set Final time\s*=\s*)(\d+\.?\d*)", f"set Final time = {final_time}", content)
    content = re.sub(r"(set Output frequency\s*=\s*)(\d+)", "set Output frequency = 100000", content)
    
    with open(TEMP_PARAM_FILE, 'w') as f:
        f.write(content)

def parse_time(output):
    match = re.search(r"\| Total wallclock time elapsed since start\s*\|\s*([\d\.]+)s", output)
    return float(match.group(1)) if match else None

def run_tests():
    results = {}
    for p, ref in SCALING_CASES:
        print(f"Running {p} procs (Ref {ref})...", end=" ", flush=True)
        update_parameter_file(ref)
        cmd = ["mpirun", "-n", str(p), "./" + EXECUTABLE, TEMP_PARAM_FILE]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            elapsed = parse_time(result.stdout)
            if elapsed:
                print(f"Done: {elapsed}s")
                results[p] = elapsed
            else:
                print("FAILED (Could not parse time)")
        else:
            print(f"FAILED (Code: {result.returncode})")
            print(f"--- STDERR ---\n{result.stderr.strip()}")
            
    with open(DATA_FILE, 'w') as f:
        json.dump(results, f)

def plot_results():
    if not os.path.exists(DATA_FILE):
        print(f"Error: {DATA_FILE} not found.")
        return

    with open(DATA_FILE, 'r') as f:
        results = json.load(f)
        
    procs = [1, 4, 16, 64]
    valid_procs = [p for p in procs if str(p) in results]
    
    if not valid_procs:
        print("No valid data to plot.")
        return

    times = [results[str(p)] for p in valid_procs]
    efficiencies = [times[0] / t for t in times]
    
    plt.figure(figsize=(8, 6))
    plt.plot(valid_procs, efficiencies, "o-", color="green", label="Parallel Efficiency")
    plt.axhline(y=1.0, color="k", linestyle="--", label="Ideal Efficiency")
    plt.ylim(0, 1.2)
    plt.xlabel("Processors")
    plt.ylabel("Efficiency")
    plt.title("Weak Scaling")
    plt.grid(True)
    plt.legend()
    plt.xticks(valid_procs)
    plt.savefig("weak_scaling.png")
    print("Grafico salvato: weak_scaling.png")

if __name__ == "__main__":
    if not HAS_PLOT:
        run_tests()
    else:
        plot_results()
