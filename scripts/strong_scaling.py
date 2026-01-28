import subprocess
import re
import matplotlib.pyplot as plt
import os

# --- CONFIGURAZIONE ---
EXECUTABLE = "../build/wave_solver"
PARAM_FILE = "../include/parameters.prm"
TEMP_PARAM_FILE = "../include/parameters_temp.prm"

# Testiamo su 1, 2, 3, 4 core
PROCS_LIST = [1, 2, 3, 4] 

# Usiamo un livello alto (8) per dare abbastanza lavoro ai core
# Se Ref 8 è troppo lento sul tuo PC, cambialo a 7.
FIXED_REFINEMENT = 8  

def update_parameter_file(refinement_level):
    with open(PARAM_FILE, 'r') as f:
        content = f.read()
    
    # Sostituisce il Refinement level
    new_content = re.sub(
        r"(set Refinement level\s*=\s*)(\d+)",
        f"\\g<1>{refinement_level}",
        content
    )
    # Riduciamo il tempo finale per rendere il test più veloce
    new_content = re.sub(
        r"(set Final time\s*=\s*)(\d+\.?\d*)", 
        "set Final time       = 0.5", 
        new_content
    )

    with open(TEMP_PARAM_FILE, 'w') as f:
        f.write(new_content)

def parse_time(output):
    # Cerca la riga: "| Total wallclock time elapsed since start |      15.9s |"
    match = re.search(r"\| Total wallclock time elapsed since start\s*\|\s*([\d\.]+)s", output)
    if match:
        return float(match.group(1))
    return None

def run_strong_scaling():
    print(f"=== STRONG SCALING TEST (Refinement {FIXED_REFINEMENT}) ===")
    print("Obiettivo: Vedere se il tempo diminuisce aumentando i processori.\n")
    
    update_parameter_file(FIXED_REFINEMENT)
    
    times = []
    
    for p in PROCS_LIST:
        print(f"Running with {p} processors...", end=" ", flush=True)
        
        cmd = ["mpirun", "-n", str(p), EXECUTABLE, TEMP_PARAM_FILE]
        # Eseguiamo ignorando output standard per pulizia
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print("FAILED")
            print(result.stderr)
            times.append(None)
            continue
            
        elapsed = parse_time(result.stdout)
        if elapsed:
            print(f"Done. Time = {elapsed:.4f}s")
            times.append(elapsed)
        else:
            print("Error parsing time.")
            times.append(None)

    # Calcolo Speedup (T_1 / T_p)
    if times[0] is not None:
        t1 = times[0]
        speedups = []
        valid_procs = []
        for i, t in enumerate(times):
            if t is not None:
                speedups.append(t1 / t)
                valid_procs.append(PROCS_LIST[i])
        
        # Plot
        plt.figure(figsize=(8, 6))
        plt.plot(valid_procs, speedups, 'o-', label='Measured Speedup', color='blue')
        plt.plot(valid_procs, valid_procs, 'k--', label='Ideal Speedup (Linear)', alpha=0.5)
        plt.xlabel('Number of Processors')
        plt.ylabel('Speedup')
        plt.title(f'Strong Scaling (Refinement {FIXED_REFINEMENT})')
        plt.grid(True)
        plt.legend()
        plt.xticks(PROCS_LIST)
        plt.savefig("strong_scaling.png")
        print("\nGrafico salvato come: strong_scaling.png")
    
    # Cleanup
    if os.path.exists(TEMP_PARAM_FILE):
        os.remove(TEMP_PARAM_FILE)

if __name__ == "__main__":
    run_strong_scaling()
