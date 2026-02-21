import subprocess
import re
import matplotlib.pyplot as plt
import os

# --- CONFIGURAZIONE ---
EXECUTABLE = "../build/wave_solver"
PARAM_FILE = "../include/parameters.prm"
TEMP_PARAM_FILE = "../include/parameters_temp.prm"

# COPPIE (Processori, Refinement Level)
# Nota: Saltiamo 2 e 3 perché in 2D la griglia quadruplica, non raddoppia.
SCALING_CASES = [(1, 7), (4, 8)]


def update_parameter_file(refinement_level):
    with open(PARAM_FILE, "r") as f:
        content = f.read()

    new_content = re.sub(
        r"(set Refinement level\s*=\s*)(\d+)", f"\\g<1>{refinement_level}", content
    )
    # Tempo breve per il test
    new_content = re.sub(
        r"(set Final time\s*=\s*)(\d+\.?\d*)", "set Final time       = 0.5", new_content
    )

    with open(TEMP_PARAM_FILE, "w") as f:
        f.write(new_content)


def parse_time(output):
    match = re.search(
        r"\| Total wallclock time elapsed since start\s*\|\s*([\d\.]+)s", output
    )
    if match:
        return float(match.group(1))
    return None


def run_weak_scaling():
    print("=== WEAK SCALING TEST ===")
    print(
        "Obiettivo: Vedere se il tempo rimane costante aumentando problema e risorse.\n"
    )
    print(
        "Nota: Eseguiamo solo 1 e 4 core perché il raffinamento quadruplica la griglia."
    )

    times = []
    procs = []

    for p, ref in SCALING_CASES:
        print(f"Running with {p} processors (Refinement {ref})...", end=" ", flush=True)

        update_parameter_file(ref)

        cmd = ["mpirun", "-n", str(p), EXECUTABLE, TEMP_PARAM_FILE]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print("FAILED")
            continue

        elapsed = parse_time(result.stdout)
        if elapsed:
            print(f"Done. Time = {elapsed:.4f}s")
            times.append(elapsed)
            procs.append(p)
        else:
            print("Error parsing time.")

    # Calcolo Efficienza (Idealmente dovrebbe essere 1.0)
    if len(times) > 0:
        t1 = times[0]
        efficiencies = [
            t1 / t for t in times
        ]  # Efficiency = Time(1) / Time(N) in weak scaling ideal scenario

        plt.figure(figsize=(8, 6))
        plt.plot(procs, efficiencies, "o-", color="green", label="Parallel Efficiency")
        plt.axhline(y=1.0, color="k", linestyle="--", label="Ideal Efficiency")
        plt.xlabel("Number of Processors")
        plt.ylabel("Efficiency")
        plt.title("Weak Scaling")
        plt.ylim(0, 1.2)  # L'efficienza non supera mai 1 di molto
        plt.grid(True)
        plt.legend()
        plt.xticks([1, 4])
        plt.savefig("weak_scaling.png")
        print("\nGrafico salvato come: weak_scaling.png")

    if os.path.exists(TEMP_PARAM_FILE):
        os.remove(TEMP_PARAM_FILE)


if __name__ == "__main__":
    run_weak_scaling()
