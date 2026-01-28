import pandas as pd
import matplotlib.pyplot as plt

try:
    df = pd.read_csv("../output_results/dispersion.csv")
    plt.figure(figsize=(10, 6))
    
    # Plot numerico
    plt.plot(df["time"], df["u_numerical"], 'r-', label="Numerical (Center)", linewidth=1)
    
    # Plot esatto (sarà 0 se non sei nello scenario 1, ma va bene)
    if df["u_exact"].abs().max() > 1e-10: 
        plt.plot(df["time"], df["u_exact"], 'k--', label="Exact", alpha=0.5)
        
    plt.xlabel("Time")
    plt.ylabel("u(0.5, 0.5)")
    plt.title("Wave Signal at Domain Center")
    plt.legend()
    plt.grid(True)
    plt.savefig("dispersion_plot.png")
    print("Dispersion plot saved to dispersion_plot.png")
except Exception as e:
    print(f"Error: {e}")
