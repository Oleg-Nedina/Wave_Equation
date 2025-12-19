import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("build/output_energy/energy.csv")

t = df.iloc[:, 0].to_numpy()
E = df.iloc[:, 1].to_numpy()

idx0 = np.argmax(E > 0) 
E0 = E[idx0] if E[idx0] > 0 else 1.0

En = E / E0

fig, ax1 = plt.subplots()

# Energia non normalizzata (rosso)
ax1.plot(t, E, "r-", label="E(t)")
ax1.set_xlabel("t")
ax1.set_ylabel("E(t)", color="r")
ax1.tick_params(axis="y", labelcolor="r")
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(t, En, "b-", label="E(t)/E0")
ax2.set_ylabel("E(t)/E0", color="b")
ax2.tick_params(axis="y", labelcolor="b")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")

plt.tight_layout()
plt.show()

