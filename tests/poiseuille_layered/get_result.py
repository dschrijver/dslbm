import h5py
from pathlib import Path
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.interpolate import interp1d

def analytical_solution(y, H):
    mu_1 = 0.05826 
    mu_2 = 0.05826 / 40
    G = 1.5e-8

    N = len(y)
    result = np.empty(N)

    for j in range(N):
        r = y[j]-0.5*H
        if (r < 0):
            result[j] = -G/(2*mu_1)*r**2 + (G*0.5*H*(mu_1 - mu_2))/(2*mu_1*(mu_1 + mu_2))*r + G*0.25*H**2/(mu_1 + mu_2)
        else:
            result[j] = -G/(2*mu_2)*r**2 + (G*0.5*H*(mu_1 - mu_2))/(2*mu_2*(mu_1 + mu_2))*r + G*0.25*H**2/(mu_1 + mu_2)
    return result

testcase_name = os.path.basename(os.path.dirname(__file__))

data_path = Path(".")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[-1], "r") as data_file:
    data = np.array(data_file["hydro/v"])
    t = np.array(data_file["t"])

NX, NY, NZ = data.shape
v_simulated = data[:, 0, 0]

# HWBB
x = np.arange(0.5, NX)
H = NX

v_analytical = analytical_solution(x, H)

Q = np.sqrt(np.sum((v_analytical-v_simulated)**2))/np.sqrt(np.sum(v_analytical**2))
print("Total deviation of numerical result from analytical result:", Q)

plt.plot(x, v_simulated)
plt.plot(x, v_analytical, ls="--", c="red", label="$Q = %.3e$"%Q)
plt.ylabel("$v$")
plt.xlabel("t")
plt.legend()
plt.savefig(testcase_name+".pdf")
plt.close()