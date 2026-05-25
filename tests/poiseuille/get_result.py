import h5py
from pathlib import Path
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def analytical_solution(y, H):
    mu = 0.1
    dpdx = -1e-5
    u_particular = -1/(2*mu)*dpdx*y*(H-y)
    return u_particular

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
x = np.arange(0, NX-0.5)
H = NX - 1

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