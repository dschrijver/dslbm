import h5py
from pathlib import Path
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def analytical_solution(y, t, H, n_max=1000):
    mu = 0.1
    rho = 1
    dpdx = -1e-5
    u_particular = -1/(2*mu)*dpdx*y*(H-y)
    u_homogeneous = 0
    for n in range(1, n_max+1):
        A_n = 1/(2*mu)*dpdx*(-2*(-2 + 2*np.cos(n*np.pi) + n*np.pi*np.sin(n*np.pi))/((n*np.pi)**3))
        u_homogeneous += A_n * H**2 * np.exp(-n**2 * np.pi**2 / H**2 * mu/rho * t)*np.sin(n*np.pi*y/H)
    return u_particular + u_homogeneous

testcase_name = os.path.basename(os.path.dirname(__file__))

data_path = Path(".")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

u_max_simulated = np.empty(n_files, dtype=np.float64)
u_max_analytical = np.empty(n_files, dtype=np.float64)
t_list = np.empty(n_files, dtype=np.int32)

for n in range(n_files):
    with h5py.File(data_files[n], "r") as data_file:
        data = np.array(data_file["hydro/v"])
        t = np.array(data_file["t"])
    t_list[n] = t
    NX, NY, NZ = data.shape
    data = data[:, NY//2, 0]

    # NEBB
    x = np.arange(0, NX-0.5)
    H = NX - 1

    u_max_simulated[n] = np.max(data)
    u_max_analytical[n] = np.max(analytical_solution(x, t, H))

Q = np.sqrt(np.sum((u_max_analytical-u_max_simulated)**2))/np.sqrt(np.sum(u_max_analytical**2))
print("Total deviation of numerical result from analytical result:", Q)

plt.plot(t_list, u_max_simulated)
plt.plot(t_list, u_max_analytical, ls="--", c="red")
plt.ylabel("$u_{\\text{max}}$")
plt.xlabel("t")
plt.savefig(testcase_name+".pdf")
plt.close()

plt.plot(t_list[1:], ((u_max_simulated - u_max_analytical)/u_max_analytical * 100)[1:])
plt.xlabel("t")
plt.ylabel("Error (%)")
plt.suptitle("Error in simulated velocity compared to analytical velocity")
plt.savefig(testcase_name+"_error"+".pdf")
plt.close()
