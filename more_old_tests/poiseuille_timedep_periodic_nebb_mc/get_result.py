import h5py
from pathlib import Path
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def analytical_solution(y, t, H, n_max=1000):
    nu = 0.1
    rho = 1000
    dpdx = -1e-5
    u_particular = -1/(2*nu*rho)*dpdx*y*(H-y)
    u_homogeneous = 0
    for n in range(1, n_max+1):
        A_n = 1/(2*nu*rho)*dpdx*(-2*(-2 + 2*np.cos(n*np.pi) + n*np.pi*np.sin(n*np.pi))/((n*np.pi)**3))
        u_homogeneous += A_n * H**2 * np.exp(-n**2 * np.pi**2 / H**2 * nu * t)*np.sin(n*np.pi*y/H)
    return u_particular + u_homogeneous

store_data = False
if not os.path.isdir("stored"):
    os.makedirs("stored")
    store_data = True

testcase_name = os.path.basename(os.path.dirname(__file__))

data_path = Path("../../")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

if not store_data:
    stored_path = Path("stored")
    stored_files = stored_path.glob("data_*.dat")
    stored_files = natsorted(stored_files, key=lambda x: x.name)

u_max_simulated = np.empty(n_files, dtype=np.float64)
u_max_analytical = np.empty(n_files, dtype=np.float64)
t_list = np.empty(n_files, dtype=np.int32)

if not store_data:
    u_max_stored = np.empty(n_files, dtype=np.float64)

for n in range(n_files):
    with h5py.File(data_files[n], "r") as data_file:
        data = np.array(data_file["u"])
        t = np.array(data_file["t"])
    t_list[n] = t
    NX, NY, NZ = data.shape
    data = data[NX//2, :, 0]

    if not store_data:
        with h5py.File(stored_files[n], "r") as stored_file:
            data_stored = np.array(stored_file["u"])
        data_stored = data_stored[NX//2, :, 0]

    # WETNODE
    y = np.arange(0, NY-0.5)
    H = NY-1

    u_max_simulated[n] = np.max(data)
    u_max_analytical[n] = np.max(analytical_solution(y, t, H))

    if not store_data:
        u_max_stored[n] = np.max(data_stored)

Q = np.sqrt(np.sum((u_max_analytical-u_max_simulated)**2))/np.sqrt(np.sum(u_max_analytical**2))
print("Total deviation of numerical result from analytical result:", Q)

if store_data:
    os.system("cp ../../*.h5 stored")
    os.system("cd stored; rename 's/.h5/.dat/' *.h5")
    print("\033[34mStoring data...\033[0m")
else:
    Q_stored = np.sqrt(np.sum((u_max_analytical-u_max_stored)**2))/np.sqrt(np.sum(u_max_analytical**2))
    print("Expected total deviation of numerical result from analytical result:", Q_stored)

    tolerance = 1e-8
    if (Q > Q_stored + tolerance*Q_stored): 
        print("\033[91m" + "Test Failed! Error is more than "+str(tolerance)+" percent bigger than expected!" + "\033[0m")
        sys.exit(1)
    else:
        print("\033[92mTest Passed!\033[0m")

plt.plot(t_list, u_max_simulated)
plt.plot(t_list, u_max_analytical, ls="--", c="red")
plt.ylabel("$u_{\\text{max}}$")
plt.xlabel("t")
plt.savefig(testcase_name+".png")
plt.close()

plt.plot(t_list[1:], ((u_max_simulated - u_max_analytical)/u_max_analytical * 100)[1:])
plt.xlabel("t")
plt.ylabel("Error (%)")
plt.suptitle("Error in simulated velocity compared to analytical velocity")
plt.savefig(testcase_name+"_error"+".png")
plt.close()
