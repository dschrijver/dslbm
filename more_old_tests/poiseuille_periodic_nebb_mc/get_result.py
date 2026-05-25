import h5py
from pathlib import Path
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def analytical_solution(y, H):
    nu = 0.1
    rho = 1
    dpdx = -1e-5
    return -1/(2*nu*rho)*dpdx*y*(H-y)

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

with h5py.File(data_files[-1], "r") as data_file:
    data = np.array(data_file["u"])
NX, NY, NZ = data.shape
data = data[NX//2, :, 0]

if not store_data:
    with h5py.File(stored_files[-1], "r", swmr=True) as stored_file:
        data_stored = np.array(stored_file["u"])
    data_stored = data_stored[NX//2, :, 0]

# WETNODE
y = np.arange(0, NY-0.5)
H = NY-1

u_simulated = data
u_analytical = analytical_solution(y, H)

if not store_data:
    u_stored = data_stored
Q = np.sqrt(np.sum((u_analytical-u_simulated)**2))/np.sqrt(np.sum(u_analytical**2))
print("Total deviation of numerical result from analytical result:", Q)

if store_data:
    os.system("cp ../../*.h5 stored")
    os.system("cd stored; rename 's/.h5/.dat/' *.h5")
    print("\033[34mStoring data...\033[0m")
else:
    Q_stored = np.sqrt(np.sum((u_analytical-u_stored)**2))/np.sqrt(np.sum(u_analytical**2))
    print("Expected total deviation of numerical result from analytical result:", Q_stored)

    tolerance = 1e-8
    if (Q > Q_stored + tolerance*Q_stored): 
        print("\033[91m" + "Test Failed! Error is more than "+str(tolerance)+" percent bigger than expected!" + "\033[0m")
        sys.exit(1)
    else:
        print("\033[92mTest Passed!\033[0m")

y_analytical = np.linspace(0, H, 100)

plt.plot(analytical_solution(y_analytical, H), y_analytical, ls="--", c="red", zorder=0)
plt.scatter(u_simulated, y, label="Q = %.3e"%Q, zorder=1)
plt.ylabel("$y$")
plt.xlabel("$u$")
plt.legend()
plt.savefig(testcase_name+".png")
plt.close()
