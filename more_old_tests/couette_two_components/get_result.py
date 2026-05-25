import numpy as np
import matplotlib.pyplot as plt
import h5py
from pathlib import Path
from natsort import natsorted
import sys
import os


def analytical_solution(x):
    a_RED, b_RED, a_BLUE, b_BLUE = params[:,0]

    result = np.zeros(x.shape)

    for i in range(len(x)):
        x_i = x[i]
        if x_i < (NX-1)/2:
            result[i] = rho_RED*(a_RED*x_i + b_RED)
        else:
            result[i] = rho_BLUE*(a_BLUE*x_i + b_BLUE)
    
    return result

store_data = False

if not os.path.isdir("stored"):
    os.makedirs("stored")
    store_data = True

NX = 200
rho_RED = 1
rho_BLUE = 0.5
nu_RED = 1/6
nu_BLUE = 1/6
mu_RED = rho_RED * nu_RED
mu_BLUE = rho_BLUE * nu_BLUE

A = np.array(
    [[0, 1, 0, 0], 
     [mu_RED, 0, -mu_BLUE, 0],
     [(NX-1)/2, 1, -(NX-1)/2, -1],
     [0, 0, (NX-1), 1]]
)

b = np.array([[1e-2, 0, 0, 0]]).T

A_inv = np.linalg.inv(A)
params = A_inv @ b

data_path = Path("../../")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

if not store_data:
    stored_path = Path("stored")
    stored_files = stored_path.glob("data_*.dat")
    stored_files = natsorted(stored_files, key=lambda x: x.name)

with h5py.File(data_files[-1], "r") as data_file:
    rho = np.array(data_file["rho"])
    v = np.array(data_file["v"])

if not store_data:
    with h5py.File(stored_files[-1], "r") as stored_file:
        stored = np.array(stored_file["v"])

x = np.arange(0, NX-0.5)

rho_analytical = np.empty(NX)
rho_analytical[:NX//2] = rho_RED
rho_analytical[NX//2:] = rho_BLUE

analytical = analytical_solution(x)/rho_analytical
Q = np.sqrt(np.sum((analytical-v[:,0,0])**2))/np.sqrt(np.sum(analytical**2))
print("Total deviation of numerical result from analytical result:", Q)

if store_data:
    os.system("cp ../../*.h5 " + "stored")
    os.system("cd "+"stored"+"; rename 's/.h5/.dat/' *.h5")
    print("\033[34mStoring data...\033[0m")
else:
    Q_stored = np.sqrt(np.sum((analytical-stored[:,0,0])**2))/np.sqrt(np.sum(analytical**2))

    print("Expected total deviation of numerical result from analytical result:", Q_stored)

    tolerance = 1e-8
    if (Q > Q_stored + tolerance*Q_stored): 
        print("\033[91m" + "Test Failed! Error is more than "+str(tolerance)+" percent bigger than expected!" + "\033[0m")
        sys.exit(1)
    else:
        print("\033[92mTest Passed!\033[0m")

plt.plot(x, rho[:,0,0]*v[:,0,0], label="LBM")
plt.plot(x, analytical_solution(x), ls="--", c="red", label="analytical")
plt.xlabel("$x$")
plt.ylabel("$p_y$")
plt.ylim(0, 0.01)
plt.xlim(0, NX-1)
plt.xticks([0, (NX-1)/2, NX-1])
plt.legend()
plt.savefig("Couette_momentum.png")
plt.close()

rho_analytical = np.empty(NX)
rho_analytical[:NX//2] = rho_RED
rho_analytical[NX//2:] = rho_BLUE
plt.axvline((NX-1)/2, ls="--", c="gray")
plt.plot(x, v[:,0,0], label="LBM, Q = %.3e"%Q)
plt.plot(x, analytical_solution(x)/rho_analytical, ls="--", c="red", label="analytical")
plt.xlabel("$x$")
plt.ylabel("$v$")
plt.ylim(0, 0.01)
plt.xlim(0, NX-1)
plt.xticks([0, (NX-1)/2, NX-1])
plt.legend()
plt.savefig("Couette_velocity.png")
plt.close()

print("Generated Couette_momentum.png and Couette_velocity.png!")
