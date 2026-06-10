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

NX = 200
rho_BLUE = 1
rho_RED = rho_BLUE * 1000
nu_RED = 0.5
nu_BLUE = nu_RED / 100
mu_RED = rho_RED * nu_RED
mu_BLUE = rho_BLUE * nu_BLUE
v_left = 1e-4
v_right = -1e-2
p_left = v_left * rho_RED
p_right = v_right * rho_BLUE

A = np.array(
    [[0, 1, 0, 0], 
     [mu_RED, 0, -mu_BLUE, 0],
     [(NX-1)/2, 1, -(NX-1)/2, -1],
     [0, 0, (NX-1), 1]]
)

b = np.array([[v_left, 0, 0, v_right]]).T

A_inv = np.linalg.inv(A)
params = A_inv @ b

data_path = Path(".")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[-1], "r") as data_file:
    rho = np.array(data_file["hydro/rho"])
    v = np.array(data_file["hydro/v"])

x = np.arange(0, NX-0.5)

rho_analytical = np.empty(NX)
rho_analytical[:NX//2] = rho_RED
rho_analytical[NX//2:] = rho_BLUE

analytical = analytical_solution(x)/rho_analytical
Q = np.sqrt(np.sum((analytical-v[:,0,0])**2))/np.sqrt(np.sum(analytical**2))
print("Total deviation of numerical result from analytical result:", Q)

plt.plot(x, rho[:,0,0]*v[:,0,0], label="LBM")
plt.plot(x, analytical_solution(x), ls="--", c="red", label="analytical")
plt.xlabel("$x$")
plt.ylabel("$p_y$")
plt.ylim(min([p_left, p_right]) - 0.1*max([np.abs(p_left), np.abs(p_right)]), max([p_left, p_right]) + 0.1*max([np.abs(p_left), np.abs(p_right)]))
plt.xlim(0, NX-1)
plt.xticks([0, (NX-1)/2, NX-1])
plt.legend()
plt.savefig("Couette_momentum.pdf")
plt.close()

rho_analytical = np.empty(NX)
rho_analytical[:NX//2] = rho_RED
rho_analytical[NX//2:] = rho_BLUE
plt.axvline((NX-1)/2, ls="--", c="gray")
plt.plot(x, v[:,0,0], label="LBM, Q = %.3e"%Q)
plt.plot(x, analytical_solution(x)/rho_analytical, ls="--", c="red", label="analytical")
plt.xlabel("$x$")
plt.ylabel("$v$")
plt.ylim(min([v_left, v_right]) - 0.1*max([np.abs(v_left), np.abs(v_right)]), max([v_left, v_right]) + 0.1*max([np.abs(v_left), np.abs(v_right)]))
plt.xlim(0, NX-1)
plt.xticks([0, (NX-1)/2, NX-1])
plt.legend()
plt.savefig("Couette_velocity.pdf")
plt.close()
