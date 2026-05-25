import h5py
import numpy as np
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import os
import sys

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble="\\usepackage{amsmath}")

store_data = False
if not os.path.isdir("stored"):
    os.makedirs("stored")
    store_data = True

data_path = Path("../../")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

if not store_data:
    stored_path = Path("stored")
    stored_files = stored_path.glob("data_*.dat")
    stored_files = natsorted(stored_files, key=lambda x: x.name)

with h5py.File(data_files[0], "r", swmr=True) as data_file:
    rho_ini = np.array(data_file["rho"])[:,:,0]

with h5py.File(data_files[-1], "r", swmr=True) as data_file:
    rho_fin = np.array(data_file["rho"])[:,:,0]
    u = np.array(data_file["u"])[:,:,0]
    v = np.array(data_file["v"])[:,:,0]

if not store_data:
    with h5py.File(stored_files[-1], "r", swmr=True) as stored_file:
        u_stored = np.array(stored_file["u"])[:,:,0]
        v_stored = np.array(stored_file["v"])[:,:,0]

NX, NY = u.shape
x = np.arange(0.5, NX, dtype=np.float64)
y = np.arange(0, NY-0.5, dtype=np.float64)

L = NX

y_exact = np.array([1, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0])*L
u_exact = np.array([1, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0])*u[NX//2, -1]

x_exact = np.array([1, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0])*L
v_exact = np.array([0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000])*u[NX//2, -1]

u_interp = CubicSpline(y, u[NX//2,:])
u_interp = u_interp(y_exact)

if not store_data:
    u_stored_interp = CubicSpline(y, u_stored[NX//2,:])
    u_stored_interp = u_stored_interp(y_exact)

Q = np.sqrt(np.sum((u_exact-u_interp)**2))/np.sqrt(np.sum(u_exact**2))
print("Total deviation of numerical result from analytical result:", Q)

if store_data:
    os.system("cp ../../*.h5 stored")
    os.system("cd stored; rename 's/.h5/.dat/' *.h5")
    print("\033[34mStoring data...\033[0m")
else:
    Q_stored = np.sqrt(np.sum((u_exact-u_stored_interp)**2))/np.sqrt(np.sum(u_exact**2))
    print("Expected total deviation of numerical result from analytical result:", Q_stored)

    tolerance = 1e-8
    if (Q > Q_stored + tolerance*Q_stored): 
        print("\033[91m" + "Test Failed! Error is more than "+str(tolerance)+" percent bigger than expected!" + "\033[0m")
        sys.exit(1)
    else:
        print("\033[92mTest Passed!\033[0m")

M_ini = np.sum(rho_ini)
M_fin = np.sum(rho_fin)
delta_M = M_fin - M_ini

M_ini_NS = 0.5*np.sum(rho_ini[:,0]) + np.sum(rho_ini[:,1:-1]) + 0.5*np.sum(rho_ini[:,-1])
M_fin_NS = 0.5*np.sum(rho_fin[:,0]) + np.sum(rho_fin[:,1:-1]) + 0.5*np.sum(rho_fin[:,-1])
delta_M_NS = M_fin_NS - M_ini_NS

fig, ax = plt.subplots(figsize=(7, 6.0))
ax.plot(y, u[NX//2,:], label="LBM", zorder=0)
ax.scatter(y_exact, u_exact, label="Ghia", color="red", zorder=1)
ax.set_xlabel("$y$")
ax.set_ylabel("$u$")
plt.legend()
plt.suptitle("$Q = $ %.3e, "%Q+"$\\Delta M^\\text{kin} = $"+" %.3e"%(delta_M/M_ini*100) + "$\\%$, $\\Delta M^\\text{hyd} = $"+" %.3e"%(delta_M_NS/M_ini_NS*100) + "$\\%$")
plt.savefig("result_1.pdf")
plt.close()

fig, ax = plt.subplots(figsize=(7, 6.0))
ax.streamplot(x, y, u.T, v.T, density=2, color="white", linewidth=0.5)
im = ax.imshow(np.sqrt(u*u+v*v).T, origin="lower", cmap="jet")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
plt.colorbar(im)
plt.savefig("figure.pdf")
plt.close()
