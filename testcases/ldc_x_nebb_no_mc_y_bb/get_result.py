import h5py
import numpy as np
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import root, minimize
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble="\\usepackage{amsmath}")

data_path = Path("../../")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)

n_files = len(data_files)

with h5py.File(data_files[-1], "r", swmr=True) as data_file:
    u = np.array(data_file["u"])[:,:,0]
    v = np.array(data_file["v"])[:,:,0]

NX, NY = u.shape
x = np.arange(0, NX-0.5, dtype=np.float64)
y = np.arange(0.5, NY, dtype=np.float64)

y_exact = np.array([1, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0])*y[-1]
u_exact = np.array([1, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0])*u[NX//2, -1]

fig, ax = plt.subplots(figsize=(7, 6.0))
ax.plot(y, u[NX//2,:], label="LBM", zorder=0)
ax.scatter(y_exact, u_exact, label="Ghia", color="red", zorder=1)
ax.set_xlabel("$y$")
ax.set_ylabel("$u$")
plt.legend()
plt.savefig("result_1.pdf")
plt.close()

fig, ax = plt.subplots(figsize=(7, 6.0))
ax.streamplot(x, y, u.T, v.T, density=2, color="white", linewidth=0.5)
im = ax.imshow(np.sqrt(u*u+v*v).T, origin="lower", cmap="jet")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
plt.colorbar(im)
plt.savefig("result_2.pdf")
plt.close()
