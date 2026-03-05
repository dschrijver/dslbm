import h5py
import numpy as np
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from skimage import measure
from scipy.optimize import least_squares
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble="\\usepackage{amsmath}")

def residuals(params, x, y):
    xc, yc, R = params
    return np.sqrt((x - xc)**2 + (y - yc)**2) - R

store_data = False
if not os.path.isdir("stored"):
    os.makedirs("stored")
    store_data = True

testcase_name = os.path.basename(os.path.dirname(__file__))

data_path = Path("../../")
data_files = data_path.glob("data_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[0], "r", swmr=True) as data_file:
    rho_RED = np.array(data_file["rho_RED"])
    rho_BLUE = np.array(data_file["rho_BLUE"])
NX, NY, NZ = rho_RED.shape
rho_lava_0 = np.max(rho_RED)
rho_air_0 = np.max(rho_BLUE)

if not store_data:
    stored_path = Path("stored")
    stored_files = stored_path.glob("data_*.dat")
    stored_files = natsorted(stored_files, key=lambda x: x.name)

x0 = NX/2
y0 = 0
R0 = 40
params0 = [x0, y0, R0]

theta_deg_list = np.empty(n_files)
t_list = np.empty(n_files)

for n in range(n_files):
    with h5py.File(data_files[n], "r", swmr=True) as data_file:
        rho_RED = np.array(data_file["rho_RED"])
        rho_BLUE = np.array(data_file["rho_BLUE"])
        u = np.array(data_file["u"])
        v = np.array(data_file["v"])
        w = np.array(data_file["w"])
        t_list[n] = np.array(data_file["t"])
    phi = rho_RED - rho_BLUE
    vel = np.sqrt(u*u + v*v + w*w)

    contour = measure.find_contours(phi[:,:,0], 0.5*(rho_lava_0 - rho_air_0))[0]
    contour[:,0] += 0.5
    res = least_squares(residuals, params0, args=(contour[:,0], contour[:,1]))
    xc, yc, R = res.x
    H = R + yc
    L = 2*np.sqrt(R**2 - yc**2)
    theta = np.arctan((R-H)/(0.5*L)) + 0.5*np.pi
    theta_deg = theta/(2*np.pi)*360
    theta_deg_list[n] = theta_deg

print("Computed angle:", theta_deg, "degrees")

if not store_data:
    with h5py.File(stored_files[-1], "r", swmr=True) as stored_file:
        rho_RED_stored = np.array(stored_file["rho_RED"])
        rho_BLUE_stored = np.array(stored_file["rho_BLUE"])
    phi_stored = rho_RED_stored - rho_BLUE_stored
    contour = measure.find_contours(phi_stored[:,:,0], 0.5*(rho_lava_0 - rho_air_0))[0]
    contour[:,0] += 0.5
    res = least_squares(residuals, params0, args=(contour[:,0], contour[:,1]))
    xc, yc, R = res.x
    H = R + yc
    L = 2*np.sqrt(R**2 - yc**2)
    theta = np.arctan((R-H)/(0.5*L)) + 0.5*np.pi
    theta_deg_stored = theta/(2*np.pi)*360
    print("Expected angle:", theta_deg_stored, "degrees")

if store_data:
    os.system("cp ../../*.h5 stored")
    os.system("cd stored; rename 's/.h5/.dat/' *.h5")
    print("\033[34mStoring data...\033[0m")


fig, ax = plt.subplots(figsize=(7, 6.0))
im = ax.imshow(vel[:,:,0].T, origin="lower", extent=[0, NX, -0.5, NY-0.5], cmap="RdBu")
ax.set_xticks([0, NX])
ax.set_yticks([0, NY-1])
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_xlim([0, NX])
ax.set_ylim([0, NY-1])
ax.set_title("$|\\boldsymbol{u}|_{\\textrm{max}} = %.5e$"%(np.max(np.sqrt(u*u + v*v + w*w))))
fig.colorbar(im)
plt.savefig("spurious_currents.pdf")
plt.close()

fig, ax = plt.subplots(figsize=(7, 6.0))
ax.plot(t_list[1:], theta_deg_list[1:])
ax.axhline(160, ls="--", color="black")
ax.set_xlabel("$t$")
ax.set_ylabel("$\\theta_c$")
ax.set_title("$\\theta_c^\\text{final} = %.3f^{\\circ}$"%(theta_deg_list[-1]))
plt.savefig("contact_angles.pdf")
plt.close()

