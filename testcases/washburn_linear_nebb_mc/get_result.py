import h5py
from pathlib import Path
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, least_squares
from skimage import measure
import sys
import os

def analytical_solution(t, theta, L, W):
    sigma = 1e-2
    rho = 1
    tau = 0.8
    theta = theta / 360 * (2*np.pi)
    nu = 1/3*(tau - 0.5)
    mu = nu*rho
    V_cap = sigma/mu
    t_d = rho * W**2 / (12 * mu)
    return V_cap * W * np.cos(theta) / (6 * L) * t_d * (np.exp(-t / t_d) + t / t_d - 1)

def circle_fit(p, x, y):
    R, x0, y0 = p
    return (x - x0)**2 + (y - y0)**2 - R**2

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

t_list = np.empty(n_files, dtype=np.int32)
h_list = np.empty(n_files, dtype=np.float64)
h_list_stored = np.empty(n_files, dtype=np.float64)
theta_list = np.zeros(n_files, dtype=np.float64)

for j in range(n_files):
    with h5py.File(data_files[j], "r") as data_file:
        rho_RED = np.array(data_file["rho_RED"])
        rho_BLUE = np.array(data_file["rho_BLUE"])
        pressure = np.array(data_file["pressure"])
        t = np.array(data_file["t"])
    t_list[j] = t
    phi = rho_RED - rho_BLUE
    NX, NY, NZ = phi.shape

    y = np.arange(0.5, NY)

    h_list[j] = 0
    for x in range(NX):
        interpolated = interp1d(phi[x, :, 0], y)
        h_list[j] += interpolated(0)
    h_list[j] /= NX

    if not store_data:
        with h5py.File(stored_files[j], "r") as stored_file:
            rho_RED_stored = np.array(stored_file["rho_RED"])
            rho_BLUE_stored = np.array(stored_file["rho_BLUE"])
        phi_stored = rho_RED_stored - rho_BLUE_stored
        h_list_stored[j] = 0
        for x in range(NX):
            interpolated = interp1d(phi_stored[x, :, 0], y)
            h_list_stored[j] += interpolated(0)
        h_list_stored[j] /= NX

    contour = measure.find_contours(phi[:,:,0], 0)[0]
    if j > 0:
        R, x0, y0 = least_squares(circle_fit, (NX/2, NX/2, h_list[j]), args=(contour[:,0], contour[:,1] + 0.5)).x
        L = x0
        H = R - np.sqrt(R*R - L*L)
        theta = np.arctan((R - H) / L) + np.pi / 2
        theta = theta / (2*np.pi) * 360 - 90
        theta_list[j] = theta
    else:
        theta = 90

    theta_list[j] = theta
    print(theta)

h_list[0] = 75
h_list_stored[0] = 75
length = NY
width = NX-1
params, _ = curve_fit(lambda t, theta: analytical_solution(t, theta, length, width), t_list, h_list - h_list[0], 80)
fitted_angle = params[0]
measured_angle = np.mean(theta_list[10:])

h_analytical = analytical_solution(t_list, measured_angle, length, width) + h_list[0]

Q = np.sqrt(np.sum((h_analytical-h_list)**2))/np.sqrt(np.sum(h_analytical**2))
print("Total deviation of numerical result from analytical result:", Q)

if store_data:
    os.system("cp ../../*.h5 stored")
    os.system("cd stored; rename 's/.h5/.dat/' *.h5")
    print("\033[34mStoring data...\033[0m")
else:
    Q_stored = np.sqrt(np.sum((h_analytical-h_list_stored)**2))/np.sqrt(np.sum(h_analytical**2))
    print("Expected total deviation of numerical result from analytical result:", Q_stored)

    tolerance = 1e-8
    if (Q > Q_stored + tolerance*Q_stored): 
        print("\033[91m" + "Test Failed! Error is more than "+str(tolerance)+" percent bigger than expected!" + "\033[0m")
        sys.exit(1)
    else:
        print("\033[92mTest Passed!\033[0m")

plt.plot(t_list, h_list, lw=2, label="Measured height")
plt.plot(t_list, analytical_solution(t_list, fitted_angle, length, width) + h_list[0], ls="--", c="red", label="Analytical with fitted contact angle")
plt.plot(t_list, analytical_solution(t_list, measured_angle, length, width) + h_list[0], ls="--", c="green", label="Analytical with measured contact angle")
plt.xlabel("t")
plt.ylabel("h")
plt.legend()
plt.savefig("washburn_heights_L=%d.png"%int(length))
plt.close()

plt.plot(t_list[10:], theta_list[10:], label="Measured contact angle")
plt.axhline(fitted_angle, ls="--", c="red", label="$\\theta_{\\text{fit}} = %.2f$"%(fitted_angle))
plt.legend()
plt.xlabel("t")
plt.ylabel("$\\theta_c$")
plt.suptitle("Final contact angle difference: %.2f degrees"%np.abs(measured_angle - fitted_angle))
plt.savefig("washburn_contact_angles_L=%d.png"%int(length))
plt.close()