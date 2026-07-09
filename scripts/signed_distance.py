import h5py
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
from skimage.measure import find_contours

@njit
def A_2(h, h_min=0.5):
    NX, NY = h.shape
    result = np.zeros(h.shape, dtype=np.float64)

    for i in range(NX):
        for j in range(NY):
            if h[i,j] < h_min:
                result[i,j] = 1.0
            else:
                result[i,j] = np.pow(h[i,j]-h_min+1.0, -3.0)
    return result

@njit
def A(h, h_min=0.5):
    if h < h_min:
        return 1.0
    else:
        return np.pow(h-h_min+1, -3)

@njit
def mod(i, n):
    if i < 0:
        return i + n
    elif i >= n:
        return i - n
    else:
        return i

@njit
def interpolate(rho_N, x, y):
    NX, NY = rho_N.shape
    i_min = int(np.floor(x - 0.5))
    i_max = i_min + 1
    j_min = int(np.floor(y - 0.5))
    j_max = j_min + 1

    rho_N_bottom = (i_max + 0.5 - x)*rho_N[mod(i_min, NX), mod(j_min, NY)] + (x - (i_min + 0.5))*rho_N[mod(i_max, NX), mod(j_min, NY)]
    rho_N_top = (i_max + 0.5 - x)*rho_N[mod(i_min, NX), mod(j_max, NY)] + (x - (i_min + 0.5))*rho_N[mod(i_max, NX), mod(j_max, NY)]

    return (j_max + 0.5 - y)*rho_N_bottom + (y - (j_min + 0.5))*rho_N_top

# @njit
def find_distance(rho_N, nx, ny, i_start, j_start, n_max=30):
    rho_N_i = rho_N[i_start, j_start]

    x_start = i_start + 0.5
    y_start = j_start + 0.5

    n = 1

    # Start in the droplet
    if rho_N_i > 0:
        while (n < n_max):
            x = x_start - n*nx
            y = y_start - n*ny
            rho_N_i = interpolate(rho_N, x, y)
            n += 1
            if (rho_N_i < 0):
                break
        else:
            return 1000

    while (n < n_max):
        x = x_start - n*nx
        y = y_start - n*ny
        rho_N_old = rho_N_i
        rho_N_i = interpolate(rho_N, x, y)
        if (rho_N_i > 0):
            x_old = x_start - (n-1)*nx
            y_old = y_start - (n-1)*ny
            x_final = rho_N_i / (rho_N_i - rho_N_old) * x_old - rho_N_old / (rho_N_i - rho_N_old) * x
            y_final = rho_N_i / (rho_N_i - rho_N_old) * y_old - rho_N_old / (rho_N_i - rho_N_old) * y
            return np.sqrt((x_start - x_final)**2 + (y_start - y_final)**2)
        n += 1

    return 1000
    

def main():
    data_path = Path(".")
    data_files = data_path.glob("data_*.h5")
    data_files = natsorted(data_files, key=lambda x: x.name)
    n_files = len(data_files)

    with h5py.File(data_files[-1], "r") as data_file:
        NX, NY, _ = data_file["cg/rho_N"].shape
        rho_N = data_file["cg/rho_N"][:,:,0]
        nx = data_file["cg/nx"][:,:,0]
        ny = data_file["cg/ny"][:,:,0]
        Gx = data_file["cg/Gx"][:,:,0]
        Gy = data_file["cg/Gy"][:,:,0]
        t = data_file["t"][()]

    G_norm = np.sqrt(Gx*Gx + Gy*Gy)

    delta_I = 0.5*G_norm
    distance = np.zeros((NX, NY), dtype=np.float64)

    for i in range(NX):
        for j in range(NY):
            distance[i,j] = find_distance(rho_N, nx[i,j], ny[i,j], i, j)

    distance_2 = distance.copy()
    distance_2[distance>999] = 0

    plt.imshow(distance_2.T, origin="lower", extent=[0, NX, 0, NY], vmin=np.min(distance))
    contours = find_contours(rho_N, 0)
    plt.xlim(150, 250)
    plt.ylim(150, 250)
    for contour in contours:
        contour += 0.5
        plt.plot(contour[:,0], contour[:,1], c="red")
    plt.colorbar(ticks=[np.min(distance), np.max(distance_2)])
    plt.savefig("distance.pdf")
    plt.close()

    force = delta_I*A_2(distance)
    plt.imshow(force.T, origin="lower", extent=[0, NX, 0, NY], vmin=0)
    contours = find_contours(rho_N, 0)
    plt.xlim(150, 250)
    plt.ylim(150, 250)
    for contour in contours:
        contour += 0.5
        plt.plot(contour[:,0], contour[:,1], c="red")
    plt.colorbar()
    plt.savefig("force.pdf")
    plt.close()

if __name__ == "__main__":
    main()