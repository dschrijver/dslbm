from read_data import read_data
import matplotlib.pyplot as plt
import os
import h5py
import numpy as np

data_files = read_data()

if not os.path.isdir("frames"):
    os.mkdir("frames")

dpi = 100
width = 1200
height = 900

fig = plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)

ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.set_xlabel("x")
ax.set_ylabel("y")

with h5py.File(data_files[0], "r") as data_file:
    NX, NY, NZ = data_file["main/rho"].shape

placeholder = np.zeros((NX, NY))
im = ax.imshow(placeholder.T, origin="lower", vmin=0.425, vmax=0.575)
cbar = plt.colorbar(im, ax=ax, ticks=[0.425, 0.5, 0.575], label="$\\phi$")

for n in range(len(data_files)):
    with h5py.File(data_files[n], "r") as data_file:
        phi = data_file["main/phi"][:,:,0]
    im.set_data(phi.T)
    fig.savefig("frames/%5d.jpg"%(n+1), dpi=dpi)
    print(n)

os.system("ffmpeg -framerate 50 -pattern_type glob -i 'frames/*.jpg' -c:v libx264 -pix_fmt yuv420p out.mp4")