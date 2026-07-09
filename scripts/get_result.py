import h5py
from pathlib import Path
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.interpolate import interp1d

def analytical_solution(y, H):
    mu_out = 0.05826 
    mu_in = 0.05826 / 40
    dpdx = -1.5e-8

    
    d = H/2
    a = H/4

    N = len(y)
    result = np.empty(N)

    for j in range(N):
        r = np.abs(y[j]-d)
        if (r < a):
            result[j] = -1/(2*mu_out)*dpdx*(d*d - a*a) - 1/(2*mu_in)*dpdx*(a*a - r*r)
        else:
            result[j] = -1/(2*mu_out)*dpdx*(d*d - r*r)
    return result

testcase_name = os.path.basename(os.path.dirname(__file__))

# HWBB
NX = 100
x_dslbm = np.arange(0.5, NX)
x_dslbm_extended = np.zeros(NX+2)
x_dslbm_extended[1:-1] = x_dslbm
x_dslbm_extended[-1] = NX
H = NX

v_theoretical_dslbm = analytical_solution(x_dslbm_extended, H)

data = np.genfromtxt("Case_C.csv", delimiter=",", skip_header=True)

x_theoretical = data[:,0]
v_theoretical = data[:,1]

# plt.plot(x_dslbm-50, v_theoretical_dslbm, label="theoretical W = 100")
# plt.plot(x_theoretical, v_theoretical, label="theoretical Saito")
# plt.ylabel("$v$")
# plt.xlabel("x")
# plt.axvline(-50, ls="--", c="gray")
# plt.axvline(50, ls="--", c="gray")
# plt.xlim(-55, 55)
# plt.legend()
# plt.savefig("theoretical_comparison.pdf")
# plt.close()

plt.scatter(x_dslbm_extended-50, v_theoretical_dslbm, label="theoretical W = 100", zorder=1)
plt.scatter(x_theoretical, v_theoretical, label="theoretical Saito", zorder=1)
plt.plot(x_dslbm_extended-50, v_theoretical_dslbm, zorder=1)
plt.plot(x_theoretical, v_theoretical, zorder=1)
plt.ylabel("$v$")
plt.xlabel("x")
plt.axvline(-50, ls="--", c="gray", zorder=0)
plt.axvline(50, ls="--", c="gray", zorder=0)
plt.axhline(0, ls="--", c="gray", zorder=0)
plt.xlim(-53, -45)
plt.ylim(-0.5e-5, 5e-5)
plt.legend()
plt.savefig("theoretical_comparison_zoomed.pdf")
plt.close()

x = data[:,0]
x2 = data[:,2]
v_simulated = data[:,-1]
v_analytical_saito = data[:,1]

plt.scatter(x_dslbm_extended-50, v_theoretical_dslbm, label="theoretical W = 100", zorder=1)
plt.scatter(x2, v_simulated, label="Saito Imp 6th", zorder=1)
plt.plot(x_dslbm_extended-50, v_theoretical_dslbm, zorder=1)
plt.plot(x2, v_simulated, zorder=1)
plt.ylabel("$v$")
plt.xlabel("x")
plt.axvline(-50, ls="--", c="gray", zorder=0)
plt.axvline(50, ls="--", c="gray", zorder=0)
plt.axhline(0, ls="--", c="gray", zorder=0)
plt.xlim(-53, -45)
plt.ylim(-0.5e-5, 5e-5)
plt.legend()
plt.savefig("data_saito_vs_our_analytical_x2_left_zoomed.pdf")
plt.close()

plt.plot(x_dslbm_extended-50, v_theoretical_dslbm, label="theoretical W = 100", zorder=1)
plt.plot(x2, v_simulated, label="Saito Imp 6th", zorder=1)
plt.ylabel("$v$")
plt.xlabel("x")
plt.axvline(-50, ls="--", c="gray", zorder=0)
plt.axvline(50, ls="--", c="gray", zorder=0)
plt.axhline(0, ls="--", c="gray", zorder=0)
plt.legend()
plt.savefig("data_saito_vs_our_analytical_x2.pdf")
plt.close()

plt.plot(x_dslbm_extended-50, v_theoretical_dslbm, label="theoretical W = 100", zorder=1)
plt.plot(np.arange(0.5, NX+1)-50, v_simulated, label="Saito Imp 6th", zorder=1)
plt.ylabel("$v$")
plt.xlabel("x")
plt.axvline(-50, ls="--", c="gray", zorder=0)
plt.axvline(50, ls="--", c="gray", zorder=0)
plt.axhline(0, ls="--", c="gray", zorder=0)
plt.legend()
plt.savefig("data_saito_vs_our_analytical_x.pdf")
plt.close()

plt.scatter(x_dslbm_extended-50, v_theoretical_dslbm, label="theoretical W = 100", zorder=1)
plt.scatter(np.arange(0.5, NX+1)-50, v_simulated, label="Saito Imp 6th", zorder=1)
plt.plot(x_dslbm_extended-50, v_theoretical_dslbm, zorder=1)
plt.plot(np.arange(0.5, NX+1)-50, v_simulated, zorder=1)
plt.ylabel("$v$")
plt.xlabel("x")
plt.axvline(-50, ls="--", c="gray", zorder=0)
plt.axvline(50, ls="--", c="gray", zorder=0)
plt.axhline(0, ls="--", c="gray", zorder=0)
plt.xlim(45, 53)
plt.ylim(-0.5e-5, 5e-5)
plt.legend()
plt.savefig("data_saito_vs_our_analytical_x_right_zoomed.pdf")
plt.close()

# plt.scatter(x_dslbm_extended-50, v_theoretical_dslbm, label="theoretical W = 100", zorder=1)
# plt.scatter(x2, v_simulated, label="Saito", zorder=1)
# plt.plot(x_dslbm_extended-50, v_theoretical_dslbm, zorder=1)
# plt.plot(x2, v_simulated, zorder=1)
# plt.ylabel("$v$")
# plt.xlabel("x")
# plt.axvline(-50, ls="--", c="gray", zorder=0)
# plt.axvline(50, ls="--", c="gray", zorder=0)
# plt.axhline(0, ls="--", c="gray", zorder=0)
# plt.xlim(-53, -45)
# plt.ylim(-0.5e-5, 5e-5)
# plt.legend()
# plt.savefig("data_saito_vs_our_analytical_zoomed.pdf")
# plt.close()