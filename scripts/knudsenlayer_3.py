import numpy as np
from scipy.optimize import fsolve
from numba import njit
import matplotlib.pyplot as plt


def correction_boundary(i, params_Kn, geq):
    rho_Kn = params_Kn[0]
    u_Kn = params_Kn[1]
    v_Kn = params_Kn[2]
    w_Kn = params_Kn[3]
    
    rho = rho_Kn + rho_NS
    u = 1/rho*(rho_Kn*u_Kn + rho_NS*u_NS)
    v = 1/rho*(rho_Kn*v_Kn + rho_NS*v_NS)
    w = 1/rho*(rho_Kn*w_Kn + rho_NS*w_NS)

    rho = rho[0]

    vel = np.empty(3, dtype=np.float64)
    vel[0] = u[0]
    vel[1] = v[0]
    vel[2] = w[0]
    c_double = c.astype(np.float64)

    uc = np.dot(vel, c_double[i,:])

    Nzh = np.zeros(3, dtype=np.float64)
    for alpha in range(3):
        if n_hat[alpha] > 0:
            Nzh[alpha] = 1/cs2 * Force[alpha,0]
        else:
            A = 1/3
            sum_parallel = 0
            for i2 in range(NP):
                cn = c[i2,0]*n_hat[0] + c[i2,1]*n_hat[1] + c[i2,2]*n_hat[2]
                if cn == 0:
                    sum_parallel += c_double[i2,alpha]*geq[i2,0]
            Nzh[alpha] = 2/(cs2*A) * (sum_parallel + (A - 1)*rho*vel[alpha] + 1/2*Force[alpha,0])

    return 2*wp[i]*rho*uc/cs2 - wp[i] * np.dot(c_double[i,:], Nzh)


def compute_eq(i, rho, vel):
    uc = np.zeros(N, dtype=np.float64)
    u2 = np.zeros(N, dtype=np.float64)
    c_double = c.astype(np.float64)
    for alpha in range(3):
        uc += vel[alpha]*c_double[i,alpha]
        u2 += vel[alpha]*vel[alpha]
    return wp[i] * rho * (1 + uc/cs2 + uc*uc/(2*cs2*cs2) - u2/(2*cs2))


def system_of_equations(f):

    # Compute macroscopic quantities
    rho_Kn = np.sum(f, axis=0)
    u_Kn = np.zeros(N, dtype=np.float64)
    v_Kn = np.zeros(N, dtype=np.float64)
    w_Kn = np.zeros(N, dtype=np.float64)
    c_double = c.astype(np.float64)

    for j in range(N):
        if rho_Kn[j] > 1e-15:
            for i in range(NP):
                u_Kn[j] += c_double[i,0] * f[i,j]
                v_Kn[j] += c_double[i,1] * f[i,j]
                w_Kn[j] += c_double[i,2] * f[i,j]
            u_Kn[j] /= rho_Kn[j]
            v_Kn[j] /= rho_Kn[j]
            w_Kn[j] /= rho_Kn[j]

    rho = rho_Kn + rho_NS
    u = 1/rho * (rho_Kn*u_Kn + rho_NS*u_NS)
    v = 1/rho * (rho_Kn*v_Kn + rho_NS*v_NS)
    w = 1/rho * (rho_Kn*w_Kn + rho_NS*w_NS)
    vel = np.empty((3, N), dtype=np.float64)
    vel[0] = u
    vel[1] = v
    vel[2] = w

    params_Kn = np.zeros((4, N), dtype=np.float64)
    params_Kn[0] = rho_Kn
    params_Kn[1] = u_Kn
    params_Kn[2] = v_Kn
    params_Kn[3] = w_Kn

    # Compute geq
    geq = np.zeros((NP, N), dtype=np.float64)
    for i in range(NP):
        geq[i] = compute_eq(i, rho, vel) - compute_eq(i, rho_NS, vel_NS)
    
    f_result = np.zeros((NP, N), dtype=np.float64)
    for i in range(NP):
        cn = c[i,0]*n_hat[0] + c[i,1]*n_hat[1] + c[i,2]*n_hat[2]

        # Zero-velocity and parallel populations
        if cn == 0:
            f_result[i] = geq[i]

        # Outgoing populations
        elif cn < 0:
            for j in range(N):
                f_result[i,j] = 0
                for k in range(j+1,N):
                    f_result[i,j] += (1-omega)**(np.abs(k-j))*geq[i,k]
            f_result[i] *= omega/(1-omega)

    for i in range(NP):
        cn = c[i,0]*n_hat[0] + c[i,1]*n_hat[1] + c[i,2]*n_hat[2]

        # Incoming populations
        if cn > 0:
            # Boundary condition at zero
            f_outgoing = f_result[p_bounceback[i],0]
            correction = correction_boundary(i, params_Kn, geq)
            f_result[i,0] = f_outgoing + correction

            for j in range(1,N):
                f_result[i,j] = 0
                for k in range(j):
                    f_result[i,j] += (1-omega)**(np.abs(k-j))*geq[i,k]
            f_result[i,1:] *= omega/(1-omega)
            for j in range(1,N):
                f_result[i,j] += (1-omega)**j*f_result[i,0]


    # Add mass correction
    correction = 0
    for i in range(NP):
        cn = c[i,0]*n_hat[0] + c[i,1]*n_hat[1] + c[i,2]*n_hat[2]

        # Outgoing populations
        if cn < 0:
            # Boundary condition at infinity
            correction += -omega*(f[i,0]-geq[i,0])
    f_result[0,0] += -1/omega*correction

    return f_result


# Stencil
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1]
cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1]
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1]
p_bounceback = np.array([0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15, 26, 25, 24, 23, 22, 21, 20, 19])
wp = np.array([8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0])

c = np.array([cx, cy, cz]).T
NP = len(wp)
cs2 = 1/3
N = 10

# Choices
omega = 1/1.2

n_hat = np.array([0.0, 1.0, 0.0])
rho_NS = np.full(N, 1)
u_NS = np.zeros(N)
v_NS = np.zeros(N)
w_NS = np.zeros(N)
vel_NS = np.empty((3, N), dtype=np.float64)
vel_NS[0] = u_NS
vel_NS[1] = v_NS
vel_NS[2] = w_NS


Force_x = np.full(N, 0)
Force_y = np.full(N, 0)
Force_y[0] = 1
Force_z = np.full(N, 0)
Force = np.array([Force_x, Force_y, Force_z])

# Initial condition
rho_Kn_0_ini = np.zeros(N)
u_Kn_0_ini = np.zeros(N)
v_Kn_0_ini = np.zeros(N)
w_Kn_0_ini = np.zeros(N)
f = np.full((NP, N), 0, dtype=np.float64)
for i in range(NP):
    f[i] = compute_eq(i, rho_Kn_0_ini, [u_Kn_0_ini, v_Kn_0_ini, w_Kn_0_ini])

n_iters = 0
while True:
    f_previous = f.copy()
    f = system_of_equations(f)
    n_iters += 1
    residual = np.sum((f_previous - f)**2)
    if (n_iters%100==0):
        print(n_iters, residual)
    if (residual < 1e-15):
        break

print(n_iters)

rho_Kn = np.sum(f, axis=0)
u_Kn = np.zeros(N, dtype=np.float64)
v_Kn = np.zeros(N, dtype=np.float64)
w_Kn = np.zeros(N, dtype=np.float64)
c_double = c.astype(np.float64)

for j in range(N):
    if rho_Kn[j] > 1e-15:
        for i in range(NP):
            u_Kn[j] += c_double[i,0] * f[i,j]
            v_Kn[j] += c_double[i,1] * f[i,j]
            w_Kn[j] += c_double[i,2] * f[i,j]
        u_Kn[j] /= rho_Kn[j]
        v_Kn[j] /= rho_Kn[j]
        w_Kn[j] /= rho_Kn[j]

plt.plot(rho_Kn)
plt.show()
plt.close()

plt.plot(u_Kn)
plt.show()

# plt.plot(rho_Kn*u_Kn)
# plt.show()

# plt.plot(rho_Kn*v_Kn)
# plt.show()

# plt.plot(rho_Kn*w_Kn)
# plt.show()

# eta = np.arange(N)
# plt.plot(eta, macro_Kn[0])
# plt.xlabel("$\\eta$")
# plt.ylabel("$\\rho^{\\text{Kn},(0)}$")
# plt.savefig("rho_N="+str(N)+".pdf")
# plt.close()

# plt.plot(eta, macro_Kn[1])
# plt.xlabel("$\\eta$")
# plt.ylabel("$u^{\\text{Kn},(0)}$")
# plt.savefig("u_N="+str(N)+".pdf")
# plt.close()

# plt.plot(eta, macro_Kn[2])
# plt.xlabel("$\\eta$")
# plt.ylabel("$v^{\\text{Kn},(0)}$")
# plt.savefig("v_N="+str(N)+".pdf")
# plt.close()

# plt.plot(eta, macro_Kn[3])
# plt.xlabel("$\\eta$")
# plt.ylabel("$w^{\\text{Kn},(0)}$")
# plt.savefig("w_N="+str(N)+".pdf")
# plt.close()