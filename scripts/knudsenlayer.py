import numpy as np
from scipy.optimize import fsolve
from numba import njit
import matplotlib.pyplot as plt

@njit
def compute_eq(i, rho, vel):
    uc = np.zeros(N, dtype=np.float64)
    u2 = np.zeros(N, dtype=np.float64)
    for alpha in range(3):
        uc += vel[alpha]*c[i,alpha]
        u2 += vel[alpha]*vel[alpha]
    return wp[i] * rho * (1 + uc/cs2 + uc*uc/(2*cs2*cs2) - u2/(2*cs2))

@njit
def compute_geq(params_Kn, params_NS):
    rho_Kn = params_Kn[0]
    u_Kn = params_Kn[1]
    v_Kn = params_Kn[2]
    w_Kn = params_Kn[3]

    rho_NS = params_NS[0]
    u_NS = params_NS[1]
    v_NS = params_NS[2]
    w_NS = params_NS[3]

    vel_NS = np.empty((3, N), dtype=np.float64)
    vel_NS[0] = u_NS
    vel_NS[1] = v_NS
    vel_NS[2] = w_NS
    
    rho = rho_Kn + rho_NS
    u = 1/rho*(rho_Kn*u_Kn + rho_NS*u_NS)
    v = 1/rho*(rho_Kn*v_Kn + rho_NS*v_NS)
    w = 1/rho*(rho_Kn*w_Kn + rho_NS*w_NS)

    vel = np.empty((3, N), dtype=np.float64)
    vel[0] = u
    vel[1] = v
    vel[2] = w

    result = np.zeros((NP, N), dtype=np.float64)
    for i in range(NP):
        result[i] = compute_eq(i, rho, vel) - compute_eq(i, rho_NS, vel_NS)

    return result

@njit
def outgoing_0(i, geq):
    result = np.zeros(N, dtype=np.float64)

    geq_i = geq[i]
    for j in range(N):
        for k in range(j+1,N):
            result[j] += (1 - omega)**np.abs(k-j) * geq_i[k]
    return result * omega/(1-omega)

@njit
def outgoing_0_wall(i, geq):
    result = 0.0

    geq_i = geq[i]
    for k in range(1,N):
        result += (1 - omega)**k * geq_i[k]
    return result * omega/(1-omega)

@njit
def correction_boundary_0(i, params_Kn, params_NS, geq):
    rho_Kn = params_Kn[0]
    u_Kn = params_Kn[1]
    v_Kn = params_Kn[2]
    w_Kn = params_Kn[3]

    rho_NS = params_NS[0]
    u_NS = params_NS[1]
    v_NS = params_NS[2]
    w_NS = params_NS[3]
    
    rho = rho_Kn + rho_NS
    u = 1/rho*(rho_Kn*u_Kn + rho_NS*u_NS)
    v = 1/rho*(rho_Kn*v_Kn + rho_NS*v_NS)
    w = 1/rho*(rho_Kn*w_Kn + rho_NS*w_NS)

    rho = rho[0]

    vel = np.empty(3, dtype=np.float64)
    vel[0] = u[0]
    vel[1] = v[0]
    vel[2] = w[0]

    uc = np.dot(vel, c[i,:].astype(np.float64))

    N = np.zeros(3, dtype=np.float64)
    for alpha in range(3):
        if n_hat[alpha] > 0:
            N[alpha] = 1/cs2 * Force_0[alpha,0]
        else:
            A = 1/3
            sum_parallel = 0
            for i2 in range(NP):
                cn = np.dot(c[i2,:].astype(np.float64), n_hat)
                if cn == 0:
                    sum_parallel += c[i2,alpha]*geq[i2,0]
            N[alpha] = 2/(cs2*A) * (sum_parallel + (A - 1)*rho*vel[alpha] + 1/2*Force_0[alpha,0])

    return 2*wp[i]*rho*uc/cs2 - wp[i] * np.dot(c[i,:].astype(np.float64), N)

@njit
def incoming_0(i, params_Kn, params_NS, geq):
    result = np.zeros(N, dtype=np.float64)

    geq_i = geq[i]
    for j in range(N):
        for k in range(j):
            result[j] += (1 - omega)**np.abs(k-j) * geq_i[k]
        result[j] *= omega/(1-omega)

        # Boundary condition
        f_outgoing = outgoing_0_wall(p_bounceback[i], geq)
        correction_boundary = correction_boundary_0(i, params_Kn, params_NS, geq)

        result[j] += (1 - omega)**j * (f_outgoing + correction_boundary)
    return result

@njit
def system_of_equations_predictor(params_Kn):
    params_NS = macro_NS_0

    f = np.zeros((NP,N), dtype=np.float64)
    geq = compute_geq(params_Kn, params_NS)

    # zero-velocity population
    f[0] = geq[0]

    # mass conservation
    # f[0,0] = 1/2*Force_y[0]

    for i in range(1, NP):
        cn = np.dot(c[i,:].astype(np.float64), n_hat)

        # parallel populations
        if cn == 0:
            f[i] = geq[i]

        # outgoing populations
        elif cn < 0:
            f[i] = outgoing_0(i, geq)

        # incoming populations
        else:
            f[i] = incoming_0(i, params_Kn, params_NS, geq)

    zeroth_moment = np.zeros(N, dtype=np.float64)
    for i in range(NP):
        zeroth_moment += f[i]

    first_moments = np.zeros((3, N), dtype=np.float64)
    for i in range(NP):
        for alpha in range(3):
            first_moments[alpha] += c[i,alpha]*f[i]

    result = np.zeros((4, N), dtype=np.float64)
    result[0] = zeroth_moment
    result[1] = first_moments[0]
    result[2] = first_moments[1]
    result[3] = first_moments[2]

    return result

@njit
def system_of_equations(params_Kn):
    params_Kn = params_Kn.reshape((4,N))

    params_NS = macro_NS_0

    f = np.zeros((NP,N), dtype=np.float64)
    geq = compute_geq(params_Kn, params_NS)

    # zero-velocity population
    f[0] = geq[0]

    # mass conservation
    # f[0,0] = 1/2*Force_y[0]

    for i in range(1, NP):
        cn = np.dot(c[i,:].astype(np.float64), n_hat)

        # parallel populations
        if cn == 0:
            f[i] = geq[i]

        # outgoing populations
        elif cn < 0:
            f[i] = outgoing_0(i, geq)

        # incoming populations
        else:
            f[i] = incoming_0(i, params_Kn, params_NS, geq)

    zeroth_moment = np.zeros(N, dtype=np.float64)
    for i in range(NP):
        zeroth_moment += f[i]

    first_moments = np.zeros((3, N), dtype=np.float64)
    for i in range(NP):
        for alpha in range(3):
            first_moments[alpha] += c[i,alpha]*f[i]

    result = np.zeros((4, N), dtype=np.float64)
    result[0] = params_Kn[0] - zeroth_moment
    result[1] = params_Kn[1] - first_moments[0]
    result[2] = params_Kn[2] - first_moments[1]
    result[3] = params_Kn[3] - first_moments[2]

    return result.flatten()


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
rho_NS_0 = np.full(N, 1)
u_NS_0 = np.zeros(N)
v_NS_0 = np.zeros(N)
w_NS_0 = np.zeros(N)
macro_NS_0 = np.array([rho_NS_0, u_NS_0, v_NS_0, w_NS_0])

Force_x = np.full(N, 0)
Force_y = np.full(N, 0.001)
Force_z = np.full(N, 0)
Force_0 = np.array([Force_x, Force_y, Force_z])

# Initial condition
rho_Kn_0_ini = np.zeros(N)
rho_Kn_0_ini[0] = 0.1
u_Kn_0_ini = np.zeros(N)
v_Kn_0_ini = np.zeros(N)
w_Kn_0_ini = np.zeros(N)
macro_Kn_ini = np.array([rho_Kn_0_ini, u_Kn_0_ini, v_Kn_0_ini, w_Kn_0_ini])

solution = fsolve(system_of_equations, macro_Kn_ini.flatten()).reshape(4,N)
macro_Kn = solution

# macro_Kn = macro_Kn_ini
# macro_Kn_previous = macro_Kn

# l = 0.1
# n_iters = 0
# while n_iters < 1000000:
#     macro_Kn = l*system_of_equations_predictor(macro_Kn) + (l-1)*macro_Kn
#     residual = np.sum((macro_Kn - macro_Kn_previous)**2)
#     macro_Kn_previous = macro_Kn
#     n_iters += 1
#     print(n_iters, residual)
#     if (residual < 1e-10):
#         break


eta = np.arange(N)
plt.plot(eta, macro_Kn[0])
plt.xlabel("$\\eta$")
plt.ylabel("$\\rho^{\\text{Kn},(0)}$")
plt.savefig("rho_N="+str(N)+".pdf")
plt.close()

plt.plot(eta, macro_Kn[1])
plt.xlabel("$\\eta$")
plt.ylabel("$u^{\\text{Kn},(0)}$")
plt.savefig("u_N="+str(N)+".pdf")
plt.close()

plt.plot(eta, macro_Kn[2])
plt.xlabel("$\\eta$")
plt.ylabel("$v^{\\text{Kn},(0)}$")
plt.savefig("v_N="+str(N)+".pdf")
plt.close()

plt.plot(eta, macro_Kn[3])
plt.xlabel("$\\eta$")
plt.ylabel("$w^{\\text{Kn},(0)}$")
plt.savefig("w_N="+str(N)+".pdf")
plt.close()