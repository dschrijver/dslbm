import numpy as np
from sympy import *

class DotZeroPrinter(StrPrinter):
    def _print_Integer(self, expr):
        return f"{expr}.0"

def m(a, b, c):
    return cx[i]**a * cy[i]**b * cz[i]**c

def t(a, b, c):
    return CX**a * CY**b * CZ**c

U = Symbol("u_i")
U2 = Symbol("u2_i")
U3 = Symbol("u3_i")
U4 = Symbol("u4_i")
V = Symbol("v_i")
V2 = Symbol("v2_i")
V3 = Symbol("v3_i")
V4 = Symbol("v4_i")
W = Symbol("w_i")
W2 = Symbol("w2_i")
W3 = Symbol("w3_i")
W4 = Symbol("w4_i")
rho = Symbol("rho_i")
Fx = Symbol("Fx_i")
Fy = Symbol("Fy_i")
Fz = Symbol("Fz_i")
omega = Symbol("omega")
nu = Symbol("nu")
cs2 = Symbol("cs2_i")
cs4 = Symbol("cs4_i")
cs6 = Symbol("cs6_i")
cs8 = Symbol("cs8_i")
G =  Symbol("G")
G2 =  Symbol("G2")
Gx = Symbol("Gx_i")
G2x = Symbol("G2x_i")
Gy = Symbol("Gy_i")
G2y = Symbol("G2y_i")
Gz = Symbol("Gz_i")
G2z = Symbol("G2z_i")
sigma = Symbol("sigma")
Qx = Symbol("Qx_i")
Qy = Symbol("Qy_i")
Qz = Symbol("Qz_i")
px = Symbol("px")
py = Symbol("py")
pz = Symbol("pz")
pressure = Symbol("pressure")

NP = 27

# Saito 2023, 10.1103/PhysRevE.108.065305
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1]
cy = [0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1]
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1]
p_bounceback = [0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25]
wp = [8 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216]
B = [-16 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648]

freq = diag(*[1, 1, 1, 1, omega, omega, omega, omega, omega, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

T = zeros(NP, NP)
M = zeros(NP, NP)
for i in range(NP):
    # Saito 2023, 10.1103/PhysRevE.108.065305
    M[0,i] = 1

    M[1,i] = m(1,0,0)
    M[2,i] = m(0,1,0)
    M[3,i] = m(0,0,1)

    M[4,i] = m(1,1,0)
    M[5,i] = m(1,0,1)
    M[6,i] = m(0,1,1)
    M[7,i] = m(2,0,0) - m(0,2,0)
    M[8,i] = m(2,0,0) - m(0,0,2)
    M[9,i] = m(2,0,0) + m(0,2,0) + m(0,0,2)

    M[10,i] = m(1,2,0)
    M[11,i] = m(1,0,2)
    M[12,i] = m(0,1,2)
    M[13,i] = m(2,1,0)
    M[14,i] = m(2,0,1)
    M[15,i] = m(0,2,1)
    M[16,i] = m(1,1,1)

    M[17,i] = m(2,2,0)
    M[18,i] = m(2,0,2)
    M[19,i] = m(0,2,2)
    M[20,i] = m(2,1,1)
    M[21,i] = m(1,2,1)
    M[22,i] = m(1,1,2)

    M[23,i] = m(1,2,2)
    M[24,i] = m(2,1,2)
    M[25,i] = m(2,2,1)

    M[26,i] = m(2,2,2) 

    # Saito 2023, 10.1103/PhysRevE.108.065305
    CX = cx[i] - U
    CY = cy[i] - V
    CZ = cz[i] - W

    T[0,i] = 1

    T[1,i] = t(1,0,0)
    T[2,i] = t(0,1,0)
    T[3,i] = t(0,0,1)

    T[4,i] = t(1,1,0)
    T[5,i] = t(1,0,1)
    T[6,i] = t(0,1,1)
    T[7,i] = t(2,0,0) - t(0,2,0)
    T[8,i] = t(2,0,0) - t(0,0,2)
    T[9,i] = t(2,0,0) + t(0,2,0) + t(0,0,2)

    T[10,i] = t(1,2,0)
    T[11,i] = t(1,0,2)
    T[12,i] = t(0,1,2)
    T[13,i] = t(2,1,0)
    T[14,i] = t(2,0,1)
    T[15,i] = t(0,2,1)
    T[16,i] = t(1,1,1)

    T[17,i] = t(2,2,0)
    T[18,i] = t(2,0,2)
    T[19,i] = t(0,2,2)
    T[20,i] = t(2,1,1)
    T[21,i] = t(1,2,1)
    T[22,i] = t(1,1,2)

    T[23,i] = t(1,2,2)
    T[24,i] = t(2,1,2)
    T[25,i] = t(2,2,1)

    T[26,i] = t(2,2,2) 

printer = DotZeroPrinter()

T = simplify(T)
print("simplified T")

M = simplify(M)
print("simplified M")

N_shift = simplify(T*M.inv())
print("simplified N_shift")

Omega = N_shift.inv() * freq * N_shift
Omega = nsimplify(simplify(Omega), tolerance=1e-12)
Omega = Omega.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
Omega = nsimplify(simplify(Omega), tolerance=1e-12)
for i in range(NP):
    print(Omega[i,:])

# Saito 2023, 10.1103/PhysRevE.108.065305
teq = zeros(NP, 1)
teq[0,0] = rho
teq[9,0] = 3*pressure
teq[17,0] = 1/3*pressure
teq[18,0] = 1/3*pressure
teq[19,0] = 1/3*pressure
teq[26,0] = 1/9*pressure

meq = N_shift.inv() * teq
meq = nsimplify(simplify(meq), tolerance=1e-12)
meq = meq.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
meq = nsimplify(simplify(meq), tolerance=1e-12)
for i in range(NP):
    print("meq[%d] = "%(i) + str(printer.doprint(meq[i,0])) + ";")

eq = T.inv() * teq
eq = nsimplify(simplify(eq), tolerance=1e-12)
eq = eq.subs([(V, 0), (W, 0), (cs2, 9/19*(1-0.9992))])
# eq = eq.subs([(U, 0), (V, 0), (W, 0), (U2, 0), (V2, 0), (W2, 0)])
eq = nsimplify(simplify(eq), tolerance=1e-12)
for i in range(NP):
    sols = solve(eq[i,0], U)
    real_sols = [s for s in sols if s.is_real]
    print(real_sols)
    print("eq[%d] = "%(i) + str(printer.doprint(eq[i,0])) + ";")

# Saito 2023, 10.1103/PhysRevE.108.065305
ts = zeros(NP, 1)
ts[1,0] = Fx
ts[2,0] = Fy
ts[3,0] = Fz
ts[10,0] = Fx/3
ts[11,0] = Fx/3
ts[12,0] = Fy/3
ts[13,0] = Fy/3
ts[14,0] = Fz/3
ts[15,0] = Fz/3
ts[23,0] = Fx/9
ts[24,0] = Fy/9
ts[25,0] = Fz/9

# Saito 2023, 10.1103/PhysRevE.108.065305
tc = zeros(NP, 1)
tc[7,0] = Qx - Qy
tc[8,0] = Qx - Qz
tc[9,0] = Qx + Qy + Qz

f = Matrix(symbols("f((0:27))"))
mf = M*f
mf = nsimplify(simplify(mf), tolerance=1e-12)
for i in range(NP):
    print("mf[%d] = "%(i) + str(printer.doprint(mf[i,0])) + ";")

mf = Matrix(symbols("mf[0:27]"))
tf = N_shift*mf
tf = nsimplify(simplify(tf), tolerance=1e-12)
tf = tf.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
tf = nsimplify(simplify(tf), tolerance=1e-12)
for i in range(NP):
    print("tf[%d] = "%(i) + str(printer.doprint(tf[i,0])) + ";")

feq = T.inv() * teq
feq_stationary = feq.subs([(U, 0), (V, 0), (W, 0)])
feq_stationary = nsimplify(simplify(feq_stationary), tolerance=1e-12)
rec = zeros(NP, 1)
for i in range(NP):
    rec[i,0] = Symbol("A")*(Gx*cx[i] + Gy*cy[i] + Gz*cz[i])/G*feq_stationary[i,0]

mrec = zeros(NP, 1)
mrec[1] = px
mrec[2] = py
mrec[3] = pz

trec = N_shift * mrec
trec = nsimplify(simplify(trec), tolerance=1e-12)
trec = trec.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
trec = nsimplify(simplify(trec), tolerance=1e-12)
for i in range(NP):
    print("trec[%d] = "%(i) + str(printer.doprint(trec[i,0])) + ";")

trec = zeros(NP, 1)
trec[1] = px
trec[2] = py
trec[3] = pz

mrec = N_shift.inv() * trec
mrec = nsimplify(simplify(mrec), tolerance=1e-12)
mrec = mrec.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
mrec = nsimplify(simplify(mrec), tolerance=1e-12)
for i in range(NP):
    print("mrec[%d] = "%(i) + str(printer.doprint(mrec[i,0])) + ";")

mrec = Matrix(symbols("mrec[0:27]"))
rec = M.inv() * mrec
rec = nsimplify(simplify(rec), tolerance=1e-12)
for i in range(NP):
    print("rec[%d] = "%(i) + str(printer.doprint(rec[i,0])) + ";")

# krec = zeros(NP, 1)
# krec[1,0] = Symbol("px")
# krec[2,0] = Symbol("py")
# krec[3,0] = Symbol("pz")
# rec = T.inv() * krec
# rec = nsimplify(simplify(rec), tolerance=1e-12)
# rec = rec.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
# rec = nsimplify(simplify(rec), tolerance=1e-12)
# for p in range(NP):
#     print("rec[%d] = "%(p) + str(printer.doprint(rec[p,0])) + ";")

# kf = T*Force
# kf = nsimplify(simplify(kf), tolerance=1e-12)
# kf = kf.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
# kf = nsimplify(simplify(kf), tolerance=1e-12)

# pert = zeros(NP, 1)
# for p in range(NP):
#     Gc = Gx*cx[p] + Gy*cy[p] + Gz*cz[p]
#     pert[p,0] = 9/8*sigma * G * (wp[p]*Gc*Gc/(G*G) - B[p])

# C = zeros(NP, 1)
# C[7,0] = Symbol("Qx") - Symbol("Qy")
# C[8,0] = Symbol("Qx") - Symbol("Qz")
# C[9,0] = Symbol("Qx") + Symbol("Qy") + Symbol("Qz")

# k_pert = T*pert
# k_pert = nsimplify(simplify(k_pert), tolerance=1e-12)
# k_pert = k_pert.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8), (G**2, G2), (Gx**2, G2x), (Gy**2, G2y), (Gz**2, G2z)])
# k_pert = nsimplify(simplify(k_pert), tolerance=1e-12)
# for p in range(NP):
#     print("k_pert[%d] = "%(p) + str(printer.doprint(apart(k_pert[p,0], G))) + ";")

# k = Matrix(symbols('k[0:27]'))
# k_pert = Matrix(symbols('k_pert[0:27]'))
# k_star = (eye(NP) - freq)*k + freq*keq + (eye(NP) - 1/2*freq)*kf + (eye(NP) - 1/2*freq)*C + freq*k_pert
# k_star = k_star.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
# k_star = nsimplify(simplify(k_star), tolerance=1e-12)

# for p in range(NP):
#     print("k_star[%d] = "%(p) + str(printer.doprint(apart(k_star[p,0], G))) + ";")

# k_star = Matrix(symbols('k_star[0:27]'))

# N_inv = M * T.inv()
# raw = N_inv * k_star
# print("computed raw")
# raw = raw.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
# raw = nsimplify(simplify(raw), tolerance=1e-12)
# for p in range(NP):
#     print("raw[%d] = "%(p) + str(printer.doprint(collect(raw[p,0], k_star))) + ";")

# raw_eq = N_inv * keq
# raw_eq = raw_eq.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
# raw_eq = nsimplify(simplify(raw_eq), tolerance=1e-12)
# for p in range(NP):
#     print("raw[%d] = "%(p) + str(printer.doprint(collect(raw_eq[p,0], k_star))) + ";")

# raw = Matrix(symbols("raw[0:27]"))
# f = M.inv()*raw
# f = nsimplify(simplify(f), tolerance=1e-12)
# for p in range(NP):
#     print("f2[INDEX_F(i, j, k, %d, n)] = "%(p) + str(printer.doprint(f[p,0])) + ";")