import numpy as np
from sympy import *

class DotZeroPrinter(StrPrinter):
    def _print_Integer(self, expr):
        return f"{expr}.0"

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
G =  Symbol("G")
G2 =  Symbol("G2")
Gx = Symbol("Gx_i")
G2x = Symbol("G2x_i")
Gy = Symbol("Gy_i")
G2y = Symbol("G2y_i")
Gz = Symbol("Gz_i")
G2z = Symbol("G2z_i")
sigma = Symbol("sigma")
gamma_sym = Symbol("gamma")

NP = 27

# De Rosis 2019, 10.1063/1.5124719
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1]  
cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1] 
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1]
wp = [8 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216]
B = [-16 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648]

freq = diag(*[1, 1, 1, 1, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega, omega])

indices_k = [0, 1, 2, 3]
# Bösch 2015, 10.1103/PhysRevE.92.043309, KBC-C4
indices_s = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] 
indices_h = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26]
for p in range(NP):
    if p in indices_h:
        freq[p,p] *= gamma_sym

T = zeros(NP, NP)
M = zeros(NP, NP)

for p in range(NP):
    # M Matrix
    M[p,0] = 1
    M[p,1] = cx[p]
    M[p,2] = cy[p]
    M[p,3] = cz[p]
    M[p,4] = cx[p]*cy[p]
    M[p,5] = cx[p]*cz[p]
    M[p,6] = cy[p]*cz[p]
    M[p,7] = cx[p]*cx[p] - cy[p]*cy[p]
    M[p,8] = cx[p]*cx[p] - cz[p]*cz[p]
    M[p,9] = cx[p]*cx[p] + cy[p]*cy[p] + cz[p]*cz[p]
    M[p,10] = cx[p]*cy[p]*cy[p]+cx[p]*cz[p]*cz[p]
    M[p,11] = cx[p]*cx[p]*cy[p]+cy[p]*cz[p]*cz[p]
    M[p,12] = cx[p]*cx[p]*cz[p]+cy[p]*cy[p]*cz[p]
    M[p,13] = cx[p]*cy[p]*cy[p]-cx[p]*cz[p]*cz[p]
    M[p,14] = cx[p]*cx[p]*cy[p]-cy[p]*cz[p]*cz[p]
    M[p,15] = cx[p]*cx[p]*cz[p]-cy[p]*cy[p]*cz[p]
    M[p,16] = cx[p]*cy[p]*cz[p]
    M[p,17] = cx[p]*cx[p]*cy[p]*cy[p]+cx[p]*cx[p]*cz[p]*cz[p]+cy[p]*cy[p]*cz[p]*cz[p]
    M[p,18] = cx[p]*cx[p]*cy[p]*cy[p]+cx[p]*cx[p]*cz[p]*cz[p]-cy[p]*cy[p]*cz[p]*cz[p]
    M[p,19] = cx[p]*cx[p]*cy[p]*cy[p]-cx[p]*cx[p]*cz[p]*cz[p]
    M[p,20] = cx[p]*cx[p]*cy[p]*cz[p]
    M[p,21] = cx[p]*cy[p]*cy[p]*cz[p]
    M[p,22] = cx[p]*cy[p]*cz[p]*cz[p]
    M[p,23] = cx[p]*cy[p]*cy[p]*cz[p]*cz[p]
    M[p,24] = cx[p]*cx[p]*cy[p]*cz[p]*cz[p]
    M[p,25] = cx[p]*cx[p]*cy[p]*cy[p]*cz[p]
    M[p,26] = cx[p]*cx[p]*cy[p]*cy[p]*cz[p]*cz[p]

    # T Matrix
    CX = cx[p] - U
    CY = cy[p] - V
    CZ = cz[p] - W

    CX2 = CX*CX
    CY2 = CY*CY
    CZ2 = CZ*CZ

    T[p,0] = 1
    T[p,1] = CX
    T[p,2] = CY
    T[p,3] = CZ
    T[p,4] = CX*CY
    T[p,5] = CX*CZ
    T[p,6] = CY*CZ
    T[p,7] = CX2 - CY2
    T[p,8] = CX2 - CZ2
    T[p,9] = CX2 + CY2 + CZ2
    T[p,10] = CX*CY2+CX*CZ2
    T[p,11] = CX2*CY+CY*CZ2
    T[p,12] = CX2*CZ+CY2*CZ
    T[p,13] = CX*CY2-CX*CZ2
    T[p,14] = CX2*CY-CY*CZ2
    T[p,15] = CX2*CZ-CY2*CZ
    T[p,16] = CX*CY*CZ
    T[p,17] = CX2*CY2+CX2*CZ2+CY2*CZ2
    T[p,18] = CX2*CY2+CX2*CZ2-CY2*CZ2
    T[p,19] = CX2*CY2-CX2*CZ2
    T[p,20] = CX2*CY*CZ
    T[p,21] = CX*CY2*CZ
    T[p,22] = CX*CY*CZ2
    T[p,23] = CX*CY2*CZ2
    T[p,24] = CX2*CY*CZ2
    T[p,25] = CX2*CY2*CZ
    T[p,26] = CX2*CY2*CZ2

printer = DotZeroPrinter()

T = simplify(T)
print("simplified T")

M = simplify(M)
print("simplified M")

keq = zeros(NP, 1)
keq[0,0] = rho
keq[9,0] = 3*rho*cs2
keq[17,0] = rho*cs2
keq[18,0] = rho*cs2**2
keq[26,0] = rho*cs2**4

eq = (T.T).inv()*keq
eq = eq.subs(U,0).subs(V,0).subs(W,0)
eq = simplify(eq)
for p in range(NP):
    print("eq[%d] = "%(p) + str(printer.doprint(eq[p,0])) + ";")

kf = zeros(NP, 1)
kf[1,0] = Fx
kf[2,0] = Fy
kf[3,0] = Fz
kf[10,0] = 2*Fx*cs2
kf[11,0] = 2*Fy*cs2
kf[12,0] = 2*Fz*cs2
kf[23,0] = Fx*cs2**2
kf[24,0] = Fy*cs2**2
kf[25,0] = Fz*cs2**2

pert = zeros(NP, 1)
for p in range(NP):
    Gc = Gx*cx[p] + Gy*cy[p] + Gz*cz[p]
    pert[p,0] = 9/8*sigma * G * (wp[p]*Gc*Gc/(G*G) - B[p])

k_pert = T.T*pert
k_pert = nsimplify(simplify(k_pert), tolerance=1e-12)
k_pert = k_pert.subs(G*G, G2)
k_pert = k_pert.subs(Gx*Gx, G2x)
k_pert = k_pert.subs(Gy*Gy, G2y)
k_pert = k_pert.subs(Gz*Gz, G2z)
k_pert = k_pert.subs(U*U, U2)
k_pert = k_pert.subs(V*V, V2)
k_pert = k_pert.subs(W*W, W2)
k_pert = nsimplify(simplify(k_pert), tolerance=1e-12)
for p in range(NP):
    print("k_pert[%d] = "%(p) + str(printer.doprint(apart(k_pert[p,0], G))) + ";")

k = Matrix(symbols('k[0:27]'))
k_pert = Matrix(symbols('k_pert[0:27]'))
k_star = (eye(NP) - freq)*k + freq*keq + (eye(NP) - 1/2*freq)*kf
k_star = k_star.subs(U*U, U2)
k_star = k_star.subs(V*V, V2)
k_star = k_star.subs(W*W, W2)
k_star = nsimplify(simplify(k_star), tolerance=1e-12)

for p in range(NP):
    print("k_star[%d] = "%(p) + str(printer.doprint(apart(k_star[p,0], G))) + ";")

k_star_hat = k_star.subs(gamma_sym, 0)
k_star_tilde = (k_star - k_star_hat)/gamma_sym
k_star_tilde = nsimplify(simplify(k_star_tilde), tolerance=1e-12)

for p in range(NP):
    print("k_star_hat[%d] = "%(p) + str(printer.doprint(apart(k_star_hat[p,0], G))) + ";")

for p in range(NP):
    print("k_star_tilde[%d] = "%(p) + str(printer.doprint(apart(k_star_tilde[p,0], G))) + ";")

k_star = Matrix(symbols('k_star[0:27]'))

N_inv = M.T * (T.T).inv()
raw = N_inv * k_star
print("computed raw")
raw = raw.subs(U*U, U2)
raw = raw.subs(V*V, V2)
raw = raw.subs(W*W, W2)
raw = nsimplify(simplify(raw), tolerance=1e-12)
for p in range(NP):
    print("raw[%d] = "%(p) + str(printer.doprint(collect(raw[p,0], k_star))) + ";")

raw_eq = N_inv * keq
raw_eq = raw_eq.subs(U*U, U2)
raw_eq = raw_eq.subs(V*V, V2)
raw_eq = raw_eq.subs(W*W, W2)
raw_eq = nsimplify(simplify(raw_eq), tolerance=1e-12)
for p in range(NP):
    print("raw[%d] = "%(p) + str(printer.doprint(collect(raw_eq[p,0], k_star))) + ";")

raw = Matrix(symbols("raw[0:27]"))
f = (M.T).inv()*raw
f = nsimplify(simplify(f), tolerance=1e-12)
for p in range(NP):
    print("f2[INDEX_F(i, j, k, %d, n)] = "%(p) + str(printer.doprint(f[p,0])) + ";")