import numpy as np
from sympy import *

class DotZeroPrinter(StrPrinter):
    def _print_Integer(self, expr):
        return f"{expr}.0"

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
alpha = 1-19/9*cs2

NP = 27

# De Rosis 2019, 10.1063/1.5124719
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1]  
cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1] 
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1]
wp = [8 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216]
B = [-16 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648]
phi = [alpha, 2/19*(1-alpha), 2/19*(1-alpha), 2/19*(1-alpha), 2/19*(1-alpha), 2/19*(1-alpha), 2/19*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/38*(1-alpha), 1/152*(1-alpha), 1/152*(1-alpha), 1/152*(1-alpha), 1/152*(1-alpha), 1/152*(1-alpha), 1/152*(1-alpha), 1/152*(1-alpha), 1/152*(1-alpha)]

freq = diag(*[1, 1, 1, 1, omega, omega, omega, omega, omega, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

T = zeros(NP, NP)
T_natural = zeros(NP, NP)
M = zeros(NP, NP)
feq_wen = zeros(NP, 1)

for p in range(NP):
    # Wen 2019, 10.1103/PhysRevE.100.023301
    uc = U*cx[p] + V*cy[p] + W*cz[p]
    u2 = U*U + V*V + W*W
    c2 = cx[p]*cx[p] + cy[p]*cy[p] + cz[p]*cz[p]

    first_order = (U*cx[p]+V*cy[p]+W*cz[p])/(1/3)
    second_order = (U*cx[p]+V*cy[p]+W*cz[p])**2/(2*(1/9)) - (U*U+V*V+W*W)/(2*(1/3))
    third_order = ((cx[p]**2-(1/3))*cy[p]*U*U*V + (cx[p]**2-(1/3))*cz[p]*U*U*W + 
                   (cy[p]**2-(1/3))*cx[p]*U*V*V + (cz[p]**2-(1/3))*cx[p]*U*W*W + 
                   (cz[p]**2-(1/3))*cy[p]*V*W*W + (cy[p]**2-(1/3))*cz[p]*V*V*W + 
                   2*( cx[p]*cy[p]*cz[p]*U*V*W)) / (2*(1/27))
    fourth_order = ((cx[p]**2-(1/3))*(cy[p]**2-(1/3))*U*U*V*V + 
                    (cx[p]**2-(1/3))*(cz[p]**2-(1/3))*U*U*W*W + 
                    (cy[p]**2-(1/3))*(cz[p]**2-(1/3))*V*V*W*W + 
                  2*(cx[p]*cy[p]*(cz[p]**2-(1/3))*U*V*W*W + 
                     cx[p]*(cy[p]**2-(1/3))*cz[p]*U*V*V*W + 
                     (cx[p]**2-(1/3))*cy[p]*cz[p]*U*U*V*W))/(4*(1/81))
    fifth_order = ((cx[p]**2-(1/3))*cy[p]*(cz[p]**2-(1/3))*U*U*V*W*W + 
                   (cx[p]**2-(1/3))*(cy[p]**2-(1/3))*cz[p]*U*U*V*V*W + 
                   cx[p]*(cy[p]**2-(1/3))*(cz[p]**2-(1/3))*U*V*V*W*W) / (4*(1/243))
    sixth_order = ((cx[p]**2-(1/3))*(cy[p]**2-(1/3))*(cz[p]**2-(1/3))*U*U*V*V*W*W) / (8*(1/729))

    feq_wen[p,0] = rho*(phi[p] + wp[p]*(first_order + second_order + third_order + fourth_order + fifth_order + sixth_order + 3/2*uc*(3*cs2 - 1)*(3*c2 - 5)))

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

    T_natural[p,0] = 1

    T_natural[p,1] = t(1,0,0)
    T_natural[p,2] = t(0,1,0)
    T_natural[p,3] = t(0,0,1)

    T_natural[p,4] = t(1,1,0)
    T_natural[p,5] = t(0,1,1)
    T_natural[p,6] = t(1,0,1)
    T_natural[p,7] = t(2,0,0)
    T_natural[p,8] = t(0,2,0)
    T_natural[p,9] = t(0,0,2)

    T_natural[p,10] = t(1,2,0)
    T_natural[p,11] = t(1,0,2)
    T_natural[p,12] = t(0,1,2)
    T_natural[p,13] = t(2,1,0)
    T_natural[p,14] = t(2,0,1)
    T_natural[p,15] = t(0,2,1)
    T_natural[p,16] = t(1,1,1)

    T_natural[p,17] = t(2,2,0)
    T_natural[p,18] = t(2,0,2)
    T_natural[p,19] = t(0,2,2)
    T_natural[p,20] = t(2,1,1)
    T_natural[p,21] = t(1,2,1)
    T_natural[p,22] = t(1,1,2)

    T_natural[p,23] = t(1,2,2)
    T_natural[p,24] = t(2,1,2)
    T_natural[p,25] = t(2,2,1)

    T_natural[p,26] = t(2,2,2) 


printer = DotZeroPrinter()

T = simplify(T)
print("simplified T")

M = simplify(M)
print("simplified M")

fneq = Matrix(symbols('fneq[0:27]'))
mneq = M.T*fneq
mneq = nsimplify(simplify(mneq), tolerance=1e-12)
mhneq = zeros(NP, 1)
for p in range(NP):
    if p in [4, 5, 6, 7, 8, 9]:
        mhneq[p,0] = mneq[p,0]
for p in range(NP):
    print("mhneq[%d] = "%(p) + str(printer.doprint(mhneq[p,0])) + ";")
sneq = (M.T).inv()*mhneq
sneq = nsimplify(simplify(sneq), tolerance=1e-12)
for p in range(NP):
    print("sneq[%d] = "%(p) + str(printer.doprint(sneq[p,0])) + ";")
for p in range(NP):
    if mhneq[p] != 0:
        mhneq[p] = Symbol("mhneq[%d]"%p)
sneq = (M.T).inv()*mhneq
sneq = nsimplify(simplify(sneq), tolerance=1e-12)
for p in range(NP):
    print("sneq[%d] = "%(p) + str(printer.doprint(sneq[p,0])) + ";")

# My equilibrium
keq = zeros(NP, 1)
keq[0,0] = rho
keq[9,0] = 3*rho*cs2
keq[17,0] = rho*cs2
keq[18,0] = rho*cs2**2
keq[26,0] = rho*cs2**3

# Wen 2019, 10.1103/PhysRevE.100.023301
feq_wen = nsimplify(simplify(feq_wen), tolerance=1e-12)
keq_wen = T_natural.T*feq_wen
keq_wen = nsimplify(simplify(keq_wen), tolerance=1e-12)
keq_wen = keq_wen.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
for p in range(NP):
    print("keq_wen[%d] = "%(p) + str(printer.doprint(keq_wen[p,0])) + ";")

eq = (T.T).inv()*keq
eq = eq.subs([(U, 0), (V, 0), (W, 0), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
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
k_pert = k_pert.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8), (G**2, G2), (Gx**2, G2x), (Gy**2, G2y), (Gz**2, G2z)])
k_pert = nsimplify(simplify(k_pert), tolerance=1e-12)
for p in range(NP):
    print("k_pert[%d] = "%(p) + str(printer.doprint(apart(k_pert[p,0], G))) + ";")

k = Matrix(symbols('k[0:27]'))
k_pert = Matrix(symbols('k_pert[0:27]'))
k_star = (eye(NP) - freq)*k + freq*keq + (eye(NP) - 1/2*freq)*kf + freq*k_pert
k_star = k_star.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
k_star = nsimplify(simplify(k_star), tolerance=1e-12)

for p in range(NP):
    print("k_star[%d] = "%(p) + str(printer.doprint(apart(k_star[p,0], G))) + ";")

k_star = Matrix(symbols('k_star[0:27]'))

N_inv = M.T * (T.T).inv()
raw = N_inv * k_star
print("computed raw")
raw = raw.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
raw = nsimplify(simplify(raw), tolerance=1e-12)
for p in range(NP):
    print("raw[%d] = "%(p) + str(printer.doprint(collect(raw[p,0], k_star))) + ";")

raw_eq = N_inv * keq
raw_eq = raw_eq.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
raw_eq = nsimplify(simplify(raw_eq), tolerance=1e-12)
for p in range(NP):
    print("raw[%d] = "%(p) + str(printer.doprint(collect(raw_eq[p,0], k_star))) + ";")

raw = Matrix(symbols("raw[0:27]"))
f = (M.T).inv()*raw
f = nsimplify(simplify(f), tolerance=1e-12)
for p in range(NP):
    print("f2[INDEX_F(i, j, k, %d, n)] = "%(p) + str(printer.doprint(f[p,0])) + ";")