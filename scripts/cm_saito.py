import numpy as np
from sympy import *

class DotZeroPrinter(StrPrinter):
    def _print_Integer(self, expr):
        return f"{expr}.0"

def m(a, b, c):
    return cx[p]**a * cy[p]**b * cz[p]**c

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
pressure = rho*cs2

NP = 27

# De Rosis 2019, 10.1063/1.5124719
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1]  
cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1] 
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1]
wp = [8 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216]
B = [-16 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 2 / 81, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648, 7 / 648]

freq = diag(*[1, 1, 1, 1, omega, omega, omega, omega, omega, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

T = zeros(NP, NP)
M = zeros(NP, NP)
Force = zeros(NP, 1)

for p in range(NP):
    # Saito 2023, 10.1103/PhysRevE.108.065305
    M[0,p] = 1

    M[1,p] = m(1,0,0)
    M[2,p] = m(0,1,0)
    M[3,p] = m(0,0,1)

    M[4,p] = m(1,1,0)
    M[5,p] = m(1,0,1)
    M[6,p] = m(0,1,1)
    M[7,p] = m(2,0,0) - m(0,2,0)
    M[8,p] = m(2,0,0) - m(0,0,2)
    M[9,p] = m(2,0,0) + m(0,2,0) + m(0,0,2)

    M[10,p] = m(1,2,0)
    M[11,p] = m(1,0,2)
    M[12,p] = m(0,1,2)
    M[13,p] = m(2,1,0)
    M[14,p] = m(2,0,1)
    M[15,p] = m(0,2,1)
    M[16,p] = m(1,1,1)

    M[17,p] = m(2,2,0)
    M[18,p] = m(2,0,2)
    M[19,p] = m(0,2,2)
    M[20,p] = m(2,1,1)
    M[21,p] = m(1,2,1)
    M[22,p] = m(1,1,2)

    M[23,p] = m(1,2,2)
    M[24,p] = m(2,1,2)
    M[25,p] = m(2,2,1)

    M[26,p] = m(2,2,2) 

    # Saito 2023, 10.1103/PhysRevE.108.065305
    CX = cx[p] - U
    CY = cy[p] - V
    CZ = cz[p] - W

    T[0,p] = 1

    T[1,p] = t(1,0,0)
    T[2,p] = t(0,1,0)
    T[3,p] = t(0,0,1)

    T[4,p] = t(1,1,0)
    T[5,p] = t(1,0,1)
    T[6,p] = t(0,1,1)
    T[7,p] = t(2,0,0) - t(0,2,0)
    T[8,p] = t(2,0,0) - t(0,0,2)
    T[9,p] = t(2,0,0) + t(0,2,0) + t(0,0,2)

    T[10,p] = t(1,2,0)
    T[11,p] = t(1,0,2)
    T[12,p] = t(0,1,2)
    T[13,p] = t(2,1,0)
    T[14,p] = t(2,0,1)
    T[15,p] = t(0,2,1)
    T[16,p] = t(1,1,1)

    T[17,p] = t(2,2,0)
    T[18,p] = t(2,0,2)
    T[19,p] = t(0,2,2)
    T[20,p] = t(2,1,1)
    T[21,p] = t(1,2,1)
    T[22,p] = t(1,1,2)

    T[23,p] = t(1,2,2)
    T[24,p] = t(2,1,2)
    T[25,p] = t(2,2,1)

    T[26,p] = t(2,2,2) 

    # Forcing
    hat_cx = cx[p]/sqrt(1/3)
    hat_cy = cy[p]/sqrt(1/3)
    hat_cz = cz[p]/sqrt(1/3)
    H1_1 = hat_cx
    H1_2 = hat_cy
    H1_3 = hat_cz
    H2_1 = hat_cx**2-1
    H2_2 = hat_cx*hat_cy
    H2_3 = hat_cx*hat_cz
    H2_4 = hat_cy*hat_cx
    H2_5 = hat_cy**2-1
    H2_6 = hat_cy*hat_cz
    H2_7 = hat_cz*hat_cx
    H2_8 = hat_cz*hat_cy
    H2_9 = hat_cz**2-1  
    H3_1 = (hat_cx**2-1)*hat_cy
    H3_2 = (hat_cx**2-1)*hat_cz
    H3_3 = hat_cx*(hat_cy**2-1)
    H3_4 = hat_cx*(hat_cz**2-1)
    H3_5 = hat_cy*(hat_cz**2-1)
    H3_6 = (hat_cy**2-1)*hat_cz
    H3_7 = hat_cx*hat_cy*hat_cz
    H3_8 = H3_7
    H3_9 = H3_7
    H4_1 = (hat_cx**2-1)*(hat_cy**2-1)
    H4_2 = (hat_cx**2-1)*(hat_cz**2-1)
    H4_3 = (hat_cy**2-1)*(hat_cz**2-1)
    H4_4 = hat_cx*hat_cy*(hat_cz**2-1)
    H4_5 = hat_cx*(hat_cy**2-1)*hat_cz
    H4_6 = (hat_cx**2-1)*hat_cy*hat_cz
    H5_1 = (hat_cx**2-1)*hat_cy*(hat_cz**2-1)
    H5_2 = (hat_cx**2-1)*(hat_cy**2-1)*hat_cz
    H5_3 = hat_cx*(hat_cy**2-1)*(hat_cz**2-1)
    H6_1 = (hat_cx**2-1)*(hat_cy**2-1)*(hat_cz**2-1)

    first_order = (Fx*H1_1+Fy*H1_2+Fz*H1_3)/sqrt(1/3)
    second_order = 1/(2*1/3)*(H2_1*(Fx*U+U*Fx)+
                              H2_5*(Fy*V+V*Fy)+ 
                              H2_9*(Fz*W+W*Fz)+                                    
                            2*H2_2*(Fx*V+U*Fy)+
                            2*H2_3*(Fx*W+U*Fz)+
                            2*H2_6*(Fy*W+V*Fz) )
    third_order =  1/(6*sqrt(1/3)**3)*(3*H3_1*(Fx*U*V+U*Fx*V+U*U*Fy)+
                              3*H3_2*(Fx*U*W+U*Fx*W+U*U*Fz)+                         
                              3*H3_3*(Fx*V*V+U*Fy*V+U*V*Fy)+ 
                              3*H3_4*(Fx*W*W+U*Fz*W+U*W*Fz)+ 
                              3*H3_5*(Fy*W*W+V*Fz*W+V*W*Fz)+ 
                              3*H3_6*(Fy*V*W+V*Fy*W+V*V*Fz)+
                              6*H3_7*(Fx*V*W+U*Fy*W+U*V*Fz) )
    fourth_order = 1./(24*1/9)*(6*H4_1*(Fx*U*V*V+U*Fx*V*V+U*U*Fy*V+U*U*V*Fy)+
                                6*H4_2*(Fx*U*W*W+U*Fx*W*W+U*U*Fz*W+U*U*W*Fz)+
                                6*H4_3*(Fy*V*W*W+V*Fy*W*W+V*V*Fz*W+V*V*W*Fz)+     
                               12*H4_4*(Fx*V*W*W+U*Fy*W*W+U*V*Fz*W+U*V*W*Fz)+
                               12*H4_5*(Fx*V*V*W+U*Fy*V*W+U*V*Fy*W+U*V*V*Fz)+
                               12*H4_6*(Fx*U*V*W+U*Fx*V*W+U*U*Fy*W+U*U*V*Fz) )
    fifth_order = 1./(120*sqrt(1/3)**5)*(30*H5_1*(Fx*U*V*W*W+U*Fx*V*W*W+U*U*Fy*W*W+U*U*V*Fz*W+U*U*V*W*Fz)+
                                30*H5_2*(Fx*U*V*V*W+U*Fx*V*V*W+U*U*Fy*V*W+U*U*V*Fy*W+U*U*V*V*Fz)+
                                30*H5_3*(Fx*V*V*W*W+U*Fy*V*W*W+U*V*Fy*W*W+U*V*V*Fz*W+U*V*V*W*Fz) )
    sixth_order = 1./(720*1/27)*90*H6_1*(Fx*U*V*V*W*W+U*Fx*V*V*W*W+U*U*Fy*V*W*W+U*U*V*Fy*W*W+U*U*V*V*Fz*W+U*U*V*V*W*Fz)

    Force[p,0] = wp[p]*(first_order + second_order + third_order + fourth_order + fifth_order + sixth_order)

printer = DotZeroPrinter()

T = simplify(T)
print("simplified T")

M = simplify(M)
print("simplified M")

# Saito 2023, 10.1103/PhysRevE.108.065305
keq = zeros(NP, 1)
keq[0,0] = rho
keq[9,0] = 3*pressure
keq[17,0] = 1/3*pressure
keq[18,0] = 1/3*pressure
keq[19,0] = 1/3*pressure
keq[26,0] = 1/9*pressure

eq = T.inv()*keq
eq = eq.subs([(U, 0), (V, 0), (W, 0), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
eq = nsimplify(simplify(eq), tolerance=1e-12)
for p in range(NP):
    print("eq[%d] = "%(p) + str(printer.doprint(eq[p,0])) + ";")

kf = T*Force
kf = nsimplify(simplify(kf), tolerance=1e-12)
kf = kf.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
kf = nsimplify(simplify(kf), tolerance=1e-12)

pert = zeros(NP, 1)
for p in range(NP):
    Gc = Gx*cx[p] + Gy*cy[p] + Gz*cz[p]
    pert[p,0] = 9/8*sigma * G * (wp[p]*Gc*Gc/(G*G) - B[p])

C = zeros(NP, 1)
C[7,0] = Symbol("Qx") - Symbol("Qy")
C[8,0] = Symbol("Qx") - Symbol("Qz")
C[9,0] = Symbol("Qx") + Symbol("Qy") + Symbol("Qz")

k_pert = T*pert
k_pert = nsimplify(simplify(k_pert), tolerance=1e-12)
k_pert = k_pert.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8), (G**2, G2), (Gx**2, G2x), (Gy**2, G2y), (Gz**2, G2z)])
k_pert = nsimplify(simplify(k_pert), tolerance=1e-12)
for p in range(NP):
    print("k_pert[%d] = "%(p) + str(printer.doprint(apart(k_pert[p,0], G))) + ";")

k = Matrix(symbols('k[0:27]'))
k_pert = Matrix(symbols('k_pert[0:27]'))
k_star = (eye(NP) - freq)*k + freq*keq + (eye(NP) - 1/2*freq)*kf + (eye(NP) - 1/2*freq)*C + freq*k_pert
k_star = k_star.subs([(U**2, U2), (V**2, V2), (W**2, W2), (cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8)])
k_star = nsimplify(simplify(k_star), tolerance=1e-12)

for p in range(NP):
    print("k_star[%d] = "%(p) + str(printer.doprint(apart(k_star[p,0], G))) + ";")

k_star = Matrix(symbols('k_star[0:27]'))

N_inv = M * T.inv()
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
f = M.inv()*raw
f = nsimplify(simplify(f), tolerance=1e-12)
for p in range(NP):
    print("f2[INDEX_F(i, j, k, %d, n)] = "%(p) + str(printer.doprint(f[p,0])) + ";")