from sympy import *

# De Rosis 2019, 10.1063/1.5124719
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1]
cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1]
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1]
wp = [8 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216]
p_bounceback = [0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15, 26, 25, 24, 23, 22, 21, 20, 19]
NP = 27


rho = Symbol("rho_i")
U = Symbol("u_i")
U2 = Symbol("u2_i")
V = Symbol("v_i")
V2 = Symbol("v2_i")
W = Symbol("w_i")
W2 = Symbol("w2_i")
cs2 = Symbol("cs2_i")
cs4 = Symbol("cs4_i")
cs6 = Symbol("cs6_i")
cs8 = Symbol("cs8_i")

T = zeros(NP, NP)
M = zeros(NP, NP)

for i in range(NP):
    # M Matrix
    M[i,0] = 1
    M[i,1] = cx[i]
    M[i,2] = cy[i]
    M[i,3] = cz[i]
    M[i,4] = cx[i]*cy[i]
    M[i,5] = cx[i]*cz[i]
    M[i,6] = cy[i]*cz[i]
    M[i,7] = cx[i]*cx[i] - cy[i]*cy[i]
    M[i,8] = cx[i]*cx[i] - cz[i]*cz[i]
    M[i,9] = cx[i]*cx[i] + cy[i]*cy[i] + cz[i]*cz[i]
    M[i,10] = cx[i]*cy[i]*cy[i]+cx[i]*cz[i]*cz[i]
    M[i,11] = cx[i]*cx[i]*cy[i]+cy[i]*cz[i]*cz[i]
    M[i,12] = cx[i]*cx[i]*cz[i]+cy[i]*cy[i]*cz[i]
    M[i,13] = cx[i]*cy[i]*cy[i]-cx[i]*cz[i]*cz[i]
    M[i,14] = cx[i]*cx[i]*cy[i]-cy[i]*cz[i]*cz[i]
    M[i,15] = cx[i]*cx[i]*cz[i]-cy[i]*cy[i]*cz[i]
    M[i,16] = cx[i]*cy[i]*cz[i]
    M[i,17] = cx[i]*cx[i]*cy[i]*cy[i]+cx[i]*cx[i]*cz[i]*cz[i]+cy[i]*cy[i]*cz[i]*cz[i]
    M[i,18] = cx[i]*cx[i]*cy[i]*cy[i]+cx[i]*cx[i]*cz[i]*cz[i]-cy[i]*cy[i]*cz[i]*cz[i]
    M[i,19] = cx[i]*cx[i]*cy[i]*cy[i]-cx[i]*cx[i]*cz[i]*cz[i]
    M[i,20] = cx[i]*cx[i]*cy[i]*cz[i]
    M[i,21] = cx[i]*cy[i]*cy[i]*cz[i]
    M[i,22] = cx[i]*cy[i]*cz[i]*cz[i]
    M[i,23] = cx[i]*cy[i]*cy[i]*cz[i]*cz[i]
    M[i,24] = cx[i]*cx[i]*cy[i]*cz[i]*cz[i]
    M[i,25] = cx[i]*cx[i]*cy[i]*cy[i]*cz[i]
    M[i,26] = cx[i]*cx[i]*cy[i]*cy[i]*cz[i]*cz[i]

    # T Matrix
    CX = cx[i] - U
    CY = cy[i] - V
    CZ = cz[i] - W

    CX2 = CX*CX
    CY2 = CY*CY
    CZ2 = CZ*CZ

    T[i,0] = 1
    T[i,1] = CX
    T[i,2] = CY
    T[i,3] = CZ
    T[i,4] = CX*CY
    T[i,5] = CX*CZ
    T[i,6] = CY*CZ
    T[i,7] = CX2 - CY2
    T[i,8] = CX2 - CZ2
    T[i,9] = CX2 + CY2 + CZ2
    T[i,10] = CX*CY2+CX*CZ2
    T[i,11] = CX2*CY+CY*CZ2
    T[i,12] = CX2*CZ+CY2*CZ
    T[i,13] = CX*CY2-CX*CZ2
    T[i,14] = CX2*CY-CY*CZ2
    T[i,15] = CX2*CZ-CY2*CZ
    T[i,16] = CX*CY*CZ
    T[i,17] = CX2*CY2+CX2*CZ2+CY2*CZ2
    T[i,18] = CX2*CY2+CX2*CZ2-CY2*CZ2
    T[i,19] = CX2*CY2-CX2*CZ2
    T[i,20] = CX2*CY*CZ
    T[i,21] = CX*CY2*CZ
    T[i,22] = CX*CY*CZ2
    T[i,23] = CX*CY2*CZ2
    T[i,24] = CX2*CY*CZ2
    T[i,25] = CX2*CY2*CZ
    T[i,26] = CX2*CY2*CZ2

keq = zeros(NP, 1)
keq[0] = rho
keq[9] = 3*rho*cs2
keq[17] = rho*cs2
keq[18] = rho*cs2**2
keq[26] =  rho*cs2**3

eq = (T.T).inv()*keq
eq = eq.subs([(cs2**2, cs4), (cs2**3, cs6), (cs2**4, cs8), (U**2, U2), (V**2, V2), (W**2, W2)])
eq = nsimplify(simplify(eq), tolerance=1e-12)

sum = 0
for i in range(NP):
    if cx[i] > 0:
        sum += eq[i]-eq[p_bounceback[i]]
print(nsimplify(simplify(sum), tolerance=1e-12))

sum = 0
for i in range(NP):
    if cx[i] > 0:
        sum += (eq[i]-eq[p_bounceback[i]])*cy[i]
print(nsimplify(simplify(sum), tolerance=1e-12))