import numpy as np
from sympy import *

# cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0]
# cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1]
# cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1]
# wp = [1 / 3, 1 / 18, 1 / 18, 1 / 18, 1 / 18, 1 / 18, 1 / 18, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36]

cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1]
cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1]
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1]

wp = [8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0]

c = np.array([cx, cy, cz])

# NP = 19
NP = len(cx)

Q = np.zeros((3,3,3,3))
Q_matrix = np.zeros((6,6))

q = np.zeros((3,3,3))

for a in range(3):
    for b in range(3):
        for mu in range(3):
            for nu in range(3):
                for i in range(NP):
                    if cy[i] > 0:
                        Q[a, b, mu, nu] += -wp[i]*c[a,i]*c[b,i]*c[mu,i]*c[nu,i]

for a in range(3):
    for mu in range(3):
        for nu in range(3):
            for i in range(NP):
                if cy[i] > 0:
                    q[a, mu, nu] += -wp[i]*c[a,i]*c[mu,i]*c[nu,i]

print(q)

indices = [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)]

for i in range(6):
    for j in range(6):
        Q_matrix[i,j] = Q[*(indices[i]), *(indices[j])]

b = symbols('b_00 b_11 b_22 b_01 b_02 b_12')
B = Matrix(6, 1, b)

M_flat = nsimplify(Matrix(Q_matrix).inv() * B, tolerance=1e-12)
print(M_flat)
M = zeros(3,3)

for k in range(6):
    M[*(indices[k])] = M_flat[k]

for j in range(3):
    for i in range(j,3):
        M[i,j] = M[j,i]

print(M)

mass_loss = 0
for i in range(NP):
    for mu in range(3):
            for nu in range(3):
                if cy[i] > 0:
                    mass_loss += M[mu, nu]*c[mu,i]*c[nu,i]*wp[i]
print(nsimplify(mass_loss, tolerance=1e-12))

m = zeros(3,1)
for a in range(3):
    for mu in range(3):
        for nu in range(3):
            for i in range(NP):
                if cy[i] > 0:
                    m[a] += -wp[i]*c[a,i]*c[mu,i]*c[nu,i]*M[mu,nu]

print(nsimplify(m, tolerance=1e-12))

n = zeros(3,3)
N_vec = Matrix(symbols('n0:3'))
for a in range(3):
    for b in range(3):
        for mu in range(3):
            for i in range(NP):
                if cy[i] > 0:
                    n[a,b] += -wp[i]*c[a,i]*c[b,i]*c[mu,i]*N_vec[mu]

print(nsimplify(n, tolerance=1e-12))