import numpy as np
from sympy import symbols, Matrix, nsimplify, zeros

# cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0]
# cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1]
# cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1]
# wp = [1 / 3, 1 / 18, 1 / 18, 1 / 18, 1 / 18, 1 / 18, 1 / 18, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36, 1 / 36]

# cx = [0, 1, -1, 0, 0, 1, -1, 1, -1]
# cy = [0, 0, 0, 1, -1, 1, 1, -1, -1]
# cz = [0, 0, 0, 0, 0, 0, 0, 0, 0]
# wp = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]

cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1]
cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1]
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1]

wp = [8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0]

c = np.array([cx, cy, cz])

NP = len(wp)

result = 0
for i in range(NP):
    if cy[i] > 0:
        result += wp[i]*c[0,i]*c[1,i]
print(result)

# # NP = 19


# Q = np.zeros((3,3,3,3))
# Q_matrix = np.zeros((6,6))

# for alpha in range(3):
#     for beta in range(3):
#         for mu in range(3):
#             for nu in range(3):
#                 for i in range(NP):
#                     if cy[i] > 0:
#                         Q[alpha, beta, mu, nu] += -wp[i]*c[alpha,i]*c[beta,i]*c[mu,i]*c[nu,i]

# indices = [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)]

# for i in range(6):
#     for j in range(6):
#         Q_matrix[i,j] = Q[*(indices[i]), *(indices[j])]

# a = symbols('a_00 a_11 a_22 a_01 a_02 a_12')
# A = Matrix(6, 1, a)

# M_flat = nsimplify(Matrix(Q_matrix).inv() * A, tolerance=1e-12)
# M = zeros(3,3)

# for k in range(6):
#     M[*(indices[k])] = M_flat[k]

# for j in range(3):
#     for i in range(j,3):
#         M[i,j] = M[j,i]

# print("M = \n", M)

# mass_loss = 0
# for i in range(NP):
#     for mu in range(3):
#             for nu in range(3):
#                 if cy[i] > 0:
#                     mass_loss += M[mu, nu]*c[mu,i]*c[nu,i]*wp[i]
# print("Mass loss =", nsimplify(mass_loss, tolerance=1e-12))

# m = zeros(3,1)
# for a in range(3):
#     for mu in range(3):
#         for nu in range(3):
#             for i in range(NP):
#                 if cy[i] > 0:
#                     m[a] += -wp[i]*c[a,i]*c[mu,i]*c[nu,i]*M[mu,nu]

# print(nsimplify(m, tolerance=1e-12))

# n = zeros(3,3)
# N_vec = Matrix(symbols('n0:3'))
# for a in range(3):
#     for b in range(3):
#         for mu in range(3):
#             for i in range(NP):
#                 if cy[i] > 0:
#                     n[a,b] += -wp[i]*c[a,i]*c[b,i]*c[mu,i]*N_vec[mu]

# print(nsimplify(n, tolerance=1e-12))

# n_vec = zeros(3,1)
# for a in range(3):
#     for b in range(3):
#         for i in range(NP):
#             if cy[i] > 0:
#                 n_vec[a] += -wp[i]*c[a,i]*c[b,i]*N_vec[b]

# print(nsimplify(n_vec, tolerance=1e-12))

# vec1 = zeros(3,1)
# for a in range(3):
#     for i in range(NP):
#         if cy[i] > 0:
#             vec1[a] += wp[i]*c[a,i]
# print(nsimplify(vec1, tolerance=1e-12))

# mass_loss = 0
# for i in range(NP):
#     for alpha in range(3):
#         if cy[i] > 0:
#             mass_loss += N_vec[alpha]*c[alpha,i]*wp[i]
# print("Mass loss =", nsimplify(mass_loss, tolerance=1e-12))