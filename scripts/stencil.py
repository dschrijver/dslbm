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
    if cy[i] == 0:
        result += wp[i]
print(result)

N = np.zeros(3)
for i in range(NP):
    for alpha in range(3):
        if cy[i] > 0:
            N[alpha] += wp[i]*c[alpha,i]
print(N)

M = np.zeros((3,3))
for i in range(NP):
    for alpha in range(3):
        for beta in range(3):
            if cy[i] > 0:
                M[alpha,beta] += wp[i]*c[alpha,i]*c[beta,i]
print(M)

Q = np.zeros((3,3,3))
for i in range(NP):
    for alpha in range(3):
        for beta in range(3):
            for gamma in range(3):
                if cy[i] == 0:
                    Q[alpha,beta,gamma] += wp[i]*c[alpha,i]*c[beta,i]*c[gamma,i]

print(Q)