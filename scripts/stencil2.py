from fractions import Fraction
import numpy as np

# cx = [0, 1, -1, 0, 0, 1, -1, 1, -1]
# cy = [0, 0, 0, 1, -1, 1, 1, -1, -1]
# cz = [0, 0, 0, 0, 0, 0, 0, 0, 0]
# wp = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
# NP = 9

# De Rosis 2019, 10.1063/1.5124719
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1]  
cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1] 
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1]
wp = [8 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216]
NP = 27
c = np.array([cx, cy, cz]).T

sum = 0
for i in range(NP):
    if (cx[i] > 0):
        sum += wp[i]
print("Incoming: ", Fraction(sum).limit_denominator(100))

sum = 0
for i in range(NP):
    if (cx[i] == 0):
        sum += wp[i]
print("Parallel + zero:", Fraction(sum).limit_denominator(100))

N = np.zeros(3)
for i in range(NP):
    if (cx[i] > 0):
        N += wp[i]*c[i]
print("Incoming Nx: ", Fraction(N[0]).limit_denominator(100))
print("Incoming Ny: ", Fraction(N[1]).limit_denominator(100))
print("Incoming Nz: ", Fraction(N[2]).limit_denominator(100))