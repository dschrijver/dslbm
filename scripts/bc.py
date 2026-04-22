from sympy import *

class DotZeroPrinter(StrPrinter):
    def _print_Integer(self, expr):
        return f"{expr}.0"

# De Rosis 2019, 10.1063/1.5124719
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1]  
cy = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1] 
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1]
wp = [8 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216]
NP = 27
cs2 = 1/3

x = Symbol("x")
y = Symbol("y")
z = Symbol("z")
t = Symbol("t")
rho0 = Function('rho')(x, y, z, t)
rho1 = Function('rho1')(x, y, z, t)
u1 = Function('u1')(x, y, z, t)
v1 = Function('v1')(x, y, z, t)
w1 = Function('w1')(x, y, z, t)
tau = Symbol("tau")
Delta_t = Symbol("Delta_t")
Delta_x = Symbol("Delta_x")

f0 = zeros(NP, 1)
for i in range(NP):
    f0[i] = wp[i]*rho0
f0 = nsimplify(simplify(f0), tolerance=1e-12)

f1 = zeros(NP, 1)
for i in range(NP):
    f1[i] = wp[i]*(rho1 + rho0*(u1*cx[i] + v1*cy[i] + w1*cz[i])/cs2) - tau/Delta_t*(Delta_x*(cx[i]*Derivative(f0, x) + cy[i]*Derivative(f0, y) + cz[i]*Derivative(f0, z)))
f1 = nsimplify(simplify(f1), tolerance=1e-12)

print(f1)