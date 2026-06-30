from sympy import *

def P_RED(rho, T_ratio):
    a = 2.0/49.0
    b = 2.0/21.0
    R = 1.0

    T_c = 0.0778/0.45724 * a/b * 1.0/R

    T = T_ratio*T_c

    omega = 0.344 # For water
    alpha = 1.0 + (0.37464 + 1.54226*omega - 0.26992*omega*omega)*(1.0 - sqrt(T_ratio))
    alpha *= alpha

    return rho*R*T/(1.0 - b*rho) - a*alpha*rho*rho/(1.0 + 2.0*b*rho - b*b*rho*rho)

rho_0_RED = Symbol("rho_0_RED")
T_ratio = Symbol("T_ratio")

rho_0_BLUE = Symbol("")

print(P_RED(rho_0_RED, T_ratio))

# solve()