from sympy import *

class DotZeroPrinter(StrPrinter):
    def _print_Integer(self, expr):
        return f"{expr}.0"

def m(a, b, c):
    return cx[i]**a * cy[i]**b * cz[i]**c

def t(a, b, c):
    return CX**a * CY**b * CZ**c

U = Symbol("u_i")
U2 = Symbol("u2_i")
V = Symbol("v_i")
V2 = Symbol("v2_i")
W = Symbol("w_i")
W2 = Symbol("w2_i")
rho = Symbol("rho_i")
Fx = Symbol("Fx_i")
Fy = Symbol("Fy_i")
Fz = Symbol("Fz_i")
omega = Symbol("omega")
Qx = Symbol("Qx_i")
Qy = Symbol("Qy_i")
Qz = Symbol("Qz_i")
pressure = Symbol("pressure_i")

NP = 27

# Saito 2023, 10.1103/PhysRevE.108.065305
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1]
cy = [0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1]
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1]
p_bounceback = [0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25]
wp = [8 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216]

freq = diag(*[1, 1, 1, 1, omega, omega, omega, omega, omega, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

T = zeros(NP, NP)
M = zeros(NP, NP)
for i in range(NP):
    # Saito 2023, 10.1103/PhysRevE.108.065305
    M[0,i] = 1

    M[1,i] = m(1,0,0)
    M[2,i] = m(0,1,0)
    M[3,i] = m(0,0,1)

    M[4,i] = m(1,1,0)
    M[5,i] = m(1,0,1)
    M[6,i] = m(0,1,1)
    M[7,i] = m(2,0,0) - m(0,2,0)
    M[8,i] = m(2,0,0) - m(0,0,2)
    M[9,i] = m(2,0,0) + m(0,2,0) + m(0,0,2)

    M[10,i] = m(1,2,0)
    M[11,i] = m(1,0,2)
    M[12,i] = m(0,1,2)
    M[13,i] = m(2,1,0)
    M[14,i] = m(2,0,1)
    M[15,i] = m(0,2,1)
    M[16,i] = m(1,1,1)

    M[17,i] = m(2,2,0)
    M[18,i] = m(2,0,2)
    M[19,i] = m(0,2,2)
    M[20,i] = m(2,1,1)
    M[21,i] = m(1,2,1)
    M[22,i] = m(1,1,2)

    M[23,i] = m(1,2,2)
    M[24,i] = m(2,1,2)
    M[25,i] = m(2,2,1)

    M[26,i] = m(2,2,2) 

    # Saito 2023, 10.1103/PhysRevE.108.065305
    CX = cx[i] - U
    CY = cy[i] - V
    CZ = cz[i] - W

    T[0,i] = 1

    T[1,i] = t(1,0,0)
    T[2,i] = t(0,1,0)
    T[3,i] = t(0,0,1)

    T[4,i] = t(1,1,0)
    T[5,i] = t(1,0,1)
    T[6,i] = t(0,1,1)
    T[7,i] = t(2,0,0) - t(0,2,0)
    T[8,i] = t(2,0,0) - t(0,0,2)
    T[9,i] = t(2,0,0) + t(0,2,0) + t(0,0,2)

    T[10,i] = t(1,2,0)
    T[11,i] = t(1,0,2)
    T[12,i] = t(0,1,2)
    T[13,i] = t(2,1,0)
    T[14,i] = t(2,0,1)
    T[15,i] = t(0,2,1)
    T[16,i] = t(1,1,1)

    T[17,i] = t(2,2,0)
    T[18,i] = t(2,0,2)
    T[19,i] = t(0,2,2)
    T[20,i] = t(2,1,1)
    T[21,i] = t(1,2,1)
    T[22,i] = t(1,1,2)

    T[23,i] = t(1,2,2)
    T[24,i] = t(2,1,2)
    T[25,i] = t(2,2,1)

    T[26,i] = t(2,2,2) 

printer = DotZeroPrinter()

T = simplify(T)
print("simplified T")

M = simplify(M)
print("simplified M")

N_shift = simplify(T*M.inv())
print("simplified N_shift")

# Saito 2023, 10.1103/PhysRevE.108.065305
teq = zeros(NP, 1)
teq[0,0] = rho
teq[9,0] = 3*pressure
teq[17,0] = 1/3*pressure
teq[18,0] = 1/3*pressure
teq[19,0] = 1/3*pressure
teq[26,0] = 1/9*pressure

# Saito 2023, 10.1103/PhysRevE.108.065305
tphi = zeros(NP, 1)
tphi[1,0] = Fx
tphi[2,0] = Fy
tphi[3,0] = Fz
tphi[10,0] = Fx/3
tphi[11,0] = Fx/3
tphi[12,0] = Fy/3
tphi[13,0] = Fy/3
tphi[14,0] = Fz/3
tphi[15,0] = Fz/3
tphi[23,0] = Fx/9
tphi[24,0] = Fy/9
tphi[25,0] = Fz/9

# Saito 2023, 10.1103/PhysRevE.108.065305
tc = zeros(NP, 1)
tc[7,0] = Qx - Qy
tc[8,0] = Qx - Qz
tc[9,0] = Qx + Qy + Qz

meq = N_shift.inv() * teq
meq = nsimplify(simplify(meq), tolerance=1e-12)
meq = meq.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
meq = nsimplify(simplify(meq), tolerance=1e-12)
for i in range(NP):
    print("mf[%d] = "%(i) + str(printer.doprint(meq[i,0])) + ";")

raw = Matrix(symbols("raw[0:27]"))
pop = M.inv() * raw
pop = nsimplify(simplify(pop), tolerance=1e-12)
for i in range(NP):
    print("pop[%d] = "%(i) + str(printer.doprint(pop[i,0])) + ";")

tf = Matrix(symbols('t(0:27)'))
t_star = (eye(NP) - freq)*tf + freq*teq + (eye(NP) - 1/2*freq)*tphi + (eye(NP) - 1/2*freq)*tc
t_star = nsimplify(simplify(t_star), tolerance=1e-12)
t_star = t_star.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
t_star = nsimplify(simplify(t_star), tolerance=1e-12)
for i in range(NP):
    print("t_star[%d] = "%(i) + str(printer.doprint(t_star[i,0])) + ";")

central = Matrix(symbols("central[0:27]"))
raw = N_shift.inv() * central
raw = nsimplify(simplify(raw), tolerance=1e-12)
raw = raw.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
raw = nsimplify(simplify(raw), tolerance=1e-12)
for i in range(NP):
    print("raw[%d] = "%(i) + str(printer.doprint(raw[i,0])) + ";")