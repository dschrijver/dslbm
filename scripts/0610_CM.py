from sympy import StrPrinter, Symbol, zeros, simplify, nsimplify, sqrt, diag, pi, exp, diff, tensorproduct, Array, Matrix, symbols, eye, solve

class DotZeroPrinter(StrPrinter):
    def _print_Integer(self, expr):
        return f"{expr}.0"

def Hf(nx, ny, nz):
    n = nx + ny + nz

    F = Array([Fx, Fy, Fz])
    vel = Array([U, V, W])

    Fu = tensorproduct(F, *((n-1)*[vel]))
    for k in range(1, n):
        Fu += tensorproduct(*(k*[vel]), F, *((n-1-k)*[vel]))

    index = tuple([0]*nx + [1]*ny + [2]*nz)

    return H(nx, ny, nz) * Fu[index]

def Hu(nx, ny, nz):
    return H(nx, ny, nz) * U**nx * V**ny * W**nz

def H(nx, ny, nz):
    return H_1D(nx, cx[i])*H_1D(ny, cy[i])*H_1D(nz, cz[i])

def H_1D(n, c_i):
    x = Symbol("x")
    return (-1/3)**n / weight(c_i) * diff(weight(x), x, n).subs(x, c_i)

def weight(x):
    return 1/sqrt(2*pi*1/3)*exp(-3/2*x**2)

def m(a, b, c):
    return cx[i]**a * cy[i]**b * cz[i]**c

def t(a, b, c):
    return CX**a * CY**b * CZ**c

rho = Symbol("rho")
U = Symbol("u")
U2 = Symbol("u2")
V = Symbol("v")
V2 = Symbol("v2")
W = Symbol("w")
W2 = Symbol("w2")
pressure = Symbol("pressure")
Fx = Symbol("Fx")
Fy = Symbol("Fy")
Fz = Symbol("Fz")
omega = Symbol("omega")
Qx = Symbol("Qx")
Qy = Symbol("Qy")
Qz = Symbol("Qz")

NP = 27
cx = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1]
cy = [0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1]
cz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1]
wp = [8 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 2 / 27, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 54, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216, 1 / 216]
Omega = diag(*[1, 1, 1, 1, omega, omega, omega, omega, omega, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

M = zeros(NP, NP)
T = zeros(NP, NP)
feq = zeros(NP, 1)
S_guo = zeros(NP, 1)
for i in range(NP):
    # De Rosis 2019, 10.1063/1.5124719
    first_order = 3 * (Hu(1,0,0) + Hu(0,1,0) + Hu(0,0,1))
    second_order = 9/2 * (Hu(2,0,0) + Hu(0,2,0) + Hu(0,0,2) + 2 * (Hu(1,1,0) + Hu(1,0,1) + Hu(0,1,1)))
    third_order = 27/2 * (Hu(2,1,0) + Hu(2,0,1) + Hu(1,2,0) + Hu(1,0,2) + Hu(0,1,2) + Hu(0,2,1) + 2 * Hu(1,1,1))
    fourth_order = 81/4 * (Hu(2,2,0) + Hu(2,0,2) + Hu(0,2,2) + 2 * (Hu(1,1,2) + Hu(1,2,1) + Hu(2,1,1)))
    fifth_order = 243/4 * (Hu(2,1,2) + Hu(2,2,1) + Hu(1,2,2))
    sixth_order = 729/8 * Hu(2,2,2)
    geq = wp[i]*rho*(1 + first_order + second_order + third_order + fourth_order + fifth_order + sixth_order)

    # Saito 2023, 10.1103/PhysRevE.108.065305
    second_order = 9/2 * (H(2,0,0) + H(0,2,0) + H(0,0,2))
    fourth_order = -27/4 * (H(2,2,0) + H(0,2,2) + H(2,0,2))
    sixth_order = 81/8 * H(2,2,2)
    E = wp[i] * (second_order + fourth_order + sixth_order)

    # Saito 2023, 10.1103/PhysRevE.108.065305
    third_order = 27/2 * (U*(H(1,2,0) + H(1,0,2)) + V*(H(2,1,0) + H(0,1,2)) + W*(H(2,0,1) + H(0,2,1)))
    fourth_order = 81/4 * ((U*U + V*V)*H(2,2,0) + (V*V + W*W)*H(0,2,2) + (U*U + W*W)*H(2,0,2) + 2*(V*W*H(2,1,1) + U*W*H(1,2,1) + U*V*H(1,1,2)))
    fifth_order = 243/4 * (U*(V*V + W*W - 1/3)*H(1,2,2) + V*(U*U + W*W - 1/3)*H(2,1,2) + W*(U*U + V*V - 1/3)*H(2,2,1))
    sixth_order = 729/8 * (U*U*V*V + V*V*W*W + U*U*W*W - 1/3*(U*U + V*V + W*W))*H(2,2,2)
    Phi = wp[i] * (third_order + fourth_order + fifth_order + sixth_order)

    feq[i,0] = geq + (pressure - rho / 3)*(E + Phi)

    # De Rosis 2019, 10.1063/1.5124719
    first_order = 3 * (Hf(1,0,0) + Hf(0,1,0) + Hf(0,0,1))
    second_order = 9/2 * (Hf(2,0,0) + Hf(0,2,0) + Hf(0,0,2) + 2 * (Hf(1,1,0) + Hf(1,0,1) + Hf(0,1,1)))
    third_order = 27/2 * (Hf(2,1,0) + Hf(2,0,1) + Hf(1,2,0) + Hf(1,0,2) + Hf(0,1,2) + Hf(0,2,1) + 2 * Hf(1,1,1))
    fourth_order = 81/4 * (Hf(2,2,0) + Hf(2,0,2) + Hf(0,2,2) + 2 * (Hf(1,1,2) + Hf(1,2,1) + Hf(2,1,1)))
    fifth_order = 243/4 * (Hf(2,1,2) + Hf(2,2,1) + Hf(1,2,2))
    sixth_order = 729/8 * Hf(2,2,2)
    # Removed rho from Eq.(9) from De Rosis 2019, 10.1063/1.5124719
    S_guo[i,0] = wp[i]*(first_order + second_order + third_order + fourth_order + fifth_order + sixth_order)

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
        
M = simplify(M)
T = simplify(T)
N_shift = T * M.inv()
N_shift = simplify(N_shift)

print("Simplifying equilibrium and forcing... ", end="", flush=True)
feq = simplify(feq)
S_guo = simplify(S_guo)
print("Done!")

sum = 0.0
for i in range(NP):
    if (cx[i] == 0):
        sum += S_guo[i]*cy[i]
print(nsimplify(simplify(sum), tolerance=1e-12))

sum = 0.0
for i in range(NP):
    # if (cx[i] == 0):
    sum += S_guo[i]*cy[i]
print(nsimplify(simplify(sum), tolerance=1e-12))

printer = DotZeroPrinter()

teq = T*feq
teq = nsimplify(simplify(teq), tolerance=1e-12)
teq = teq.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
teq = nsimplify(simplify(teq), tolerance=1e-12)
print("\nEquilibrium in central moment space:\n")
for i in range(NP):
    print("teq[%d] = "%(i) + str(printer.doprint(teq[i,0])) + ";")

t_s = T * S_guo
t_s = nsimplify(simplify(t_s), tolerance=1e-12)
t_s = t_s.subs([(U**2, U2), (V**2, V2), (W**2, W2)])
t_s = nsimplify(simplify(t_s), tolerance=1e-12)
print("\nForce source in central space:\n")
for i in range(NP):
    print("t_s[%d] = "%(i) + str(printer.doprint(t_s[i,0])) + ";")

t_c = zeros(NP, 1)
t_c[7,0] = Qx - Qy
t_c[8,0] = Qx - Qz
t_c[9,0] = Qx + Qy + Qz

tf = Matrix(symbols("t(0:27)"))
t_star = tf - Omega*(tf - teq) + (eye(NP) - Omega/2)*t_s + (eye(NP) - Omega/2)*t_c
t_star = nsimplify(simplify(t_star), tolerance=1e-12)
print("\nCollision in central space:\n")
for i in range(NP):
    print("t_star[%d] = "%(i) + str(printer.doprint(t_star[i,0])) + ";")