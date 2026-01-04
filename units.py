import numpy as np

# The RED fluid is denser than the BLUE fluid. 
# Reynolds number and Weber number are based on the RED fluid. 
# Ratios are RED/BLUE fluid.

# Dimensionless quantities
Re = 10         # = U L / nu_RED
We = 10         # = rho_RED U^2 L / sigma
Fr = 10         # = U / sqrt(g L)
rho_ratio = 10  # = rho_RED / rho_BLUE
nu_ratio = 10   # = nu_RED / nu_BLUE

# Choose length scale. 
L = 15 # Lattice Units

# LBM choices
tau_RED = 1 # Lattice Units
rho_RED = 1 # Lattice Units
print("Chosen lattice values:")
print("L =", L)
print("tau_RED =", tau_RED)
print("rho_RED =", rho_RED)


# Computed lattice values
# --- DON'T CHANGE THIS ---
print("\nComputed lattice values:")
rho_BLUE = rho_RED / rho_ratio
print("rho_BLUE =", rho_BLUE)
nu_RED = 1.0/3.0 * (tau_RED - 0.5)
nu_BLUE = nu_RED / nu_ratio
tau_BLUE = 3.0*nu_BLUE + 0.5
print("tau_BLUE =", tau_BLUE)
U = Re * nu_RED / L
sigma = rho_RED * U**2 * L / We
print("sigma = %.2e"%sigma)
g = U**2 / (Fr**2 * L)
print("g = %.2e"%g)
# -------------------------

# Physical choices
nu_RED_phys = 1.004e-6  # m^2/s
rho_RED_phys = 1000   # kg/m^3
print("\nChosen physical values:")
print("nu_RED_phys =", nu_RED_phys, "m^2/s")
print("rho_RED_phys =", rho_RED_phys, "kg/m^3")

# Computed physical values
# --- DON'T CHANGE THIS ---
print("\nComputed physical values:")
rho_BLUE_phys = rho_RED_phys / rho_ratio
print("rho_BLUE_phys =", rho_BLUE_phys, "kg/m^3")
nu_BLUE_phys = nu_RED_phys / nu_ratio
print("nu_BLUE_phys =", nu_BLUE_phys, "m^2/s")
# -------------------------