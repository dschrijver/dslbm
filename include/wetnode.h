#ifndef WETNODE_H
#define WETNODE_H

#include "datatypes.h"

void wetnode_boundary_conditions(SimulationBag *sim);
void wetnode_compute_density(int i, int j, int k, int nx, int ny, int nz, SimulationBag *sim);
void wetnode_compute_velocity(int i, int j, int k, int nx, int ny, int nz, SimulationBag *sim);
void wetnode_mass_conservation_streaming(int i, int j, int k, int nx, int ny, int nz, SimulationBag *sim);
void non_equilibrium_bounce_back_x(int i, int nx, SimulationBag *sim);
void non_equilibrium_bounce_back_y(int j, int ny, SimulationBag *sim);
void non_equilibrium_bounce_back_z(int k, int nz, SimulationBag *sim);

#endif