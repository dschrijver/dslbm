#ifndef WETNODE_H
#define WETNODE_H

#include "datatypes.h"

void wetnode_macroscopic_fields(SimulationBag *sim);
void wetnode_distributions(SimulationBag *sim);
void wetnode_mass_conservation_streaming(const int i, const int j, const int k, const int nx, const int ny, const int nz, SimulationBag *sim);
void wetnode_compute_density(const int i, const int j, const int k, const int nx, const int ny, const int nz, SimulationBag *sim);
void wetnode_compute_velocity(const int i, const int j, const int k, const int nx, const int ny, const int nz, SimulationBag *sim);
void non_equilibrium_bounce_back_x(int i, int nx, SimulationBag *sim);
void non_equilibrium_bounce_back_y(int j, int ny, SimulationBag *sim);
void non_equilibrium_bounce_back_z(int k, int nz, SimulationBag *sim);

#endif