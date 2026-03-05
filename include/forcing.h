#ifndef FORCING_H
#define FORCING_H

#include "datatypes.h"

void evaluate_forces(SimulationBag *sim);
void evaluate_force(int i, int j, int k, SimulationBag *sim);
void evaluate_shan_chen_force(int i, int j, int k, SimulationBag *sim);
void set_wall_densities_angle_y(int i, int j, int k, int ny, double theta_c, double *rho_RED, double *rho_BLUE, SimulationBag *sim);

#endif