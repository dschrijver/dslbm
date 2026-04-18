#ifndef COLLIDE_H
#define COLLIDE_H

#include "datatypes.h"

void compute_equilibrium(double rho, double u, double v, double w, double cs2, double *feq, SimulationBag *sim);
void evaluate_color_gradients(SimulationBag *sim);
double extrapolate_wall_rho_N(int i, int j, int k, SimulationBag *sim);
void collide(SimulationBag *sim);
void compute_stationary_equilibrium(double rho, double cs2, double *feq);

#endif