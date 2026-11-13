#ifndef COLLIDE_H
#define COLLIDE_H

#include "datatypes.h"

void compute_equilibrium(double rho, double u, double v, double w, double *feq, SimulationBag *sim);
void collide(SimulationBag *sim);

#endif