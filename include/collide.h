#ifndef COLLIDE_H
#define COLLIDE_H

#include "datatypes.h"

void collide_distributions_MRT(SimulationBag *sim);
void collide_distributions_CGM(SimulationBag *sim);
void extrapolate_wall_density(int i, int j, int k, double *rho_RED, double *rho_BLUE, SimulationBag *sim);

#endif