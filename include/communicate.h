#ifndef COMMUNICATE_H
#define COMMUNICATE_H

#include "datatypes.h"

void communicate_field(double *field, SimulationBag *sim);
void communicate_dist(double *dist, SimulationBag *sim);
void communicate_matrix(matrix *field, SimulationBag *sim);

#endif