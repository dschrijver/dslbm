#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "datatypes.h"

// fields.c
void evaluate_pressure(SimulationBag *sim);

void dims_create(int dims[3], ParamBag *params);
int block_start(int coord, int n, int p);

#endif