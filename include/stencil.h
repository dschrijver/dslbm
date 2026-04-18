#ifndef STENCIL_H
#define STENCIL_H

#include <string.h>
#include <math.h>

#include "datatypes.h"
#include "memory.h"

static inline void initialize_stencil(SimulationBag *sim)
{
    Stencil *stencil = sim->stencil;
    ParamBag *params = sim->params;

    stencil->NP = 27;

    allocate_stencil(sim);

    stencil->zeta = 9.0 / 19.0;
    stencil->cs2[RED] = stencil->zeta * (1.0 - params->alpha_RED);
    stencil->cs2[BLUE] = stencil->zeta * (1.0 - params->alpha_BLUE);
    stencil->tau[RED] = 1.0 / stencil->cs2[RED] * params->nu_RED + 0.5;
    stencil->tau[BLUE] = 1.0 / stencil->cs2[BLUE] * params->nu_BLUE + 0.5;

    stencil->C_norm = 1.0 / 6.0;
    stencil->C_par = 1.0 / 18.0;

    // De Rosis 2019, 10.1063/1.5124719
    int cx[] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1};
    int cy[] = {0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1};
    int cz[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1};
    double wp[] = {8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0};
    int p_bounceback[] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15, 26, 25, 24, 23, 22, 21, 20, 19};

    memcpy(stencil->cx, cx, stencil->NP * sizeof(int));
    memcpy(stencil->cy, cy, stencil->NP * sizeof(int));
    memcpy(stencil->cz, cz, stencil->NP * sizeof(int));
    memcpy(stencil->wp, wp, stencil->NP * sizeof(double));
    memcpy(stencil->p_bounceback, p_bounceback, stencil->NP * sizeof(int));
}

#endif