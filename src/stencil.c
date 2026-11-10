#include <string.h>

#include "../definitions.h"
#include "../include/stencil.h"

void initialize_stencil(SimulationBag *sim)
{
    Stencil *stencil = sim->stencil;

#ifdef   DIM_2D 
    stencil->NP = 9;
#else // DIM_3D
    stencil->NP = 27;
#endif

    allocate_stencil(sim);

#ifdef   DIM_2D 
    int cx[] = {0, 1, -1, 0, 0, 1, -1, -1, 1};
    int cy[] = {0, 0, 0, 1, -1, 1, -1, 1, -1};
    int cz[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double wp[] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    int p_bounceback[] = {0, 2, 1, 4, 3, 6, 5, 8, 7};
#else // DIM_3D
    int cx[] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1};
    int cy[] = {0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1};
    int cz[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1};
    double wp[] = {8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0};
    int p_bounceback[] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
#endif

    memcpy(stencil->cx, cx, stencil->NP * sizeof(int));
    memcpy(stencil->cy, cy, stencil->NP * sizeof(int));
    memcpy(stencil->cz, cz, stencil->NP * sizeof(int));
    memcpy(stencil->wp, wp, stencil->NP * sizeof(double));
    memcpy(stencil->p_bounceback, p_bounceback, stencil->NP * sizeof(int));
}