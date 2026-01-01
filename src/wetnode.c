#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/wetnode.h"

void wetnode_boundary_conditions(SimulationBag *sim)
{
    DistributionBag *dists = sim->dists;
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

#ifdef BOTTOM_NEBB_NOSLIP
    // Set bottom velocities
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            u[INDEX_GLOB(i, 0, k)] = 0.0;
            v[INDEX_GLOB(i, 0, k)] = 0.0;
            w[INDEX_GLOB(i, 0, k)] = 0.0;
        }
    }
#endif

#ifdef TOP_NEBB_NOSLIP
    // Set bottom velocities
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            u[INDEX_GLOB(i, NY - 1, k)] = 0.0;
            v[INDEX_GLOB(i, NY - 1, k)] = 0.0;
            w[INDEX_GLOB(i, NY - 1, k)] = 0.0;
        }
    }
#endif
}