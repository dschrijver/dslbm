#include <math.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"
#include "../include/forcing.h"

void evaluate_forces(SimulationBag *sim)
{
    ParamBag *params = sim->params;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    FOR_DOMAIN
    {
        evaluate_force(i, j, k, sim);
    }
}

void evaluate_force(int i, int j, int k, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    double rho_i;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;

    double gx = params->gx;
    double gy = params->gy;
    double gz = params->gz;

    double *rho_comp = comp_fields->rho_comp;

    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    for (int n = 0; n < NCOMP; n++)
    {
        rho_i = rho_comp[INDEX(i, j, k, n)];

        // Fx[INDEX(i, j, k, n)] = rho_i * gx;
        // Fy[INDEX(i, j, k, n)] = rho_i * gy;
        // Fz[INDEX(i, j, k, n)] = rho_i * gz;

        Fx[INDEX(i, j, k, n)] = gx;
        Fy[INDEX(i, j, k, n)] = gy;
        Fz[INDEX(i, j, k, n)] = gz;
    }
}