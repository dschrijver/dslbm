#include "../include/datatypes.h"
#include "../include/forcing.h"

void evaluate_forces(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    FOR_DOMAIN
    {
        for (int n = 0; n < NCOMP; n++)
        {
            Fx[INDEX(i, j, k, n)] = 0.0;
            Fy[INDEX(i, j, k, n)] = 0.0;
            Fz[INDEX(i, j, k, n)] = 0.0;
        }
    }
}