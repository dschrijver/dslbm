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

    double Fx_ext = params->Fx_ext;
    double Fy_ext = params->Fy_ext;
    double Fz_ext = params->Fz_ext;

    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    FOR_DOMAIN
    {
        for (int n = 0; n < NCOMP; n++)
        {
            Fx[INDEX(i, j, k, n)] = Fx_ext;
            Fy[INDEX(i, j, k, n)] = Fy_ext;
            Fz[INDEX(i, j, k, n)] = Fz_ext;
        }
    }
}

void evaluate_force(int i, int j, int k, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;

    double Fx_ext = params->Fx_ext;
    double Fy_ext = params->Fy_ext;
    double Fz_ext = params->Fz_ext;

    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    for (int n = 0; n < NCOMP; n++)
    {
        Fx[INDEX(i, j, k, n)] = Fx_ext;
        Fy[INDEX(i, j, k, n)] = Fy_ext;
        Fz[INDEX(i, j, k, n)] = Fz_ext;
    }
}