#include "../include/datatypes.h"
#include "../include/forcing.h"

void evaluate_forces(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    double rho_i;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double gx = params->gx;
    double gy = params->gy;
    double gz = params->gz;

    double *rho_comp = comp_fields->rho_comp;

    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    FOR_DOMAIN
    {
        for (int n = 0; n < NCOMP; n++)
        {
            rho_i = rho_comp[INDEX(i, j, k, n)];

            Fx[INDEX(i, j, k, n)] = rho_i * gx;
            Fy[INDEX(i, j, k, n)] = rho_i * gy;
            Fz[INDEX(i, j, k, n)] = rho_i * gz;
        }
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

        Fx[INDEX(i, j, k, n)] = rho_i * gx;
        Fy[INDEX(i, j, k, n)] = rho_i * gy;
        Fz[INDEX(i, j, k, n)] = rho_i * gz;
    }
}