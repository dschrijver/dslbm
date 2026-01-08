#include "../include/datatypes.h"
#include "../include/fields.h"

void extract_moments(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    double rho_RED_i, rho_BLUE_i, rho_i, u_i, v_i, w_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    double zeta = stencil->zeta;
    double alpha_RED = params->alpha_RED;
    double alpha_BLUE = params->alpha_BLUE;

    double *rho = glob_fields->rho;
    double *pressure = glob_fields->pressure;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *rho_comp = comp_fields->rho_comp;

    double *f1 = dists->f1;

    FOR_DOMAIN
    {
        rho_RED_i = 0.0;
        rho_BLUE_i = 0.0;
        u_i = 0.0;
        v_i = 0.0;
        w_i = 0.0;

        for (int p = 0; p < NP; p++)
        {
            rho_RED_i += f1[INDEX_F(i, j, k, p, RED)];
            rho_BLUE_i += f1[INDEX_F(i, j, k, p, BLUE)];
            u_i += (f1[INDEX_F(i, j, k, p, RED)] + f1[INDEX_F(i, j, k, p, BLUE)]) * (double)cx[p];
            v_i += (f1[INDEX_F(i, j, k, p, RED)] + f1[INDEX_F(i, j, k, p, BLUE)]) * (double)cy[p];
            w_i += (f1[INDEX_F(i, j, k, p, RED)] + f1[INDEX_F(i, j, k, p, BLUE)]) * (double)cz[p];
        }

        rho_i = rho_RED_i + rho_BLUE_i;
        rho[INDEX_GLOB(i, j, k)] = rho_i;
        rho_comp[INDEX(i, j, k, RED)] = rho_RED_i;
        rho_comp[INDEX(i, j, k, BLUE)] = rho_BLUE_i;

        pressure[INDEX_GLOB(i, j, k)] = rho_comp[INDEX(i, j, k, RED)] * zeta * (1.0 - alpha_RED) + rho_comp[INDEX(i, j, k, BLUE)] * zeta * (1.0 - alpha_BLUE);

        u[INDEX_GLOB(i, j, k)] = u_i / rho_i;
        v[INDEX_GLOB(i, j, k)] = v_i / rho_i;
        w[INDEX_GLOB(i, j, k)] = w_i / rho_i;
    }
}

void update_final_velocity(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    double rho_i, Fx_i, Fy_i, Fz_i;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *rho = glob_fields->rho;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX_GLOB(i, j, k)];

        Fx_i = Fx[INDEX(i, j, k, RED)] + Fx[INDEX(i, j, k, BLUE)];
        Fy_i = Fy[INDEX(i, j, k, RED)] + Fy[INDEX(i, j, k, BLUE)];
        Fz_i = Fz[INDEX(i, j, k, RED)] + Fz[INDEX(i, j, k, BLUE)];

        u[INDEX_GLOB(i, j, k)] += 1.0 / (2.0 * rho_i) * Fx_i;
        v[INDEX_GLOB(i, j, k)] += 1.0 / (2.0 * rho_i) * Fy_i;
        w[INDEX_GLOB(i, j, k)] += 1.0 / (2.0 * rho_i) * Fz_i;
    }
}

void evaluate_density(int i, int j, int k, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    double rho_RED_i, rho_BLUE_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;

    double *rho_comp = comp_fields->rho_comp;

    double *f1 = dists->f1;

    rho_RED_i = 0.0;
    rho_BLUE_i = 0.0;
    for (int p = 0; p < NP; p++)
    {
        rho_RED_i += f1[INDEX_F(i, j, k, p, RED)];
        rho_BLUE_i += f1[INDEX_F(i, j, k, p, BLUE)];
    }
    rho_comp[INDEX(i, j, k, RED)] = rho_RED_i;
    rho_comp[INDEX(i, j, k, BLUE)] = rho_BLUE_i;
}