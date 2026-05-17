#include <math.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/forcing.h"
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

    double rho_0_RED = params->rho_0_RED;
    double rho_0_BLUE = params->rho_0_BLUE;

    double *rho = glob_fields->rho;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;
    double *rho_N = glob_fields->rho_N;

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

        u[INDEX_GLOB(i, j, k)] = u_i / rho_i;
        v[INDEX_GLOB(i, j, k)] = v_i / rho_i;
        w[INDEX_GLOB(i, j, k)] = w_i / rho_i;

        // De Rosis 2019, 10.1063/1.5124719
        rho_N[INDEX_GLOB(i, j, k)] = (rho_RED_i / rho_0_RED - rho_BLUE_i / rho_0_BLUE) / (rho_RED_i / rho_0_RED + rho_BLUE_i / rho_0_BLUE);
    }
}

void update_final_velocity(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;

    double rho_i, Fx_i, Fy_i, Fz_i;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *rho = glob_fields->rho;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *Fx = glob_fields->Fx;
    double *Fy = glob_fields->Fy;
    double *Fz = glob_fields->Fz;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX_GLOB(i, j, k)];

        Fx_i = Fx[INDEX_GLOB(i, j, k)];
        Fy_i = Fy[INDEX_GLOB(i, j, k)];
        Fz_i = Fz[INDEX_GLOB(i, j, k)];

        u[INDEX_GLOB(i, j, k)] += 0.5 * Fx_i / rho_i;
        v[INDEX_GLOB(i, j, k)] += 0.5 * Fy_i / rho_i;
        w[INDEX_GLOB(i, j, k)] += 0.5 * Fz_i / rho_i;
    }
}

double evaluate_mass(int n, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    double M_local, M_total;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int NY = params->NY;
    int NZ = params->NZ;

    double *rho_comp = comp_fields->rho_comp;

    M_local = 0.0;
    FOR_DOMAIN
    {
        M_local += rho_comp[INDEX(i, j, k, n)];
    }

    if (params->comm_xslices == MPI_COMM_NULL)
        return 0.0;

    MPI_Allreduce(&M_local, &M_total, 1, MPI_DOUBLE, MPI_SUM, params->comm_xslices);

    return M_total;
}

void evaluate_density(int i, int j, int k, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    Stencil *stencil = sim->stencil;

    double rho_RED_i, rho_BLUE_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;

    double *rho_comp = comp_fields->rho_comp;
    double *rho = glob_fields->rho;

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
    rho[INDEX_GLOB(i, j, k)] = rho_RED_i + rho_BLUE_i;
}

void evaluate_velocity(int i, int j, int k, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    Stencil *stencil = sim->stencil;

    double rho_i, u_i, v_i, w_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *Fx = glob_fields->Fx;
    double *Fy = glob_fields->Fy;
    double *Fz = glob_fields->Fz;

    double *f1 = dists->f1;

    rho_i = 0.0;
    u_i = 0.0;
    v_i = 0.0;
    w_i = 0.0;
    for (int p = 0; p < NP; p++)
    {
        rho_i += f1[INDEX_F(i, j, k, p, RED)] + f1[INDEX_F(i, j, k, p, BLUE)];

        u_i += (f1[INDEX_F(i, j, k, p, RED)] + f1[INDEX_F(i, j, k, p, BLUE)]) * (double)cx[p];
        v_i += (f1[INDEX_F(i, j, k, p, RED)] + f1[INDEX_F(i, j, k, p, BLUE)]) * (double)cy[p];
        w_i += (f1[INDEX_F(i, j, k, p, RED)] + f1[INDEX_F(i, j, k, p, BLUE)]) * (double)cz[p];
    }

    evaluate_force(i, j, k, sim);

    u[INDEX_GLOB(i, j, k)] = (u_i + 0.5 * Fx[INDEX_GLOB(i, j, k)]) / rho_i;
    v[INDEX_GLOB(i, j, k)] = (v_i + 0.5 * Fy[INDEX_GLOB(i, j, k)]) / rho_i;
    w[INDEX_GLOB(i, j, k)] = (w_i + 0.5 * Fz[INDEX_GLOB(i, j, k)]) / rho_i;
}

void evaluate_pressure(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    double rho_RED_i, rho_BLUE_i;
    double b;

    int i_start = params->i_start;
    int i_end = params->i_end;
    int NY = params->NY;
    int NZ = params->NZ;

    double cs2_RED = params->cs2_RED;
    double cs2_BLUE = params->cs2_BLUE;

    double rho_0_RED = params->rho_0_RED;
    double rho_0_BLUE = params->rho_0_BLUE;
    double p_star_RED = rho_0_RED * cs2_RED - rho_0_BLUE * cs2_BLUE;

    double *pressure = glob_fields->pressure;
    double *rho_comp = comp_fields->rho_comp;

    FOR_DOMAIN
    {
        rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
        rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];

        b = rho_RED_i * cs2_RED + rho_BLUE_i * cs2_BLUE - p_star_RED;

        pressure[INDEX_GLOB(i, j, k)] = 0.5 * (b + sqrt(b * b + 4.0 * rho_BLUE_i * cs2_BLUE * p_star_RED));
    }
}