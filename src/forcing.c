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
    GlobalFieldBag *glob_fields = sim->glob_fields;

    double rho_i;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;

    double gx = params->gx;
    double gy = params->gy;
    double gz = params->gz;

    double Fb_x = params->Fb_x;
    double Fb_y = params->Fb_y;
    double Fb_z = params->Fb_z;

    double *rho = glob_fields->rho;

    double *Fx = glob_fields->Fx;
    double *Fy = glob_fields->Fy;
    double *Fz = glob_fields->Fz;

    rho_i = rho[INDEX_GLOB(i, j, k)];

    Fx[INDEX_GLOB(i, j, k)] = rho_i * gx + Fb_x;
    Fy[INDEX_GLOB(i, j, k)] = rho_i * gy + Fb_y;
    Fz[INDEX_GLOB(i, j, k)] = rho_i * gz + Fb_z;

    evaluate_surface_force(i, j, k, sim);
}

void evaluate_surface_force(int i, int j, int k, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    Stencil *stencil = sim->stencil;

    double kappa;
    double n_i[3];
    double n_local[3];
    double dn[3][3];
    int ic, jc, kc;

    int i_start = params->i_start;
    int NY = params->NY;
    int NZ = params->NZ;

    double NP = stencil->NP;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;
    int *c_vec[3] = {cx, cy, cz};

    double sigma = params->sigma;

    double *Gx = glob_fields->Gx;
    double *Gy = glob_fields->Gy;
    double *Gz = glob_fields->Gz;

    double *nx = glob_fields->nx;
    double *ny = glob_fields->ny;
    double *nz = glob_fields->nz;

    double *Fx = glob_fields->Fx;
    double *Fy = glob_fields->Fy;
    double *Fz = glob_fields->Fz;

    int *flag = glob_fields->flag;

    for (int alpha = 0; alpha < 3; alpha++)
    {
        for (int beta = 0; beta < 3; beta++)
        {
            dn[alpha][beta] = 0.0;
        }
    }

    for (int p = 0; p < NP; p++)
    {
        ic = i + cx[p];
        jc = j + cy[p];
        kc = k + cz[p];

#ifdef YPERIODIC
        jc = mod(jc, NY);
#endif
#ifdef ZPERIODIC
        kc = mod(kc, NZ);
#endif

        if (flag[INDEX_FLAG(ic, jc, kc)] > 0)
        {
            n_local[0] = extrapolate_wall_n(ic, jc, kc, 0, sim);
            n_local[1] = extrapolate_wall_n(ic, jc, kc, 1, sim);
            n_local[2] = extrapolate_wall_n(ic, jc, kc, 2, sim);
            goto skip;
        }

        n_local[0] = nx[INDEX_GLOB(ic, jc, kc)];
        n_local[1] = ny[INDEX_GLOB(ic, jc, kc)];
        n_local[2] = nz[INDEX_GLOB(ic, jc, kc)];

    skip:

        for (int alpha = 0; alpha < 3; alpha++)
        {
            for (int beta = 0; beta < 3; beta++)
            {
                dn[alpha][beta] += 3.0 * wp[p] * n_local[beta] * (double)c_vec[alpha][p];
            }
        }
    }

    n_i[0] = nx[INDEX_GLOB(i, j, k)];
    n_i[1] = ny[INDEX_GLOB(i, j, k)];
    n_i[2] = nz[INDEX_GLOB(i, j, k)];

    kappa = 0.0;
    for (int alpha = 0; alpha < 3; alpha++)
    {
        kappa += -dn[alpha][alpha];
        for (int beta = 0; beta < 3; beta++)
        {
            kappa += n_i[alpha] * n_i[beta] * dn[beta][alpha];
        }
    }

    Fx[INDEX_GLOB(i, j, k)] += 0.5 * sigma * kappa * Gx[INDEX_GLOB(i, j, k)];
    Fy[INDEX_GLOB(i, j, k)] += 0.5 * sigma * kappa * Gy[INDEX_GLOB(i, j, k)];
    Fz[INDEX_GLOB(i, j, k)] += 0.5 * sigma * kappa * Gz[INDEX_GLOB(i, j, k)];
}

double extrapolate_wall_n(int i, int j, int k, int alpha, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    Stencil *stencil = sim->stencil;
    GlobalFieldBag *glob_fields = sim->glob_fields;

    int ic, jc, kc;
    double sum_n, sum_wp;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;

    int NP = stencil->NP;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double *nx = glob_fields->nx;
    double *ny = glob_fields->ny;
    double *nz = glob_fields->nz;
    double *n_vec[3] = {nx, ny, nz};

    sum_n = 0.0;
    sum_wp = 0.0;

    for (int p = 1; p < NP; p++)
    {
        ic = i + cx[p];
        jc = j + cy[p];
        kc = k + cz[p];

#ifndef XPERIODIC
        if ((ic < 0) || (ic > params->NX - 1))
            continue;
#endif
#ifndef YPERIODIC
        if ((jc < 0) || (jc > NY - 1))
            continue;
#else
        jc = mod(jc, NY);
#endif
#ifndef ZPERIODIC
        if ((kc < 0) || (kc > NZ - 1))
            continue;
#else
        kc = mod(kc, NZ);
#endif

        sum_n += wp[p] * n_vec[alpha][INDEX_GLOB(ic, jc, kc)];
        sum_wp += wp[p];
    }

    if (sum_wp == 0.0)
        return 0.0;
    else
        return sum_n / sum_wp;
}