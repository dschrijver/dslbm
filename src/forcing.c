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

        Fx[INDEX(i, j, k, n)] = rho_i * gx;
        Fy[INDEX(i, j, k, n)] = rho_i * gy;
        Fz[INDEX(i, j, k, n)] = rho_i * gz;
    }

#ifdef SHAN_CHEN
    evaluate_shan_chen_force(i, j, k, sim);
#endif
}

void evaluate_shan_chen_force(int i, int j, int k, SimulationBag *sim)
{
#ifndef SHAN_CHEN
    (void)i;
    (void)j;
    (void)k;
    (void)sim;
#else
    ParamBag *params = sim->params;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    double rho_RED_i, rho_BLUE_i, rho_RED_local, rho_BLUE_local;
    double Fx_RED_i, Fy_RED_i, Fz_RED_i;
    double Fx_BLUE_i, Fy_BLUE_i, Fz_BLUE_i;
    int ic, jc, kc;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double G_SC = params->G_SC;

    double *rho_comp = comp_fields->rho_comp;

    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
    rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];

    Fx_RED_i = 0.0;
    Fy_RED_i = 0.0;
    Fz_RED_i = 0.0;

    Fx_BLUE_i = 0.0;
    Fy_BLUE_i = 0.0;
    Fz_BLUE_i = 0.0;

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

#if defined(LEFT_BOUNCEBACK_NOSLIP) || defined(LEFT_NEBB_NOSLIP) || defined(LEFT_BOUNCEBACK_PRESSURE) || defined(LEFT_NEBB_PRESSURE)
        if (ic < 0)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(LEFT_BOUNCEBACK_NOSLIP) || defined(LEFT_NEBB_NOSLIP)
            rho_RED_local *= (1.0 - XI_LEFT);
            rho_BLUE_local *= (1.0 + XI_LEFT);
#endif
            goto skip;
        }
#endif

#if defined(RIGHT_BOUNCEBACK_NOSLIP) || defined(RIGHT_NEBB_NOSLIP) || defined(RIGHT_BOUNCEBACK_PRESSURE) || defined(RIGHT_NEBB_PRESSURE)
        if (ic > params->NX - 1)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(RIGHT_BOUNCEBACK_NOSLIP) || defined(RIGHT_NEBB_NOSLIP)
            rho_RED_local *= (1.0 - XI_RIGHT);
            rho_BLUE_local *= (1.0 + XI_RIGHT);
#endif
            goto skip;
        }
#endif

#if defined(BOTTOM_BOUNCEBACK_NOSLIP) || defined(BOTTOM_NEBB_NOSLIP) || defined(BOTTOM_BOUNCEBACK_PRESSURE) || defined(BOTTOM_NEBB_PRESSURE)
        if (jc < 0)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(BOTTOM_BOUNCEBACK_NOSLIP) || defined(BOTTOM_NEBB_NOSLIP)
            rho_RED_local *= (1.0 - XI_BOTTOM);
            rho_BLUE_local *= (1.0 + XI_BOTTOM);
#endif
            goto skip;
        }
#endif

#if defined(TOP_BOUNCEBACK_NOSLIP) || defined(TOP_NEBB_NOSLIP) || defined(TOP_BOUNCEBACK_PRESSURE) || defined(TOP_NEBB_PRESSURE)
        if (jc > NY - 1)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(TOP_BOUNCEBACK_NOSLIP) || defined(TOP_NEBB_NOSLIP)
            rho_RED_local *= (1.0 - XI_TOP);
            rho_BLUE_local *= (1.0 + XI_TOP);
#endif
            goto skip;
        }
#endif

#if defined(BACK_BOUNCEBACK_NOSLIP) || defined(BACK_NEBB_NOSLIP) || defined(BACK_BOUNCEBACK_PRESSURE) || defined(BACK_NEBB_PRESSURE)
        if (kc < 0)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(BACK_BOUNCEBACK_NOSLIP) || defined(BACK_NEBB_NOSLIP)
            rho_RED_local *= (1.0 - BACK);
            rho_BLUE_local *= (1.0 + BACK);
#endif
            goto skip;
        }
#endif

#if defined(FRONT_BOUNCEBACK_NOSLIP) || defined(FRONT_NEBB_NOSLIP) || defined(FRONT_BOUNCEBACK_PRESSURE) || defined(FRONT_NEBB_PRESSURE)
        if (kc > NZ - 1)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(FRONT_BOUNCEBACK_NOSLIP) || defined(FRONT_NEBB_NOSLIP)
            rho_RED_local *= (1.0 - FRONT);
            rho_BLUE_local *= (1.0 + FRONT);
#endif
            goto skip;
        }
#endif

#ifdef YPERIODIC_FLIP
        if ((jc < 0) || (jc > NY - 1))
        {
            jc = mod(jc, NY);
            rho_RED_local = rho_comp[INDEX(ic, jc, kc, BLUE)];
            rho_BLUE_local = rho_comp[INDEX(ic, jc, kc, RED)];
            goto skip;
        }
#endif

        rho_RED_local = rho_comp[INDEX(ic, jc, kc, RED)];
        rho_BLUE_local = rho_comp[INDEX(ic, jc, kc, BLUE)];

    skip:

        Fx_RED_i += wp[p] * rho_BLUE_local * (double)cx[p];
        Fy_RED_i += wp[p] * rho_BLUE_local * (double)cy[p];
        Fz_RED_i += wp[p] * rho_BLUE_local * (double)cz[p];

        Fx_BLUE_i += wp[p] * rho_RED_local * (double)cx[p];
        Fy_BLUE_i += wp[p] * rho_RED_local * (double)cy[p];
        Fz_BLUE_i += wp[p] * rho_RED_local * (double)cz[p];
    }

    Fx[INDEX(i, j, k, RED)] += -G_SC * rho_RED_i * Fx_RED_i;
    Fy[INDEX(i, j, k, RED)] += -G_SC * rho_RED_i * Fy_RED_i;
    Fz[INDEX(i, j, k, RED)] += -G_SC * rho_RED_i * Fz_RED_i;

    Fx[INDEX(i, j, k, BLUE)] += -G_SC * rho_BLUE_i * Fx_BLUE_i;
    Fy[INDEX(i, j, k, BLUE)] += -G_SC * rho_BLUE_i * Fy_BLUE_i;
    Fz[INDEX(i, j, k, BLUE)] += -G_SC * rho_BLUE_i * Fz_BLUE_i;
#endif
}