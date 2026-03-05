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

#if defined(LEFT_BOUNCEBACK_VELOCITY) || defined(LEFT_NEBB_VELOCITY) || defined(LEFT_BOUNCEBACK_PRESSURE) || defined(LEFT_NEBB_PRESSURE)
        if (ic < 0)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(LEFT_BOUNCEBACK_VELOCITY) || defined(LEFT_NEBB_VELOCITY)
            rho_RED_local *= (1.0 - XI_LEFT);
            rho_BLUE_local *= (1.0 + XI_LEFT);
#endif
            goto skip;
        }
#endif

#if defined(RIGHT_BOUNCEBACK_VELOCITY) || defined(RIGHT_NEBB_VELOCITY) || defined(RIGHT_BOUNCEBACK_PRESSURE) || defined(RIGHT_NEBB_PRESSURE)
        if (ic > params->NX - 1)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(RIGHT_BOUNCEBACK_VELOCITY) || defined(RIGHT_NEBB_VELOCITY)
            rho_RED_local *= (1.0 - XI_RIGHT);
            rho_BLUE_local *= (1.0 + XI_RIGHT);
#endif
            goto skip;
        }
#endif

#if defined(BOTTOM_BOUNCEBACK_VELOCITY) || defined(BOTTOM_NEBB_VELOCITY) || defined(BOTTOM_BOUNCEBACK_PRESSURE) || defined(BOTTOM_NEBB_PRESSURE)
        if (jc < 0)
        {
#ifdef THETA_C_BOTTOM
            set_wall_densities_angle_y(ic, jc, kc, 1, THETA_C_BOTTOM, &rho_RED_local, &rho_BLUE_local, sim);
#else
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
    #if defined(BOTTOM_BOUNCEBACK_VELOCITY) || defined(BOTTOM_NEBB_VELOCITY)
            rho_RED_local *= (1.0 - XI_BOTTOM);
            rho_BLUE_local *= (1.0 + XI_BOTTOM);
    #endif
#endif
            goto skip;
        }
#endif

#if defined(TOP_BOUNCEBACK_VELOCITY) || defined(TOP_NEBB_VELOCITY) || defined(TOP_BOUNCEBACK_PRESSURE) || defined(TOP_NEBB_PRESSURE)
        if (jc > NY - 1)
        {
#ifdef THETA_C_TOP
            set_wall_densities_angle_y(ic, jc, kc, -1, THETA_C_TOP, &rho_RED_local, &rho_BLUE_local, sim);
#else
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
    #if defined(TOP_BOUNCEBACK_VELOCITY) || defined(TOP_NEBB_VELOCITY)
            rho_RED_local *= (1.0 - XI_TOP);
            rho_BLUE_local *= (1.0 + XI_TOP);
    #endif
#endif
            goto skip;
        }
#endif

#if defined(BACK_BOUNCEBACK_VELOCITY) || defined(BACK_NEBB_VELOCITY) || defined(BACK_BOUNCEBACK_PRESSURE) || defined(BACK_NEBB_PRESSURE)
        if (kc < 0)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(BACK_BOUNCEBACK_VELOCITY) || defined(BACK_NEBB_VELOCITY)
            rho_RED_local *= (1.0 - BACK);
            rho_BLUE_local *= (1.0 + BACK);
#endif
            goto skip;
        }
#endif

#if defined(FRONT_BOUNCEBACK_VELOCITY) || defined(FRONT_NEBB_VELOCITY) || defined(FRONT_BOUNCEBACK_PRESSURE) || defined(FRONT_NEBB_PRESSURE)
        if (kc > NZ - 1)
        {
            extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
#if defined(FRONT_BOUNCEBACK_VELOCITY) || defined(FRONT_NEBB_VELOCITY)
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

void set_wall_densities_angle_y(int i, int j, int k, int ny, double theta_c, double *rho_RED, double *rho_BLUE, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    double grad_x_RED, grad_x_BLUE, grad_z_RED, grad_z_BLUE, grad_par_RED, grad_par_BLUE;

    int i_start = params->i_start;

    int NY = params->NY;
    int NZ = params->NZ;

    double *rho_comp = comp_fields->rho_comp;

#ifndef XPERIODIC
    if ((i < 0) || (i > params->NX - 1))
    {
        *rho_RED = 0.0;
        *rho_BLUE = 0.0;
        return;
    }
#endif

#ifndef ZPERIODIC
    if ((z < 0) || (z > NZ - 1))
    {
        *rho_RED = 0.0;
        *rho_BLUE = 0.0;
        return;
    }
#endif

#ifdef XPERIODIC
    grad_x_RED = 0.5 * (rho_comp[INDEX(i + 1, j + ny, k, RED)] - rho_comp[INDEX(i - 1, j + ny, k, RED)]);
    grad_x_BLUE = 0.5 * (rho_comp[INDEX(i + 1, j + ny, k, BLUE)] - rho_comp[INDEX(i - 1, j + ny, k, BLUE)]);
#else
    if (i == 0)
    {
        grad_x_RED = rho_comp[INDEX(i + 1, j + ny, k, RED)] - rho_comp[INDEX(i, j + ny, k, RED)];
        grad_x_BLUE = rho_comp[INDEX(i + 1, j + ny, k, BLUE)] - rho_comp[INDEX(i, j + ny, k, BLUE)];
    }
    else if (i == params->NX - 1)
    {
        grad_x_RED = rho_comp[INDEX(i, j + ny, k, RED)] - rho_comp[INDEX(i - 1, j + ny, k, RED)];
        grad_x_BLUE = rho_comp[INDEX(i, j + ny, k, BLUE)] - rho_comp[INDEX(i - 1, j + ny, k, BLUE)];
    }
    else 
    {
        grad_x_RED = 0.5 * (rho_comp[INDEX(i + 1, j + ny, k, RED)] - rho_comp[INDEX(i - 1, j + ny, k, RED)]);
        grad_x_BLUE = 0.5 * (rho_comp[INDEX(i + 1, j + ny, k, BLUE)] - rho_comp[INDEX(i - 1, j + ny, k, BLUE)]);
    }
#endif

#ifdef ZPERIODIC
    grad_z_RED = 0.5 * (rho_comp[INDEX(i, j + ny, mod(k + 1, NZ), RED)] - rho_comp[INDEX(i, j + ny, mod(k - 1, NZ), RED)]);
    grad_z_BLUE = 0.5 * (rho_comp[INDEX(i, j + ny, mod(k + 1, NZ), BLUE)] - rho_comp[INDEX(i, j + ny, mod(k - 1, NZ), BLUE)]);
#else
    if (k == 0)
    {
        grad_z_RED = rho_comp[INDEX(i, j + ny, k + 1, RED)] - rho_comp[INDEX(i, j + ny, k, RED)];
        grad_z_BLUE = rho_comp[INDEX(i, j + ny, k + 1, BLUE)] - rho_comp[INDEX(i, j + ny, k, BLUE)];
    }
    else if (k == NZ - 1)
    {
        grad_z_RED = rho_comp[INDEX(i, j + ny, k, RED)] - rho_comp[INDEX(i, j + ny, k - 1, RED)];
        grad_z_BLUE = rho_comp[INDEX(i, j + ny, k, BLUE)] - rho_comp[INDEX(i, j + ny, k - 1, BLUE)];
    }
    else 
    {
        grad_z_RED = 0.5 * (rho_comp[INDEX(i, j + ny, k + 1, RED)] - rho_comp[INDEX(i, j + ny, k - 1, RED)]);
        grad_z_BLUE = 0.5 * (rho_comp[INDEX(i, j + ny, k + 1, BLUE)] - rho_comp[INDEX(i, j + ny, k - 1, BLUE)]);
    }
#endif

    grad_par_RED = sqrt(grad_x_RED*grad_x_RED + grad_z_RED*grad_z_RED);
    grad_par_BLUE = sqrt(grad_x_BLUE*grad_x_BLUE + grad_z_BLUE*grad_z_BLUE);

    *rho_RED = 2.0*rho_comp[INDEX(i, j + 2*ny, k, RED)] - 2.0*tan(DS_PI / 2.0 - (180.0 - theta_c / 360.0 * 2.0 * DS_PI))*grad_par_RED;
    *rho_BLUE = 2.0*rho_comp[INDEX(i, j + 2*ny, k, BLUE)] + 2.0*tan(DS_PI / 2.0 - (180.0 - theta_c / 360.0 * 2.0 * DS_PI))*grad_par_BLUE;
}