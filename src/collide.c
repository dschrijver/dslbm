#include <math.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"

void collide_distributions_MRT(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    double rho_i, rho_RED_i, rho_BLUE_i;
    double u_i, v_i, w_i, u2, uc;
    double feq, S;
    double Fx_RED_i, Fy_RED_i, Fz_RED_i;
    double Fx_BLUE_i, Fy_BLUE_i, Fz_BLUE_i;
    double rho_N, tau;
    double f_RED_i, f_BLUE_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double tau_RED = params->tau_RED;
    double tau_BLUE = params->tau_BLUE;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double *M = stencil->M;
    double *M_inv = stencil->M_inv;
    double **omega = stencil->omega;
    double **phi_eq = stencil->phi_eq;

    double *rho = glob_fields->rho;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *rho_comp = comp_fields->rho_comp;
    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX_GLOB(i, j, k)];
        u_i = u[INDEX_GLOB(i, j, k)];
        v_i = v[INDEX_GLOB(i, j, k)];
        w_i = w[INDEX_GLOB(i, j, k)];
        u2 = u_i * u_i + v_i * v_i + w_i * w_i;

        rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
        rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];
        Fx_RED_i = Fx[INDEX(i, j, k, RED)];
        Fx_BLUE_i = Fx[INDEX(i, j, k, BLUE)];
        Fy_RED_i = Fy[INDEX(i, j, k, RED)];
        Fy_BLUE_i = Fy[INDEX(i, j, k, BLUE)];
        Fz_RED_i = Fz[INDEX(i, j, k, RED)];
        Fz_BLUE_i = Fz[INDEX(i, j, k, BLUE)];

        // Average viscosity
        rho_N = (rho_RED_i - rho_BLUE_i) / rho_i;
        tau = 0.5 * (1.0 + rho_N) * tau_RED + 0.5 * (1.0 - rho_N) * tau_BLUE;

        omega[RED][9] = 1.0 / tau;
        omega[BLUE][9] = 1.0 / tau;
        omega[RED][11] = 1.0 / tau;
        omega[BLUE][11] = 1.0 / tau;
        omega[RED][13] = 1.0 / tau;
        omega[BLUE][13] = 1.0 / tau;
        omega[RED][14] = 1.0 / tau;
        omega[BLUE][14] = 1.0 / tau;
        omega[RED][15] = 1.0 / tau;
        omega[BLUE][15] = 1.0 / tau;

        // Collision in moment space
        for (int p = 0; p < NP; p++)
        {
            f2[INDEX_F(i, j, k, p, RED)] = 0.0;
            f2[INDEX_F(i, j, k, p, BLUE)] = 0.0;

            for (int pk = 0; pk < NP; pk++)
            {
                uc = u_i * (double)cx[pk] + v_i * (double)cy[pk] + w_i * (double)cz[pk];

                feq = wp[pk] * (uc / cs2 + (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                S = wp[pk] * ((((double)cx[pk] - u_i) / cs2 + uc / (cs2 * cs2) * (double)cx[pk]) * Fx_RED_i + (((double)cy[pk] - v_i) / cs2 + uc / (cs2 * cs2) * (double)cy[pk]) * Fy_RED_i + (((double)cz[pk] - w_i) / cs2 + uc / (cs2 * cs2) * (double)cz[pk]) * Fz_RED_i);

                f2[INDEX_F(i, j, k, p, RED)] += -omega[RED][p] * M[INDEX_MRT(p, pk)] * (f1[INDEX_F(i, j, k, pk, RED)] - rho_RED_i * (feq + phi_eq[RED][pk])) + (1.0 - 0.5 * omega[RED][p]) * M[INDEX_MRT(p, pk)] * S;

                S = wp[pk] * ((((double)cx[pk] - u_i) / cs2 + uc / (cs2 * cs2) * (double)cx[pk]) * Fx_BLUE_i + (((double)cy[pk] - v_i) / cs2 + uc / (cs2 * cs2) * (double)cy[pk]) * Fy_BLUE_i + (((double)cz[pk] - w_i) / cs2 + uc / (cs2 * cs2) * (double)cz[pk]) * Fz_BLUE_i);

                f2[INDEX_F(i, j, k, p, BLUE)] += -omega[BLUE][p] * M[INDEX_MRT(p, pk)] * (f1[INDEX_F(i, j, k, pk, BLUE)] - rho_BLUE_i * (feq + phi_eq[BLUE][pk])) + (1.0 - 0.5 * omega[BLUE][p]) * M[INDEX_MRT(p, pk)] * S;
            }
        }

        // Transform back to population space
        for (int p = 0; p < NP; p++)
        {
            f_RED_i = 0.0;
            f_BLUE_i = 0.0;
            for (int pk = 0; pk < NP; pk++)
            {
                f_RED_i += M_inv[INDEX_MRT(p, pk)] * f2[INDEX_F(i, j, k, pk, RED)];
                f_BLUE_i += M_inv[INDEX_MRT(p, pk)] * f2[INDEX_F(i, j, k, pk, BLUE)];
            }
            f1[INDEX_F(i, j, k, p, RED)] += f_RED_i;
            f1[INDEX_F(i, j, k, p, BLUE)] += f_BLUE_i;
        }

        // Swap populations
        for (int p = 0; p < NP; p++)
        {
            f2[INDEX_F(i, j, k, p, RED)] = f1[INDEX_F(i, j, k, p, RED)];
            f2[INDEX_F(i, j, k, p, BLUE)] = f1[INDEX_F(i, j, k, p, BLUE)];
        }
    }
}

void collide_distributions_CGM(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    double Gx, Gy, Gz, G, Gc;
    double rho_i, rho_RED_i, rho_BLUE_i;
    double rho_RED_local, rho_BLUE_local, rho_local;
    double rho_N, tau, A, fstar, feq, cos_phi;
    int ic, jc, kc;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double *B = stencil->B;
    double **phi_eq = stencil->phi_eq;

    double tau_RED = params->tau_RED;
    double tau_BLUE = params->tau_BLUE;
    double sigma = params->sigma;
    double beta = params->beta;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *rho = glob_fields->rho;
    double *rho_comp = comp_fields->rho_comp;

    double *f2 = dists->f2;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX_GLOB(i, j, k)];
        rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
        rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];

        // Compute Color Gradient
        Gx = 0.0;
        Gy = 0.0;
        Gz = 0.0;
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

#ifdef LEFT_BOUNCEBACK
            if (ic < 0)
            {
                extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
                goto skip;
            }
#endif

#ifdef RIGHT_BOUNCEBACK
            if (ic > params->NX - 1)
            {
                extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
                goto skip;
            }
#endif

#ifdef BOTTOM_BOUNCEBACK
            if (jc < 0)
            {
                extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
                goto skip;
            }
#endif

#ifdef TOP_BOUNCEBACK
            if (jc > NY - 1)
            {
                extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
                goto skip;
            }
#endif

#ifdef BACK_BOUNCEBACK
            if (kc < 0)
            {
                extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
                goto skip;
            }
#endif

#ifdef FRONT_BOUNCEBACK
            if (kc > NZ - 1)
            {
                extrapolate_wall_density(ic, jc, kc, &rho_RED_local, &rho_BLUE_local, sim);
                goto skip;
            }
#endif

            rho_RED_local = rho_comp[INDEX(ic, jc, kc, RED)];
            rho_BLUE_local = rho_comp[INDEX(ic, jc, kc, BLUE)];

        skip:

            rho_local = rho_RED_local + rho_BLUE_local;

            Gx += wp[p] * (rho_RED_local - rho_BLUE_local) / rho_local * (double)cx[p];
            Gy += wp[p] * (rho_RED_local - rho_BLUE_local) / rho_local * (double)cy[p];
            Gz += wp[p] * (rho_RED_local - rho_BLUE_local) / rho_local * (double)cz[p];
        }
        Gx /= cs2;
        Gy /= cs2;
        Gz /= cs2;

#ifdef LEFT_BOUNCEBACK
        if (i == 0)
        {
            Gx = -tan(DS_PI / 2.0 - (180.0 - THETA_C_LEFT) / 360.0 * 2.0 * DS_PI) * sqrt(Gy * Gy + Gz * Gz);
        }
#endif

#ifdef RIGHT_BOUNCEBACK
        if (i == params->NX - 1)
        {
            Gx = tan(DS_PI / 2.0 - (180.0 - THETA_C_RIGHT) / 360.0 * 2.0 * DS_PI) * sqrt(Gy * Gy + Gz * Gz);
        }
#endif

#ifdef BOTTOM_BOUNCEBACK
        if (j == 0)
        {
            Gy = -tan(DS_PI / 2.0 - (180.0 - THETA_C_BOTTOM) / 360.0 * 2.0 * DS_PI) * sqrt(Gx * Gx + Gz * Gz);
        }
#endif

#ifdef TOP_BOUNCEBACK
        if (j == NY - 1)
        {
            Gy = tan(DS_PI / 2.0 - (180.0 - THETA_C_TOP) / 360.0 * 2.0 * DS_PI) * sqrt(Gx * Gx + Gz * Gz);
        }
#endif

#ifdef BACK_BOUNCEBACK
        if (k == 0)
        {
            Gz = -tan(DS_PI / 2.0 - (180.0 - THETA_C_LEFT) / 360.0 * 2.0 * DS_PI) * sqrt(Gx * Gx + Gy * Gy);
        }
#endif

#ifdef FRONT_BOUNCEBACK
        if (k == NZ - 1)
        {
            Gy = tan(DS_PI / 2.0 - (180.0 - THETA_C_RIGHT) / 360.0 * 2.0 * DS_PI) * sqrt(Gx * Gx + Gy * Gy);
        }
#endif

        G = sqrt(Gx * Gx + Gy * Gy + Gz * Gz);

        if (G < 1e-15)
            continue;

        // Average viscosity
        rho_N = (rho_RED_i - rho_BLUE_i) / rho_i;
        tau = 0.5 * (1.0 + rho_N) * tau_RED + 0.5 * (1.0 - rho_N) * tau_BLUE;

        // Perturbation
        A = 9.0 / 8.0 * sigma / tau;
        for (int p = 0; p < NP; p++)
        {
            Gc = Gx * (double)cx[p] + Gy * (double)cy[p] + Gz * (double)cz[p];
            f2[INDEX_F(i, j, k, p, RED)] += A * G * (wp[p] * Gc * Gc / (G * G) - B[p]);
            f2[INDEX_F(i, j, k, p, BLUE)] += A * G * (wp[p] * Gc * Gc / (G * G) - B[p]);
        }

        // Recoloring
        for (int p = 0; p < NP; p++)
        {
            fstar = f2[INDEX_F(i, j, k, p, RED)] + f2[INDEX_F(i, j, k, p, BLUE)];
            f2[INDEX_F(i, j, k, p, RED)] = rho_RED_i / rho_i * fstar;
            f2[INDEX_F(i, j, k, p, BLUE)] = rho_BLUE_i / rho_i * fstar;
        }
        for (int p = 1; p < NP; p++)
        {
            feq = rho_RED_i * phi_eq[RED][p] + rho_BLUE_i * phi_eq[BLUE][p];
            Gc = Gx * (double)cx[p] + Gy * (double)cy[p] + Gz * (double)cz[p];
            cos_phi = Gc / (sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]) * G);
            f2[INDEX_F(i, j, k, p, RED)] += beta * rho_RED_i * rho_BLUE_i / (rho_i * rho_i) * cos_phi * feq;
            f2[INDEX_F(i, j, k, p, BLUE)] -= beta * rho_RED_i * rho_BLUE_i / (rho_i * rho_i) * cos_phi * feq;
        }
    }
}

void extrapolate_wall_density(int i, int j, int k, double *rho_RED, double *rho_BLUE, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    Stencil *stencil = sim->stencil;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    int ic, jc, kc;
    double sum_rho_RED, sum_rho_BLUE, sum_wp;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;

    int NP = stencil->NP;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double *rho_comp = comp_fields->rho_comp;

    sum_rho_RED = 0.0;
    sum_rho_BLUE = 0.0;
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

        sum_rho_RED += wp[p] * rho_comp[INDEX(ic, jc, kc, RED)];
        sum_rho_BLUE += wp[p] * rho_comp[INDEX(ic, jc, kc, BLUE)];
        sum_wp += wp[p];
    }

    *rho_RED = sum_rho_RED / sum_wp;
    *rho_BLUE = sum_rho_BLUE / sum_wp;
}