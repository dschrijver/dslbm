#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/stream.h"
#include "../include/wetnode.h"

void stream_distributions(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    int ic, jc, kc;

    int p_bb;
    double rho_i, u_i, v_i, w_i, uc, u2;
    double cs2 = stencil->cs2;
    double *wp = stencil->wp;
    int *p_bounceback = stencil->p_bounceback;
    double **phi_eq = stencil->phi_eq;

    (void)glob_fields;
    (void)comp_fields;
    (void)p_bb; 
    (void)rho_i; (void)u_i; (void)v_i; (void)w_i; (void)uc; (void)u2;
    (void)cs2;
    (void)wp;
    (void)p_bounceback;
    (void)phi_eq;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    FOR_DOMAIN
    {
        for (int p = 0; p < NP; p++)
        {
            ic = i - cx[p];
            jc = j - cy[p];
            kc = k - cz[p];

#ifdef YPERIODIC
            jc = mod(jc, NY);
#endif
#ifdef ZPERIODIC
            kc = mod(kc, NZ);
#endif

#if defined(LEFT_NEBB_VELOCITY) || defined(LEFT_NEBB_PRESSURE)
            if (ic < 0)
                continue;
#endif

#if defined(RIGHT_NEBB_VELOCITY) || defined(RIGHT_NEBB_PRESSURE)
            if (ic > params->NX - 1)
                continue;
#endif

#if defined(BOTTOM_NEBB_VELOCITY) || defined(BOTTOM_NEBB_PRESSURE)
            if (jc < 0)
                continue;
#endif

#if defined(TOP_NEBB_VELOCITY) || defined(TOP_NEBB_PRESSURE)
            if (jc > NY - 1)
                continue;
#endif

#if defined(BACK_NEBB_VELOCITY) || defined(BACK_NEBB_PRESSURE)
            if (kc < 0)
                continue;
#endif

#if defined(FRONT_NEBB_VELOCITY) || defined(FRONT_NEBB_PRESSURE)
            if (kc > NZ - 1)
                continue;
#endif

#ifdef LEFT_BOUNCEBACK_VELOCITY
            if (ic < 0)
            {
                p_bb = p_bounceback[p];
                u_i = LEFT_U_VELOCITY;
                v_i = LEFT_V_VELOCITY;
                w_i = LEFT_W_VELOCITY;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], RED)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, RED)] * uc / cs2;
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], BLUE)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, BLUE)] * uc / cs2;
                continue;
            }
#endif

#ifdef RIGHT_BOUNCEBACK_VELOCITY
            if (ic > params->NX - 1)
            {
                p_bb = p_bounceback[p];
                u_i = RIGHT_U_VELOCITY;
                v_i = RIGHT_V_VELOCITY;
                w_i = RIGHT_W_VELOCITY;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], RED)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, RED)] * uc / cs2;
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], BLUE)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, BLUE)] * uc / cs2;
                continue;
            }
#endif

#ifdef BOTTOM_BOUNCEBACK_VELOCITY
            if (jc < 0)
            {
                p_bb = p_bounceback[p];
                u_i = BOTTOM_U_VELOCITY;
                v_i = BOTTOM_V_VELOCITY;
                w_i = BOTTOM_W_VELOCITY;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], RED)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, RED)] * uc / cs2;
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], BLUE)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, BLUE)] * uc / cs2;
                continue;
            }
#endif

#ifdef TOP_BOUNCEBACK_VELOCITY
            if (jc > NY - 1)
            {
                p_bb = p_bounceback[p];
                u_i = TOP_U_VELOCITY;
                v_i = TOP_V_VELOCITY;
                w_i = TOP_W_VELOCITY;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], RED)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, RED)] * uc / cs2;
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], BLUE)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, BLUE)] * uc / cs2;
                continue;
            }
#endif

#ifdef BACK_BOUNCEBACK_VELOCITY
            if (kc < 0)
            {
                p_bb = p_bounceback[p];
                u_i = BACK_U_VELOCITY;
                v_i = BACK_V_VELOCITY;
                w_i = BACK_W_VELOCITY;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], RED)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, RED)] * uc / cs2;
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], BLUE)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, BLUE)] * uc / cs2;
                continue;
            }
#endif

#ifdef FRONT_BOUNCEBACK_VELOCITY
            if (kc > NZ - 1)
            {
                p_bb = p_bounceback[p];
                u_i = FRONT_U_VELOCITY;
                v_i = FRONT_V_VELOCITY;
                w_i = FRONT_W_VELOCITY;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], RED)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, RED)] * uc / cs2;
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], BLUE)] - 2.0 * wp[p_bb] * comp_fields->rho_comp[INDEX(i, j, k, BLUE)] * uc / cs2;
                continue;
            }
#endif

#ifdef LEFT_BOUNCEBACK_PRESSURE
            if (ic < 0)
            {
                p_bb = p_bounceback[p];
                u_i = 1.5 * glob_fields->u[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->u[INDEX_GLOB(i + 1, j, k)];
                v_i = 1.5 * glob_fields->v[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->v[INDEX_GLOB(i + 1, j, k)];
                w_i = 1.5 * glob_fields->w[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->w[INDEX_GLOB(i + 1, j, k)];
                u2 = u_i * u_i + v_i * v_i + w_i * w_i;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];

                rho_i = LEFT_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                f1[INDEX_F(i, j, k, p, RED)] = -f2[INDEX_F(i, j, k, p_bb, RED)] + 2.0 * rho_i * (phi_eq[RED][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                rho_i = LEFT_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
                f1[INDEX_F(i, j, k, p, BLUE)] = -f2[INDEX_F(i, j, k, p_bb, BLUE)] + 2.0 * rho_i * (phi_eq[BLUE][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                continue;
            }
#endif

#ifdef RIGHT_BOUNCEBACK_PRESSURE
            if (ic > params->NX - 1)
            {
                p_bb = p_bounceback[p];
                u_i = 1.5 * glob_fields->u[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->u[INDEX_GLOB(i - 1, j, k)];
                v_i = 1.5 * glob_fields->v[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->v[INDEX_GLOB(i - 1, j, k)];
                w_i = 1.5 * glob_fields->w[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->w[INDEX_GLOB(i - 1, j, k)];
                u2 = u_i * u_i + v_i * v_i + w_i * w_i;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];

                rho_i = RIGHT_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                f1[INDEX_F(i, j, k, p, RED)] = -f2[INDEX_F(i, j, k, p_bb, RED)] + 2.0 * rho_i * (phi_eq[RED][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                rho_i = RIGHT_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
                f1[INDEX_F(i, j, k, p, BLUE)] = -f2[INDEX_F(i, j, k, p_bb, BLUE)] + 2.0 * rho_i * (phi_eq[BLUE][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                continue;
            }
#endif

#ifdef BOTTOM_BOUNCEBACK_PRESSURE
            if (jc < 0)
            {
                p_bb = p_bounceback[p];
                u_i = 1.5 * glob_fields->u[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->u[INDEX_GLOB(i, j + 1, k)];
                v_i = 1.5 * glob_fields->v[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->v[INDEX_GLOB(i, j + 1, k)];
                w_i = 1.5 * glob_fields->w[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->w[INDEX_GLOB(i, j + 1, k)];
                u2 = u_i * u_i + v_i * v_i + w_i * w_i;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];

                rho_i = BOTTOM_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                f1[INDEX_F(i, j, k, p, RED)] = -f2[INDEX_F(i, j, k, p_bb, RED)] + 2.0 * rho_i * (phi_eq[RED][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                rho_i = BOTTOM_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
                f1[INDEX_F(i, j, k, p, BLUE)] = -f2[INDEX_F(i, j, k, p_bb, BLUE)] + 2.0 * rho_i * (phi_eq[BLUE][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                continue;
            }
#endif

#ifdef TOP_BOUNCEBACK_PRESSURE
            if (jc > NY - 1)
            {
                p_bb = p_bounceback[p];
                u_i = 1.5 * glob_fields->u[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->u[INDEX_GLOB(i, j - 1, k)];
                v_i = 1.5 * glob_fields->v[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->v[INDEX_GLOB(i, j - 1, k)];
                w_i = 1.5 * glob_fields->w[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->w[INDEX_GLOB(i, j - 1, k)];
                u2 = u_i * u_i + v_i * v_i + w_i * w_i;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];

                rho_i = TOP_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                f1[INDEX_F(i, j, k, p, RED)] = -f2[INDEX_F(i, j, k, p_bb, RED)] + 2.0 * rho_i * (phi_eq[RED][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                rho_i = TOP_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
                f1[INDEX_F(i, j, k, p, BLUE)] = -f2[INDEX_F(i, j, k, p_bb, BLUE)] + 2.0 * rho_i * (phi_eq[BLUE][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                continue;
            }
#endif

#ifdef BACK_BOUNCEBACK_PRESSURE
            if (kc < 0)
            {
                p_bb = p_bounceback[p];
                u_i = 1.5 * glob_fields->u[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->u[INDEX_GLOB(i, j, k + 1)];
                v_i = 1.5 * glob_fields->v[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->v[INDEX_GLOB(i, j, k + 1)];
                w_i = 1.5 * glob_fields->w[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->w[INDEX_GLOB(i, j, k + 1)];
                u2 = u_i * u_i + v_i * v_i + w_i * w_i;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];

                rho_i = BACK_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                f1[INDEX_F(i, j, k, p, RED)] = -f2[INDEX_F(i, j, k, p_bb, RED)] + 2.0 * rho_i * (phi_eq[RED][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                rho_i = BACK_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
                f1[INDEX_F(i, j, k, p, BLUE)] = -f2[INDEX_F(i, j, k, p_bb, BLUE)] + 2.0 * rho_i * (phi_eq[BLUE][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                continue;
            }
#endif

#ifdef FRONT_BOUNCEBACK_PRESSURE
            if (kc > NZ - 1)
            {
                p_bb = p_bounceback[p];
                u_i = 1.5 * glob_fields->u[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->u[INDEX_GLOB(i, j, k - 1)];
                v_i = 1.5 * glob_fields->v[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->v[INDEX_GLOB(i, j, k - 1)];
                w_i = 1.5 * glob_fields->w[INDEX_GLOB(i, j, k)] - 0.5 * glob_fields->w[INDEX_GLOB(i, j, k - 1)];
                u2 = u_i * u_i + v_i * v_i + w_i * w_i;
                uc = u_i * (double)cx[p_bb] + v_i * (double)cy[p_bb] + w_i * (double)cz[p_bb];

                rho_i = FRONT_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                f1[INDEX_F(i, j, k, p, RED)] = -f2[INDEX_F(i, j, k, p_bb, RED)] + 2.0 * rho_i * (phi_eq[RED][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                rho_i = FRONT_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
                f1[INDEX_F(i, j, k, p, BLUE)] = -f2[INDEX_F(i, j, k, p_bb, BLUE)] + 2.0 * rho_i * (phi_eq[BLUE][p_bb] + wp[p_bb] * (uc * uc) / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
                continue;
            }
#endif

#ifdef YPERIODIC_FLIP
            if ((jc < 0) || (jc > NY - 1))
            {
                jc = mod(jc, NY);
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(ic, jc, kc, p, BLUE)];
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(ic, jc, kc, p, RED)];
                continue;
            }
#endif

            f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(ic, jc, kc, p, RED)];
            f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(ic, jc, kc, p, BLUE)];
        }
    }
}