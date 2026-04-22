#include <math.h>
#include <string.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"

void compute_equilibrium(double rho, double u, double v, double w, double cs2, double *feq, SimulationBag *sim)
{
    DistributionBag *dists = sim->dists;
    double *raw = dists->raw;
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;
    double cs4 = cs2 * cs2;
    double cs6 = cs2 * cs4;

    raw[0] = rho;
    raw[1] = rho * u;
    raw[2] = rho * v;
    raw[3] = rho * w;
    raw[4] = rho * u * v;
    raw[5] = rho * u * w;
    raw[6] = rho * v * w;
    raw[7] = rho * (u2 - v2);
    raw[8] = rho * (u2 - w2);
    raw[9] = rho * (3.0 * cs2 + u2 + v2 + w2);
    raw[10] = rho * u * (2.0 * cs2 + v2 + w2);
    raw[11] = rho * v * (2.0 * cs2 + u2 + w2);
    raw[12] = rho * w * (2.0 * cs2 + u2 + v2);
    raw[13] = rho * u * (v2 - w2);
    raw[14] = rho * v * (u2 - w2);
    raw[15] = rho * w * (u2 - v2);
    raw[16] = rho * u * v * w;
    raw[17] = rho * (2.0 * cs2 * (u2 + v2 + w2) + cs2 + u2 * v2 + u2 * w2 + v2 * w2);
    raw[18] = rho * (2.0 * cs2 * u2 + cs4 + u2 * v2 + u2 * w2 - v2 * w2);
    raw[19] = rho * (cs2 + u2) * (v2 - w2);
    raw[20] = rho * v * w * (cs2 + u2);
    raw[21] = rho * u * w * (cs2 + v2);
    raw[22] = rho * u * v * (cs2 + w2);
    raw[23] = rho * u * (2.0 * cs2 * (v2 + w2) + cs2 - cs4 + 2.0 * v2 * w2) / 2.0;
    raw[24] = rho * v * (4.0 * cs2 * (u2 + w2) + cs2 + cs4 + 4.0 * u2 * w2) / 4.0;
    raw[25] = rho * w * (4.0 * cs2 * (u2 + v2) + cs2 + cs4 + 4.0 * u2 * v2) / 4.0;
    raw[26] = rho * (cs2 * (2.0 * u2 + v2 + w2) + 4.0 * cs2 * (u2 * v2 + u2 * w2 + v2 * w2) + cs4 * (-2.0 * u2 + v2 + w2) + 4.0 * cs6 + 4.0 * u2 * v2 * w2) / 4.0;

    feq[0] = raw[0] + raw[17] - raw[26] - raw[9];
    feq[1] = -raw[10] / 2.0 - raw[17] / 4.0 - raw[18] / 4.0 + raw[1] / 2.0 + raw[23] / 2.0 + raw[26] / 2.0 + raw[7] / 6.0 + raw[8] / 6.0 + raw[9] / 6.0;
    feq[2] = raw[10] / 2.0 - raw[17] / 4.0 - raw[18] / 4.0 - raw[1] / 2.0 - raw[23] / 2.0 + raw[26] / 2.0 + raw[7] / 6.0 + raw[8] / 6.0 + raw[9] / 6.0;
    feq[3] = -raw[11] / 2.0 - 3.0 * raw[17] / 8.0 + raw[18] / 8.0 - raw[19] / 4.0 + raw[24] / 2.0 + raw[26] / 2.0 + raw[2] / 2.0 - raw[7] / 3.0 + raw[8] / 6.0 + raw[9] / 6.0;
    feq[4] = raw[11] / 2.0 - 3.0 * raw[17] / 8.0 + raw[18] / 8.0 - raw[19] / 4.0 - raw[24] / 2.0 + raw[26] / 2.0 - raw[2] / 2.0 - raw[7] / 3.0 + raw[8] / 6.0 + raw[9] / 6.0;
    feq[5] = -raw[12] / 2.0 - 3.0 * raw[17] / 8.0 + raw[18] / 8.0 + raw[19] / 4.0 + raw[25] / 2.0 + raw[26] / 2.0 + raw[3] / 2.0 + raw[7] / 6.0 - raw[8] / 3.0 + raw[9] / 6.0;
    feq[6] = raw[12] / 2.0 - 3.0 * raw[17] / 8.0 + raw[18] / 8.0 + raw[19] / 4.0 - raw[25] / 2.0 + raw[26] / 2.0 - raw[3] / 2.0 + raw[7] / 6.0 - raw[8] / 3.0 + raw[9] / 6.0;
    feq[7] = raw[10] / 8.0 + raw[11] / 8.0 + raw[13] / 8.0 + raw[14] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 + raw[19] / 8.0 - raw[22] / 4.0 - raw[23] / 4.0 - raw[24] / 4.0 - raw[26] / 4.0 + raw[4] / 4.0;
    feq[8] = -raw[10] / 8.0 + raw[11] / 8.0 - raw[13] / 8.0 + raw[14] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 + raw[19] / 8.0 + raw[22] / 4.0 + raw[23] / 4.0 - raw[24] / 4.0 - raw[26] / 4.0 - raw[4] / 4.0;
    feq[9] = raw[10] / 8.0 - raw[11] / 8.0 + raw[13] / 8.0 - raw[14] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 + raw[19] / 8.0 + raw[22] / 4.0 - raw[23] / 4.0 + raw[24] / 4.0 - raw[26] / 4.0 - raw[4] / 4.0;
    feq[10] = -raw[10] / 8.0 - raw[11] / 8.0 - raw[13] / 8.0 - raw[14] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 + raw[19] / 8.0 - raw[22] / 4.0 + raw[23] / 4.0 + raw[24] / 4.0 - raw[26] / 4.0 + raw[4] / 4.0;
    feq[11] = raw[10] / 8.0 + raw[12] / 8.0 - raw[13] / 8.0 + raw[15] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 - raw[19] / 8.0 - raw[21] / 4.0 - raw[23] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 + raw[5] / 4.0;
    feq[12] = -raw[10] / 8.0 + raw[12] / 8.0 + raw[13] / 8.0 + raw[15] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 - raw[19] / 8.0 + raw[21] / 4.0 + raw[23] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 - raw[5] / 4.0;
    feq[13] = raw[10] / 8.0 - raw[12] / 8.0 - raw[13] / 8.0 - raw[15] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 - raw[19] / 8.0 + raw[21] / 4.0 - raw[23] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 - raw[5] / 4.0;
    feq[14] = -raw[10] / 8.0 - raw[12] / 8.0 + raw[13] / 8.0 - raw[15] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 - raw[19] / 8.0 - raw[21] / 4.0 + raw[23] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 + raw[5] / 4.0;
    feq[15] = raw[11] / 8.0 + raw[12] / 8.0 - raw[14] / 8.0 - raw[15] / 8.0 + raw[17] / 8.0 - raw[18] / 8.0 - raw[20] / 4.0 - raw[24] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 + raw[6] / 4.0;
    feq[16] = -raw[11] / 8.0 + raw[12] / 8.0 + raw[14] / 8.0 - raw[15] / 8.0 + raw[17] / 8.0 - raw[18] / 8.0 + raw[20] / 4.0 + raw[24] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 - raw[6] / 4.0;
    feq[17] = raw[11] / 8.0 - raw[12] / 8.0 - raw[14] / 8.0 + raw[15] / 8.0 + raw[17] / 8.0 - raw[18] / 8.0 + raw[20] / 4.0 - raw[24] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 - raw[6] / 4.0;
    feq[18] = -raw[11] / 8.0 - raw[12] / 8.0 + raw[14] / 8.0 + raw[15] / 8.0 + raw[17] / 8.0 - raw[18] / 8.0 - raw[20] / 4.0 + raw[24] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 + raw[6] / 4.0;
    feq[19] = raw[16] / 8.0 + raw[20] / 8.0 + raw[21] / 8.0 + raw[22] / 8.0 + raw[23] / 8.0 + raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
    feq[20] = -raw[16] / 8.0 + raw[20] / 8.0 - raw[21] / 8.0 - raw[22] / 8.0 - raw[23] / 8.0 + raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
    feq[21] = -raw[16] / 8.0 - raw[20] / 8.0 + raw[21] / 8.0 - raw[22] / 8.0 + raw[23] / 8.0 - raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
    feq[22] = raw[16] / 8.0 - raw[20] / 8.0 - raw[21] / 8.0 + raw[22] / 8.0 - raw[23] / 8.0 - raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
    feq[23] = -raw[16] / 8.0 - raw[20] / 8.0 - raw[21] / 8.0 + raw[22] / 8.0 + raw[23] / 8.0 + raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
    feq[24] = raw[16] / 8.0 - raw[20] / 8.0 + raw[21] / 8.0 - raw[22] / 8.0 - raw[23] / 8.0 + raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
    feq[25] = raw[16] / 8.0 + raw[20] / 8.0 - raw[21] / 8.0 - raw[22] / 8.0 + raw[23] / 8.0 - raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
    feq[26] = -raw[16] / 8.0 + raw[20] / 8.0 + raw[21] / 8.0 + raw[22] / 8.0 - raw[23] / 8.0 - raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
}

void evaluate_color_gradients(SimulationBag *sim)
{

    ParamBag *params = sim->params;
    Stencil *stencil = sim->stencil;
    GlobalFieldBag *glob_fields = sim->glob_fields;

    double Gx_i, Gy_i, Gz_i, G_i;
    double nhat_x, nhat_y, nhat_z;
    double rho_N_local;
    int ic, jc, kc;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

#ifdef THETA_C
    double nhat_norm, Gn;
    double theta = THETA_C / 360.0 * 2.0 * DS_PI;
#endif

    int *flag = glob_fields->flag;

    double *rho_N = glob_fields->rho_N;
    double *G_norm = glob_fields->G_norm;
    double *Gx = glob_fields->Gx;
    double *Gy = glob_fields->Gy;
    double *Gz = glob_fields->Gz;

    FOR_DOMAIN
    {
        // Compute Color Gradient
        Gx_i = 0.0;
        Gy_i = 0.0;
        Gz_i = 0.0;
        nhat_x = 0.0;
        nhat_y = 0.0;
        nhat_z = 0.0;
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
                rho_N_local = extrapolate_wall_rho_N(ic, jc, kc, sim);
                nhat_x += -wp[p] * (double)cx[p];
                nhat_y += -wp[p] * (double)cy[p];
                nhat_z += -wp[p] * (double)cz[p];
                goto skip;
            }

            rho_N_local = rho_N[INDEX_GLOB(ic, jc, kc)];

        skip:

            Gx_i += wp[p] * rho_N_local * (double)cx[p];
            Gy_i += wp[p] * rho_N_local * (double)cy[p];
            Gz_i += wp[p] * rho_N_local * (double)cz[p];
        }
        Gx_i *= 3.0;
        Gy_i *= 3.0;
        Gz_i *= 3.0;

#ifdef THETA_C
        nhat_norm = sqrt(nhat_x * nhat_x + nhat_y * nhat_y + nhat_z * nhat_z);

        if (nhat_norm > 0.0)
        {
            nhat_x /= nhat_norm;
            nhat_y /= nhat_norm;
            nhat_z /= nhat_norm;

            Gn = nhat_x * Gx_i + nhat_y * Gy_i + nhat_z * Gz_i;
            Gx_i -= Gn * nhat_x;
            Gy_i -= Gn * nhat_y;
            Gz_i -= Gn * nhat_z;

            Gn = tan(theta - 0.5 * DS_PI) * sqrt(Gx_i * Gx_i + Gy_i * Gy_i + Gz_i * Gz_i);

            Gx_i += Gn * nhat_x;
            Gy_i += Gn * nhat_y;
            Gz_i += Gn * nhat_z;
        }
#endif

        G_i = sqrt(Gx_i * Gx_i + Gy_i * Gy_i + Gz_i * Gz_i);

        G_norm[INDEX_GLOB(i, j, k)] = G_i;
        Gx[INDEX_GLOB(i, j, k)] = Gx_i;
        Gy[INDEX_GLOB(i, j, k)] = Gy_i;
        Gz[INDEX_GLOB(i, j, k)] = Gz_i;
    }
}

double extrapolate_wall_rho_N(int i, int j, int k, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    Stencil *stencil = sim->stencil;
    GlobalFieldBag *glob_fields = sim->glob_fields;

    int ic, jc, kc;
    double sum_rho_N, sum_wp;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;

    int NP = stencil->NP;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double *rho_N = glob_fields->rho_N;

    sum_rho_N = 0.0;
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

        sum_rho_N += wp[p] * rho_N[INDEX_GLOB(ic, jc, kc)];
        sum_wp += wp[p];
    }

    if (sum_wp == 0.0)
        return 0.0;
    else
        return sum_rho_N / sum_wp;
}

void collide(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    double G, G2, Gx_i, Gy_i, Gz_i, G2x_i, G2y_i, G2z_i;
    double rho_i;
    double u_i, v_i, w_i;
    double u2_i, v2_i, w2_i;
    double Fx_i, Fy_i, Fz_i;
    double k4, k5, k6, k7, k8;
    double cx_bar, cy_bar, cz_bar;
    double cs2_i, cs4_i, cs6_i;
    double rho_N_i, omega;
    double f_star, rho_RED_i, rho_BLUE_i;
    double Gc, cos_phi, p_rec;

    int i_start = params->i_start;
    int i_end = params->i_end;
    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *cs2 = stencil->cs2;
    double *tau = stencil->tau;

    double sigma = params->sigma;
    double beta = params->beta;

    // Wen 2019, 10.1103/PhysRevE.100.023301
    double delta_tau = 0.98;

    // Grunau 1993, 10.1063/1.858769
    double alpha_tau = 2.0 * tau[RED] * tau[BLUE] / (tau[RED] + tau[BLUE]);
    double beta_tau = 2.0 * (tau[RED] - alpha_tau) / delta_tau;
    double kappa_tau = -beta_tau / (2.0 * delta_tau);
    double eta_tau = 2.0 * (alpha_tau - tau[BLUE]) / delta_tau;
    double xi_tau = eta_tau / (2.0 * delta_tau);

    double *G_norm = glob_fields->G_norm;
    double *Gx = glob_fields->Gx;
    double *Gy = glob_fields->Gy;
    double *Gz = glob_fields->Gz;

    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *rho_N = glob_fields->rho_N;

    double *rho_comp = comp_fields->rho_comp;
    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    double *k_pert = dists->k_pert;
    double *k_star = dists->k_star;
    double *raw = dists->raw;
    double *feq = dists->feq;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    FOR_DOMAIN
    {
        u_i = u[INDEX_GLOB(i, j, k)];
        v_i = v[INDEX_GLOB(i, j, k)];
        w_i = w[INDEX_GLOB(i, j, k)];

        u2_i = u_i * u_i;
        v2_i = v_i * v_i;
        w2_i = w_i * w_i;

        G = G_norm[INDEX_GLOB(i, j, k)];
        Gx_i = Gx[INDEX_GLOB(i, j, k)];
        Gy_i = Gy[INDEX_GLOB(i, j, k)];
        Gz_i = Gz[INDEX_GLOB(i, j, k)];

        G2 = G * G;
        G2x_i = Gx_i * Gx_i;
        G2y_i = Gy_i * Gy_i;
        G2z_i = Gz_i * Gz_i;

        rho_N_i = rho_N[INDEX_GLOB(i, j, k)];
        if (rho_N_i > delta_tau)
        {
            omega = 1.0 / tau[RED];
        }
        else if (rho_N_i > 0)
        {
            omega = 1.0 / (alpha_tau + beta_tau * rho_N_i + kappa_tau * rho_N_i * rho_N_i);
        }
        else if (rho_N_i >= -delta_tau)
        {
            omega = 1.0 / (alpha_tau + eta_tau * rho_N_i + xi_tau * rho_N_i * rho_N_i);
        }
        else
        {
            omega = 1.0 / tau[BLUE];
        }

        if (G < 1e-15)
        {
            memset(k_pert, 0, NP * sizeof(double));
        }
        else
        {
            k_pert[0] = -3.0 * sigma * (G2 - G2x_i - G2y_i - G2z_i) / (8.0 * G);
            k_pert[1] = 3.0 * sigma * u_i * (G2 - G2x_i - G2y_i - G2z_i) / (8.0 * G);
            k_pert[2] = 3.0 * sigma * v_i * (G2 - G2x_i - G2y_i - G2z_i) / (8.0 * G);
            k_pert[3] = 3.0 * sigma * w_i * (G2 - G2x_i - G2y_i - G2z_i) / (8.0 * G);
            k_pert[4] = -sigma * (3.0 * G2 * u_i * v_i - 3.0 * G2x_i * u_i * v_i - 3.0 * G2y_i * u_i * v_i - 3.0 * G2z_i * u_i * v_i - 2.0 * Gx_i * Gy_i) / (8.0 * G);
            k_pert[5] = -sigma * (3.0 * G2 * u_i * w_i - 3.0 * G2x_i * u_i * w_i - 3.0 * G2y_i * u_i * w_i - 3.0 * G2z_i * u_i * w_i - 2.0 * Gx_i * Gz_i) / (8.0 * G);
            k_pert[6] = -sigma * (3.0 * G2 * v_i * w_i - 3.0 * G2x_i * v_i * w_i - 3.0 * G2y_i * v_i * w_i - 3.0 * G2z_i * v_i * w_i - 2.0 * Gy_i * Gz_i) / (8.0 * G);
            k_pert[7] = -sigma * (3.0 * G2 * u2_i - 3.0 * G2 * v2_i - 3.0 * G2x_i * u2_i + 3.0 * G2x_i * v2_i - 2.0 * G2x_i - 3.0 * G2y_i * u2_i + 3.0 * G2y_i * v2_i + 2.0 * G2y_i - 3.0 * G2z_i * u2_i + 3.0 * G2z_i * v2_i) / (8.0 * G);
            k_pert[8] = -sigma * (3.0 * G2 * u2_i - 3.0 * G2 * w2_i - 3.0 * G2x_i * u2_i + 3.0 * G2x_i * w2_i - 2.0 * G2x_i - 3.0 * G2y_i * u2_i + 3.0 * G2y_i * w2_i - 3.0 * G2z_i * u2_i + 3.0 * G2z_i * w2_i + 2.0 * G2z_i) / (8.0 * G);
            k_pert[9] = -sigma * (3.0 * G2 * u2_i + 3.0 * G2 * v2_i + 3.0 * G2 * w2_i + 9.0 * G2 - 3.0 * G2x_i * u2_i - 3.0 * G2x_i * v2_i - 3.0 * G2x_i * w2_i - 5.0 * G2x_i - 3.0 * G2y_i * u2_i - 3.0 * G2y_i * v2_i - 3.0 * G2y_i * w2_i - 5.0 * G2y_i - 3.0 * G2z_i * u2_i - 3.0 * G2z_i * v2_i - 3.0 * G2z_i * w2_i - 5.0 * G2z_i) / (8.0 * G);
            k_pert[10] = sigma * (3.0 * G2 * u_i * v2_i + 3.0 * G2 * u_i * w2_i + 6.0 * G2 * u_i - 3.0 * G2x_i * u_i * v2_i - 3.0 * G2x_i * u_i * w2_i - 2.0 * G2x_i * u_i - 3.0 * G2y_i * u_i * v2_i - 3.0 * G2y_i * u_i * w2_i - 4.0 * G2y_i * u_i - 3.0 * G2z_i * u_i * v2_i - 3.0 * G2z_i * u_i * w2_i - 4.0 * G2z_i * u_i - 4.0 * Gx_i * Gy_i * v_i - 4.0 * Gx_i * Gz_i * w_i) / (8.0 * G);
            k_pert[11] = sigma * (3.0 * G2 * u2_i * v_i + 3.0 * G2 * v_i * w2_i + 6.0 * G2 * v_i - 3.0 * G2x_i * u2_i * v_i - 3.0 * G2x_i * v_i * w2_i - 4.0 * G2x_i * v_i - 3.0 * G2y_i * u2_i * v_i - 3.0 * G2y_i * v_i * w2_i - 2.0 * G2y_i * v_i - 3.0 * G2z_i * u2_i * v_i - 3.0 * G2z_i * v_i * w2_i - 4.0 * G2z_i * v_i - 4.0 * Gx_i * Gy_i * u_i - 4.0 * Gy_i * Gz_i * w_i) / (8.0 * G);
            k_pert[12] = sigma * (3.0 * G2 * u2_i * w_i + 3.0 * G2 * v2_i * w_i + 6.0 * G2 * w_i - 3.0 * G2x_i * u2_i * w_i - 3.0 * G2x_i * v2_i * w_i - 4.0 * G2x_i * w_i - 3.0 * G2y_i * u2_i * w_i - 3.0 * G2y_i * v2_i * w_i - 4.0 * G2y_i * w_i - 3.0 * G2z_i * u2_i * w_i - 3.0 * G2z_i * v2_i * w_i - 2.0 * G2z_i * w_i - 4.0 * Gx_i * Gz_i * u_i - 4.0 * Gy_i * Gz_i * v_i) / (8.0 * G);
            k_pert[13] = sigma * (3.0 * G2 * u_i * v2_i - 3.0 * G2 * u_i * w2_i - 3.0 * G2x_i * u_i * v2_i + 3.0 * G2x_i * u_i * w2_i - 3.0 * G2y_i * u_i * v2_i + 3.0 * G2y_i * u_i * w2_i - 2.0 * G2y_i * u_i - 3.0 * G2z_i * u_i * v2_i + 3.0 * G2z_i * u_i * w2_i + 2.0 * G2z_i * u_i - 4.0 * Gx_i * Gy_i * v_i + 4.0 * Gx_i * Gz_i * w_i) / (8.0 * G);
            k_pert[14] = sigma * (3.0 * G2 * u2_i * v_i - 3.0 * G2 * v_i * w2_i - 3.0 * G2x_i * u2_i * v_i + 3.0 * G2x_i * v_i * w2_i - 2.0 * G2x_i * v_i - 3.0 * G2y_i * u2_i * v_i + 3.0 * G2y_i * v_i * w2_i - 3.0 * G2z_i * u2_i * v_i + 3.0 * G2z_i * v_i * w2_i + 2.0 * G2z_i * v_i - 4.0 * Gx_i * Gy_i * u_i + 4.0 * Gy_i * Gz_i * w_i) / (8.0 * G);
            k_pert[15] = sigma * (3.0 * G2 * u2_i * w_i - 3.0 * G2 * v2_i * w_i - 3.0 * G2x_i * u2_i * w_i + 3.0 * G2x_i * v2_i * w_i - 2.0 * G2x_i * w_i - 3.0 * G2y_i * u2_i * w_i + 3.0 * G2y_i * v2_i * w_i + 2.0 * G2y_i * w_i - 3.0 * G2z_i * u2_i * w_i + 3.0 * G2z_i * v2_i * w_i - 4.0 * Gx_i * Gz_i * u_i + 4.0 * Gy_i * Gz_i * v_i) / (8.0 * G);
            k_pert[16] = sigma * (3.0 * G2 * u_i * v_i * w_i - 3.0 * G2x_i * u_i * v_i * w_i - 3.0 * G2y_i * u_i * v_i * w_i - 3.0 * G2z_i * u_i * v_i * w_i - 2.0 * Gx_i * Gy_i * w_i - 2.0 * Gx_i * Gz_i * v_i - 2.0 * Gy_i * Gz_i * u_i) / (8.0 * G);
            k_pert[17] = -sigma * (9.0 * G2 * u2_i * v2_i + 9.0 * G2 * u2_i * w2_i + 18.0 * G2 * u2_i + 9.0 * G2 * v2_i * w2_i + 18.0 * G2 * v2_i + 18.0 * G2 * w2_i + 15.0 * G2 - 9.0 * G2x_i * u2_i * v2_i - 9.0 * G2x_i * u2_i * w2_i - 6.0 * G2x_i * u2_i - 9.0 * G2x_i * v2_i * w2_i - 12.0 * G2x_i * v2_i - 12.0 * G2x_i * w2_i - 7.0 * G2x_i - 9.0 * G2y_i * u2_i * v2_i - 9.0 * G2y_i * u2_i * w2_i - 12.0 * G2y_i * u2_i - 9.0 * G2y_i * v2_i * w2_i - 6.0 * G2y_i * v2_i - 12.0 * G2y_i * w2_i - 7.0 * G2y_i - 9.0 * G2z_i * u2_i * v2_i - 9.0 * G2z_i * u2_i * w2_i - 12.0 * G2z_i * u2_i - 9.0 * G2z_i * v2_i * w2_i - 12.0 * G2z_i * v2_i - 6.0 * G2z_i * w2_i - 7.0 * G2z_i - 24.0 * Gx_i * Gy_i * u_i * v_i - 24.0 * Gx_i * Gz_i * u_i * w_i - 24.0 * Gy_i * Gz_i * v_i * w_i) / (24.0 * G);
            k_pert[18] = -sigma * (9.0 * G2 * u2_i * v2_i + 9.0 * G2 * u2_i * w2_i + 18.0 * G2 * u2_i - 9.0 * G2 * v2_i * w2_i + 5.0 * G2 - 9.0 * G2x_i * u2_i * v2_i - 9.0 * G2x_i * u2_i * w2_i - 6.0 * G2x_i * u2_i + 9.0 * G2x_i * v2_i * w2_i - 6.0 * G2x_i * v2_i - 6.0 * G2x_i * w2_i - 5.0 * G2x_i - 9.0 * G2y_i * u2_i * v2_i - 9.0 * G2y_i * u2_i * w2_i - 12.0 * G2y_i * u2_i + 9.0 * G2y_i * v2_i * w2_i + 6.0 * G2y_i * w2_i - G2y_i - 9.0 * G2z_i * u2_i * v2_i - 9.0 * G2z_i * u2_i * w2_i - 12.0 * G2z_i * u2_i + 9.0 * G2z_i * v2_i * w2_i + 6.0 * G2z_i * v2_i - G2z_i - 24.0 * Gx_i * Gy_i * u_i * v_i - 24.0 * Gx_i * Gz_i * u_i * w_i + 24.0 * Gy_i * Gz_i * v_i * w_i) / (24.0 * G);
            k_pert[19] = -sigma * (9.0 * G2 * u2_i * v2_i - 9.0 * G2 * u2_i * w2_i + 9.0 * G2 * v2_i - 9.0 * G2 * w2_i - 9.0 * G2x_i * u2_i * v2_i + 9.0 * G2x_i * u2_i * w2_i - 9.0 * G2x_i * v2_i + 9.0 * G2x_i * w2_i - 9.0 * G2y_i * u2_i * v2_i + 9.0 * G2y_i * u2_i * w2_i - 6.0 * G2y_i * u2_i - 3.0 * G2y_i * v2_i + 3.0 * G2y_i * w2_i - 2.0 * G2y_i - 9.0 * G2z_i * u2_i * v2_i + 9.0 * G2z_i * u2_i * w2_i + 6.0 * G2z_i * u2_i - 3.0 * G2z_i * v2_i + 3.0 * G2z_i * w2_i + 2.0 * G2z_i - 24.0 * Gx_i * Gy_i * u_i * v_i + 24.0 * Gx_i * Gz_i * u_i * w_i) / (24.0 * G);
            k_pert[20] = -sigma * (9.0 * G2 * u2_i * v_i * w_i + 9.0 * G2 * v_i * w_i - 9.0 * G2x_i * u2_i * v_i * w_i - 9.0 * G2x_i * v_i * w_i - 9.0 * G2y_i * u2_i * v_i * w_i - 3.0 * G2y_i * v_i * w_i - 9.0 * G2z_i * u2_i * v_i * w_i - 3.0 * G2z_i * v_i * w_i - 12.0 * Gx_i * Gy_i * u_i * w_i - 12.0 * Gx_i * Gz_i * u_i * v_i - 6.0 * Gy_i * Gz_i * u2_i - 2.0 * Gy_i * Gz_i) / (24.0 * G);
            k_pert[21] = -sigma * (9.0 * G2 * u_i * v2_i * w_i + 9.0 * G2 * u_i * w_i - 9.0 * G2x_i * u_i * v2_i * w_i - 3.0 * G2x_i * u_i * w_i - 9.0 * G2y_i * u_i * v2_i * w_i - 9.0 * G2y_i * u_i * w_i - 9.0 * G2z_i * u_i * v2_i * w_i - 3.0 * G2z_i * u_i * w_i - 12.0 * Gx_i * Gy_i * v_i * w_i - 6.0 * Gx_i * Gz_i * v2_i - 2.0 * Gx_i * Gz_i - 12.0 * Gy_i * Gz_i * u_i * v_i) / (24.0 * G);
            k_pert[22] = -sigma * (9.0 * G2 * u_i * v_i * w2_i + 9.0 * G2 * u_i * v_i - 9.0 * G2x_i * u_i * v_i * w2_i - 3.0 * G2x_i * u_i * v_i - 9.0 * G2y_i * u_i * v_i * w2_i - 3.0 * G2y_i * u_i * v_i - 9.0 * G2z_i * u_i * v_i * w2_i - 9.0 * G2z_i * u_i * v_i - 6.0 * Gx_i * Gy_i * w2_i - 2.0 * Gx_i * Gy_i - 12.0 * Gx_i * Gz_i * v_i * w_i - 12.0 * Gy_i * Gz_i * u_i * w_i) / (24.0 * G);
            k_pert[23] = sigma * (9.0 * G2 * u_i * v2_i * w2_i + 9.0 * G2 * u_i * v2_i + 9.0 * G2 * u_i * w2_i + 5.0 * G2 * u_i - 9.0 * G2x_i * u_i * v2_i * w2_i - 3.0 * G2x_i * u_i * v2_i - 3.0 * G2x_i * u_i * w2_i - G2x_i * u_i - 9.0 * G2y_i * u_i * v2_i * w2_i - 3.0 * G2y_i * u_i * v2_i - 9.0 * G2y_i * u_i * w2_i - 3.0 * G2y_i * u_i - 9.0 * G2z_i * u_i * v2_i * w2_i - 9.0 * G2z_i * u_i * v2_i - 3.0 * G2z_i * u_i * w2_i - 3.0 * G2z_i * u_i - 12.0 * Gx_i * Gy_i * v_i * w2_i - 4.0 * Gx_i * Gy_i * v_i - 12.0 * Gx_i * Gz_i * v2_i * w_i - 4.0 * Gx_i * Gz_i * w_i - 24.0 * Gy_i * Gz_i * u_i * v_i * w_i) / (24.0 * G);
            k_pert[24] = sigma * (9.0 * G2 * u2_i * v_i * w2_i + 9.0 * G2 * u2_i * v_i + 9.0 * G2 * v_i * w2_i + 5.0 * G2 * v_i - 9.0 * G2x_i * u2_i * v_i * w2_i - 3.0 * G2x_i * u2_i * v_i - 9.0 * G2x_i * v_i * w2_i - 3.0 * G2x_i * v_i - 9.0 * G2y_i * u2_i * v_i * w2_i - 3.0 * G2y_i * u2_i * v_i - 3.0 * G2y_i * v_i * w2_i - G2y_i * v_i - 9.0 * G2z_i * u2_i * v_i * w2_i - 9.0 * G2z_i * u2_i * v_i - 3.0 * G2z_i * v_i * w2_i - 3.0 * G2z_i * v_i - 12.0 * Gx_i * Gy_i * u_i * w2_i - 4.0 * Gx_i * Gy_i * u_i - 24.0 * Gx_i * Gz_i * u_i * v_i * w_i - 12.0 * Gy_i * Gz_i * u2_i * w_i - 4.0 * Gy_i * Gz_i * w_i) / (24.0 * G);
            k_pert[25] = sigma * (9.0 * G2 * u2_i * v2_i * w_i + 9.0 * G2 * u2_i * w_i + 9.0 * G2 * v2_i * w_i + 5.0 * G2 * w_i - 9.0 * G2x_i * u2_i * v2_i * w_i - 3.0 * G2x_i * u2_i * w_i - 9.0 * G2x_i * v2_i * w_i - 3.0 * G2x_i * w_i - 9.0 * G2y_i * u2_i * v2_i * w_i - 9.0 * G2y_i * u2_i * w_i - 3.0 * G2y_i * v2_i * w_i - 3.0 * G2y_i * w_i - 9.0 * G2z_i * u2_i * v2_i * w_i - 3.0 * G2z_i * u2_i * w_i - 3.0 * G2z_i * v2_i * w_i - G2z_i * w_i - 24.0 * Gx_i * Gy_i * u_i * v_i * w_i - 12.0 * Gx_i * Gz_i * u_i * v2_i - 4.0 * Gx_i * Gz_i * u_i - 12.0 * Gy_i * Gz_i * u2_i * v_i - 4.0 * Gy_i * Gz_i * v_i) / (24.0 * G);
            k_pert[26] = -sigma * (27.0 * G2 * u2_i * v2_i * w2_i + 27.0 * G2 * u2_i * v2_i + 27.0 * G2 * u2_i * w2_i + 15.0 * G2 * u2_i + 27.0 * G2 * v2_i * w2_i + 15.0 * G2 * v2_i + 15.0 * G2 * w2_i + 7.0 * G2 - 27.0 * G2x_i * u2_i * v2_i * w2_i - 9.0 * G2x_i * u2_i * v2_i - 9.0 * G2x_i * u2_i * w2_i - 3.0 * G2x_i * u2_i - 27.0 * G2x_i * v2_i * w2_i - 9.0 * G2x_i * v2_i - 9.0 * G2x_i * w2_i - 3.0 * G2x_i - 27.0 * G2y_i * u2_i * v2_i * w2_i - 9.0 * G2y_i * u2_i * v2_i - 27.0 * G2y_i * u2_i * w2_i - 9.0 * G2y_i * u2_i - 9.0 * G2y_i * v2_i * w2_i - 3.0 * G2y_i * v2_i - 9.0 * G2y_i * w2_i - 3.0 * G2y_i - 27.0 * G2z_i * u2_i * v2_i * w2_i - 27.0 * G2z_i * u2_i * v2_i - 9.0 * G2z_i * u2_i * w2_i - 9.0 * G2z_i * u2_i - 9.0 * G2z_i * v2_i * w2_i - 9.0 * G2z_i * v2_i - 3.0 * G2z_i * w2_i - 3.0 * G2z_i - 72.0 * Gx_i * Gy_i * u_i * v_i * w2_i - 24.0 * Gx_i * Gy_i * u_i * v_i - 72.0 * Gx_i * Gz_i * u_i * v2_i * w_i - 24.0 * Gx_i * Gz_i * u_i * w_i - 72.0 * Gy_i * Gz_i * u2_i * v_i * w_i - 24.0 * Gy_i * Gz_i * v_i * w_i) / (72.0 * G);
        }

        for (int n = 0; n < NCOMP; n++)
        {
            rho_i = rho_comp[INDEX(i, j, k, n)];

            Fx_i = Fx[INDEX(i, j, k, n)];
            Fy_i = Fy[INDEX(i, j, k, n)];
            Fz_i = Fz[INDEX(i, j, k, n)];

            cs2_i = cs2[n];
            cs4_i = cs2_i * cs2_i;
            cs6_i = cs2_i * cs4_i;

            k4 = 0.0;
            k5 = 0.0;
            k6 = 0.0;
            k7 = 0.0;
            k8 = 0.0;

            for (int p = 0; p < NP; p++)
            {
                cx_bar = (double)cx[p] - u_i;
                cy_bar = (double)cy[p] - v_i;
                cz_bar = (double)cz[p] - w_i;

                k4 += f(p) * cx_bar * cy_bar;
                k5 += f(p) * cx_bar * cz_bar;
                k6 += f(p) * cy_bar * cz_bar;
                k7 += f(p) * (cx_bar * cx_bar - cy_bar * cy_bar);
                k8 += f(p) * (cx_bar * cx_bar - cz_bar * cz_bar);
            }

            k_star[0] = k_pert[0] + rho_i;
            k_star[1] = Fx_i / 2.0 + k_pert[1];
            k_star[2] = Fy_i / 2.0 + k_pert[2];
            k_star[3] = Fz_i / 2.0 + k_pert[3];
            k_star[4] = (1.0 - omega) * k4 + k_pert[4] * omega;
            k_star[5] = (1.0 - omega) * k5 + k_pert[5] * omega;
            k_star[6] = (1.0 - omega) * k6 + k_pert[6] * omega;
            k_star[7] = (1.0 - omega) * k7 + k_pert[7] * omega;
            k_star[8] = (1.0 - omega) * k8 + k_pert[8] * omega;
            k_star[9] = 3.0 * cs2_i * rho_i + k_pert[9];
            k_star[10] = Fx_i * cs2_i + k_pert[10];
            k_star[11] = Fy_i * cs2_i + k_pert[11];
            k_star[12] = Fz_i * cs2_i + k_pert[12];
            k_star[13] = k_pert[13];
            k_star[14] = k_pert[14];
            k_star[15] = k_pert[15];
            k_star[16] = k_pert[16];
            k_star[17] = cs2_i * rho_i + k_pert[17];
            k_star[18] = cs4_i * rho_i + k_pert[18];
            k_star[19] = k_pert[19];
            k_star[20] = k_pert[20];
            k_star[21] = k_pert[21];
            k_star[22] = k_pert[22];
            k_star[23] = Fx_i * cs4_i / 2.0 + k_pert[23];
            k_star[24] = Fy_i * cs4_i / 2.0 + k_pert[24];
            k_star[25] = Fz_i * cs4_i / 2.0 + k_pert[25];
            k_star[26] = cs6_i * rho_i + k_pert[26];

            raw[0] = k_star[0];
            raw[1] = k_star[0] * u_i + k_star[1];
            raw[2] = k_star[0] * v_i + k_star[2];
            raw[3] = k_star[0] * w_i + k_star[3];
            raw[4] = k_star[0] * u_i * v_i + k_star[1] * v_i + k_star[2] * u_i + k_star[4];
            raw[5] = k_star[0] * u_i * w_i + k_star[1] * w_i + k_star[3] * u_i + k_star[5];
            raw[6] = k_star[0] * v_i * w_i + k_star[2] * w_i + k_star[3] * v_i + k_star[6];
            raw[7] = k_star[0] * (u2_i - v2_i) + 2.0 * k_star[1] * u_i - 2.0 * k_star[2] * v_i + k_star[7];
            raw[8] = k_star[0] * (u2_i - w2_i) + 2.0 * k_star[1] * u_i - 2.0 * k_star[3] * w_i + k_star[8];
            raw[9] = k_star[0] * (u2_i + v2_i + w2_i) + 2.0 * k_star[1] * u_i + 2.0 * k_star[2] * v_i + 2.0 * k_star[3] * w_i + k_star[9];
            raw[10] = k_star[0] * u_i * (v2_i + w2_i) + k_star[10] + k_star[1] * (v2_i + w2_i) + 2.0 * k_star[2] * u_i * v_i + 2.0 * k_star[3] * u_i * w_i + 2.0 * k_star[4] * v_i + 2.0 * k_star[5] * w_i - k_star[7] * u_i / 3.0 - k_star[8] * u_i / 3.0 + 2.0 * k_star[9] * u_i / 3.0;
            raw[11] = k_star[0] * v_i * (u2_i + w2_i) + k_star[11] + 2.0 * k_star[1] * u_i * v_i + k_star[2] * (u2_i + w2_i) + 2.0 * k_star[3] * v_i * w_i + 2.0 * k_star[4] * u_i + 2.0 * k_star[6] * w_i + 2.0 * k_star[7] * v_i / 3.0 - k_star[8] * v_i / 3.0 + 2.0 * k_star[9] * v_i / 3.0;
            raw[12] = k_star[0] * w_i * (u2_i + v2_i) + k_star[12] + 2.0 * k_star[1] * u_i * w_i + 2.0 * k_star[2] * v_i * w_i + k_star[3] * (u2_i + v2_i) + 2.0 * k_star[5] * u_i + 2.0 * k_star[6] * v_i - k_star[7] * w_i / 3.0 + 2.0 * k_star[8] * w_i / 3.0 + 2.0 * k_star[9] * w_i / 3.0;
            raw[13] = k_star[0] * u_i * (v2_i - w2_i) + k_star[13] + k_star[1] * (v2_i - w2_i) + 2.0 * k_star[2] * u_i * v_i - 2.0 * k_star[3] * u_i * w_i + 2.0 * k_star[4] * v_i - 2.0 * k_star[5] * w_i - k_star[7] * u_i + k_star[8] * u_i;
            raw[14] = k_star[0] * v_i * (u2_i - w2_i) + k_star[14] + 2.0 * k_star[1] * u_i * v_i + k_star[2] * (u2_i - w2_i) - 2.0 * k_star[3] * v_i * w_i + 2.0 * k_star[4] * u_i - 2.0 * k_star[6] * w_i + k_star[8] * v_i;
            raw[15] = k_star[0] * w_i * (u2_i - v2_i) + k_star[15] + 2.0 * k_star[1] * u_i * w_i - 2.0 * k_star[2] * v_i * w_i + k_star[3] * (u2_i - v2_i) + 2.0 * k_star[5] * u_i - 2.0 * k_star[6] * v_i + k_star[7] * w_i;
            raw[16] = k_star[0] * u_i * v_i * w_i + k_star[16] + k_star[1] * v_i * w_i + k_star[2] * u_i * w_i + k_star[3] * u_i * v_i + k_star[4] * w_i + k_star[5] * v_i + k_star[6] * u_i;
            raw[17] = k_star[0] * (u2_i * v2_i + u2_i * w2_i + v2_i * w2_i) + 2.0 * k_star[10] * u_i + 2.0 * k_star[11] * v_i + 2.0 * k_star[12] * w_i + k_star[17] + 2.0 * k_star[1] * u_i * (v2_i + w2_i) + 2.0 * k_star[2] * v_i * (u2_i + w2_i) + 2.0 * k_star[3] * w_i * (u2_i + v2_i) + 4.0 * k_star[4] * u_i * v_i + 4.0 * k_star[5] * u_i * w_i + 4.0 * k_star[6] * v_i * w_i + k_star[7] * (-u2_i / 3.0 + 2.0 * v2_i / 3.0 - w2_i / 3.0) + k_star[8] * (-u2_i / 3.0 - v2_i / 3.0 + 2.0 * w2_i / 3.0) + k_star[9] * (2.0 * u2_i / 3.0 + 2.0 * v2_i / 3.0 + 2.0 * w2_i / 3.0);
            raw[18] = k_star[0] * (u2_i * v2_i + u2_i * w2_i - v2_i * w2_i) + 2.0 * k_star[10] * u_i + 2.0 * k_star[14] * v_i + 2.0 * k_star[15] * w_i + k_star[18] + 2.0 * k_star[1] * u_i * (v2_i + w2_i) + 2.0 * k_star[2] * v_i * (u2_i - w2_i) + 2.0 * k_star[3] * w_i * (u2_i - v2_i) + 4.0 * k_star[4] * u_i * v_i + 4.0 * k_star[5] * u_i * w_i - 4.0 * k_star[6] * v_i * w_i + k_star[7] * (-u2_i / 3.0 + w2_i) + k_star[8] * (-u2_i / 3.0 + v2_i) + 2.0 * k_star[9] * u2_i / 3.0;
            raw[19] = k_star[0] * u2_i * (v2_i - w2_i) + k_star[11] * v_i - k_star[12] * w_i + 2.0 * k_star[13] * u_i + k_star[14] * v_i - k_star[15] * w_i + k_star[19] + 2.0 * k_star[1] * u_i * (v2_i - w2_i) + 2.0 * k_star[2] * u2_i * v_i - 2.0 * k_star[3] * u2_i * w_i + 4.0 * k_star[4] * u_i * v_i - 4.0 * k_star[5] * u_i * w_i + k_star[7] * (-u2_i + v2_i / 3.0 - w2_i / 3.0) + k_star[8] * (u2_i + v2_i / 3.0 - w2_i / 3.0) + k_star[9] * (v2_i / 3.0 - w2_i / 3.0);
            raw[20] = k_star[0] * u2_i * v_i * w_i + k_star[11] * w_i / 2.0 + k_star[12] * v_i / 2.0 + k_star[14] * w_i / 2.0 + k_star[15] * v_i / 2.0 + 2.0 * k_star[16] * u_i + 2.0 * k_star[1] * u_i * v_i * w_i + k_star[20] + k_star[2] * u2_i * w_i + k_star[3] * u2_i * v_i + 2.0 * k_star[4] * u_i * w_i + 2.0 * k_star[5] * u_i * v_i + k_star[6] * u2_i + k_star[7] * v_i * w_i / 3.0 + k_star[8] * v_i * w_i / 3.0 + k_star[9] * v_i * w_i / 3.0;
            raw[21] = k_star[0] * u_i * v2_i * w_i + k_star[10] * w_i / 2.0 + k_star[12] * u_i / 2.0 + k_star[13] * w_i / 2.0 - k_star[15] * u_i / 2.0 + 2.0 * k_star[16] * v_i + k_star[1] * v2_i * w_i + k_star[21] + 2.0 * k_star[2] * u_i * v_i * w_i + k_star[3] * u_i * v2_i + 2.0 * k_star[4] * v_i * w_i + k_star[5] * v2_i + 2.0 * k_star[6] * u_i * v_i - 2.0 * k_star[7] * u_i * w_i / 3.0 + k_star[8] * u_i * w_i / 3.0 + k_star[9] * u_i * w_i / 3.0;
            raw[22] = k_star[0] * u_i * v_i * w2_i + k_star[10] * v_i / 2.0 + k_star[11] * u_i / 2.0 - k_star[13] * v_i / 2.0 - k_star[14] * u_i / 2.0 + 2.0 * k_star[16] * w_i + k_star[1] * v_i * w2_i + k_star[22] + k_star[2] * u_i * w2_i + 2.0 * k_star[3] * u_i * v_i * w_i + k_star[4] * w2_i + 2.0 * k_star[5] * v_i * w_i + 2.0 * k_star[6] * u_i * w_i + k_star[7] * u_i * v_i / 3.0 - 2.0 * k_star[8] * u_i * v_i / 3.0 + k_star[9] * u_i * v_i / 3.0;
            raw[23] = k_star[0] * u_i * v2_i * w2_i + k_star[10] * (v2_i / 2.0 + w2_i / 2.0) + k_star[11] * u_i * v_i + k_star[12] * u_i * w_i + k_star[13] * (-v2_i / 2.0 + w2_i / 2.0) - k_star[14] * u_i * v_i - k_star[15] * u_i * w_i + 4.0 * k_star[16] * v_i * w_i + k_star[17] * u_i / 2.0 - k_star[18] * u_i / 2.0 + k_star[1] * v2_i * w2_i + 2.0 * k_star[21] * w_i + 2.0 * k_star[22] * v_i + k_star[23] + 2.0 * k_star[2] * u_i * v_i * w2_i + 2.0 * k_star[3] * u_i * v2_i * w_i + 2.0 * k_star[4] * v_i * w2_i + 2.0 * k_star[5] * v2_i * w_i + 4.0 * k_star[6] * u_i * v_i * w_i + k_star[7] * u_i * (v2_i - 2.0 * w2_i) / 3.0 - k_star[8] * u_i * (2.0 * v2_i - w2_i) / 3.0 + k_star[9] * u_i * (v2_i + w2_i) / 3.0;
            raw[24] = k_star[0] * u2_i * v_i * w2_i + k_star[10] * u_i * v_i + k_star[11] * (u2_i / 2.0 + w2_i / 2.0) + k_star[12] * v_i * w_i - k_star[13] * u_i * v_i + k_star[14] * (-u2_i / 2.0 + w2_i / 2.0) + k_star[15] * v_i * w_i + 4.0 * k_star[16] * u_i * w_i + k_star[17] * v_i / 4.0 + k_star[18] * v_i / 4.0 - k_star[19] * v_i / 2.0 + 2.0 * k_star[1] * u_i * v_i * w2_i + 2.0 * k_star[20] * w_i + 2.0 * k_star[22] * u_i + k_star[24] + k_star[2] * u2_i * w2_i + 2.0 * k_star[3] * u2_i * v_i * w_i + 2.0 * k_star[4] * u_i * w2_i + 4.0 * k_star[5] * u_i * v_i * w_i + 2.0 * k_star[6] * u2_i * w_i + k_star[7] * v_i * (u2_i + w2_i) / 3.0 - k_star[8] * v_i * (2.0 * u2_i - w2_i) / 3.0 + k_star[9] * v_i * (u2_i + w2_i) / 3.0;
            raw[25] = k_star[0] * u2_i * v2_i * w_i + k_star[10] * u_i * w_i + k_star[11] * v_i * w_i + k_star[12] * (u2_i / 2.0 + v2_i / 2.0) + k_star[13] * u_i * w_i + k_star[14] * v_i * w_i + k_star[15] * (-u2_i / 2.0 + v2_i / 2.0) + 4.0 * k_star[16] * u_i * v_i + k_star[17] * w_i / 4.0 + k_star[18] * w_i / 4.0 + k_star[19] * w_i / 2.0 + 2.0 * k_star[1] * u_i * v2_i * w_i + 2.0 * k_star[20] * v_i + 2.0 * k_star[21] * u_i + k_star[25] + 2.0 * k_star[2] * u2_i * v_i * w_i + k_star[3] * u2_i * v2_i + 4.0 * k_star[4] * u_i * v_i * w_i + 2.0 * k_star[5] * u_i * v2_i + 2.0 * k_star[6] * u2_i * v_i - k_star[7] * w_i * (2.0 * u2_i - v2_i) / 3.0 + k_star[8] * w_i * (u2_i + v2_i) / 3.0 + k_star[9] * w_i * (u2_i + v2_i) / 3.0;
            raw[26] = k_star[0] * u2_i * v2_i * w2_i + k_star[10] * u_i * (v2_i + w2_i) + k_star[11] * v_i * (u2_i + w2_i) + k_star[12] * w_i * (u2_i + v2_i) - k_star[13] * u_i * (v2_i - w2_i) - k_star[14] * v_i * (u2_i - w2_i) - k_star[15] * w_i * (u2_i - v2_i) + 8.0 * k_star[16] * u_i * v_i * w_i + k_star[17] * (u2_i / 2.0 + v2_i / 4.0 + w2_i / 4.0) + k_star[18] * (-u2_i / 2.0 + v2_i / 4.0 + w2_i / 4.0) + k_star[19] * (-v2_i / 2.0 + w2_i / 2.0) + 2.0 * k_star[1] * u_i * v2_i * w2_i + 4.0 * k_star[20] * v_i * w_i + 4.0 * k_star[21] * u_i * w_i + 4.0 * k_star[22] * u_i * v_i + 2.0 * k_star[23] * u_i + 2.0 * k_star[24] * v_i + 2.0 * k_star[25] * w_i + k_star[26] + 2.0 * k_star[2] * u2_i * v_i * w2_i + 2.0 * k_star[3] * u2_i * v2_i * w_i + 4.0 * k_star[4] * u_i * v_i * w2_i + 4.0 * k_star[5] * u_i * v2_i * w_i + 4.0 * k_star[6] * u2_i * v_i * w_i + k_star[7] * (u2_i * v2_i / 3.0 - 2.0 * u2_i * w2_i / 3.0 + v2_i * w2_i / 3.0) + k_star[8] * (-2.0 * u2_i * v2_i / 3.0 + u2_i * w2_i / 3.0 + v2_i * w2_i / 3.0) + k_star[9] * (u2_i * v2_i / 3.0 + u2_i * w2_i / 3.0 + v2_i * w2_i / 3.0);

            f2[INDEX_F(i, j, k, 0, n)] = raw[0] + raw[17] - raw[26] - raw[9];
            f2[INDEX_F(i, j, k, 1, n)] = -raw[10] / 2.0 - raw[17] / 4.0 - raw[18] / 4.0 + raw[1] / 2.0 + raw[23] / 2.0 + raw[26] / 2.0 + raw[7] / 6.0 + raw[8] / 6.0 + raw[9] / 6.0;
            f2[INDEX_F(i, j, k, 2, n)] = raw[10] / 2.0 - raw[17] / 4.0 - raw[18] / 4.0 - raw[1] / 2.0 - raw[23] / 2.0 + raw[26] / 2.0 + raw[7] / 6.0 + raw[8] / 6.0 + raw[9] / 6.0;
            f2[INDEX_F(i, j, k, 3, n)] = -raw[11] / 2.0 - 3.0 * raw[17] / 8.0 + raw[18] / 8.0 - raw[19] / 4.0 + raw[24] / 2.0 + raw[26] / 2.0 + raw[2] / 2.0 - raw[7] / 3.0 + raw[8] / 6.0 + raw[9] / 6.0;
            f2[INDEX_F(i, j, k, 4, n)] = raw[11] / 2.0 - 3.0 * raw[17] / 8.0 + raw[18] / 8.0 - raw[19] / 4.0 - raw[24] / 2.0 + raw[26] / 2.0 - raw[2] / 2.0 - raw[7] / 3.0 + raw[8] / 6.0 + raw[9] / 6.0;
            f2[INDEX_F(i, j, k, 5, n)] = -raw[12] / 2.0 - 3.0 * raw[17] / 8.0 + raw[18] / 8.0 + raw[19] / 4.0 + raw[25] / 2.0 + raw[26] / 2.0 + raw[3] / 2.0 + raw[7] / 6.0 - raw[8] / 3.0 + raw[9] / 6.0;
            f2[INDEX_F(i, j, k, 6, n)] = raw[12] / 2.0 - 3.0 * raw[17] / 8.0 + raw[18] / 8.0 + raw[19] / 4.0 - raw[25] / 2.0 + raw[26] / 2.0 - raw[3] / 2.0 + raw[7] / 6.0 - raw[8] / 3.0 + raw[9] / 6.0;
            f2[INDEX_F(i, j, k, 7, n)] = raw[10] / 8.0 + raw[11] / 8.0 + raw[13] / 8.0 + raw[14] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 + raw[19] / 8.0 - raw[22] / 4.0 - raw[23] / 4.0 - raw[24] / 4.0 - raw[26] / 4.0 + raw[4] / 4.0;
            f2[INDEX_F(i, j, k, 8, n)] = -raw[10] / 8.0 + raw[11] / 8.0 - raw[13] / 8.0 + raw[14] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 + raw[19] / 8.0 + raw[22] / 4.0 + raw[23] / 4.0 - raw[24] / 4.0 - raw[26] / 4.0 - raw[4] / 4.0;
            f2[INDEX_F(i, j, k, 9, n)] = raw[10] / 8.0 - raw[11] / 8.0 + raw[13] / 8.0 - raw[14] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 + raw[19] / 8.0 + raw[22] / 4.0 - raw[23] / 4.0 + raw[24] / 4.0 - raw[26] / 4.0 - raw[4] / 4.0;
            f2[INDEX_F(i, j, k, 10, n)] = -raw[10] / 8.0 - raw[11] / 8.0 - raw[13] / 8.0 - raw[14] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 + raw[19] / 8.0 - raw[22] / 4.0 + raw[23] / 4.0 + raw[24] / 4.0 - raw[26] / 4.0 + raw[4] / 4.0;
            f2[INDEX_F(i, j, k, 11, n)] = raw[10] / 8.0 + raw[12] / 8.0 - raw[13] / 8.0 + raw[15] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 - raw[19] / 8.0 - raw[21] / 4.0 - raw[23] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 + raw[5] / 4.0;
            f2[INDEX_F(i, j, k, 12, n)] = -raw[10] / 8.0 + raw[12] / 8.0 + raw[13] / 8.0 + raw[15] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 - raw[19] / 8.0 + raw[21] / 4.0 + raw[23] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 - raw[5] / 4.0;
            f2[INDEX_F(i, j, k, 13, n)] = raw[10] / 8.0 - raw[12] / 8.0 - raw[13] / 8.0 - raw[15] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 - raw[19] / 8.0 + raw[21] / 4.0 - raw[23] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 - raw[5] / 4.0;
            f2[INDEX_F(i, j, k, 14, n)] = -raw[10] / 8.0 - raw[12] / 8.0 + raw[13] / 8.0 - raw[15] / 8.0 + raw[17] / 16.0 + raw[18] / 16.0 - raw[19] / 8.0 - raw[21] / 4.0 + raw[23] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 + raw[5] / 4.0;
            f2[INDEX_F(i, j, k, 15, n)] = raw[11] / 8.0 + raw[12] / 8.0 - raw[14] / 8.0 - raw[15] / 8.0 + raw[17] / 8.0 - raw[18] / 8.0 - raw[20] / 4.0 - raw[24] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 + raw[6] / 4.0;
            f2[INDEX_F(i, j, k, 16, n)] = -raw[11] / 8.0 + raw[12] / 8.0 + raw[14] / 8.0 - raw[15] / 8.0 + raw[17] / 8.0 - raw[18] / 8.0 + raw[20] / 4.0 + raw[24] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 - raw[6] / 4.0;
            f2[INDEX_F(i, j, k, 17, n)] = raw[11] / 8.0 - raw[12] / 8.0 - raw[14] / 8.0 + raw[15] / 8.0 + raw[17] / 8.0 - raw[18] / 8.0 + raw[20] / 4.0 - raw[24] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 - raw[6] / 4.0;
            f2[INDEX_F(i, j, k, 18, n)] = -raw[11] / 8.0 - raw[12] / 8.0 + raw[14] / 8.0 + raw[15] / 8.0 + raw[17] / 8.0 - raw[18] / 8.0 - raw[20] / 4.0 + raw[24] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 + raw[6] / 4.0;
            f2[INDEX_F(i, j, k, 19, n)] = raw[16] / 8.0 + raw[20] / 8.0 + raw[21] / 8.0 + raw[22] / 8.0 + raw[23] / 8.0 + raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
            f2[INDEX_F(i, j, k, 20, n)] = -raw[16] / 8.0 + raw[20] / 8.0 - raw[21] / 8.0 - raw[22] / 8.0 - raw[23] / 8.0 + raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
            f2[INDEX_F(i, j, k, 21, n)] = -raw[16] / 8.0 - raw[20] / 8.0 + raw[21] / 8.0 - raw[22] / 8.0 + raw[23] / 8.0 - raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
            f2[INDEX_F(i, j, k, 22, n)] = raw[16] / 8.0 - raw[20] / 8.0 - raw[21] / 8.0 + raw[22] / 8.0 - raw[23] / 8.0 - raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
            f2[INDEX_F(i, j, k, 23, n)] = -raw[16] / 8.0 - raw[20] / 8.0 - raw[21] / 8.0 + raw[22] / 8.0 + raw[23] / 8.0 + raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
            f2[INDEX_F(i, j, k, 24, n)] = raw[16] / 8.0 - raw[20] / 8.0 + raw[21] / 8.0 - raw[22] / 8.0 - raw[23] / 8.0 + raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
            f2[INDEX_F(i, j, k, 25, n)] = raw[16] / 8.0 + raw[20] / 8.0 - raw[21] / 8.0 - raw[22] / 8.0 + raw[23] / 8.0 - raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
            f2[INDEX_F(i, j, k, 26, n)] = -raw[16] / 8.0 + raw[20] / 8.0 + raw[21] / 8.0 + raw[22] / 8.0 - raw[23] / 8.0 - raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
        }

        if (G > 1e-15)
        {
            rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
            rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];
            rho_i = rho_RED_i + rho_BLUE_i;

            for (int p = 0; p < NP; p++)
            {
                f_star = f2[INDEX_F(i, j, k, p, RED)] + f2[INDEX_F(i, j, k, p, BLUE)];
                f2[INDEX_F(i, j, k, p, RED)] = rho_RED_i / rho_i * f_star;
                f2[INDEX_F(i, j, k, p, BLUE)] = rho_BLUE_i / rho_i * f_star;
            }

            // Compute total stationary equilibrium
            compute_stationary_equilibrium(rho_RED_i, cs2[RED], raw);
            compute_stationary_equilibrium(rho_BLUE_i, cs2[BLUE], feq);
            for (int p = 1; p < NP; p++)
            {
                feq[p] += raw[p];
            }

            for (int p = 1; p < NP; p++)
            {
                Gc = Gx_i * (double)cx[p] + Gy_i * (double)cy[p] + Gz_i * (double)cz[p];
                cos_phi = Gc / (sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]) * G);

                p_rec = beta * rho_RED_i * rho_BLUE_i / (rho_i * rho_i) * cos_phi * feq[p];

                f2[INDEX_F(i, j, k, p, RED)] += p_rec;
                f2[INDEX_F(i, j, k, p, BLUE)] -= p_rec;
            }
        }
    }
}

void compute_stationary_equilibrium(double rho, double cs2, double *feq)
{
    double cs4 = cs2 * cs2;
    double cs6 = cs2 * cs4;

    feq[0] = rho * (-2.0 * cs2 - cs6 + 1.0);
    feq[1] = rho * (cs2 - cs4 + 2.0 * cs6) / 4.0;
    feq[2] = rho * (cs2 - cs4 + 2.0 * cs6) / 4.0;
    feq[3] = rho * (cs2 + cs4 + 4.0 * cs6) / 8.0;
    feq[4] = rho * (cs2 + cs4 + 4.0 * cs6) / 8.0;
    feq[5] = rho * (cs2 + cs4 + 4.0 * cs6) / 8.0;
    feq[6] = rho * (cs2 + cs4 + 4.0 * cs6) / 8.0;
    feq[7] = rho * (cs2 + cs4 - 4.0 * cs6) / 16.0;
    feq[8] = rho * (cs2 + cs4 - 4.0 * cs6) / 16.0;
    feq[9] = rho * (cs2 + cs4 - 4.0 * cs6) / 16.0;
    feq[10] = rho * (cs2 + cs4 - 4.0 * cs6) / 16.0;
    feq[11] = rho * (cs2 + cs4 - 4.0 * cs6) / 16.0;
    feq[12] = rho * (cs2 + cs4 - 4.0 * cs6) / 16.0;
    feq[13] = rho * (cs2 + cs4 - 4.0 * cs6) / 16.0;
    feq[14] = rho * (cs2 + cs4 - 4.0 * cs6) / 16.0;
    feq[15] = rho * (cs2 - cs4 - 2.0 * cs6) / 8.0;
    feq[16] = rho * (cs2 - cs4 - 2.0 * cs6) / 8.0;
    feq[17] = rho * (cs2 - cs4 - 2.0 * cs6) / 8.0;
    feq[18] = rho * (cs2 - cs4 - 2.0 * cs6) / 8.0;
    feq[19] = cs6 * rho / 8.0;
    feq[20] = cs6 * rho / 8.0;
    feq[21] = cs6 * rho / 8.0;
    feq[22] = cs6 * rho / 8.0;
    feq[23] = cs6 * rho / 8.0;
    feq[24] = cs6 * rho / 8.0;
    feq[25] = cs6 * rho / 8.0;
    feq[26] = cs6 * rho / 8.0;
}