#include <math.h>
#include <string.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"

void compute_equilibrium(double rho, double u, double v, double w, double P, double *feq, SimulationBag *sim)
{
    DistributionBag *dists = sim->dists;
    double *meq = dists->meq;
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;

    meq[0] = rho;
    meq[1] = rho * u;
    meq[2] = rho * v;
    meq[3] = rho * w;
    meq[4] = rho * u * v;
    meq[5] = rho * u * w;
    meq[6] = rho * v * w;
    meq[7] = rho * (u2 - v2);
    meq[8] = rho * (u2 - w2);
    meq[9] = 3.0 * P + rho * (u2 + v2 + w2);
    meq[10] = u * (P + rho * v2);
    meq[11] = u * (P + rho * w2);
    meq[12] = v * (P + rho * w2);
    meq[13] = v * (P + rho * u2);
    meq[14] = w * (P + rho * u2);
    meq[15] = w * (P + rho * v2);
    meq[16] = rho * u * v * w;
    meq[17] = P * (u2 + v2) + P / 3.0 + rho * u2 * v2;
    meq[18] = P * (u2 + w2) + P / 3.0 + rho * u2 * w2;
    meq[19] = P * (v2 + w2) + P / 3.0 + rho * v2 * w2;
    meq[20] = v * w * (P + rho * u2);
    meq[21] = u * w * (P + rho * v2);
    meq[22] = u * v * (P + rho * w2);
    meq[23] = u * (3.0 * P * (v2 + w2) + P + 3.0 * rho * v2 * w2) / 3.0;
    meq[24] = v * (3.0 * P * (u2 + w2) + P + 3.0 * rho * u2 * w2) / 3.0;
    meq[25] = w * (3.0 * P * (u2 + v2) + P + 3.0 * rho * u2 * v2) / 3.0;
    meq[26] = P * u2 / 3.0 + P * v2 / 3.0 + P * w2 / 3.0 + P * (u2 * v2 + u2 * w2 + v2 * w2) + P / 9.0 + rho * u2 * v2 * w2;

    raw_to_pop(meq, feq);
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
    double n_star_x, n_star_y, n_star_z;
    double n_plus_x, n_plus_y, n_plus_z;
    double n_min_x, n_min_y, n_min_z;
    double theta_prime;
    double dnx, dny, dnz, D_plus, D_min;
#endif

    int *flag = glob_fields->flag;

    double *rho_N = glob_fields->rho_N;
    double *G_norm = glob_fields->G_norm;
    double *Gx = glob_fields->Gx;
    double *Gy = glob_fields->Gy;
    double *Gz = glob_fields->Gz;
    double *nx = glob_fields->nx;
    double *ny = glob_fields->ny;
    double *nz = glob_fields->nz;

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

        G_i = sqrt(Gx_i * Gx_i + Gy_i * Gy_i + Gz_i * Gz_i);

#ifdef THETA_C
        nhat_norm = sqrt(nhat_x * nhat_x + nhat_y * nhat_y + nhat_z * nhat_z);

        if (nhat_norm > 0.0)
        {
            nhat_x /= nhat_norm;
            nhat_y /= nhat_norm;
            nhat_z /= nhat_norm;

            // n_star_x = -Gx_i / G_i;
            // n_star_y = -Gy_i / G_i;
            // n_star_z = -Gz_i / G_i;

            // theta_prime = acos(nhat_x * n_star_x + nhat_y * n_star_y + nhat_z * n_star_z);

            // n_plus_x = (cos(theta) - sin(theta) * cos(theta_prime) / sin(theta_prime)) * nhat_x + sin(theta) / sin(theta_prime) * n_star_x;
            // n_plus_y = (cos(theta) - sin(theta) * cos(theta_prime) / sin(theta_prime)) * nhat_y + sin(theta) / sin(theta_prime) * n_star_y;
            // n_plus_z = (cos(theta) - sin(theta) * cos(theta_prime) / sin(theta_prime)) * nhat_z + sin(theta) / sin(theta_prime) * n_star_z;

            // n_min_x = (cos(-theta) - sin(-theta) * cos(theta_prime) / sin(theta_prime)) * nhat_x + sin(-theta) / sin(theta_prime) * n_star_x;
            // n_min_y = (cos(-theta) - sin(-theta) * cos(theta_prime) / sin(theta_prime)) * nhat_y + sin(-theta) / sin(theta_prime) * n_star_y;
            // n_min_z = (cos(-theta) - sin(-theta) * cos(theta_prime) / sin(theta_prime)) * nhat_z + sin(-theta) / sin(theta_prime) * n_star_z;

            // dnx = n_plus_x - nhat_x;
            // dny = n_plus_y - nhat_y;
            // dnz = n_plus_z - nhat_z;
            // D_plus = sqrt(dnx * dnx + dny * dny + dnz * dnz);

            // dnx = n_min_x - nhat_x;
            // dny = n_min_y - nhat_y;
            // dnz = n_min_z - nhat_z;
            // D_min = sqrt(dnx * dnx + dny * dny + dnz * dnz);

            // if (D_plus <= D_min)
            // {
            //     Gx_i = -G_i * n_plus_x;
            //     Gy_i = -G_i * n_plus_y;
            //     Gz_i = -G_i * n_plus_z;
            // }
            // else
            // {
            //     Gx_i = -G_i * n_min_x;
            //     Gy_i = -G_i * n_min_y;
            //     Gz_i = -G_i * n_min_z;
            // }

            Gn = nhat_x * Gx_i + nhat_y * Gy_i + nhat_z * Gz_i;
            Gx_i -= Gn * nhat_x;
            Gy_i -= Gn * nhat_y;
            Gz_i -= Gn * nhat_z;

            Gn = tan(theta - 0.5 * DS_PI) * sqrt(Gx_i * Gx_i + Gy_i * Gy_i + Gz_i * Gz_i);

            Gx_i += Gn * nhat_x;
            Gy_i += Gn * nhat_y;
            Gz_i += Gn * nhat_z;

            G_i = sqrt(Gx_i * Gx_i + Gy_i * Gy_i + Gz_i * Gz_i);
        }
#endif

        G_norm[INDEX_GLOB(i, j, k)] = G_i;
        Gx[INDEX_GLOB(i, j, k)] = Gx_i;
        Gy[INDEX_GLOB(i, j, k)] = Gy_i;
        Gz[INDEX_GLOB(i, j, k)] = Gz_i;
        nx[INDEX_GLOB(i, j, k)] = Gx_i / (G_i + 1e-15);
        ny[INDEX_GLOB(i, j, k)] = Gy_i / (G_i + 1e-15);
        nz[INDEX_GLOB(i, j, k)] = Gz_i / (G_i + 1e-15);
    }
}

void compute_Q_corrections(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    Stencil *stencil = sim->stencil;

    double Qx_i, Qy_i, Qz_i;
    double qx, qy, qz;

    int ic, jc, kc;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int *flag = glob_fields->flag;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double *Qx = glob_fields->Qx;
    double *Qy = glob_fields->Qy;
    double *Qz = glob_fields->Qz;

    double *rho = glob_fields->rho;
    double *pressure = glob_fields->pressure;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    FOR_DOMAIN
    {
        Qx_i = 0.0;
        Qy_i = 0.0;
        Qz_i = 0.0;

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
                qx = extrapolate_wall_q(ic, jc, kc, 0, sim);
                qy = extrapolate_wall_q(ic, jc, kc, 1, sim);
                qz = extrapolate_wall_q(ic, jc, kc, 2, sim);
                goto skip;
            }

            qx = (pressure[INDEX_GLOB(ic, jc, kc)] - rho[INDEX_GLOB(ic, jc, kc)] / 3.0) * u[INDEX_GLOB(ic, jc, kc)];
            qy = (pressure[INDEX_GLOB(ic, jc, kc)] - rho[INDEX_GLOB(ic, jc, kc)] / 3.0) * v[INDEX_GLOB(ic, jc, kc)];
            qz = (pressure[INDEX_GLOB(ic, jc, kc)] - rho[INDEX_GLOB(ic, jc, kc)] / 3.0) * w[INDEX_GLOB(ic, jc, kc)];

        skip:

            Qx_i += wp[p] * qx * (double)cx[p];
            Qy_i += wp[p] * qy * (double)cy[p];
            Qz_i += wp[p] * qz * (double)cz[p];
        }

        Qx[INDEX_GLOB(i, j, k)] = -9.0 * Qx_i;
        Qy[INDEX_GLOB(i, j, k)] = -9.0 * Qy_i;
        Qz[INDEX_GLOB(i, j, k)] = -9.0 * Qz_i;
    }
}

void collide(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    double t4, t5, t6, t7, t8;
    double cx_bar, cy_bar, cz_bar;
    double u_i, v_i, w_i;
    double Fx_i, Fy_i, Fz_i;
    double Qx_i, Qy_i, Qz_i;
    double rho_i, pressure_i, rho_RED_i, rho_BLUE_i;
    double rho_N_i, nu;
    double omega;
    double G_i, nx_i, ny_i, nz_i;
    double prefac, px, py, pz;

    int i_start = params->i_start;
    int i_end = params->i_end;
    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    double nu_RED = params->nu_RED;
    double nu_BLUE = params->nu_BLUE;

    double beta = params->beta;

    double *rho_comp = comp_fields->rho_comp;

    double *rho = glob_fields->rho;
    double *pressure = glob_fields->pressure;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;
    double *Fx = glob_fields->Fx;
    double *Fy = glob_fields->Fy;
    double *Fz = glob_fields->Fz;
    double *Qx = glob_fields->Qx;
    double *Qy = glob_fields->Qy;
    double *Qz = glob_fields->Qz;
    double *rho_N = glob_fields->rho_N;
    double *nx = glob_fields->nx;
    double *ny = glob_fields->ny;
    double *nz = glob_fields->nz;
    double *G_norm = glob_fields->G_norm;

    double *t_star = dists->t_star;
    double *m_star = dists->m_star;
    double *f_star = dists->f_star;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX_GLOB(i, j, k)];
        pressure_i = pressure[INDEX_GLOB(i, j, k)];
        rho_N_i = rho_N[INDEX_GLOB(i, j, k)];

        rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
        rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];

        u_i = u[INDEX_GLOB(i, j, k)];
        v_i = v[INDEX_GLOB(i, j, k)];
        w_i = w[INDEX_GLOB(i, j, k)];

        Fx_i = Fx[INDEX_GLOB(i, j, k)];
        Fy_i = Fy[INDEX_GLOB(i, j, k)];
        Fz_i = Fz[INDEX_GLOB(i, j, k)];

        Qx_i = Qx[INDEX_GLOB(i, j, k)];
        Qy_i = Qy[INDEX_GLOB(i, j, k)];
        Qz_i = Qz[INDEX_GLOB(i, j, k)];

        G_i = G_norm[INDEX_GLOB(i, j, k)];

        // Saito 2023
        nu = 0.5 * (rho_N_i + 1.0) * nu_RED + 0.5 * (rho_N_i - 1.0) * nu_BLUE;
        omega = 1.0 / (rho_i / pressure_i * nu + 0.5);

        t4 = 0.0;
        t5 = 0.0;
        t6 = 0.0;
        t7 = 0.0;
        t8 = 0.0;

        for (int p = 0; p < NP; p++)
        {
            cx_bar = (double)cx[p] - u_i;
            cy_bar = (double)cy[p] - v_i;
            cz_bar = (double)cz[p] - w_i;

            t4 += f(p) * cx_bar * cy_bar;
            t5 += f(p) * cx_bar * cz_bar;
            t6 += f(p) * cy_bar * cz_bar;
            t7 += f(p) * (cx_bar * cx_bar - cy_bar * cy_bar);
            t8 += f(p) * (cx_bar * cx_bar - cz_bar * cz_bar);
        }

        t_star[0] = rho_i;
        t_star[1] = Fx_i / 2.0;
        t_star[2] = Fy_i / 2.0;
        t_star[3] = Fz_i / 2.0;
        t_star[4] = t4 * (1.0 - omega);
        t_star[5] = t5 * (1.0 - omega);
        t_star[6] = t6 * (1.0 - omega);
        t_star[7] = -t7 * (omega - 1.0) - (Qx_i - Qy_i) * (omega - 2.0) / 2.0;
        t_star[8] = -t8 * (omega - 1.0) - (Qx_i - Qz_i) * (omega - 2.0) / 2.0;
        t_star[9] = Qx_i / 2.0 + Qy_i / 2.0 + Qz_i / 2.0 + 3.0 * pressure_i;
        t_star[10] = Fx_i / 6.0;
        t_star[11] = Fx_i / 6.0;
        t_star[12] = Fy_i / 6.0;
        t_star[13] = Fy_i / 6.0;
        t_star[14] = Fz_i / 6.0;
        t_star[15] = Fz_i / 6.0;
        t_star[16] = 0.0;
        t_star[17] = pressure_i / 3.0;
        t_star[18] = pressure_i / 3.0;
        t_star[19] = pressure_i / 3.0;
        t_star[20] = 0.0;
        t_star[21] = 0.0;
        t_star[22] = 0.0;
        t_star[23] = Fx_i / 18.0;
        t_star[24] = Fy_i / 18.0;
        t_star[25] = Fz_i / 18.0;
        t_star[26] = pressure_i / 9.0;

        central_to_raw(t_star, m_star, u_i, v_i, w_i);

        raw_to_pop(m_star, f_star);

        // Mixing
        for (int p = 0; p < NP; p++)
        {
            f2[INDEX_F(i, j, k, p, RED)] = rho_RED_i / rho_i * f_star[p];
            f2[INDEX_F(i, j, k, p, BLUE)] = rho_BLUE_i / rho_i * f_star[p];
        }

        if (G_i > 1e-15)
        {
            nx_i = nx[INDEX_GLOB(i, j, k)];
            ny_i = ny[INDEX_GLOB(i, j, k)];
            nz_i = nz[INDEX_GLOB(i, j, k)];

            prefac = beta * rho_RED_i * rho_BLUE_i / (rho_i * rho_i) * pressure_i;
            px = prefac * nx_i;
            py = prefac * ny_i;
            pz = prefac * nz_i;

            f2[INDEX_F(i, j, k, 1, RED)] += px / 2.0;
            f2[INDEX_F(i, j, k, 2, RED)] += -px / 2.0;
            f2[INDEX_F(i, j, k, 3, RED)] += py / 2.0;
            f2[INDEX_F(i, j, k, 4, RED)] += -py / 2.0;
            f2[INDEX_F(i, j, k, 5, RED)] += pz / 2.0;
            f2[INDEX_F(i, j, k, 6, RED)] += -pz / 2.0;

            f2[INDEX_F(i, j, k, 1, BLUE)] -= px / 2.0;
            f2[INDEX_F(i, j, k, 2, BLUE)] -= -px / 2.0;
            f2[INDEX_F(i, j, k, 3, BLUE)] -= py / 2.0;
            f2[INDEX_F(i, j, k, 4, BLUE)] -= -py / 2.0;
            f2[INDEX_F(i, j, k, 5, BLUE)] -= pz / 2.0;
            f2[INDEX_F(i, j, k, 6, BLUE)] -= -pz / 2.0;
        }
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

double extrapolate_wall_q(int i, int j, int k, int alpha, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    Stencil *stencil = sim->stencil;
    GlobalFieldBag *glob_fields = sim->glob_fields;

    int ic, jc, kc;
    double sum_q, sum_wp;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;

    int NP = stencil->NP;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    double *rho = glob_fields->rho;
    double *pressure = glob_fields->pressure;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;
    double *u_vec[3] = {u, v, w};

    sum_q = 0.0;
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

        sum_q += wp[p] * (pressure[INDEX_GLOB(ic, jc, kc)] - rho[INDEX_GLOB(ic, jc, kc)] / 3.0) * u_vec[alpha][INDEX_GLOB(ic, jc, kc)];
        sum_wp += wp[p];
    }

    if (sum_wp == 0.0)
        return 0.0;
    else
        return sum_q / sum_wp;
}