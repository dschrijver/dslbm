#include <math.h>
#include <string.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"

void compute_equilibrium(const double rho, const double u, const double v, const double w, const double P, double *restrict const feq, SimulationBag *sim)
{
    UNPACK_BAGS

    WRITE_DIST(meq)

    const double u2 = u * u;
    const double v2 = v * v;
    const double w2 = w * w;

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

void compute_equilibrium_comp(const double rho_comp, const double rho_tot, const double u, const double v, const double w, const double P_tot, double *restrict const feq, SimulationBag *sim)
{
    UNPACK_BAGS

    WRITE_DIST(meq)

    const double u2 = u * u;
    const double v2 = v * v;
    const double w2 = w * w;

    const double rho = rho_comp;
    const double P = rho_comp / rho_tot * P_tot;

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
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho_N)

    WRITE_FIELD(G_norm)
    WRITE_FIELD(Gx)
    WRITE_FIELD(Gy)
    WRITE_FIELD(Gz)
    WRITE_FIELD(nx)
    WRITE_FIELD(ny)
    WRITE_FIELD(nz)

    const int *restrict const flag = fields->flag;

#ifdef THETA_C
    const double theta = THETA_C / 360.0 * 2.0 * DS_PI;
#endif

    FOR_DOMAIN
    {
        // Compute Color Gradient
        double Gx_i = 0.0;
        double Gy_i = 0.0;
        double Gz_i = 0.0;
        double nhat_x = 0.0;
        double nhat_y = 0.0;
        double nhat_z = 0.0;
        for (int p = 0; p < NP; p++)
        {
            int ic = i + cx[p];
            int jc = j + cy[p];
            int kc = k + cz[p];

#ifdef YPERIODIC
            jc = mod(jc, NY);
#endif
#ifdef ZPERIODIC
            kc = mod(kc, NZ);
#endif

            double rho_N_local;

            if (flag[INDEX_FLAG(ic, jc, kc)] > 0)
            {
                rho_N_local = extrapolate_wall_rho_N(ic, jc, kc, sim);
                nhat_x += -wp[p] * (double)cx[p];
                nhat_y += -wp[p] * (double)cy[p];
                nhat_z += -wp[p] * (double)cz[p];
                goto skip;
            }

            rho_N_local = rho_N[INDEX(ic, jc, kc)];

        skip:

            Gx_i += wp[p] * rho_N_local * (double)cx[p];
            Gy_i += wp[p] * rho_N_local * (double)cy[p];
            Gz_i += wp[p] * rho_N_local * (double)cz[p];
        }
        Gx_i *= 3.0;
        Gy_i *= 3.0;
        Gz_i *= 3.0;

        double G_i = sqrt(Gx_i * Gx_i + Gy_i * Gy_i + Gz_i * Gz_i);

#ifdef THETA_C
        const double nhat_norm = sqrt(nhat_x * nhat_x + nhat_y * nhat_y + nhat_z * nhat_z);

        if ((nhat_norm > 0.0) && (G_i > 1e-15))
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

            double Gn = nhat_x * Gx_i + nhat_y * Gy_i + nhat_z * Gz_i;
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

        const int idx = INDEX(i, j, k);

        G_norm[idx] = G_i;
        Gx[idx] = Gx_i;
        Gy[idx] = Gy_i;
        Gz[idx] = Gz_i;
        nx[idx] = Gx_i / (G_i + 1e-15);
        ny[idx] = Gy_i / (G_i + 1e-15);
        nz[idx] = Gz_i / (G_i + 1e-15);
    }
}

void compute_Q_corrections(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho)
    READ_FIELD(pressure)
    READ_FIELD(u)
    READ_FIELD(v)
    READ_FIELD(w)

    WRITE_FIELD(Qx)
    WRITE_FIELD(Qy)
    WRITE_FIELD(Qz)

    const int *restrict const flag = fields->flag;

    FOR_DOMAIN
    {
        double Qx_i = 0.0;
        double Qy_i = 0.0;
        double Qz_i = 0.0;

        for (int p = 0; p < NP; p++)
        {
            int ic = i + cx[p];
            int jc = j + cy[p];
            int kc = k + cz[p];

#ifdef YPERIODIC
            jc = mod(jc, NY);
#endif
#ifdef ZPERIODIC
            kc = mod(kc, NZ);
#endif

            double qx;
            double qy;
            double qz;

            if (flag[INDEX_FLAG(ic, jc, kc)] > 0)
            {
                qx = extrapolate_wall_q(ic, jc, kc, 0, sim);
                qy = extrapolate_wall_q(ic, jc, kc, 1, sim);
                qz = extrapolate_wall_q(ic, jc, kc, 2, sim);
                goto skip;
            }

            const int idxl = INDEX(ic, jc, kc);

            qx = (pressure[idxl] - rho[idxl] / 3.0) * u[idxl];
            qy = (pressure[idxl] - rho[idxl] / 3.0) * v[idxl];
            qz = (pressure[idxl] - rho[idxl] / 3.0) * w[idxl];

        skip:

            Qx_i += wp[p] * qx * (double)cx[p];
            Qy_i += wp[p] * qy * (double)cy[p];
            Qz_i += wp[p] * qz * (double)cz[p];
        }

        const int idx = INDEX(i, j, k);

        Qx[idx] = -9.0 * Qx_i;
        Qy[idx] = -9.0 * Qy_i;
        Qz[idx] = -9.0 * Qz_i;
    }
}

void collide(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho)
    READ_FIELD(rho_RED)
    READ_FIELD(rho_BLUE)
    READ_FIELD(pressure)
    READ_FIELD(u)
    READ_FIELD(v)
    READ_FIELD(w)
    READ_FIELD(Fx)
    READ_FIELD(Fy)
    READ_FIELD(Fz)
    READ_FIELD(Qx)
    READ_FIELD(Qy)
    READ_FIELD(Qz)
    READ_FIELD(rho_N)
    READ_FIELD(nx)
    READ_FIELD(ny)
    READ_FIELD(nz)
    READ_FIELD(G_norm)

    WRITE_DIST(m_star)
    WRITE_DIST(t_star)
    WRITE_DIST(f_star)
    WRITE_DIST(f1_RED)
    WRITE_DIST(f2_RED)
    WRITE_DIST(f1_BLUE)
    WRITE_DIST(f2_BLUE)

    PARAM(mu_RED)
    PARAM(mu_BLUE)
    PARAM(rho_0_BLUE)
    PARAM(cs2_BLUE)
    PARAM(sigma)
    PARAM(beta)

    const double p_0 = cs2_BLUE * rho_0_BLUE;
    const double omega_RED = 1.0 / (mu_RED / p_0 + 0.5);
    const double omega_BLUE = 1.0 / (mu_BLUE / p_0 + 0.5);

    FOR_DOMAIN
    {
        const int idx = INDEX(i, j, k);

        const double rho_i = rho[idx];
        const double pressure_i = pressure[idx];
        const double rho_N_i = rho_N[idx];

        const double rho_RED_i = rho_RED[idx];
        const double rho_BLUE_i = rho_BLUE[idx];

        const double u_i = u[idx];
        const double v_i = v[idx];
        const double w_i = w[idx];

        const double Fx_i = Fx[idx];
        const double Fy_i = Fy[idx];
        const double Fz_i = Fz[idx];

        const double Qx_i = Qx[idx];
        const double Qy_i = Qy[idx];
        const double Qz_i = Qz[idx];

        const double G_i = G_norm[idx];

        // const double omega = interpolate_omega_wen(rho_N_i, params);
        // const double omega = interpolate_omega_ba(rho_N_i, params);
        // const double nu = interpolate_nu_saito(rho_N_i, params);
        // const double omega = 1.0 / (rho_i / pressure_i * nu + 0.5);
        const double omega = 0.5*(1.0 + rho_N_i)*omega_RED + 0.5*(1.0 - rho_N_i)*omega_BLUE;

        double t4 = 0.0;
        double t5 = 0.0;
        double t6 = 0.0;
        double t7 = 0.0;
        double t8 = 0.0;

        for (int p = 0; p < NP; p++)
        {
            double cx_bar = (double)cx[p] - u_i;
            double cy_bar = (double)cy[p] - v_i;
            double cz_bar = (double)cz[p] - w_i;

            const int idxf = INDEX_F(i, j, k, p);
            double f = f1_RED[idxf] + f1_BLUE[idxf];

            t4 += f * cx_bar * cy_bar;
            t5 += f * cx_bar * cz_bar;
            t6 += f * cy_bar * cz_bar;
            t7 += f * (cx_bar * cx_bar - cy_bar * cy_bar);
            t8 += f * (cx_bar * cx_bar - cz_bar * cz_bar);
        }

        t_star[0] = rho_i;
        t_star[1] = Fx_i / 2.0;
        t_star[2] = Fy_i / 2.0;
        t_star[3] = Fz_i / 2.0;
        t_star[4] = t4 * (1.0 - omega);
        t_star[5] = t5 * (1.0 - omega);
        t_star[6] = t6 * (1.0 - omega);
        t_star[7] = t7 * (1.0 - omega) - (Qx_i - Qy_i) * (omega - 2.0) / 2.0;
        t_star[8] = t8 * (1.0 - omega) - (Qx_i - Qz_i) * (omega - 2.0) / 2.0;
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
            const int idxf = INDEX_F(i, j, k, p);

            f2_RED[idxf] = rho_RED_i / rho_i * f_star[p];
            f2_BLUE[idxf] = rho_BLUE_i / rho_i * f_star[p];
        }

        if (G_i > 1e-15)
        {
            const double nx_i = nx[idx];
            const double ny_i = ny[idx];
            const double nz_i = nz[idx];

            const double prefac = beta * rho_RED_i * rho_BLUE_i / (rho_i * rho_i) * pressure_i;
            const double px = prefac * nx_i;
            const double py = prefac * ny_i;
            const double pz = prefac * nz_i;

            f2_RED[INDEX_F(i, j, k, 1)] += px / 2.0;
            f2_RED[INDEX_F(i, j, k, 2)] += -px / 2.0;
            f2_RED[INDEX_F(i, j, k, 3)] += py / 2.0;
            f2_RED[INDEX_F(i, j, k, 4)] += -py / 2.0;
            f2_RED[INDEX_F(i, j, k, 5)] += pz / 2.0;
            f2_RED[INDEX_F(i, j, k, 6)] += -pz / 2.0;

            f2_BLUE[INDEX_F(i, j, k, 1)] -= px / 2.0;
            f2_BLUE[INDEX_F(i, j, k, 2)] -= -px / 2.0;
            f2_BLUE[INDEX_F(i, j, k, 3)] -= py / 2.0;
            f2_BLUE[INDEX_F(i, j, k, 4)] -= -py / 2.0;
            f2_BLUE[INDEX_F(i, j, k, 5)] -= pz / 2.0;
            f2_BLUE[INDEX_F(i, j, k, 6)] -= -pz / 2.0;
        }

        // if (G_i > 1e-15)
        // {
        //     const double nx_i = nx[idx];
        //     const double ny_i = ny[idx];
        //     const double nz_i = nz[idx];

        //     for (int p = 1; p < NP; p++)
        //     {
        //         const double feq_i = 3.0*wp[p]*pressure_i;
        //         const double nc = nx_i * (double)cx[p] + ny_i * (double)cy[p] + nz_i * (double)cz[p];
        //         const double cos_phi = nc / sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);

        //         const double mom_exchange = beta * rho_RED_i * rho_BLUE_i / (rho_i * rho_i) * cos_phi * feq_i;
        //         f2_RED[INDEX_F(i, j, k, p)] += mom_exchange;
        //         f2_BLUE[INDEX_F(i, j, k, p)] -= mom_exchange;
        //     }
        // }
    }
}

double extrapolate_wall_rho_N(const int i, const int j, const int k, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho_N)

    double sum_rho_N = 0.0;
    double sum_wp = 0.0;

    for (int p = 1; p < NP; p++)
    {
        int ic = i + cx[p];
        int jc = j + cy[p];
        int kc = k + cz[p];

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

        sum_rho_N += wp[p] * rho_N[INDEX(ic, jc, kc)];
        sum_wp += wp[p];
    }

    if (sum_wp == 0.0)
        return 0.0;
    else
        return sum_rho_N / sum_wp;
}

double extrapolate_wall_q(const int i, const int j, const int k, const int alpha, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho)
    READ_FIELD(pressure)

    double *restrict u = fields->u;
    double *restrict v = fields->v;
    double *restrict w = fields->w;

    double *u_vec[3] = {u, v, w};

    double sum_q = 0.0;
    double sum_wp = 0.0;

    for (int p = 1; p < NP; p++)
    {
        int ic = i + cx[p];
        int jc = j + cy[p];
        int kc = k + cz[p];

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

        sum_q += wp[p] * (pressure[INDEX(ic, jc, kc)] - rho[INDEX(ic, jc, kc)] / 3.0) * u_vec[alpha][INDEX(ic, jc, kc)];
        sum_wp += wp[p];
    }

    if (sum_wp == 0.0)
        return 0.0;
    else
        return sum_q / sum_wp;
}