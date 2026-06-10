#ifndef COLLIDE_H
#define COLLIDE_H

#include "datatypes.h"

void compute_equilibrium(const double rho, const double u, const double v, const double w, const double P, double * restrict const feq, SimulationBag *sim);
void compute_equilibrium_comp(const double rho_comp, const double rho_tot, const double u, const double v, const double w, const double P_tot, double * restrict const feq, SimulationBag *sim);
void evaluate_color_gradients(SimulationBag *sim);
void compute_Q_corrections(SimulationBag *sim);

void collide(SimulationBag *sim);

double extrapolate_wall_rho_N(const int i, const int j, const int k, SimulationBag *sim);
double extrapolate_wall_q(const int i, const int j, const int k, const int alpha, SimulationBag *sim);

static inline double interpolate_omega_wen(const double phi, ParamBag *params)
{
    PARAM(mu_RED)
    PARAM(mu_BLUE)
    PARAM(cs2_BLUE)
    PARAM(rho_0_BLUE)

    const double p_0 = cs2_BLUE * rho_0_BLUE;

    const double tau_RED = mu_RED / p_0 + 0.5;
    const double tau_BLUE = mu_BLUE / p_0 + 0.5;

    // Wen 2019, 10.1103/PhysRevE.100.023301
    const double delta = 0.98;

    // Grunau 1993, 10.1063/1.858769
    const double alpha = 2.0 * tau_RED * tau_BLUE / (tau_RED + tau_BLUE);
    const double beta = 2.0 * (tau_RED - alpha) / delta;
    const double kappa = -beta / (2.0 * delta);
    const double eta = 2.0 * (alpha - tau_BLUE) / delta;
    const double xi = eta / (2.0 * delta);

    double tau;
    if (phi > delta)
    {
        tau = tau_RED;
    }
    else if (phi > 0)
    {
        tau = alpha + beta * phi + kappa * phi * phi;
    }
    else if (phi >= -delta)
    {
        tau = alpha + eta * phi + xi * phi * phi;
    }
    else
    {
        tau = tau_BLUE;
    }

    return 1.0 / tau;
}

static inline double interpolate_omega_ba(const double phi, ParamBag *params)
{
    PARAM(mu_RED)
    PARAM(mu_BLUE)
    PARAM(cs2_BLUE)
    PARAM(rho_0_BLUE)

    const double p_0 = cs2_BLUE * rho_0_BLUE;

    const double omega_RED = 1.0 / (mu_RED / p_0 + 0.5);
    const double omega_BLUE = 1.0 / (mu_BLUE / p_0 + 0.5);

    // Wen 2019, 10.1103/PhysRevE.100.023301
    const double delta = 0.1;

    // Grunau 1993, 10.1063/1.858769
    const double alpha = 2.0 * omega_RED * omega_BLUE / (omega_RED + omega_BLUE);
    const double beta = 2.0 * (omega_RED - alpha) / delta;
    const double kappa = -beta / (2.0 * delta);
    const double eta = 2.0 * (alpha - omega_BLUE) / delta;
    const double xi = eta / (2.0 * delta);

    double omega;
    if (phi > delta)
    {
        omega = omega_RED;
    }
    else if (phi > 0)
    {
        omega = alpha + beta * phi + kappa * phi * phi;
    }
    else if (phi >= -delta)
    {
        omega = alpha + eta * phi + xi * phi * phi;
    }
    else
    {
        omega = omega_BLUE;
    }

    return omega;
}

static inline double interpolate_nu_saito(const double phi, ParamBag *params)
{
    PARAM(mu_RED)
    PARAM(mu_BLUE)
    PARAM(rho_0_RED)
    PARAM(rho_0_BLUE)

    const double nu_RED = mu_RED / rho_0_RED;
    const double nu_BLUE = mu_BLUE / rho_0_BLUE;

    const double nu = 0.5*(1.0 + phi)*nu_RED + 0.5*(1.0 - phi)*nu_BLUE;
    return nu;
}

static inline void raw_to_pop(const double *restrict raw, double *restrict pop)
{
    pop[0] = raw[0] + raw[17] + raw[18] + raw[19] - raw[26] - raw[9];
    pop[1] = -raw[10] / 2.0 - raw[11] / 2.0 - raw[17] / 2.0 - raw[18] / 2.0 + raw[1] / 2.0 + raw[23] / 2.0 + raw[26] / 2.0 + raw[7] / 6.0 + raw[8] / 6.0 + raw[9] / 6.0;
    pop[2] = raw[10] / 2.0 + raw[11] / 2.0 - raw[17] / 2.0 - raw[18] / 2.0 - raw[1] / 2.0 - raw[23] / 2.0 + raw[26] / 2.0 + raw[7] / 6.0 + raw[8] / 6.0 + raw[9] / 6.0;
    pop[3] = -raw[12] / 2.0 - raw[13] / 2.0 - raw[17] / 2.0 - raw[19] / 2.0 + raw[24] / 2.0 + raw[26] / 2.0 + raw[2] / 2.0 - raw[7] / 3.0 + raw[8] / 6.0 + raw[9] / 6.0;
    pop[4] = raw[12] / 2.0 + raw[13] / 2.0 - raw[17] / 2.0 - raw[19] / 2.0 - raw[24] / 2.0 + raw[26] / 2.0 - raw[2] / 2.0 - raw[7] / 3.0 + raw[8] / 6.0 + raw[9] / 6.0;
    pop[5] = -raw[14] / 2.0 - raw[15] / 2.0 - raw[18] / 2.0 - raw[19] / 2.0 + raw[25] / 2.0 + raw[26] / 2.0 + raw[3] / 2.0 + raw[7] / 6.0 - raw[8] / 3.0 + raw[9] / 6.0;
    pop[6] = raw[14] / 2.0 + raw[15] / 2.0 - raw[18] / 2.0 - raw[19] / 2.0 - raw[25] / 2.0 + raw[26] / 2.0 - raw[3] / 2.0 + raw[7] / 6.0 - raw[8] / 3.0 + raw[9] / 6.0;
    pop[7] = raw[10] / 4.0 + raw[13] / 4.0 + raw[17] / 4.0 - raw[22] / 4.0 - raw[23] / 4.0 - raw[24] / 4.0 - raw[26] / 4.0 + raw[4] / 4.0;
    pop[8] = -raw[10] / 4.0 - raw[13] / 4.0 + raw[17] / 4.0 - raw[22] / 4.0 + raw[23] / 4.0 + raw[24] / 4.0 - raw[26] / 4.0 + raw[4] / 4.0;
    pop[9] = raw[10] / 4.0 - raw[13] / 4.0 + raw[17] / 4.0 + raw[22] / 4.0 - raw[23] / 4.0 + raw[24] / 4.0 - raw[26] / 4.0 - raw[4] / 4.0;
    pop[10] = -raw[10] / 4.0 + raw[13] / 4.0 + raw[17] / 4.0 + raw[22] / 4.0 + raw[23] / 4.0 - raw[24] / 4.0 - raw[26] / 4.0 - raw[4] / 4.0;
    pop[11] = raw[12] / 4.0 + raw[15] / 4.0 + raw[19] / 4.0 - raw[20] / 4.0 - raw[24] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 + raw[6] / 4.0;
    pop[12] = -raw[12] / 4.0 - raw[15] / 4.0 + raw[19] / 4.0 - raw[20] / 4.0 + raw[24] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 + raw[6] / 4.0;
    pop[13] = raw[12] / 4.0 - raw[15] / 4.0 + raw[19] / 4.0 + raw[20] / 4.0 - raw[24] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 - raw[6] / 4.0;
    pop[14] = -raw[12] / 4.0 + raw[15] / 4.0 + raw[19] / 4.0 + raw[20] / 4.0 + raw[24] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 - raw[6] / 4.0;
    pop[15] = raw[11] / 4.0 + raw[14] / 4.0 + raw[18] / 4.0 - raw[21] / 4.0 - raw[23] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 + raw[5] / 4.0;
    pop[16] = -raw[11] / 4.0 - raw[14] / 4.0 + raw[18] / 4.0 - raw[21] / 4.0 + raw[23] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 + raw[5] / 4.0;
    pop[17] = raw[11] / 4.0 - raw[14] / 4.0 + raw[18] / 4.0 + raw[21] / 4.0 - raw[23] / 4.0 + raw[25] / 4.0 - raw[26] / 4.0 - raw[5] / 4.0;
    pop[18] = -raw[11] / 4.0 + raw[14] / 4.0 + raw[18] / 4.0 + raw[21] / 4.0 + raw[23] / 4.0 - raw[25] / 4.0 - raw[26] / 4.0 - raw[5] / 4.0;
    pop[19] = raw[16] / 8.0 + raw[20] / 8.0 + raw[21] / 8.0 + raw[22] / 8.0 + raw[23] / 8.0 + raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
    pop[20] = -raw[16] / 8.0 + raw[20] / 8.0 + raw[21] / 8.0 + raw[22] / 8.0 - raw[23] / 8.0 - raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
    pop[21] = -raw[16] / 8.0 - raw[20] / 8.0 - raw[21] / 8.0 + raw[22] / 8.0 + raw[23] / 8.0 + raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
    pop[22] = raw[16] / 8.0 - raw[20] / 8.0 - raw[21] / 8.0 + raw[22] / 8.0 - raw[23] / 8.0 - raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
    pop[23] = -raw[16] / 8.0 - raw[20] / 8.0 + raw[21] / 8.0 - raw[22] / 8.0 + raw[23] / 8.0 - raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
    pop[24] = raw[16] / 8.0 - raw[20] / 8.0 + raw[21] / 8.0 - raw[22] / 8.0 - raw[23] / 8.0 + raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
    pop[25] = -raw[16] / 8.0 + raw[20] / 8.0 - raw[21] / 8.0 - raw[22] / 8.0 - raw[23] / 8.0 + raw[24] / 8.0 + raw[25] / 8.0 + raw[26] / 8.0;
    pop[26] = raw[16] / 8.0 + raw[20] / 8.0 - raw[21] / 8.0 - raw[22] / 8.0 + raw[23] / 8.0 - raw[24] / 8.0 - raw[25] / 8.0 + raw[26] / 8.0;
}

static inline void central_to_raw(const double *restrict central, double *restrict raw, const double u, const double v, const double w)
{
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;

    raw[0] = central[0];
    raw[1] = central[0] * u + central[1];
    raw[2] = central[0] * v + central[2];
    raw[3] = central[0] * w + central[3];
    raw[4] = central[0] * u * v + central[1] * v + central[2] * u + central[4];
    raw[5] = central[0] * u * w + central[1] * w + central[3] * u + central[5];
    raw[6] = central[0] * v * w + central[2] * w + central[3] * v + central[6];
    raw[7] = central[0] * (u2 - v2) + 2.0 * central[1] * u - 2.0 * central[2] * v + central[7];
    raw[8] = central[0] * (u2 - w2) + 2.0 * central[1] * u - 2.0 * central[3] * w + central[8];
    raw[9] = central[0] * (u2 + v2 + w2) + 2.0 * central[1] * u + 2.0 * central[2] * v + 2.0 * central[3] * w + central[9];
    raw[10] = central[0] * u * v2 + central[10] + central[1] * v2 + 2.0 * central[2] * u * v + 2.0 * central[4] * v - 2.0 * central[7] * u / 3.0 + central[8] * u / 3.0 + central[9] * u / 3.0;
    raw[11] = central[0] * u * w2 + central[11] + central[1] * w2 + 2.0 * central[3] * u * w + 2.0 * central[5] * w + central[7] * u / 3.0 - 2.0 * central[8] * u / 3.0 + central[9] * u / 3.0;
    raw[12] = central[0] * v * w2 + central[12] + central[2] * w2 + 2.0 * central[3] * v * w + 2.0 * central[6] * w + central[7] * v / 3.0 - 2.0 * central[8] * v / 3.0 + central[9] * v / 3.0;
    raw[13] = central[0] * u2 * v + central[13] + 2.0 * central[1] * u * v + central[2] * u2 + 2.0 * central[4] * u + central[7] * v / 3.0 + central[8] * v / 3.0 + central[9] * v / 3.0;
    raw[14] = central[0] * u2 * w + central[14] + 2.0 * central[1] * u * w + central[3] * u2 + 2.0 * central[5] * u + central[7] * w / 3.0 + central[8] * w / 3.0 + central[9] * w / 3.0;
    raw[15] = central[0] * v2 * w + central[15] + 2.0 * central[2] * v * w + central[3] * v2 + 2.0 * central[6] * v - 2.0 * central[7] * w / 3.0 + central[8] * w / 3.0 + central[9] * w / 3.0;
    raw[16] = central[0] * u * v * w + central[16] + central[1] * v * w + central[2] * u * w + central[3] * u * v + central[4] * w + central[5] * v + central[6] * u;
    raw[17] = central[0] * u2 * v2 + 2.0 * central[10] * u + 2.0 * central[13] * v + central[17] + 2.0 * central[1] * u * v2 + 2.0 * central[2] * u2 * v + 4.0 * central[4] * u * v - central[7] * (2.0 * u2 - v2) / 3.0 + central[8] * (u2 + v2) / 3.0 + central[9] * (u2 + v2) / 3.0;
    raw[18] = central[0] * u2 * w2 + 2.0 * central[11] * u + 2.0 * central[14] * w + central[18] + 2.0 * central[1] * u * w2 + 2.0 * central[3] * u2 * w + 4.0 * central[5] * u * w + central[7] * (u2 + w2) / 3.0 - central[8] * (2.0 * u2 - w2) / 3.0 + central[9] * (u2 + w2) / 3.0;
    raw[19] = central[0] * v2 * w2 + 2.0 * central[12] * v + 2.0 * central[15] * w + central[19] + 2.0 * central[2] * v * w2 + 2.0 * central[3] * v2 * w + 4.0 * central[6] * v * w + central[7] * (v2 - 2.0 * w2) / 3.0 - central[8] * (2.0 * v2 - w2) / 3.0 + central[9] * (v2 + w2) / 3.0;
    raw[20] = central[0] * u2 * v * w + central[13] * w + central[14] * v + 2.0 * central[16] * u + 2.0 * central[1] * u * v * w + central[20] + central[2] * u2 * w + central[3] * u2 * v + 2.0 * central[4] * u * w + 2.0 * central[5] * u * v + central[6] * u2 + central[7] * v * w / 3.0 + central[8] * v * w / 3.0 + central[9] * v * w / 3.0;
    raw[21] = central[0] * u * v2 * w + central[10] * w + central[15] * u + 2.0 * central[16] * v + central[1] * v2 * w + central[21] + 2.0 * central[2] * u * v * w + central[3] * u * v2 + 2.0 * central[4] * v * w + central[5] * v2 + 2.0 * central[6] * u * v - 2.0 * central[7] * u * w / 3.0 + central[8] * u * w / 3.0 + central[9] * u * w / 3.0;
    raw[22] = central[0] * u * v * w2 + central[11] * v + central[12] * u + 2.0 * central[16] * w + central[1] * v * w2 + central[22] + central[2] * u * w2 + 2.0 * central[3] * u * v * w + central[4] * w2 + 2.0 * central[5] * v * w + 2.0 * central[6] * u * w + central[7] * u * v / 3.0 - 2.0 * central[8] * u * v / 3.0 + central[9] * u * v / 3.0;
    raw[23] = central[0] * u * v2 * w2 + central[10] * w2 + central[11] * v2 + 2.0 * central[12] * u * v + 2.0 * central[15] * u * w + 4.0 * central[16] * v * w + central[19] * u + central[1] * v2 * w2 + 2.0 * central[21] * w + 2.0 * central[22] * v + central[23] + 2.0 * central[2] * u * v * w2 + 2.0 * central[3] * u * v2 * w + 2.0 * central[4] * v * w2 + 2.0 * central[5] * v2 * w + 4.0 * central[6] * u * v * w + central[7] * u * (v2 - 2.0 * w2) / 3.0 - central[8] * u * (2.0 * v2 - w2) / 3.0 + central[9] * u * (v2 + w2) / 3.0;
    raw[24] = central[0] * u2 * v * w2 + 2.0 * central[11] * u * v + central[12] * u2 + central[13] * w2 + 2.0 * central[14] * v * w + 4.0 * central[16] * u * w + central[18] * v + 2.0 * central[1] * u * v * w2 + 2.0 * central[20] * w + 2.0 * central[22] * u + central[24] + central[2] * u2 * w2 + 2.0 * central[3] * u2 * v * w + 2.0 * central[4] * u * w2 + 4.0 * central[5] * u * v * w + 2.0 * central[6] * u2 * w + central[7] * v * (u2 + w2) / 3.0 - central[8] * v * (2.0 * u2 - w2) / 3.0 + central[9] * v * (u2 + w2) / 3.0;
    raw[25] = central[0] * u2 * v2 * w + 2.0 * central[10] * u * w + 2.0 * central[13] * v * w + central[14] * v2 + central[15] * u2 + 4.0 * central[16] * u * v + central[17] * w + 2.0 * central[1] * u * v2 * w + 2.0 * central[20] * v + 2.0 * central[21] * u + central[25] + 2.0 * central[2] * u2 * v * w + central[3] * u2 * v2 + 4.0 * central[4] * u * v * w + 2.0 * central[5] * u * v2 + 2.0 * central[6] * u2 * v - central[7] * w * (2.0 * u2 - v2) / 3.0 + central[8] * w * (u2 + v2) / 3.0 + central[9] * w * (u2 + v2) / 3.0;
    raw[26] = central[0] * u2 * v2 * w2 + 2.0 * central[10] * u * w2 + 2.0 * central[11] * u * v2 + 2.0 * central[12] * u2 * v + 2.0 * central[13] * v * w2 + 2.0 * central[14] * v2 * w + 2.0 * central[15] * u2 * w + 8.0 * central[16] * u * v * w + central[17] * w2 + central[18] * v2 + central[19] * u2 + 2.0 * central[1] * u * v2 * w2 + 4.0 * central[20] * v * w + 4.0 * central[21] * u * w + 4.0 * central[22] * u * v + 2.0 * central[23] * u + 2.0 * central[24] * v + 2.0 * central[25] * w + central[26] + 2.0 * central[2] * u2 * v * w2 + 2.0 * central[3] * u2 * v2 * w + 4.0 * central[4] * u * v * w2 + 4.0 * central[5] * u * v2 * w + 4.0 * central[6] * u2 * v * w + central[7] * (u2 * v2 - 2.0 * u2 * w2 + v2 * w2) / 3.0 + central[8] * (-2.0 * u2 * v2 + u2 * w2 + v2 * w2) / 3.0 + central[9] * (u2 * v2 + u2 * w2 + v2 * w2) / 3.0;
}

#endif