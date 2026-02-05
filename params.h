#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 25000;
    params->NSTORE = 250;
    params->NLOG = 250;
    params->NX = 201;
    params->NY = 101;
    params->NZ = 1;

    // Relaxation times
    params->tau_RED = 1.0;
    params->tau_BLUE = 0.505;

    // MRT parameters
    params->s_rho_RED = 1.0;
    params->s_rho_BLUE = 1.0;
    params->s_j_RED = 1.0;
    params->s_j_BLUE = 1.0;
    params->s_e_RED = 1.1;
    params->s_e_BLUE = 1.1;
    params->s_eps_RED = 1.1;
    params->s_eps_BLUE = 1.1;
    params->s_q_RED = 1.1;
    params->s_q_BLUE = 1.1;
    params->s_pi_RED = 1.1;
    params->s_pi_BLUE = 1.1;
    params->s_m_RED = 1.2;
    params->s_m_BLUE = 1.2;

    // Starting densities
    params->rho_0_RED = 1.0;
    params->rho_0_BLUE = 1.0;

    params->sigma = 1e-4;
    params->beta = 0.7;
    params->alpha_BLUE = 0.2;
    // params->alpha_BLUE = 1.0/3.0;
    double gamma = params->rho_0_RED / params->rho_0_BLUE;
    params->alpha_RED = 1.0 - (1.0 - params->alpha_BLUE) / gamma;

    params->gx = 0.0;
    params->gy = 0.0;
    params->gz = 0.0;

    params->G_SC = 4.0;
}

#endif