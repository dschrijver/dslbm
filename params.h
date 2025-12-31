#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 1000;
    params->NSTORE = 100;
    params->NX = 60;
    params->NY = 60;
    params->NZ = 60;

    // Relaxation times
    params->tau_RED = 1.0;
    params->tau_BLUE = 0.625;

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

    params->sigma = 1e-3;
    params->beta = 0.7;
    params->alpha_BLUE = 0.2;
    double gamma = params->rho_0_RED / params->rho_0_BLUE;
    params->alpha_RED = 1.0 - (1.0 - params->alpha_BLUE) / gamma;

    params->Fx_ext = 0.0;
    params->Fy_ext = 0.0;
    params->Fz_ext = 0.0;
}

#endif