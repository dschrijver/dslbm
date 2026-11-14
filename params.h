#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 10000;
    params->NSTORE = 100;
    params->NLOG = 1;
    params->NX = 64;
    params->NY = 64;
    params->NZ = 64;

    // Initial denities
    params->rho_0_a   = 20.0;
    params->rho_0_p  = 20.0;

    // Kinematic viscosities
    const double mu_a = 10.0 / 3.0;
    const double mu_p = 10.0 / 3.0;

    params->nu_a  = mu_a / params->rho_0_a;
    params->nu_p = mu_p / params->rho_0_p;

    params->kappa = 5.0;
    params->a = 0.0025;
    params->gamma = 4.0;

    params->S_nem = 1.0;
    params->C = 0.1;
    params->K_LC = 0.15;
    params->lambda = 0.0;
    params->zeta = 0.1;
    params->Gamma_Q = 0.1;

    params->N_active_substeps = 5;
}

#endif