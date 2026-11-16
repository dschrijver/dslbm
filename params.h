#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 6000;
    params->NSTORE = 60;
    params->NLOG = 6;
    params->NX = 1;
    params->NY = 33;
    params->NZ = 1;

    // Initial denities
    params->rho_0_a = 0.5;
    params->rho_0_p = 0.5;

    // Dynamic viscosities
    params->eta_a = 0.075;
    params->eta_p = 0.075;

    params->kappa = 5.0;
    params->a = 0.0025;
    params->gamma = 4.0;

    params->S_nem = 1.0;
    params->C = 0.1;
    params->K_LC = 0.15;
    params->lambda = 0.0;
    params->zeta = 0.1;
    params->Gamma_Q = 0.1;

    params->Fx_ext_a = 0.5e-5;
    params->Fy_ext_a = 0.0;
    params->Fz_ext_a = 0.0;
    params->Fx_ext_p = 0.5e-5;
    params->Fy_ext_p = 0.0;
    params->Fz_ext_p = 0.0;

    params->N_active_substeps = 5;
}

#endif