#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 400;
    params->NSTORE = 400;
    params->NLOG = 5;
    params->NX = 200;
    params->NY = 200;
    params->NZ = 1;

    // Initial denities
    params->rho_0_a = 20.0;
    params->rho_0_p = 20.0;

    // Kinematic viscosities
    params->eta_a = 10.0 / 3.0;
    params->eta_p = 10.0 / 3.0;

    params->kappa = 5.0;
    params->a = 0.0025;
    params->gamma = 4.0;

    params->S_nem = 1.0;
    params->C = 0.1;
    params->K_LC = 0.15;
    params->lambda = 0.0;
    params->zeta = 0.1;
    params->Gamma_Q = 0.1;

    params->Fx_ext_a = 0.0;
    params->Fy_ext_a = 0.0;
    params->Fz_ext_a = 0.0;
    params->Fx_ext_p = 0.0;
    params->Fy_ext_p = 0.0;
    params->Fz_ext_p = 0.0;

    params->N_active_substeps = 5;
}

#endif