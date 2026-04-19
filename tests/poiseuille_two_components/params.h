#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 1000000;
    params->NSTORE = 10000;
    params->NLOG = 1000;
    params->NX = 100;
    params->NY = 1;
    params->NZ = 1;

    // Initial denities
    params->rho_0_RED   = 8.0;
    params->rho_0_BLUE  = 0.008;

    // Kinematic viscosities
    params->nu_BLUE  = 1.0/6.0;
    params->nu_RED = params->nu_BLUE/25.0;

    // Speed of sound parameters
    double gamma = params->rho_0_RED / params->rho_0_BLUE;
    params->alpha_BLUE  = 0.2;
    params->alpha_RED   = 1.0 - (1.0 - params->alpha_BLUE) / gamma;

    // Surface tension
    params->sigma = 0.0;

    // Recoloring parameter
    params->beta = 0.7;

    params->gx = 0.0;
    params->gy = 1.5e-8;
    params->gz = 0.0;
}

#endif