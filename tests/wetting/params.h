#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME   = 10;
    params->NSTORE  = 1;
    params->NLOG    = 1;
    params->NX      = 100;
    params->NY      = 100;
    params->NZ      = 1;

    // Initial denities
    params->rho_0_RED   = 1.0;
    params->rho_0_BLUE  = 1.0;

    // Kinematic viscosities
    params->nu_RED  = 1.0/6.0;
    params->nu_BLUE = 1.0/6.0;

    // Speed of sound parameters
    double gamma = params->rho_0_RED / params->rho_0_BLUE;
    params->alpha_BLUE  = 8.0/27.0;
    params->alpha_RED   = 1.0 - (1.0 - params->alpha_BLUE) / gamma;

    // Surface tension
    params->sigma = 1e-2;

    // Recoloring parameter
    params->beta = 0.7;

    // Gravitational accelerations
    params->gx = 0.0;
    params->gy = 0.0;
    params->gz = 0.0;
}

#endif