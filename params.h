#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME   = 10000;
    params->NSTORE  = 100;
    params->NLOG    = 10;
    params->NX      = 61;
    params->NY      = 61;
    params->NZ      = 61;

    // Initial denities
    params->rho_0_RED   = 1000.0;
    params->rho_0_BLUE  = 1.0;

    // Kinematic viscosities
    params->nu_RED  = 0.00001;
    params->nu_BLUE = 0.0001;

    // Speed of sound parameters
    double gamma = params->rho_0_RED / params->rho_0_BLUE;
    params->alpha_BLUE  = 0.2;
    params->alpha_RED   = 1.0 - (1.0 - params->alpha_BLUE) / gamma;

    // Surface tension
    params->sigma = 1e-4;

    // Recoloring parameter
    params->beta = 0.7;

    // Gravitational accelerations
    params->gx = 0.0;
    params->gy = 0.0;
    params->gz = 0.0;
}

#endif