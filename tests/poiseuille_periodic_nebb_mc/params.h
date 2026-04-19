#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 100000;
    params->NSTORE = 1000;
    params->NLOG = 100;
    params->NX = 2;
    params->NY = 16;
    params->NZ = 1;

    // Initial denities
    params->rho_0_RED   = 1.0;
    params->rho_0_BLUE  = 0.001;

    // Kinematic viscosities
    params->nu_RED  = 0.1;
    params->nu_BLUE = 0.1;

    // Speed of sound parameters
    double gamma = params->rho_0_RED / params->rho_0_BLUE;
    params->alpha_BLUE  = 0.2;
    params->alpha_RED   = 1.0 - (1.0 - params->alpha_BLUE) / gamma;

    // Surface tension
    params->sigma = 1e-4;

    // Recoloring parameter
    params->beta = 0.2;

    params->gx = 1e-5;
    params->gy = 0.0;
    params->gz = 0.0;
}

#endif