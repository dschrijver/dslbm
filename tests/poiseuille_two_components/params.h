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
    params->rho_0_RED   = 1.0;
    params->rho_0_BLUE  = 0.001;

    // Kinematic viscosities
    params->mu_RED  = 0.1;
    params->mu_BLUE = 0.0025;

    // Speeds of sound squared
    params->cs2_RED = 1.0/3.0/1000.0;
    params->cs2_BLUE = 1.0/3.0;

    // Surface tension
    params->sigma = 1e-4;

    // Recoloring parameter
    params->beta = 0.7;

    // Gravitational accelerations
    params->gx = 0.0;
    params->gy = 0.0;
    params->gz = 0.0;

    // Constant body forces
    params->Fb_x = 0.0;
    params->Fb_y = 1.5e-8;
    params->Fb_z = 0.0;
}

#endif