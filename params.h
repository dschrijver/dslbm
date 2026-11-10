#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 100000;
    params->NSTORE = 10000;
    params->NLOG = 1000;
    params->NX = 100;
    params->NY = 2;
    params->NZ = 2;
}

#endif