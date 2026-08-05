#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "datatypes.h"

// fields.c
void evaluate_pressure(SimulationBag *sim);

static inline double physx(int i)
{
#if defined(LEFT_NEBB_VELOCITY) || defined(LEFT_NEBB_PRESSURE)
    return (double)i;
#else
    return (double)i + 0.5;
#endif
}

static inline double physy(int j)
{
#if defined(BOTTOM_NEBB_VELOCITY) || defined(BOTTOM_NEBB_PRESSURE)
    return (double)j;
#else
    return (double)j + 0.5;
#endif
}

static inline double physz(int k)
{
#if defined(BACK_NEBB_VELOCITY) || defined(BACK_NEBB_PRESSURE)
    return (double)k;
#else
    return (double)k + 0.5;
#endif
}

static inline double physlx(ParamBag *params)
{
    double result = (double)params->NX;
#if defined(LEFT_NEBB_VELOCITY) || defined(LEFT_NEBB_PRESSURE)
    result -= 0.5;
#endif

#if defined(RIGHT_NEBB_VELOCITY) || defined(RIGHT_NEBB_PRESSURE)
    result -= 0.5;
#endif
    return result;
}

static inline double physly(ParamBag *params)
{
    double result = (double)params->NY;
#if defined(BOTTOM_NEBB_VELOCITY) || defined(BOTTOM_NEBB_PRESSURE)
    result -= 0.5;
#endif

#if defined(TOP_NEBB_VELOCITY) || defined(TOP_NEBB_PRESSURE)
    result -= 0.5;
#endif
    return result;
}

static inline double physlz(ParamBag *params)
{
    double result = (double)params->NZ;
#if defined(BACK_NEBB_VELOCITY) || defined(BACK_NEBB_PRESSURE)
    result -= 0.5;
#endif

#if defined(FRONT_NEBB_VELOCITY) || defined(FRONT_NEBB_PRESSURE)
    result -= 0.5;
#endif
    return result;
}

#endif