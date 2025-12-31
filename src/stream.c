#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/stream.h"

void stream_distributions(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    int ic, jc, kc;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    FOR_DOMAIN
    {
        for (int p = 0; p < NP; p++)
        {
            ic = i - cx[p];
            jc = j - cy[p];
            kc = k - cz[p];

#ifdef YPERIODIC
            jc = mod(jc, NY);
#endif
#ifdef ZPERIODIC
            kc = mod(kc, NZ);
#endif

#ifdef BOTTOM_BOUNCEBACK
            if (jc < 0)
            {
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], RED)];
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], BLUE)];
                continue;
            }
#endif

#ifdef TOP_BOUNCEBACK
            if (jc > NY-1)
            {
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], RED)];
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(i, j, k, stencil->p_bounceback[p], BLUE)];
                continue;
            }
#endif

            f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(ic, jc, kc, p, RED)];
            f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(ic, jc, kc, p, BLUE)];
        }
    }
}