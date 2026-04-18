#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/stream.h"
#include "../include/wetnode.h"

void stream_distributions(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    int ic, jc, kc, p_bb;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int *flag = sim->glob_fields->flag;

    int i_start = params->i_start;
    int i_end = params->i_end;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    int *p_bounceback = stencil->p_bounceback;

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

            if (flag[INDEX_FLAG(ic, jc, kc)] == WETNODE)
                continue;

            if (flag[INDEX_FLAG(ic, jc, kc)] == BOUNDARY)
            {
                p_bb = p_bounceback[p];
                f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(i, j, k, p_bb, RED)];
                f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(i, j, k, p_bb, BLUE)];
                continue;
            }

            f1[INDEX_F(i, j, k, p, RED)] = f2[INDEX_F(ic, jc, kc, p, RED)];
            f1[INDEX_F(i, j, k, p, BLUE)] = f2[INDEX_F(ic, jc, kc, p, BLUE)];
        }
    }
}