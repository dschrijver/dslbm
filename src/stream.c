#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/stream.h"
#include "../include/wetnode.h"

void stream_distributions(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_DIST(f2_RED)
    READ_DIST(f2_BLUE)

    WRITE_DIST(f1_RED)
    WRITE_DIST(f1_BLUE)

    const int * restrict const flag = fields->flag;

    FOR_DOMAIN
    {
        for (int p = 0; p < NP; p++)
        {
            int ic = i - cx[p];
            int jc = j - cy[p];
            int kc = k - cz[p];

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
                const int p_bb = p_bounceback[p];

                f1_RED[INDEX_F(i, j, k, p)] = f2_RED[INDEX_F(i, j, k, p_bb)];
                f1_BLUE[INDEX_F(i, j, k, p)] = f2_BLUE[INDEX_F(i, j, k, p_bb)];
                continue;
            }

            f1_RED[INDEX_F(i, j, k, p)] = f2_RED[INDEX_F(ic, jc, kc, p)];
            f1_BLUE[INDEX_F(i, j, k, p)] = f2_BLUE[INDEX_F(ic, jc, kc, p)];
        }
    }
}