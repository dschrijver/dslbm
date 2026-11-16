#include "../include/datatypes.h"
#include "../definitions.h"

void stream_distributions(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_DIST(f2_a)
    READ_DIST(f2_p)

    WRITE_DIST(f1_a)
    WRITE_DIST(f1_p)

    FOR_DOMAIN
    {
        for (int p = 0; p < NP; p++)
        {
            int ic = i - cx[p];
            int jc = j - cy[p];
            int kc = k - cz[p];

            f1_a[INDEX_F(i, j, k, p)] = f2_a[INDEX_F(ic, jc, kc, p)];
            f1_p[INDEX_F(i, j, k, p)] = f2_p[INDEX_F(ic, jc, kc, p)];
        }
    }
}