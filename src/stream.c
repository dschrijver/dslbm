#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/stream.h"
#include "../include/wetnode.h"

void stream_distributions(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_DIST(f2)
    READ_DIST(g2)

    WRITE_DIST(f1)
    WRITE_DIST(g1)

    FOR_DOMAIN
    {
        for (int p = 0; p < NP; p++)
        {
            int ic = i - cx[p];
            int jc = j - cy[p];
            int kc = k - cz[p];

            f1[INDEX_F(i, j, k, p)] = f2[INDEX_F(ic, jc, kc, p)];
            g1[INDEX_F(i, j, k, p)] = g2[INDEX_F(ic, jc, kc, p)];
        }
    }
}