#include "../include/datatypes.h"
#include "../definitions.h"

void bc_distributions(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    WRITE_DIST(f2_a)
    WRITE_DIST(f2_p)

#ifdef LEFT_NOSLIP
    if (i_start == 0)
    {
        for (int j = j_start; j < j_end; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                for (int p = 0; p < NP; p++)
                {
                    if (cx[p] < 0)
                    {
                        const int p_bb = p_bounceback[p];
                        const int jc = j + cy[p];
                        const int kc = k + cz[p];

                        f2_a[INDEX_F(i_start - 1, jc, kc, p_bb)] = f2_a[INDEX_F(i_start, j, k, p)];
                        f2_p[INDEX_F(i_start - 1, jc, kc, p_bb)] = f2_p[INDEX_F(i_start, j, k, p)];
                    }
                }
            }
        }
    }
#endif

#ifdef RIGHT_NOSLIP
    if (i_end == NX)
    {
        for (int j = j_start; j < j_end; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                for (int p = 0; p < NP; p++)
                {
                    if (cx[p] > 0)
                    {
                        const int p_bb = p_bounceback[p];
                        const int jc = j + cy[p];
                        const int kc = k + cz[p];

                        f2_a[INDEX_F(i_end, jc, kc, p_bb)] = f2_a[INDEX_F(i_end - 1, j, k, p)];
                        f2_p[INDEX_F(i_end, jc, kc, p_bb)] = f2_p[INDEX_F(i_end - 1, j, k, p)];
                    }
                }
            }
        }
    }
#endif

#ifdef BOTTOM_NOSLIP
    if (j_start == 0)
    {
        for (int i = i_start; i < i_end; i++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                for (int p = 0; p < NP; p++)
                {
                    if (cy[p] < 0)
                    {
                        const int p_bb = p_bounceback[p];
                        const int ic = i + cx[p];
                        const int kc = k + cz[p];

                        f2_a[INDEX_F(ic, j_start - 1, kc, p_bb)] = f2_a[INDEX_F(i, j_start, k, p)];
                        f2_p[INDEX_F(ic, j_start - 1, kc, p_bb)] = f2_p[INDEX_F(i, j_start, k, p)];
                    }
                }
            }
        }
    }
#endif

#ifdef TOP_NOSLIP
    if (j_end == NY)
    {
        for (int i = i_start; i < i_end; i++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                for (int p = 0; p < NP; p++)
                {
                    if (cy[p] > 0)
                    {
                        const int p_bb = p_bounceback[p];
                        const int ic = i + cx[p];
                        const int kc = k + cz[p];

                        f2_a[INDEX_F(ic, j_end, kc, p_bb)] = f2_a[INDEX_F(i, j_end - 1, k, p)];
                        f2_p[INDEX_F(ic, j_end, kc, p_bb)] = f2_p[INDEX_F(i, j_end - 1, k, p)];
                    }
                }
            }
        }
    }
#endif

#ifdef BACK_NOSLIP
    if (k_start == 0)
    {
        for (int i = i_start; i < i_end; i++)
        {
            for (int j = j_start; j < j_end; j++)
            {
                for (int p = 0; p < NP; p++)
                {
                    if (cz[p] < 0)
                    {
                        const int p_bb = p_bounceback[p];
                        const int ic = i + cx[p];
                        const int jc = j + cy[p];

                        f2_a[INDEX_F(ic, jc, k_start - 1, p_bb)] = f2_a[INDEX_F(i, j, k_start, p)];
                        f2_p[INDEX_F(ic, jc, k_start - 1, p_bb)] = f2_p[INDEX_F(i, j, k_start, p)];
                    }
                }
            }
        }
    }
#endif

#ifdef FRONT_NOSLIP
    if (k_end == NZ)
    {
        for (int i = i_start; i < i_end; i++)
        {
            for (int j = j_start; j < j_end; j++)
            {
                for (int p = 0; p < NP; p++)
                {
                    if (cz[p] > 0)
                    {
                        const int p_bb = p_bounceback[p];
                        const int ic = i + cx[p];
                        const int jc = j + cy[p];

                        f2_a[INDEX_F(ic, jc, k_end, p_bb)] = f2_a[INDEX_F(i, j, k_end - 1, p)];
                        f2_p[INDEX_F(ic, jc, k_end, p_bb)] = f2_p[INDEX_F(i, j, k_end - 1, p)];
                    }
                }
            }
        }
    }
#endif

}