#include <math.h>
#include <string.h>

#include "../definitions.h"
#include "../include/collide.h"

void f_eq(const double rho, const double u, const double v, const double w, const int p, SimulationBag *sim)
{
    UNPACK_STENCIL

    const double u2 = u*u + v*v + w*w;
    const double uc = u*(double)cx[p] + v*(double)cy[p] + w*(double)cz[p];

    return wp[p] * rho * (1.0 + 3.0*uc + 4.5*uc*uc - 1.5*u2);
}

void g_eq(SimulationBag *sim)
{

}

void collide(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho)
    READ_FIELD(phi)
    READ_FIELD(u)
    READ_FIELD(v)
    READ_FIELD(w)
    READ_FIELD(Fx)
    READ_FIELD(Fy)
    READ_FIELD(Fz)
    
    WRITE_DIST(f1)
    WRITE_DIST(f2)
    WRITE_DIST(g1)
    WRITE_DIST(g2)

    FOR_DOMAIN
    {
        const int idx = INDEX(i, j, k);

        // TODO
    }
}