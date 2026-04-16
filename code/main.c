#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Definitions
#define IDX(i,j) (NP*j + i)

// Stencil
static int cx[] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1};
static int cy[] = {0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1};
static int cz[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1};
static int p_bounceback[] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15, 26, 25, 24, 23, 22, 21, 20, 19};
static double wp[] = {8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0};
static int NP = 27;
static double cs2 = 1.0/3.0;

// Simulation
#define N 50
static int nhat[] = {0, 1, 0};

// Choices 
static double omega = 1.0/1.2;

// Global variables;
double *rhobulk;
double *ubulk;
double *vbulk;
double *wbulk;
double *rhokn;
double *rhokn_old;
double *ukn;
double *vkn;
double *wkn;
double *Fx;
double *Fy;
double *Fz;

double *F[3];

double *fkn;
double *geq;

double residual()
{
    double item;
    double result = 0.0;
    for (int j = 0; j < N; j++)
    {
        item = fabs((rhokn[j] - rhokn_old[j]));
        if (item > result) 
            result = item;
    }
    return result;
}

double eq(int i, double rho, double u2, double u, double v, double w)
{
    double uc = (double)cx[i]*u + (double)cy[i]*v + (double)cz[i]*w;
    return wp[i]*rho*(1.0 + uc/cs2 + (uc*uc)/(2.0*cs2*cs2) - u2/(2.0*cs2));
}

void compute_geq()
{
    double rho;
    double u2;
    double u;
    double v;
    double w;
    double u2bulk;

    for (int j = 0; j < N; j++)
    {
        rho = rhobulk[j] + rhokn[j];
        u = (rhobulk[j]*ubulk[j] + rhokn[j]*ukn[j])/rho;
        v = (rhobulk[j]*vbulk[j] + rhokn[j]*vkn[j])/rho;
        w = (rhobulk[j]*wbulk[j] + rhokn[j]*wkn[j])/rho;
        u2 = u*u + v*v + w*w;
        u2bulk = ubulk[j]*ubulk[j] + vbulk[j]*vbulk[j] + wbulk[j]*wbulk[j];
        for (int i = 0; i < NP; i++)
        {
            geq[IDX(i,j)] = eq(i, rho, u2, u, v, w) - eq(i, rhobulk[j], u2bulk, ubulk[j], vbulk[j], wbulk[j]);
        }
    }
}

void momentum_corrections(double *Nmom)
{
    for (int alpha = 0; alpha < 3; alpha++)
    {
        if (nhat[alpha] > 0.0)
        {
            Nmom[alpha] = 1.0/cs2*F[alpha][0];
        }
        else 
        {
            Nmom[alpha] = 1.0/(cs2*cs2)*F[alpha][0];
        }
    }
}

void macroscopic()
{
    // Compute macroscopic fields
    for (int j = 0; j < N; j++)
    {
        rhokn[j] = 0.0;
        ukn[j] = 0.0;
        vkn[j] = 0.0;
        wkn[j] = 0.0;
        for (int i = 0; i < NP; i++)
        {
            rhokn[j] += fkn[IDX(i,j)];
            ukn[j] += (double)cx[i] * fkn[IDX(i,j)];
            vkn[j] += (double)cy[i] * fkn[IDX(i,j)];
            wkn[j] += (double)cz[i] * fkn[IDX(i,j)];
        }
        if (rhokn[j] == 0.0)
        {
            ukn[j] = 0.0;
            vkn[j] = 0.0;
            wkn[j] = 0.0;
        }
        else
        {
            ukn[j] /= rhokn[j];
            vkn[j] /= rhokn[j];
            wkn[j] /= rhokn[j];
        }
    }
}

void update()
{
    int cn;
    double uc, cN;
    double correction;
    double Nmom[3];

    compute_geq();

    // Outgoing, parallel and zero-velocity populations
    for (int i = 0; i < NP; i++)
    {
        cn = cx[i]*nhat[0] + cy[i]*nhat[1] + cz[i]*nhat[2];
        
        if (cn == 0)
        {
            for (int j = 0; j < N; j++)
            {
                fkn[IDX(i,j)] = geq[IDX(i,j)];
            }
        }
        else if (cn < 0)
        {
            for (int j = 0; j < N; j++)
            {
                fkn[IDX(i,j)] = 0.0;
                for (int k = j+1; k < N; k++)
                {
                    fkn[IDX(i,j)] += pow((1.0 - omega), fabs((double)(k-j)))*geq[IDX(i,k)];
                }
                fkn[IDX(i,j)] *= omega/(1.0-omega);
            }
        }
    }

    // Boundary condition 
    momentum_corrections(Nmom);
    for (int i = 0; i < NP; i++)
    {
        cn = cx[i]*nhat[0] + cy[i]*nhat[1] + cz[i]*nhat[2];
        if (cn > 0)
        {
            uc = (double)cx[i]*ukn[0] + (double)cy[i]*vkn[0] + (double)cz[i]*wkn[0];
            cN = (double)cx[i]*Nmom[0] + (double)cy[i]*Nmom[1] + (double)cz[i]*Nmom[2];
            fkn[IDX(i,0)] = fkn[IDX(p_bounceback[i],0)] + 2.0*wp[i]*rhokn[0]*uc/cs2 - wp[i]*cN;
        }
    }

    // Incoming populations
    for (int i = 0; i < NP; i++)
    {
        cn = cx[i]*nhat[0] + cy[i]*nhat[1] + cz[i]*nhat[2];
        
        if (cn > 0)
        {
            for (int j = 1; j < N; j++)
            {
                fkn[IDX(i,j)] = 0.0;
                for (int k = 0; k < j; k++)
                {
                    fkn[IDX(i,j)] += pow((1.0 - omega), fabs((double)(k-j)))*geq[IDX(i,k)];
                }
                fkn[IDX(i,j)] = omega/(1.0-omega)*fkn[IDX(i,j)] + pow((1.0-omega), (double)j)*fkn[IDX(i,0)];
            }
        }
    }

    // Mass correction to zero-velocity population
    correction = 0.5*Fy[0];
    for (int i = 0; i < NP; i++)
    {
        cn = cx[i]*nhat[0] + cy[i]*nhat[1] + cz[i]*nhat[2];
        if (cn < 0){
            correction += -omega*(fkn[IDX(i,0)]-geq[IDX(i,0)]);
        }
    }
    fkn[IDX(0,0)] += correction;
}

int main(void)
{
    // Choices 
    rhobulk = (double*) malloc(N*sizeof(double));
    ubulk = (double*) malloc(N*sizeof(double));
    vbulk = (double*) malloc(N*sizeof(double));
    wbulk = (double*) malloc(N*sizeof(double));
    Fx = (double*) malloc(N*sizeof(double));
    Fy = (double*) malloc(N*sizeof(double));
    Fz = (double*) malloc(N*sizeof(double));
    for (int j = 0; j < N; j++)
    {
        rhobulk[j] = 1.0;
        ubulk[j] = 0.0;
        vbulk[j] = 0.0;
        wbulk[j] = 0.0;
        Fx[j] = 0.0;
        Fy[j] = 0.0;
        Fz[j] = 0.0;
    }
    Fy[0] = 2.0;

    // Allocate
    fkn = (double*) malloc(N*NP*sizeof(double));
    geq = (double*) malloc(N*NP*sizeof(double));
    rhokn = (double*) malloc(N*sizeof(double));
    rhokn_old = (double*) malloc(N*sizeof(double));
    ukn = (double*) malloc(N*sizeof(double));
    vkn = (double*) malloc(N*sizeof(double));
    wkn = (double*) malloc(N*sizeof(double));
    F[0] = Fx;
    F[1] = Fy;
    F[2] = Fz;

    // Initial condition
    for (int j = 0; j < N; j++)
    {
        rhokn[j] = 1.0;
        ukn[j] = 0.0;
        vkn[j] = 0.0;
        wkn[j] = 0.0;
    }

    int niters = 0;
    double res;
    double lambda = 1e-6;
    while (niters < 10000000)
    {
        update();

        memcpy(rhokn_old, rhokn, N*sizeof(double));
        macroscopic();

        for (int j = 0; j < N; j++)
        {
            rhokn[j] = lambda*rhokn[j] + (1.0 - lambda)*rhokn_old[j];
        }

        res = residual();
        if (niters % 100 == 0)
        {
            printf("\rIters = %d, Res = %.3e", niters+1, res);
            fflush(stdout);
        }

        if (niters % 5000 == 0)
        {
            printf("\n");
            for (int j = 0; j < N; j++)
            {
                printf("j = %d: %.3e %.3e %.3e %.3e\n", j, rhokn[j], rhokn[j]*ukn[j], rhokn[j]*vkn[j], rhokn[j]*wkn[j]);
            }
        }

        niters++;
        if (res < 1e-12)
            break;
    }
    printf("\rTotal iters = %d                         \n", niters);
    printf("\n");
    for (int j = 0; j < N; j++)
    {
        printf("j = %d: %.3e %.3e %.3e %.3e\n", j, rhokn[j], rhokn[j]*ukn[j], rhokn[j]*vkn[j], rhokn[j]*wkn[j]);
    }

    

    return 0;
}