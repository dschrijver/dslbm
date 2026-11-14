#include "../include/datatypes.h"

void matmul(matrix A, matrix B, matrix C)
{
    for (int a = 0; a < 3; a++)
    {
        for (int b = 0; b < 3; b++)
        {
            double result = 0.0;
            for (int c = 0; c < 3; c++)
            {
                result += A[a][c] * B[c][b];
            }
            C[a][b] = result;
        }
    }
}

void matsca(matrix A, double s, matrix C)
{
    for (int a = 0; a < 3; a++)
    {
        for (int b = 0; b < 3; b++)
        {
            C[a][b] = s * A[a][b];
        }
    }
}

double matcontraction(matrix A, matrix B)
{
    double result = 0.0;
    for (int a = 0; a < 3; a++)
    {
        for (int b = 0; b < 3; b++)
        {
            result += A[a][b]*B[a][b];
        }
    }
    return result;
}