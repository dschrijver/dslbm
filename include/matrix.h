#ifndef MATRIX_H
#define MATRIX_H

#include "datatypes.h"

void matmul(matrix A, matrix B, matrix C);
void matsca(matrix A, double s, matrix C);
double matcontraction(matrix A, matrix B);

#endif