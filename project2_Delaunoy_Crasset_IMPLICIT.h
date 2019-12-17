#ifndef IMPLICIT_H_
#define IMPLICIT_H_

#include "project2_Delaunoy_Crasset_SPARSE.h"

double dotProduct(double* x, double* y, unsigned int size);
double vectorNorm(double* x, unsigned int size);
double* conjugateGradient(double** A, double* b, unsigned int size, double rThresh);
double* sparseConjugateGradient(SparseMatrix* A, double* b, unsigned int size, double rThresh);

#endif