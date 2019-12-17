#ifndef IMPLICIT_H_
#define IMPLICIT_H_

double dotProduct(double* x, double* y, unsigned int size);
double vectorNorm(double* x, unsigned int size);
double* conjugateGradient(double** A, double* b, unsigned int size, double rThresh);

#endif