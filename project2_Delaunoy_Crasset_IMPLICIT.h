#ifndef IMPLICIT_H_
#define IMPLICIT_H_

#include "project2_Delaunoy_Crasset_SPARSE.h"
#include "project2_Delaunoy_Crasset_IO.h"

double dotProduct(double* x, double* y, unsigned int size);
double vectorNorm(double* x, unsigned int size);
double* conjugateGradient(double** A, double* b, unsigned int size, double rThresh);
double* sparseConjugateGradient(SparseMatrix* A, double* b, unsigned int size, double rThresh);
double* MPISparseConjugateGradient(SparseMatrix* A, double* b, unsigned int size, 
                                   double rThresh, int nbproc, int myrank, unsigned int startIndex, 
                                   unsigned int endIndex, int* recvcounts, int* displs);
void initAb(SparseMatrix* A, double* b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t);
void initAtb(SparseMatrix* A, double* b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t);
void makeDefinitePositive(SparseMatrix** At, double** b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t);
void BuildSystemMatrix(SparseMatrix** A, double** b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t);
void save_inputs(double* x, int xSize, int ySize);
int eulerImplicitMPI(Map* map, Parameters* params, double** eta, double** u, double** v, int debug, int debug_rank);

#endif