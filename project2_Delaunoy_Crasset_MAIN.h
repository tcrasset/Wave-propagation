#ifndef MAIN_H_
#define MAIN_H_

#include "project2_Delaunoy_Crasset_IO.h"

typedef struct SparseMatrix {
    unsigned int nonzeroSize;
    unsigned int* m;
    unsigned int* n;
    unsigned int* v;
} SparseMatrix;

SparseMatrix* toSparseMatrix(double** matrix, int xSize, int ySize);
void printUsefulMapInformation(Map* map);
double* transformMatrixToArray(double** matrix, int x, int y);
void get_array_sizes(int rank, int nbproc, int xSize, int* size_X, int* size_X_u, int* size_X_h, int* startval_X_h, int* endval_X_h);
int gather_and_save(double** eta, double**  u, double**  v, int xSize, int ySize,  int debug, unsigned int iteration, Parameters* params);
int eulerExplicitMPI(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug, int debug_rank);
#endif  // MAIN_H_