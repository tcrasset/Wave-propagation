#ifndef EXPLICIT_H_
#define EXPLICIT_H_

#define MAX_FILENAME_SIZE 500

#include "project2_Delaunoy_Crasset_IO.h"

double bilinearInterpolation(Map* map, double x, double y);
double getGridValueAtDomainCoordinates(Map* map, double x, double y);
double** allocateDoubleMatrix(int x, int y);
void freeDoubleMatrix(double** matrix, int x, int debug);
double* transformMatrixToArray(double** matrix, int x, int y);
void get_array_sizes(int rank, int nbproc, int xSize, int* size_X, int* size_X_u, int* size_X_h, int* startval_X_h, int* endval_X_h);
void gather_and_save(double** eta, double**  u, double**  v, int xSize, int ySize,  int debug, unsigned int iteration, Parameters* params);
int eulerExplicitMPI(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug, int debug_rank);

#endif  // EXPLICIT_H_