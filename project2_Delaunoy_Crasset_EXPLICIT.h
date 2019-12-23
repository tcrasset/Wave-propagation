#ifndef EXPLICIT_H_
#define EXPLICIT_H_

#define MAX_FILENAME_SIZE 500

#include "project2_Delaunoy_Crasset_IO.h"

/**
 * Compute the size of the arrays this process is responsible for
 *
 * Parameters:
 * rank: The rank of the calling process
 * nbproc: The number of processes
 * xSize: The discretization along the x axis
 * size_X: A pointer to an integer that will be set to the x size of eta and v
 * size_X_u: A pointer to an integer that will be set to the x size of u
 * size_X: A pointer to an integer that will be set to the x size of h
 * startval_X_h: A pointer to an integer that will be set to the starting value of h
 * endval_X_h: A pointer to an integer that will be set to the ending value of h
 */
void get_array_sizes(int rank, int nbproc, int xSize, int* size_X, int* size_X_u, int* size_X_h, int* startval_X_h, int* endval_X_h);

/**
 * Gather results from all process and save to disk
 *
 * Parameters:
 * eta: The eta array of the calling process
 * u: The u array of the calling process
 * v: The v array of the calling process
 * xSize: The discretization size along the x axis
 * ySize: The discretization size along the y axis
 * iteration: The iteration at which the save is performed
 * params: The structure holding the parameters of the run
 */
void gather_and_save(double** eta, double**  u, double**  v, int xSize, int ySize, unsigned int iteration, Parameters* params);

/**
 * Solve the Navier-Stockes equations using explicit Euler method
 *
 * Parameters:
 * map: A structure containing the map infos
 * params: The parameters of the run
 * eta: A pointer to a matrix that will be set to the result of eta
 * u: A pointer to a matrix that will be set to the result of u
 * v: A pointer to a matrix that will be set to the result of v
 *
 * Returns:
 * An integer indicating whether the algorithm run with success or not
 */
int eulerExplicitMPI(Map* map, Parameters* params, double*** eta, double*** u, double*** v);

#endif  // EXPLICIT_H_