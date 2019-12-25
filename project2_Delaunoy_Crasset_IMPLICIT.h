#ifndef IMPLICIT_H_
#define IMPLICIT_H_

#include "project2_Delaunoy_Crasset_SPARSE.h"
#include "project2_Delaunoy_Crasset_IO.h"

/**
 * Compute the indices of the A matrix this process is responsible for.
 *
 * Parameters:
 * nbproc: The number of process
 * myrank: The rank of this process
 * size: The number of lines in A
 * start: Pointer to an unsigned int that will be set to the starting index
 * start: Pointer to an unsigned int that will be set to the ending index
 */
void get_indices(int nbproc, int myrank, unsigned int size, unsigned int* start, unsigned int* end);

/**
 * Perform a linear system resolution using the conjugate gradient algorithm
 *
 * Parameters:
 * A: A sparse matrix representing the variables coefficients in the equations
 *    the process calling the function is reponsible for
 * b: An array representing the offsets in the equations
 * size: The number of variables in the system
 * rThresh: The threshold at which to stop the algorithm
 * nbproc: The number of processors
 * myrank: The rank of the processor running the function
 * startIndex: The starting index in the A matrix the calling process in responsible for.
 * endIndex: The ending index in the A matrix the calling process in responsible for.
 * recvcounts: The number of elements received by each process on gather calls.
 * dspls: The offsets at which to place the elemnts received by each process.
 * 
 * Returns:
 * A pointer to an allocated array containing the solution of the system
 */
double* MPISparseConjugateGradient(SparseMatrix* A, double* b, unsigned int size, 
                                   double rThresh, int nbproc, int myrank, unsigned int startIndex, 
                                   unsigned int endIndex, int* recvcounts, int* displs);


/**
 * Transform 2d indices into the corresponding 1d index of a flattened array
 *
 * Parameters:
 * x: The value of the first index
 * y: The value of the second index
 * ySize: The size of the second axis of the array
 * 
 * Returns:
 * The equivalent 1d index
 */
int get_1d(int x, int y, int ySize);

/**
 * Transform 1d index into the corresponding 2d indices of a reshaped array
 *
 * Parameters:
 * i: The 1d index
 * ySize: The size of the second axis of the array.
 * x: A pointer to the first index that will be set
 * y: A pointer to the second index that will be set
 */
void get_2d(int i, int ySize, int* x, int* y);

/**
 * Compute the equation matrix of the system
 *
 * Parameters:
 * A: A pointer to a pointer of SparseMatrix that will be set with a SparseMatrix containing the
 *    elements of the equation matrix the calling process is responsible for
 * b: A pointer to an array of double that will be set with the offsets of the equations
 * result: The result of the previous iteration that will be used to fill the matrix A and vector b.
 * start: The starting index the calling process is reponsible for
 * end: The ending index the calling process is repsonsible for
 * xSize: The x axis size of the discretization
 * ySize: The y axis size of the discretization
 * h: The filled h matrix with the depth at rest
 * params: The parameters of the algorithm
 * t: The current iteration
 */
void BuildSystemMatrix(SparseMatrix** A, double** b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t);

/**
 * Init the transposed equation coefficients matrix and offset matrix
 *
 * Parameters:
 * A: A pointed to an allocated equation matrix which elements will be set
 * b: A pointer to an allocated array of double which will be filled with offsets.
 * result: The result of the previous iteration that will be used to fill the matrix A and vector b.
 * start: The starting index the calling process is reponsible for
 * end: The ending index the calling process is repsonsible for
 * xSize: The x axis size of the discretization
 * ySize: The y axis size of the discretization
 * h: The filled h matrix with the depth at rest
 * params: The parameters of the algorithm
 * t: The current iteration
 */
void initAtb(SparseMatrix* A, double* b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t);

/**
 * Make the A matrix semi definite positive and update b accordingly.
 * The elements pointed by At and b will be updated.
 *
 * Parameters:
 * A: A pointer to a pointer of SparseMatrix containing the elements of the equation matrix 
 * the calling process is responsible for
 * b: A pointer to an array of double containing the offsets of the equations
 * result: The result of the previous iteration that will be used to fill the matrix A and vector b.
 * start: The starting index the calling process is reponsible for
 * end: The ending index the calling process is repsonsible for
 * xSize: The x axis size of the discretization
 * ySize: The y axis size of the discretization
 * h: The filled h matrix with the depth at rest
 * params: The parameters of the algorithm
 * t: The current iteration
 */
void makeDefinitePositive(SparseMatrix** At, double** b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t);

/**
 * Save the results.
 *
 * Parameters:
 * x: A pointer to the results
 * xSize: The x axis size of the discretization
 * ySize: The y axis size of the discretization
 * iteration: The current iteration
 * params: The parameters of the algorithm
 * nbproc: The number of processes
 * nbThread: The number of threads
 */
void save_inputs(double* x, int xSize, int ySize);

/**
 * Solve Navier-Stokes equations using the implicit Euler method
 *
 * Parameters:
 * map: A pointer to a map structure
 * params: The parameters with which to solve the problem
 * eta: A pointer to an array of double that will be set with the final iteration results of eta
 * u: A pointer to an array of double that will be set with the final iteration results of u
 * v: A pointer to an array of double that will be set with the final iteration results of v
 * 
 * Returns:
 * A integer indicating whether the function executed successfully
 */
int eulerImplicitMPI(Map* map, Parameters* params, double** eta, double** u, double** v);

#endif