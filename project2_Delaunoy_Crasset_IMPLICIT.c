#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "project2_Delaunoy_Crasset_IMPLICIT.h"
#include "project2_Delaunoy_Crasset_SPARSE.h"
#include "project2_Delaunoy_Crasset_IO.h"

#define M_PI 3.14159265358979323846

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
void get_indices(int nbproc, int myrank, unsigned int size, unsigned int* start, unsigned int* end){
    
    // Compute number of indices a process is responsible for
    int nbIndex = size/nbproc;

    // Compute The remaining indices
    int remaining = size%nbproc;

    // Assign indices to processes
    if(myrank < remaining){
        *start = (nbIndex+1) * myrank;
        *end = (nbIndex+1) * (myrank + 1) - 1;
    }
    else{
        *start = nbIndex * myrank + remaining;
        *end = nbIndex * (myrank + 1) - 1 + remaining;
    }
}

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
                                   unsigned int endIndex, int* recvcounts, int* displs){

    // Allocate memory
    double* xCurr = malloc(size * sizeof(double));
    if(!xCurr)
        return NULL;

    double* xNext = malloc(size * sizeof(double));
    if(!xNext){
        free(xCurr);
        return NULL;
    }

    double* rCurr = malloc(size * sizeof(double));
    if(!rCurr){
        free(xCurr);
        free(xNext);
        return NULL;
    }

    double* rNext = malloc(size * sizeof(double));
    if(!rNext){
        free(xCurr);
        free(xNext);
        free(rCurr);
        return NULL;
    }

    double* pCurr = malloc(size * sizeof(double));
    if(!pCurr){
        free(xCurr);
        free(xNext);
        free(rCurr);
        free(rNext);
        return NULL;
    }

    double* pNext = malloc(size * sizeof(double));
    if(!pNext){
        free(xCurr);
        free(xNext);
        free(rCurr);
        free(rNext);
        free(pCurr);
        return NULL;
    }

    double* Ap = malloc(size * sizeof(double));
    if(!Ap){
        free(xCurr);
        free(xNext);
        free(rCurr);
        free(rNext);
        free(pCurr);
        free(pNext);
        return NULL;
    }

    double* tmpBuff = malloc((endIndex - startIndex + 1) *sizeof(double));
    //double* tmpBuff = malloc(3 *sizeof(double));
    if(!tmpBuff){
        free(xCurr);
        free(xNext);
        free(rCurr);
        free(rNext);
        free(pCurr);
        free(pNext);
        free(Ap);
        return NULL;
    }

    // Initialise x
    for(unsigned int i = 0; i < size; i++){
        xCurr[i] = 0;
    }

    // Initialise r
    for(unsigned int i = 0; i < size; i++){
        rCurr[i] = b[i]; // No need to add Ax_0 term as x_0 = 0
    }

    // Initialise p
    for(unsigned int i = 0; i < size; i++){
        pCurr[i] = rCurr[i];
    }

    // Compute r_0
    double rBaseNorm = MPIDotProduct(rCurr, rCurr, size, startIndex, endIndex, myrank, nbproc);

    double num, denum, alpha, beta, rNextNorm;
    double rCurrNorm = rBaseNorm;
    
    // Loop until convergence reached
    while(rCurrNorm / rBaseNorm >= rThresh){
    
        // Compute alpha

        MPIMatVecMul(A, pCurr, tmpBuff, Ap, startIndex, endIndex, myrank, nbproc, recvcounts, displs);

        alpha = rCurrNorm/MPIDotProduct(pCurr, Ap, size, startIndex, endIndex, myrank, nbproc);

        // Compute xNext
        #pragma omp parallel default(shared)
        { 

            #pragma omp for schedule(static)
            for(unsigned int i = 0; i < size; i++)
                xNext[i] = xCurr[i] + alpha * pCurr[i];

            // Compute rNext
            #pragma omp for schedule(static)
            for(unsigned int i = 0; i < size; i++)
                rNext[i] = rCurr[i] - alpha * Ap[i];
        }

        // Compute next rNorm
        rNextNorm = MPIDotProduct(rNext, rNext, size, startIndex, endIndex, myrank, nbproc);
        
        // Compute beta
        beta = rNextNorm/rCurrNorm;

        // Compute pNext
        #pragma omp parallel default(shared)
        { 

            #pragma omp for schedule(static)
            for(unsigned int i = 0; i < size; i++)
                pNext[i] = rNext[i] + beta * pCurr[i];
        }
          
        // Go to next iteration
        double* buf;

        buf = xNext;
        xNext = xCurr;
        xCurr = buf;

        buf = rNext;
        rNext = rCurr;
        rCurr = buf;

        buf = pNext;
        pNext = pCurr;
        pCurr = buf;

        rCurrNorm = rNextNorm;
    }

    // Free memory
    free(tmpBuff);
    free(xNext);
    free(rCurr);
    free(rNext);
    free(pCurr);
    free(pNext);
    free(Ap);

    return xCurr;
}

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
int get_1d(int x, int y, int ySize){
    return x * ySize + y;
}

/**
 * Transform 1d index into the corresponding 2d indices of a reshaped array
 *
 * Parameters:
 * i: The 1d index
 * ySize: The size of the second axis of the array.
 * x: A pointer to the first index that will be set
 * y: A pointer to the second index that will be set
 */
void get_2d(int i, int ySize, int* x, int* y){
    *x = i/ySize;
    *y = i%ySize;
}


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
void BuildSystemMatrix(SparseMatrix** A, double** b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t){
    
    // Compute the indices of each part of the result vector
    int uBegin = (xSize + 1) * (ySize + 1);
    int vBegin = uBegin + (xSize + 2) * (ySize + 1);
    int vEnd = vBegin + (xSize + 1) * (ySize + 2) - 1;

    // Allocate memory
    *A = createSparseMatrix(start, end, 5);
    *b = malloc((vEnd + 1) * sizeof(double));

    // Initialize the matrices
    initAtb(*A, *b, result, start, end, xSize, ySize, h, params, t);
   
    // Make those semi-definite positive
    makeDefinitePositive(A, b, result, start, end, xSize, ySize, h, params, t);
}

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
void makeDefinitePositive(SparseMatrix** At, double** b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t){
    
    // Allocate new A
    SparseMatrix* AtA = createSparseMatrix(start, end, 25);

    // Compute the indices of each part of the result vector
    int uBegin = (xSize + 1) * (ySize + 1);
    int vBegin = uBegin + (xSize + 2) * (ySize + 1);
    int vEnd = vBegin + (xSize + 1) * (ySize + 2) - 1;

    // Allocate new b
    double* Atb = malloc((vEnd + 1) * sizeof(double));
    
    // Construct full At and mult by this
    #pragma omp parallel default(shared)
    {   
        // Temporary variable used for computation
        int x, y;
        SparseVector* vec = createSparseVector(5);

        // Loop on each line of A 
        #pragma omp for schedule(static)
        for(int i = 0; i <= vEnd; i++){
            
            // Set vec with the indices associated with the current line performed

            //eta variable
            if(i < uBegin){
                get_2d(i, ySize + 1, &x, &y);
                sparseVecInsertElement(vec, i, 1.0/params->deltaT);

                if(x != 0)
                    sparseVecInsertElement(vec, uBegin + get_1d(x, y, ySize + 1), params->g/params->deltaX);
                if(x != xSize)
                    sparseVecInsertElement(vec, uBegin + get_1d(x + 1, y, ySize + 1), -params->g/params->deltaX);

                if(y != 0)
                    sparseVecInsertElement(vec, vBegin + get_1d(x, y, ySize + 2), params->g/params->deltaY);
                
                if(y != ySize)
                    sparseVecInsertElement(vec, vBegin + get_1d(x, y + 1, ySize + 2), -params->g/params->deltaY);
            }

            //u variable
            else if(i < vBegin){
                get_2d(i - uBegin, ySize + 1, &x, &y);

                if(x != 0)
                    sparseVecInsertElement(vec, get_1d(x - 1, y, ySize + 1), h[2*x-1][2*y]/params->deltaX);

                if(x == 0)
                    sparseVecInsertElement(vec, get_1d(x, y, ySize + 1), -h[0][2*y]/params->deltaX);
                else if(x != xSize + 1)
                    sparseVecInsertElement(vec, get_1d(x, y, ySize + 1), -h[2*x-1][2*y]/params->deltaX);

                // At the border
                if(x == 0 || x == xSize + 1){
                    sparseVecInsertElement(vec, i, 1.0);
                }
                else{
                    sparseVecInsertElement(vec, i, 1.0/params->deltaT + params->gamma);
                }

            }

            //v variable
            else{
                get_2d(i - vBegin, ySize + 2, &x, &y);

                if(y != 0)
                    sparseVecInsertElement(vec, get_1d(x, y - 1, ySize + 1), h[2*x][2*y-1]/params->deltaY);

                if(y == 0)
                    sparseVecInsertElement(vec, get_1d(x, y, ySize + 1), -h[2*x][0]/params->deltaY);
                else if (y != ySize + 1)
                    sparseVecInsertElement(vec, get_1d(x, y, ySize + 1), -h[2*x][2*y-1]/params->deltaY);
                
                // At the border
                if(y == 0 || y == ySize + 1){
                    sparseVecInsertElement(vec, i, 1.0);
                }

                else{
                    sparseVecInsertElement(vec, i, 1.0/params->deltaT + params->gamma);
                }
            }

            // Mult A by this vec and set the corresponding elements of AtA
            double result;
            for(int j = start; j <= end; j++){
                result = sparseMatVecDotProduct(*At, j, vec);
                if(result != 0.0){
                    #pragma omp critical(fill_AtA)
                        sparseInsertElement(AtA, j, i, result);
                }
            }

            // Mult b by this vec and set the corresponding element of Atb
            Atb[i] = vecSparseDotProduct(vec, *b);

            // Reset buffer vector
            resetSparseVector(vec);
        }

        // Free space of tempary vector
        freeSparseVector(vec);
    }

    // Free old matrices
    freeSparseMatrix(*At);
    free(*b);

    // Replace old matrices with the new ones
    *At = AtA;
    *b = Atb;

}

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
void initAtb(SparseMatrix* A, double* b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t){
    
    // Reset the matrix
    resetSparseMatrix(A);

    // Compute the indices of each part of the result vector
    int uBegin = (xSize + 1) * (ySize + 1);
    int vBegin = uBegin + (xSize + 2) * (ySize + 1);
    int vEnd = vBegin + (xSize + 1) * (ySize + 2) - 1;

    // Fill At
    #pragma omp parallel default(shared)
    {   
        int x, y;

        // Loop on lines of At
        #pragma omp for schedule(static)
        for(int i = start; i <= end; i++){

            //eta variable
            if(i < uBegin){
                get_2d(i, ySize + 1, &x, &y);
                sparseInsertElement(A, i, i, 1.0/params->deltaT);

                if(x != 0)
                    sparseInsertElement(A, i, uBegin + get_1d(x, y, ySize + 1), params->g/params->deltaX);
                if(x != xSize)
                    sparseInsertElement(A, i, uBegin + get_1d(x + 1, y, ySize + 1), -params->g/params->deltaX);

                if(y != 0)
                    sparseInsertElement(A, i, vBegin + get_1d(x, y, ySize + 2), params->g/params->deltaY);
                
                if(y != ySize)
                    sparseInsertElement(A, i, vBegin + get_1d(x, y + 1, ySize + 2), -params->g/params->deltaY);
            }

            //u variable
            else if(i < vBegin){
                get_2d(i - uBegin, ySize + 1, &x, &y);

                if(x != 0)
                    sparseInsertElement(A, i, get_1d(x - 1, y, ySize + 1), h[2*x-1][2*y]/params->deltaX);

                if(x == 0)
                    sparseInsertElement(A, i, get_1d(x, y, ySize + 1), -h[0][2*y]/params->deltaX);
                else if(x != xSize + 1)
                    sparseInsertElement(A, i, get_1d(x, y, ySize + 1), -h[2*x-1][2*y]/params->deltaX);

                // At the border
                if(x == 0 || x == xSize + 1){
                    sparseInsertElement(A, i, i, 1.0);
                }
                else{
                    sparseInsertElement(A, i, i, 1.0/params->deltaT + params->gamma);
                }

            }

            //v variable
            else{
                get_2d(i - vBegin, ySize + 2, &x, &y);

                if(y != 0)
                    sparseInsertElement(A, i, get_1d(x, y - 1, ySize + 1), h[2*x][2*y-1]/params->deltaY);

                if(y == 0)
                    sparseInsertElement(A, i, get_1d(x, y, ySize + 1), -h[2*x][0]/params->deltaY);
                else if (y != ySize + 1)
                    sparseInsertElement(A, i, get_1d(x, y, ySize + 1), -h[2*x][2*y-1]/params->deltaY);
                
                // At the border
                if(y == 0 || y == ySize + 1){
                    sparseInsertElement(A, i, i, 1.0);
                }

                else{
                    sparseInsertElement(A, i, i, 1.0/params->deltaT + params->gamma);
                }
            }
        }

        // Fill b
        #pragma omp for schedule(static)
        for(int i = 0; i <= vEnd; i++){
            
            // eta constraint
            if(i < uBegin){
                b[i] = result[i]/params->deltaT;
            }

            // u constraint
            else if(i < vBegin){
                get_2d(i - uBegin, ySize + 1, &x, &y);

                // At the border
                if(x == 0 || x == xSize + 1){
                    b[i] = 0;
                }

                // In the middle
                else{
                    b[i] = result[i]/params->deltaT;
                }
            }

            // v constraint
            else{
                get_2d(i - vBegin, ySize + 2, &x, &y);
                // At the bottom border
                if(y == 0){
                    b[i] = 0;
                }

                // At the top border
                else if(y == ySize + 1){
                    if(params->s == 0)
                        b[i] = params->A * sin(2 * M_PI * params->f * t * params->deltaT);
                    else
                        b[i] = params->A * sin(2 * M_PI * params->f * t * params->deltaT) * exp(- t * params->deltaT / 500);
                }

                // In the middle
                else{
                    b[i] = result[i]/params->deltaT;
                }
            }
        }
    }
}

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
void saveImplicit(double* x, int xSize, int ySize, unsigned int iteration, Parameters* params, int nbproc, int nbThread){

    // Compute the indices of each part of the result vector
    int uBegin = (xSize + 1) * (ySize + 1);
    int vBegin = uBegin + (xSize + 2) * (ySize + 1);

    // Save results
    saveToDisk(x, x + uBegin, x + vBegin, xSize, ySize, iteration, params, nbproc, nbThread);
}

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
int eulerImplicitMPI(Map* map, Parameters* params, double** eta, double** u, double** v, int debug, int debug_rank){
    
    // Assert arugments
    assert(map);
    assert(params);

    // Initialise the number of processes and process num
    int nbproc, myrank ;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    int openMP_nbthreads = omp_get_num_threads();

    // Compute the size of the discretization
    int xSize = (int)(map->a / params->deltaX);
    int ySize = (int)(map->b / params->deltaY);

    // Allocate memory
    // x = [eta u v]
    int inputSize = (xSize + 1) * (ySize + 1) + (xSize + 2) * (ySize + 1) + (xSize + 1) * (ySize + 2);

    double* x = calloc(inputSize, sizeof(double));
    if(!x)
        return -1;

    double** h = allocateDoubleMatrix(2 * xSize + 3, 2 * ySize + 3);
    if(!h){
        free(x);
        return -1;
    }

    // Fill depth at rest matrix
    for(int i = 0; i < 2 * xSize + 3; i++){
        for(int j = 0; j < 2 * ySize + 3; j++){
            h[i][j] = getGridValueAtDomainCoordinates(map, ((float)(i * xSize)/(xSize + 1)) * (params->deltaX / 2), ((float)(j * ySize)/(ySize + 1)) * (params->deltaY / 2));
        }
    }

    // Compute the indices this process will be responsible for
    unsigned int start, end;
    get_indices(nbproc, myrank, inputSize, &start, &end);

    // Loop over ietrations to perform
    for(unsigned int t = 1; t <= params->TMax/params->deltaT; t++){

        if(myrank == 0)
            fprintf(stderr, "t = %d\n", t);
        
        // Build the equation matrix
        SparseMatrix* A;
        double* b;
        BuildSystemMatrix(&A, &b, x, start, end, xSize, ySize, h, params, t);
        
        // Free no longer used results of previous iteration
        free(x);

        // Compute received counts and displacements
        unsigned int tmpStart, tmpEnd;
        int* recvcounts = malloc(nbproc * sizeof(int));
        int* displs = malloc(nbproc * sizeof(int));
        
        displs[0] = 0;
        for(unsigned int k = 0; k < nbproc; k++){
            get_indices(nbproc, k, inputSize, &tmpStart, &tmpEnd);
            int tmpSize = tmpEnd - tmpStart + 1;
            recvcounts[k] = tmpSize;
            if(k < nbproc - 1)
                displs[k+1] = displs[k] + tmpSize;
        }
    
        // Get threshold       
        double rThresh = params->r_threshold;
        
        // Solve system using conjugate gradient method
        x = MPISparseConjugateGradient(A, b, inputSize, rThresh, nbproc, myrank, start, end, recvcounts, displs);

        // Process 0 saves arrays to disk
        if(params->S != 0 && t % params->S == 0 && myrank == 0){
            // Save to disk
            saveImplicit(x, xSize, ySize, t, params, nbproc, openMP_nbthreads);
        }

        // Free memory unused
        free(recvcounts);
        free(displs);
        freeSparseMatrix(A);
        free(b);
    }

    // free depth at rest matrix
    freeDoubleMatrix(h, 2*xSize+3);

    // Save last iteration solution
    if(params->S != 0 && ((int) (params->TMax/params->deltaT) % params->S) == 0 && myrank == 0){
        saveImplicit(x, xSize, ySize, params->TMax/params->deltaT, params, nbproc, openMP_nbthreads);
    }

    // Set the solution pointers
    int uBegin = (xSize + 1) * (ySize + 1);
    int vBegin = uBegin + (xSize + 2) * (ySize + 1);

    *eta = x;
    *u = x + uBegin;
    *v = x + vBegin;
    
    return 0;
}