#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>

#include "project2_Delaunoy_Crasset_IMPLICIT.h"
#include "project2_Delaunoy_Crasset_SPARSE.h"
#include "project2_Delaunoy_Crasset_IO.h"

#define M_PI 3.14159265358979323846

double dotProduct(double* x, double* y, unsigned int size){
    //fprintf(stderr, "dot product \n");
    double result = 0.0;

    for(unsigned int i = 0; i < size; i++){
        //fprintf(stderr, "x[%d] = %lf, y[%d] = %lf\n", i, x[i], i, y[i]);
        result += x[i] * y[i];
    }

    return result;
}

double vectorNorm(double* x, unsigned int size){
    double norm = 0.0;

    for(unsigned int i = 0; i < size; i++){
        norm += x[i] * x[i];
    }

    return norm;
}

double* conjugateGradient(double** A, double* b, unsigned int size, double rThresh){
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

    // Initialise x
    for(unsigned int i = 0; i < size; i++){
        xCurr[i] = 0;
    }

    // Initialise r
    for(unsigned int i = 0; i < size; i++){
        rCurr[i] = b[i]; // No need to add Ax_0 term as x_0 = 0
    }

    double rBaseNorm = vectorNorm(rCurr, size);

    // Initialise p
    for(unsigned int i = 0; i < size; i++){
        pCurr[i] = rCurr[i];
    }

    double num, denum, alpha, beta, rNextNorm;
    double rCurrNorm = rBaseNorm;
    
    while(rCurrNorm / rBaseNorm >= rThresh){
        fprintf(stderr, "loop \n");
    
        // Compute alpha

        fprintf(stderr, "pCurr = ");
        for(unsigned int i = 0; i < size; i++)
            fprintf(stderr, "%lf ", pCurr[i]);
        fprintf(stderr, "\n");

        fprintf(stderr, "Ap =");
        for(unsigned int i = 0; i < size; i++){
            Ap[i] = dotProduct(A[i], pCurr, size);
            fprintf(stderr, " %lf", Ap[i]);
        }
        fprintf(stderr, "\n");

        alpha = rCurrNorm/dotProduct(pCurr, Ap, size);

        fprintf(stderr, "alpha = %lf\n", alpha);

        // Compute xNext
        fprintf(stderr, "xNext =");
        for(unsigned int i = 0; i < size; i++){
            xNext[i] = xCurr[i] + alpha * pCurr[i];
            fprintf(stderr, " %lf", xNext[i]);
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "rNext =");
        // Compute rNext
        for(unsigned int i = 0; i < size; i++){
            rNext[i] = rCurr[i] - alpha * Ap[i];
            fprintf(stderr, " %lf", rNext[i]);
        }
        fprintf(stderr, "\n");

        rNextNorm = vectorNorm(rNext, size);
        // Compute beta
        beta = rNextNorm/rCurrNorm;
        fprintf(stderr, "beta = %lf\n", beta);

        fprintf(stderr, "pNext =");
        // Compute pNext
        for(unsigned int i = 0; i < size; i++){
            pNext[i] = rNext[i] + beta * pCurr[i];
            fprintf(stderr, " %lf", pNext[i]);
        }
        fprintf(stderr, "\n");

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

    free(xNext);
    free(rCurr);
    free(rNext);
    free(pCurr);
    free(pNext);
    free(Ap);

    return xCurr;
}

double* sparseConjugateGradient(SparseMatrix* A, double* b, unsigned int size, double rThresh){
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

    // Initialise x
    for(unsigned int i = 0; i < size; i++){
        xCurr[i] = 0;
    }

    // Initialise r
    for(unsigned int i = 0; i < size; i++){
        rCurr[i] = b[i]; // No need to add Ax_0 term as x_0 = 0
    }

    double rBaseNorm = vectorNorm(rCurr, size);

    // Initialise p
    for(unsigned int i = 0; i < size; i++){
        pCurr[i] = rCurr[i];
    }

    double num, denum, alpha, beta, rNextNorm;
    double rCurrNorm = rBaseNorm;
    
    while(rCurrNorm / rBaseNorm >= rThresh){
        fprintf(stderr, "loop \n");
    
        // Compute alpha

        fprintf(stderr, "pCurr = ");
        for(unsigned int i = 0; i < size; i++)
            fprintf(stderr, "%lf ", pCurr[i]);
        fprintf(stderr, "\n");

        fprintf(stderr, "Ap =");
        for(unsigned int i = 0; i < size; i++){
            Ap[i] = sparseDotProduct(A, i, pCurr);
            fprintf(stderr, " %lf", Ap[i]);
        }
        fprintf(stderr, "\n");

        alpha = rCurrNorm/dotProduct(pCurr, Ap, size);

        fprintf(stderr, "alpha = %lf\n", alpha);

        // Compute xNext
        fprintf(stderr, "xNext =");
        for(unsigned int i = 0; i < size; i++){
            xNext[i] = xCurr[i] + alpha * pCurr[i];
            fprintf(stderr, " %lf", xNext[i]);
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "rNext =");
        // Compute rNext
        for(unsigned int i = 0; i < size; i++){
            rNext[i] = rCurr[i] - alpha * Ap[i];
            fprintf(stderr, " %lf", rNext[i]);
        }
        fprintf(stderr, "\n");

        rNextNorm = vectorNorm(rNext, size);
        // Compute beta
        beta = rNextNorm/rCurrNorm;
        fprintf(stderr, "beta = %lf\n", beta);

        fprintf(stderr, "pNext =");
        // Compute pNext
        for(unsigned int i = 0; i < size; i++){
            pNext[i] = rNext[i] + beta * pCurr[i];
            fprintf(stderr, " %lf", pNext[i]);
        }
        fprintf(stderr, "\n");

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

    free(xNext);
    free(rCurr);
    free(rNext);
    free(pCurr);
    free(pNext);
    free(Ap);

    return xCurr;
}

void get_indices(int nbproc, int myrank, unsigned int size, unsigned int* start, unsigned int* end){
    int nbIndex = size/nbproc;
    int remaining = size%nbproc;

    if(myrank < remaining){
        *start = (nbIndex+1) * myrank;
        *end = (nbIndex+1) * (myrank + 1) - 1;
    }
    else{
        *start = nbIndex * myrank + remaining;
        *end = nbIndex * (myrank + 1) - 1 + remaining;
    }
}

double* MPISparseConjugateGradient(SparseMatrix* A, double* b, unsigned int size, unsigned int procSize, 
                                   double rThresh, int nbproc, int myrank, unsigned int startIndex, 
                                   unsigned int endIndex, int* recvcounts, int* displs){

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

    double* tmpBuff = malloc((startIndex - endIndex + 1) *sizeof(double));
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

    double rBaseNorm = MPIDotProduct(rCurr, rCurr, size, startIndex, endIndex, myrank, nbproc);

    //MPI_Allreduce(&rBaseNorm_, &rBaseNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD)

    double num, denum, alpha, beta, rNextNorm;
    double rCurrNorm = rBaseNorm;
    
    while(rCurrNorm / rBaseNorm >= rThresh){
    
        // Compute alpha

        MPIMatVecMul(A, pCurr, tmpBuff, Ap, startIndex, endIndex, myrank, nbproc, recvcounts, displs);

        alpha = rCurrNorm/MPIDotProduct(pCurr, Ap, size, startIndex, endIndex, myrank, nbproc);

        // Compute xNext

        for(unsigned int i = 0; i < size; i++)
            xNext[i] = xCurr[i] + alpha * pCurr[i];

        // Compute rNext
        for(unsigned int i = 0; i < size; i++)
            rNext[i] = rCurr[i] - alpha * Ap[i];

        // Compute next rNorm
        rNextNorm = MPIDotProduct(rNext, rNext, size, startIndex, endIndex, myrank, nbproc);
        
        // Compute beta
        beta = rNextNorm/rCurrNorm;

        // Compute pNext
        for(unsigned int i = 0; i < size; i++)
            pNext[i] = rNext[i] + beta * pCurr[i];
      
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

    free(xNext);
    free(rCurr);
    free(rNext);
    free(pCurr);
    free(pNext);
    free(Ap);

    return xCurr;
}

int get_1d(int x, int y, int ySize){
    return x * ySize + y;
}

//mettre inline
void get_2d(int i, int ySize, int* x, int* y){
    *x = i/ySize;
    *y = i%ySize;
}

void setSystemMatrixElements(SparseMatrix* A, double* b, double* result, unsigned int start, unsigned int end, int xSize, int ySize, double** h, Parameters* params, double t){
    resetSparseMatrix(A);

    int uBegin = (xSize + 1) * (ySize + 1);
    int vBegin = uBegin + (xSize + 2) * (ySize + 1);
    int vEnd = vBegin + (xSize + 1) * (ySize + 2) - 1;

    // Fill A
    int x, y;
    for(int i = start; i <= end; i++){
        
        // eta constraint
        if(i < uBegin){
            get_2d(i, ySize + 1, &x, &y);

            sparseInsertElement(A, i, i, 1.0/params->deltaT);

            sparseInsertElement(A, i, i + uBegin + ySize + 1, h[2*x+1][2*y]/params->deltaX);
            if(x == 0)
                sparseInsertElement(A, i, i + uBegin, -h[0][2*y]/params->deltaX);
            else
                sparseInsertElement(A, i, i + uBegin, -h[2*x-1][2*y]/params->deltaX);

            sparseInsertElement(A, i, i + vBegin + 1, h[2*x][2*y+1]/params->deltaY);
            if(y == 0)
                sparseInsertElement(A, i, i + uBegin, -h[2*x][0]/params->deltaY);
            else
                sparseInsertElement(A, i, i + uBegin, -h[2*x][2*y-1]/params->deltaY);

        }

        // u constraint
        else if(i < vBegin){
            get_2d(i - uBegin, ySize + 1, &x, &y);

            // At the border
            if(x == 0 || x == xSize + 1){
                sparseInsertElement(A, i, i, 1.0);
            }

            // In the middle
            else{
                sparseInsertElement(A, i, i, 1.0/params->deltaT + params->gamma);
                sparseInsertElement(A, i, get_1d(x, y, ySize + 1), params->g/params->deltaX);
                sparseInsertElement(A, i, get_1d(x - 1, y, ySize + 1), -params->g/params->deltaX);
            }
        }

        // v constraint
        else{
            get_2d(i - vBegin, ySize + 2, &x, &y);
            // At the border
            if(y == 0 || y == ySize + 1){
                sparseInsertElement(A, i, i, 1.0);
            }

            // In the middle
            else{
                sparseInsertElement(A, i, i, 1.0/params->deltaT + params->gamma);
                sparseInsertElement(A, i, get_1d(x, y, ySize + 1), params->g/params->deltaY);
                sparseInsertElement(A, i, get_1d(x, y - 1, ySize + 1), -params->g/params->deltaY);
            }
        }
    }

    // Fill b
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

void save_inputs(double* x, int xSize, int ySize){

}

int eulerImplicitMPI(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug, int debug_rank){
    assert(map);
    assert(params);

    int nbproc, myrank ;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

    /*
    fprintf(stderr, "map->a = %lf\n", map->a);
    fprintf(stderr, "map->b = %lf\n", map->b);
    */

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

    SparseMatrix* A = createSparseMatrix(inputSize, inputSize, 5*inputSize);
    if(!A){
        free(x);
        freeDoubleMatrix(h, 2 * xSize + 3, 0);
        return -1;
    }
    
    double* b = malloc(inputSize * sizeof(double));
    if(!b){
        free(x);
        freeDoubleMatrix(h, 2 * xSize + 3, 0);
        freeSparseMatrix(A);
        return -1;
    }

    for(int i = 0; i < 2 * xSize + 3; i++){
        for(int j = 0; j < 2 * ySize + 3; j++){
            h[i][j] = getGridValueAtDomainCoordinates(map, ((float)(i * xSize)/(xSize + 1)) * (params->deltaX / 2), ((float)(j * ySize)/(ySize + 1)) * (params->deltaY / 2));
        }
    }

    unsigned int start, end;

    get_indices(nbproc, myrank, inputSize, &start, &end);

    for(unsigned int t = 1; t <= params->TMax/params->deltaT; t++){
        //NULL params
        setSystemMatrixElements(A, b, x, start, end, xSize, ySize, h, NULL, t);
        free(x);
        //calculer procsize et les recvcounts et displacements
        // choper rThresh des params
        int procSize = 0;
        int* recvcounts = NULL;
        int* displs = NULL;
        double rThresh = 0.0;
        x = MPISparseConjugateGradient(A, b, inputSize, procSize, rThresh, nbproc, myrank, start, end, recvcounts, displs);

        // faire les saves
    }

    // assigner eta, u et v
}