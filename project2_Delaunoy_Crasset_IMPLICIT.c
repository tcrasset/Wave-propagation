#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "project2_Delaunoy_Crasset_IMPLICIT.h"
#include "project2_Delaunoy_Crasset_SPARSE.h"

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