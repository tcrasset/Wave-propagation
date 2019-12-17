#include <stdlib.h>
#include <stdio.h>

#include "project2_Delaunoy_Crasset_IMPLICIT.h"

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

    double* tmp = malloc(size * sizeof(double));
    if(!tmp){
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

    double rNorm = vectorNorm(rCurr, size);

    // Initialise p
    for(unsigned int i = 0; i < size; i++){
        pCurr[i] = rCurr[i];
    }

    double num;
    double denum;
    
    while(vectorNorm(rCurr, size) / rNorm >= rThresh){
        fprintf(stderr, "loop \n");
    
        // Compute alpha
        num = dotProduct(rCurr, rCurr, size);

        fprintf(stderr, "pCurr = ");
        for(unsigned int i = 0; i < size; i++)
            fprintf(stderr, "%lf ", pCurr[i]);
        fprintf(stderr, "\n");

        fprintf(stderr, "tmp =");
        for(unsigned int i = 0; i < size; i++){
            tmp[i] = dotProduct(A[i], pCurr, size);
            fprintf(stderr, " %lf", tmp[i]);
        }
        fprintf(stderr, "\n");

        denum = dotProduct(pCurr, tmp, size);
        fprintf(stderr, "denum = %lf\n", denum);
        double alpha = num/denum;

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
            rNext[i] = rCurr[i] - alpha * dotProduct(A[i], pCurr, size);
            fprintf(stderr, " %lf", rNext[i]);
        }
        fprintf(stderr, "\n");

        // Compute beta
        double beta = dotProduct(rNext, rNext, size)/dotProduct(rCurr, rCurr, size);
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
    }

    free(xNext);
    free(rCurr);
    free(rNext);
    free(pCurr);
    free(pNext);
    free(tmp);

    return xCurr;
}