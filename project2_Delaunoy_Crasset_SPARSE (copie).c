#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

#include "project2_Delaunoy_Crasset_SPARSE.h"

struct SparseMatrix_t{
	double* A;
	unsigned int* IA;
	unsigned int* JA;
	unsigned int x;
	unsigned int y;
	unsigned int maxNbElements;
	unsigned int currNbElement;
};

struct SparseVector_t{
	double* A;
	double* index;
	unsigned int maxNbElements;
	unsigned int currNbElement;
}

SparseMatrix* createSparseMatrix(unsigned int x, unsigned int y, unsigned int maxNbElements){
	SparseMatrix* mat = malloc(sizeof(SparseMatrix));
	if(!mat)
		return NULL;

	mat->A = malloc(maxNbElements * sizeof(double));
	if(!mat->A){
		free(mat);
		return NULL;
	}

	mat->IA = calloc(x + 1, sizeof(unsigned int));
	if(!mat->IA){
		free(mat->A);
		free(mat);
		return NULL;
	}

	mat->JA = malloc(maxNbElements * sizeof(unsigned int));
	if(!mat->JA){
		free(mat->IA);
		free(mat->A);
		free(mat);
		return NULL;
	}

	mat->x = x;
	mat->y = y;
	mat->maxNbElements = maxNbElements;
	mat->currNbElement = 0;

	return mat;
}

void freeSparseMatrix(SparseMatrix* mat){
	free(mat->A);
	free(mat->IA);
	free(mat->JA);
	free(mat);
}

void printSparseMatrix(SparseMatrix* mat){
	for(int i = 0; i < mat->currNbElement; i++){
		int j = 0;
		for(j = 0; j < mat->x + 1 && mat->IA[j] <= i; j++);
		fprintf(stderr, "(%d, %d) = %lf\n", j-1, mat->JA[i], mat->A[i]);
	}
}

double sparseDotProduct(SparseMatrix* mat, unsigned int row, double* vector){
	double result = 0.0;

	for(unsigned int i = mat->IA[row]; i < mat->IA[row+1]; i++){
		result += mat->A[i] * vector[mat->JA[i]];
	}

	return result;
}

void sparseInsertElement(SparseMatrix* mat, unsigned int i, unsigned int j, double elem){

	mat->A[mat->currNbElement] = elem;

	for(unsigned int k = i+1; k < mat->x + 1; k++)
		mat->IA[k]++;

	mat->JA[mat->currNbElement] = j;

	mat->currNbElement++;
}

void resetSparseMatrix(SparseMatrix * mat){
	mat->currNbElement = 0;

	for(unsigned int i = 0; i < mat->x + 1; i++){
		mat->IA[i] = 0;
	}
}

double MPIDotProduct(double* x, double* y, unsigned int size, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc){
    double myResult = 0.0;
    for(unsigned int i = startIndex; i <= endIndex; i++){
        myResult += x[i] * y[i];
    }

    double totResult;

    MPI_Allreduce(&myResult, &totResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return totResult;
}

void MPIMatVecMul(SparseMatrix* A, double* x, double* tmpBuff, double* result, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc, int* recvcounts, int* displs){

    for(unsigned int i = startIndex ; i <= endIndex; i++){
        tmpBuff[i-startIndex] = sparseDotProduct(A, i, x);
    }

    MPI_Allgatherv(tmpBuff, (endIndex - startIndex + 1), MPI_DOUBLE, result, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
}

void MPIMatSparseVecMul(SparseMatrix* A, SparseMatrix* vec, unsigned int vecLine, double* tmpBuff, double* result, unsigned int size, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc, int* recvcounts, int* displs){
	
}

void MPISelfMatMul(SparseMatrix* A, SparseMatrix* result, unsigned int size, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc, int* recvcounts, int* displs){
	
}
