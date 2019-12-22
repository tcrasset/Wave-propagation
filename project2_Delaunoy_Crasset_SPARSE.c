#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

#include "project2_Delaunoy_Crasset_SPARSE.h"

struct SparseMatrix_t{
	SparseVector** vectors;
	unsigned int begin;
	unsigned int end;
};

struct SparseVector_t{
	double* A;
	int* indices;
	unsigned int maxNbElements;
	unsigned int currNbElement;
};

SparseVector* createSparseVector(unsigned int maxNbElements){
	SparseVector* vec = malloc(sizeof(SparseVector));
	if(!vec)
		return NULL;

	vec->A = malloc(maxNbElements * sizeof(double));
	if(!vec->A){
		free(vec);
		return NULL;
	}

	vec->indices = malloc(maxNbElements * sizeof(int));
	if(!vec->indices){
		free(vec->A);
		free(vec);
		return NULL;
	}

	vec->maxNbElements = maxNbElements;
	vec->currNbElement = 0;

	return vec;
}

void freeSparseVector(SparseVector* vec){
	free(vec->A);
	free(vec->indices);
	free(vec);
}

SparseMatrix* createSparseMatrix(unsigned int begin, unsigned int end, unsigned int maxNbElements){
	SparseMatrix* mat = malloc(sizeof(SparseMatrix));
	if(!mat)
		return NULL;

	mat->vectors = malloc((end-begin+1) * sizeof(SparseVector*));
	if(!mat->vectors){
		free(mat);
		return NULL;
	}

	for(unsigned int i = 0; i < end-begin+1; i++){
		mat->vectors[i] = createSparseVector(maxNbElements);
		if(!mat->vectors[i]){
			for(int j = i-1; j >= 0; j--){
				freeSparseVector(mat->vectors[j]);
				free(mat->vectors);
				free(mat);
				return NULL;
			}
		}
	}

	mat->begin = begin;
	mat->end = end;

	return mat;
}

void freeSparseMatrix(SparseMatrix* mat){
	for(int i = 0; i < mat->end - mat->begin + 1; i++){
		freeSparseVector(mat->vectors[i]);
	}
	free(mat->vectors);
	free(mat);
}

void printSparseVector(SparseVector* vec, int row){
	for(int i = 0; i < vec->currNbElement; i++)
		fprintf(stderr, "(%d, %d) = %lf\n", row, vec->indices[i], vec->A[i]);
}

void printSparseMatrix(SparseMatrix* mat){
	for(int i = mat->begin; i < mat->end + 1; i++)
		printSparseVector(mat->vectors[i - mat->begin], i);
}

double vecSparseDotProduct(SparseVector* vec1, double* vec2){
	double result = 0.0;
	#pragma omp parallel reduction(+: result) private(vec1, vec2) default(shared)
	{
		#pragma omp for schedule(static)
		for(unsigned int i = 0; i < vec1->currNbElement; i++){
			result += vec1->A[i] * vec2[vec1->indices[i]];
		}
	}
	return result;
}

double sparseDotProduct(SparseMatrix* mat, unsigned int row, double* vector){
	return vecSparseDotProduct(mat->vectors[row - mat->begin], vector);
}

void sparseVecInsertElement(SparseVector* vec, unsigned int j, double elem){
	vec->A[vec->currNbElement] = elem;
	vec->indices[vec->currNbElement] = j;

	vec->currNbElement++;
}

void sparseInsertElement(SparseMatrix* mat, unsigned int i, unsigned int j, double elem){
	sparseVecInsertElement(mat->vectors[i - mat->begin], j, elem);
}

void resetSparseVector(SparseVector* vec){
	vec->currNbElement = 0;
}

void resetSparseMatrix(SparseMatrix * mat){
	for(int i = 0; i < mat->end - mat->begin + 1; i++)
		resetSparseVector(mat->vectors[i]);
}

double MPIDotProduct(double* x, double* y, unsigned int size, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc){
    double myResult = 0.0;

    #pragma omp parallel reduction(+: myResult) private(x, y) default(shared)
	{	
		#pragma omp for schedule(static)
	    for(unsigned int i = startIndex; i <= endIndex; i++){
	        myResult += x[i] * y[i];
	    }
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

double sparseVecVecDotProduct(SparseVector* vec1, SparseVector* vec2){
	double result = 0.0;
	int i1 = 0;
	int i2 = 0;
	while(i1 < vec1->currNbElement && i2 < vec2->currNbElement){
		if(vec1->indices[i1] < vec2->indices[i2]){
			i1++;
		}
		else if(vec1->indices[i1] == vec2->indices[i2]){
			result += vec1->A[i1] * vec2->A[i2];
			i1++;
			i2++;
		}
		else{
			i2++;
		}
	}
	return result;
}

double sparseMatVecDotProduct(SparseMatrix* mat, unsigned int row, SparseVector* vec){
	return sparseVecVecDotProduct(mat->vectors[row - mat->begin], vec);
}