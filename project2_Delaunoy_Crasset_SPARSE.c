#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <omp.h>

#include "project2_Delaunoy_Crasset_SPARSE.h"

/**
 * Structure representing a sparse matrix
 */
struct SparseMatrix_t{
	SparseVector** vectors;
	unsigned int begin;
	unsigned int end;
};

/**
 * Structure representing a sparse vector
 */
struct SparseVector_t{
	double* A;
	int* indices;
	unsigned int maxNbElements;
	unsigned int currNbElement;
};

/**
 * Allocate a SparseVector structure
 *
 * Parameters:
 * maxNbElements: The maximum number of element this vector can contain
 * 
 * Returns:
 * A pointer to the allocated SparseVector structure
 */
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

/**
 * Free a SparseVector structure
 *
 * Parameters:
 * vec: A pointer to the structure to free
 */
void freeSparseVector(SparseVector* vec){
	free(vec->A);
	free(vec->indices);
	free(vec);
}

/**
 * Allocate a SparseMatrix structure
 *
 * Parameters:
 * begin: The starting row index of the matrix
 * end: The ending row index of the matrix
 * maxNbElements: The maximum number of elements each of its vectors can contain
 * 
 * Returns:
 * A pointer to the allocated structure
 */
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

/**
 * Free a SparseMatrix structure
 *
 * Parameters:
 * mat: A pointer to the structure to free
 */
void freeSparseMatrix(SparseMatrix* mat){
	for(int i = 0; i < mat->end - mat->begin + 1; i++){
		freeSparseVector(mat->vectors[i]);
	}
	free(mat->vectors);
	free(mat);
}

/**
 * Print a SparseVector structure
 *
 * Parameters:
 * vec: A pointer to the SparseVector to print
 * row: The row this vector is associated with
 */
void printSparseVector(const SparseVector* vec, int row){
	for(int i = 0; i < vec->currNbElement; i++)
		fprintf(stderr, "(%d, %d) = %lf\n", row, vec->indices[i], vec->A[i]);
}

/**
 * Print a SparseMatrix structure
 *
 * Parameters:
 * mat: A pointer to the SparseMatrix to print
 */
void printSparseMatrix(const SparseMatrix* mat){
	for(int i = mat->begin; i < mat->end + 1; i++)
		printSparseVector(mat->vectors[i - mat->begin], i);
}

/**
 * Perform dot product between a SparseVector and an array
 *
 * Parameters:
 * vec1: A pointer to the SparseVector
 * vec2: A pointer to the array
 * 
 * Returns:
 * The dot product of vec1 and vec2
 */
double vecSparseDotProduct(const SparseVector* vec1, const double* vec2){
	double result = 0.0;
	
	for(unsigned int i = 0; i < vec1->currNbElement; i++){
		result += vec1->A[i] * vec2[vec1->indices[i]];
	}

	return result;
}

/**
 * Perform dot product between a row of a SparseMatrix and a SparseVector
 *
 * Parameters:
 * mat: A pointer to the sparseMatrix
 * row: The row of the SparseMatrix to use
 * vector: A pointer to the SparseVector
 * 
 * Returns:
 * The dot product between the row 'row' of mat and vector
 */
double sparseDotProduct(const SparseMatrix* mat, unsigned int row, const double* vector){
	return vecSparseDotProduct(mat->vectors[row - mat->begin], vector);
}

/**
 * Insert an element in a sparse vector
 *
 * Parameters:
 * vec: A pointer to the SparseVector
 * j: The indices at which to insert the element
 * elem: The value to insert
 */
void sparseVecInsertElement(SparseVector* vec, unsigned int j, double elem){
	vec->currNbElement++;
	int i;

	// Loop while element is high that the one to insert such as to keep the vector sorted
	for (i = vec->currNbElement - 1; i >= 1 && vec->indices[i-1] > j; i--){
		vec->indices[i] = vec->indices[i-1];
		vec->A[i] = vec->A[i-1];
	}

	vec->indices[i] = j;
	vec->A[i] = elem;
}

/**
 * Insert an element in a sparse matrix
 *
 * Parameters:
 * mat: A pointer to a SparseMatrix
 * i: The row at which to insert the element
 * j: The column at which to insert the element
 * elem: The value to insert
 */
void sparseInsertElement(SparseMatrix* mat, unsigned int i, unsigned int j, double elem){
	sparseVecInsertElement(mat->vectors[i - mat->begin], j, elem);
}

/**
 * Reset a SparseVector removing all its elements.
 *
 * Parameters:
 * vec: A pointer to the SparseVector to reset
 */
void resetSparseVector(SparseVector* vec){
	vec->currNbElement = 0;
}

/**
 * Reset a SparseMatrix removing all its elements.
 *
 * Parameters:
 * mat: A pointer to the SparseMatrix to reset
 */
void resetSparseMatrix(SparseMatrix * mat){
	for(int i = 0; i < mat->end - mat->begin + 1; i++)
		resetSparseVector(mat->vectors[i]);
}

/**
 * Perform dot product between to vectors in parallel
 *
 * Parameters:
 * x: A pointer to the first vector
 * y: A pointer to the second vector
 * size: The size of the vectors
 * startIndex: The starting index the calling process is responsible for
 * endIndex: The ending index the calling process is responsible for
 * myrank: The rank of the calling process
 * nbproc: The number of processes
 * 
 * Returns:
 * The dot product between the row 'row' of mat and vector
 */
double MPIDotProduct(const double* x, const double* y, unsigned int size, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc){
    double myResult = 0.0;

    // Compute the dot product of the elements this process is responsible for
    // Share the computations between each thread of this process
    #pragma omp parallel reduction(+: myResult) default(shared)
	{	
		#pragma omp for schedule(static)
	    for(unsigned int i = startIndex; i <= endIndex; i++){
	        myResult += x[i] * y[i];
	    }
	}

    double totResult;

    // Sum the results of all processes
    MPI_Allreduce(&myResult, &totResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return totResult;
}

/**
 * Perform A matrix multiplication between a SparseMatrix and an array in parallel
 *
 * Parameters:
 * A: A pointer to the SparseMatrix
 * x: A pointer to the array
 * tmpBuff: A pointer to a buffer used for computation
 * result: A pointer to an array which elements will be set to the result.
 * startIndex: The starting index the calling process is responsible for
 * endIndex: The ending index the calling process is responsible for
 * myrank: The rank of the calling process
 * nbproc: The number of processes
 * 
 * Returns:
 * The dot product between the row 'row' of mat and vector
 */
void MPIMatVecMul(const SparseMatrix* A, const double* x, double* tmpBuff, double* result, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc, int* recvcounts, int* displs){

	// Compute the dot product with each line this process is responsible for
	// Share the line betweens the thread of this process
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
	    for(unsigned int i = startIndex ; i <= endIndex; i++){
	        tmpBuff[i-startIndex] = sparseDotProduct(A, i, x);
	    }
	}

	// Gather the result of all processes
    MPI_Allgatherv(tmpBuff, (endIndex - startIndex + 1), MPI_DOUBLE, result, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
}

/**
 * Perform dot product between two SparseVectors
 *
 * Parameters:
 * vec1: A pointer to the first SparseVector
 * vec2: A pointer to the second SparseVector
 * 
 * Returns:
 * The dot product between the two SparseVectors
 */
double sparseVecVecDotProduct(const SparseVector* vec1, const SparseVector* vec2){
	
	// Initialise variables
	double result = 0.0;
	int i1 = 0;
	int i2 = 0;

	// While there are still element in the two vectors
	while(i1 < vec1->currNbElement && i2 < vec2->currNbElement){
		
		// If the current index of vec1 is lower than the one of vec2, go to next index of vec1
		if(vec1->indices[i1] < vec2->indices[i2]){
			i1++;
		}
		// If the current index of both vectors are the same, multiply both values, add to result and go to next index
		else if(vec1->indices[i1] == vec2->indices[i2]){
			result += vec1->A[i1] * vec2->A[i2];
			i1++;
			i2++;
		}
		// If the current index of vec2 is lower than the one of vec1, go to next index of vec2
		else{
			i2++;
		}
	}
	return result;
}

/**
 * Perform dot product between a row of a SparseMatrix and a SparseVector
 *
 * Parameters:
 * mat: A pointer to the SparseMatrix
 * row: The row of the SparseMatrix to use
 * vec: A pointer to the SparseVector
 * 
 * Returns:
 * The dot product between the row 'row' of mat and vec
 */
double sparseMatVecDotProduct(const SparseMatrix* mat, unsigned int row, const SparseVector* vec){
	return sparseVecVecDotProduct(mat->vectors[row - mat->begin], vec);
}