#ifndef SPARSE_H_
#define SPARSE_H_

/**
 * Structure representing a sparse matrix
 */
typedef struct SparseMatrix_t SparseMatrix;

/**
 * Structure representing a sparse vector
 */
typedef struct SparseVector_t SparseVector;

/**
 * Allocate a SparseVector structure
 *
 * Parameters:
 * maxNbElements: The maximum number of element this vector can contain
 * 
 * Returns:
 * A pointer to the allocated SparseVector structure
 */
SparseVector* createSparseVector(unsigned int maxNbElements);

/**
 * Free a SparseVector structure
 *
 * Parameters:
 * vec: A pointer to the structure to free
 */
void freeSparseVector(SparseVector* vec);

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
SparseMatrix* createSparseMatrix(unsigned int begin, unsigned int end, unsigned int maxNbElements);

/**
 * Free a SparseMatrix structure
 *
 * Parameters:
 * mat: A pointer to the structure to free
 */
void freeSparseMatrix(SparseMatrix* mat);

/**
 * Print a SparseVector structure
 *
 * Parameters:
 * vec: A pointer to the SparseVector to print
 * row: The row this vector is associated with
 */
void printSparseVector(const SparseVector* vec, int row);

/**
 * Print a SparseMatrix structure
 *
 * Parameters:
 * mat: A pointer to the SparseMatrix to print
 */
void printSparseMatrix(const SparseMatrix* mat);

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
double vecSparseDotProduct(const SparseVector* vec1, const double* vec2);

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
double sparseDotProduct(const SparseMatrix* mat, unsigned int row, const double* vector);

/**
 * Insert an element in a sparse vector
 *
 * Parameters:
 * vec: A pointer to the SparseVector
 * j: The indices at which to insert the element
 * elem: The value to insert
 */
void sparseVecInsertElement(SparseVector* vec, unsigned int j, double elem);

/**
 * Insert an element in a sparse matrix
 *
 * Parameters:
 * mat: A pointer to a SparseMatrix
 * i: The row at which to insert the element
 * j: The column at which to insert the element
 * elem: The value to insert
 */
void sparseInsertElement(SparseMatrix* mat, unsigned int i, unsigned int j, double elem);

/**
 * Reset a SparseVector removing all its elements.
 *
 * Parameters:
 * vec: A pointer to the SparseVector to reset
 */
void resetSparseVector(SparseVector* vec);

/**
 * Reset a SparseMatrix removing all its elements.
 *
 * Parameters:
 * mat: A pointer to the SparseMatrix to reset
 */
void resetSparseMatrix(SparseMatrix * mat);

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
double MPIDotProduct(const double* x, const double* y, unsigned int size, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc);

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
void MPIMatVecMul(const SparseMatrix* A, const double* x, double* tmpBuff, double* result, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc, int* recvcounts, int* displs);

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
double sparseVecVecDotProduct(const SparseVector* vec1, const SparseVector* vec2);

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
double sparseMatVecDotProduct(const SparseMatrix* mat, unsigned int row, const SparseVector* vec);

#endif