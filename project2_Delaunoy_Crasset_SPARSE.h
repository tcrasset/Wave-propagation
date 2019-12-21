#ifndef SPARSE_H_
#define SPARSE_H_

typedef struct SparseMatrix_t SparseMatrix;
typedef struct SparseVector_t SparseVector;

SparseVector* createSparseVector(unsigned int maxNbElements);
void freeSparseVector(SparseVector* vec);
SparseMatrix* createSparseMatrix(unsigned int begin, unsigned int end, unsigned int maxNbElements);
void freeSparseMatrix(SparseMatrix* mat);
void printSparseVector(SparseVector* vec, int row);
void printSparseMatrix(SparseMatrix* mat);
double vecSparseDotProduct(SparseVector* vec1, double* vec2);
double sparseDotProduct(SparseMatrix* mat, unsigned int row, double* vector);
void sparseVecInsertElement(SparseVector* vec, unsigned int j, double elem);
void sparseInsertElement(SparseMatrix* mat, unsigned int i, unsigned int j, double elem);
void resetSparseVector(SparseVector* vec);
void resetSparseMatrix(SparseMatrix * mat);
void MPIMatVecMul(SparseMatrix* A, double* x, double* tmpBuff, double* result, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc, int* recvcounts, int* displs);
double MPIDotProduct(double* x, double* y, unsigned int size, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc);
double sparseVecVecDotProduct(SparseVector* vec1, SparseVector* vec2);
double sparseMatVecDotProduct(SparseMatrix* mat, unsigned int row, SparseVector* vec);

#endif