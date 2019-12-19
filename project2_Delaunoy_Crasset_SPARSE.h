#ifndef SPARSE_H_
#define SPARSE_H_

typedef struct SparseMatrix_t SparseMatrix;

SparseMatrix* createSparseMatrix(unsigned int x, unsigned int y, unsigned int maxNbElements);
void freeSparseMatrix(SparseMatrix* mat);
void printSparseMatrix(SparseMatrix* mat);
void sparseInsertElement(SparseMatrix* mat, unsigned int i, unsigned int j, double elem);
double sparseDotProduct(SparseMatrix* mat, unsigned int row, double* vector);
void resetSparseMatrix(SparseMatrix * mat);
void MPIMatVecMul(SparseMatrix* A, double* x, double* tmpBuff, double* result, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc, int* recvcounts, int* displs);
double MPIDotProduct(double* x, double* y, unsigned int size, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc);

#endif