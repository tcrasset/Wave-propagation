#ifndef SPARSE_H_
#define SPARSE_H_


typedef struct SparseMatrix_t SparseMatrix;
typedef struct SparseVector_t SparseVector;

SparseVector* createSparseVector(unsigned int maxNbElements);
void freeSparseVector(SparseVector* vec);
SparseMatrix* createSparseMatrix(unsigned int begin, unsigned int end, unsigned int maxNbElements);
void freeSparseMatrix(SparseMatrix* mat);
void printSparseVector(const SparseVector* vec, int row);
void printSparseMatrix(const SparseMatrix* mat);
double vecSparseDotProduct(const SparseVector* vec1, const double* vec2);
double sparseDotProduct(const SparseMatrix* mat, unsigned int row, const double* vector);
void sparseVecInsertElement(SparseVector* vec, unsigned int j, double elem);
void sparseInsertElement(SparseMatrix* mat, unsigned int i, unsigned int j, double elem);
void resetSparseVector(SparseVector* vec);
void resetSparseMatrix(SparseMatrix * mat);
void MPIMatVecMul(const SparseMatrix* A, const double* x, double* tmpBuff, double* result, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc, int* recvcounts, int* displs);
double MPIDotProduct(const double* x, const double* y, unsigned int size, unsigned int startIndex, unsigned int endIndex, int myrank, int nbproc);
double sparseVecVecDotProduct(const SparseVector* vec1, const SparseVector* vec2);
double sparseMatVecDotProduct(const SparseMatrix* mat, unsigned int row, const SparseVector* vec);

#endif