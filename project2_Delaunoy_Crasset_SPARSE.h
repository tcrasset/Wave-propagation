#ifndef SPARSE_H_
#define SPARSE_H_

typedef struct SparseMatrix_t SparseMatrix;

SparseMatrix* createSparseMatrix(unsigned int x, unsigned int y, unsigned int maxNbElements);
void freeSparseMatrix(SparseMatrix* mat);
void sparseInsertElement(SparseMatrix* mat, unsigned int i, unsigned int j, double elem);
double sparseDotProduct(SparseMatrix* mat, unsigned int row, double* vector);

#endif