#include <stdlib.h>
#include <stdio.h>

#include "project2_Delaunoy_Crasset_IMPLICIT.h"
#include "project2_Delaunoy_Crasset_SPARSE.h"

void test_conjugate(){
	double** A = malloc(3*sizeof(double*));
	for(int i = 0; i < 3; i++){
		A[i] = malloc(3*sizeof(double));
	}

	double* b = malloc(3*sizeof(double));

	A[0][0] = 1;
	A[0][1] = 2;
	A[0][2] = -1;
	A[1][0] = 2;
	A[1][1] = -1;
	A[1][2] = 3;
	A[2][0] = -1;
	A[2][1] = 3;
	A[2][2] = -1;

	b[0] = 2;
	b[1] = 9;
	b[2] = 2;

	double* x = conjugateGradient(A, b, 3, 0.00001);

	fprintf(stderr, "x = %lf\n y = %lf\n z = %lf\n", x[0], x[1], x[2]);
}

void test_sparse_conjugate(){
	SparseMatrix* A = createSparseMatrix(3, 3, 9);

	double* b = malloc(3*sizeof(double));

	sparseInsertElement(A, 0, 0, 1);
	sparseInsertElement(A, 0, 1, 2);
	sparseInsertElement(A, 0, 2, -1);
	sparseInsertElement(A, 1, 0, 2);
	sparseInsertElement(A, 1, 1, -1);
	sparseInsertElement(A, 1, 2, 3);
	sparseInsertElement(A, 2, 0, -1);
	sparseInsertElement(A, 2, 1, 3);
	sparseInsertElement(A, 2, 2, -1);

	b[0] = 2;
	b[1] = 9;
	b[2] = 2;

	double* x = sparseConjugateGradient(A, b, 3, 0.00001);
	fprintf(stderr, "x = %lf\n y = %lf\n z = %lf\n", x[0], x[1], x[2]);
}

int main(){
	//test_conjugate();
	test_sparse_conjugate();

	return 0;
}