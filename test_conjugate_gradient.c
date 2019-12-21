#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

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

void test_mpi_dot_product(int argc, char* argv[]){
	MPI_Init(&argc,&argv);

    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	double* x = malloc(5 * sizeof(double));
	x[0] = 3;
	x[1] = 2;
	x[2] = 4;
	x[3] = 1;
	x[4] = 5;

	double* y = malloc(5 * sizeof(double));
	y[0] = 6;
	y[1] = 3;
	y[2] = 8;
	y[3] = 4;
	y[4] = 1;

	double result;

	if(myrank == 0){
		result = MPIDotProduct(x, y, 5, 0, 2, myrank, nbproc);
	}

	if(myrank == 1){
		result = MPIDotProduct(x, y, 5, 3, 4, myrank, nbproc);
	}

	fprintf(stderr, "result proc %d = %lf\n", myrank, result);

	MPI_Finalize();

}

void test_mpi_mat_vec_mul(int argc, char* argv[]){
	MPI_Init(&argc,&argv);

    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	double* x = malloc(3 * sizeof(double));
	x[0] = 3.0;
	x[1] = 2.0;
	x[2] = 4.0;

	double* result = malloc(3 * sizeof(double));

	int* recvcounts = malloc(2 * sizeof(int));
	int* displs = malloc(2 * sizeof(int));

	recvcounts[0] = 2;
	recvcounts[1] = 1;

	displs[0] = 0;
	displs[1] = 2;

	if(myrank == 0){
		SparseMatrix* A = createSparseMatrix(0, 1, 5);

		sparseInsertElement(A, 0, 0, 1);
		sparseInsertElement(A, 0, 1, 2);
		sparseInsertElement(A, 0, 2, -1);
		
		sparseInsertElement(A, 1, 0, 2);
		sparseInsertElement(A, 1, 1, -1);
		sparseInsertElement(A, 1, 2, 3);

		double* tmpBuff = malloc(2 *sizeof(double));
		MPIMatVecMul(A, x, tmpBuff, result, 0, 1, myrank, nbproc, recvcounts, displs);
	}
	
	if(myrank == 1){
		SparseMatrix* A = createSparseMatrix(2, 2, 5);

		sparseInsertElement(A, 2, 0, -1);
		sparseInsertElement(A, 2, 1, 3);
		sparseInsertElement(A, 2, 2, -1);

		double* tmpBuff = malloc(sizeof(double));
		MPIMatVecMul(A, x, tmpBuff, result, 2, 2, myrank, nbproc, recvcounts, displs);
	}

	for(unsigned int i = 0; i < 3; i++){
		fprintf(stderr, "process %d result[%d] = %lf\n", myrank, i, result[i]);
	}

	MPI_Finalize();

}

void test_sparse_vec_dot_product(){
	SparseVector* vec1 = createSparseVector(10);
	SparseVector* vec2 = createSparseVector(10);

	sparseVecInsertElement(vec1, 4, 4.0);
	sparseVecInsertElement(vec1, 5, 1.0);
	sparseVecInsertElement(vec1, 9, 5.0);
	sparseVecInsertElement(vec1, 500, 2.0);
	sparseVecInsertElement(vec1, 1000, 3.0);

	sparseVecInsertElement(vec2, 1, 8.0);
	sparseVecInsertElement(vec2, 2, 9.0);
	sparseVecInsertElement(vec2, 5, 6.0);
	sparseVecInsertElement(vec2, 9, 10.0);
	sparseVecInsertElement(vec2, 500, 7.0);

	fprintf(stderr, "result = %lf\n", sparseVecVecDotProduct(vec1, vec2));
	
}

void test_sparse_mat_vec_dot_product(){
	SparseMatrix* mat = createSparseMatrix(4, 7, 10);
	SparseVector* vec = createSparseVector(10);

	sparseInsertElement(mat, 5, 4, 4.0);
	sparseInsertElement(mat, 5, 5, 1.0);
	sparseInsertElement(mat, 4, 5, 8.0);
	sparseInsertElement(mat, 5, 9, 5.0);
	sparseInsertElement(mat, 6, 9, 10.0);
	sparseInsertElement(mat, 5, 500, 2.0);
	sparseInsertElement(mat, 5, 1000, 3.0);

	sparseVecInsertElement(vec, 1, 8.0);
	sparseVecInsertElement(vec, 2, 9.0);
	sparseVecInsertElement(vec, 5, 6.0);
	sparseVecInsertElement(vec, 9, 10.0);
	sparseVecInsertElement(vec, 500, 7.0);

	fprintf(stderr, "result = %lf\n", sparseMatVecDotProduct(mat, 5, vec));
	
}

void test_mpi_conjugate(int argc, char* argv[]){
	MPI_Init(&argc,&argv);

    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int* recvcounts = malloc(2 * sizeof(int));
	int* displs = malloc(2 * sizeof(int));

	recvcounts[0] = 2;
	recvcounts[1] = 1;

	displs[0] = 0;
	displs[1] = 2;

	double* b = malloc(3*sizeof(double));
	b[0] = 2;
	b[1] = 9;
	b[2] = 2;

    if(myrank == 0){
    	SparseMatrix* A = createSparseMatrix(0, 1, 6);

		sparseInsertElement(A, 0, 0, 1);
		sparseInsertElement(A, 0, 1, 2);
		sparseInsertElement(A, 0, 2, -1);
		
		sparseInsertElement(A, 1, 0, 2);
		sparseInsertElement(A, 1, 1, -1);
		sparseInsertElement(A, 1, 2, 3);
		
		double* x = MPISparseConjugateGradient(A, b, 3, 2, 0.00001, nbproc, myrank, 0, 1, recvcounts, displs);
		fprintf(stderr, "proc %d:\n x = %lf\n y = %lf\n z = %lf\n", myrank, x[0], x[1], x[2]);
    }
    if(myrank == 1){
    	SparseMatrix* A = createSparseMatrix(2, 2, 3);
		
		sparseInsertElement(A, 2, 0, -1);
		sparseInsertElement(A, 2, 1, 3);
		sparseInsertElement(A, 2, 2, -1);

		double* x = MPISparseConjugateGradient(A, b, 3, 1, 0.00001, nbproc, myrank, 2, 2, recvcounts, displs);
		fprintf(stderr, "proc %d:\n x = %lf\n y = %lf\n z = %lf\n", myrank, x[0], x[1], x[2]);
    }

    MPI_Finalize();
    //fprintf(stderr, "x = %lf\n y = %lf\n z = %lf\n", x[0], x[1], x[2]);
}

void test_initAb(int argc, char* argv[]){
	MPI_Init(&argc,&argv);

    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double* b = malloc(16 * sizeof(double));

    double* x = malloc(16 * sizeof(double));

    for(int i = 0; i < 16; i++){
    	x[i] = i+1;
    }

    double** h = malloc(7*sizeof(double*));
    for(int i = 0; i < 7; i++){
    	h[i] = malloc(7*sizeof(double));
    }

    for(int i = 0; i < 7; i++){
    	for(int j = 0; j < 7; j++)
    		h[i][j] = i*7+j;
    }

    Parameters* params = malloc(sizeof(Parameters));
    params->deltaX = 5.0;
    params->deltaY = 10.0;
    params->deltaT = 0.1;
    params->gamma = 1.6785;
    params->f = 0.87368;
    params->s = 0;
    params->g = 9.789;
    params->A = 1;

    if(myrank == 0){
    	SparseMatrix* A = createSparseMatrix(0, 7, 5);
	    sparseInsertElement(A, 0, 0, 999999);
	    sparseInsertElement(A, 1, 12, 999999);
	    sparseInsertElement(A, 3, 5, 999999);
	    sparseInsertElement(A, 7, 7, 999999);
    	initAb(A, b, x, 0, 7, 1, 1, h, params, 1);
    	fprintf(stderr, "process 0\n");
    	fprintf(stderr, "A:\n");
    	printSparseMatrix(A);
    	fprintf(stderr, "b:");
    	for(int i = 0; i < 16; i++){
    		fprintf(stderr, " %lf", b[i]);
    	}
    	fprintf(stderr, "\n");
    }
    if(myrank == 0){
    	SparseMatrix* A = createSparseMatrix(8, 15, 5);
	    sparseInsertElement(A, 8, 5, 999999);
	    sparseInsertElement(A, 15, 3, 999999);
	    sparseInsertElement(A, 10, 5, 999999);
    	initAb(A, b, x, 8, 15, 1, 1, h, params, 1);
    	fprintf(stderr, "process 1\n");
    	fprintf(stderr, "A:\n");
    	printSparseMatrix(A);
    	fprintf(stderr, "b:");
    	for(int i = 0; i < 16; i++){
    		fprintf(stderr, " %lf", b[i]);
    	}
    	fprintf(stderr, "\n");
    }

	MPI_Finalize();	
}

void test_initAtb(int argc, char* argv[]){
	MPI_Init(&argc,&argv);

    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double* b = malloc(16 * sizeof(double));

    double* x = malloc(16 * sizeof(double));

    for(int i = 0; i < 16; i++){
    	x[i] = i+1;
    }

    double** h = malloc(7*sizeof(double*));
    for(int i = 0; i < 7; i++){
    	h[i] = malloc(7*sizeof(double));
    }

    for(int i = 0; i < 7; i++){
    	for(int j = 0; j < 7; j++)
    		h[i][j] = i*7+j;
    }

    Parameters* params = malloc(sizeof(Parameters));
    params->deltaX = 5.0;
    params->deltaY = 10.0;
    params->deltaT = 0.1;
    params->gamma = 1.6785;
    params->f = 0.87368;
    params->s = 0;
    params->g = 9.789;
    params->A = 1;

    if(myrank == 0){
    	SparseMatrix* A = createSparseMatrix(0, 7, 5);
	    sparseInsertElement(A, 0, 0, 999999);
	    sparseInsertElement(A, 1, 12, 999999);
	    sparseInsertElement(A, 3, 5, 999999);
	    sparseInsertElement(A, 7, 7, 999999);
    	initAb(A, b, x, 0, 7, 1, 1, h, params, 1);
    	fprintf(stderr, "process 0\n");
    	fprintf(stderr, "A:\n");
    	printSparseMatrix(A);
    	fprintf(stderr, "b:");
    	for(int i = 0; i < 16; i++){
    		fprintf(stderr, " %lf", b[i]);
    	}
    	fprintf(stderr, "\n");
    	initAtb(A, b, x, 0, 7, 1, 1, h, params, 1);
    	fprintf(stderr, "process 0\n");
    	fprintf(stderr, "A:\n");
    	printSparseMatrix(A);
    	fprintf(stderr, "b:");
    	for(int i = 0; i < 16; i++){
    		fprintf(stderr, " %lf", b[i]);
    	}
    	fprintf(stderr, "\n");
    }
    if(myrank == 0){
    	SparseMatrix* A = createSparseMatrix(8, 15, 5);
	    sparseInsertElement(A, 8, 5, 999999);
	    sparseInsertElement(A, 15, 3, 999999);
	    sparseInsertElement(A, 10, 5, 999999);
    	initAb(A, b, x, 8, 15, 1, 1, h, params, 1);
    	fprintf(stderr, "process 1\n");
    	fprintf(stderr, "A:\n");
    	printSparseMatrix(A);
    	fprintf(stderr, "b:");
    	for(int i = 0; i < 16; i++){
    		fprintf(stderr, " %lf", b[i]);
    	}
    	fprintf(stderr, "\n");
    	initAtb(A, b, x, 8, 15, 1, 1, h, params, 1);
    	fprintf(stderr, "process 1\n");
    	fprintf(stderr, "A:\n");
    	printSparseMatrix(A);
    	fprintf(stderr, "b:");
    	for(int i = 0; i < 16; i++){
    		fprintf(stderr, " %lf", b[i]);
    	}
    	fprintf(stderr, "\n");
    }

	MPI_Finalize();	
}

void test_build_system_matrix(int argc, char* argv[]){
	MPI_Init(&argc,&argv);

    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	Parameters* params = malloc(sizeof(Parameters));
    params->deltaX = 5.0;
    params->deltaY = 10.0;
    params->deltaT = 0.1;
    params->gamma = 1.6785;
    params->f = 0.87368;
    params->s = 0;
    params->g = 9.789;
    params->A = 1;

    double* x = malloc(16 * sizeof(double));

    for(int i = 0; i < 16; i++){
    	x[i] = i+1;
    }

    double** h = malloc(7*sizeof(double*));
    for(int i = 0; i < 7; i++){
    	h[i] = malloc(7*sizeof(double));
    }

    for(int i = 0; i < 7; i++){
    	for(int j = 0; j < 7; j++)
    		h[i][j] = i*7+j;
    }

    fprintf(stderr, "before build system\n");

    if(myrank == 0){
    	SparseMatrix* A;
		double* b;
    	BuildSystemMatrix(&A, &b, x, 0, 7, 1, 1, h, params, 1);
    	fprintf(stderr, "process 0\n");
    	fprintf(stderr, "A:\n");
    	printSparseMatrix(A);
    	fprintf(stderr, "b:");
    	for(int i = 0; i < 16; i++){
    		fprintf(stderr, " %lf", b[i]);
    	}
    }

    if(myrank == 0){
    	SparseMatrix* A;
		double* b;
    	BuildSystemMatrix(&A, &b, x, 8, 15, 1, 1, h, params, 1);
    	fprintf(stderr, "process 0\n");
    	fprintf(stderr, "A:\n");
    	printSparseMatrix(A);
    	fprintf(stderr, "b:");
    	for(int i = 0; i < 16; i++){
    		fprintf(stderr, " %lf", b[i]);
    	}
    }

    MPI_Finalize();	
}

void test_print_sparse_matrix(){
	SparseMatrix* A = createSparseMatrix(0, 2, 5);
	sparseInsertElement(A, 0, 1, 5);
	sparseInsertElement(A, 2, 0, 3);
	sparseInsertElement(A, 2, 2, 1);

	printSparseMatrix(A);
}

int main(int argc, char* argv[]){
	// /!\ MPI are to be run with 2 processes


	//test_conjugate();
	//test_sparse_conjugate();
	//test_mpi_dot_product(argc, argv);
	//test_mpi_mat_vec_mul(argc, argv);
	//test_sparse_vec_dot_product();
	//test_sparse_mat_vec_dot_product();
	//test_mpi_conjugate(argc, argv);
	//test_print_sparse_matrix();
	//test_initAb(argc, argv);
	//test_initAtb(argc, argv);
	test_build_system_matrix(argc, argv);


	return 0;
}