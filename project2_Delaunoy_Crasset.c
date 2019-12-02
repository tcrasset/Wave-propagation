#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include "project2_Delaunoy_Crasset.h"

#define M_PI 3.14159265358979323846

void writeTestMap(char* filename, int debug){

    FILE* fp;

    fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    double a = 5;
    double b = 10;
    int X = 20;
    int Y = 10;

    assert(a > 0);
    assert(b > 0);
    assert(X > 0);
    assert(Y > 0);

    fwrite(&a, sizeof(a),1,fp);
    fwrite(&b, sizeof(b),1,fp);
    fwrite(&X, sizeof(X),1,fp);
    fwrite(&Y, sizeof(Y),1,fp);

    for (int row = 0; row < Y; row++){   
        for (int col = 0; col < X; col++){
            double value =(double) (Y - row - 1)* X + col; 
            if(debug == 1)
                printf("%lf \n", value);
            fwrite(&value, sizeof(value),1,fp); 
        }
    }    

    fclose(fp);
}


void writeResultMatrix(char* filename, int xsize, int ysize,
                         double** matrix, int debug){
    
    FILE* fp;

    fp = fopen(filename, "wb");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fwrite(&xsize, sizeof(xsize),1,fp);
    fwrite(&ysize, sizeof(ysize),1,fp);

    for(int row = ysize -1; row >= 0; row--){
        for(int col = 0; col < xsize;col++){
            if(debug == 1){
                printf("%lf \n", matrix[col][row]);
            }
            fwrite(&matrix[col][row], 8,1,fp);
        }
    }

    fclose(fp);
}

SparseMatrix* toSparseMatrix(double** matrix, int xSize, int ySize){

    for(int i = 0; i < xSize; i++){
        for (int j = 0; j < ySize; j++){
            if(matrix[i][j] != 0){

            }
        }
        
    }

}


Parameters* readParameterFile(const char* filename) {
    FILE* fp;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    Parameters* params = malloc(sizeof(Parameters));

    if (params == NULL) {
        fprintf(stderr, "Unable to allocate memory for parameters\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    fscanf(fp, "%lf", &params->g);
    fscanf(fp, "%lf", &params->gamma);
    fscanf(fp, "%lf", &params->deltaX);
    fscanf(fp, "%lf", &params->deltaY);
    fscanf(fp, "%lf", &params->deltaT);
    fscanf(fp, "%u", &params->TMax);
    fscanf(fp, "%lf", &params->A);
    fscanf(fp, "%lf", &params->f);
    fscanf(fp, "%u", &params->S);
    fscanf(fp, "%u", &params->s);
    fscanf(fp, "%lf", &params->r_threshold);

    fclose(fp);

    return params;
}

void printUsefulMapInformation(Map* map){

    printf("X: %d Y: %d\n", map->X, map->Y);
    printf("a : %lf, b : %lf \n", map->a, map->b);
    printf("Sampling steps: dx = %lf, dy = %lf\n", map->dx, map->dy);
}

double bilinearInterpolation(Map* map, double x, double y){

    assert(0 <= x);
    assert(0 <= y);
    assert(x <= map->a);
    assert(y <= map->b);

    // Sampling coordinates
    // NB: trunc comes from the math library, which is not included in base gcc, you have to 
    // compile your code using the flag -lm to add it at compile time like this:
    // gcc -g yourfile.c -lm -o yourOutFile
    // i.e. the flag should come after the c code
    int k = trunc(x/map->dx);
    int l = trunc(y/map->dy);

    // printUsefulMapInformation(map);
    // printf("x : %lf, y : %lf \n", x, y);
    // printf("k: %d, l: %d \n", k, l);

    double x_k = k * map->dx;
    double x_k1 = (k+1) * map->dx;
    double y_l = l * map->dy;
    double y_l1 = (l+1) * map->dy;

    double prod1 = (x_k1 - x) * (y_l1 - y);
    double prod2 = (x_k1 - x) * (y - y_l);
    double prod3 = (x - x_k) * (y_l1 - y);
    double prod4 = (x - x_k) * (y - y_l);

    double return_value = prod1 * map->grid[k][l];  
    double epsilon = 10e-6;

    // Robust != statement
    if(fabs(x - map->a) > epsilon)
        return_value += prod3 * map->grid[k+1][l];

    if(fabs(y - map->b) > epsilon)
        return_value += prod2 * map->grid[k][l+1];

    if(fabs(x - map->a) > epsilon && fabs(y - map->b) > epsilon)
        return_value += prod4 * map->grid[k+1][l+1];

    return_value /= map->dx*map->dy;

    return return_value;

    /*
    return (prod1 * map->grid[k][l]
            + prod2 * map->grid[k][l+1]
            + prod3 * map->grid[k+1][l]
            + prod4 * map->grid[k+1][l+1])/(map->dx*map->dy);
    */
}

double getGridValueAtDomainCoordinates(Map* map, double x, double y){
    assert(x >= 0);
    assert(y >= 0);
    double epsilon = 10e-6;
    // Sampling step

    // If value already in the grid, use that instead of interpolating
    if(fmod(x, map->dx) < epsilon && fmod(y, map->dy)  < epsilon){
        return map->grid[(int) trunc(x/map->dx)][(int) trunc(y/map->dy)]; 
    } else {
        return bilinearInterpolation(map, x, y);
    }
}

void printGrid(Map* map){
    assert(map);
    for(int i = 0; i < map->X; i++){
        for(int j = 0; j < map->Y; j++){
            printf("%lf ", map->grid[i][j]);
        }
        printf("\n");
    }
}


double** allocateDoubleMatrix(int x, int y){
    assert(x > 0);
    assert(y > 0);
    double** matrix = malloc(x * sizeof(double*));
    if(!matrix)
        return NULL;

    for(int i = 0; i < x; i++){
        matrix[i] = malloc(y * sizeof(double));
        if(!matrix[i]){
            for(int j = i-1; j >= 0; j--)
                free(matrix[j]);
            free(matrix);
            return NULL;
        }
    }

    return matrix;
}

void freeDoubleMatrix(double** matrix, int x, int debug){
    assert(x > 0);
    assert(matrix);


    for(int i = 0; i < x; i++){
        if(debug == 1){
            int nbproc, myrank ;

            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
            MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

            if(myrank == 1){

                int startval_Y = (10 * myrank) / (nbproc);
                int endval_Y = (10 * (myrank+1)) / (nbproc);

                int start = (myrank == 0) ? 0 : 2*startval_Y + 3;

                for(int j = start; j < 2 * endval_Y +3 ; j++){
                    // printf("P%d (%d,%d) = %lf\n",myrank, i,j,*(matrix[i] + j));
                    fprintf(stderr, "P%d (%d,%d) = %lf\n",myrank, i,j,*(matrix[i] + j));
                }
            }
        }

        if(matrix[i] != NULL){
            free(matrix[i]);
        }
    }
    free(matrix);
}

Map* readMapFile(const char* filename, int debug) {
    FILE* fp;
    char buffer[8];

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Read the parameters at the top of the map file
    // Assumption : there are no spaces and no end of lines, just
    // contiguous bytes in double precision (8 bytes per unit)

    Map* map = malloc(sizeof(Map));

    if (map == NULL) {
        fprintf(stderr, "Unable to allocate memory for the map\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }


    // Read constants from the map file
    fread(buffer, 8, 1, fp);
    map->a = *((double*)buffer);

    fread(buffer, 8, 1, fp);
    map->b = *((double*)buffer);

    fread(buffer, 4, 1, fp);
    map->X = *((int*)buffer);

    fread(buffer, 4, 1, fp);
    map->Y = *((int*)buffer);

    // Sampling step
    map->dx = map->a/(map->X -1);
    map->dy = map->b/(map->Y -1);

    map->grid = allocateDoubleMatrix(map->X, map->Y);

    if (map->grid == NULL) {
        fprintf(stderr, "Unable to allocate memory for the grid\n");
        free(map);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    // Read the bathymetry depth grid from the map file

    long long i = 0;
    fread(buffer, 8, 1, fp);
    for(int row = map->Y-1; row >= 0; row--){
        for(int col = 0; col < map->X;fread(buffer, 8, 1, fp), col++){
            if(debug == 1)
                printf("(%d,%d)=%lf\n",col,row,*((double*)buffer));
            map->grid[col][row] = *((double*)buffer);
        }
    }
    
    fclose(fp);
    return map;
}

void printDoubleMatrix(double** matrix, int x, int y, int process_rank){
    assert(x > 0);
    assert(y > 0);
    assert(matrix != NULL);

    for(int i = 0; i < x; i++){
        for(int j = 0; j < y; j++){
        fprintf(stderr,"\tP%d\t", process_rank);
            fprintf(stderr, "%lf ", matrix[i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void printLinearArray(double* array, int x, int y){
    for(int i = 0; i < x * y; i ++){
        fprintf(stderr,"%lf ", array[i]);
        if(i != 0 && (i+1) % (x) == 0){
            fprintf(stderr,"\n");
        }
    }
    fprintf(stderr,"\n");
}

double* transformMatrixToArray(double** matrix, int x, int y){
    double * array = calloc(x*y, sizeof(double));
   
    int cnt = 0;
    for(int i = 0; i < x; i++){
        for(int j = 0; j < y; j++){
            array[cnt++] = matrix[i][j];
        }
    }

    return array;
}


int eulerExplicitMPI(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug, int debug_rank){

    // int i = 0;
    // char hostname[256];
    // gethostname(hostname, sizeof(hostname));
    // printf("PID %d on %s ready for attach\n", getpid(), hostname);
    // fflush(stderr);
    // while (0 == i)
    // sleep(5);
    
    assert(map);
    assert(params);

    int nbproc, myrank ;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

    int xSize = (int)(map->a / params->deltaX);
    int ySize = (int)(map->b / params->deltaY);

    int mpi_xsize = xSize/nbproc;

    int startval_X, endval_X;
    int startval_X_h, endval_X_h;
    int startval_X_u, endval_X_u;
    if(myrank == 0){
        startval_X = 0;
        endval_X = mpi_xsize;
        startval_X_h = 0;
        endval_X_h = 2*mpi_xsize + 2;
        startval_X_u = 0;
        endval_X_u = mpi_xsize;
    }else if(myrank == nbproc -1){
        startval_X = myrank * mpi_xsize + 1;
        endval_X = (myrank+1) * mpi_xsize;
        // startval_X_h = 2 * myrank * mpi_xsize + 3;
        startval_X_h = 2 * myrank * mpi_xsize + 2; //Include this line so as to not transfer h
        endval_X_h = 2 * (myrank+1) * mpi_xsize + 2;
        startval_X_u = myrank * mpi_xsize + 1;
        endval_X_u = (myrank+1) * mpi_xsize;
    }else{
        startval_X = myrank * mpi_xsize + 1;
        endval_X = (myrank+1) * mpi_xsize; 
        // startval_X_h = 2 * myrank * mpi_xsize + 3;
        startval_X_h = 2 * myrank * mpi_xsize + 2;
        endval_X_h = 2 * (myrank+1) * mpi_xsize + 2;
        startval_X_u = myrank * mpi_xsize + 1;
        endval_X_u = (myrank+1) * mpi_xsize + 1;
    }

    if(debug == 1 && myrank == debug_rank){
        fprintf(stderr, "Process %d (mpi_xsize, ySize) = (%d,%d)\n", myrank, mpi_xsize, ySize);
        fprintf(stderr, "Process %d  (start,end) = (%d,%d)\n", myrank, startval_X, endval_X);
    }

    int size_X = endval_X - startval_X + 1;
    int size_X_u = endval_X_u - startval_X_u + 1;
    int size_X_h = endval_X_h - startval_X_h + 1;

    fprintf(stderr,"Processs %d \n xSize %d \n Size X u %d\n Size X_h %d\n Size_X %d \n", myrank, xSize, size_X_u, size_X_h, size_X);

    // Allocate memory
    // eta in {0, 1, ..., a/dx}X{0, 1, ..., b/dy}
    double** etaCurr = allocateDoubleMatrix(size_X, ySize + 1);
    if(!etaCurr){
        return -1;
    }

    double** etaNext = allocateDoubleMatrix(size_X, ySize + 1);
    if(!etaNext){
        freeDoubleMatrix(etaCurr, size_X,0);
        return -1;
    }

    // u in {-1/2, 1/2, ..., a/dx + 1/2}X{0, 1, ..., b/dy}
    double** uCurr = allocateDoubleMatrix(size_X_u, ySize + 1);
    if(!uCurr){
        freeDoubleMatrix(etaCurr,size_X,0);
        freeDoubleMatrix(etaNext,size_X,0);
        return -1;
    }

    double** uNext = allocateDoubleMatrix(size_X_u, ySize + 1);
    if(!uNext){
        freeDoubleMatrix(etaCurr,size_X,0);
        freeDoubleMatrix(etaNext,size_X,0);
        freeDoubleMatrix(uCurr, size_X_u,0);
        return -1;
    }

    // v in {0, 1, .., a/dx}X{-1/2, 1/2, ..., b/dy + 1/2}
    double** vCurr = allocateDoubleMatrix(size_X, ySize + 2);
    if(!vCurr){
        freeDoubleMatrix(etaCurr, size_X,0);
        freeDoubleMatrix(etaNext, size_X,0);
        freeDoubleMatrix(uCurr, size_X_u,0);
        freeDoubleMatrix(uNext, size_X_u,0);
        return -1;
    }

    double** vNext = allocateDoubleMatrix(size_X, ySize + 2);
    if(!vNext){
        freeDoubleMatrix(etaCurr, size_X,0);
        freeDoubleMatrix(etaNext, size_X,0);
        freeDoubleMatrix(uCurr, size_X_u,0);
        freeDoubleMatrix(uNext, size_X_u,0);
        freeDoubleMatrix(vCurr, size_X,0);
        return -1;
    }

    // h in {-1/2, 0, 1/2, ..., a/dx, a/dx + 1/2}X{-1/2, 0, 1/2, ..., b/dy, b/dy + 1/2}
    double** h = allocateDoubleMatrix(size_X_h, 2 * ySize + 3);
    if(!h){
        freeDoubleMatrix(etaCurr, size_X,0);
        freeDoubleMatrix(etaNext, size_X,0);
        freeDoubleMatrix(uCurr, size_X_u,0);
        freeDoubleMatrix(uNext, size_X_u,0);
        freeDoubleMatrix(vCurr, size_X,0);
        freeDoubleMatrix(vNext, size_X,0);
        return -1;
    }

    // fprintf(stderr, "%d\n",startval_X_h);
    // fprintf(stderr, "%d\n",endval_X_h);
    for(int i = startval_X_h; i <= endval_X_h; i++){
        for(int j = 0; j < 2 * ySize + 3; j++){
            h[i-startval_X_h][j] = getGridValueAtDomainCoordinates(map, ((float)(i * xSize)/(xSize + 1)) * (params->deltaX / 2), ((float)(j * ySize)/(ySize + 1)) * (params->deltaY / 2));
        }
    }
    fprintf(stderr, "Does not fail at h fill\n");

    // if(debug == 1 && myrank == debug_rank){
    //     printf("*************Process %d *******************\n", myrank);
    //     printf("h\n");
    //     for(int i = startval_X_h; i < endval_X_h; i++){
    //         for(int j = 0; j < 2 * ySize + 3; j++){
    //             fprintf(stderr, "%lf ",h[i][j]);
    //         }
    //         fprintf(stderr, "\n");
    //     }
    //     printf("apres h\n");
    // }


    for(int i = 0; i < size_X; i++){
        for(int j = 0; j < ySize; j++){
            etaCurr[i][j] = 0;
        }
    }

    for(int i = 0; i < size_X_u; i++){
        for(int j = 0; j < ySize; j++){
            uCurr[i][j] = 0;
        }
    }

    for(int i = 0; i < size_X; i++){
        for(int j = 0; j < ySize; j++)
            vCurr[i][j] = 0;
    }


    for(unsigned int t = 1; t <= params->TMax; t++){
        if(debug == 1 && myrank == debug_rank){
            printf("************* Process %d *******************\n", myrank);
            printf("t = %u\n", t);
        }

        // Compute etaNext
        // Separate for loop for cache optimization
        /*
        for(int i = 0; i < xSize + 1; i++)
            etaNext[i][0] = 0;

        for(int i = 0; i < xSize + 1; i++)
            etaNext[i][ySize] = 0;

        for(int i = 0; i < ySize + 1; i++)
            etaNext[0][i] = 0;

        for(int i = 0; i < ySize + 1; i++)
            etaNext[xSize][i] = 0;
        */

        //Receive etaCurr[mpi_xSize] or etaCurr[(myrank+1)*mpi_xsize]



        fprintf(stderr, "Fails before etaNext\n");
        /*
        double* uReceived = malloc((ySize + 1) * sizeof(double));
        fprintf(stderr, "allocated uReceived");
        if(myrank == 0){
            MPI_Send(uCurr[size_X-1], ySize + 1, MPI_DOUBLE, 1, 62, MPI_COMM_WORLD); //Tag 62 is for u
        }else if (myrank = nbproc -1){
            MPI_Recv(uReceived, ySize + 1, MPI_DOUBLE, myrank - 1, 62, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }else{
            MPI_Sendrecv(uCurr[size_X-1], ySize + 1, MPI_DOUBLE, myrank + 1, 62,
                            uReceived, ySize + 1, MPI_DOUBLE, myrank, 62, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        fprintf(stderr,"*****************PROCESS %d *****************",myrank);
        for(int j = 0; j < ySize + 1; j++){
            fprintf(stderr, "%lf ", uReceived[j]);
        }
        */

        // Process etaNext in one block
        for(int i = 0; i < size_X; i++){
            for(int j = 0; j < ySize + 1; j++){
                etaNext[i][j] = (-(h[2*i+2][2*j+1] * uCurr[i+1][j] - h[2*i][2*j+1] * uCurr[i][j]) / params->deltaX 
                                -(h[2*i+1][2*j+2] * vCurr[i][j+1] - h[2*i+1][2*j] * vCurr[i][j]) / params->deltaY)
                                * params->deltaT + etaCurr[i][j];
            }
        }
        /*
        if(myrank == 0){ 
            
        }
        
        //Process etaNext[0] alone because of uReceived
        else { 
            for(int j = 0; j < ySize + 1; j++){
                etaNext[0][j] = (-(h[2][2*j+1] * uCurr[0][j] - h[0][2*j+1] * uReceived[j]) / params->deltaX 
                                    -(h[1][2*j+2] * vCurr[0][j+1] - h[1][2*j] * vCurr[0][j]) / params->deltaY)
                                    * params->deltaT + etaCurr[0][j];
            }

            for(int i = 1; i < size_X +1; i++){
                for(int j = 0; j < ySize + 1; j++){
                    etaNext[i][j] = (-(h[2*i+2][2*j+1] * uCurr[i+1][j] - h[2*i][2*j+1] * uCurr[i][j]) / params->deltaX 
                                    -(h[2*i+1][2*j+2] * vCurr[i][j+1] - h[2*i+1][2*j] * vCurr[i][j]) / params->deltaY)
                                    * params->deltaT + etaCurr[i][j];
                }
            }
        }
        */

        fprintf(stderr, "Does not fail at etaNext\n");

        // Compute uNext
        if(myrank == 0){
            for(int i = 0; i < ySize + 1; i++){
                uNext[0][i] = 0;
            }
        }
        else if(myrank == nbproc -1){
            for(int i = 0; i < ySize + 1; i++){
                uNext[xSize-1][i] = 0;
            }
        }

        double* etaReceived = malloc(ySize + 1 * sizeof(double));
        fprintf(stderr, "allocated etaReceived");
        if(myrank == 0){
            MPI_Send(etaCurr[size_X-1], ySize + 1, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD); //Tag 42 is for eta
        }else if (myrank = nbproc -1){
            MPI_Recv(etaReceived, ySize + 1, MPI_DOUBLE, myrank - 1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }else{
            MPI_Sendrecv(etaCurr[size_X-1], ySize + 1, MPI_DOUBLE, myrank + 1, 42,
                            etaReceived, ySize + 1, MPI_DOUBLE, myrank, 42,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        fprintf(stderr,"*****************PROCESS %d *****************",myrank);
        for(int j = 0; j < ySize + 1; j++){
            fprintf(stderr, "%lf ", etaReceived[j]);
        }

        // Process uNext in one block
        if(myrank == 0){
            for(int i = 1; i < size_X_u; i++){ // Shouldn't that be size_X_u ? Or is it because one starts at 1 and not 0
                for(int j = 0; j < ySize + 1; j++){
                    uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                }
            }
        }
        else if(myrank == nbproc-1){
            for(int j = 0; j < ySize + 1; j++){
                uNext[0][j] = (-params->g * (etaCurr[0][j] - etaReceived[j]) / params->deltaX
                               -params->gamma * uCurr[0][j]) * params->deltaT + uCurr[0][j];
            }
            for(int i = 1; i < size_X_u-1; i++){ // Shouldn't that be size_X_u ? Or is it because one starts at 1 and not 0
                for(int j = 0; j < ySize + 1; j++){
                    uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                }
            }
        }

        // Process uNext[0] alone because of etaReceived in one block
        else{
            
            for(int j = 0; j < ySize + 1; j++){
                uNext[0][j] = (-params->g * (etaCurr[0][j] - etaReceived[j]) / params->deltaX
                               -params->gamma * uCurr[0][j]) * params->deltaT + uCurr[0][j];
            }
            for(int i = 1; i < size_X_u; i++){ // Shouldn't that be size_X_u ? Or is it because one starts at 1 and not 0
                for(int j = 0; j < ySize + 1; j++){
                    uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                }
            }
        }

        for(int i = 0; i < size_X; i++)
            vNext[i][0] = 0;

        for(int i = 0; i < size_X; i++){
            if(params->s == 0)
                vNext[i][ySize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT);
            else
                vNext[i][ySize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT) * exp(- t * params->deltaT / 500);
        }
        
        for(int i = 0; i < size_X+1; i++){
            for(int j = 1; j < ySize + 1; j++){
                vNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i][j-1]) / params->deltaY
                               -params->gamma * vCurr[i][j]) * params->deltaT + vCurr[i][j];
            }
        }

        if(debug == 1 && myrank == debug_rank){
            printf("\n\n\n*************Process %d *******************\n\n\n\n", myrank);
            printf("etaCurr\n");
            printDoubleMatrix(etaCurr, size_X+1, ySize + 1, myrank);
            printf("etaNext\n");
            printDoubleMatrix(etaNext, size_X+1, ySize + 1,myrank);
            printf("uCurr\n");
            printDoubleMatrix(uCurr, size_X_u+1, ySize + 1,myrank);
            printf("uNext\n");
            printDoubleMatrix(uNext, size_X_u+1, ySize + 1,myrank);
            printf("vCurr\n");
            printDoubleMatrix(vCurr, size_X+1, ySize + 2,myrank);
            printf("vNext\n");
            printDoubleMatrix(vNext, size_X+1, ySize + 2,myrank);
        }

        // Go to next step
        double** tmp;
        
        tmp = etaCurr;
        etaCurr = etaNext;
        etaNext = tmp;

        tmp = uCurr;
        uCurr = uNext;
        uNext = tmp;

        tmp = vCurr;
        vCurr = vNext;
        vNext = tmp;

    }

    *eta = etaCurr;
    *u = uCurr;
    *v = vCurr;
    
    freeDoubleMatrix(etaNext, size_X+1,0);
    freeDoubleMatrix(uNext, size_X_u+1,0);
    freeDoubleMatrix(vNext, size_X+1,0);
    // freeDoubleMatrix(h, size_X_h,1);

    if(myrank == 1) {

        // printf("*************Process %d *******************\n", myrank);
        // printf("h\n");
        // int row = (myrank == 0) ? 0 : 2*startval_Y + 3;
        // for(; row < 2 * endval_Y + 3 ; row++){
        //     for(int col = 0; col <  2 * xSize + 3; col++){
        //         fprintf(stderr, "%lf ",h[col][row]);
        //     }
        //     fprintf(stderr, "\n");
        // }
        // printf("apres h\n");

        // free(h[0]);
        // free(h[1]);
        // free(h[2]);
        // free(h[3]);
        // free(h[4]);
        // free(h[5]);
        // free(h[6]);
        // free(h[7]);
        // free(h[8]);
        // free(h[9]);
        // free(h[10]);
        // free(h[11]);
        // free(h[12]);
        // free(h[13]);
    }

    return 0;
}

int main(int argc, char* argv[]) {
    // Check number of arguments
    int debug, debug_rank;
    assert(argc >= 4);
    const char* parameter_file = argv[1];
    const char* map_file = argv[2];
    const unsigned int scheme = atoi(argv[3]);

    if(argc == 4)
        debug = 0;
    else if(argc == 5){
        fprintf(stderr, "Call with rank of process to debug as 6th argument \n");
        return -1;
    }else{
        assert(argc == 6);
        debug = atoi(argv[4]);
        debug_rank = atoi(argv[5]);
    }

    //Check argument validity
    assert((scheme == 0) || (scheme == 1));

    MPI_Init(&argc,&argv);
    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // writeTestMap("test_map.dat", debug);

    Parameters* params = readParameterFile(parameter_file);
    // Map* map = readMapFile("serverFiles/refraction.dat", 0);
    Map* map = readMapFile("test_map.dat", 0);
    if(debug == 1)
        printDoubleMatrix(map->grid, map->X, map->Y, myrank);

    // Explicit
    if (scheme == 0) {
        double** eta;
        double** u;
        double** v;

        // if(eulerExplicit(map, params, &eta, &u, &v, debug) == -1){
        //     fprintf(stderr, "error in euler function\n");
        //     free(params);
        //     free(map->grid);
        //     free(map);
        //     //exit(EXIT_FAILURE)
        // }


        if(eulerExplicitMPI(map, params, &eta, &u, &v, debug, debug_rank) == -1){
            fprintf(stderr, "error in euler function\n");
            free(params);
            free(map->grid);
            free(map);
            //exit(EXIT_FAILURE)
        }
        

        int nbproc, myrank ;

        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

        int xSize = (int)(map->a / params->deltaX);
        int ySize = (int)(map->b / params->deltaY);

        int mpi_xsize = xSize/nbproc;

        int startval_X, endval_X;
        int startval_X_h, endval_X_h;
        int startval_X_u, endval_X_u;
        if(myrank == 0){
            startval_X = 0;
            endval_X = mpi_xsize; 
            startval_X_h = 0;
            endval_X_h = 2*mpi_xsize + 3;
            startval_X_u = 0;
            endval_X_u = mpi_xsize + 1;
        }else if(myrank == nbproc -1 && nbproc != 2){
            startval_X = myrank * mpi_xsize + 1;
            endval_X = (myrank+1) * mpi_xsize;
            // startval_X_h = 2 * myrank * mpi_xsize + 3;
            startval_X_h = 2 * myrank * mpi_xsize + 2; //Include this line so as to not transfer h
            endval_X_h = 2 * (myrank+1) * mpi_xsize + 2;
            startval_X_u = myrank * mpi_xsize + 1;
            endval_X_u = (myrank+1) * mpi_xsize;
        }else{
            startval_X = myrank * mpi_xsize + 1;
            endval_X = (myrank+1) * mpi_xsize; 
            // startval_X_h = 2 * myrank * mpi_xsize + 3;
            startval_X_h = 2 * myrank * mpi_xsize + 2;
            endval_X_h = 2 * (myrank+1) * mpi_xsize + 3;
            startval_X_u = myrank * mpi_xsize + 1;
            endval_X_u = (myrank+1) * mpi_xsize + 1;
        }

        int size_X = endval_X - startval_X;
        int size_X_u = endval_X_u - startval_X_u;
        int size_X_h = endval_X_h - startval_X_h;


        //MPI
        if(debug == 0){
            fprintf(stderr,"***********ETA**************\n");
            printDoubleMatrix(eta, size_X, ySize + 1,myrank);
            fprintf(stderr,"***********U**************\n");
            printDoubleMatrix(u, size_X_u, ySize + 1,myrank);
            fprintf(stderr,"***********V**************\n");
            printDoubleMatrix(v, size_X, ySize + 2,myrank);
        }

        double* etaPartial = transformMatrixToArray(eta, xSize+1, ySize +1);
        double* uPartial = transformMatrixToArray(u, xSize+2, ySize +1);
        double* vPartial = transformMatrixToArray(v, xSize+1, ySize +2);

        double * etaTotal = malloc((xSize + 1) * (ySize  + 1)* sizeof(double));
        double * uTotal = malloc((size_X_u) * (ySize  + 1)* sizeof(double));
        double * vTotal = malloc((size_X_u) * (ySize  + 2)* sizeof(double)); // had to increase xsize by 1 to not get segfault
        MPI_Gather(etaPartial, (size_X) * (ySize + 1) , MPI_DOUBLE, etaTotal, (size_X) * (ySize  + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(uPartial, (size_X_u) * (ySize + 1) , MPI_DOUBLE, uTotal, (size_X_u) * (ySize  + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(vPartial, (size_X) * (ySize + 2) , MPI_DOUBLE, vTotal, (size_X) * (ySize  + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(myrank == 0){
            fprintf(stderr,"***********ETA TOTAL**************\n");
            printLinearArray(etaTotal, size_X, ySize +1);
            fprintf(stderr,"***********U TOTAL**************\n");
            printLinearArray(uTotal, size_X_u, ySize +1);
            fprintf(stderr,"***********V TOTAL**************\n"); // had to increase xsize by 1 to not get segfault
            printLinearArray(vTotal, size_X_u, ySize +2);
        }
        free(etaTotal);
        free(uTotal);
        free(vTotal);

        // writeResultMatrix("eta_test.dat", xSize+1, ySize+1, eta, debug);

        freeDoubleMatrix(eta, size_X, 0);
        freeDoubleMatrix(u, size_X_u, 0);
        freeDoubleMatrix(v, size_X, 0);

    }
    // Implicit
    else {
        printf("Implicit ");
        printf("%s %s %u", parameter_file, map_file, scheme);
    }

    MPI_Finalize();

    free(params);
    freeDoubleMatrix(map->grid, map->X,0);
    free(map);
    /* code */
    return 0;
}
