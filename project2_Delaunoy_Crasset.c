#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
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

void freeDoubleMatrix(double** matrix, int x){
    assert(x > 0);
    assert(matrix);
    for(int i = 0; i < x; i++)
        free(matrix[i]);
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

void printDoubleMatrix(double** matrix, int x, int y){
    assert(x > 0);
    assert(y > 0);
    assert(matrix != NULL);

    for(int i = 0; i < y; i++){
        for(int j = 0; j < x; j++){
            fprintf(stdout, "%lf ", matrix[j][i]);
        }
        fprintf(stdout, "\n");
    }
}

int eulerExplicit(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug){

    assert(map);
    assert(params);

    int xSize = (int)(map->a / params->deltaX);
    int ySize = (int)(map->b / params->deltaY);

    if(debug == 1){
        fprintf(stdout, "xSize = %d ySize = %d \n", xSize, ySize);
    }

    // Allocate memory
    // eta in {0, 1, ..., a/dx}X{0, 1, ..., b/dy}
    double** etaCurr = allocateDoubleMatrix(xSize + 1, ySize + 1);
    if(!etaCurr){
        return -1;
    }

    double** etaNext = allocateDoubleMatrix(xSize + 1, ySize + 1);
    if(!etaNext){
        freeDoubleMatrix(etaCurr, xSize + 1);
        return -1;
    }

    // u in {-1/2, 1/2, ..., a/dx + 1/2}X{0, 1, ..., b/dy}
    double** uCurr = allocateDoubleMatrix(xSize + 2, ySize + 1);
    if(!uCurr){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        return -1;
    }

    double** uNext = allocateDoubleMatrix(xSize + 2, ySize + 1);
    if(!uNext){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        return -1;
    }

    // v in {0, 1, .., a/dx}X{-1/2, 1/2, ..., b/dy + 1/2}
    double** vCurr = allocateDoubleMatrix(xSize + 1, ySize + 2);
    if(!vCurr){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        freeDoubleMatrix(uNext, xSize + 2);
        return -1;
    }

    double** vNext = allocateDoubleMatrix(xSize + 1, ySize + 2);
    if(!vNext){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        freeDoubleMatrix(uNext, xSize + 2);
        freeDoubleMatrix(vCurr, xSize + 1);
        return -1;
    }

    // h in {-1/2, 0, 1/2, ..., a/dx, a/dx + 1/2}X{-1/2, 0, 1/2, ..., b/dy, b/dy + 1/2}
    double** h = allocateDoubleMatrix(2 * xSize + 3, 2 * ySize + 3);
    if(!h){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        freeDoubleMatrix(uNext, xSize + 2);
        freeDoubleMatrix(vCurr, xSize + 1);
        freeDoubleMatrix(vNext, xSize + 1);
        return -1;
    }

    // Initialise matrices

    for(int i = 0; i < 2 * xSize + 3; i++){
        for(int j = 0; j < 2 * ySize + 3; j++)
            h[i][j] = getGridValueAtDomainCoordinates(map, ((float)(i * xSize)/(xSize + 1)) * (params->deltaX / 2), ((float)(j * ySize)/(ySize + 1)) * (params->deltaY / 2));
    }

    if(debug == 1){
        printf("h\n");
        printDoubleMatrix(h, 2*xSize+3, 2*ySize+3);
        printf("apres h\n");
    }

    for(int i = 0; i < xSize + 1; i++){
        for(int j = 0; j < ySize + 1; j++)
            etaCurr[i][j] = 0;
    }

    for(int i = 0; i < xSize + 2; i++){
        for(int j = 0; j < ySize + 1; j++)
            uCurr[i][j] = 0;
    }

    for(int i = 0; i < xSize + 1; i++){
        for(int j = 0; j < ySize + 1; j++)
            vCurr[i][j] = 0;
    }

    for(unsigned int t = 1; t <= params->TMax; t++){
        if(debug == 1){   
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

        // For h, 1 step in theory is 2 step in the matrix (because we go up by 1/2)
        // Moreover, full index j in theory are 2*j+1 in the matrix
        // Thus, j - 1/2 in theory is 2*j and
        // j + 1/2 in theory is 2*j + 2
        // Except that in theory, the indices start at -1/2
        // and here they start at 0
        // Example: [... j - 1/2,   j,  j + 1/2] in theory
        //          [... 2j,     2j + 1,   2j + 2]

        // For u and v, 1 step in theory is 1 here in the matrix
        // Except that in theory, the indices start at -1/2
        // and here they start at 0.
        // As u is in the x direction, we have indices j that
        // Example: [... i - 1/2,   /,  i + 1/2] in theory
        //          [... i,     /,   i + 1]
        for(int i = 0; i < xSize + 1; i++){
            for(int j = 0; j < ySize + 1; j++){
                etaNext[i][j] = (-(h[2*i+2][2*j+1] * uCurr[i+1][j] - h[2*i][2*j+1] * uCurr[i][j]) / params->deltaX 
                                -(h[2*i+1][2*j+2] * vCurr[i][j+1] - h[2*i+1][2*j] * vCurr[i][j]) / params->deltaY)
                                * params->deltaT + etaCurr[i][j];
            }
        }

        // Compute uNext
        for(int i = 0; i < ySize + 1; i++)
            uNext[0][i] = 0;

        for(int i = 0; i < ySize + 1; i++)
            uNext[xSize+1][i] = 0;

        for(int i = 1; i < xSize + 1; i++){
            for(int j = 0; j < ySize + 1; j++){
                uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                               -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
            }
        }

        // Compute vNext
        for(int i = 0; i < xSize + 1; i++)
            vNext[i][0] = 0;

        for(int i = 0; i < xSize + 1; i++){
            if(params->s == 0)
                vNext[i][ySize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT);
            else
                vNext[i][ySize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT) * exp(- t * params->deltaT / 500);
        }

        for(int i = 0; i < xSize + 1; i++){
            for(int j = 1; j < ySize + 1; j++){
                vNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i][j-1]) / params->deltaY
                               -params->gamma * vCurr[i][j]) * params->deltaT + vCurr[i][j];
            }
        }

        if(debug == 1){
            printf("etaCurr\n");
            printDoubleMatrix(etaCurr, xSize + 1, ySize + 1);
            printf("etaNext\n");
            printDoubleMatrix(etaNext, xSize + 1, ySize + 1);
            printf("uCurr\n");
            printDoubleMatrix(uCurr, xSize + 2, ySize + 1);
            printf("uNext\n");
            printDoubleMatrix(uNext, xSize + 2, ySize + 1);
            printf("vCurr\n");
            printDoubleMatrix(vCurr, xSize + 1, ySize + 2);
            printf("vNext\n");
            printDoubleMatrix(vNext, xSize + 1, ySize + 2);
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
    
    freeDoubleMatrix(etaNext, xSize + 1);
    freeDoubleMatrix(uNext, xSize + 2);
    freeDoubleMatrix(vNext, xSize + 1);
    freeDoubleMatrix(h, 2 * xSize + 3);

    return 0;
}

int eulerExplicitMPI(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug, int debug_rank){

    assert(map);
    assert(params);

    int nbproc, myrank ;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

    int xSize = (int)(map->a / params->deltaX);
    int ySize = (int)(map->b / params->deltaY);

    int mpi_ysize = ySize /nbproc;
    int startval_Y = (ySize * myrank) / (nbproc);
    int endval_Y = (ySize * (myrank+1)) / (nbproc);

    if(debug == 1 && myrank == debug_rank){
        fprintf(stdout, "Process %d :xSize = %d mpi_ysize = %d, endval_y - startval_Y = %d  \n", myrank, xSize, mpi_ysize, endval_Y - startval_Y);
        fprintf(stdout, "Process %d :endval_Y = %d startval_Y = %d\n", myrank, endval_Y, startval_Y);
    }

    // Allocate memory
    // eta in {0, 1, ..., a/dx}X{0, 1, ..., b/dy}
    double** etaCurr = allocateDoubleMatrix(xSize + 1, mpi_ysize + 1);
    if(!etaCurr){
        return -1;
    }

    double** etaNext = allocateDoubleMatrix(xSize + 1, mpi_ysize + 1);
    if(!etaNext){
        freeDoubleMatrix(etaCurr, xSize + 1);
        return -1;
    }

    // u in {-1/2, 1/2, ..., a/dx + 1/2}X{0, 1, ..., b/dy}
    double** uCurr = allocateDoubleMatrix(xSize + 2, mpi_ysize + 1);
    if(!uCurr){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        return -1;
    }

    double** uNext = allocateDoubleMatrix(xSize + 2, mpi_ysize + 1);
    if(!uNext){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        return -1;
    }

    // v in {0, 1, .., a/dx}X{-1/2, 1/2, ..., b/dy + 1/2}
    double** vCurr = allocateDoubleMatrix(xSize + 1, mpi_ysize + 2);
    if(!vCurr){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        freeDoubleMatrix(uNext, xSize + 2);
        return -1;
    }

    double** vNext = allocateDoubleMatrix(xSize + 1, mpi_ysize + 2);
    if(!vNext){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        freeDoubleMatrix(uNext, xSize + 2);
        freeDoubleMatrix(vCurr, xSize + 1);
        return -1;
    }

    // h in {-1/2, 0, 1/2, ..., a/dx, a/dx + 1/2}X{-1/2, 0, 1/2, ..., b/dy, b/dy + 1/2}
    double** h = allocateDoubleMatrix(2 * xSize + 3, 2 * mpi_ysize + 3);
    if(!h){
        freeDoubleMatrix(etaCurr, xSize + 1);
        freeDoubleMatrix(etaNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        freeDoubleMatrix(uNext, xSize + 2);
        freeDoubleMatrix(vCurr, xSize + 1);
        freeDoubleMatrix(vNext, xSize + 1);
        return -1;
    }

    // Initialise matrices
    if(debug == 1 && myrank == debug_rank){
        printf("********Process %d **********\n",myrank);
    }
    for(int i = 0; i < 2 * xSize + 3; i++){
        int j;
        if(myrank == 0){
            j= 0;
        }
        else {
            j = 2*startval_Y + 3;
        }

        for(; j < 2 * endval_Y + 3; j++){
            h[i][j] = getGridValueAtDomainCoordinates(map, ((float)(i * xSize)/(xSize + 1)) * (params->deltaX / 2), ((float)(j * ySize)/(ySize + 1)) * (params->deltaY / 2));
        }
    }

    if(debug == 1 && myrank == debug_rank){
        printf("*************Process %d *******************\n", myrank);
        printf("h\n");
        int row;
        if(myrank == 0)row= 0;
        else row = 2*startval_Y + 3;

        if(myrank == 0)row= 0;
        else row = 2*startval_Y + 3;
        for(; row < 2 * endval_Y + 3 ; row++){
            for(int col = 0; col <  2 * xSize + 3; col++){
                // fprintf(stdout, "(col, row): (%d,%d) = %lf", col, row, h[col][row]);
                fprintf(stdout, "%lf ",h[col][row]);
            }
            fprintf(stdout, "\n");
        }
        printf("apres h\n");
    }

    return 0;

    for(int i = 0; i < xSize + 1; i++){
        for(int j = 0; j < mpi_ysize + 1; j++)
            etaCurr[i][j] = 0;
    }

    for(int i = 0; i < xSize + 2; i++){
        for(int j = 0; j < mpi_ysize + 1; j++)
            uCurr[i][j] = 0;
    }

    for(int i = 0; i < xSize + 1; i++){
        for(int j = 0; j < mpi_ysize + 1; j++)
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
            etaNext[i][mpi_ysize] = 0;

        for(int i = 0; i < mpi_ysize + 1; i++)
            etaNext[0][i] = 0;

        for(int i = 0; i < mpi_ysize + 1; i++)
            etaNext[xSize][i] = 0;
        */

        // For h, 1 step in theory is 2 step in the matrix (because we go up by 1/2)
        // Moreover, full index j in theory are 2*j+1 in the matrix
        // Thus, j - 1/2 in theory is 2*j and
        // j + 1/2 in theory is 2*j + 2
        // Except that in theory, the indices start at -1/2
        // and here they start at 0
        // Example: [... j - 1/2,   j,  j + 1/2] in theory
        //          [... 2j,     2j + 1,   2j + 2]

        // For u and v, 1 step in theory is 1 here in the matrix
        // Except that in theory, the indices start at -1/2
        // and here they start at 0.
        // As u is in the x direction, we have indices j that
        // Example: [... i - 1/2,   /,  i + 1/2] in theory
        //          [... i,     /,   i + 1]
        for(int i = 0; i < xSize + 1; i++){
            for(int j = 0; j < mpi_ysize + 1; j++){
                etaNext[i][j] = (-(h[2*i+2][2*j+1] * uCurr[i+1][j] - h[2*i][2*j+1] * uCurr[i][j]) / params->deltaX 
                                -(h[2*i+1][2*j+2] * vCurr[i][j+1] - h[2*i+1][2*j] * vCurr[i][j]) / params->deltaY)
                                * params->deltaT + etaCurr[i][j];
            }
        }

        // Compute uNext
        for(int i = 0; i < mpi_ysize + 1; i++)
            uNext[0][i] = 0;

        for(int i = 0; i < mpi_ysize + 1; i++)
            uNext[xSize+1][i] = 0;

        for(int i = 1; i < xSize + 1; i++){
            for(int j = 0; j < mpi_ysize + 1; j++){
                uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                               -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
            }
        }

        // Compute vNext
        for(int i = 0; i < xSize + 1; i++)
            vNext[i][0] = 0;

        for(int i = 0; i < xSize + 1; i++){
            if(params->s == 0)
                vNext[i][mpi_ysize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT);
            else
                vNext[i][mpi_ysize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT) * exp(- t * params->deltaT / 500);
        }

        for(int i = 0; i < xSize + 1; i++){
            for(int j = 1; j < mpi_ysize + 1; j++){
                vNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i][j-1]) / params->deltaY
                               -params->gamma * vCurr[i][j]) * params->deltaT + vCurr[i][j];
            }
        }

        if(debug == 1 && myrank == debug_rank){
            printf("\n\n\n*************Process %d *******************\n\n\n\n", myrank);
            printf("etaCurr\n");
            printDoubleMatrix(etaCurr, xSize + 1, mpi_ysize + 1);
            printf("etaNext\n");
            printDoubleMatrix(etaNext, xSize + 1, mpi_ysize + 1);
            printf("uCurr\n");
            printDoubleMatrix(uCurr, xSize + 2, mpi_ysize + 1);
            printf("uNext\n");
            printDoubleMatrix(uNext, xSize + 2, mpi_ysize + 1);
            printf("vCurr\n");
            printDoubleMatrix(vCurr, xSize + 1, mpi_ysize + 2);
            printf("vNext\n");
            printDoubleMatrix(vNext, xSize + 1, mpi_ysize + 2);
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
    
    freeDoubleMatrix(etaNext, xSize + 1);
    freeDoubleMatrix(uNext, xSize + 2);
    freeDoubleMatrix(vNext, xSize + 1);
    freeDoubleMatrix(h, 2 * xSize + 3);

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

    // writeTestMap("test_map.dat", debug);

    Parameters* param = readParameterFile(parameter_file);
    // Map* map = readMapFile("serverFiles/refraction.dat", 0);
    Map* map = readMapFile("test_map.dat", 0);
    if(debug == 1)
        printDoubleMatrix(map->grid, map->X, map->Y);

    // Explicit
    if (scheme == 0) {
        double** eta;
        double** u;
        double** v;

        // if(eulerExplicit(map, param, &eta, &u, &v, debug) == -1){
        //     fprintf(stderr, "error in euler function\n");
        //     free(param);
        //     free(map->grid);
        //     free(map);
        //     //exit(EXIT_FAILURE)
        // }


        if(eulerExplicitMPI(map, param, &eta, &u, &v, debug, debug_rank) == -1){
            fprintf(stderr, "error in euler function\n");
            free(param);
            free(map->grid);
            free(map);
            //exit(EXIT_FAILURE)
        }
        

        int xSize = (int)(map->a / param->deltaX);
        int ySize = (int)(map->b / param->deltaY);
        int nbproc;
        MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
        int mpi_ysize = (int)(map->b / param->deltaY)/nbproc;

        //MPI
        if(debug == 1){
            printf("eta\n");
            printDoubleMatrix(eta, xSize + 1, mpi_ysize + 1);
            printf("u\n");
            printDoubleMatrix(u, xSize + 2, mpi_ysize + 1);
            printf("v\n");
            printDoubleMatrix(v, xSize + 1, mpi_ysize + 2);
        }

        // Non MPI
        // if(debug == 1){
        //     printf("eta\n");
        //     printDoubleMatrix(eta, xSize + 1, ySize + 1);
        //     printf("u\n");
        //     printDoubleMatrix(u, xSize + 2, ySize + 1);
        //     printf("v\n");
        //     printDoubleMatrix(v, xSize + 1, ySize + 2);
        // }

        // writeResultMatrix("eta_test.dat", xSize+1, ySize+1, eta, debug);

        freeDoubleMatrix(eta, xSize + 1);
        freeDoubleMatrix(u, xSize + 2);
        freeDoubleMatrix(v, xSize + 1);

    }
    // Implicit
    else {
        printf("Implicit ");
        printf("%s %s %u", parameter_file, map_file, scheme);
    }

    MPI_Finalize();

    free(param);
    freeDoubleMatrix(map->grid, map->X);
    free(map);
    /* code */
    return 0;
}
