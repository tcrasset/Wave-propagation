#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>

#include "project2_Delaunoy_Crasset.h"

#define M_PI 3.14159265358979323846

void writeTestMap(char* filename){

    FILE* fp;

    fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    double a = 5;
    double b = 10;
    int X = 10;
    int Y = 20;

    fwrite(&a, sizeof(a),1,fp);
    fwrite(&b, sizeof(b),1,fp);
    fwrite(&X, sizeof(X),1,fp);
    fwrite(&Y, sizeof(Y),1,fp);

    for(int i =0; i < X * Y; i++){
        double value = (double) i;
        fwrite(&value, sizeof(value),1,fp); 
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
    double dx = map->a/map->X;
    double dy = map->b/map->Y;
    printf("X: %lld Y: %lld\n", map->X, map->Y);
    printf("a : %lf, b : %lf \n", map->a, map->b);
    printf("Sampling steps: dx = %lf, dy = %lf\n", dx, dy);
}

double getGridValueAtSamplingCoordinates(Map* map, int x, int y){
    return map->grid[map->Y * y + x];
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
    // printf("k: %lld, l: %lld \n", k, l);

    double x_k = k * map->dx;
    double x_k1 = (k+1) * map->dx;
    double y_l = l * map->dy;
    double y_l1 = (l+1) * map->dy;

    double prod1 = (x_k1 - x) * (y_l1 - y);
    double prod2 = (x_k1 - x) * (y - y_l);
    double prod3 = (x - x_k) * (y_l1 - y);
    double prod4 = (x - x_k) * (y - y_l);

    return (prod1 * getGridValueAtSamplingCoordinates(map, k, l)
            + prod2 * getGridValueAtSamplingCoordinates(map, k, l+1)
            + prod3 * getGridValueAtSamplingCoordinates(map, k+1, l)
            + prod4 * getGridValueAtSamplingCoordinates(map, k+1, l+1))/(map->dx*map->dy);
}

double getGridValueAtDomainCoordinates(Map* map, double x, double y){
    fprintf(stderr, "x = %lf, y = %lf\n", x, y);
    double epsilon = 10e-6;
    // Sampling step
    double dx = map->a/map->X;
    double dy = map->b/map->Y;

    // If value already in the grid, use that instead of interpolating
    if(fmod(x, dx) < epsilon && fmod(y, dy)  < epsilon){
        return getGridValueAtSamplingCoordinates(map, trunc(x/dx), trunc(y/dy)); 
    } else {
        return bilinearInterpolation(map, x, y);
    }
}

void printGrid(Map* map){
    for(int i = 0; i < map->X * map->Y; i++){
        printf("%lf ", map->grid[i]);
        if(i != 0 && i % map-> X == 0){
            printf("\n");
        }
    }
}

Map* readMapFile(const char* filename) {
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
    map->dx = map->a/map->X;
    map->dy = map->b/map->Y;

    map->grid = malloc((map->X * map->Y) * sizeof(double));

    if (map->grid == NULL) {
        fprintf(stderr, "Unable to allocate memory for the grid\n");
        free(map);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    // Read the barymetric depth grid from the map file

    long long i = 0;
    while (fread(buffer, 8, 1, fp) == 1) {
    // While we still read 8 contiguous bytes, fill the array
        map->grid[i] = *((double*)buffer);
        i++;
    }

    fclose(fp);
    return map;
}

double** allocateDoubleMatrix(int x, int y){
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
    for(int i = 0; i < x; i++)
        free(matrix[i]);

    free(matrix);
}

void printDoubleMatrix(double** matrix, int x, int y){
    for(int j = 0; j < y; j++){
        for(int i = 0; i < x; i++){
            fprintf(stderr, "%lf ", matrix[i][j]);
        }
        fprintf(stderr, "\n");
    }
}

int eulerExplicit(Map* map, Parameters* params, double*** nu, double*** u, double*** v){
    int xSize = (int)(map->a / params->deltaX);
    int ySize = (int)(map->b / params->deltaY);
    fprintf(stderr, "xSize = %d ySize = %d \n", xSize, ySize);

    // Allocate memory
    // nu in {0, 1, ..., a/dx}X{0, 1, ..., b/dy}
    double** nuCurr = allocateDoubleMatrix(xSize + 1, ySize + 1);
    if(!nuCurr){
        return -1;
    }

    double** nuNext = allocateDoubleMatrix(xSize + 1, ySize + 1);
    if(!nuNext){
        freeDoubleMatrix(nuCurr, xSize + 1);
        return -1;
    }

    // u in {-1/2, 1/2, ..., a/dx + 1/2}X{0, 1, ..., b/dy}
    double** uCurr = allocateDoubleMatrix(xSize + 2, ySize + 1);
    if(!uCurr){
        freeDoubleMatrix(nuCurr, xSize + 1);
        freeDoubleMatrix(nuNext, xSize + 1);
        return -1;
    }

    double** uNext = allocateDoubleMatrix(xSize + 2, ySize + 1);
    if(!uNext){
        freeDoubleMatrix(nuCurr, xSize + 1);
        freeDoubleMatrix(nuNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        return -1;
    }

    // v in {0, 1, .., a/dx}X{-1/2, 1/2, ..., b/dy + 1/2}
    double** vCurr = allocateDoubleMatrix(xSize + 1, ySize + 2);
    if(!vCurr){
        freeDoubleMatrix(nuCurr, xSize + 1);
        freeDoubleMatrix(nuNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        freeDoubleMatrix(uNext, xSize + 2);
        return -1;
    }

    double** vNext = allocateDoubleMatrix(xSize + 1, ySize + 2);
    if(!vNext){
        freeDoubleMatrix(nuCurr, xSize + 1);
        freeDoubleMatrix(nuNext, xSize + 1);
        freeDoubleMatrix(uCurr, xSize + 2);
        freeDoubleMatrix(uNext, xSize + 2);
        freeDoubleMatrix(vCurr, xSize + 1);
        return -1;
    }

    // h in {-1/2, 0, 1/2, ..., a/dx, a/dx + 1/2}X{-1/2, 0, 1/2, ..., b/dy, b/dy + 1/2}
    double** h = allocateDoubleMatrix(2 * xSize + 3, 2 * ySize + 3);
    if(!h){
        freeDoubleMatrix(nuCurr, xSize + 1);
        freeDoubleMatrix(nuNext, xSize + 1);
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

    printf("h\n");
    printDoubleMatrix(h, 2*xSize+3, 2*ySize+3);
    printf("apres h\n");

    for(int i = 0; i < xSize + 1; i++){
        for(int j = 0; j < ySize + 1; j++)
            nuCurr[i][j] = 0;
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
        printf("t = %u\n", t);

        // Compute nuNext
        // Separate for loop for cache optimization
        /*
        for(int i = 0; i < xSize + 1; i++)
            nuNext[i][0] = 0;

        for(int i = 0; i < xSize + 1; i++)
            nuNext[i][ySize] = 0;

        for(int i = 0; i < ySize + 1; i++)
            nuNext[0][i] = 0;

        for(int i = 0; i < ySize + 1; i++)
            nuNext[xSize][i] = 0;
        */

        for(int i = 0; i < xSize + 1; i++){
            for(int j = 0; j < ySize + 1; j++){
                nuNext[i][j] = (-(h[2*i+2][2*j+1] * uCurr[i+1][j] - h[2*i][2*j+1] * uCurr[i][j]) / params->deltaX 
                                -(h[2*i+1][2*j+2] * vCurr[i][j+1] - h[2*i+1][2*j] * vCurr[i][j]) / params->deltaY)
                                * params->deltaT + nuCurr[i][j];
            }
        }

        // Compute uNext
        for(int i = 0; i < ySize + 1; i++)
            uNext[0][i] = 0;

        for(int i = 0; i < ySize + 1; i++)
            uNext[xSize+1][i] = 0;

        for(int i = 1; i < xSize + 1; i++){
            for(int j = 0; j < ySize + 1; j++){
                uNext[i][j] = (-params->g * (nuCurr[i][j] - nuCurr[i-1][j]) / params->deltaX
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
                vNext[i][j] = (-params->g * (nuCurr[i][j] - nuCurr[i][j-1]) / params->deltaY
                               -params->gamma * vCurr[i][j]) * params->deltaT + vCurr[i][j];
            }
        }

        printf("nuCurr\n");
        printDoubleMatrix(nuCurr, xSize + 1, ySize + 1);
        printf("nuNext\n");
        printDoubleMatrix(nuNext, xSize + 1, ySize + 1);
        printf("uCurr\n");
        printDoubleMatrix(uCurr, xSize + 2, ySize + 1);
        printf("uNext\n");
        printDoubleMatrix(uNext, xSize + 2, ySize + 1);
        printf("vCurr\n");
        printDoubleMatrix(vCurr, xSize + 1, ySize + 2);
        printf("vNext\n");
        printDoubleMatrix(vNext, xSize + 1, ySize + 2);

        // Go to next step
        double** tmp;
        
        tmp = nuCurr;
        nuCurr = nuNext;
        nuNext = tmp;

        tmp = uCurr;
        uCurr = uNext;
        uNext = tmp;

        tmp = vCurr;
        vCurr = vNext;
        vNext = tmp;

    }

    *nu = nuCurr;
    *u = uCurr;
    *v = vCurr;
    
    freeDoubleMatrix(nuNext, xSize + 1);
    freeDoubleMatrix(uNext, xSize + 2);
    freeDoubleMatrix(vNext, xSize + 1);
    freeDoubleMatrix(h, 2 * xSize + 3);

    return 0;
}

int main(int argc, char const* argv[]) {
    // Check number of arguments
    assert(argc == 4);
    const char* parameter_file = argv[1];
    const char* map_file = argv[2];
    const unsigned int scheme = atoi(argv[3]);

    //Check argument validity
    assert((scheme == 0) || (scheme == 1));

    writeTestMap("test_map.dat");

    Parameters* param = readParameterFile(parameter_file);
    // Map* map = readMapFile("serverFiles/sriLanka.dat");
    Map* map = readMapFile("test_map.dat");

    // printGrid(map);
    // Explicit
    if (scheme == 0) {
        printf("Explicit ");
        printf("%s %s %u \n", parameter_file, map_file, scheme);
        
        double** nu;
        double** u;
        double** v;

        if(eulerExplicit(map, param, &nu, &u, &v) == -1){
            fprintf(stderr, "error in euler function\n");
            free(param);
            free(map->grid);
            free(map);
        }

        int xSize = (int)(map->a / param->deltaX);
        int ySize = (int)(map->b / param->deltaY);

        printf("nu\n");
        printDoubleMatrix(nu, xSize + 1, ySize + 1);
        printf("u\n");
        printDoubleMatrix(u, xSize + 2, xSize + 1);
        printf("v\n");
        printDoubleMatrix(v, xSize + 1, ySize + 2);

        freeDoubleMatrix(nu, xSize + 1);
        freeDoubleMatrix(u, xSize + 2);
        freeDoubleMatrix(v, xSize + 1);

    }
    // Implicit
    else {
        printf("Implicit ");
        printf("%s %s %u", parameter_file, map_file, scheme);
    }

    free(param);
    free(map->grid);
    free(map);
    /* code */
    return 0;
}
