#include <assert.h>
#include <libgen.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <mpi.h>

#include "project2_Delaunoy_Crasset_IO.h"

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
    double epsilon = 10e-6;
    assert(x >= 0);
    assert(y >= 0);
    // Sampling step

    // If value already in the grid, use that instead of interpolating
    if(fmod(x, map->dx) < epsilon && fmod(y, map->dy)  < epsilon){
        return map->grid[(int) trunc(x/map->dx)][(int) trunc(y/map->dy)]; 
    } else {
        return bilinearInterpolation(map, x, y);
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

void printDoubleMatrix(double** matrix, int x, int y, int process_rank) {
    assert(x > 0);
    assert(y > 0);
    assert(matrix != NULL);

    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            fprintf(stderr, "\tP%d\t", process_rank);
            fprintf(stderr, "%lf ", matrix[i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void printLinearArray(double* array, int x, int y) {
    for (int i = 0; i < x * y; i++) {
        fprintf(stderr, "%lf ", array[i]);
        if (i != 0 && (i + 1) % (y) == 0) {
            fprintf(stderr, "\n");
        }
    }
    fprintf(stderr, "\n");
}

Map* readMapFile(const char* filename, int debug) {
    FILE* fp;
    char buffer[8];

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Read the parameters at the top of the map file
    // Assumption : there are no spaces and no end of lines, just
    // contiguous bytes in double precision (8 bytes per unit)

    Map* map = malloc(sizeof(Map));

    if (map == NULL) {
        fprintf(stderr, "Unable to allocate memory for the map\n");
        fclose(fp);
        MPI_Finalize();
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
    map->dx = map->a / (map->X - 1);
    map->dy = map->b / (map->Y - 1);

    map->grid = allocateDoubleMatrix(map->X, map->Y);

    if (map->grid == NULL) {
        fprintf(stderr, "Unable to allocate memory for the grid\n");
        free(map);
        fclose(fp);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Read the bathymetry depth grid from the map file

    long long i = 0;
    fread(buffer, 8, 1, fp);
    for (int row = map->Y - 1; row >= 0; row--) {
        for (int col = 0; col < map->X; fread(buffer, 8, 1, fp), col++) {
            if (debug == 1)
                printf("(%d,%d)=%lf\n", col, row, *((double*)buffer));
            map->grid[col][row] = *((double*)buffer);
        }
    }

    fclose(fp);
    return map;
}

Parameters* readParameterFile(const char* filename) {
    FILE* fp;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    Parameters* params = malloc(sizeof(Parameters));

    if (params == NULL) {
        fprintf(stderr, "Unable to allocate memory for parameters\n");
        fclose(fp);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    params->filename = filename;
    fscanf(fp, "%lf", &params->g);
    fscanf(fp, "%lf", &params->gamma);
    fscanf(fp, "%lf", &params->deltaX);
    fscanf(fp, "%lf", &params->deltaY);
    fscanf(fp, "%lf", &params->deltaT);
    fscanf(fp, "%lf", &params->TMax);
    fscanf(fp, "%lf", &params->A);
    fscanf(fp, "%lf", &params->f);
    fscanf(fp, "%u", &params->S);
    fscanf(fp, "%u", &params->s);
    fscanf(fp, "%lf", &params->r_threshold);

    fclose(fp);

    return params;
}

void printUsefulMapInformation(Map* map) {
    printf("X: %d Y: %d\n", map->X, map->Y);
    printf("a : %lf, b : %lf \n", map->a, map->b);
    printf("Sampling steps: dx = %lf, dy = %lf\n", map->dx, map->dy);
}

void writeResultMatrix(char* filename, int xsize, int ysize, double** matrix, int debug) {
    FILE* fp;

    fp = fopen(filename, "wb");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    fwrite(&xsize, sizeof(xsize), 1, fp);
    fwrite(&ysize, sizeof(ysize), 1, fp);

    for (int row = ysize - 1; row >= 0; row--) {
        for (int col = 0; col < xsize; col++) {
            if (debug == 1) {
                printf("%lf \n", matrix[col][row]);
            }
            fwrite(&matrix[col][row], 8, 1, fp);
        }
    }

    fclose(fp);
}

void writeResultArray(char* filename, int xsize, int ysize, double* array, int debug) {
    
    FILE* fp = fopen(filename, "wb");

    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    fwrite(&xsize, sizeof(xsize), 1, fp);
    fwrite(&ysize, sizeof(ysize), 1, fp);

    for (int row = ysize - 1; row >= 0; row--) {
        for (int col = 0; col < xsize; col++) {
            int index = col * ysize + row;
            fwrite(&array[index], 8, 1, fp);
        }
    }

    fclose(fp);
}

void writeTestMap(char* filename, int debug) {
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

    fwrite(&a, sizeof(a), 1, fp);
    fwrite(&b, sizeof(b), 1, fp);
    fwrite(&X, sizeof(X), 1, fp);
    fwrite(&Y, sizeof(Y), 1, fp);

    for (int row = 0; row < Y; row++) {
        for (int col = 0; col < X; col++) {
            double value = (double)(Y - row - 1) * X + col;
            if (debug == 1)
                printf("%lf \n", value);
            fwrite(&value, sizeof(value), 1, fp);
        }
    }

    fclose(fp);
}

void printGrid(Map* map) {
    assert(map);
    for (int i = 0; i < map->X; i++) {
        for (int j = 0; j < map->Y; j++) {
            printf("%lf ", map->grid[i][j]);
        }
        printf("\n");
    }
}

void getFileNames(char* etaName, char* uName, char* vName, char* dir_name, unsigned int iteration) {
    char* etaPrefix = "eta";
    char* uPrefix = "u";
    char* vPrefix = "v";

    char file_suffix[MAX_FILE_SIZE];
    snprintf(file_suffix, MAX_FILE_SIZE, "_%u", iteration);

    strncpy(etaName, dir_name, MAX_FILE_SIZE);
    strncat(etaName, etaPrefix, MAX_FILE_SIZE);
    strncat(etaName, file_suffix, MAX_FILE_SIZE);
    strncat(etaName, ".dat", MAX_FILE_SIZE);

    strncpy(uName, dir_name, MAX_FILE_SIZE);
    strncat(uName, uPrefix, MAX_FILE_SIZE);
    strncat(uName, file_suffix, MAX_FILE_SIZE);
    strncat(uName, ".dat", MAX_FILE_SIZE);

    strncpy(vName, dir_name, MAX_FILE_SIZE);
    strncat(vName, vPrefix, MAX_FILE_SIZE);
    strncat(vName, file_suffix, MAX_FILE_SIZE);
    strncat(vName, ".dat", MAX_FILE_SIZE);
}

int saveToDisk(double* etaTotal, double* uTotal, double* vTotal, unsigned int xSize, unsigned int ySize, unsigned int iteration, Parameters* params) {
    static int createDirectory = 0;
    static char full_path[MAX_FILE_SIZE];
    int status = 0;

    // Attempt creating a directory when calling this function for the first time
    if (createDirectory == 0) {
        createDirectory = 1;

        //Get current working directory
        char current_dir[MAX_FILE_SIZE];
        getcwd(current_dir, MAX_FILE_SIZE);

        // Get parameter file and remove '.txt' extension
        char* parameter_file = basename((char*)params->filename);
        parameter_file[strlen(parameter_file)-4] = 0; 
        
        // Create output directory
        char new_dir[MAX_FILE_SIZE];
        snprintf(new_dir, MAX_FILE_SIZE, "/Results/matrices_of_%s/", parameter_file);
        strncpy(full_path, current_dir, MAX_FILE_SIZE);
        strncat(full_path, new_dir, MAX_FILE_SIZE);

        // Check if file exists, if not, create it.
        if(access(full_path, F_OK) == -1) { 
            status = mkdir(full_path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (status == -1) {
                fprintf(stderr, "Error in saveToDisk = %s\n" , strerror(errno));
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        }

    }
    
    // Create file names
    char etaFilename[MAX_FILE_SIZE];
    char uFilename[MAX_FILE_SIZE];
    char vFilename[MAX_FILE_SIZE];
    getFileNames(etaFilename, uFilename, vFilename, full_path, iteration);

    writeResultArray(etaFilename, xSize + 1, ySize + 1, etaTotal, 0);
    writeResultArray(uFilename, xSize + 2, ySize + 1, uTotal, 0);
    writeResultArray(vFilename, xSize + 1, ySize + 2, vTotal, 0);

}