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
#include "project2_Delaunoy_Crasset_EXPLICIT.h"

/**
 * Compute the value of the map at a given point using bilinear interpolation
 *
 * Parameters:
 * map: The map used
 * x: The x coordinate at which to evaluate the map
 * y: The y coordinate at which to evaluate the map
 * 
 * Returns:
 * The evaluation of map a point x, y
 */
double bilinearInterpolation(Map* map, double x, double y){

    assert(0 <= x);
    assert(0 <= y);
    assert(x <= map->a);
    assert(y <= map->b);

    // Sampling coordinates
    int k = trunc(x/map->dx);
    int l = trunc(y/map->dy);

    // Compute border values
    double x_k = k * map->dx;
    double x_k1 = (k+1) * map->dx;
    double y_l = l * map->dy;
    double y_l1 = (l+1) * map->dy;

    // Compute coeficients
    double prod1 = (x_k1 - x) * (y_l1 - y);
    double prod2 = (x_k1 - x) * (y - y_l);
    double prod3 = (x - x_k) * (y_l1 - y);
    double prod4 = (x - x_k) * (y - y_l);

    // Interpolate
    double return_value = prod1 * map->grid[k][l];  
    double epsilon = 10e-6;

    // Robust implementation of the statement
    if(fabs(x - map->a) > epsilon)
        return_value += prod3 * map->grid[k+1][l];

    if(fabs(y - map->b) > epsilon)
        return_value += prod2 * map->grid[k][l+1];

    if(fabs(x - map->a) > epsilon && fabs(y - map->b) > epsilon)
        return_value += prod4 * map->grid[k+1][l+1];

    return_value /= map->dx*map->dy;

    return return_value;

}

/**
 * Compute the value of the map at a given point
 *
 * Parameters:
 * map: The map used
 * x: The x coordinate at which to evaluate the map
 * y: The y coordinate at which to evaluate the map
 * 
 * Returns:
 * The evaluation of map a point x, y
 */
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

/**
 * Allocate a double matrix
 *
 * Parameters:
 * x: The first index size
 * y: The second index size
 * 
 * Returns:
 * The allocated matrix
 */
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

/**
 * Free a double matrix
 *
 * Parameters:
 * matrix: The matrix to free
 * x: The first index size
 */
void freeDoubleMatrix(double** matrix, int x){
    assert(x > 0);
    assert(matrix);

    for(int i = 0; i < x; i++){
        if(matrix[i] != NULL){
            free(matrix[i]);
        }
    }
    free(matrix);
}

/**
 * Tansform a matrix to a 1d flattenned array
 *
 * Parameters:
 * matrix: The matrix to flatten
 * x: The first index size
 * y: The second index size
 * 
 * Returns:
 * The flattenned array
 */
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

/**
 * Read a map file and load it in a Map structure
 *
 * Parameters:
 * filename: The name of the file to load
 * 
 * Returns:
 * A pointer to a Map structure containing the info in the file
 */
Map* readMapFile(const char* filename) {
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
            map->grid[col][row] = *((double*)buffer);
        }
    }

    fclose(fp);
    return map;
}

/**
 * Read a parameter file and load it in a Parameters structure
 *
 * Parameters:
 * filename: The name of the file to load
 * 
 * Returns:
 * A pointer to a Parameters structure containing the info in the file
 */
Parameters* readParameterFile(const char* filename) {
    FILE* fp;

    // Open file
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Allocate structure
    Parameters* params = malloc(sizeof(Parameters));

    if (params == NULL) {
        fprintf(stderr, "Unable to allocate memory for parameters\n");
        fclose(fp);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Fille structure with info from the file
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

/**
 * Write an array into a file
 *
 * Parameters:
 * filename: The name of the file in which to write the array
 * xsize: The size of the discretization along the x axis
 * ysize: The size of the discretization along the y axis
 * array: The array to write in the file
 */
void writeResultArray(char* filename, int xsize, int ysize, double* array) {
    
    // Open file
    FILE* fp = fopen(filename, "wb");

    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Write sizes
    fwrite(&xsize, sizeof(xsize), 1, fp);
    fwrite(&ysize, sizeof(ysize), 1, fp);

    // Write array
    for (int row = ysize - 1; row >= 0; row--) {
        for (int col = 0; col < xsize; col++) {
            int index = col * ysize + row;
            fwrite(&array[index], 8, 1, fp);
        }
    }

    fclose(fp);
}

/**
 * Build file names for each variable to save
 *
 * Parameters:
 * etaName: A pointer to char that will be set to the name of eta
 * uName: A pointer to char that will be set to the name of u
 * vName: A pointer to char that will be set to the name of v
 * dir_name: A string containg the name of the directory to include in the filename
 * iteration: The iteration at which to make the save
 */
void getFileNames(char* etaName, char* uName, char* vName, char* dir_name, unsigned int iteration) {
    
    // Base prefixes
    char* etaPrefix = "eta";
    char* uPrefix = "u";
    char* vPrefix = "v";

    //Build suffix based on ieration
    char file_suffix[MAX_FILENAME_SIZE];
    snprintf(file_suffix, MAX_FILENAME_SIZE, "_%u", iteration);

    // Build eta name
    strncpy(etaName, dir_name, MAX_FILENAME_SIZE);
    strncat(etaName, etaPrefix, MAX_FILENAME_SIZE);
    strncat(etaName, file_suffix, MAX_FILENAME_SIZE);
    strncat(etaName, ".dat", MAX_FILENAME_SIZE);

    // Build u name
    strncpy(uName, dir_name, MAX_FILENAME_SIZE);
    strncat(uName, uPrefix, MAX_FILENAME_SIZE);
    strncat(uName, file_suffix, MAX_FILENAME_SIZE);
    strncat(uName, ".dat", MAX_FILENAME_SIZE);

    // Build v name
    strncpy(vName, dir_name, MAX_FILENAME_SIZE);
    strncat(vName, vPrefix, MAX_FILENAME_SIZE);
    strncat(vName, file_suffix, MAX_FILENAME_SIZE);
    strncat(vName, ".dat", MAX_FILENAME_SIZE);
}

/**
 * Save system state to disk
 *
 * Parameters:
 * etaTotal: An array containing the values of eta
 * uTotal: An array containing the values of u
 * vTotal: An array containing the values of v
 * xSize: The size of the discretization along the x axis
 * ySize: The size of the discretization along the y axis
 * iteration: The iteration at which the save is performed
 * params: The parameters of the run
 * nbproc: The number of processes
 * nbthreads: The number of threads
 * 
 * Returns:
 * A pointer to a Map structure containing the info in the file
 */
void saveToDisk(double* etaTotal, double* uTotal, double* vTotal, unsigned int xSize, 
    unsigned int ySize, unsigned int iteration, Parameters* params, int nbproc, int nbthreads) {

    static int createDirectory = 0;
    static char full_path[MAX_FILENAME_SIZE];
    int status = 0;

    // Attempt creating a directory when calling this function for the first time
    if(createDirectory == 0){
        createDirectory = 1;

        //Get current working directory
        char current_dir[MAX_FILENAME_SIZE];
        getcwd(current_dir, MAX_FILENAME_SIZE);

        // Get parameter file and remove '.txt' extension
        char* parameter_file = basename((char*)params->filename);
        parameter_file[strlen(parameter_file)-4] = 0; 
        
        // Create output directory
        char new_dir[MAX_FILENAME_SIZE];
        snprintf(new_dir, MAX_FILENAME_SIZE, "/Results/matrices_of_%s_%d_%d/", parameter_file, nbproc, nbthreads);
        strncpy(full_path, current_dir, MAX_FILENAME_SIZE);
        strncat(full_path, new_dir, MAX_FILENAME_SIZE);

        // Check if file exists, if not, create it.
        if(access(full_path, F_OK) == -1) { 
            status = mkdir(full_path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if(status == -1){
                fprintf(stderr, "Error in saveToDisk = %s\n" , strerror(errno));
                fprintf(stderr, "Attempted to create: %s\n" , full_path);
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        }
    }

    // Create file names
    char etaFilename[MAX_FILENAME_SIZE];
    char uFilename[MAX_FILENAME_SIZE];
    char vFilename[MAX_FILENAME_SIZE];
    getFileNames(etaFilename, uFilename, vFilename, full_path, iteration);

    writeResultArray(etaFilename, xSize + 1, ySize + 1, etaTotal);
    writeResultArray(uFilename, xSize + 2, ySize + 1, uTotal);
    writeResultArray(vFilename, xSize + 1, ySize + 2, vTotal);

}