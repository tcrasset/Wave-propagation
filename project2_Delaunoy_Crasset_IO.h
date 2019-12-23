#ifndef IO_H_
#define IO_H_

#define MAX_FILENAME_SIZE 500

/**
 * Structure representing the parameters with which
 * to make the computations
 */
typedef struct Parameters {
    const char* filename;
    double g;
    double gamma;
    double deltaX;
    double deltaY;
    double deltaT;
    double TMax;
    double A;
    double f;
    unsigned int S;
    unsigned int s;
    double r_threshold;
} Parameters;

/**
 * Structure representing a map
 */
typedef struct Map {
    double a;
    double b;
    int X;
    int Y;
    double dx;
    double dy;
    double** grid;
} Map;

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
double bilinearInterpolation(Map* map, double x, double y);

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
double getGridValueAtDomainCoordinates(Map* map, double x, double y);

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
double** allocateDoubleMatrix(int x, int y);

/**
 * Free a double matrix
 *
 * Parameters:
 * matrix: The matrix to free
 * x: The first index size
 */
void freeDoubleMatrix(double** matrix, int x);

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
double* transformMatrixToArray(double** matrix, int x, int y);

/**
 * Read a map file and load it in a Map structure
 *
 * Parameters:
 * filename: The name of the file to load
 * 
 * Returns:
 * A pointer to a Map structure containing the info in the file
 */
Map* readMapFile(const char* filename);

/**
 * Read a parameter file and load it in a Parameters structure
 *
 * Parameters:
 * filename: The name of the file to load
 * 
 * Returns:
 * A pointer to a Parameters structure containing the info in the file
 */
Parameters* readParameterFile(const char* filename);

/**
 * Write an array into a file
 *
 * Parameters:
 * filename: The name of the file in which to write the array
 * xsize: The size of the discretization along the x axis
 * ysize: The size of the discretization along the y axis
 * array: The array to write in the file
 */
void writeResultArray(char* filename, int xsize, int ysize,double* array);

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
void getFileNames(char* etaName, char* uName, char* vName, char* dir_name, unsigned int iteration);

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
    unsigned int ySize, unsigned int iteration, Parameters* params, int nbproc, int nbthreads);
    
#endif  // IO_H_