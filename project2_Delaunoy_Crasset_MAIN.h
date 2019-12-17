#ifndef MAIN_H_
#define MAIN_H_

#define MAX_FILENAME_SIZE 500

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

typedef struct Map {
    double a;
    double b;
    int X;
    int Y;
    double dx;
    double dy;
    double** grid;
} Map;

typedef struct SparseMatrix {
    unsigned int nonzeroSize;
    unsigned int* m;
    unsigned int* n;
    unsigned int* v;
} SparseMatrix;

SparseMatrix* toSparseMatrix(double** matrix, int xSize, int ySize);
void printUsefulMapInformation(Map* map);
double bilinearInterpolation(Map* map, double x, double y);
double getGridValueAtDomainCoordinates(Map* map, double x, double y);
double** allocateDoubleMatrix(int x, int y);
void freeDoubleMatrix(double** matrix, int x, int debug);
double* transformMatrixToArray(double** matrix, int x, int y);
void get_array_sizes(int rank, int nbproc, int xSize, int* size_X, int* size_X_u, int* size_X_h, int* startval_X_h, int* endval_X_h);
void gather_and_save(double** eta, double**  u, double**  v, int xSize, int ySize,  int debug, unsigned int iteration, Parameters* params);
int eulerExplicitMPI(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug, int debug_rank);
#endif  // MAIN_H_