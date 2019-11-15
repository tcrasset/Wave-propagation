#ifndef MAIN_H_
#define MAIN_H_

typedef struct Parameters {
    double g;
    double gamma;
    double deltaX;
    double deltaY;
    double deltaT;
    unsigned int TMax;
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

void writeTestMap(char* filename, int debug);
void writeResultMatrix(char* filename, int xsize, int ysize, double** matrix, int debug);

Parameters* readParameterFile(const char* filename);
void printUsefulMapInformation(Map* map);
double bilinearInterpolation(Map* map, double x, double y);
double getGridValueAtDomainCoordinates(Map* map, double x, double y);
void printGrid(Map* map);
Map* readMapFile(const char* filename, int debug);
int eulerExplicit(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug);

#endif  // MAIN_H_