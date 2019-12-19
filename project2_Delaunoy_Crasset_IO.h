#ifndef IO_H_
#define IO_H_

#define MAX_FILE_SIZE 500

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

void printDoubleMatrix(double** matrix, int x, int y, int process_rank);
void printLinearArray(double* array, int x, int y);
Map* readMapFile(const char* filename, int debug);
Parameters* readParameterFile(const char* filename);
double bilinearInterpolation(Map* map, double x, double y);
double getGridValueAtDomainCoordinates(Map* map, double x, double y);
double** allocateDoubleMatrix(int x, int y);
void freeDoubleMatrix(double** matrix, int x, int debug);
void printUsefulMapInformation(Map* map);
void writeResultMatrix(char* filename, int xsize, int ysize,double** matrix, int debug);
void writeResultArray(char* filename, int xsize, int ysize,double* array, int debug);
void writeTestMap(char* filename, int debug);
void printGrid(Map* map);
void getFileNames(char* etaName, char* uName, char* vName, char* dir_name, unsigned int iteration);
int saveToDisk(double* etaTotal, double* uTotal, double* vTotal, unsigned int xSize, unsigned int ySize, unsigned int iteration, Parameters* params);
#endif  // IO_H_