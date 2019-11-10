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
    long long X;
    long long Y;
    double dx;
    double dy;
    double* grid;
} Map;

void writeTestMap(char* filename);
Parameters* readParameterFile(const char* filename);
void printUsefulMapInformation(Map* map);
double getGridValueAtSamplingCoordinates(Map* map, long long x, long long y);
double bilinearInterpolation(Map* map, double x, double y);
double getGridValueAtDomainCoordinates(Map* map, double x, double y);
void printGrid(Map* map);
Map* readMapFile(const char* filename);



#endif // MAIN_H_