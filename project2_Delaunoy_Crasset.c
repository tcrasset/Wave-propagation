#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>

typedef struct Parameters {
    double g;
    double gamma;
    double deltaX;
    double deltaY;
    double deltaT;
    double A;
    double func;
    unsigned int S;
    unsigned int s;
    double r_threshold;
} Parameters;

typedef struct Map {
    double a;
    double b;
    long long X;
    long long Y;
    double* grid;
} Map;


void writeTestMap(char* filename){

    FILE* fp;

    fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    double a = 5;
    double b = 10;
    long long X = 10;
    long long Y = 20;

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
        free(fp);
        exit(EXIT_FAILURE);
    }

    fscanf(fp, "%lf", &params->g);
    fscanf(fp, "%lf", &params->gamma);
    fscanf(fp, "%lf", &params->deltaX);
    fscanf(fp, "%lf", &params->deltaY);
    fscanf(fp, "%lf", &params->deltaT);
    fscanf(fp, "%lf", &params->A);
    fscanf(fp, "%lf", &params->func);
    fscanf(fp, "%u", &params->S);
    fscanf(fp, "%u", &params->s);
    fscanf(fp, "%lf", &params->r_threshold);

    fclose(fp);

    return params;
}


long long getIndex(Map* map, long long x, long long y){
    return map->X * x +  y;
}

void printCoordinate(Map* map, long long x, long long y){
    
    long long index = getIndex(map, x, y);
    printf("map(%lld,%lld) = %lf \t map[%lld] \n", x, y, map->grid[index], index);

}

double getGridValueAtSampling(Map* map, long long x, long long y){
    return map->grid[map->X * x +  y];
}

double bilinearInterpolation(Map* map, double x, double y){

    double dx = map->a/map->X;
    double dy = map->b/map->Y;

    long long k = trunc(x/dx);
    long long l = trunc(y/dy);

    double x_k = k * dx;
    double x_k1 = (k+1) * dx;
    double y_l = l * dy;
    double y_l1 = (l+1) * dy;

    double prod1 = (x_k1 - x) * (y_l1 - y);
    double prod2 = (x_k1 - x) * (y - y_l);
    double prod3 = (x - x_k) * (y_l1 - y);
    double prod4 = (x - x_k) * (y - y_l);

    return (prod1 * getGridValueAtSampling(map, k, l)
            + prod2 * getGridValueAtSampling(map, k, l+1)
            + prod3 * getGridValueAtSampling(map, k+1, l)
            + prod4 * getGridValueAtSampling(map, k+1, l+1))/(dx*dy);
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
        free(fp);
        exit(EXIT_FAILURE);
    }

    char buffer[8];
    fread(buffer, 8, 1, fp);
    map->a = *((double*)buffer);

    fread(buffer, 8, 1, fp);
    map->b = *((double*)buffer);

    fread(buffer, 8, 1, fp);
    map->X = *((long long*)buffer);

    fread(buffer, 8, 1, fp);
    map->Y = *((long long*)buffer);

    // //Override the parameters (temporary)
    // map->X = 8000;
    // map->Y = 8000;

    map->grid = malloc((map->X * map->Y) * sizeof(double));

    if (map->grid == NULL) {
        fprintf(stderr, "Unable to allocate memory for the grid\n");
        free(map);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    // While we still read 8 contiguous bytes,
    // fill the array
    long long i = 0;
    while (fread(buffer, 8, 1, fp) == 1) {
        map->grid[i] = *((double*)buffer);
        i++;
    }

    fclose(fp);
    return map;
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
    Map* map = readMapFile("test_map.dat");

    printCoordinate(map ,0, 0);
    printf("Bilinear interp : %lf\n", bilinearInterpolation(map, 2,0.2));
    printGrid(map);

    free(param);
    free(map->grid);
    free(map);

    // Explicit
    if (scheme == 0) {
        printf("Explicit ");
        printf("%s %s %u", parameter_file, map_file, scheme);

        return 0;
    }
    // Implicit
    else {
        printf("Implicit ");
        printf("%s %s %u", parameter_file, map_file, scheme);
        return 0;
    }
    /* code */
    return 0;
}
