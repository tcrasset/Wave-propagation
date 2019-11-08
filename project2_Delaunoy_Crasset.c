#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>


void readParameterFile(const char* filename, double* g, double* gamma, 
                        double* deltaX, double* deltaY, double* deltaT, 
                        double* A, double* func, unsigned int* S, 
                        unsigned int* s, double* r_threshold){
    FILE * fp;

    fp = fopen(filename, "r");
    if (fp == NULL){
        exit(EXIT_FAILURE);
    }

    fscanf(fp,"%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%u\n%u\n%lf", g, gamma, deltaX, deltaY, deltaT, A, func, S, s, r_threshold);
    fclose(fp);

    return;
}

double* readMapFile(const char* filename, double* a, double* b, 
                        long long* X,  long long* Y){
    FILE * fp;

    fp = fopen(filename, "rb");
    if (fp == NULL){
        exit(EXIT_FAILURE);
    }

    // Read the parameters at the top of the map file
    // Assumption : there are no spaces and no end of lines, just
    // contiguous bytes in double precision (8 bytes per unit)
    char buffer[8];
    fread(buffer, 8, 1, fp);
    a = ((double*)buffer);

    fread(buffer, 8, 1, fp);
    b = ((double*)buffer);

    fread(buffer, 8, 1, fp);
    X = ((long long*)buffer);

    fread(buffer, 8, 1, fp);
    Y = ((long long*)buffer);
    
    double* grid = malloc((*X) * (*Y) * sizeof(double));
    
    if (grid == NULL){
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    // While we still read 8 contiguous bytes,
    // fill the array
    long long i = 0;
    while(fread(buffer, 8, 1, fp) == 8) {
        grid[i] = *((double*)buffer);
    }

    fclose(fp);
    return grid;
}

int main(int argc, char const* argv[]) {

    // Check number of arguments
    assert(argc == 4);  
    const char* parameter_file = argv[1];
    const char* map_file = argv[2];
    const unsigned int scheme = atoi(argv[3]);
    
    //Check argument validity
    assert((scheme == 0) || (scheme == 1));

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

    //readParameterFile(parameter_file, &g, &gamma, &deltaX, &deltaY, &deltaT, &A, &func, &S, &s, &r_threshold);

    double a; 
    double b; 
    long long X;  
    long long Y;
    double* grid = readMapFile(map_file, &a, &b, &X, &Y);

    // Explicit
    if (scheme == 0) {
        printf("Explicit");
        printf("%s %s %u", parameter_file, map_file, scheme);

        return 0;
    }
    // Implicit
    else {
        printf("Implicit");
        printf("%s %s %u", parameter_file, map_file, scheme);
        return 0;
    }
    /* code */
    return 0;
}
