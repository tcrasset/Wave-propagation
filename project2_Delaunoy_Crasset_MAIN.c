#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>

#include "project2_Delaunoy_Crasset_IO.h"
#include "project2_Delaunoy_Crasset_IMPLICIT.h"
#include "project2_Delaunoy_Crasset_EXPLICIT.h"

#define M_PI 3.14159265358979323846

int main(int argc, char* argv[]) {

    // Check number of arguments
    assert(argc >= 4);
    const char* parameter_file = argv[1];
    const char* map_file = argv[2];
    const unsigned int scheme = atoi(argv[3]);

    //Check argument validity
    assert((scheme == 0) || (scheme == 1));

    // Init MPI
    int nbproc, myrank;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if(provided != MPI_THREAD_FUNNELED){
        fprintf(stderr, "wrong provided = %d", provided);
    }   
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Measure execution time
    double startTime = MPI_Wtime();

    // Read parameters and map
    Parameters* params = readParameterFile(parameter_file);
    Map* map = readMapFile(map_file);

    // Explicit
    if (scheme == 0) {
        double** eta;
        double** u;
        double** v;

        int status = eulerExplicitMPI(map, params, &eta, &u, &v);
        if(status == -1){
            fprintf(stderr, "Error in euler function\n");
            free(params);
            free(map->grid);
            free(map);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }

    // Implicit
    else {
        double* eta;
        double* u;
        double* v;
        if(eulerImplicitMPI(map, params, &eta, &u, &v) == -1){
            fprintf(stderr, "error in euler function\n");
            free(params);
            free(map->grid);
            free(map);
            MPI_Finalize();
            exit(EXIT_FAILURE);
            
        }
    }

    // Mesure execution time
    double endTime = MPI_Wtime();
    double executionTime = endTime - startTime;
    char* openMP_nbthreads = getenv("OMP_NUM_THREADS");

    // Print statistics to standard output (for later analysis)
    fprintf(stdout,"%d,%d,%d,%s,%lf,%lf,%lf,%lf,%u,%lf\n", scheme, myrank, nbproc, \
                openMP_nbthreads, executionTime, params->deltaX, params->deltaY, \
                params->deltaT, params->s, params->r_threshold);
    
    MPI_Finalize();

    return 0;
}
