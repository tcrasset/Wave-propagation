#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>

#include "project2_Delaunoy_Crasset_EXPLICIT.h"
#include "project2_Delaunoy_Crasset_IO.h"

#define M_PI 3.14159265358979323846

void get_array_sizes(int rank, int nbproc, int xSize, int* size_X, int* size_X_u, int* size_X_h, int* startval_X_h, int* endval_X_h){
    int mpi_xsize = xSize/nbproc;

    int startval_X, endval_X;
    int startval_X_u, endval_X_u;
    if(nbproc == 1){//Only 1 process
        startval_X = 0;
        endval_X = xSize;
        *startval_X_h = 0;
        *endval_X_h = 2*xSize + 2;
        startval_X_u = 0;
        endval_X_u = xSize+1;
    }
    else if(rank == 0){//Multiprocess
        startval_X = 0;
        endval_X = mpi_xsize;
        *startval_X_h = 0;
        *endval_X_h = 2*mpi_xsize + 2;
        startval_X_u = 0;
        endval_X_u = mpi_xsize;
    }else if(rank == nbproc -1){
        startval_X = rank * mpi_xsize + 1;
        endval_X = (rank+1) * mpi_xsize;
        *startval_X_h = 2 * rank * mpi_xsize + 2;
        *endval_X_h = 2 * (rank+1) * mpi_xsize + 2;
        startval_X_u = rank * mpi_xsize + 1;
        endval_X_u = (rank+1) * mpi_xsize + 1;
    }else{
        startval_X = rank * mpi_xsize + 1;
        endval_X = (rank+1) * mpi_xsize; 
        *startval_X_h = 2 * rank * mpi_xsize + 2;
        *endval_X_h = 2 * (rank+1) * mpi_xsize + 2;
        startval_X_u = rank * mpi_xsize + 1;
        endval_X_u = (rank+1) * mpi_xsize;
    }

    int remaining = xSize%nbproc;
    if(rank < remaining){
        startval_X += rank;
        endval_X += rank + 1;
        startval_X_u += rank;
        endval_X_u += rank + 1;
        *startval_X_h += rank * 2;
        *endval_X_h += (rank + 1) * 2;
    }
    else{
        *startval_X_h += remaining * 2;
        *endval_X_h += remaining * 2;
    }

    *size_X = endval_X - startval_X + 1;
    *size_X_u = endval_X_u - startval_X_u + 1;
    *size_X_h = *endval_X_h - *startval_X_h + 1;
}



void gather_and_save(double** eta, double**  u, double**  v, int xSize, int ySize,  int debug, unsigned int iteration, Parameters* params){

    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int size_X, size_X_u, size_X_h, startval_X_h, endval_X_h;
    get_array_sizes(myrank, nbproc, xSize, &size_X, &size_X_u, &size_X_h, &startval_X_h, &endval_X_h);
    
    int openMP_nbthreads = atoi(getenv("OMP_NUM_THREADS"));

    double* etaTotal;
    double* uTotal;
    double* vTotal;

    double* etaPartial = transformMatrixToArray(eta, size_X, ySize +1);
    double* uPartial = transformMatrixToArray(u, size_X_u, ySize +1);
    double* vPartial = transformMatrixToArray(v, size_X, ySize +2);

    if(nbproc != 1){

        int tmp_size_X;
        int tmp_size_X_u;
        int tmp_size_X_h;
        int tmp_startval_X_h;
        int tmp_endval_X_h;

        int* recvcounts_eta = malloc(nbproc * sizeof(int));
        int* recvcounts_u = malloc(nbproc * sizeof(int));
        int* recvcounts_v = malloc(nbproc * sizeof(int));
        int* disp_eta = malloc(nbproc * sizeof(int));
        int* disp_u = malloc(nbproc * sizeof(int));
        int* disp_v = malloc(nbproc * sizeof(int));

        if(!recvcounts_eta || !recvcounts_u || !recvcounts_v || !disp_eta || !disp_u || !disp_v){
            fprintf(stderr, "error malloc recvcounts\n");
            MPI_Finalize();
            exit(-1);
        }

        for(int i = 0; i < nbproc; i++){
            get_array_sizes(i, nbproc, xSize, &tmp_size_X, &tmp_size_X_u, &tmp_size_X_h, &tmp_startval_X_h, &tmp_endval_X_h);
            recvcounts_eta[i] = tmp_size_X * (ySize + 1);
            recvcounts_u[i] = tmp_size_X_u * (ySize + 1);
            recvcounts_v[i] = tmp_size_X * (ySize + 2);

            if(i == 0){
                disp_eta[0] = 0;
                disp_u[0] = 0;
                disp_v[0] = 0;
            }
            if (i < nbproc - 1){
                disp_eta[i + 1] = disp_eta[i] + tmp_size_X * (ySize + 1);
                disp_u[i + 1] = disp_u[i] + tmp_size_X_u * (ySize + 1);
                disp_v[i + 1] = disp_v[i] + tmp_size_X * (ySize + 2);
            }
        }

        etaTotal = malloc((xSize + 1) * (ySize  + 1)* sizeof(double));
        uTotal = malloc((xSize + 2) * (ySize  + 1)* sizeof(double));
        vTotal = malloc((xSize + 1) * (ySize  + 2)* sizeof(double));

        MPI_Gatherv(etaPartial, (size_X) * (ySize + 1) , MPI_DOUBLE, etaTotal, recvcounts_eta, disp_eta, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gatherv(uPartial, (size_X_u) * (ySize + 1) , MPI_DOUBLE, uTotal, recvcounts_u, disp_u, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gatherv(vPartial, (size_X) * (ySize + 2) , MPI_DOUBLE, vTotal, recvcounts_v, disp_v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        

        free(etaPartial);
        free(uPartial);
        free(vPartial);
        free(recvcounts_eta);
        free(recvcounts_u);
        free(recvcounts_v);
        free(disp_eta);
        free(disp_u);
        free(disp_v);

        if(myrank == 0){
            saveToDisk(etaTotal, uTotal, vTotal, xSize, ySize, iteration, params, nbproc, openMP_nbthreads);
        }
    }
    else{
        etaTotal = transformMatrixToArray(eta, xSize + 1, ySize +1);
        uTotal = transformMatrixToArray(u, xSize + 2, ySize +1);
        vTotal= transformMatrixToArray(v, xSize + 1, ySize +2);

        saveToDisk(etaTotal, uTotal, vTotal, xSize, ySize, iteration, params, nbproc, openMP_nbthreads);
    }

    free(etaTotal);
    free(uTotal);
    free(vTotal);
}

int eulerExplicitMPI(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug, int debug_rank){
    

    assert(map);
    assert(params);

    int nbproc, myrank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

    int xSize = (int)(map->a / params->deltaX);
    int ySize = (int)(map->b / params->deltaY);

    int size_X;
    int size_X_u;
    int size_X_h;
    int startval_X_h;
    int endval_X_h;
    get_array_sizes(myrank, nbproc, xSize, &size_X, &size_X_u, &size_X_h, &startval_X_h, &endval_X_h);

    // Allocate memory
    // eta in {0, 1, ..., a/dx}X{0, 1, ..., b/dy}
    double** etaCurr_ = allocateDoubleMatrix(size_X, ySize + 1);
    if(!etaCurr_){
        return -1;
    }

    double** etaNext = allocateDoubleMatrix(size_X, ySize + 1);
    if(!etaNext){
        freeDoubleMatrix(etaCurr_, size_X,0);
        return -1;
    }

    // u in {-1/2, 1/2, ..., a/dx + 1/2}X{0, 1, ..., b/dy}
    double** uCurr_ = allocateDoubleMatrix(size_X_u, ySize + 1);
    if(!uCurr_){
        freeDoubleMatrix(etaCurr_,size_X,0);
        freeDoubleMatrix(etaNext,size_X,0);
        return -1;
    }

    double** uNext = allocateDoubleMatrix(size_X_u, ySize + 1);
    if(!uNext){
        freeDoubleMatrix(etaCurr_,size_X,0);
        freeDoubleMatrix(etaNext,size_X,0);
        freeDoubleMatrix(uCurr_, size_X_u,0);
        return -1;
    }

    // v in {0, 1, .., a/dx}X{-1/2, 1/2, ..., b/dy + 1/2}
    double** vCurr_ = allocateDoubleMatrix(size_X, ySize + 2);
    if(!vCurr_){
        freeDoubleMatrix(etaCurr_, size_X,0);
        freeDoubleMatrix(etaNext, size_X,0);
        freeDoubleMatrix(uCurr_, size_X_u,0);
        freeDoubleMatrix(uNext, size_X_u,0);
        return -1;
    }

    double** vNext = allocateDoubleMatrix(size_X, ySize + 2);
    if(!vNext){
        freeDoubleMatrix(etaCurr_, size_X,0);
        freeDoubleMatrix(etaNext, size_X,0);
        freeDoubleMatrix(uCurr_, size_X_u,0);
        freeDoubleMatrix(uNext, size_X_u,0);
        freeDoubleMatrix(vCurr_, size_X,0);
        return -1;
    }
    
    // h in {-1/2, 0, 1/2, ..., a/dx, a/dx + 1/2}X{-1/2, 0, 1/2, ..., b/dy, b/dy + 1/2}
    double** h_ = allocateDoubleMatrix(size_X_h, 2 * ySize + 3);
    if(!h_){
        freeDoubleMatrix(etaCurr_, size_X,0);
        freeDoubleMatrix(etaNext, size_X,0);
        freeDoubleMatrix(uCurr_, size_X_u,0);
        freeDoubleMatrix(uNext, size_X_u,0);
        freeDoubleMatrix(vCurr_, size_X,0);
        freeDoubleMatrix(vNext, size_X,0);
        return -1;
    }

    // Compute h from the provided map file
    for(int i = startval_X_h; i <= endval_X_h; i++){
        for(int j = 0; j < 2 * ySize + 3; j++){
            h_[i-startval_X_h][j] = getGridValueAtDomainCoordinates(map, ((float)(i * xSize)/(xSize + 1)) * (params->deltaX / 2), ((float)(j * ySize)/(ySize + 1)) * (params->deltaY / 2));
        }
    }

    #pragma omp parallel default(shared)
    {   
        #pragma omp for schedule(static)
        for(int i = 0; i < size_X; i++){
            for(int j = 0; j < ySize; j++){
                etaCurr_[i][j] = 0;
            }
        }

        #pragma omp for schedule(static)
        for(int i = 0; i < size_X_u; i++){
            for(int j = 0; j < ySize; j++){
                uCurr_[i][j] = 0;
            }
        }

        #pragma omp for schedule(static)
        for(int i = 0; i < size_X; i++){
            for(int j = 0; j < ySize; j++)
                vCurr_[i][j] = 0;
        }
    }

    const double** etaCurr = (const double**) etaCurr_;
    const double** uCurr = (const double**) uCurr_;
    const double** vCurr = (const double**) vCurr_;
    const double** h = (const double**) h_;

    // Alocate arrays for receiving data from other process
    double* uReceived = malloc((ySize + 1) * sizeof(double));
    double* etaReceived = malloc((ySize + 1) * sizeof(double));

    // Starting time loop
    for(unsigned int t = 1; t <= params->TMax/params->deltaT; t++){
        
        if(myrank == 0){
            fprintf(stderr, "in loop t = %u\n", t);
        }

        // In a multiprocess environment, sending the leftmost column of u of the domain controlled
        // by the current process to the process with the previous rank 
        if(nbproc != 1){
            if(myrank == nbproc-1){
                MPI_Send(uCurr[0], ySize + 1, MPI_DOUBLE, myrank - 1, 62, MPI_COMM_WORLD); //Tag 62 is for u
            }else if (myrank == 0){
                MPI_Recv(uReceived, ySize + 1, MPI_DOUBLE, 1, 62, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }else{
                MPI_Sendrecv(uCurr[0], ySize + 1, MPI_DOUBLE, myrank - 1, 62,
                            uReceived, ySize + 1, MPI_DOUBLE, myrank + 1, 62,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // Compute the next value of eta
        #pragma omp parallel default(shared)
        {   
            // Process etaNext in one block
            if(myrank == nbproc-1 || nbproc == 1){

                #pragma omp for schedule(static)
                for(int i = 0; i < size_X; i++){
                    for(int j = 0; j < ySize + 1; j++){
                        etaNext[i][j] = (-(h[2*i+2][2*j+1] * uCurr[i+1][j] - h[2*i][2*j+1] * uCurr[i][j]) / params->deltaX 
                                        -(h[2*i+1][2*j+2] * vCurr[i][j+1] - h[2*i+1][2*j] * vCurr[i][j]) / params->deltaY)
                                        * params->deltaT + etaCurr[i][j];
                    }
                }
            }
            else{ // Process the last column separately from the rest because we need to use uReceived from the
                // the process with higher rank

                #pragma omp for schedule(static)
                for(int i = 0; i < size_X - 1; i++){
                    for(int j = 0; j < ySize + 1; j++){
                        etaNext[i][j] = (-(h[2*i+2][2*j+1] * uCurr[i+1][j] - h[2*i][2*j+1] * uCurr[i][j]) / params->deltaX 
                                        -(h[2*i+1][2*j+2] * vCurr[i][j+1] - h[2*i+1][2*j] * vCurr[i][j]) / params->deltaY)
                                        * params->deltaT + etaCurr[i][j];
                    }
                }
                #pragma omp for schedule(static)
                for(int j = 0; j < ySize + 1; j++){
                    etaNext[size_X-1][j] = (-(h[2*(size_X-1)+2][2*j+1] * uReceived[j] - h[2*(size_X-1)][2*j+1] * uCurr[size_X-1][j]) / params->deltaX 
                                    -(h[2*(size_X-1)+1][2*j+2] * vCurr[size_X-1][j+1] - h[2*(size_X-1)+1][2*j] * vCurr[size_X-1][j]) / params->deltaY)
                                    * params->deltaT + etaCurr[size_X-1][j];
                }
            }
        }

        // In a multiprocess environment, sending the rightmost column of eta of the domain controlled
        // by the current process to the process with the previous rank 
        if(nbproc != 1){
            if(myrank == 0){
                MPI_Send(etaCurr[size_X-1], ySize + 1, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD); //Tag 42 is for eta
            }else if (myrank == nbproc -1){
                MPI_Recv(etaReceived, ySize + 1, MPI_DOUBLE, myrank - 1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }else{
                MPI_Sendrecv(etaCurr[size_X-1], ySize + 1, MPI_DOUBLE, myrank + 1, 42,
                            etaReceived, ySize + 1, MPI_DOUBLE, myrank - 1, 42,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // uNext Boundary conditions
        if(myrank == 0 || nbproc == 1){
            for(int i = 0; i < ySize + 1; i++){
                uNext[0][i] = 0;
            }
        }
        if(myrank == nbproc -1 || nbproc == 1){
            for(int i = 0; i < ySize + 1; i++){
                uNext[size_X_u - 1][i] = 0;
            }
        }
    
        // Compute the next value of u
        #pragma omp parallel default(shared)
        {
            // Process uNext in one block
            if(nbproc == 1){
                #pragma omp for schedule(static)
                for(int i = 1; i < size_X_u-1; i++){
                    for(int j = 0; j < ySize + 1; j++){
                        uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                    -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                    }
                }
            }
            else if(myrank == 0){
                #pragma omp for schedule(static)
                for(int i = 1; i < size_X_u; i++){
                    for(int j = 0; j < ySize + 1; j++){
                        uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                    -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                    }
                }
            }
            else if(myrank == nbproc-1){
                // Process the first column separately from the rest because we need to use etaReceived from the
                // the process with lower rank
                // The last process has a smaller size along the x axis
                #pragma omp for schedule(static)
                for(int j = 0; j < ySize + 1; j++){
                    uNext[0][j] = (-params->g * (etaCurr[0][j] - etaReceived[j]) / params->deltaX
                                -params->gamma * uCurr[0][j]) * params->deltaT + uCurr[0][j];
                }
                #pragma omp for schedule(static)
                for(int i = 1; i < size_X_u-1; i++){ 
                    for(int j = 0; j < ySize + 1; j++){
                        uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                                        -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                    }
                }
            }
            else{
                // Process the first column separately from the rest because we need to use etaReceived from the
                // the process with lower rank
                #pragma omp for schedule(static)
                for(int j = 0; j < ySize + 1; j++){
                    uNext[0][j] = (-params->g * (etaCurr[0][j] - etaReceived[j]) / params->deltaX
                                -params->gamma * uCurr[0][j]) * params->deltaT + uCurr[0][j];
                }
                #pragma omp for schedule(static)
                for(int i = 1; i < size_X_u; i++){
                    for(int j = 0; j < ySize + 1; j++){
                        uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                    -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                    }
                }
            }
        }
        
        // Boundary conditions for v
        for(int i = 0; i < size_X; i++)
            vNext[i][0] = 0;

        // Setting the excitation on the rightmost column of the whole domain space
        for(int i = 0; i < size_X; i++){
            if(params->s == 0) //Sinusoidal excitation
                vNext[i][ySize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT);
            else // Exponentially decaying excitation
                vNext[i][ySize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT) * exp(- t * params->deltaT / 500);
        }


        // Compute the next value of v
        #pragma omp parallel default(shared)
        {
            #pragma omp for schedule(static)
            for(int i = 0; i < size_X; i++){
                for(int j = 1; j < ySize + 1; j++){
                    vNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i][j-1]) / params->deltaY
                                -params->gamma * vCurr[i][j]) * params->deltaT + vCurr[i][j];
                }
            }
        }

        // Process 0 gathers the sub-matrices of the processes and saves them to disk
        if(params->S != 0 && t % params->S == 0){
            gather_and_save(etaNext,uNext,vNext, xSize,ySize, debug, t, params);
        }

        // Go to next step
        double** tmp;
        
        tmp = (double**) etaCurr;
        etaCurr = (const double**) etaNext;
        etaNext = tmp;

        tmp = (double**) uCurr;
        uCurr = (const double**) uNext;
        uNext = tmp;

        tmp = (double**) vCurr;
        vCurr = (const double**) vNext;
        vNext = tmp;
    }   

    // Return values
    *eta = (double**) etaCurr;
    *u = (double**) uCurr;
    *v = (double**) vCurr;
    
    freeDoubleMatrix(etaNext, size_X,0);
    freeDoubleMatrix(uNext, size_X_u,0);
    freeDoubleMatrix(vNext, size_X,0);
    freeDoubleMatrix((double**) h, size_X_h,0);
    
    free(uReceived);
    free(etaReceived);

    return 0;
}

