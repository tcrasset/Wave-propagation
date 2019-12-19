#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>

#include "project2_Delaunoy_Crasset_MAIN.h"
#include "project2_Delaunoy_Crasset_IO.h"

#define M_PI 3.14159265358979323846

SparseMatrix* toSparseMatrix(double** matrix, int xSize, int ySize){

    for(int i = 0; i < xSize; i++){
        for (int j = 0; j < ySize; j++){
            if(matrix[i][j] != 0){

            }
        }
        
    }

}


void get_array_sizes(int rank, int nbproc, int xSize, int* size_X, int* size_X_u, int* size_X_h, int* startval_X_h, int* endval_X_h){
    int mpi_xsize = xSize/nbproc;

    int startval_X, endval_X;
    int startval_X_u, endval_X_u;
    if(rank == 0){
        startval_X = 0;
        endval_X = mpi_xsize;
        *startval_X_h = 0;
        *endval_X_h = 2*mpi_xsize + 2;
        startval_X_u = 0;
        endval_X_u = mpi_xsize;
    }else if(rank == nbproc -1){
        startval_X = rank * mpi_xsize + 1;
        endval_X = (rank+1) * mpi_xsize;
        // startval_X_h = 2 * rank * mpi_xsize + 3;
        *startval_X_h = 2 * rank * mpi_xsize + 2; //Include this line so as to not transfer h
        *endval_X_h = 2 * (rank+1) * mpi_xsize + 2;
        startval_X_u = rank * mpi_xsize + 1;
        endval_X_u = (rank+1) * mpi_xsize + 1;
    }else{
        startval_X = rank * mpi_xsize + 1;
        endval_X = (rank+1) * mpi_xsize; 
        // startval_X_h = 2 * rank * mpi_xsize + 3;
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



int gather_and_save(double** eta, double**  u, double**  v, int xSize, int ySize,  int debug, unsigned int iteration, Parameters* params){

    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int size_X, size_X_u, size_X_h, startval_X_h, endval_X_h;
    get_array_sizes(myrank, nbproc, xSize, &size_X, &size_X_u, &size_X_h, &startval_X_h, &endval_X_h);


    double* etaPartial = transformMatrixToArray(eta, size_X, ySize +1);
    double* uPartial = transformMatrixToArray(u, size_X_u, ySize +1);
    double* vPartial = transformMatrixToArray(v, size_X, ySize +2);

    if(debug == 1){

        fprintf(stderr, "process %d size_X = %d\n", myrank, size_X);
        fprintf(stderr, "process %d size_X_u = %d\n", myrank, size_X_u);
        fprintf(stderr, "process %d xSize = %d\n", myrank, xSize);
        fprintf(stderr, "process %d ySize = %d\n", myrank, ySize);

        fprintf(stderr, "process %d begin test etaPartial\n", myrank);
        fprintf(stderr, "process %d etaPartial = %lf\n", myrank, etaPartial[size_X * (ySize + 1) - 1]);
        fprintf(stderr, "process %d end test etaPartial\n", myrank);
        fprintf(stderr, "process %d begin test uPartial\n", myrank);
        fprintf(stderr, "process %d uPartial = %lf\n", myrank, uPartial[size_X_u * (ySize + 1) - 1]);
        fprintf(stderr, "process %d end test uPartial\n", myrank);
        fprintf(stderr, "process %d begin test vPartial\n", myrank);
        fprintf(stderr, "process %d vPartial = %lf\n", myrank, vPartial[size_X * (ySize + 2) - 1]);
        fprintf(stderr, "process %d end test vPartial\n", myrank);
    }
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

    /*
    for (int i = 0; i < nbproc; i++){
        fprintf(stderr, "recvcount % d = %d \n", i, recvcounts_v[i]);
        fprintf(stderr, "disp % d = %d \n", i, disp_v[i]);
    }
    */

    double * etaTotal = malloc((xSize + 1) * (ySize  + 1)* sizeof(double));
    double * uTotal = malloc((xSize + 2) * (ySize  + 1)* sizeof(double));
    double * vTotal = malloc((xSize + 1) * (ySize  + 2)* sizeof(double));
    MPI_Gatherv(etaPartial, (size_X) * (ySize + 1) , MPI_DOUBLE, etaTotal, recvcounts_eta, disp_eta, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(uPartial, (size_X_u) * (ySize + 1) , MPI_DOUBLE, uTotal, recvcounts_u, disp_u, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(vPartial, (size_X) * (ySize + 2) , MPI_DOUBLE, vTotal, recvcounts_v, disp_v, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(recvcounts_eta);
    free(recvcounts_u);
    free(recvcounts_v);
    free(disp_eta);
    free(disp_u);
    free(disp_v);

    /*
    if(myrank == 0){
        fprintf(stderr,"***********ETA TOTAL**************\n");
        printLinearArray(etaTotal, xSize + 1, ySize +1);
        fprintf(stderr,"***********U TOTAL**************\n");
        printLinearArray(uTotal, xSize + 2, ySize +1);
        fprintf(stderr,"***********V TOTAL**************\n");
        printLinearArray(vTotal, xSize + 1, ySize +2);
    }
    */

    if(myrank == 0){
        for(int i = 0; i < (xSize + 1) * (ySize + 1); i++){
            if(etaTotal[i] != 0.0)
                fprintf(stderr, "etaTotal = %lf\n", etaTotal[i]);
        }
    }

    if(debug == 1 && myrank == 0){
        fprintf(stderr,"***********ETA TOTAL**************\n");
        printLinearArray(etaTotal, xSize + 1, ySize +1);
        fprintf(stderr,"***********U TOTAL**************\n");
        printLinearArray(uTotal, xSize + 2, ySize +1);
        fprintf(stderr,"***********V TOTAL**************\n");
        printLinearArray(vTotal, xSize + 1, ySize +2);
    }

    if(myrank == 0){
        saveToDisk(etaTotal, uTotal, vTotal, xSize, ySize, iteration, params);
    }

    free(etaTotal);
    free(uTotal);
    free(vTotal);
}

int eulerExplicitMPI(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug, int debug_rank){
    
    // int i = 0;
    // char hostname[256];
    // gethostname(hostname, sizeof(hostname));
    // printf("PID %d on %s ready for attach\n", getpid(), hostname);
    // fflush(stderr);
    // while (0 == i)
    // sleep(5);

    assert(map);
    assert(params);

    int nbproc, myrank ;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);

    /*
    fprintf(stderr, "map->a = %lf\n", map->a);
    fprintf(stderr, "map->b = %lf\n", map->b);
    */

    int xSize = (int)(map->a / params->deltaX);
    int ySize = (int)(map->b / params->deltaY);


    int size_X;
    int size_X_u;
    int size_X_h;
    int startval_X_h;
    int endval_X_h;
    get_array_sizes(myrank, nbproc, xSize, &size_X, &size_X_u, &size_X_h, &startval_X_h, &endval_X_h);

    /*
    fprintf(stderr, "process %d size_X = %d\n", myrank, size_X);
    fprintf(stderr, "process %d size_X_u = %d\n", myrank, size_X_u);
    fprintf(stderr, "process %d size_X_h = %d\n", myrank, size_X_h);
    */

    if(debug == 1){
        fprintf(stderr, "process %d size_X = %d\n", myrank, size_X);
        fprintf(stderr, "process %d size_X_u = %d\n", myrank, size_X_u);
        fprintf(stderr, "process %d size_X_h = %d\n", myrank, size_X_h);
    }

    if(debug == 1){
        fprintf(stderr,"Processs %d \n xSize %d \n Size X u %d\n Size X_h %d\n Size_X %d \n", myrank, xSize, size_X_u, size_X_h, size_X);
    }

    // Allocate memory

    // eta in {0, 1, ..., a/dx}X{0, 1, ..., b/dy}
    double** etaCurr = allocateDoubleMatrix(size_X, ySize + 1);
    if(!etaCurr){
        return -1;
    }

    double** etaNext = allocateDoubleMatrix(size_X, ySize + 1);
    if(!etaNext){
        freeDoubleMatrix(etaCurr, size_X,0);
        return -1;
    }

    // u in {-1/2, 1/2, ..., a/dx + 1/2}X{0, 1, ..., b/dy}
    double** uCurr = allocateDoubleMatrix(size_X_u, ySize + 1);
    if(!uCurr){
        freeDoubleMatrix(etaCurr,size_X,0);
        freeDoubleMatrix(etaNext,size_X,0);
        return -1;
    }

    double** uNext = allocateDoubleMatrix(size_X_u, ySize + 1);
    if(!uNext){
        freeDoubleMatrix(etaCurr,size_X,0);
        freeDoubleMatrix(etaNext,size_X,0);
        freeDoubleMatrix(uCurr, size_X_u,0);
        return -1;
    }

    // v in {0, 1, .., a/dx}X{-1/2, 1/2, ..., b/dy + 1/2}
    double** vCurr = allocateDoubleMatrix(size_X, ySize + 2);
    if(!vCurr){
        freeDoubleMatrix(etaCurr, size_X,0);
        freeDoubleMatrix(etaNext, size_X,0);
        freeDoubleMatrix(uCurr, size_X_u,0);
        freeDoubleMatrix(uNext, size_X_u,0);
        return -1;
    }

    double** vNext = allocateDoubleMatrix(size_X, ySize + 2);
    if(!vNext){
        freeDoubleMatrix(etaCurr, size_X,0);
        freeDoubleMatrix(etaNext, size_X,0);
        freeDoubleMatrix(uCurr, size_X_u,0);
        freeDoubleMatrix(uNext, size_X_u,0);
        freeDoubleMatrix(vCurr, size_X,0);
        return -1;
    }

    if(debug == 1){
        fprintf(stderr, "Process %d before vnext\n", myrank);
    }
    
    if(debug == 1){
        fprintf(stderr, "Process %d after vnext\n", myrank);
    }
    // h in {-1/2, 0, 1/2, ..., a/dx, a/dx + 1/2}X{-1/2, 0, 1/2, ..., b/dy, b/dy + 1/2}
    double** h = allocateDoubleMatrix(size_X_h, 2 * ySize + 3);
    if(!h){
        freeDoubleMatrix(etaCurr, size_X,0);
        freeDoubleMatrix(etaNext, size_X,0);
        freeDoubleMatrix(uCurr, size_X_u,0);
        freeDoubleMatrix(uNext, size_X_u,0);
        freeDoubleMatrix(vCurr, size_X,0);
        freeDoubleMatrix(vNext, size_X,0);
        return -1;
    }

    for(int i = startval_X_h; i <= endval_X_h; i++){
        for(int j = 0; j < 2 * ySize + 3; j++){
            h[i-startval_X_h][j] = getGridValueAtDomainCoordinates(map, ((float)(i * xSize)/(xSize + 1)) * (params->deltaX / 2), ((float)(j * ySize)/(ySize + 1)) * (params->deltaY / 2));
        }
    }

    /*
    if(myrank == 0){
        fprintf(stderr, "Process %i \n", myrank);
        printDoubleMatrix(h, endval_X_h - startval_X_h + 1, 2*ySize + 3, myrank);
    }
    */

    if(debug == 1){
        fprintf(stderr, "Process %d Does not fail at h fill\n", myrank);
    }

    // if(debug == 1 && myrank == debug_rank){
    //     printf("*************Process %d *******************\n", myrank);
    //     printf("h\n");
    //     for(int i = startval_X_h; i < endval_X_h; i++){
    //         for(int j = 0; j < 2 * ySize + 3; j++){
    //             fprintf(stderr, "%lf ",h[i][j]);
    //         }
    //         fprintf(stderr, "\n");
    //     }
    //     printf("apres h\n");
    // }


    for(int i = 0; i < size_X; i++){
        for(int j = 0; j < ySize + 1; j++){
            etaCurr[i][j] = 0;
        }
    }

    for(int i = 0; i < size_X_u; i++){
        for(int j = 0; j < ySize + 1; j++){
            uCurr[i][j] = 0;
        }
    }

    for(int i = 0; i < size_X; i++){
        for(int j = 0; j < ySize + 2; j++)
            vCurr[i][j] = 0;
    }

    double* uReceived = malloc((ySize + 1) * sizeof(double));
    double* etaReceived = malloc((ySize + 1) * sizeof(double));

    for(unsigned int t = 1; t <= params->TMax/params->deltaT; t++){

        //fprintf(stderr, "in loop t = %u\n", t);

        if(debug == 1){
            fprintf(stderr, "Process%d begin loop %d/%f\n", myrank, t, params->TMax/params->deltaT);
        }

        if(debug == 1 && myrank == debug_rank){
            printf("************* Process %d *******************\n", myrank);
            printf("t = %u\n", t);
        }

        if(debug == 1){
            fprintf(stderr, "Process %d Fails before etaNext\n", myrank);
        }

        if(debug == 1){
            fprintf(stderr, "Process %d allocated uReceived\n", myrank);
        }

        if(myrank == nbproc-1){
            MPI_Send(uCurr[0], ySize + 1, MPI_DOUBLE, myrank - 1, 42, MPI_COMM_WORLD); //Tag 42 is for eta
        }else if (myrank == 0){
            MPI_Recv(uReceived, ySize + 1, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }else{
            MPI_Sendrecv(uCurr[0], ySize + 1, MPI_DOUBLE, myrank - 1, 42,
                         uReceived, ySize + 1, MPI_DOUBLE, myrank + 1, 42,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if(myrank == nbproc-1){
            // Process etaNext in one block
            for(int i = 0; i < size_X; i++){
                for(int j = 0; j < ySize + 1; j++){
                    etaNext[i][j] = (-(h[2*i+2][2*j+1] * uCurr[i+1][j] - h[2*i][2*j+1] * uCurr[i][j]) / params->deltaX 
                                    -(h[2*i+1][2*j+2] * vCurr[i][j+1] - h[2*i+1][2*j] * vCurr[i][j]) / params->deltaY)
                                    * params->deltaT + etaCurr[i][j];
                }
            }
        }
        else{
            for(int i = 0; i < size_X - 1; i++){
                for(int j = 0; j < ySize + 1; j++){
                    etaNext[i][j] = (-(h[2*i+2][2*j+1] * uCurr[i+1][j] - h[2*i][2*j+1] * uCurr[i][j]) / params->deltaX 
                                    -(h[2*i+1][2*j+2] * vCurr[i][j+1] - h[2*i+1][2*j] * vCurr[i][j]) / params->deltaY)
                                    * params->deltaT + etaCurr[i][j];
                }
            }
            for(int j = 0; j < ySize + 1; j++){
                etaNext[size_X-1][j] = (-(h[2*(size_X-1)+2][2*j+1] * uReceived[j] - h[2*(size_X-1)][2*j+1] * uCurr[size_X-1][j]) / params->deltaX 
                                -(h[2*(size_X-1)+1][2*j+2] * vCurr[size_X-1][j+1] - h[2*(size_X-1)+1][2*j] * vCurr[size_X-1][j]) / params->deltaY)
                                * params->deltaT + etaCurr[size_X-1][j];
            }
        }

        if(debug == 1){
            fprintf(stderr, "Process %d Does not fail at etaNext\n", myrank);
        }
        // Compute uNext
        if(myrank == 0){
            for(int i = 0; i < ySize + 1; i++){
                uNext[0][i] = 0;
            }
        }
        else if(myrank == nbproc -1){
            for(int i = 0; i < ySize + 1; i++){
                uNext[size_X_u - 1][i] = 0;
            }
        }

        if(debug == 1){
            fprintf(stderr, "Process %d allocated etaReceived\n", myrank);
        }

        if(myrank == 0){
            MPI_Send(etaCurr[size_X-1], ySize + 1, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD); //Tag 42 is for eta
        }else if (myrank == nbproc -1){
            MPI_Recv(etaReceived, ySize + 1, MPI_DOUBLE, myrank - 1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }else{
            MPI_Sendrecv(etaCurr[size_X-1], ySize + 1, MPI_DOUBLE, myrank + 1, 42,
                         etaReceived, ySize + 1, MPI_DOUBLE, myrank - 1, 42,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if(debug == 1){

            fprintf(stderr,"*****************PROCESS %d *****************\n",myrank);
            /*
            for(int j = 0; j < ySize + 1; j++){
                fprintf(stderr, "%lf ", etaReceived[j]);
            }
            */
            fprintf(stderr, "\n");
        }

        if(debug == 1){
            fprintf(stderr, "Process %d start unext\n", myrank);
        }
        
        // Process uNext in one block
        if(myrank == 0){
            for(int i = 1; i < size_X_u; i++){
                for(int j = 0; j < ySize + 1; j++){
                    uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                }
            }
        }
        else if(myrank == nbproc-1){
            for(int j = 0; j < ySize + 1; j++){
                uNext[0][j] = (-params->g * (etaCurr[0][j] - etaReceived[j]) / params->deltaX
                               -params->gamma * uCurr[0][j]) * params->deltaT + uCurr[0][j];
            }
            for(int i = 1; i < size_X_u-1; i++){ 
                for(int j = 0; j < ySize + 1; j++){
                    /*
                    fprintf(stderr, "eta curr = %lf\n", etaCurr[i][j]);
                    fprintf(stderr, "eta curr -1 = %lf\n", etaCurr[i-1][j]);
                    fprintf(stderr, "eta curr = %lf\n", uCurr[i][j]);
                    */
                    uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                   -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                    //fprintf(stderr, "uNext = %lf\n", uNext[i][j]);
                }
            }
        }

        // Process uNext[0] alone because of etaReceived in one block
        else{
            
            for(int j = 0; j < ySize + 1; j++){
                uNext[0][j] = (-params->g * (etaCurr[0][j] - etaReceived[j]) / params->deltaX
                               -params->gamma * uCurr[0][j]) * params->deltaT + uCurr[0][j];
            }
            for(int i = 1; i < size_X_u; i++){ // Shouldn't that be size_X_u ? Or is it because one starts at 1 and not 0
                for(int j = 0; j < ySize + 1; j++){
                    uNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i-1][j]) / params->deltaX
                                -params->gamma * uCurr[i][j]) * params->deltaT + uCurr[i][j];
                }
            }
        }
        if(debug == 1){
            fprintf(stderr, "Process %d end unext\n", myrank);
            fprintf(stderr, "Process %d start vnext\n", myrank);
        }

        for(int i = 0; i < size_X; i++)
            vNext[i][0] = 0;

        if(debug == 1)
            fprintf(stderr, "Process %d after first v loop\n", myrank);

        for(int i = 0; i < size_X; i++){
            if(params->s == 0)
                vNext[i][ySize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT);
            else
                vNext[i][ySize+1] = params->A * sin(2 * M_PI * params->f * t * params->deltaT) * exp(- t * params->deltaT / 500);

            //fprintf(stderr, "process %d vnext = %lf", myrank, vNext[i][ySize + 1]);
        }
        if(debug == 1)
            fprintf(stderr, "Process %d before v loop \n", myrank);
        for(int i = 0; i < size_X; i++){
            for(int j = 1; j < ySize + 1; j++){
                vNext[i][j] = (-params->g * (etaCurr[i][j] - etaCurr[i][j-1]) / params->deltaY
                               -params->gamma * vCurr[i][j]) * params->deltaT + vCurr[i][j];
            }
        }

        if(debug == 1)
            fprintf(stderr, "Process %d end vnext\n", myrank);
        

        if(debug == 1 && myrank == debug_rank){
            printf("\n\n\n*************Process %d *******************\n\n\n\n", myrank);
            printf("etaCurr\n");
            printDoubleMatrix(etaCurr, size_X+1, ySize + 1, myrank);
            printf("etaNext\n");
            printDoubleMatrix(etaNext, size_X+1, ySize + 1,myrank);
            printf("uCurr\n");
            printDoubleMatrix(uCurr, size_X_u+1, ySize + 1,myrank);
            printf("uNext\n");
            printDoubleMatrix(uNext, size_X_u+1, ySize + 1,myrank);
            printf("vCurr\n");
            printDoubleMatrix(vCurr, size_X+1, ySize + 2,myrank);
            printf("vNext\n");
            printDoubleMatrix(vNext, size_X+1, ySize + 2,myrank);
        }

        /*
        printf("\n\n\n*************Process %d *******************\n\n\n\n", myrank);
        printf("etaCurr\n");
        printDoubleMatrix(etaCurr, size_X, ySize + 1, myrank);
        printf("etaNext\n");
        printDoubleMatrix(etaNext, size_X, ySize + 1,myrank);
        printf("uCurr\n");
        printDoubleMatrix(uCurr, size_X_u, ySize + 1,myrank);
        printf("uNext\n");
        printDoubleMatrix(uNext, size_X_u, ySize + 1,myrank);
        printf("vCurr\n");
        printDoubleMatrix(vCurr, size_X, ySize + 2,myrank);
        printf("vNext\n");
        printDoubleMatrix(vNext, size_X, ySize + 2,myrank);
        */


        // Process 0 saves arrays to disk
        if(params->S != 0 && t % params->S == 0){
            // Gather the matrices and save to disk
            gather_and_save(etaNext,uNext,vNext, xSize,ySize, debug, t, params);
        }

        // Go to next step
        double** tmp;
        
        tmp = etaCurr;
        etaCurr = etaNext;
        etaNext = tmp;

        tmp = uCurr;
        uCurr = uNext;
        uNext = tmp;

        tmp = vCurr;
        vCurr = vNext;
        vNext = tmp;

        if(debug == 1){
            fprintf(stderr, "Process%d end loop %d/%lf\n", myrank, t, params->TMax);
        }
    }   

    *eta = etaCurr;
    *u = uCurr;
    *v = vCurr;
    
    freeDoubleMatrix(etaNext, size_X,0);
    freeDoubleMatrix(uNext, size_X_u,0);
    freeDoubleMatrix(vNext, size_X,0);
    freeDoubleMatrix(h, size_X_h,0);
    
    free(uReceived);
    free(etaReceived);

    return 0;
}

inline double dotProduct(double* x, double* y, unsigned int size){
    double result = 0.0;

    for(unsigned int i = 0; i < size; i++){
        result += x[i] * y[i];
    }

    return result;
}

inline double vectorNorm(double* x, unsigned int size){
    double norm = 0.0;

    for(unsigned int i = 0; i < size; i++){
        norm += x[i] * x[i];
    }

    return norm;
}

int EulerImplicit(Map* map, Parameters* params, double*** eta, double*** u, double*** v, int debug, int debug_rank){

}

int main(int argc, char* argv[]) {

    // Check number of arguments
    int debug, debug_rank;
    assert(argc >= 4);
    const char* parameter_file = argv[1];
    const char* map_file = argv[2];
    const unsigned int scheme = atoi(argv[3]);

    if(argc == 4)
        debug = 0;
    else if(argc == 5){
        fprintf(stderr, "Call with rank of process to debug as 6th argument \n");
        return -1;
    }else{
        assert(argc == 6);
        debug = atoi(argv[4]);
        debug_rank = atoi(argv[5]);
    }

    //Check argument validity
    assert((scheme == 0) || (scheme == 1));

    MPI_Init(&argc,&argv);
    
    // Measure execution time
    double startTime = MPI_Wtime();

    int nbproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    Parameters* params = readParameterFile(parameter_file);
    Map* map = readMapFile(map_file, 0);
    // Map* map = readMapFile("test_map.dat", 0);
    if(debug == 1)
        printDoubleMatrix(map->grid, map->X, map->Y, myrank);

    // Explicit
    if (scheme == 0) {
        double** eta;
        double** u;
        double** v;


        if(eulerExplicitMPI(map, params, &eta, &u, &v, debug, debug_rank) == -1){

            fprintf(stderr, "error in euler function\n");
            free(params);
            free(map->grid);
            free(map);
            MPI_Finalize();
            exit(EXIT_FAILURE);
            
        }

        int xSize = (int)(map->a / params->deltaX);
        int ySize = (int)(map->b / params->deltaY);

        int size_X, size_X_u, size_X_h, startval_X_h, endval_X_h;

        get_array_sizes(myrank, nbproc, xSize, &size_X, &size_X_u, &size_X_h, &startval_X_h, &endval_X_h);

        //MPI
        if(debug == 1){
            fprintf(stderr,"***********ETA**************\n");
            printDoubleMatrix(eta, size_X, ySize + 1,myrank);
            fprintf(stderr,"***********U**************\n");
            printDoubleMatrix(u, size_X_u, ySize + 1,myrank);
            fprintf(stderr,"***********V**************\n");
            printDoubleMatrix(v, size_X, ySize + 2,myrank);
        }

        // Gather the matrices and save to disk
        gather_and_save(eta,u,v, xSize,ySize, debug, params->TMax/params->deltaT + 1, params);
        

        if(debug == 1){
            fprintf(stderr, "process %d before free matrix\n", myrank);
            freeDoubleMatrix(eta, size_X, 0);
            freeDoubleMatrix(u, size_X_u, 0);
            freeDoubleMatrix(v, size_X, 0);
            fprintf(stderr, "process %d after free matrix\n", myrank);
        }

    }
    // Implicit
    else {
        printf("Implicit ");
        printf("%s %s %u", parameter_file, map_file, scheme);
    }


    if(debug == 1){
        fprintf(stderr, "before free params\n");
    }
    free(params);
    freeDoubleMatrix(map->grid, map->X,0);
    free(map);
    if(debug == 1){
        fprintf(stderr, "after free params\n");
    }    

    // Mesure execution time
    double endTime = MPI_Wtime();
    double executionTime = endTime - startTime;
    char* openMP_nbthreads = getenv("OMP_NUM_THREADS");

    // Print statistics to standard output
    fprintf(stdout,"%d,%d,%d,%s,%lf,%lf,%lf,%lf,%u,%lf\n", scheme, myrank, nbproc, \
                openMP_nbthreads, executionTime, params->deltaX, params->deltaY, \
                params->deltaT, params->s, params->r_threshold);
    MPI_Finalize();

    return 0;

}
