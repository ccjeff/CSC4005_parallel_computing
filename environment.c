#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main(int argc, char** argv){
    // int* unsorted = generateArray(1000);
    MPI_Comm comm;
    int rank;
    int pnum;


    clock_t startTime, endTime,duration;
    startTime = clock(); 
    MPI_Init(&argc, &argv);


    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &pnum);

    MPI_Finalize();

    endTime = clock();
    duration = endTime - startTime; 
    double time_taken = ((double)duration)/CLOCKS_PER_SEC; // in seconds 
    
    printf("start up time cost: %f\n", time_taken);

    return 0;
}