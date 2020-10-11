#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int send(int localRank, int pnum,MPI_Comm comm, int size){
    int myRight = (localRank + 1) % pnum;
    int myLeft = (localRank + pnum - 1) % pnum;
    int* toSend = malloc(sizeof(int) * size);
    int* toRecv = malloc(sizeof(int)*size);
    for (int i = 0; i < size; i++){
        toSend[i] = 0;
        toRecv[i] = 0;
    }
    if (localRank != 0){
        MPI_Send(toSend, size, MPI_INT, myLeft, 0, comm);
        MPI_Recv(toRecv, size, MPI_INT, myLeft, 0, comm, MPI_STATUS_IGNORE);
    }
    if (localRank != pnum - 1){
        MPI_Recv(toRecv, size, MPI_INT, myRight, 0, comm, MPI_STATUS_IGNORE);  
        MPI_Send(toSend, size, MPI_INT, myRight, 0, comm);
    }

    return 0;
}



int main(int argc, char** argv){
    int size = atoi(argv[1]);

    MPI_Comm comm;
	int rank;
	int pnum;
	
    MPI_Init(&argc, &argv);

	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &pnum);

   clock_t startTime, endTime,duration; 
    startTime = clock(); 
    
    send(rank, pnum, comm, size);

    endTime = clock();
    duration = endTime - startTime; 
    double time_taken = ((double)duration)/CLOCKS_PER_SEC; // in seconds 
    
    printf("communication took %f seconds to execute \n", time_taken); 
    
    MPI_Finalize();

    return 0;
}