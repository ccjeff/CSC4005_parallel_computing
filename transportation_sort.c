#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//helper functions
int generateArray(int* target,int size){
    srand(100);   // Initialization
    for (int i = 0; i < size; i++){
        target[i] = rand()%1000;
        // printf("the number of index %d is: %d \n",i,target[i]);
    }
    return 0;
}

int printArray(int *target, int size){
    for (int i = 0; i < size; i++){
        printf("%d,",target[i]);
    }
    return 0;
}

int swap(int* array, int startPos, int endPos)
{
    int temp = array[startPos];
    array[startPos] = array[endPos];
    array[endPos] = temp;

    return 0;
}


int oddEvenSort(int* localArray,const int chunkSize,const int localRank, const int globalSize,int pnum, MPI_Comm comm){
    
	int toSend = 0;
	int toRecv = 0;
    int myRight = (localRank + 1) % pnum;
    int myLeft = (localRank + pnum - 1) % pnum;

    printf("ODD_EVEN_SORT: the array to sort is: \n");
    printArray(localArray,chunkSize);
    printf("\n");

    for (int i = 0; i < globalSize; i++){
        if (i % 2 == 0){
            for (int j = chunkSize - 1; j > 0; j -= 2) {
                if (localArray[j-1] > localArray[j]){
                    swap(localArray, j, j-1);
                }
            }
        } else {
            for (int j = chunkSize - 2; j > 0; j -= 2) {
                if (localArray[j-1] > localArray[j]){
                    swap(localArray, j, j-1);
                }
            }
            if (localRank != 0){
                toSend = localArray[0];
				MPI_Send(&toSend, 1, MPI_INT, myLeft, 0, comm);
				MPI_Recv(&toRecv, 1, MPI_INT, myLeft, 0, comm, MPI_STATUS_IGNORE);
				if (toRecv > localArray[0]){
                    localArray[0] = toRecv;
                }
            }
            if (localRank != pnum - 1){
                toSend = localArray[chunkSize-1];
                MPI_Recv(&toRecv, 1, MPI_INT, myRight, 0, comm, MPI_STATUS_IGNORE);  
				MPI_Send(&toSend, 1, MPI_INT, myRight, 0, comm);
                if (toRecv < localArray[chunkSize-1]){
                    localArray[chunkSize-1] = toRecv;
                }
            }
        }
    }
	
    return 0;
}

int main(int argc, char** argv){
    // int* unsorted = generateArray(1000);
    if (argc != 2){
        printf("Please enter the array size in command line arguments!\n");
        return 0;
    }
    MPI_Comm comm;
	int rank;
	int pnum;
	int allSize = atoi(argv[1]);
	int chunkSize = 0;
    int isOdd;

    MPI_Init(&argc, &argv);

	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &pnum);

    int* globalArray = (int*) malloc(sizeof(int) * allSize);
    
    //divisability check
    int divisable = allSize%pnum;
	chunkSize = allSize / pnum;
    
    if (rank == 0){
        generateArray(globalArray,allSize);
        printf("Parallel Odd-Even Transportation Sort Demo, ID: %d \n", 117010008);
        printf("the array to sort is: \n");
        printArray(globalArray,allSize);
        printf("\n");
    }
    // int* localArray = (int*) malloc(sizeof(int) * chunkSize);
    int* localArray;
    
    clock_t startTime, endTime,duration; 
    startTime = clock(); 
    if (divisable == 0){
        localArray = (int*) malloc(sizeof(int) * chunkSize);
        MPI_Scatter(globalArray, chunkSize, MPI_INT, localArray, chunkSize, MPI_INT, 0, comm);
	    oddEvenSort(localArray,chunkSize, rank, allSize,pnum, comm);
	    MPI_Gather(localArray, chunkSize, MPI_INT, globalArray, chunkSize, MPI_INT, 0, comm);
    } else {
        int* displs = malloc(sizeof(int) * pnum);
        int* sendCounts = malloc(sizeof(int) * pnum);
        for (int i = 0; i < pnum - 1; i++){
            sendCounts[i] = chunkSize;
        }
        sendCounts[pnum - 1] = chunkSize + divisable;
        displs[0] = 0;
        for (int i = 1; i < pnum; i++){
            displs[i] = displs[i-1] + sendCounts[i-1];
        }
        if (rank == 0){
            printf("aaa\n");
            printArray(displs,pnum);
            printf("aaa\n");
            printArray(sendCounts,pnum);
            printf("aaa\n");

        }
        
        int recvCount = (rank == pnum - 1) ? chunkSize + divisable : chunkSize;
        localArray = (int*) malloc(sizeof(int) * recvCount);

        MPI_Scatterv(globalArray,sendCounts, displs, MPI_INT, localArray, recvCount, MPI_INT, 0, comm);
        oddEvenSort(localArray,sendCounts[rank], rank, allSize,pnum, comm);
        MPI_Gatherv(localArray, sendCounts[rank], MPI_INT, globalArray, sendCounts, displs, MPI_INT, 0, comm);

        // MPI_Scatterv(globalArray, sendCount,0, MPI_INT, localArray, recvCount, MPI_INT, 0, comm);
        // printf("sendCount for rank %d is: %d \n", rank, sendCount);
    }
    
    endTime = clock();
    duration = endTime - startTime; 
    double time_taken = ((double)duration)/CLOCKS_PER_SEC; // in seconds 
    
    if (rank == 0){
        printf("the array after sort is: \n");
        printArray(globalArray,allSize);
        printf("the parallel sort took %f seconds to execute \n", time_taken); 
    }
    free(localArray);

    MPI_Finalize();

    free(globalArray);

    return 0;
}