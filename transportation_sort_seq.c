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
    printf("\n");
    return 0;
}

int swap(int* array, int startPos, int endPos)
{
    int temp = array[startPos];
    array[startPos] = array[endPos];
    array[endPos] = temp;
    return 0;
}

int oddEvenSort(int* arr, int n){
    int isSorted = 0;

    while(!isSorted){
        isSorted = 1;
        for(int i = 1; i <= n-2; i = i+2){
            if(arr[i] > arr[i+1]){
                swap(arr, i, i+1);
                isSorted = 0;
            }
            if(arr[n-1] < arr[n-2]){
                swap(arr, n-1, n-2);
                isSorted = 0;
            }
        }

        for(int i = 0; i < n-2; i = i+2){
            if(arr[i] > arr[i+1]){
                swap(arr, i, i+1);
                isSorted = 0;
            }
        }
        
    }

    return isSorted;
}


int main(int argc, char** argv){
    if (argc != 2){
        return 0;
    }
    int allSize = atoi(argv[1]);
    int* allArray = malloc(sizeof(int) * allSize);
    generateArray(allArray,allSize);
    printf("Sequential Odd-Even Transportation Sort Demo, ID: %d \n", 117010008);
    printf("the array to sort is: \n");
    printArray(allArray,allSize);

    clock_t startTime, endTime,duration; 
    startTime = clock(); 
    oddEvenSort(allArray, allSize);
    endTime = clock();
    duration = endTime - startTime; 
    double time_taken = ((double)duration)/CLOCKS_PER_SEC; // in seconds 
    
    printf("the array after sort is: \n");
    printArray(allArray,allSize);

    printf("the sequential sort took %f seconds to execute \n", time_taken); 

    return 0;
}