#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "ants.h"

#define MAX_ANTS 16

struct ant_t{
 float tourLength;
 int curNode, nextNode, pathIndex;
 //int path[MAX_NODES];
 //int tabu[MAX_NODES];
};




int colony_size;
int *res;

int *getCalcs(){
  srand (getpid());
  int *result = new int[colony_size];
  for (int i=0; i<colony_size; i++){
    result[i] = rand()%10;
    
  }
  return result;
}

int main(int argc, char* argv[]){
  MPI_Init(NULL, NULL);
 
  //Create  Datatype ant_type
  MPI_Datatype ant_type;
  int lengths[4] = {1, 1, 1, 1};
  const MPI_Aint displacements[6]={0, sizeof(int), sizeof(int)+sizeof(int), sizeof(int)*3};
  MPI_Datatype types[4] = {MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Type_create_struct(4, lengths, displacements, types, &ant_type);
  MPI_Type_commit(&ant_type);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  colony_size = MAX_ANTS/world_size;

  res = getCalcs();  

  if (world_rank != 0){
   struct ant_t buffer;
   buffer.tourLength = 2.0;
   buffer.curNode = 10;
   buffer.nextNode = 20;
   buffer.pathIndex = 1;
   
   // int hormigas[colony_size];
   //MPI_Send(res, colony_size, MPI_INT, 0, 1, MPI_COMM_WORLD);
   MPI_Send(&buffer, 1, ant_type, 0,1, MPI_COMM_WORLD);
   printf("Se envia: ");/*
   /
   for (int i = 0; i<colony_size; i++){
    printf("%d, ", res[i]);
   }
   printf("\n");*/
  }

  else {
   int container[world_size-1][colony_size];
   struct ant_t buffer_r;
   for (int i=1; i<world_size; i++){
     	    
   //MPI_Recv(container[i-1], colony_size, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&buffer_r, 1, ant_type, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }
  
   for (int j=1; j<world_size; j++){
    printf("Se recibieron estos datos del process: %d: \n",j);
    printf("Data --> %d", buffer_r.curNode);}

   /*
    for (int k=0; k<colony_size; k++){
      printf("%d ", container[j-1][k]);
   } 
   
  }*/
  }  
  MPI_Finalize();


}
