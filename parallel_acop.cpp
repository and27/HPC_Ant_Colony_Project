#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string>
#include "ants.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
//#define MAX_ANTS 16

using namespace std; 

struct ant_t{
 float tourLength;
 int curNode, nextNode, pathIndex;
 int path[MAX_NODES];
 int tabu[MAX_NODES];
};


int colony_size;
int *res;
ant_t *antColony;
float pheromone[MAX_NODES][MAX_NODES];
float best = (float)MAX_TOUR;
int bestIndex;
//int MAX_ANTS;
cityType cities[MAX_NODES];
string graph;
EdgeMatrix *dist;

//struct ant_t container_buffer_r[][];

//this can be implemented as a function member or constructor if we construct the ant class
void initAnt(){
  for (int from = 0; from < MAX_NODES; from ++){
   for (int to=0; to < MAX_NODES; to ++){
     pheromone[from][to] = INIT_PHER;
   }
  }

  //Initialize the departure node on all the antColony
  for (int ant=0; ant<colony_size; ant++){
    antColony[ant].curNode = rand()%MAX_NODES;
    //Put all the cities available (tabu = false (0))
    for (int to = 0; to <MAX_NODES; to++){
      antColony[ant].tabu[to]=0;
      }
    antColony[ant].pathIndex = 0; 
    antColony[ant].path[0] = antColony[ant].curNode;
    antColony[ant].nextNode = -1;
    antColony[ant].tourLength=0;
    //put the selected city as visited (tabu = true (1))
    antColony[ant].tabu[antColony[ant].curNode] = 1; 
  }
}


void restartAnts() {
  for (int ant = 0; ant < colony_size; ant++) {
    //At the begining we put the nextNode as disabled (-1)	  
    antColony[ant].nextNode = -1;
    antColony[ant].tourLength = 0.0;
    for (int i = 0; i < MAX_NODES; i++) {
      antColony[ant].tabu[i] = 0;
      antColony[ant].path[i] = -1;
    }
    //Select the first node randomly
    antColony[ant].curNode = rand() % MAX_NODES;
    antColony[ant].pathIndex = 1;
    antColony[ant].path[0] = antColony[ant].curNode;
    antColony[ant].tabu[antColony[ant].curNode] = 1;
  }
}

//Calculate node selection probability
float calculateProbability(int from, int to){
   
   //cout << " from "<<from<<" to " <<to <<"there is " <<(*dist)[from][to]<<endl;
   //cout << (pow(pheromone[from][to], ALPHA)* pow((1.0/(*dist)[from][to]), BETA));
   return (pow(pheromone[from][to], ALPHA) * pow((1.0/(*dist)[from][to]), BETA));
}

//Select the next node using the calculateProbability function (p)
int selectNextCity(int ant){
  int from = antColony[ant].curNode;
  float sum = 0.0;
  for (int to = 0; to <MAX_NODES; to++){
    if (antColony[ant].tabu[to]==0){
	    
      sum+=calculateProbability(from, to);
    }
   }
  //cout << "the sum is : "<< sum <<endl;

  int lastBestIndex = 0.0; 
  float luckyNumber = (float)rand()/RAND_MAX;
  //cout << "random: " << luckyNumber<<endl; 

  for (int to = 0; to < MAX_NODES; to++) {
    if (antColony[ant].tabu[to]==0){

      float product = calculateProbability(from, to) /sum;
      sleep(0.5);
      if(product > 0 ){
       // cout << "probability: "  << product <<endl; 
       luckyNumber-= product; 
       lastBestIndex=to;
       	
       if(luckyNumber<=0.0){
	  //cout << "Selected node "<< to<<endl;
        return to;
       }
      }
    }
   }
  cout << "not reached" << endl; 
  return lastBestIndex;
 }


void updatePheromone(ant_t *antColony){
 int from, to, i, ant;
 for (from = 0; from <MAX_NODES; from ++){
  for (to=0; to<from; to++){
   pheromone[from][to]*=1.0-RHO;
   
   if (pheromone[from][to]<0.0){
    pheromone[from][to]=INIT_PHER;
   }
   pheromone[to][from] = pheromone[from][to];
  }
  
 }
 

  for (ant = 0; ant < colony_size; ant++) {
    for (i = 0; i < MAX_NODES; i++) {
	    
     from = antColony[ant].path[i];
      //cout << "from is : " << from;

      if (i < MAX_NODES - 1) {
      to = antColony[ant].path[i+1];
      //cout << " to is " <<to<<endl;
      } else {
      to = antColony[ant].path[0];
      //cout << " to is " <<to<<endl;
       }
      pheromone[from][to] += (QVAL / antColony[ant].tourLength);
      pheromone[to][from] = pheromone[from][to];
    }
     }
} 

float euclideanDistance(int x1, int x2, int y1, int y2) {
  int xd = x1 - x2;
  int yd = y1 - y2;
  return (int) (sqrt(xd * xd + yd * yd) + 0.5);
}

float pseudoEuclideanDistance(int x1, int x2, int y1, int y2) {
  int xd = x1 - x2;
  int yd = y1 - y2;
  float rij = sqrt((xd * xd + yd * yd) / 10.0); 
  return ceil(rij);
}


void constructTSP(string graph, cityType *cities, EdgeMatrix *dist) {
  // Load cities from file
  ifstream infile(("instances/"+graph + ".tsp").c_str());
  string line;
  bool euclidean = true; // whether to run EUC_2D or ATT distance metric

  int city;
  int x, y;
  bool reading_nodes = false;
  while (std::getline(infile, line)) {
    istringstream iss(line);
    string word;
    if (!reading_nodes) {
      cout << "Reading lines " <<endl;
      iss >> word;
      if (word.compare("EDGE_WEIGHT_TYPE") == 0) {
        iss >> word >> word;
        std::cout << "edge type: " << word << std::endl;
        euclidean = !word.compare("EUC_2D");
      } else if (word.compare("NODE_COORD_SECTION") == 0) {
        reading_nodes = true;
	cout << "yeah baby " <<endl;
      }
    } else if (iss >> city >> x >> y) {
      cout << " x value is: " << x << " " << endl;
      cities[city-1].x = x;
      cities[city-1].y = y;
    }
  }
  infile.close();

  // Compute distances between cities (edge weights)
  for (int from = 0; from < MAX_NODES; from++) {
    (*dist)[from][from] = 0.0;

    for (int to = from + 1; to < MAX_NODES; to++) {
      float edge_dist;
      if (euclidean) {
	//cout << "cities from x " << cities[from].x << "." <<endl;
        edge_dist = euclideanDistance(cities[from].x, cities[to].x, cities[from].y, cities[to].y);
      } else {
        edge_dist = pseudoEuclideanDistance(cities[from].x, cities[to].x, cities[from].y, cities[to].y);
      }
      if (edge_dist == 0) {
        edge_dist = 1.0;
      }
      //printf("edge[%d][%d] = %f\n", from, to, edge_dist);
      (*dist)[from][to] = edge_dist;
      (*dist)[to][from] = edge_dist;
    }
  }
  printf("Graph constructed!\n");
}
  
//Select the best found solution
float bestSolution(ant_t *container){	 
  for (int ant = 0; ant < colony_size; ant++) {
    if (container[ant].tourLength < best) {
      best = container[ant].tourLength;
      printf("new best: %1.f\n", best);
      bestIndex = ant;
    }
  }
 
  return best, bestIndex;
}
 
int solveTSP(){
  cout << "Initializing MPI ... " <<endl;
  sleep(1);
  //Create  Datatype ant_type
  MPI_Datatype ant_type;
  int lengths[6] = {1, 1, 1, 1, MAX_NODES, MAX_NODES};
  const MPI_Aint displacements[6]={0, sizeof(float), sizeof(float)+sizeof(int), sizeof(int)*2+sizeof(float), sizeof(int)*3+sizeof(float), sizeof(int)*MAX_NODES+sizeof(int)*3+sizeof(float)};
  MPI_Datatype types[6] = {MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Type_create_struct(6, lengths, displacements, types, &ant_type);
  MPI_Type_commit(&ant_type);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  colony_size = MAX_ANTS/world_size;

  antColony = new ant_t[colony_size];
  initAnt();
 
 cout << "Initializing ants and running 5 generations" <<endl;

 for (int a = 0; a<5; a++){
  if (world_rank != 0){

     cout << "------------"<<endl;
     cout << "Running iteration: "<<a<<endl;
     cout << "------------"<<endl;

     cout << "Pheromone value at the iteration start: " << pheromone[0][1]<<endl;
     sleep(3);

   //cout << "node length from other " <<(*dist)[0][1]<<endl;
   //printf("Process %d is executing aco algorithm with colony size: %d\n", world_rank, colony_size);

   struct ant_t ants[colony_size];

   for (int p = 0; p< sizeof pheromone /sizeof pheromone[0]; p++){
	   //for (int pc = 0; pc <sizeof pheromone[0] / sizeof pheromone[0][0]; pc++){
	  // cout<<setprecision(3)<< pheromone[p][pc]<<"  " ;  
	   //}
       	   cout <<endl;
    }


   for (int ant = 0; ant < colony_size; ant++){
      cout << "\nAnt "<< ant << " walking"<<endl;
      while (antColony[ant].pathIndex < MAX_NODES){
      //cout << "Lets select the next city " << endl;
        antColony[ant].nextNode = selectNextCity(ant);
	//cout << antColony[ant].tabu[antColony[ant].nextNode]<<" ";
        antColony[ant].tabu[antColony[ant].nextNode]=1;
	
        //cout << "current -> " << antColony[ant].curNode << endl;
        //cout << "next -> " << antColony[ant].nextNode << endl;
        antColony[ant].path[antColony[ant].pathIndex] = antColony[ant].nextNode;

	//cout << "City selected " << endl;
	cout << antColony[ant].nextNode << " ";
	//cout << antColony[ant].path[antColony[ant].pathIndex]<<endl;

	antColony[ant].pathIndex = antColony[ant].pathIndex+1;
        antColony[ant].tourLength += (*dist)[antColony[ant].curNode][antColony[ant].nextNode];   
        //Now the current node is the node selected (next node)
        //printf("Tour Length %f\n", antColony[ant].tourLength);
        antColony[ant].curNode = antColony[ant].nextNode;
      }
     antColony[ant].tourLength += (*dist)[antColony[ant].path[MAX_NODES-1]][antColony[ant].path[0]]; 
     ants[ant].tourLength = antColony[ant].tourLength;
     //printf("Se ha enviado %f\n", ants[ant].tourLength);
     ants[ant] = antColony[ant];
   }

  MPI_Send(&ants, colony_size, ant_type, 0,1, MPI_COMM_WORLD);
  }
  

  else if (world_rank == 0) {
    	  
   ant_t container_buffer_r[world_size][colony_size];	  
   ant_t *container_buff[world_size];
   for (int i=1; i<world_size; i++){
   // MPI_Recv(container[i-1], colony_size, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&container_buffer_r[i-1], colony_size, ant_type, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
   }
    for (int colony=0; colony<world_size-1; colony++){
      updatePheromone(container_buffer_r[colony]);
      int bestSol, bestant;
      bestSol, bestant = bestSolution(container_buffer_r[colony]);
      //for (int j1=0; j1<MAX_NODES; j1++){
       // cout << container_buffer_r[colony][bestant].path[j1]<<" "; 
      }
    //}

    //cout <<" Phermone is : " << pheromone[0][2]<<endl;
 

   //float b = bestSolution(container_buff, world_size);
   int contsize = sizeof container_buffer_r / sizeof container_buffer_r[0] ;
   //cout << "Container buffer rows: " << contsize <<endl;
   int colonys = sizeof container_buffer_r[0] / sizeof (ant_t);
   //cout << "Container buffer ants: " << colonys <<endl;

   //for (int j=1; j<world_size; j++){
    
   // printf("Se recibieron estos datos del process: %d: \n",j);
    //printf("Tour Length of the first ant: %f\n", container_buffer_r[j-1][0].tourLength);
    //printf("Tour Length of the last ant: %f\n", container_buffer_r[j-1][colony_size-1].tourLength);
   // for (int k =0; k<MAX_NODES; k++){
    //  cout << container_buffer_r[j-1][0].path[k] << " ";
    //}
   }

  	
  cout << "Generation best tour is : " << best << endl;
  //}

  //lets broadcast the pheromone values;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&pheromone[0][0], MAX_NODES*MAX_NODES, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank!=0){
  cout << "Pheromone values at the iteration end: " << pheromone[0][1] << endl;
  restartAnts();
  }
  }
 }


int main(int argc, char* argv[]){
  cout << "--------------------------------------"<<endl;
  cout <<" PARLLEL ACO WITH MPI"<<endl;
  cout << "--------------------------------------"<<endl;

  if (argc < 3){
    cerr << "Usage: " <<argv[0] << " instance_name " << "number of ants" <<endl;
  }
  else{ 
  graph = argv[1];
  //MAX_ANTS = atoi(argv[2]);
  dist = new EdgeMatrix();
  constructTSP(graph, cities, dist);

  MPI_Init(NULL, NULL);
   //for (int i =0; i<6; i++){
	   
    solveTSP();
    
  //}

  MPI_Finalize();

 }

}
