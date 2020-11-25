// Traveling Antsmen constants and struct declarations

#define MAX_NODES 280 // number of vertices
#define MAX_DIST 100000
#define MAX_TOUR (MAX_NODES * MAX_DIST)
//#define MAX_ANTS MAX_NODES
#define MAX_ANTS 100 
#define NUM_EDGES ((MAX_NODES * MAX_NODES - MAX_NODES) / 2)

struct cityType {
  int x, y;
};

// allows edge queries as a 2D array
class EdgeMatrix {
  float *dist;
public:
  EdgeMatrix() {
    dist = new float[MAX_NODES * MAX_NODES];
  }
  ~EdgeMatrix() {
    delete dist;
  }
  float* operator[](unsigned int i) {
   return &dist[MAX_NODES * i];
 }

  float *get_array(){
    return dist;
  }
};

//Ant algorithm parameters
#define ALPHA 7.0 // alpha controls the pheromone influence
#define BETA 1.0 // this parameter raises the weight of distance over pheromone
#define RHO 0.1 // evaporation rate
#define QVAL 1  // amount of pheromone that will be incremented in the update trails processs
#define MAX_ITERATIONS 10 
#define MAX_TIME (MAX_ITERATIONS * MAX_NODES)
#define INIT_PHER (1.0 / MAX_NODES)
