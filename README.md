## To use mpicc and mpirun we have use conda to get mpi library
`conda activate acoenv`

## To compile the parallel program use: 
`mpicc -lstdc++ -lm -lpthread parallel_acop.cpp`

## Then run 

`mpirun -n 4 ./a.out rl1889`

To select a different tsp instance, modify the MAX_NODES number according to the instance size in ants.h.
The algorithm parameters (alpha, beta, rho, qval, iterations) can be modified in ants.h.
Depending on the coordinate format of the instance you can configure the type of coordinate (x,y) in cityType in ants.h and in parllel_aco.cpp as int or double. 

