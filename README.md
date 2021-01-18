## To use mpicc and mpirun we have use conda to get mpi library
`conda activate acoenv`

## To compile the parallel program use: 
`mpicc -lstdc++ -lm -lpthread parallel_acop.cpp`

## Then run 

`mpirun -n 4 ./a.out rl1889`

If you want to run with a different instance your first need to configure the ants.h file with the max nodes = to the instance size
According to the instance coordinates, if you have int then you need to configure the int x and int y in the ants.h and in the parallel_aco construct tsp function
