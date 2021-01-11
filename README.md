## To execute the parallel program use: 
`mpicc -lstdc++ -lm -lpthread parallel_acop.cpp`

## Then run 

`mpirun -n 4 ./a.out wi29 100`

If you want to run with a different instance your first need to configure the ants.h file with the max nodes = to the instance size
