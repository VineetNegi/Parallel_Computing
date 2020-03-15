
/**********************************************************************
                  Simulation Specifications (constants)
By: VINEET NEGI

Description -
This is header file which contains all the simulation specifications as
constants

History -
created on 4/5/2019

***********************************************************************/
// DEFINE MACROS AND CONSTANTS

#define NODAL_Z 4

#define NUM_NODAL_DATA 5

#define DX 0.004
#define DY 0.004
#define DT 1e-6
#define NUM_ITER 10000

#define GRAPH_COMM 0

#define NUM_THREADS 1

/**********************************************************************/

// #define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#define BGQ 0 // BGQ = 1 when running on BG/Q
#if BGQ == 1
    #include<hwi/include/bqc/A2_inlines.h>
#else
    #define GetTimeBase MPI_Wtime
#endif


/***********************************************************************/
