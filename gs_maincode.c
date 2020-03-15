#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <pthread.h>
/**********************************************************************/
// Include primitive dtypes
#include "gs_pdtypes.h"

// Include constants
#include "gs_constants.h"

// Include primitive dtypes
#include "gs_sdtypes.h"

// Include utilities functions
#include "gs_utils.h"

/**********************************************************************
                         PARALLEL GRIDSOLVER - 2D
By: VINEET NEGI

Description -
This code performs explicit time marching on a 2D grid. Each nodal points
on the grid has three properties -
1) disp - a array of dof at a node
2) forces - an array of forces at a node
3) field - an array of fields defined at a node

Only disps are updated in solver. Forces and fields may be updated but their
update will have to be defined explicitly

* USAGE :
* mpirun -np #tasks executable  g_parts g_xadj g_adjncy g_global_vert_ids g_nodal_data

History -
02/27/2019: Serial version of the code finished. Some touchup is still required
***********************************************************************/


/**********************************************************************/
// DEFINE GLOBAL DATA STRUCTURES

// total number of domain nodes
int_t num_nodes;

// total number of ghost nodes
int_t num_ghost_nodes;

// array of all (local domain) nodes
struct node *nodes;

// array of ghost nodes (nodes lying on adjacent processes)
struct node *ghost_nodes;

// define physics
struct physics phy;

// define solution
struct solution sol;

// define communicating processes called as adj_processes
struct process *adj_processes;

// number of adj processes
int_t num_adj_processes;

// node offset
int_t node_offset;

//define main input data structures
/*Input GLOBAL datastructure ... start with 'g' prefix ... 
 * This DS is only required at rank 0 */
int_t nvtxs; // number of global vertices
int_t *g_parts; //list of partitions
int_t *g_xadj; //list of node adjacency pointer
int_t *g_adjncy; //node adjacencies
int_t *g_global_vert_ids; //a list of global ids of the vertices
real_t **g_nodal_data; // a double 2D array of nodal data

/*Input LOCAL datastructures ... starts with 'l' prefix ...
 * They are required at al ranks
 */
int_t l_nvtxs; // number of vertices on a partition
int_t *l_xadj; // list of node adjacency pointer
int_t *l_adjncy; // list of node adjacencies 
int_t *l_global_vert_ids; //a list of global ids of the vertices
real_t **l_nodal_data; //a double 2D array of nodal data
 
int_t *g_n2r_map; // a map of global node id to process or rank


/* ptheards related global vars*/
pthread_t p_threads[NUM_THREADS - 1]; 
struct thread_args t_args[NUM_THREADS];
pthread_attr_t attr; 
pthread_barrier_t barr; //barrier variable

/**********************************************************************/
//inlude key functions for the -> grid/physics/solution <-
#include "gs_keyFuncs.h"

/**********************************************************************/
// MPI RELATED GLOBAL FUNCTION

void get_graph_comm(MPI_Comm old_comm, MPI_Comm *new_comm)
{
	/* This function changes the old communicator (Default)
	 * to a new optimized one (graph) based
	 */
	 int_t i;
	 int_t *adj_ranks = calloc(num_adj_processes, sizeof(*adj_ranks)); 
	 int_t *adj_weights = calloc(num_adj_processes, sizeof(*adj_weights));
	 
	 // get adj_ranks and adj_weights from adj_processes list
	 for(i = 0; i < num_adj_processes; i++)
	 {
		 adj_ranks[i] = adj_processes[i].rank;
		 adj_weights[i] = adj_processes[i].num_recv_nodes; //num_recv_nodes = num_send_nodes (symmetric communication)
	 }
	 
	 //create graph topology
	 if(MPI_Dist_graph_create_adjacent(old_comm, num_adj_processes,
							adj_ranks, adj_weights, num_adj_processes,
							adj_ranks, adj_weights, MPI_INFO_NULL, 0, 
							new_comm) != MPI_SUCCESS)
		{
			printf("ERROR: graph communicator creation failed\n");
		}
	
	//free function vars on heap
	free(adj_ranks);
	free(adj_weights);
	 
	 return;
 }

void request_ghost_nodes(MPI_Request *recv_requests, MPI_Comm comm)
{
    /*
    This function goes over all the processes and requests the data
    from adjacent processes
    */

    int_t i;
    for(i = 0; i < num_adj_processes; i++)
    {
        MPI_Irecv(adj_processes[i].recv_buf, adj_processes[i].num_recv_nodes,
        mpi_real_t, adj_processes[i].rank, 0, comm, &recv_requests[i]); // post the receive request for a rank
    }

    return;
}

void send_ghosted_nodes(MPI_Request *send_requests, MPI_Comm comm)
{
    /*
    This function goes over all the processes and sends the data to
    adjacent processes
    */

    int_t i;
    for(i = 0; i < num_adj_processes; i++)
    {
        MPI_Isend(adj_processes[i].send_buf, adj_processes[i].num_send_nodes,
        mpi_real_t, adj_processes[i].rank, 0, comm, &send_requests[i]); // send the data to adj processes
    }

    return;
}

void update_ghost_nodes(int_t it, MPI_Request *recv_requests, MPI_Comm comm)
{
    /*
    This function updates the ghost nodes once the recv communication is finished
    */
    int_t i, j, idx;
    idx = it%2; //determine if the current iteration is odd or even
    
    for(i = 0; i <  num_adj_processes; i++)
    {
        MPI_Wait(&recv_requests[i], MPI_STATUS_IGNORE); // wait for the Irecv to finish

        for(j = 0; j < adj_processes[i].num_recv_nodes; j++)
        {
			adj_processes[i].recv_locs[j]->iter = it; //update the iteration counter (not required)
            adj_processes[i].recv_locs[j]->disp[idx] = adj_processes[i].recv_buf[j]; // copy from recv buffer
        }
    }

    return;

}

void copy_ghosted_nodes(int_t it)
{
    /*
    This function copies the ghosted nodes to the send_buf
    */
    int_t i, j, idx;
    idx = it%2; //determine if the current iteration is odd or even
    
    for(i = 0; i <  num_adj_processes; i++)
    {
        for(j = 0; j < adj_processes[i].num_send_nodes; j++)
        {
             adj_processes[i].send_buf[j] = adj_processes[i].send_locs[j]->disp[idx]; // copy to send buffer
        }
    }

    return;
}

/**********************************************************************/

/**********************************************************************/
// DEFINE SOLUTION RELATED GLOBAL FUNCTIONS

/****************** functions for data input at rank 0 ****************/

void read_data(char **argv, int_t npes)
{
	/*
	 * This function reads the main input files.
	 * This function is called at rank 0
	 */
	 	 
	 int_t i;
	 FILE *input_file;
	 
	 //update
     printf("Reading input files ...\n");
         
     //read g_parts
     g_parts = calloc(npes + 1, sizeof(*g_parts));
     input_file = fopen(argv[1], "rb");
     if(fread(g_parts, sizeof(*g_parts), npes + 1, input_file) != npes + 1)
     {
		 printf("ERROR: g_parts file could not be read\n");
		 exit(-1);
	 }
	 fclose(input_file);

	 nvtxs = g_parts[npes];
	 
	 //read g_xadj
     g_xadj = calloc(nvtxs + 1, sizeof(*g_xadj));
     input_file = fopen(argv[2], "rb");
     if(fread(g_xadj, sizeof(*g_xadj), nvtxs + 1, input_file) != nvtxs + 1)
     {
		 printf("ERROR: g_xadj file could not be read\n");
		 exit(-1);
	 }
	 fclose(input_file);
	 	 
	 //read g_adjncy
     g_adjncy = calloc(g_xadj[nvtxs], sizeof(*g_adjncy));
     input_file = fopen(argv[3], "rb");
     if(fread(g_adjncy, sizeof(*g_adjncy), g_xadj[nvtxs], input_file) != g_xadj[nvtxs])
     {
		 printf("ERROR: g_adjncy file could not be read\n");
		 exit(-1);
	 }
	 fclose(input_file);
	 	 
	 //read g_global_vert_ids
     g_global_vert_ids = calloc(nvtxs, sizeof(*g_global_vert_ids));
     input_file = fopen(argv[4], "rb");
     if(fread(g_global_vert_ids, sizeof(*g_global_vert_ids), nvtxs, input_file) != nvtxs)
     {
		 printf("ERROR: g_global_vert_ids file could not be read\n");
		 exit(-1);
	 }
	 fclose(input_file);	 
	 
	 //read g_nodal_data (try to make the read cache friendly)
	 
	 g_nodal_data = calloc(NUM_NODAL_DATA + 2, sizeof(*g_nodal_data));
	 
	 for(i=0; i<NUM_NODAL_DATA + 2; i++)
	 {
		 g_nodal_data[i] = calloc(nvtxs, sizeof(*g_nodal_data[i]));
	 }
     input_file = fopen(argv[5], "r");
     
	 for(i=0; i < 2+NUM_NODAL_DATA; i++)
	 {		
		if(fread(g_nodal_data[i], sizeof(*g_nodal_data[i]), nvtxs, input_file) != nvtxs)
		{
			printf("ERROR: g_nodal_data file could not be read\n");
			exit(-1);
		}
	 }
	 fclose(input_file);	 
	 		 
	 return;
 }

void set_n2r_map(int_t npes, int_t *g_partition, int_t *g_nodes, int_t *mapping)
{
	/*
	 * This function creates a mapping from nodes to the rank/process/-
	 * partions where they reside
	 */
	 
	 int_t i, j, n;

	 for(i=1; i < npes + 1; i++)
	 {
		 for(j = g_partition[i - 1]; j < g_partition[i]; j++)
		 {
			 n = g_nodes[j]; // get the global node id
			 mapping[n] = i - 1; // set the rank to the mapping at node loc
		 }
	 }
	 
	 return;
 }
	
/****** functions for data distribution during simulation setup ******/
	 
void receive_data(int_t rank, MPI_Comm comm)
{
	/* This function is called at all ranks except 0
	 * This function receives the data from rank 0
	 */
	 int_t i;
	 MPI_Request recv_request;
	 MPI_Request recv_requests[NUM_NODAL_DATA + 2]; //to receive nodal_data
	 
	 //recv node_offset
	 MPI_Irecv(&node_offset, 1, mpi_int_t, 0, 99 , comm, &recv_request); 
	 MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
	 
	 //recv l_nvtxs
	 MPI_Irecv(&l_nvtxs, 1, mpi_int_t, 0, 0 , comm, &recv_request); 
	 MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
	 
	 //recv l_xadj
	 l_xadj = calloc(l_nvtxs + 1, sizeof(*l_xadj));
	 MPI_Irecv(l_xadj, l_nvtxs + 1, mpi_int_t, 0, 1, comm, &recv_request); 
	 MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
	
	 //recv l_adjncy
	 l_adjncy = calloc(l_xadj[l_nvtxs] - l_xadj[0], sizeof(*l_adjncy));
	 MPI_Irecv(l_adjncy, l_xadj[l_nvtxs] - l_xadj[0], mpi_int_t, 0, 2, comm, &recv_request); 
	 MPI_Wait(&recv_request, MPI_STATUS_IGNORE); 
	 
	 //recv l_global_vert_ids
	 l_global_vert_ids = calloc(l_nvtxs, sizeof(*l_global_vert_ids));
	 MPI_Irecv(l_global_vert_ids, l_nvtxs, mpi_int_t, 0, 3, comm, &recv_request); 
	 MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
	 
	 //recv g_n2r_map
	 MPI_Bcast(&nvtxs, 1, mpi_int_t, 0, comm);	 
	 g_n2r_map = calloc(nvtxs, sizeof(*g_n2r_map));
	 MPI_Bcast(g_n2r_map, nvtxs,mpi_int_t, 0, comm);
	 
	 //recv l_nodal_data
	 l_nodal_data = calloc(NUM_NODAL_DATA + 2, sizeof(*l_nodal_data));
	 for(i=0; i<NUM_NODAL_DATA + 2; i++)
	 {
		 l_nodal_data[i] = calloc(l_nvtxs, sizeof(*l_nodal_data[i]));
		 MPI_Irecv(l_nodal_data[i], l_nvtxs, mpi_real_t, 0, 4 + i, comm, &recv_requests[i]); 		 
	 }
	 MPI_Waitall(NUM_NODAL_DATA + 2, recv_requests, MPI_STATUSES_IGNORE);
	 	 
	 return;
 }
	 
void send_data(int_t npes, MPI_Comm comm)
{
	/* This function sends the data to all the process from rank 0
	 */
	 
	 // At rank 0 send data internally
	 
	 int_t i, j, k, temp_int, temp_int2, *temp_arr;
	 MPI_Request *send_requests;
	 
	 temp_arr = calloc(npes, sizeof(*temp_arr));
	 send_requests = calloc(npes - 1, sizeof(*send_requests));
	 
	 //update
	 printf("Starting data distribution ...\n");
	 
	 //send node_offset (sum of nodes in processes before this)
	 for(i = 1; i < npes + 1; i++)
	 {
		 temp_int = i - 1;
		 temp_arr[temp_int] = g_parts[i-1]; //buff size
		 
		 if(temp_int > 0 )
		 {
			 MPI_Isend(&temp_arr[temp_int], 1, mpi_int_t, temp_int, 99, comm, &send_requests[temp_int - 1]); //send to other ranks
		 }
		 else
		 {
			 node_offset = temp_arr[temp_int]; // assign internally
		 }
		 
	 }
	 MPI_Waitall(npes - 1, send_requests, MPI_STATUSES_IGNORE);

	 
	 //send the l_nvtxs
	 for(i = 1; i < npes + 1; i++)
	 {
		 temp_int = i - 1;
		 temp_arr[temp_int] = g_parts[i] - g_parts[i - 1]; //buff size
		 
		 if(temp_int > 0 )
		 {
			 MPI_Isend(&temp_arr[temp_int], 1, mpi_int_t, temp_int, 0, comm, &send_requests[temp_int - 1]); //send to other ranks
		 }
		 else
		 {
			 l_nvtxs = temp_arr[temp_int]; // assign internally
		 }
		 
	 }
	 MPI_Waitall(npes - 1, send_requests, MPI_STATUSES_IGNORE);
	 
	 // send the l_xadj
	 for(i = 1; i < npes + 1; i++)
	 {
		 temp_int = i - 1; // rank of the target process
		 temp_arr[temp_int] = g_parts[i] - g_parts[i - 1] + 1; //buff size
		 
		 if(temp_int > 0 )
		 {
			 MPI_Isend(&g_xadj[g_parts[temp_int]], temp_arr[temp_int], mpi_int_t, temp_int, 1, comm, &send_requests[temp_int - 1]); //send to other ranks
		 }
		 else
		 {
			 // assign internally
			 l_xadj = calloc(l_nvtxs + 1, sizeof(*l_xadj));
			 for(j = g_parts[i - 1]; j <= g_parts[i]; j++)
			 {
				 l_xadj[j - g_parts[i - 1]] = g_xadj[j];
			 }
		 }
	 }
	 MPI_Waitall(npes - 1, send_requests, MPI_STATUSES_IGNORE);
	 
	 //send the l_adjncy
	 for(i = 1; i < npes + 1; i++)
	 {
		 temp_int = i - 1; // rank of the target process
		 temp_int2 = g_xadj[g_parts[temp_int]]; //starting point of adjncy list
		 temp_arr[temp_int] = g_xadj[g_parts[temp_int + 1]]  - g_xadj[g_parts[temp_int]]; //buff size
		 
		 if(temp_int > 0 )
		 {
			 MPI_Isend(&g_adjncy[temp_int2], temp_arr[temp_int], mpi_int_t, temp_int, 2, comm, &send_requests[temp_int - 1]); //send to other ranks
		 }
		 else
		 {
			 // assign internally
			 l_adjncy = calloc(l_xadj[l_nvtxs] - l_xadj[0], sizeof(*l_adjncy));
			 for(j = l_xadj[0]; j < l_xadj[l_nvtxs]; j++)
			 {
				 l_adjncy[j - l_xadj[0]] = g_adjncy[j];
			 }
		 }
	 }
	 MPI_Waitall(npes - 1, send_requests, MPI_STATUSES_IGNORE);
	 
	 // send the l_global_vert_ids
	 for(i = 1; i < npes + 1; i++)
	 {
		 temp_int = i - 1; // rank of the target process
		 temp_arr[temp_int] = g_parts[i] - g_parts[i - 1]; //buff size
		 
		 if(temp_int > 0 )
		 {
			 MPI_Isend(&g_global_vert_ids[g_parts[temp_int]], temp_arr[temp_int], mpi_int_t, temp_int, 3, comm, &send_requests[temp_int - 1]); //send to other ranks
		 }
		 else
		 {
			 // assign internally
			 l_global_vert_ids = calloc(l_nvtxs, sizeof(*l_global_vert_ids));
			 for(j = g_parts[i - 1]; j < g_parts[i]; j++)
			 {
				 l_global_vert_ids[j - g_parts[i - 1]] = g_global_vert_ids[j];
			 }
		 }
	 }
	 MPI_Waitall(npes - 1, send_requests, MPI_STATUSES_IGNORE);	 
	 
	 // send the g_n2r_map (broadcast)
	 MPI_Bcast(&nvtxs, 1, mpi_int_t, 0, comm); 
	 MPI_Bcast(g_n2r_map, nvtxs, mpi_int_t, 0, comm);
	 
	 //send the nodal_data
	 l_nodal_data = calloc(NUM_NODAL_DATA + 2, sizeof(*l_nodal_data));
	 for(i=0; i<NUM_NODAL_DATA + 2; i++)
	 {
		 l_nodal_data[i] = calloc(l_nvtxs, sizeof(*l_nodal_data[i]));
	 }
	 
	 for(k = 0; k < NUM_NODAL_DATA + 2; k++)
	 {
		 // send data for each field
		 for(i = 1; i < npes + 1; i++)
		 {
			 temp_int = i - 1; // rank of the target process	 
			 temp_arr[temp_int] = g_parts[i] - g_parts[i - 1]; //buff size
			 if(temp_int > 0 )
			 {
				 MPI_Isend(&g_nodal_data[k][g_parts[temp_int]], temp_arr[temp_int], mpi_real_t, temp_int, 4 + k, comm, &send_requests[temp_int - 1]); //send to other ranks
			 }
			 else
			 {
				 // assign internally
				 for(j = g_parts[i - 1]; j < g_parts[i]; j++)
				 {
					 l_nodal_data[k][j - g_parts[i - 1]] = g_nodal_data[k][j];
				 }
			 }
		 }	
		 MPI_Waitall(npes - 1, send_requests, MPI_STATUSES_IGNORE);
	 }		

	 //update
	 printf("Finished data distribution ...\n"); 
	 
	 //free function local vars on heap
	 free(temp_arr);
	 free(send_requests);
	 
	 //free global mesh data on heap
	 free(g_xadj);
	 free(g_adjncy);
	 free(g_global_vert_ids);
	 for(i=0; i<NUM_NODAL_DATA + 2; i++)
	 {
		 free(g_nodal_data[i]);
	 }
	 free(g_nodal_data);
	 free(g_parts);
	 
			 
	 return;	 
 }

/**************** function for initial setup of system ****************/

void grid_initialize(int_t rank)
{
	/* This function  sets up the grid at each rank based on the data
	 * recevied in MPI exchange
	 */
	 
	 int_t i, j;
	 int_t idx, j0, j1, n_idx, counter = 0;
	 int_t temp0;
	 int_t flag;
	 
	 //set the number of domain nodes
	 num_nodes = l_nvtxs;
	 
	 // Allocate memory for nodes
	 nodes = calloc(num_nodes, sizeof(*nodes));
	
	 //loop through the global nodes vertices and create the domain nodes
	 for(i = 0; i < num_nodes; i++)
	 {	
		 nodes[i].id = l_global_vert_ids[i];
		 nodes[i].nodetype = (int_t) l_nodal_data[2][i]; //nodetype is a nodal data
		 nodes[i].rank_status = 0;
		 nodes[i].iter = 0;
		 for(j = 0; j < 2; j++)
		 {
			 nodes[i].loc[j] = l_nodal_data[j][i];
			 nodes[i].disp[j] = 0.0;
			 nodes[i].forces[j] = 0.0;
		 }
		 for(j = 0; j < NODAL_Z; j++)
		 {
			 nodes[i].neighbors[j] = NULL;
		 }
		 for(j=0; j<4; j++)
		 {
			 nodes[i].field[j] = 0.0;
		 }
		 
		 //set values of nodal data except loc and nodetype
		 for(j = 0; j < 2; j++)
		 {
			 nodes[i].disp[j] = l_nodal_data[3][i]; //set initial disp value
			 nodes[i].forces[j] = l_nodal_data[4][i]; //set forces
			 nodes[i].field[2*j] = l_nodal_data[5][i]; //set field - x dir
			 nodes[i].field[2*j + 1] = l_nodal_data[6][i]; //set field - y dir
		 }
	 }

	 //loop through all the adjacencies and fix them on nodes
	 
	 
	 for(i = 0; i < num_nodes; i++)
	 {
		 // get the starting and end index of the adjacencies in l_adjncy (subract the starting node index)
		 j0 = l_xadj[i] - l_xadj[0];
		 j1 = l_xadj[i + 1] - l_xadj[0];
		 
		 
		 for(j = j0; j < j1; j++)
		 {
			 //get the id of the adj node
			 n_idx = l_adjncy[j];

			 if(g_n2r_map[n_idx] == rank)
			 {
				 //node is in the domain
				 
				 //get the index of the adjacent node
				 idx = binarySearch(l_global_vert_ids, 0, num_nodes - 1, n_idx);
				 
				 //attach the node pointer at suitable location
				 nodes[i].neighbors[get_rel_position_tag(nodes[i].id,nodes[idx].id)] = &nodes[idx];
				  
			 }
			 else
			 {
				 // node is outside of the domain
				 counter++;
			 } 
		 }
	 }
	 
	 //create local arrays for ghosted, ghost nodes ids and shared processors ids
	 int_t b_size = counter;
	 int_t *N_i, *N_j, *P_j, *n_j, *p_j;
	 N_i = calloc(b_size, sizeof(*N_i));
	 N_j = calloc(b_size, sizeof(*N_j));
	 P_j = calloc(b_size, sizeof(*P_j));
	 n_j = calloc(b_size, sizeof(*n_j));
	 p_j = calloc(b_size, sizeof(*p_j));
	 
	 counter = 0;
	 for(i = 0; i < l_nvtxs; i++)
	 {
		 // get the starting and end index of the adjacencies in l_adjncy (subract the starting node index)
		 j0 = l_xadj[i] - l_xadj[0];
		 j1 = l_xadj[i + 1] - l_xadj[0];
		 
		 
		 for(j = j0; j < j1; j++)
		 {
			 //get the id of the adj node
			 n_idx = l_adjncy[j];

			 if(g_n2r_map[n_idx] != rank)
			 {
				 // node is outside of the domain
				 N_i[counter] = l_global_vert_ids[i];
				 N_j[counter] = n_idx;
				 n_j[counter] = n_idx;
				 P_j[counter] = g_n2r_map[n_idx];
				 p_j[counter] = P_j[counter];
				 counter++;
			 }
		 }
	 }
	 
	 //sort and get unique no. of adj nodes
	 quickSort(n_j, 0, counter - 1); 
	 j0 = unique(n_j, counter); //this is total number of unique ghost nodes
	 num_ghost_nodes = j0;
	 
	 //alocate memory for ghost nodes and create them
	 ghost_nodes = calloc(num_ghost_nodes, sizeof(*ghost_nodes));
	 for(i=0; i<num_ghost_nodes; i++)
	 {
		 ghost_nodes[i].id = n_j[i];
		 ghost_nodes[i].nodetype = 0;
		 ghost_nodes[i].rank_status = 2;
		 ghost_nodes[i].iter = 0;
		 for(j = 0; j < NODAL_Z; j++)
		 {
			 ghost_nodes[i].neighbors[j] = NULL;
		 }
		 for(j = 0; j < 2; j++)
		 {
			 ghost_nodes[i].loc[j] = 0.0;
			 ghost_nodes[i].disp[j] = 0.0;
			 ghost_nodes[i].forces[j] = 0.0;
		 }
		 for(j = 0; j < 4; j++)
		 {
			 ghost_nodes[i].field[j] = 0.0;
		 }
	 }
	 
	 //sort and get unique no. of adj ranks
	 quickSort(p_j, 0, counter - 1);
	 j1 = unique(p_j, counter);
	 num_adj_processes = j1;
	 
	 //create adj process and their feed counter
	 int_t *proc_feed_counter = calloc(j1, sizeof(*proc_feed_counter));
	 adj_processes = calloc(num_adj_processes, sizeof(*adj_processes));
	 for(i=0; i<num_adj_processes;i++)
	 {
		 adj_processes[i].rank = p_j[i];
	 }
	 
	 //fix the ghost nodes
	 
	 //count the num_recv_nodes
	 for(i=0; i < num_ghost_nodes; i++)
	 {
		 j = binarySearch(p_j, 0, num_adj_processes - 1, g_n2r_map[n_j[i]]);
		 if(j > -1)
		 {
			 proc_feed_counter[j] += 1;
		 }
		 else
		 {
			 printf("ERROR: process not found\n");
		 }
	 }
	 
	 //allocate the memory within processes
	 for(i=0; i<num_adj_processes; i++)
	 {
		 adj_processes[i].num_recv_nodes = proc_feed_counter[i];
		 adj_processes[i].recv_buf = calloc(proc_feed_counter[i], sizeof(*adj_processes[i].recv_buf));
		 adj_processes[i].recv_locs = calloc(proc_feed_counter[i], sizeof(*adj_processes[i].recv_locs));
		 proc_feed_counter[i] = 0;
	 }
	 
	 //add the ghost node addresses to the processes
	 for(i=0; i < num_ghost_nodes; i++)
	 {
		 j = binarySearch(p_j, 0, num_adj_processes - 1, g_n2r_map[n_j[i]]);
		 if(j > -1)
		 {
			 adj_processes[j].recv_locs[proc_feed_counter[j]] = &ghost_nodes[i];
			 proc_feed_counter[j] += 1;
		 }
		 else
		 {
			 printf("ERROR: process not found\n");
		 }
	 }
	 
	 //fix the ghosted nodes
	 /*****************************************************************/
	 for(i=0; i<num_adj_processes; i++)
	 {
		 proc_feed_counter[i] = 0;
	 }
	 
	 for(i=0; i < b_size; i++)
	 {
		 idx = binarySearch(l_global_vert_ids, 0, num_nodes - 1, N_i[i]); //get the idx of n_i node (in the domain)
		 temp0 = binarySearch(n_j, 0, j0 - 1, N_j[i]); //get the idx of n_i node (in the domain)
		 
		 //modify the status of the ghosted node
		 nodes[idx].rank_status = 1;
		 
		 //attach the node pointer at suitable location
		 nodes[idx].neighbors[get_rel_position_tag(nodes[idx].id,ghost_nodes[temp0].id)] = &ghost_nodes[temp0];
		 
		 flag = 0;
		 if( i < b_size - 1)
		 {
			 j = 1;
			 while(N_i[i + j] == N_i[i])
			 {
				 if(P_j[i + j] == P_j[i])
				 {
					 flag = 1;
					 break;
				 }
				 j++;
				 
				 if((i+j) >= b_size)
				 {
					 break;
				 }
			 }
		 }
		 
		 if(flag == 0)
		 {
			temp0 = binarySearch(p_j, 0, num_adj_processes - 1, P_j[i]); // find the id of the P_j process in p_j (unique)
			proc_feed_counter[temp0] += 1;

		 }
	 }
	 
	 //allocate the memory within processes
	 
	 for(i=0; i < num_adj_processes; i++)
	 {
		 adj_processes[i].num_send_nodes = proc_feed_counter[i];
		 adj_processes[i].send_buf = calloc(proc_feed_counter[i], sizeof(*adj_processes[i].send_buf));
		 adj_processes[i].send_locs = calloc(proc_feed_counter[i], sizeof(*adj_processes[i].send_locs));
		 proc_feed_counter[i] = 0;
	 }
	 
	 //add the ghosted nodes to the adj_processes
	 
	 for(i=0; i < b_size; i++)
	 {
		 idx = binarySearch(l_global_vert_ids, 0, num_nodes - 1, N_i[i]); //get the idx of n_i node (in the domain)
		 		 
		 flag = 0;
		 if( i < b_size - 1)
		 {
			 j = 1;
			 while(N_i[i + j] == N_i[i])
			 {
				 if(P_j[i + j] == P_j[i])
				 {
					 flag = 1;
					 break;
				 }
				 
				 j++;
				 
				 if((i+j) >= b_size)
				 {
					 break;
				 }
			 }
		 }
		 
		 if(flag == 0)
		 {
			 temp0 = binarySearch(p_j, 0, num_adj_processes - 1, P_j[i]); // find the id of the P_j process in p_j (unique)
			 
			 if(temp0 > -1)
			 {
				 adj_processes[temp0].send_locs[proc_feed_counter[temp0]] = &nodes[idx];
				 proc_feed_counter[temp0] += 1;
			 }
			 else
			 {
				 printf("ERROR: process not found\n");
			 }
		 }
	 }
	 		 
	 /*****************************************************************/
	 
	//free function local vars on heap
	free(N_i);
	free(N_j);
	free(P_j);
	free(n_j);
	free(p_j);
	free(proc_feed_counter);
	
	//free local mesh data on heap
	free(l_xadj);
	free(l_adjncy);
	free(l_global_vert_ids);
	free(g_n2r_map);
	for(i=0; i<NUM_NODAL_DATA + 2; i++)
	{
		free(l_nodal_data[i]);
	}
	free(l_nodal_data);
		
	 return;
 }

void init_ghost_nodes(int rank, MPI_Comm comm)
{
	/* This function initializes the ghost nodes at the beginning of the 
	 * simulation
	 */
	 
	 int i, j, r;
	 MPI_Request *recv_requests, *send_requests;
	 
	 recv_requests = calloc(num_adj_processes, sizeof(*recv_requests));
	 send_requests = calloc(num_adj_processes, sizeof(*send_requests));
	 
	 //exchange node_loc
	 for(r = 0; r < 2 ; r++)
	 {
		 request_ghost_nodes(recv_requests, comm);
		 
		 for(i = 0; i <  num_adj_processes; i++)
		{
			for(j = 0; j < adj_processes[i].num_send_nodes; j++)
			{
				 adj_processes[i].send_buf[j] = adj_processes[i].send_locs[j]->loc[r]; // copy to send buffer
			}
		}
		
		 send_ghosted_nodes(send_requests, comm);
		 
		for(i = 0; i <  num_adj_processes; i++)
		{
			MPI_Wait(&recv_requests[i], MPI_STATUS_IGNORE); // wait for the Irecv to finish

			for(j = 0; j < adj_processes[i].num_recv_nodes; j++)
			{
				adj_processes[i].recv_locs[j]->loc[r] = adj_processes[i].recv_buf[j]; // copy from recv buffer
			}
		}
		
		MPI_Waitall(num_adj_processes, send_requests, MPI_STATUSES_IGNORE);
	}
	
	//exchange node_disp
	 for(r = 0; r < 2 ; r++)
	 {
		 request_ghost_nodes(recv_requests, comm);
		 
		 for(i = 0; i <  num_adj_processes; i++)
		{
			for(j = 0; j < adj_processes[i].num_send_nodes; j++)
			{
				 adj_processes[i].send_buf[j] = adj_processes[i].send_locs[j]->disp[r]; // copy to send buffer
			}
		}
		
		 send_ghosted_nodes(send_requests, comm);
		 
		for(i = 0; i <  num_adj_processes; i++)
		{
			MPI_Wait(&recv_requests[i], MPI_STATUS_IGNORE); // wait for the Irecv to finish

			for(j = 0; j < adj_processes[i].num_recv_nodes; j++)
			{
				adj_processes[i].recv_locs[j]->disp[r] = adj_processes[i].recv_buf[j]; // copy from recv buffer
			}
		}
		
		MPI_Waitall(num_adj_processes, send_requests, MPI_STATUSES_IGNORE);
	}
	
	//exchange node_forces
	 for(r = 0; r < 2 ; r++)
	 {
		 request_ghost_nodes(recv_requests, comm);
		 
		 for(i = 0; i <  num_adj_processes; i++)
		{
			for(j = 0; j < adj_processes[i].num_send_nodes; j++)
			{
				 adj_processes[i].send_buf[j] = adj_processes[i].send_locs[j]->forces[r]; // copy to send buffer
			}
		}
		
		 send_ghosted_nodes(send_requests, comm);
		 
		for(i = 0; i <  num_adj_processes; i++)
		{
			MPI_Wait(&recv_requests[i], MPI_STATUS_IGNORE); // wait for the Irecv to finish

			for(j = 0; j < adj_processes[i].num_recv_nodes; j++)
			{
				adj_processes[i].recv_locs[j]->forces[r] = adj_processes[i].recv_buf[j]; // copy from recv buffer
			}
		}
		
		MPI_Waitall(num_adj_processes, send_requests, MPI_STATUSES_IGNORE);
	}
	
	//exchange node_field
	 for(r = 0; r < 4 ; r++)
	 {
		 request_ghost_nodes(recv_requests, comm);
		 
		 for(i = 0; i <  num_adj_processes; i++)
		{
			for(j = 0; j < adj_processes[i].num_send_nodes; j++)
			{
				 adj_processes[i].send_buf[j] = adj_processes[i].send_locs[j]->field[r]; // copy to send buffer
			}
		}
		
		 send_ghosted_nodes(send_requests, comm);
		 
		for(i = 0; i <  num_adj_processes; i++)
		{
			MPI_Wait(&recv_requests[i], MPI_STATUS_IGNORE); // wait for the Irecv to finish

			for(j = 0; j < adj_processes[i].num_recv_nodes; j++)
			{
				adj_processes[i].recv_locs[j]->field[r] = adj_processes[i].recv_buf[j]; // copy from recv buffer
			}
		}
		
		MPI_Waitall(num_adj_processes, send_requests, MPI_STATUSES_IGNORE);
	}
	
	if(rank == 0)
	{
		printf("adjacent processes handshake complete...\n");
	}
	
	//free function vars on heap
	free(recv_requests);
	free(send_requests);
	
	 return;
 }
 
void *solve_at_thread(void *arg)
{
	/* This function solves the grid at a thread. It essentially calls grid_update() function at its core */
	int *thread_id_p = (int*) arg;
	int thread_id = *thread_id_p;
	int i;
	int barrier_return;
	
	// do time marching
	for(i = 1; i <= sol.Niter; i++)
	{
		/*************** PTHREADS (ENTRY) SYNC POINT ******************/
		
		// Pthreads synchronization point (w Barrier)
		// ...(wait for the thread 0 to cross the barrier)
		
		barrier_return = pthread_barrier_wait(&barr);
		
		if(barrier_return != 0 && barrier_return != PTHREAD_BARRIER_SERIAL_THREAD)
		{
			printf("Could not wait on barrier\n");
			exit(-1);
		}

		/******************** CALCULATE NEW STATES ********************/

		// call grid update function
		grid_update( i, t_args[thread_id].n0, t_args[thread_id].n1);
		
		/**************** PTHREADS (EXIT) SYNC POINT ******************/
		// Pthreads synchronization point (w Barrier)
		barrier_return = pthread_barrier_wait(&barr);
		
		if(barrier_return != 0 && barrier_return != PTHREAD_BARRIER_SERIAL_THREAD)
		{
			printf("Could not wait on barrier\n");
			exit(-1);
		}
	}
	
	pthread_exit(0); // terminate the thread
}

void solve(int rank, int npes, MPI_Comm comm)
{
	/*
	This function runs explicit time marching on the grid based on the physics setup
	*/
	int_t i, j, temp_int;
	int barrier_return;
	
	// time measurement
	//Declare time measurement related variables
    double time_in_secs = 0;
    double processor_frequency = 0.0;
    double start_cycles= 0.0;
    double end_cycles = 0.0;
    
    if(BGQ)
    {
        processor_frequency = 1600000000.0;
    }
    else
    {
        processor_frequency = 1.0;
    }

	MPI_Request *recv_requests, *send_requests;

	recv_requests = calloc(num_adj_processes, sizeof(*recv_requests));
	send_requests = calloc(num_adj_processes, sizeof(*send_requests));
	
	// initialize pthread attributes to default
	pthread_attr_init(&attr);

    //run the time marching scheme
    if(rank == 0) printf("running solution ...\n");
    
    /******************** SETUP PTHREAD INTERFACE ********************/

	//initialize barrier (before creating threads)
	if(pthread_barrier_init(&barr, NULL, NUM_THREADS))
	{
		printf("ERROR: Could not create a barrier\n");
		exit(-1);
	}

	//setup the thread arguments
	temp_int = num_nodes/NUM_THREADS + 1;
		
	j = 0;
	
	for(i = 0; i < NUM_THREADS; i++)
	{
		 //setup thread args
		 t_args[i].rank = rank; //set the rank
		 t_args[i].thread_id = i;//set the thread id
		 t_args[i].n0 = j; // starting node index
		 
		 /***************************/
		 if(j + temp_int >  num_nodes)
		 {
			 j = num_nodes;
		 }
		 else
		 {
			 j += temp_int;
		 } 
		 /***************************/
		 
		 t_args[i].n1 = j - 1; // ending node index	 
		 
	}

	//create threads (NUM_THREADS - 1)
	for(i = 1; i < NUM_THREADS; i++)
	{
		 if(pthread_create(&p_threads[i - 1], &attr, solve_at_thread, (void *)&t_args[i].thread_id) != 0)
		 {
			 printf("ERROR: Thread creation failed at rank %d\n", rank);
		 }
	}	 
    
    /******************************************************************/
    
    for(i = 1; i <= sol.Niter; i++)
    {
       
        //request for ghost nodes data
        request_ghost_nodes(recv_requests, comm);
        
        /*************** PTHREADS (ENTRY) SYNC POINT ******************/

		// Pthreads synchronization point (w Barrier)
		barrier_return = pthread_barrier_wait(&barr);
		
		if(barrier_return != 0 && barrier_return != PTHREAD_BARRIER_SERIAL_THREAD)
		{
			printf("Could not wait on barrier\n");
			exit(-1);
		}

		/************************* GRID UPDATE ************************/
		//start the time measurement
		start_cycles= (double)GetTimeBase();
		
		grid_update( i, t_args[0].n0, t_args[0].n1);

			
		/**************** PTHREADS (EXIT) SYNC POINT ******************/

		// Pthreads synchronization point (w Barrier)
		barrier_return = pthread_barrier_wait(&barr);
		
		if(barrier_return != 0 && barrier_return != PTHREAD_BARRIER_SERIAL_THREAD)
		{
			printf("Could not wait on barrier\n");
			exit(-1);
		}
		
		//end the time measurement
		end_cycles= (double)GetTimeBase();
		
		//calculate the execution time at each rank
		time_in_secs += (end_cycles - start_cycles)/processor_frequency;
        
        /***************** FINISH MPI DATA EXCHANGE *******************/
        
        //copy ghosted nodes to buf
        copy_ghosted_nodes(i);
        
        //send the ghosted nodes
        send_ghosted_nodes(send_requests, comm);
        
        //update the ghost nodes
        update_ghost_nodes(i, recv_requests, comm);
        
        //wait for send completion
        MPI_Waitall(num_adj_processes, send_requests, MPI_STATUSES_IGNORE);
                
    }
    
    /******************************************************************/

	//...Join threads
	for(i = 1; i < NUM_THREADS; i++)
	{
		pthread_join(p_threads[i - 1], NULL);
	}

	 /*****************************************************************/

        
    /*******************************************************************/
    
    double total_time = 0;
    
		
	MPI_Reduce( &time_in_secs, &total_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	
	if(rank == 0)
	{
		printf("Cummulative computation (only) time = %lf seconds \n", total_time);
		printf("Mean computation (only) time per process = %lf seconds \n", total_time/((double) npes));
		
	} 
    
    /*******************************************************************/
    
	//free function vars on heap
	free(recv_requests);
	free(send_requests);

    return;
}

void mpi_result_output(char *filename, int_t rank, MPI_Comm comm)
{
    /*
    This function outputs (MPI_IO) the nodal solution for visualization and post processing
    
    -- verified (5/4/2019)
    */

	int_t i, j;
	real_t *disp = calloc(num_nodes, sizeof(*disp));
	
	//copy the displacements
	j = nodes[0].iter%2; //even or odd iteration
	for(i=0; i < num_nodes; i++)
	{
		disp[i] = nodes[i].disp[j];
	}
	 
	// declare file handle and offset
	MPI_File fh;
    MPI_Offset rank_offset;
    
    //open IO file
	if( MPI_File_open(comm, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh) == MPI_SUCCESS)
	{
		//calculate the rank offset
		rank_offset = node_offset*sizeof(real_t);
		
		//perform blocking write operation
		MPI_File_write_at(fh, rank_offset, disp, num_nodes, mpi_real_t, MPI_STATUS_IGNORE);
	
		//close IO file
		MPI_File_close(&fh);
	}
	else
	{
		printf("Error opening MPI_File for grid output\n");
	}
	
	//free function vars on heap
	free(disp);
	
	return;

}

void free_heap()
{
	/* This function frees the heap variables - nodes, ghost_nodes, 
	 * adj_processes */
	 int_t i;
	 
	 //free nodes
	 free(nodes);
	 
	 //free ghost_nodes
	 free(ghost_nodes);
	 
	 
	 //free adj_processes
	 for(i=0; i < num_adj_processes; i++)
	 {
		 free(adj_processes[i].recv_buf);
		 free(adj_processes[i].send_buf);
		 free(adj_processes[i].recv_locs);
		 free(adj_processes[i].send_locs);
	 }
	 free(adj_processes);
	 
	 
	 return;
 }
 
/**********************************************************************/
//MAIN FUNCTION DEF
int main(int argc, char **argv)
{


        
    //Declare time measurement related variables
    double time_in_secs = 0;
    double processor_frequency = 0.0;
    double start_cycles= 0.0;
    double end_cycles = 0.0;

    if(BGQ)
    {
        processor_frequency = 1600000000.0;
    }
    else
    {
        processor_frequency = 1.0;
    }
    
    //Initialize MPI Interface
    MPI_Init(&argc, &argv);
    
    // Declare MPI related variables
    int_t rank, npes;
    MPI_Comm comm;

    // get the comm_size and process rank
    MPI_Comm_size(MPI_COMM_WORLD, (int*) &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, (int*) &rank);
    
    //check if the inputs are correct
    if(rank == 0 && argc < 7)
    {
        printf("ERROR: Only %d arguments given, 6 expected\n 1) g_parts filename\n 2) g_xadj filename\n 3) g_adjncy filename\n 4) g_global_vert_ids filename\n 5) g_nodal_data filename\n 6) output filename\n", argc - 1);
        exit(-2);
    }
    
    //Read input data
    if(rank == 0)
    {
		read_data(argv, npes);
	}
	else
	{
		receive_data(rank, MPI_COMM_WORLD);
	}

	//prepare the node to rank mapping
	if(rank == 0)
	{
		g_n2r_map = calloc(nvtxs, sizeof(*g_n2r_map));
		set_n2r_map(npes, g_parts, g_global_vert_ids, g_n2r_map);
	}
	
	//Distribute data to all processes from rank 0
	if(rank == 0)
	{
		send_data(npes, MPI_COMM_WORLD);
	}
	
	// Initialize the grid, Process data at each rank
	grid_initialize(rank);
	
	//set up the physics
    set_physics();

    //set up solution
    set_solution();
    
    // change the communicator to graph
    if(GRAPH_COMM == 1)
    {
		get_graph_comm(MPI_COMM_WORLD, &comm);
	}
	else
	{
		comm = MPI_COMM_WORLD;
	}
    
    //perform initial exchange between MPI ranks (initialize ghost nodes)
    init_ghost_nodes(rank, comm);
	
	/**************************** RUN SOLUTION ************************/
	
	//start the time measurement
	start_cycles= (double)GetTimeBase();
	
    //run the solution
    solve(rank, npes, comm);

	//end the time measurement
	end_cycles= (double)GetTimeBase();
	
	//calculate the execution time at each rank
	time_in_secs = (end_cycles - start_cycles)/processor_frequency;
	
	//print solution time
	if(rank == 0)
	{
		printf("Parallel execution time (computation  + communication) = %lf seconds\n", time_in_secs);
	}
	
	/**************************** OUTPUT RESULTS ***********************/
	//start the time measurement
	start_cycles= (double)GetTimeBase();
    
    //output the solution (binary)
    mpi_result_output(argv[6], rank, comm);
    
    //end the time measurement
	end_cycles= (double)GetTimeBase();
	
	//calculate the execution time at each rank
	time_in_secs = (end_cycles - start_cycles)/processor_frequency;
	
	//print IO time
	if(rank == 0)
	{
		printf("MPI IO completed in %lf seconds\n", time_in_secs);
	}
    
    /******************************************************************/

    // Close the MPI interface
    MPI_Finalize();

    //free the heap
    free_heap();
    
	return 0;
}


