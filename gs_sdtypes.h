

/**********************************************************************
                        Primary Datatypes
By: VINEET NEGI

Description -
This is header file which contains all the Secondary datatypes

* node - to stores nodal info
* process - to store info about a process (rank)
* physics - to store into about a physics
* solution - to store info about a solution

History -
created on 4/5/2019

***********************************************************************/
// DEFINE SECONDARY DATATYPES


// template of a node
struct node
{
    // This struct provides a template for storing all nodal info

    //Node identification and classification
    int_t id; //global node id
    int_t nodetype; //type of node: 0 - free, 1 - not free

    //node location
    real_t loc[2];

    //nodal solution value
    real_t disp[2];

    //nodal forces
    real_t forces[2];

    //specify the nodal vector field
    real_t field[4];

    //iteration info
    int_t iter;

    //nodal neighborhood info. This is an array of pointers to the neighboring nodes
    struct node *neighbors[NODAL_Z];

    /*
    status with reference to parallel implementation. Here,
    0 - node is owned and not shared
    1 - node is owned and shared
    2 - node is not owned
    */
    int_t rank_status;

};

//template of process (neighbor)
struct process
{
    // This struct provides a template for defining a neighboring process
    // for exchange of data

    int_t rank; //rank of the neighboring process
    int_t num_recv_nodes; // no. of nodes to be received
    int_t num_send_nodes; // no. of nodes to be sent
    real_t *recv_buf; // buffer for receive data
    real_t *send_buf; // buffer for send data
    struct node **recv_locs; // locations for received data
    struct node **send_locs; //locations for sent data
};

//template of a thread args
struct thread_args
{
	//This struct provides a template to define a composite argument
	// to be passed to threads
	
	int rank; // rank of the process
	int thread_id; //id of the thread
	int n0; // starting idx of the nodes to be processed
	int n1; // end index of the nodes to be processed
	
};

// template of a thermal physics
struct physics
{
    // This struct provides a template for describing heat conduction problem

    real_t density, heat_capacity, conductivity;

};

//template of solution
struct solution
{
    //This struct holds the setting and parameters required for solution
    real_t dt, dx, dy; //discretization in time, x, and y dimensions
    int_t Niter; //Number of iterations
    
    //miscs vars (specific to physics
    real_t sx, sy, sq, fx, fy;

};

/***********************************************************************/
