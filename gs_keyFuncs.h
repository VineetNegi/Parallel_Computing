

/**********************************************************************
                        KEU FUNCTIONS
By: VINEET NEGI

Description -
This is header file which contains all the key functions required to 
solve the physics

* grid_update - this function holds the grid update rule
* get_rel_position_tag - this function gets the position tag (loc in 
neighors array) of one node w.r.t. other
* set_physics - Initializes the physics
* set_solution - Initializes the solution
* 
History -
created on 4/5/2019

***********************************************************************/
// DEFINE DATATYPES KEY FUNCTIONS

/***************************** GRID RELATED ***************************/
int get_rel_position_tag(int_t n0, int_t n1)
{
	/* this function calclates the relative position of point p0 w.r.t. 
	 * point p1
	 * n0 = global id of p0
	 * n1 = global id of p1
	 * (modify it accordingly)
	 */
	 
	 if(abs(n1 - n0) == 1) // node along x
	 {
		 return (n1 > n0) ? 1: 0;
	 }
	 else //node along y
	 {
		 return (n1 > n0) ? 3: 2;
	 }
 }
 
/************************** PHYSICS RELATED ***************************/ 
 
void set_physics()
{
    /*
    This function reads in the physics and sets the system for solve
    */

    phy.conductivity = 1.0;
    phy.density = 1.0;
    phy.heat_capacity = 1.0;

    return;
}

/************************* SOLUTION RELATED ***************************/ 
void set_solution()
{
    /*
    This function reads in and sets the solution parameters
    */

    sol.dt = DT;
    sol.dx = DX;
    sol.dy = DY;
    sol.Niter = NUM_ITER;
    
    // calculate relevant physics and solution coefficients and setup the run
	sol.sx = (phy.conductivity/(phy.density*phy.heat_capacity))*(sol.dt/(sol.dx*sol.dx));
	sol.sy = (phy.conductivity/(phy.density*phy.heat_capacity))*(sol.dt/(sol.dy*sol.dy));
	sol.sq = sol.dt/(phy.density*phy.heat_capacity);
	sol.fx = sol.dt/(2.0*sol.dx);
	sol.fy = sol.dt/(2.0*sol.dy);

    return;
}

void grid_update(int_t it, int_t n0, int_t n1)
{
	/* This function performs the main grid update */
	
	int_t j;
	int_t it_present, it_previous;

	
	it_present = it%2; //determine if the current iteration is odd or even
    it_previous = (it - 1)%2; //determine if the previous iteration was odd or even
	
	for(j = n0; j <= n1; j++)
	{
		nodes[j].iter = it; // update the iteration number of the nodes
				
		// specify the time marching rule
		if(nodes[j].nodetype == 0)
		{
			
			nodes[j].disp[it_present] =
			(1.0 - 2.0*sol.sx - 2.0*sol.sy)*nodes[j].disp[it_previous] +
			sol.sx*(nodes[j].neighbors[0]->disp[it_previous] + nodes[j].neighbors[1]->disp[it_previous]) + /* heat conduction in X */
			sol.sy*(nodes[j].neighbors[2]->disp[it_previous] + nodes[j].neighbors[3]->disp[it_previous]) + /* heat conduction in Y */
			sol.sq*nodes[j].forces[it_previous] + /* heat generation */
			-sol.fx*nodes[j].field[2*it_previous]*(nodes[j].neighbors[1]->disp[it_previous] - nodes[j].neighbors[0]->disp[it_previous]) + /* heat convection in X */
			-sol.fy*nodes[j].field[2*it_previous + 1]*(nodes[j].neighbors[3]->disp[it_previous] - nodes[j].neighbors[2]->disp[it_previous]); /* heat convection in Y */

		}
	}
	
	return;
}

/***********************************************************************/
