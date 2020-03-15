#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <metis.h>
#include <arpa/inet.h>

/**********************************************************************
                         MESH PARTITION - 2D
By: VINEET NEGI

Description -
This code uses METIS to partition a mesh feed in the form of
Compressed Storage Format. There are three inputs -
1) xadj - an array which contains the location of adjacency in adjncy array
2) adjncy - an array which stores the neighbors of a node
3) nodal_coords - an array storing the nodal coordinates

output -
g_parts - reference to g_xadj for each partition
g_xadj - reference to g_adjncy for each node
g_adjncy - an array of vertices adjacancies
g_global_vert_ids - an array of global node ids
g_nodal_data - a array of nodal data (including the nodal positions)


History -
02/12/2019 - work on the code started
04/02/2019 - code verified
04/08/2019 - output files are changed to binaries from ASCII. Endianess taken into account
***********************************************************************/
// DEFINE MACROS
#define NUM_NODAL_DATA 5


/***********************************************************************/
// DEFINE GLOBAL DATASTRUCTURES
/*Input datastructure*/
idx_t nvtxs;
idx_t ncon;
idx_t *xadj;
idx_t *adjncy;
real_t *nodal_data;

idx_t nparts;
idx_t objval;
idx_t *part;

/*Output datastructure ... start with 'w' prefix */
int *w_parts;
int *w_xadj;
int *w_adjncy;
int *w_global_vert_ids;
double **w_nodal_data;

/***********************************************************************/
// DEFINE GLOBAL FUNCTIONS

void readinput(char *fn_xadj , char *fn_adjncy, char *fn_nodal_data)
{
    /*
    This function reads the input files
    */
    int i, j;
    FILE *input_file;
    int fsize = 0;
    int temp_int;
    double temp_dbl;

    //update
    printf("Reading Input files\n");

    // READ xadj

	if( (input_file = fopen(fn_xadj, "r")) == NULL)
    {
        printf("Failed to open xadj data file: %s \n", fn_xadj);
        exit(-1);
    }
    fscanf(input_file, "%d", &fsize); // get the number of lines to read
    xadj = calloc(fsize, sizeof(*xadj)); // allocate space to xadj array
    nvtxs = (idx_t) (fsize - 1);

    for(i=0; i<fsize; i++)
    {
        fscanf(input_file, "%d", &temp_int);
        xadj[i] = (idx_t) temp_int;
    }

    fclose(input_file);

    printf("xadj input file successfully read\n");

    // READ adjncy

	if( (input_file = fopen(fn_adjncy, "r")) == NULL)
    {
        printf("Failed to open adjncy data file: %s \n", fn_adjncy);
        exit(-1);
    }
    fscanf(input_file, "%d", &fsize); // get the number of lines to read
    adjncy = calloc(fsize, sizeof(*adjncy)); // allocate space to adjncy array

    for(i=0; i<fsize; i++)
    {
        fscanf(input_file, "%d", &temp_int);
        adjncy[i] = (idx_t) temp_int;
    }

    fclose(input_file);

    printf("adjncy input file successfully read\n");

    // READ nodal data

	if( (input_file = fopen(fn_nodal_data, "r")) == NULL)
    {
        printf("Failed to open nodal_data data file: %s \n", fn_nodal_data);
        exit(-1);
    }
    fscanf(input_file, "%d", &fsize); // get the number of lines to read
    nodal_data = calloc((2 + NUM_NODAL_DATA)*fsize, sizeof(*nodal_data)); // for 2D

    for(i=0; i<fsize; i++)
    {
		for(j = 0; j < 2 + NUM_NODAL_DATA; j++)
		{
			fscanf(input_file, "%lf", &temp_dbl); //scan the value
			nodal_data[(2 + NUM_NODAL_DATA)*i + j] = (real_t) temp_dbl;
		 }
	 }

    fclose(input_file);

    printf("nodal_coords input file successfully read\n");

    return;

}

void get_partition(int num_parts)
{
	/*
	 * This function call METIS API and gets the partition
	 */

	 int s;

	 //specify the output vars
	 ncon = 1;
	 nparts = num_parts;
	 part = calloc((int)nvtxs, sizeof(*part)); // allocate space to part array

	 printf("...calling METIS for partitioning of %d nodes in %d domains...\n", (int) nvtxs, num_parts);
	 
	 if(num_parts > 1)
	 {
		s = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts, NULL, NULL, NULL, &objval, part);
	 

	 if(s == METIS_OK)
	 {
		 printf("grid partitioning finished successfully\n");
	 }
	 else
	 {
		 printf("METIS returned with error %d\n", s);
	 }
	 
	 }

	 return;
 }
	
void get_num_part_vert(int *num_part_vert)
{
	/*
	 * This function gets the number of vertices associated per part
	 */
	 int i, j;
	 
	 for(i=0; i < (int) nvtxs; i++)
	 {
		 j = (int) part[i];
		 num_part_vert[j] += 1;
	 }
	 		
	 return;
 }
 
void prepare_results()
{
	/*
	 * This function prepares the final results. It prepares 4 DS
	 * 
	 * parts - This is an array of size - num partition. It points to 
	 * read location on xadj array
	 * 
	 * xadj - This is an array of size num vertices + 1. It points to the
	 * adjncy array
	 * 
	 * adjncy - This is an array which stores the neighboring vert. of a 
	 * vertex
	 * 
	 * nodal_data - This is an array of size num vertices x num field. 
	 * It stores all the nodal properties.
	 */
	 int i, j, p, v, temp_int;
	 int temp_int1, temp_int2;
	 int counter1, counter2;
	 int *num_part_vert = calloc(nparts, sizeof(*num_part_vert));
	 
	 int *part_feed_counter = calloc(nparts, sizeof(*part_feed_counter));;
	 int **part_vertices;
	 
	 //get the number of vertices associated with a part
	 get_num_part_vert(num_part_vert);
	 
	 //allocate space
	 part_vertices = calloc(nparts, sizeof(*part_vertices));
	 for(i=0; i<(int) nparts; i++)
	 {
		 part_vertices[i] = calloc(num_part_vert[i], sizeof(*part_vertices[i]));
	 }
	 
	 //separate the vertices as per their partition
	 for(i=0; i < (int) nvtxs; i++)
	 {
		 p = (int) part[i];
		 part_vertices[p][part_feed_counter[p]] = i;
		 part_feed_counter[p] += 1;
	 }
	 
	 //allocate memory to the output data
	 w_parts = calloc(nparts + 1, sizeof(*w_parts));
	 w_xadj = calloc(nvtxs + 1, sizeof(*w_xadj));
	 w_global_vert_ids = calloc(nvtxs, sizeof(*w_global_vert_ids));
	 w_adjncy = calloc(xadj[nvtxs], sizeof(*w_adjncy));
	 w_nodal_data = calloc(NUM_NODAL_DATA + 2, sizeof(*w_nodal_data));
	 for(i=0; i<NUM_NODAL_DATA + 2; i++)
	 {
		 w_nodal_data[i] = calloc(nvtxs, sizeof(*w_nodal_data[i]));
	 }
	 
	 //fill the output data
	 w_parts[0] = 0;
	 w_xadj[0] = 0;
	 temp_int = 0;
	 counter1 = 0;
	 counter2 = 1;
	 for(i = 1; i < nparts + 1; i++)
	 {
		 w_parts[i] = temp_int + num_part_vert[i - 1]; // add to w_parts
		 temp_int += num_part_vert[i - 1];
		 
		 for(j = 0; j < num_part_vert[i - 1]; j++)
		 {
			 v = part_vertices[i - 1][j]; // get the global vertex ID
			 w_global_vert_ids[counter2 - 1] = v; //add to w_global_vert_ids
			 
			 //fix the nodal data (add to the nodal data)
			 for(p=0; p<(2+NUM_NODAL_DATA);p++)
			 {
				w_nodal_data[p][counter2 - 1] = (double) nodal_data[(2 + NUM_NODAL_DATA)*v + p];
			}
					 
			 temp_int1 = (int) xadj[v];
			 temp_int2 = (int) xadj[v + 1];
			 
			 for(p = temp_int1; p < temp_int2; p++) //reusing p
			 {
				 w_adjncy[counter1] = adjncy[p]; // add to w_adjncy
				 counter1++;
			 }
			 w_xadj[counter2] = counter1; // add to w_xadj
			 counter2++;
		 }
		 
	 }
			 
 
	 //free memory
	 for(i=0; i<nparts;i++) // free 2D array or double pointer
	 {
		 free(part_vertices[i]);
	 }
	 free(num_part_vert);
	 free(part_feed_counter);
	 free(part_vertices);
	 
	 return;
 }
 
size_t endian_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream, int endian)
{
	int i;
	int s0, s1;
	
	unsigned char temp;

    if ( endian == 0)
    {
        /* Little endian machine, use fwrite directly */
        return fwrite(ptr, size, nmemb, stream);
    }
    else
    {
        /* Big endian machine, pre-process first */
        
		unsigned char *buffer = (unsigned char*) ptr;
		
        for(i=0; i<nmemb; i++)
        {           
			//reverse the order of bytes
			s0 = ((int) size)*i; //start idx
			s1 = ((int) size)*i + ((int) size) - 1; //end idx
			while(s1 > s0)
			{
				temp = buffer[s1];
				buffer[s1] = buffer[s0];
				buffer[s0] = temp;
				s1--;
				s0++;
			}
		}
        return fwrite(ptr, size, nmemb, stream);
    }  
}
 
 void write_results(int ENDIAN)
 {
	 /*
	  * This function writes the prepared results to a file.
	  * Note - all the files start with prefix 'g'
	  */
	  FILE *out = NULL;
	  int i, temp_int = w_xadj[nvtxs];
	  	  
	  //write the partitions

	  out = fopen("g_parts", "wb");
	  endian_fwrite((void*) w_parts, sizeof(*w_parts), nparts + 1, out, ENDIAN);
	  //fwrite(w_parts, sizeof(*w_parts), nparts + 1, out);
	  fclose(out); 
	  
	  
	  //write the xadj
	  out = fopen("g_xadj", "w");
	  endian_fwrite((void*) w_xadj, sizeof(*w_xadj), nvtxs +1, out, ENDIAN);
	  //fwrite(w_xadj, sizeof(*w_xadj), nvtxs +1, out);
	  fclose(out);
	  
	  	  
	  //write the adjncy
	  out = fopen("g_adjncy", "w");
	  endian_fwrite((void*) w_adjncy, sizeof(*w_adjncy), temp_int, out, ENDIAN);
	  //fwrite(w_adjncy, sizeof(*w_adjncy), w_xadj[nvtxs], out);
	  fclose(out);
	  	  
	  //write the w_global_vert_ids
	  out = fopen("g_global_vert_ids", "w");
	  endian_fwrite((void*) w_global_vert_ids, sizeof(*w_global_vert_ids), nvtxs, out, ENDIAN);
	  //fwrite(w_global_vert_ids, sizeof(*w_global_vert_ids), nvtxs, out);
	  fclose(out);  

	  //write the w_nodal_data
	  out = fopen("g_nodal_data", "w");
	  
	  for(i = 0; i < 2+NUM_NODAL_DATA; i++)
	  {
		endian_fwrite((void*) w_nodal_data[i], sizeof(*w_nodal_data[i]),  nvtxs, out, ENDIAN);
		//fwrite(w_nodal_data[i], sizeof(*w_nodal_data[i]), nvtxs, out);

	  }
	  fclose(out); 
  }
	  
/***********************************************************************/
// DEFINE MAIN FUNCTION

int main(int argc, char **argv)
{
	int i;
    int num_parts;
    int ENDIAN;

    //check if the inputs are correct
    if(argc < 6)
    {
        printf("Enter input files. Only %d arguments given, 5 expected -> xadj, adjncy, nodal_data, #partitions, ENDIANNESS (0 for little, 1 for big) \n", argc - 1);
        exit(-1);
    }

	// read inputs
    readinput(argv[1], argv[2], argv[3]);
    
    //set the endianness
    ENDIAN = atoi(argv[5]); // 0 for little endian & 1 for big endian
    
    if(ENDIAN  == 0)
    {
		printf("Little Endian output format selected\n");
	}
	else if(ENDIAN == 1)
	{
		printf("Big Endian output format selected\n");
	}
	else
	{
		printf("ERROR: unknown endianness format\n");
		exit(-1);
	}
	
    
    //calculate the partition
    num_parts = atoi(argv[4]);
    get_partition(num_parts);
       
    //preparea and write results
    prepare_results();
    write_results(ENDIAN);

    //free the memory
    free(xadj);
    free(adjncy);
    free(nodal_data);
    free(part);
    free(w_parts);
    free(w_xadj);
    free(w_adjncy);
    free(w_global_vert_ids);
    for(i=0; i < 2+NUM_NODAL_DATA; i++) // free 2D array or double pointer
    {
		free(w_nodal_data[i]);
	}
    free(w_nodal_data);
   
    return 0;

}
