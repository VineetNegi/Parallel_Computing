
/**********************************************************************
                         Utility functions
By: VINEET NEGI

Description -
This is header file which contains utility functions - 
* quick sort (int)
* binary search (int)
* unique (int)

*** these are copied from Geeks for Geeks website

History -

***********************************************************************/

// Define functions
int_t unique(int_t* arr, int_t n) 
{ 
	/* This function removes the duplicates in-place and return the size
	 of the new array
	*/
	int_t i, j;
	
    if (n==0 || n==1)
    { 
        return n; 
	}
  
    // To store index of next unique element 
    j = 0; 
  
    // Doing same as done in Method 1 
    // Just maintaining another updated index i.e. j 
    for (i=0; i < n-1; i++)
    {
		// If current element is not equal 
        // to next element then store that 
        // current element 
         
        if (arr[i] != arr[i+1])
        { 
            arr[j++] = arr[i]; 
		}
	}
		
    // Store the last element as whether 
    // it is unique or repeated, it hasn't 
    // stored previously  
    arr[j++] = arr[n-1]; 
  
    return j; 
}

int_t binarySearch(int_t arr[], int_t l, int_t r, int_t x) 
{ 
	/* A iterative binary search function. It returns location of x in 
	 * given array arr[l..r] if present, otherwise -1 
	 */
	 
    while (l <= r) { 
        int_t m = l + (r - l)/2; 
  
        // Check if x is present at mid 
        if (arr[m] == x) 
            return m; 
  
        // If x greater, ignore left half 
        if (arr[m] < x) 
            l = m + 1; 
  
        // If x is smaller, ignore right half 
        else
            r = m - 1; 
    } 
  
    // if we reach here, then element was 
    // not present 
    return -1; 
} 
  
void swap(int_t* a, int_t* b) 
{ 	
	/* A utility function to swap two elements  */

	int_t t = *a; 
	*a = *b; 
	*b = t; 
	return;
} 

int_t partition(int_t arr[], int_t low, int_t high) 
{ 
	/* This function takes last element as pivot, places the pivot 
	 * element at its correct position in sorted array, and places all
	 * smaller (smaller than pivot) to left of pivot and all greater 
	 * elements to right of pivot
	 */
	 
	int_t pivot = arr[high]; // pivot 
	int_t i = (low - 1); // Index of smaller element 
	int_t j;

	for (j = low; j <= high - 1; j++) 
	{ 
		// If current element is smaller than or 
		// equal to pivot 
		if (arr[j] <= pivot) 
		{ 
			i++; // increment index of smaller element 
			swap(&arr[i], &arr[j]); 
		} 
	} 
	swap(&arr[i + 1], &arr[high]); 
	return (i + 1); 
}

void quickSort(int_t arr[], int_t low, int_t high) 
{
	/* The main function that implements QuickSort
	 * arr[] --> Array to be sorted, 
	 * low --> Starting index, 
	 * high --> Ending index 
	 */
	 
	if (low < high) 
	{ 
		/* pi is partitioning index, arr[p] is now 
		at right place */
		int_t pi = partition(arr, low, high); 

		// Separately sort elements before 
		// partition and after partition 
		quickSort(arr, low, pi - 1); 
		quickSort(arr, pi + 1, high); 
	} 
	
	return;
} 

/***********************************************************************/
