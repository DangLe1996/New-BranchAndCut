

#include "headers.h"

int main(int argc, char *argv[])
{

	/*runInstance(argv[1]);*/


	return 0;
}


EXPORT int add(int i, int j) {
	return i + j;
}

EXPORT void runInstance(char* filename, __int32* sequence[]) {
	clock_t  start, end;
	double opt_value;
	double cputime;

	n_conn_comp = 0;
	old_objval = 0;
	count_same_node = 0;
	n_int_feas_subproblems = 0;
	n_frac_feas_subproblems = 0;
	n_feas_subproblems = 0;
	printf("File name is  %s " , filename);
	start = clock();
	read_INSTANCE_AMAZON(filename);
	printf("Solving TSP with B&C algorithm: \n");

	opt_value = solve_TSP(sequence);
	end = clock();
	cputime = (double)(end - start) / (double)CLOCKS_PER_SEC;   //Compute CPU time

	free_memory();


}

EXPORT void read(__int32* input, size_t size)
{
	int i;
	for (i = 0;i < size;i++)
		input[i] = i;
}

