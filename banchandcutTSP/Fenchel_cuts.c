#include "headers.h"


/// INFO ///////////////////////////////////////////////////////////////////////
//
// master.c: Master Problem Functions
//
// Author:  Carlos Luna-Mota 
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "/home/marlui/experiments/concorde/concorde.h"	


// CREATE LINE VECTOR//

LINE *create_line_vector(int cells) {
    LINE *ptr = (LINE *) calloc(cells, sizeof(LINE));
    if (ptr == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for int_vector\n");
    }
    return ptr;
}

void ///Funcion de estructura de distancias para concorde
fill_dat ( int ncount, CCdatagroup * dat, int **adj )
{

    int i, j;

    CCutil_init_datagroup ( dat );
    CCutil_dat_setnorm ( dat, CC_MATRIXNORM );

	
    dat->adj = CC_SAFE_MALLOC ( ncount, int * );
    dat->adjspace = CC_SAFE_MALLOC ( ncount * ( ncount + 1 ) / 2, int );
    if ( dat->adj == ( int ** ) NULL || dat->adjspace == ( int * ) NULL )
    {
        CC_IFFREE ( dat->adj, int * );
        CC_IFFREE ( dat->adjspace, int );
    }
    for ( i = 0, j = 0; i < ncount; i++ )
    {
        dat->adj[i] = dat->adjspace + j;
        j += ( i + 1 );
    }
    for ( i = 0; i < ncount; i++ )
    {
        for ( j = 0; j <= i; j++ )
        {
            dat->adj[i][j] = adj[i][j];
        }
    }
    return;
}


//Heuristico de Concorde
int heuristic_tsp_Concorde(GLOBAL_INFO global,double *dist) {
	//////For lazyness/////
	//clock_t inicio,final;
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	double  ovStep2;//=global.results->ovStep2;
	int status= 0;
	int kk=0;int m=0;
	int semente=rand();
    double szeit, *mybnd, *mytimebound;
    int success, foundtour, hit_timebound = 0;
    int *in_tour = (int *) NULL;
    int *out_tour = (int *) NULL;
    CCrandstate rstate;
    char *probname = (char *) NULL;
    char *probname0 = (char *) NULL;
    static int run_silently = 1;
    CCutil_sprand(semente, &rstate);
    mybnd = (double *) NULL;
   	mytimebound = (double *) NULL;
    int *elist = create_int_vector(E*2);
    int *elen =create_int_vector(E);
    int edge = 0;
    int edge_peso = 0;
    int e;
    
    int i=0;	
    for(e=0;e<E;e++){
		elist[i]=global.data->index_i[e];
		elist[i+1]=global.data->index_j[e];
	}

    		
    CCdatagroup dat;
     int **adj ;
     int rval = 0;
     double val;
     int j;  
    	      
     adj=create_int_matrix(N,N);
    	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
		  if(i!=j) adj[i][j]=(int) round(dist[global.data->index_e[i][j]]*10000);
		}
		} 	
				
	 fill_dat ( N, &dat, adj );
     out_tour= CC_SAFE_MALLOC (N, int); 
     rval=    CClinkern_tour (N, &dat, E, elist, 10000000,N,in_tour,out_tour,&val,1,-1,-1,(char *) NULL,  CC_LK_RANDOM_KICK, &rstate);  
                  
     val=val/10000;	
     printf("val %f\n",val);
 
		
     for(kk=0;kk<E;kk++){
		global.results->yStep2[kk]=0;
	 }
		
	 for(kk=0;kk<N-1;kk++){
		global.results->yStep2[global.data->index_e[out_tour[kk]][out_tour[kk+1]]]=1;
	 }
	 global.results->yStep2[global.data->index_e[out_tour[0]][out_tour[N-1]]]=1;
	 
	 szeit = CCutil_zeit();		
	 free(out_tour);
	 free(probname);
	 free(elist);
	 free(elen);
	 return status;
}





int solve_tsp_Concorde(GLOBAL_INFO global) {
	
	struct timespec time1,time2;
	
	 clock_gettime(CLOCK_MONOTONIC, &time1);
	
	//////For lazyness/////
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	double  ovStep2;//=global.results->ovStep2;
	
	///////////////////////
	
		int status= 0;
		int kk=0;int m=0;
		int semente=rand();
		semente = rand();
    	double szeit, *mybnd, *mytimebound;
    	int success, foundtour, hit_timebound = 0;
    	int *in_tour = (int *) NULL;
    	int *out_tour = (int *) NULL;
    	CCrandstate rstate;
    	char *probname = (char *) NULL;
    	char *probname0 = (char *) NULL;
    	static int run_silently = 1;
    	CCutil_sprand(semente, &rstate);
    	mybnd = (double *) NULL;
   		mytimebound = (double *) NULL;
    	int *elist = create_int_vector(E*2);
    	int *elen =create_int_vector(E);
    	int edge = 0;
    	int edge_peso = 0;
    	int e;
    	int **adj ;
		int j,i;
		CCdatagroup dat;
    	
		for (kk = 0; kk < N; kk++) {
      			for ( m = 0; m < N; m++) {
				if (kk < m) {
				  elist[edge] = kk;
	 			  elist[edge + 1] = m; 
	  			  elen[edge_peso] =-(int) round(global.results->alpha_cur[global.data->index_e[kk][m]]*10000);///-alpha_cur[index_e[kk][m]];
	  			  edge_peso++;
	 			  edge = edge + 2;
				}
     			 }
    		}
 

		int minimo=1;
		for(kk=0;kk<E;kk++){
		if(elen[kk]<minimo){
		minimo=elen[kk];
		}
		}

		for (kk = 0; kk < E; kk++) {
			elen[kk]=elen[kk]-minimo;
	  	}	
 
		
    	adj=create_int_matrix(N,N);
		kk=0;
		for(i=0;i<N-1;i++){
			for(j=i+1;j<N;j++){
				if(i<j){ 
					adj[i][j]=elen[kk];
					adj[j][i]=elen[kk];						 
				}
			kk+=1;
			}
		} 	
    
    	fill_dat(N,&dat,adj );	
    		
    	out_tour = CC_SAFE_MALLOC (N, int); 
	    probname = CCtsp_problabel(" "); 
	    
		status=CCtsp_solve_dat (N, &dat,in_tour,
         out_tour, NULL, &ovStep2, &success,
        &foundtour, probname,   mytimebound, &hit_timebound,
        run_silently, &rstate);
		ovStep2=-(ovStep2+(N)*minimo)/10000;
		global.results->ovStep2=ovStep2;
		printf("status %d",status);
		printf("success %d\n",success);
		printf("optval %f\n",ovStep2);
		printf("foundtour %d\n",foundtour);

		
		
		///Tras resolver con concorde, guardamos el tour en yStep2 [global]///
		for(kk=0;kk<E;kk++){
		global.results->yStep2[kk]=0;
		}
		
		for(kk=0;kk<N-1;kk++){
		global.results->yStep2[global.data->index_e[out_tour[kk]][out_tour[kk+1]]]=1;
		}
		global.results->yStep2[global.data->index_e[out_tour[0]][out_tour[N-1]]]=1;
	
		szeit = CCutil_zeit();
		
		
		free(out_tour);
		free(probname);
		free(elist);
		free(elen);
		
		clock_gettime(CLOCK_MONOTONIC, &time2);
	  global.results->concorde_time = (time2.tv_sec - time1.tv_sec);
	  global.results->concorde_time+= (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;
	
		return status;

}


//*Version sin concorde resolviendo con MTZ
////////////////////////////////////////////////////////////////////////////////
/*int solve_tsp(GLOBAL_INFO global) {
	int status = 0;
	int i, j,e;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	//double best_upper_bound, best_lower_bound;

	//Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *obj;    // objective function coefficients ..............................
	double    *rhs;    // right and side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zelo element ...................
	double    *matval; // coefficient values fo the non-zero elements of constraints....
	double    *lb;     // lower bounds of variables.....................................
	double    *ub;     // upper bounds of variables.....................................
//	int       status;  // optimization status......................... .................
	double    *x;      // solution vector (double, even if the problem is integer) .....
	char probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	//double    value;   // objevtive value of solution ..................................
	int       num_z_var, num_x_var, whichmodel;
	int       **pos_y;
	int       *pos_u;
	//float	  **coef;
	//Just for lazyness:
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	int		*index_i = global.data->index_i;
	int		*index_j = global.data->index_j;
	int		**index_k = global.data->index_k;
	double	**W = global.data->W;
	double	**C = global.data->C;
	double  **F = global.data->F;
	double	**C1 = global.data->C1;
	double	**C2 = global.data->C2;
	COMMODITY *Com = global.data->Com;
	char	**colname;

	for ( e = 0; e < E; e++)global.results->yStep2[e] = 0;

	//Initialize CPLEX environment
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}

	CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, 2);

	// Create the problem in CPLEX 
	strcpy(probname, "TSP");
	lp = CPXcreateprob(env, &status, probname);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}

	CPXchgobjsen(env, lp, CPX_MAX);

	pos_y = create_int_matrix(N,N);
	//coef = create_int_matrix(N, N);
	pos_u = create_int_vector(N);
	x = create_double_vector(N*N + N);
	//for (int i = 0; i < N-1; i++){
	//	for (j = i+1; j < N; j++){
	//		coef[i][j] = global.results->alpha_cur[index_e[i][j]];
			coef[j][i] = global.results->alpha_cur[index_e[i][j]];
	//	}
	//}

	
	index1 = 0;  // index of columns
	numcols = N*(N - 1);
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	colname = (char**)calloc(numcols, sizeof(char*));
	for (i = 0; i < numcols; i++){
		colname[i] = (char*)calloc(255, sizeof(char));
	}


	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			if (i != j){
				pos_y[i][j] = index1;
				obj[index1] =global.results->alpha_cur[index_e[i][j]];
				lb[index1] = 0;
				ub[index1] = 1;
				ctype[index1] = 'B';
				sprintf(colname[index1], "y_%d,%d ", i , j );
				index1++;
			}
		}
	}

	status = CPXnewcols(env,lp, index1, obj, lb, ub, ctype, colname);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(colname);

	//Define u_k variables
	index1 = 0;  // index of columns
	numcols = N;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	colname = (char**)calloc(numcols, sizeof(char*));
	for (i = 0; i < numcols; i++){
		colname[i] = (char*)calloc(255, sizeof(char));
	}

	for (i = 0; i<N; i++){
		pos_u[i] = N*(N-1) + index1;
		obj[index1] =0;
		lb[index1] = -CPX_INFBOUND;
		ub[index1] = CPX_INFBOUND;
		ctype[index1] = 'C';
		sprintf(colname[index1], "u_%d ", i);
		index1++;
	}
	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, colname);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(colname);


	//Degree constraints

	numrows = 2*N;
	numnz = 2 * N*(N-1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (i = 0; i < N; i++){
		sense[index1] = 'E';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		for (j = 0; j < N; j++){
			if (i!=j){
				matind[index] = pos_y[i][j];
				matval[index++] = 1;
			}
		}
		sense[index1] = 'E';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		for (j = 0; j < N; j++){
			if (i != j){
				matind[index] = pos_y[j][i];
				matval[index++] = 1;
			}
		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//MTZ

	numrows = N*(N-1);
	numnz = 3*N*(N - 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");


	index = 0;
	index1 = 0;
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			if (i != j && i!=1  && j!=1){
				sense[index1] = 'L';
				rhs[index1] = N - 1;
				matbeg[index1++] = index;
				matind[index] = pos_u[i];
				matval[index++] = 1;
				matind[index] = pos_u[j];
				matval[index++] = -1;
				matind[index] = pos_y[i][j];
				matval[index++] = N;
			}
			}
		}
	

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	CPXwriteprob(env, lp, "model_TSP.lp", NULL);

	CPXmipopt(env, lp);

	CPXgetmipobjval(env, lp, &global.results->ovStep2);
	printf("TSP: %.2f   ", global.results->ovStep2);


	i = CPXgetstat(env, lp);
	if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 103)
		printf(" infeasible solution\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);



	numcols = CPXgetnumcols(env, lp);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(env, lp, x, 0, numcols - 1);




	printf("solucion 'y' que se encuentra en el tsp \n");

	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			if (i != j && x[pos_y[i][j]]>0){
				printf("%d %d %f %f\n", i, j, x[pos_y[i][j]],global.results->alpha_cur[index_e[i][j]]);
				global.results->yStep2[index_e[i][j]] = 1;
			}
		}
	}



	if (lp != NULL) {
		status = CPXfreeprob(env, &lp);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	free(pos_u);
	for (i = 0; i < N; i++){
		free(pos_y[i]);

	}
	free(pos_y);

	///for (i = 0; i < N; i++){
	//	free(coef[i]);

	//}
	//free(coef);


	


	return status;
}*/
////////////////////////////////////////////////////////////////////////////////


///Crear entorno para modelo del paso 1
////////////////////////////////////////////////////////////////////////////////
int create_CPLEX_Step1_enviroment(CPXENVptr *env, GLOBAL_INFO global) {
	int status = 0;

	// Create enviroment:
	*env = CPXopenCPLEX(&status);
	assert(*env != NULL);

	// Set parameters:
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_THREADS,    1));						// Threads usados
    status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_TILIM,		global.param->MAX_CPU_TIME));		// Time Limit
	status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_EPGAP,		EPSILON*EPSILON));			// Gap de Epsilon Optimalidad
    status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_CUTSFACTOR, 1.0));						// <= 1.0 No cuts will be generated, >1: Limits the number of cuts that can be added
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPSEARCH,  CPX_MIPSEARCH_TRADITIONAL));// Turn on traditional search for use with control callbacks 
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPCBREDLP, CPX_OFF));					// Let MIP callbacks work on the original model 
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_HEURFREQ,	-1));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPINTERVAL, 1));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PRESLVND, -1));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PREIND, 0));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_REPEATPRESOLVE, 0));
	//status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_CUTUP, (global.results->UpperBound)+0.01));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPEMPHASIS, 2));		//	0	CPX_MIPEMPHASIS_BALANCED		Balance optimality and feasibility; default
																				//	1	CPX_MIPEMPHASIS_FEASIBILITY		Emphasize feasibility over optimality
																				//	2	CPX_MIPEMPHASIS_OPTIMALITY		Emphasize optimality over feasibility
																				//	3	CPX_MIPEMPHASIS_BESTBOUND		Emphasize moving best bound
																				//	4	CPX_MIPEMPHASIS_HIDDENFEAS		Emphasize finding hidden feasible solutions	
	// The Lazy Constraints Callback will switch off these anyway:
    status = MAX(status, CPXsetintparam(*env, CPX_PARAM_REDUCE,		CPX_OFF));	// 0: No primal or dual reductions,   1: Only primal reductions ???
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PRELINEAR,	CPX_OFF));	// Assure linear mappings between the presolved and original models 
	

	//status = MAX(status, CPXsetintparam(*env, CPX_PARAM_VARSEL,  3)); // 0: default, 1: max_infeas, 2: pseudo_cost, 3: strong_branching
    


	return status;
}
////////////////////////////////////////////////////////////////////////////////



////Crear modelo paso 1
////////////////////////////////////////////////////////////////////////////////
int create_CPLEX_Step1_lp(CPXENVptr env, CPXLPptr *lp_ptr, GLOBAL_INFO global) {

	clock_t	start;
	int i, j, k, l, c, e, t;
	int added_cuts, status = 0;

	// Just for lazyness:
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	int		*index_i  = global.data->index_i;
	int		*index_j  = global.data->index_j;
	int		**index_k = global.data->index_k;
	double	**W = global.data->W;
	double	**C = global.data->C;
	double  **F = global.data->F;
	COMMODITY *Com = global.data->Com;

	// New data:
	//double **Y = create_double_matrix(N,N);

    /// CPLEX PROBLEM ////////////////////////////////////////////////////
    int       numcols; // número de variables ..................................
    int       numrows; // número de restricciones sin contar las cotas .........
    int       numnz;   // número de elementos no nulos de la matriz ............
    int       objsen;  // sentido de la optimizacion (min:1, max:-1 ) ..........
    double    *obj;    // coeficientes de la función objetivo ..................
    double    *rhs;    // términos independientes de las restricciones .........
    char      *sense;  // sentido de las restricciones (<=: 'L', =:'E', >=:'G'). 
    int       *matbeg; // índice del primer coeficiente no nulo de cada columna.
    int       *matcnt; // número de elementos no nulos de cada columna .........
    int       *matind; // fila a la que corresponde cada elemento no nulo  .....
    double    *matval; // valores de los coef no nulos en las restricciones ....
    double    *lb;     // cotas inferiores de las variables ....................
    double    *ub;     // cotas superiores de las variables ....................
	int       index, index1;
	int       *pos_alpha;
	double	  *solfact;
	int		  numedges;
	char	**colname;
    ////////////////////////////////////////////////////////////////////////////
	
	
	//objsen  = CPX_MAX;		// Dirección de la optimización (minimizar!)
	numcols = E ;		// Number of variables of the Master Prob (Y's & Z's)
	//numrows = N + 1;		// Number of initial Constrains of the Master Problem (single_node cutsets & sum y = N-1)
	//numnz   = N*(N-1) + E;	// Number of Non-Zero coefs (single_node cutsets & sum y = N-1)
	
	pos_alpha = create_int_vector(E);
	solfact = create_double_vector(E);
	
	colname = (char**)calloc(numcols, sizeof(char*));
	for (i = 0; i < numcols; i++){
		colname[i] = (char*)calloc(255, sizeof(char));
	}

	// Define the linear program
	*lp_ptr = CPXcreateprob(env, &status, "Step1_FENCHEL_PROBLEM");
	assert(*lp_ptr != NULL);
	CPXchgobjsen(env, *lp_ptr, CPX_MAX);


	//Define alpha_ij variables
	index1 = 0;  // index of columns
	numcols = E;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	for (e = 0; e<E; e++){
		pos_alpha[e] = index1;
		obj[index1] = 0;//global.results->relax_solution[e];
		lb[index1] = 0;//-CPX_INFBOUND;
		ub[index1] = CPX_INFBOUND;
		sprintf(colname[index1], "alpha_%d,%d ",index_i[e],index_j[e]);
		index1++;
	}
	status = CPXnewcols(env, *lp_ptr, index1, obj, lb, ub, NULL, colname);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);


		CPXwriteprob(env, *lp_ptr, "model2.lp", NULL);

	// Remember initial number of constraints:
	global.results->initial_cons = CPXgetnumrows(env, *lp_ptr);
	
	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int solve_Fenchel_procedure( GLOBAL_INFO global,CPXENVptr env2,CPXLPptr  lp2) {
	clock_t start;
	int		status = 0;
	int		cuts_found = 0;
	int     imp_it = 0;
	int		no_imp = 0;
	double	old_LowerBound;
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		Fcuts;
	double cota; //Solución en step 1. Si es menor que 1, entonces no hay corte de Fenchel
	int e;
	struct timespec time_now;
	double elapsedn=0;
	int *indice;
	double *coef;
	int e1,i;
	double *dist_aux=create_double_vector(E);
	FILE *output_prueba2 = NULL;
   FILE *output_prueba3 = NULL;
      FILE *output_prueba4 = NULL;
	int ol=0;
	
	/*SHRINKING----------------------------------------------------------------------------------------*/
	struct LINE *lines;
	int* degree=create_int_vector(N);
	int* visited=create_int_vector(N);
	int maxim=(int) round(N/2);
	int* in_line=create_int_vector(N);
	
	
	lines=create_line_vector(maxim);
	for(i=0;i<N;i++){
	lines[i].seq=create_int_vector(N-1);
	}
	int lcount=0;
	int count=0;
	int j=0;
	
	for(i=0;i<maxim;i++){
		for(j=0;j<N-1;j++){
			lines[i].seq[j]=0;
		}
	}
	
	
	/*vector degree*/
	for(e=0;e<E;e++){
	if(global.results->relax_solution[e]==1){
		degree[global.data->index_i[e]]+=1;
		degree[global.data->index_j[e]]+=1;
	}
	}
	/////////////////*/
	
	/*Comprobacion de vector degree*/
	for(e=0;e<E;e++){
	if(global.results->relax_solution[e]>0){
	printf("%d %d value %f\n",global.data->index_i[e],global.data->index_j[e],global.results->relax_solution[e]);
	}
	}
	
	for(i=0;i<N;i++){
	printf("node %d degree %d\n",i,degree[i]);
	}
	//////////////////////////////////*/
	
	
	
	///estructura lines///
	
	/*Quitar lo que son solo aristas y no caminos*/
	for(i=0;i<N;i++){
		if(degree[i]==1 && visited[i]==0){
			for(j=0;j<N;j++){
				if(global.results->relax_solution[global.data->index_e[i][j]]==1 && visited[i]==0){
					if(degree[j]==1){
						visited[i]=1;
						visited[j]=1;
						break;
					}
				}
			}
		}
	}
	
	/*Guardar caminos*/
	for(i=0;i<N;i++){
		printf("%d\n",i);
		if(degree[i]==1 && visited[i]==0){
			visited[i]=1;
			lines[lcount].dim+=1;
			lines[lcount].start=i;
			for(j=0;j<N;j++){
				if(global.results->relax_solution[global.data->index_e[i][j]]==1 && visited[j]==0){
					visited[j]=1;
					lines[lcount].dim+=1;
					lines[lcount].seq[count]=j;
					count+=1;
					if(degree[j]==2){
						i=j;
						j=-1;
					}
					if(degree[j]==1){
						i=-1;
						count=0;
						lcount+=1;
						break;
					}
				}
			}
		}
	}
	
	for(i=0;i<lcount;i++){
	printf("linea %d start %d\n",i,lines[i].start);
	}
	
	for(i=0;i<lcount;i++){
		printf("linea %d: \n",i); 
		printf("dimension %d\n",lines[i].dim);
		printf("start %d",lines[i].start);
	for(j=0;j<N-1;j++){
	printf(" %d",lines[i].seq[j]);
	}
	printf("\n");
	}
	
		
	for(i=0;i<N;i++)in_line[i]=-1;
	
	for(i=0;i<lcount;i++){
	in_line[lines[i].start]=i;
	for(j=0;j<lines[i].dim-1;j++){
	in_line[lines[i].seq[j]]=i;
	}
	}
	
	for(i=0;i<N;i++)printf("%d in line %d\n",i,in_line[i]);
	
	//getchar();
	
	
	/*****Pruebas para generar el reducido*****/
	int newN=lcount*2+N;
	for(i=0;i<lcount;i++){
		newN-=lines[i].dim;
	}	
	printf("%d\n",newN);
	//getchar();
    int* oldnod=create_int_vector(newN);
    double* edges=create_double_vector(newN*(newN-1)/2);
    int* index_ii;
    int* index_jj;
    int** index_ee;
    
    index_ii = create_int_vector(newN*(newN-1)/2);
	index_jj = create_int_vector(newN*(newN-1)/2);
	index_ee = create_int_matrix(newN, newN);

	e=0;
	//for (i=0; i<newN-1; i++) {
		//for (j=i+1; j<newN; j++) {
			//index_ii[e] = i;
			//index_jj[e] = j;
			//index_ee[i][j] = e;
			//index_ee[j][i] = e;
			//e++;
		//}
	//}
    
    
    j=0;
    for(i=0;i<lcount;i++){
		oldnod[j++]=lines[i].start;
		oldnod[j++]=lines[i].seq[lines[i].dim-2];
	}
	for(i=0;i<N;i++){
		if(in_line[i]==-1){
		oldnod[j++]=i;
		}
	}
    
    for(i=0;i<newN;i++){
	printf("%d %d\n",i,oldnod[i]);
	}
	//getchar();
	
	
	
	int newE;
	
	for(i=0;i<newN;i++){
		for(j=0;j<newN;j++){
			index_ee[i][j]=-1;
		}
	}
	
	//Arreglar//
	e=0;
	for(i=0;i<newN-1;i++){
		if(in_line[oldnod[i]]==-1){
			for(j=i+1;j<newN;j++){
				if(global.results->relax_solution[global.data->index_e[oldnod[i]][oldnod[j]]]>0.0001){
				index_ii[e] = i;
				index_jj[e] = j;
				index_ee[i][j] = e;
				index_ee[j][i] = e;
				edges[index_ee[i][j]]=global.results->relax_solution[global.data->index_e[oldnod[i]][oldnod[j]]];
				e++;
				}
			}
		}
	}
	
	for(i=0;i<newN-1;i++){
		if(in_line[oldnod[i]]!=-1){
			 for(j=i+1;j<newN;j++){
				if(in_line[oldnod[i]]==in_line[oldnod[j]]){
				index_ii[e] = i;
				index_jj[e] = j;
				index_ee[i][j] = e;
				index_ee[j][i] = e;
				edges[index_ee[i][j]]=1;
				e++;
				}
				else if(in_line[oldnod[i]]!=in_line[oldnod[j]] && global.results->relax_solution[global.data->index_e[oldnod[i]][oldnod[j]]]>0.0001){
				edges[index_ee[i][j]]=1;
				index_ii[e] = i;
				index_jj[e] = j;
				index_ee[i][j] = e;
				index_ee[j][i] = e;
				edges[index_ee[i][j]]=global.results->relax_solution[global.data->index_e[oldnod[i]][oldnod[j]]];
				e++;
				}
			 }
		}
	}
	
	newE=e;
	
	printf("PRUEBA\n");
	for(e=0;e<newE;e++){
		printf("%d %f %d %d %d %d\n",e,edges[e],index_ii[e],index_jj[e],oldnod[index_ii[e]],oldnod[index_jj[e]]);
	}
	
	printf("\n\n");
	for(i=0;i<newN;i++){
	for(j=0;j<newN;j++){
	printf("%d %d %d\n",index_ee[i][j],i,j);
	}
	}
	
	

	
	
	printf("FUNCION ");
     solve_STEP1(newN,newE,index_ii,index_jj,index_ee,edges);
   //  getchar();
	
		
    
    /*****************************************/
    ////////////////////////////////////////////////////////////////////////////////////////////////////
	
	indice=create_int_vector(E);
	coef=create_double_vector(E);
	for( e=0;e<E;e++) indice[e]=e;
	
	
	/////////////////////////////////////////////////////////////
	/////SECCION A MODIFICAR: IMPRIMIR PARA VER LA EVOLUCION//////
	
	output_prueba2 = open_file("soluciones.txt", "a");
	output_prueba3 = open_file("comprobacion.txt", "a");
	output_prueba4 = open_file("cortes_ejemplos.txt", "a");
	
	
	printf("Solucion inicial:\n");
	
	fprintf(output_prueba2,"Solucion inicial:\n");
	fprintf(output_prueba3,"Solucion inicial:\n");
		fprintf(output_prueba4,"Solucion inicial:\n");
		
	for(e=0;e<E;e++){
	if(global.results->relax_solution[e]>0)printf("y0_%d,%d=%f\n",global.data->index_i[e],global.data->index_j[e],global.results->relax_solution[e]); 
	}	
	
	for(e=0;e<E;e++){
	if(global.results->relax_solution[e]>0)fprintf(output_prueba2,"y0_%d,%d=%f\n",global.data->index_i[e],global.data->index_j[e],global.results->relax_solution[e]); 
	}
	
	for(e=0;e<E;e++){
	if(global.results->relax_solution[e]>0)fprintf(output_prueba3,"y0_%d,%d=%f\n",global.data->index_i[e],global.data->index_j[e],global.results->relax_solution[e]); 
	}
	for(e=0;e<E;e++){
	if(global.results->relax_solution[e]>0)fprintf(output_prueba4,"y0_%d,%d=%f\n",global.data->index_i[e],global.data->index_j[e],global.results->relax_solution[e]); 
	}
	
	for(e=0;e<E;e++){
	if(global.results->relax_solution[e]>0)printf("y0_%d,%d=%f\n",global.data->index_i[e],global.data->index_j[e],global.results->relax_solution[e]); 
	}
	///////////////////////////////////////////////////////////////////////////////
	
	global.results->yStep2 = create_double_vector(E);
	for ( e = 0; e < E; e++)global.results->yStep2[e] = 0;

	if (global.param->SCREEN_OUTPUT >= 2) { printf("\t\tNumRows: %d\n", CPXgetnumrows(env2, lp2)); }


	//Inicializamos la cota=2 para que entre inicialmente en el bucl
	cota=2;
	
	for( e=0;e<E;e++){
	  coef[e]=global.results->relax_solution[e];
	}
	
	////Cambiamos la funcion objetivo del paso 1
	CPXchgobj(env2,lp2,E,indice,coef);
	
	/////Eliminar restricciones no necesarias:////////////////////////////////////////////////////////////
	
	CPXwriteprob(env2, lp2, "model2_antes_de_quitar.lp", NULL);
	//printf("\n antiguo \n");
	//getchar();
     
     double coeff=0;
     double sum=0;
     int nozeros=0;
    
     for(e=0;e<E;e++){
		global.results->const_record[e]=0; //// const_record: para guardar si las variables aparecen en alguna restriccion del paso 1.
		if(global.results->relax_solution[e]>0) nozeros+=1; ///numero de aristas con coeficientes mayor que 0 en el paso inicial
	 }
	 
	 
	////En cada restriccion del paso 1, contamos el numero de variables que forman parte de la restriccion y que tienen coeficiente mayor
	////que 0 en la funcion objetivo (sum). Si es menor que el 75% de las variables que no son 0 en la funcion objetivo, entonces eliminamos la restriccion 
	for(i=0;i<CPXgetnumrows(env2,lp2);i++){
		sum=0;
		for(e=0;e<E;e++){
			status=CPXgetcoef(env2, lp2, i, e, &coeff);
			if(coeff>0 && global.results->relax_solution[e]>0) sum+=1;
			//sum+=coeff*global.results->relax_solution[e];
		}
		//printf("%d %f\n",i,sum);
		if(sum<round(0.75*nozeros)){
			 CPXdelrows(env2, lp2, i, i);
		 }
	}
	
	
	
	for(i=0;i<CPXgetnumrows(env2,lp2);i++){
		for(e=0;e<E;e++){
			status=CPXgetcoef(env2, lp2, i, e, &coeff);
			if(coeff>0) global.results->const_record[e]=1;
		}
	}
	
	
	CPXwriteprob(env2, lp2, "model2_tras_quitar.lp", NULL);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	////-----Add constraints to avoid unboundedness-----////
	
	//////////////////////////////////////////idea 0: <=1 cada arista que aparece en la fobj///////////////////////
	int       numcols; // número de variables ..................................
    int       numrows; // número de restricciones sin contar las cotas .........
    int       numnz;   // número de elementos no nulos de la matriz ............
    int       objsen;  // sentido de la optimizacion (min:1, max:-1 ) ..........
    double    *obj;    // coeficientes de la función objetivo ..................
    double    *rhs;    // términos independientes de las restricciones .........
    char      *sense;  // sentido de las restricciones (<=: 'L', =:'E', >=:'G'). 
    int       *matbeg; // índice del primer coeficiente no nulo de cada columna.
    int       *matcnt; // número de elementos no nulos de cada columna .........
    int       *matind; // fila a la que corresponde cada elemento no nulo  .....
    double    *matval; // valores de los coef no nulos en las restricciones ....
    double    *lb;     // cotas inferiores de las variables ....................
    double    *ub;     // cotas superiores de las variables ....................
	int       index, index1;
	int       *pos_alpha;
	double	  *solfact;
	int		  numedges;
	char	**colname;
		numedges = 0;
	for (e = 0; e < E; e++){
		if (global.results->relax_solution[e] > 0 && global.results->const_record[e]<0.5){
			numedges += 1;
		}
	}
	
	numrows = numedges;
	numnz = numedges;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (e = 0; e < E; e++){
		if (global.results->relax_solution[e]>0 && global.results->const_record[e]<0.5){
		global.results->const_record[e]=1;
		sense[index1] = 'L';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		matind[index] = e;
		matval[index++] = 1;
		}
	}
	

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	//getchar();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*///////////////////////////////////////////idea2///////TOUR INICIAL///////////////////////////////////////////////
	//for(e=0;e<E;e++){
		//dist_aux[e]=-coef[e];
	//} 
	//printf("distancias\n");
		//for(e=0;e<E;e++){
		////if(dist_aux[e]>0){
		//printf("%d %f\n",e,dist_aux[e]);
		//}
		//printf("-------------------------------");
		//getchar();
    
	//heuristic_tsp_Concorde(global,dist_aux); 		
		//add_cuts_to_Step1(global,env2,lp2); 	
			//CPXwriteprob(env2, lp2, "model2.lp", NULL);
			//printf("stop");
			//getchar();
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////idea1///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for(e=0;e<E;e++){
		//if(global.results->const_record[e]==0 && global.results->relax_solution[e]>0){
	     //////Creo unas distancias que hagan que esa arista este en la solucion del tsp. Uso el heuristico en lugar del exacto para resolver
	     //dist_aux[e]=0;
			//for(e1=0;e1<E;e1++){
				//if(e1!=e){
					//dist_aux[e1]=1.1;//1+rand() / ((double) RAND_MAX);
				//}
				//if(e1!=e && global.results->relax_solution[e1]>0){
					//dist_aux[e1]=global.results->relax_solution[e1];
				//}
				//printf("%d %f %f\n",e1,dist_aux[e1], global.results->relax_solution[e1]);
				////getchar();
			//}
		
		
		//printf("e------------------------- %d----------------------- %f",e,global.results->relax_solution[e] );
		//getchar();
		
		//printf("distancias\n");
		//for(e=0;e<E;e++){
		////if(dist_aux[e]>0){
		//printf("%d %f\n",e,dist_aux[e]);
		//}
		//printf("-------------------------------");
		//getchar();
		
		//////////LLamo a concorde con el heuristico y creo la restriccion. Al crear la restriccion debo contar
		 
		//heuristic_tsp_Concorde(global,dist_aux); 		
		//add_cuts_to_Step1(global,env2,lp2); 	
			//CPXwriteprob(env2, lp2, "model2.lp", NULL);
			//printf("stop");
			//getchar();
	//}
	//}*/
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////------------------------------------------------////*/
	CPXwriteprob(env2, lp2, "model2_conbasicas.lp", NULL);

	
	////While cota>1 (Es decir, mientras que el valor objetico optimo de la funcion objetivo sea mayor que 1.	
	int it2=0;
	while(cota>1.01){
		//QUITAR EL BREAK DE ABAJO
		break;
		clock_gettime(CLOCK_MONOTONIC, &time_now);
		
		elapsedn = (time_now.tv_sec - global.results->start_total.tv_sec);
	    elapsedn += (time_now.tv_nsec - global.results->start_total.tv_nsec) / 1000000000.0;
	    
		if(elapsedn > global.param-> MAX_CPU_TIME){
			Fcuts=0;
			break;
		}
		
	
		start = clock();
		CPXwriteprob(env2, lp2, "model2.lp", NULL);
		global.results->numb_rest=CPXgetnumrows(env2, lp2);
		status = CPXlpopt(env2, lp2);
		assert(status == 0);
		global.results->master_time += elapsed(start);
   

		status = CPXgetobjval(env2, lp2, &cota);


		status = CPXgetx(env2, lp2, global.results->alpha_cur, 0, CPXgetnumcols(env2, lp2) - 1);;


		it2+=1;// Para contar iteraciones
	
	
		fprintf(output_prueba3,"------------Iteracion %d---------------------------------- \n",it2);
		fprintf(output_prueba3,"Solucion de paso 1:\n");
		for(e=0;e<E;e++){
			if(global.results->alpha_cur[e]!=0){
			fprintf(output_prueba3,"alpha %d %d %f\n",global.data->index_i[e],global.data->index_j[e],global.results->alpha_cur[e]);
			}
		}
		fprintf(output_prueba2,"Valor objetivo: %f\n", cota);
	
		
		if (cota <= 1.01){
			Fcuts=0;
			break;
		}

		solve_tsp_Concorde(global);
		global.results->concorde+=1;
		
		printf("Valor objetivo de paso 2: %f Tiempo de concorde en esta iteracion: %f \n",global.results->ovStep2,global.results->concorde_time);
		

		
		if (global.results->ovStep2 <= 1.001){
			
			printf("Corte que hay que incluir en el problema principal\n");
			for(e=0;e<E;e++){
			if(global.results->alpha_cur[e]>0){
			printf("%f y_%d,%d+",global.results->alpha_cur[e],global.data->index_i[e],global.data->index_j[e]);
			}
			}
			printf("<=1\n\n\n");
		
			printf("Valor objetivo de paso 2: %f Tiempo de concorde en esta iteracion: %f \n",global.results->ovStep2,global.results->concorde_time);
			
			fprintf(output_prueba2,"Corte que hay que incluir en el problema principal\n");
			for(e=0;e<E;e++){
			if(global.results->alpha_cur[e]>0){
			fprintf(output_prueba2,"%f y_%d,%d+",global.results->alpha_cur[e],global.data->index_i[e],global.data->index_j[e]);
			}
			}
			fprintf(output_prueba2,"<=1\n\n\n");
			
			fprintf(output_prueba3,"Corte que hay que incluir en el problema principal\n");
			for(e=0;e<E;e++){
			if(global.results->alpha_cur[e]>0){
			fprintf(output_prueba3,"%f y_%d,%d+",global.results->alpha_cur[e],global.data->index_i[e],global.data->index_j[e]);
			}
			}
			fprintf(output_prueba3,"<=1\n\n\n");
			
			fprintf(output_prueba4,"Corte que hay que incluir en el problema principal\n");
			for(e=0;e<E;e++){
			if(global.results->alpha_cur[e]>0){
			fprintf(output_prueba4,"%f y_%d,%d+",global.results->alpha_cur[e],global.data->index_i[e],global.data->index_j[e]);
			}
			}
			fprintf(output_prueba4,"<=1\n\n\n");
			
			printf("Corte que hay que incluir en el problema principal\n");
			for(e=0;e<E;e++){
			if(global.results->alpha_cur[e]>0){
			printf("%f y_%d,%d+",global.results->alpha_cur[e],global.data->index_i[e],global.data->index_j[e]);
			}
			}
			printf("<=1\n\n\n");
			
			fprintf(output_prueba4,"Solucion inicial:\n");
	
				
			fprintf(output_prueba3,"Alguno con coeficiente 0 en el paso 1\n");
			for(e=0;e<E;e++){
			if(global.results->alpha_cur[e]>0 && global.results->relax_solution[e]<=0){
		    fprintf(output_prueba3,"i %d j%d alpha %f y0 %f\n",global.data->index_i[e],global.data->index_j[e],global.results->alpha_cur[e],global.results->relax_solution[e]);
		    }
		     }
			Fcuts=1;
			break;
		}
		else if (global.results->ovStep2 >1.001){
			
		fprintf(output_prueba3,"Res. paso 1 %f Res. paso 2 %f Cota inferior actual %f\n",cota,global.results->ovStep2,global.results->LowerBound);
		fprintf(output_prueba3,"Restriccion a incluir en paso 1:\n");
	    for(e=0;e<E;e++){
		if(global.results->yStep2[e]>0.5){
		fprintf(output_prueba3,"alpha_%d_%d +",global.data->index_i[e],global.data->index_j[e]);
		}
		}
			fprintf(output_prueba3,"<=1\n\n");	
			fprintf(output_prueba3,"Alguno con coeficiente 0 en el paso 1\n");
			for(e=0;e<E;e++){
			if(global.results->alpha_cur[e]>0 && global.results->relax_solution[e]<=0){
		    fprintf(output_prueba3,"i %d j%d alpha %f y0 %f\n",global.data->index_i[e],global.data->index_j[e],global.results->alpha_cur[e],global.results->relax_solution[e]);
		    }
		     }
			
			add_cuts_to_Step1(global, env2, lp2);
		}
		


	}


	return Fcuts;
}



///////////////////////////////////////////////////////////////////////////
/*void add_basic(GLOBAL_INFO global, CPXCENVptr env2, CPXLPptr lp2) {	

	int		e, k;
	int		E = global.data->E;
	int		K = global.data->K;
	int		status = 0;
	int		numnz;
	double    *rhs;    // términos independientes de las restricciones .........
    char      *sense;  // sentido de las restricciones (<=: 'L', =:'E', >=:'G'). 
    int       *matbeg; // índice del primer coeficiente no nulo de cada columna.
    int       *matcnt; // número de elementos no nulos de cada columna .........
    int       *matind; // fila a la que corresponde cada elemento no nulo  .....
    double    *matval; // valores de los coef no nulos en las restricciones ....
	//double  *violation = create_double_vector(K);
	double  threshold;
	int     n_rows;
	int		added_cuts = 0;
	int		numrows, index1, index;


	int numedges = 0;
	for (e = 0; e < E; e++){
		if (global.results->relax_solution[e] > 0){
			numedges += 1;
		}
	}
	
	numrows = numedges;
	numnz = numedges;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (e = 0; e < E; e++){
		if (global.results->relax_solution[e]>0){
		sense[index1] = 'L';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		matind[index] = e;
		matval[index++] = 1;
		
		}
	}

	status = CPXaddrows(env2,lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);




	//free(matind);
	//free(matval);


}*/


//Funcion para annadir restricciones al paso 1
void add_cuts_to_Step1(GLOBAL_INFO global, CPXCENVptr env2, CPXLPptr lp2) {	

	int		e, k;
	int		E = global.data->E;
	int		K = global.data->K;
	int		status = 0;
	int		numnz=0;
	double  *rhs = create_double_vector(1);
	char	*sense = create_char_vector(1);
	int		*matbeg = create_int_vector(1);
	int		*matind;// = create_int_vector(E);
	double	*matval;// = create_double_vector(E);
	double  threshold;
	int     n_rows;
	int		added_cuts = 0;
	int		numrows, index1, index;

	numnz = 0;
	for (e = 0; e < E; e++){
		if (global.results->yStep2[e] > 0){
			numnz += 1;
		}
	}

	matind = create_int_vector(numnz);
	matval = create_double_vector(numnz);

	numrows = 1;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	sense[0] = 'L';
	rhs[0] = 1;
	matbeg[0] = index;
	for (e = 0; e < E; e++){
			if (global.results->yStep2[e]>0){
			global.results->const_record[e]=1;
			matind[index] = e;
			matval[index++] = global.results->yStep2[e];
		}
	}

	n_rows = CPXgetnumrows(env2, lp2);

	status = CPXaddrows(env2, lp2, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
	printf("%d ", CPXgetnumrows(env2, lp2));


	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");

	free(matind);
	free(matval);
}


////Annadir un corte de Fenchel////
void add_Fenchel_cut(GLOBAL_INFO global, CPXCENVptr env, CPXLPptr lp) {

	int		e, k;
	int		E = global.data->E;
	int		K = global.data->K;
	int		status = 0;
	int		numnz = 0;
	double  *rhs = create_double_vector(1);
	char	*sense = create_char_vector(1);
	int		*matbeg = create_int_vector(1);
	int		*matind;// = create_int_vector(E);
	double	*matval;// = create_double_vector(E);
	//double  *violation = create_double_vector(K);
	double  threshold;
	int     n_rows;
	int		added_cuts = 0;
	int		numrows, index1, index;

	numnz = 0;
	for (e = 0; e < E; e++){
		if (global.results->alpha_cur[e] != 0){
			numnz += 1;
		}
	}

	matind = create_int_vector(numnz);
	matval = create_double_vector(numnz);


	numrows = 1;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	sense[0] = 'L';
	rhs[0] = 1;
	matbeg[0] = index;
	for (e = 0; e < E; e++){
		if (global.results->alpha_cur[e]!=0){
			matind[index] = e;
			matval[index++] = global.results->alpha_cur[e];
		}
	}

	n_rows = CPXgetnumrows(env, lp);

	status = CPXaddrows(env, lp, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);

	printf("%d ", CPXgetnumrows(env, lp));


	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
free(matind);
free(matval);

}



////Funcion antigua para resolver tsp/Formulacion con MTZ

/*int solve_tsp2(GLOBAL_INFO global) {
	int status = 0;
	int i, j, e;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	//double best_upper_bound, best_lower_bound;

	//Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *obj;    // objective function coefficients ..............................
	double    *rhs;    // right and side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zelo element ...................
	double    *matval; // coefficient values fo the non-zero elements of constraints....
	double    *lb;     // lower bounds of variables.....................................
	double    *ub;     // upper bounds of variables.....................................
	//	int       status;  // optimization status......................... .................
	double    *x;      // solution vector (double, even if the problem is integer) .....
	char probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	//double    value;   // objevtive value of solution ..................................
	int       num_z_var, num_x_var, whichmodel;
	int       **pos_y;
	int       *pos_u;
	//float	  **coef;
	////////////////////////////////////////
	// Just for lazyness:
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	int		*index_i = global.data->index_i;
	int		*index_j = global.data->index_j;
	int		**index_k = global.data->index_k;
	double	**W = global.data->W;
	double	**C = global.data->C;
	double  **F = global.data->F;
	double	**C1 = global.data->C1;
	double	**C2 = global.data->C2;
	COMMODITY *Com = global.data->Com;
	char	**colname;
	//////////////////////////////////////////
//	for ( e = 0; e < E; e++)global.results->yStep2[e] = 0;


	//Initialize CPLEX environment
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}

	CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, 2);

	// Create the problem in CPLEX 
	strcpy(probname, "TSP");
	lp = CPXcreateprob(env, &status, probname);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}


	pos_y = create_int_matrix(N, N);
	//coef = create_int_matrix(N, N);
	pos_u = create_int_vector(N);
	x = create_double_vector(N*N + N);
	//for (int i = 0; i < N-1; i++){
	//for (j = i+1; j < N; j++){
	//coef[i][j] = global.results->alpha_cur[index_e[i][j]];
	//coef[j][i] = global.results->alpha_cur[index_e[i][j]];
	//}
	//}


	index1 = 0;  // index of columns
	numcols = N*(N - 1);
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	colname = (char**)calloc(numcols, sizeof(char*));
	for (i = 0; i < numcols; i++){
		colname[i] = (char*)calloc(255, sizeof(char));
	}


	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			if (i != j){
				pos_y[i][j] = index1;
				obj[index1] = global.data->C[i][j]*100;
				lb[index1] = 0;
				ub[index1] = 1;
				ctype[index1] = 'B';
				sprintf(colname[index1], "y_%d,%d ", i, j);
				index1++;
			}
		}
	}

	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, colname);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(colname);

	//Define u_k variables
	index1 = 0;  // index of columns
	numcols = N;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	colname = (char**)calloc(numcols, sizeof(char*));
	for (i = 0; i < numcols; i++){
		colname[i] = (char*)calloc(255, sizeof(char));
	}

	for (i = 0; i<N; i++){
		pos_u[i] = N*(N - 1) + index1;
		obj[index1] = 0;
		lb[index1] = -CPX_INFBOUND;
		ub[index1] = CPX_INFBOUND;
		ctype[index1] = 'C';
		sprintf(colname[index1], "u_%d ", i);
		index1++;
	}
	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, colname);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(colname);
	

	//Degree constraints

	numrows = 2 * N;
	numnz = 2 * N*(N - 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (i = 0; i < N; i++){
		sense[index1] = 'E';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		for (j = 0; j < N; j++){
			if (i != j){
				matind[index] = pos_y[i][j];
				matval[index++] = 1;
			}
		}
		sense[index1] = 'E';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		for (j = 0; j < N; j++){
			if (i != j){
				matind[index] = pos_y[j][i];
				matval[index++] = 1;
			}
		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//MTZ

	numrows = N*(N - 1);
	numnz = 3 * N*(N - 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");


	index = 0;
	index1 = 0;
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			if (i != j && i != 1 && j != 1){
				sense[index1] = 'L';
				rhs[index1] = N - 1;
				matbeg[index1++] = index;
				matind[index] = pos_u[i];
				matval[index++] = 1;
				matind[index] = pos_u[j];
				matval[index++] = -1;
				matind[index] = pos_y[i][j];
				matval[index++] = N;
			}
		}
	}


	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	CPXwriteprob(env, lp, "model_TSP2.lp", NULL);

	CPXmipopt(env, lp);
	


	CPXgetmipobjval(env, lp, &global.results->ovStep2);
	printf("TSP: %.2f   ", global.results->ovStep2);



	i = CPXgetstat(env, lp);
	if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 103)
		printf(" infeasible solution\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);



	numcols = CPXgetnumcols(env, lp);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(env, lp, x, 0, numcols - 1);




	//printf("solucion 'y' que se encuentra en el tsp \n");

	//for (i = 0; i < N; i++){
		//for (j = 0; j < N; j++){
			//if (i != j && x[pos_y[i][j]]>0){
				//printf("%d %d %f %f\n", i, j, x[pos_y[i][j]], F[i][j]);
				//global.results->yStep2[index_e[i][j]] = 1;
			//}
		//}
	//}






	if (lp != NULL) {
		status = CPXfreeprob(env, &lp);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	free(pos_u);
	for (i = 0; i < N; i++){
		free(pos_y[i]);

	}
	free(pos_y);

	//for (i = 0; i < N; i++){
	//free(coef[i]);
	//}
	//free(coef);

	return status;
}*/





//void free_line_vector(int **ptr) {
	//if (*ptr == NULL) {
		//fprintf(stderr, "ERROR: Unable to free memory from int_vector\n");
	//} else {
		//free(*ptr);
		//*ptr = NULL;
	//}
//}



//void add_cuts_to_Step12(int *tour, env, lp) {	

	//int		e, k;
	///int		E = global.data->E;
	//int		K = global.data->K;
	//int		status = 0;
	//int		numnz=0;
	//double  *rhs = create_double_vector(1);
	//char	*sense = create_char_vector(1);
	//int		*matbeg = create_int_vector(1);
	//int		*matind;// = create_int_vector(E);
	//double	*matval;// = create_double_vector(E);
	//double  threshold;
	//int     n_rows;
	//int		added_cuts = 0;
	//int		numrows, index1, index;

	//numnz = 0;
	//for (e = 0; e < E; e++){
		//if (global.results->yStep2[e] > 0){
			//numnz += 1;
		//}
	//}

	//matind = create_int_vector(numnz);
	//matval = create_double_vector(numnz);

	//numrows = 1;
	//d_vector(&rhs, numrows, "open_cplex:2");
	//c_vector(&sense, numrows, "open_cplex:3");
	//i_vector(&matbeg, numrows, "open_cplex:4");
	//i_vector(&matind, numnz, "open_cplex:6");
	//d_vector(&matval, numnz, "open_cplex:7");

	//index = 0;
	//index1 = 0;
	//sense[0] = 'L';
	//rhs[0] = 1;
	//matbeg[0] = index;
	//for (e = 0; e < E; e++){
			//if (global.results->yStep2[e]>0){
			//global.results->const_record[e]=1;
			//matind[index] = e;
			//matval[index++] = global.results->yStep2[e];
		//}
	//}

	//n_rows = CPXgetnumrows(env2, lp2);

	//status = CPXaddrows(env2, lp2, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
	//printf("%d ", CPXgetnumrows(env2, lp2));


	//if (status)
		//fprintf(stderr, "CPXaddrows failed.\n");

	//free(matind);
	//free(matval);
//}


/*Concorde para subgrafo*/

int * solve_tsp_Concorde2(int newN,int newE, int** adj) {	
	//time//
	struct timespec time1,time2;
	clock_gettime(CLOCK_MONOTONIC, &time1);
	////////
	int objConcorde;
	int status= 0;
	int kk=0;int m=0;
	int e,j,i;
	int semente=rand();
	double ovStep2;
    double szeit, *mybnd, *mytimebound;
    int success, foundtour, hit_timebound = 0;
    int *in_tour = (int *) NULL;
    int *out_tour = (int *) NULL;
    CCrandstate rstate;
    char *probname = (char *) NULL;
    static int run_silently = 1;
    CCutil_sprand(semente, &rstate);
    mybnd = (double *) NULL;
   	mytimebound = (double *) NULL;
	CCdatagroup dat;
	double concorde_time;
   
    fill_dat(newN,&dat,adj);	

    out_tour = CC_SAFE_MALLOC (newE, int); 
	probname = CCtsp_problabel(" "); 
	    
	status=CCtsp_solve_dat (newN, &dat,in_tour,out_tour, NULL, &ovStep2, &success,&foundtour, probname,   mytimebound, &hit_timebound, run_silently, &rstate);
       
	printf("status %d",status);
	printf("success %d\n",success);
	printf("optval %f\n",ovStep2);
	printf("foundtour %d\n",foundtour);
        
	szeit = CCutil_zeit();
	free(probname);
	CCutil_freedatagroup (&dat);
	//free(adj);
		
	clock_gettime(CLOCK_MONOTONIC, &time2);
	concorde_time = (time2.tv_sec - time1.tv_sec);
	concorde_time+= (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;
	  
	return out_tour;
	//free(out_tour);

}



/* FUNCION DEFINIENDO LP PARA CORTES FENCHEL */

void solve_STEP1(int newN,int newE,int* index_ii,int* index_jj,int** index_ee,double* edges){
	int i,j,e;
	int index,index1;  // auxiliar indices to fill in the constraint matrix
	//Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *obj;    // objective function coefficients ..............................
	double    *rhs;    // right and side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zelo element ...................
	double    *matval; // coefficient values fo the non-zero elements of constraints....
	double    *lb;     // lower bounds of variables.....................................
	double    *ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	double    *x;      // solution vector (double, even if the problem is integer) .....
	char probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	int       num_z_var, num_x_var, whichmodel;
	int       *pos_alpha;
	char      **colname;
	double	  objSTEP1=2;
	double 	  objSTEP2;
	int 	  iter=0;
	int       concorde_calls=0;
	int		  **adj;
	int		  *tour;
	int        k;
	int       *cut;
	int       edgecut;
	FILE *out = NULL;
	
	out= open_file("some_info.txt", "a");
  //Initialize vector with the index of the edges appearing in the cut
	cut=create_int_vector(newE);
  //initialize vector tour;
  tour=create_int_vector(newN);

  //Initialize CPLEX environment
	env = CPXopenCPLEX (&status);
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}

    // Create the problem in CPLEX 
	strcpy(probname,"STEP1");
	lp = CPXcreateprob (env, &status, probname);
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not create LP. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}
	
	///"POR AQUI///
	index1 = 0;  // index of columns
	numcols = newE;
	
	pos_alpha = create_int_vector(newE);
	colname = (char**)calloc(numcols, sizeof(char*));
	for (i = 0; i < numcols; i++){
		colname[i] = (char*)calloc(255, sizeof(char));
	}
	
	//Sentido funcion objetivo
	CPXchgobjsen(env, lp, CPX_MAX);
	
	//Define variables
	d_vector(&obj,numcols,"open_cplex:1");
	d_vector(&lb,numcols,"open_cplex:8");
	d_vector(&ub,numcols,"open_cplex:9");
	c_vector(&ctype,numcols,"open_cplex:01");
	
	
	for (e = 0; e<newE; e++){
		pos_alpha[e] = index1;
		obj[index1] = edges[e];
		lb[index1] = 0;
		ub[index1] = CPX_INFBOUND;
		sprintf(colname[index1], "alpha_%d,%d ",index_ii[e],index_jj[e]);
		index1++;
	}
	
	
	status = CPXnewcols(env,lp, index1, obj, lb, ub, NULL, colname);
	
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	
	
	//Define constraints
	numrows = newE;
	numnz = newE;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (e = 0; e < newE; e++){
		sense[index1] = 'L';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		matind[index] = e;
		matval[index++] = 1;
	}
	

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	
	CPXwriteprob(env,lp, "MODELO_NUEVO.lp", NULL);
	
	//Solve the model and add constraints if necessary (Fenchel procedure)
	while(objSTEP1>1.01){
		
		
		status = CPXlpopt(env, lp);
		assert(status == 0);
		status = CPXgetobjval(env, lp, &objSTEP1);
	
	
		//printf("objSTEP1 %f\n",objSTEP1);
		//getchar();
		d_vector(&x,numcols,"open_cplex:0");
		
	    CPXgetx(env,lp,x,0, numcols-1);  // obtain the values of the decision variables
		
		
		iter+=1;// Para contar iteraciones
		
		if (objSTEP1 <= 1.01){
			break;
		}

		//initialize_matrix_concorde
		//Initialize distance matrix for Concorde
		adj=create_int_matrix(newN,newN);
		
		for(i=0;i<newN;i++){
			for(j=0;j<newN;j++){
				adj[i][j]=0;
			}
		}
		
		
		for(e=0;e<newE;e++){
			adj[index_ii[e]][index_jj[e]]=-(int) round(10000*x[pos_alpha[e]]);
			adj[index_jj[e]][index_ii[e]]=-(int) round(10000*x[pos_alpha[e]]);
		}
		
		
		//~ for(i=0;i<newN;i++){
			//~ for(j=0;j<newN;j++){
				//~ printf("%f ",-(double) adj[i][j]/10000);
			//~ }
			//~ printf("\n");
	    //~ }
	    //~ printf("\n");
		
		tour=solve_tsp_Concorde2(newN,newE,adj);
		concorde_calls+=1;

		//for(e=0;e<newN;e++){
		 //printf("%d\n",tour[e]);
		//}

		//Calcular valor del TSP//
		
		
		for(i=0;i<newN;i++){
			for(j=0;j<newN;j++){
				adj[i][j]=0;
			}
		}
		
		for(e=0;e<newE;e++){
			adj[index_ii[e]][index_jj[e]]=(int) round(10000*x[pos_alpha[e]]);
			adj[index_jj[e]][index_ii[e]]=(int) round(10000*x[pos_alpha[e]]);
		}
		
		
		objSTEP2=0;
		for(i=0;i<newN-1;i++){
			//printf("%d %d %f\n",tour[i],tour[i+1], (double) adj[tour[i]][tour[i+1]]/10000);
				objSTEP2+= (double) adj[tour[i]][tour[i+1]]/10000;
		}
		objSTEP2+= (double) adj[tour[0]][tour[newN-1]]/10000;
		
		k=0;
		for(i=0;i<newN-1;i++){
			if(index_ee[tour[i]][tour[i+1]]>-1){
				cut[k++]=index_ee[tour[i]][tour[i+1]];
			}
		}
		if(index_ee[tour[0]][tour[newN-1]]>-1){
			cut[k++]=index_ee[tour[0]][tour[newN-1]];
		}
		
		edgecut=k;
		
		//for(k=0;k<edgecut;k++){
			//printf("cut %d: %d %d %d\n",k,cut[k],index_ii[cut[k]],index_jj[cut[k]]);
		//}
		//getchar();
			
		//printf("STEP 2 %f\n",objSTEP2);
		//getchar();
		fprintf(out,"iteration %d STEP1 %f STEP2 %f\n",iter,objSTEP1,objSTEP2);
		
		
		if (objSTEP2 <= 1.001){			
			printf("Corte que hay que incluir en el problema principal\n");
			fprintf(out, "Fenchel cut to add in the main model\n");
			for(e=0;e<newE;e++){
				if(x[pos_alpha[e]]>0){
				fprintf(out,"+%f *alpha_{%d,%d} ",x[pos_alpha[e]],index_ii[e],index_jj[e]);
			}
			}
			fprintf(out,"<=1\n");
			fclose(out);
			//Tengo que hacer el codigo para transformar el grafo y meter cortes en el problema original
			getchar();
			break;
		}
		else{		
			////add_cuts_to_Step1(tour, env, lp);
			int		numnz=0;
			double  *rhs = create_double_vector(1);
			char	*sense = create_char_vector(1);
			int		*matbeg = create_int_vector(1);
			int		*matind;// = create_int_vector(E);
			double	*matval;// = create_double_vector(E);
			////double  threshold;
			int     n_rows;
			int		numrows, index1, index;
			
			numnz=edgecut;

			matind = create_int_vector(numnz);
			matval = create_double_vector(numnz);

			numrows = 1;
			d_vector(&rhs, numrows, "open_cplex:2");
			c_vector(&sense, numrows, "open_cplex:3");
			i_vector(&matbeg, numrows, "open_cplex:4");
			i_vector(&matind, numnz, "open_cplex:6");
			d_vector(&matval, numnz, "open_cplex:7");

			index = 0;
			index1 = 0;
			sense[0] = 'L';
			rhs[0] = 1;
			matbeg[0] = index;
			for (k = 0; k < edgecut; k++){
				matind[index] = pos_alpha[cut[k]];
				matval[index++] = 1;
			}

			n_rows = CPXgetnumrows(env, lp);

			status = CPXaddrows(env, lp, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
			CPXwriteprob(env,lp, "MODELO_NUEVO.lp", NULL);
			////printf("%d ", CPXgetnumrows(env2, lp2));


			if (status)
				fprintf(stderr, "CPXaddrows failed.\n");
			free(matind);
			free(matval);
			}
		
		
		
	}
	///////////////
	if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
        fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
    }
	if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);
      if ( status ) {
        char  errmsg[1024];
        fprintf (stderr, "Could not close CPLEX environment.\n");
        CPXgeterrorstring (env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
      }
    }
   //free memory
   free(cut);
   free(tour);
   for(i=0;i<numcols;i++) free(colname[i]);
   free(colname);
   free(x);
   for(i=0;i<newN;i++)free(adj[i]);
   free(adj);
    
    
	
}


