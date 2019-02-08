#include "open_mp_functions.h"

/* Get arguments */
void get_args( int argc, char** argv, char** filename, int *rows, int *columns, char **type, int *selection, int *iterations, int *threads_num) {

	int i, j;

	if (argc == 15) {

    	for (i = 1; i < argc-1; i += 2) {

			// Check if argument is given again
			for( j = i+2; j < argc-1; j += 2 ) {

				if( !strcmp(argv[i], argv[j]) ) {
					printf("Error: Argument given again.\nPlease check README for more info.\n");
					MPI_Abort(MPI_COMM_WORLD, -1);
			        MPI_Finalize();
					exit(-1);
				}
			}

			// Get arguments
			if( !strcmp(argv[i], "-f") ) {
				*filename = malloc(strlen(argv[i+1])+1);
				strcpy(*filename, argv[i+1]);
			} else if( !strcmp(argv[i], "-r") ) {
				*rows = atoi(argv[i+1]);
			} else if( !strcmp(argv[i], "-c") ) {
				*columns = atoi(argv[i+1]);
			} else if( !strcmp(argv[i], "-t") ) {
	
				*type = malloc(strlen(argv[i+1])+1);
				strcpy(*type, argv[i+1]);

				if( strcmp(*type,RGB) && strcmp(*type,GREY) ) {
					printf("Error: Type must be 'rgb' or 'grey'.\nPlease check README for more info.\n");
					MPI_Abort(MPI_COMM_WORLD, -1);
					MPI_Finalize();
					exit(-1);
				}

			} else if( !strcmp(argv[i], "-s") ) {
				*selection = atoi(argv[i+1]);
				
				if( *selection < 1 || *selection > 10 ) {
					printf("Error: Selection must be a number from 1-10.\nPlease check README for more info.\n");
					MPI_Abort(MPI_COMM_WORLD, -1);
					MPI_Finalize();
					exit(-1);
				}

			} else if( !strcmp(argv[i], "-i") ) {
				*iterations = atoi(argv[i+1]);
			} else if( !strcmp(argv[i], "-d") ){
				*threads_num = atoi(argv[i+1]);
			} else {
				printf("Error: Wrong argument given.\nPlease check README for more info.\n");
				MPI_Abort(MPI_COMM_WORLD, -1);
				MPI_Finalize();
				exit(-1);
			}
		}

	} else {
		printf("Error: Not all arguments given.\nPlease check README for more info.\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
   		MPI_Finalize();
		exit(-1);
	}

	printf("Running with:\nfilename = %s\nRows = %d\nColumns = %d\nType = %s\nSelection = %d\nIterations = %d\nThreads = %d\n",
	*filename, *rows, *columns, *type, *selection, *iterations, *threads_num);
}

/* Return chosen filter */
 float** filter_normalizer(int selection) {

	float **h = NULL;
	
    switch(selection) {

		case 1:
			h = chosen_filter_normalizer(Blur);
			break;
		case 2:
			h = chosen_filter_normalizer(Sharpen);
			break;
		case 3:
			h = chosen_filter_normalizer(Edge);
			break;
		case 4:
			h = chosen_filter_normalizer(EdgeH);
			break;
		case 5:
			h = chosen_filter_normalizer(EdgeV);
			break;
		case 6:
			h = chosen_filter_normalizer(GradientH);
			break;
		case 7:
			h = chosen_filter_normalizer(GradientV);
			break;
		case 8:
			h = chosen_filter_normalizer(SobelH);
			break;
		case 9:
			h = chosen_filter_normalizer(SobelV);
			break;
		case 10:
			h = chosen_filter_normalizer(Emboss);
			break;

	}

	return h;
}

/* Return chosen filter normalized */
 float ** chosen_filter_normalizer(const int filter[S+2][S+2]) {

	int i, j;
	float **h = NULL, sum = 0.0;

	h = malloc(sizeof(float *)*(S+2));

	for( i = 0; i < (S+2); i++)
		h[i] = malloc(sizeof(float)*(S+2));
	
	for( i = 0; i < (S+2); i++) {
		for( j = 0; j < (S+2); j++) {
			h[i][j] = filter[i][j];
			sum += h[i][j];
		}
	}

	for( i = 0; i < (S+2); i++){
		for( j = 0; j < (S+2); j++) {
			if( sum != 0 )
				h[i][j] = h[i][j]/sum;
		}
	}

    return h;
}

/* Free chosen filter */
 void purge_filter(float ** filter) {

	int i;

    for( i = 0; i < (S+2); i++)
		free(filter[i]);

	free(filter);
}

/* Initialization of matrix */
 void init_matrix(unsigned char* matrix, unsigned char* new_matrix, unsigned char* read_buffer, int p_rows, int p_cols, int pix, int threads_num){
	int i, j;
#	pragma omp parallel for num_threads(threads_num) shared(p_rows,p_cols,matrix,pix,new_matrix,read_buffer) private(i,j) collapse(2)
	for( i = 0; i < (p_rows+2); i++) {
		for( j = 0; j < (p_cols+2*pix); j++) {
			if(i==0 || (i == p_rows + 1) || ((j>=0) && (j < pix)) || ( (j >= (p_cols + 1*pix)) && (j < (p_cols + 2*pix))) ){
				matrix[i*(p_cols+2*pix) + j]=(unsigned char) 0;
				new_matrix[i*(p_cols+2*pix) + j]=(unsigned char) 0;
			}else{
					matrix[i*(p_cols+2*pix) + j]=read_buffer[(i-1)*p_cols+(j-1*pix)];
			}
		}
	}
}

/* Initialization of neighbours */
 void init_neighbours(MPI_Comm new_comm, int * dim_size, int my_rank, int ndims, int* North, int* South, int* East, int* West, int* n_w, int* n_e, int* s_w, int* s_e){
	int coords[D], nghbr_coords[D];

	//find coordinates of the process
	MPI_Cart_coords(new_comm, my_rank, ndims, coords);
	//North-South
	MPI_Cart_shift( new_comm, 0, 1, North, South );
	//printf( "source and destination of rank %d for non-periodic shift 1 and direction 0 is %d and %d\n", my_rank, North , South);

	//West-East
	MPI_Cart_shift( new_comm, 1, 1, West, East );
	//printf( "source and destination of rank %d for non-periodic shift 1 and direction 1 is %d and %d\n", my_rank, West , East);

	//arxikopoihsh twn metavlitwn prosanatolismou
	*n_w=MPI_PROC_NULL;  *n_e=MPI_PROC_NULL;  *s_w=MPI_PROC_NULL;  *s_e=MPI_PROC_NULL;
	
	//North-West
	nghbr_coords[0] = coords[0]-1;
	nghbr_coords[1] = coords[1]-1;
	if(!((nghbr_coords[0] > (dim_size[0]-1) || nghbr_coords[0] < 0 )|| (nghbr_coords[1] < 0 || nghbr_coords[1] > (dim_size[1]-1))))
		MPI_Cart_rank(new_comm, nghbr_coords, n_w);
	
	//North-East
	nghbr_coords[0] = coords[0]-1;
	nghbr_coords[1] = coords[1]+1;
	if(!((nghbr_coords[0] > (dim_size[0]-1) || nghbr_coords[0] < 0 )|| (nghbr_coords[1] < 0 || nghbr_coords[1] > (dim_size[1]-1))))
		MPI_Cart_rank(new_comm, nghbr_coords, n_e);
	
	//South-West
	nghbr_coords[0] = coords[0]+1;
	nghbr_coords[1] = coords[1]-1;
	if(!((nghbr_coords[0] > (dim_size[0]-1) || nghbr_coords[0] < 0 )|| (nghbr_coords[1] < 0 || nghbr_coords[1] > (dim_size[1]-1))))
		MPI_Cart_rank(new_comm, nghbr_coords, s_w);
	
	//South-East
	nghbr_coords[0] = coords[0]+1;
	nghbr_coords[1] = coords[1]+1;
	if(!((nghbr_coords[0] > (dim_size[0]-1) || nghbr_coords[0] < 0 )|| (nghbr_coords[1] < 0 || nghbr_coords[1] > (dim_size[1]-1))))
		MPI_Cart_rank(new_comm, nghbr_coords, s_e);
}

/* Set output filename */
void set_output_filename(char** output_filename, char* type, int rows, int columns, int selection) {

	char filter_name[20];

	*output_filename = malloc(MAX_NUM_C*sizeof(char));

	if(*output_filename == NULL) {
		printf("Error: Wrong argument given.\nPlease check README for more info.\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
		MPI_Finalize();
		exit(-1);
	}

	switch(selection) {
		case 1:
			strcpy(filter_name,"Blur");
			break;
		case 2:
			strcpy(filter_name,"Sharpen");
			break;
		case 3:
			strcpy(filter_name,"Edge");
			break;
		case 4:
			strcpy(filter_name,"EdgeH");
			break;
		case 5:
			strcpy(filter_name,"EdgeV");
			break;
		case 6:
			strcpy(filter_name,"GradientH");
			break;
		case 7:
			strcpy(filter_name,"GradientV");
			break;
		case 8:
			strcpy(filter_name,"SobelH");
			break;
		case 9:
			strcpy(filter_name,"SobelV");
			break;
		case 10:
			strcpy(filter_name,"Emboss");
			break;
	}

	sprintf(*output_filename,"%s_image_%s_%d_%d.raw", filter_name, type, rows, columns);

	printf("The name of output file is %s\n", *output_filename);
}