#include "open_mp_functions.h"
 
int main(int argc, char *argv[]){

	int i, j;
	float **h = NULL;
	unsigned char *matrix=NULL, *new_matrix=NULL, *temp=NULL;
	unsigned char *read_buffer = NULL, *buffer = NULL;

	int my_rank, comm_sz;
	int threads_num;
	int pix;
	int rows, columns, selection, iterations;
  	int file_open_error, num_of_bytes, count;
  	int in_size, out_size, type_size;
	char *filename = NULL, *type = NULL, *output_filename = NULL;

	MPI_Status status;
    MPI_Datatype subarray;

	MPI_Offset my_offset;
	MPI_Offset total_number_of_bytes;
  	MPI_File fh, fh2;

	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if(my_rank==0){
		// Get args of program
		get_args(argc, argv, &filename, &rows, &columns, &type, &selection, &iterations, &threads_num);

		set_output_filename(&output_filename, type, rows, columns, selection);

		in_size=strlen(filename);
		out_size=strlen(output_filename);
		type_size=strlen(type); 
	}

	MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&columns, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&selection, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&threads_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&in_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&out_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&type_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(my_rank!=0){

		filename=malloc((in_size + 1)*sizeof(char));
		output_filename=malloc((out_size + 1)*sizeof(char));
		type=malloc((type_size + 1)*sizeof(char));
		if(filename == NULL || output_filename == NULL || type == NULL){
			printf("Memory allocation failed\n");
			MPI_Finalize();
			exit(-1);
		}
	}

	MPI_Bcast(filename, in_size + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(output_filename, out_size + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(type, type_size + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

	
	//open file with communicator MPI_COMM_WORLD
	file_open_error = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
 
	// https://www.hpc.ntnu.no/display/hpc/Reading+from+MPI+Files
  	if (file_open_error != MPI_SUCCESS) {
    	MPI_Abort(MPI_COMM_WORLD, file_open_error);
  	}

  	MPI_File_get_size(fh, &total_number_of_bytes);
  	//printf("%3d: total_number_of_bytes = %lld\n", my_rank, total_number_of_bytes);


  	if( !strcmp(type,RGB) )
		pix = 3;
	else
		pix = 1;

/****************For multiple of image********************/

  	if((rows*columns*pix)/total_number_of_bytes >= 2) {

  		int div_num = total_number_of_bytes/comm_sz;

		buffer = malloc(div_num*sizeof(unsigned char));

		if(buffer == NULL){
			printf("Memory allocation failed\n");
			MPI_Finalize();
			exit(-1);
		}

		MPI_File_seek(fh, my_rank * div_num, MPI_SEEK_SET);

		MPI_File_read_all(fh, buffer, div_num, MPI_BYTE, &status);

		MPI_File_close(&fh);


		//reopen with the same file header for reusability
		file_open_error = MPI_File_open(MPI_COMM_WORLD, output_filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh2);
 
  		if (file_open_error != MPI_SUCCESS) {
    		MPI_Abort(MPI_COMM_WORLD, file_open_error);
  		}

  		/*
		int k = (rows*columns*pix)/total_number_of_bytes;
		//printf("k = %d total = %d div num = %d\n", k, total_number_of_bytes, div_num);
		for( i = 0; i < div_num; i+=pix ) {

			for( j = 0; j < k; j++) {
				//printf("%d %d\n", i ,j);
				//printf("i*k+j = %d \n", i*k+j*pix);
				MPI_File_seek(fh2, k*my_rank*div_num+i*k+j*pix, MPI_SEEK_SET);
				MPI_File_write_all(fh2, &buffer[i], pix, MPI_UNSIGNED_CHAR, &status);
			}
		}
		*/

		for (i = 0 ; i < ((rows*columns*pix)/total_number_of_bytes) ; i++) {
			MPI_File_seek(fh2, my_rank*div_num+i*total_number_of_bytes, MPI_SEEK_SET);
			MPI_File_write_all(fh2, buffer, div_num, MPI_UNSIGNED_CHAR, &status);
		}

		MPI_File_close(&fh2);

		file_open_error = MPI_File_open(MPI_COMM_WORLD, output_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
 
  		if (file_open_error != MPI_SUCCESS) {
    		MPI_Abort(MPI_COMM_WORLD, file_open_error);
  		}

  		free(buffer);
  	}

/*******************************************************/


  	//check if rows or columns are divided by sqrt_comm_sz
	int sqrt_comm_sz = sqrt(comm_sz);
	if( rows%sqrt_comm_sz != 0 || columns%sqrt_comm_sz != 0 ) {
		printf("Error: rows or columns cannot be divided by comm_sz\n");
		MPI_Finalize();
		exit(-1);
	}


	// p_rows = process rows, p_cols = process columns
	int p_rows = rows/sqrt_comm_sz;
	int p_cols = pix*(columns/sqrt_comm_sz);

	// p_cols with extra cols
	int total_p_cols = p_cols+2*pix;

	num_of_bytes = p_rows * p_cols * sizeof(unsigned char);

	//create a subarray type for reading
	int ndims = D;
	int size[D] = {rows,pix*columns};
    int subsize[D] = {p_rows, p_cols};
    int start[D] = {0,0};

	MPI_Type_create_subarray(ndims, size, subsize, start, MPI_ORDER_C, MPI_CHAR, &subarray); 
	MPI_Type_commit(&subarray);

	// find box offset of current process
	//mporoume na xrisimopoihsoume kai cartesian topologies gia kalitero typo ama einai
	my_offset = (MPI_Offset) (my_rank/sqrt_comm_sz) * p_rows * p_cols * sqrt_comm_sz + (my_rank%sqrt_comm_sz) * p_cols;

	MPI_File_set_view(fh, my_offset, MPI_CHAR, subarray, "native", MPI_INFO_NULL);

	matrix = malloc((p_rows+2)*(total_p_cols)*sizeof(unsigned char));
	new_matrix = malloc((p_rows+2)*(total_p_cols)*sizeof(unsigned char));
	read_buffer = malloc(num_of_bytes);

	if(matrix == NULL || new_matrix == NULL || read_buffer == NULL){
		printf("Memory allocation failed\n");
		MPI_Finalize();
		exit(-1);
	}


	MPI_File_read_all(fh, read_buffer, num_of_bytes, MPI_UNSIGNED_CHAR, &status);
	//arxikopoihsh me mia stathera to perifragma kai gemisma tou eswterikou matrix
	init_matrix(matrix, new_matrix, read_buffer, p_rows, p_cols, pix, threads_num);
	
	double local_start_time, local_end_time, local_elapsed, elapsed;


	/**************************************************internal code*************************************************/

	
	int North, South, East, West, n_w, n_e, s_w, s_e;

	MPI_Request request_isend_n, request_irecv_n, request_isend_s, request_irecv_s,request_isend_w, request_irecv_w, request_isend_e, request_irecv_e,
    request_isend_nw, request_irecv_nw, request_isend_ne, request_irecv_ne, request_isend_sw, request_irecv_sw, request_isend_se, request_irecv_se;

	// Get filter
    h = filter_normalizer(selection);

	//Cartesian-topologies
	MPI_Comm old_comm, new_comm;
	int reorder, periods[D], dim_size[D];

	old_comm = MPI_COMM_WORLD;
	ndims = D; 							// 2D matrix/grid

	if(comm_sz < 4){
		dim_size[0] = 1; 				// rows
		dim_size[1] = comm_sz; 			// columns
	}else{
		dim_size[0] = sqrt_comm_sz; 	// rows
		dim_size[1] = sqrt_comm_sz; 	// columns
	}

	periods[0] = 0; 					// row non-periodic (each column forms a ring)
	periods[1] = 0; 					// columns non-periodic
	reorder = 1; 						// allows processes reordered for efficiency

	MPI_Cart_create(old_comm, ndims, dim_size, periods, reorder, &new_comm); 	// creates communicator


	//data-type creation
	MPI_Datatype neighbor_row, neighbor_col;

	MPI_Type_contiguous(p_cols, MPI_UNSIGNED_CHAR, &neighbor_row);
	MPI_Type_commit(&neighbor_row);

	//p_rows stoixeia kai p_cols + 2 stoixeia pou kanei skip
	MPI_Type_vector(p_rows, 1*pix, p_cols + 2*pix, MPI_UNSIGNED_CHAR, &neighbor_col);
	MPI_Type_commit(&neighbor_col);

	//find neighbours with topologies
	init_neighbours(new_comm, dim_size, my_rank, ndims, &North, &South, &East, &West, &n_w, &n_e, &s_w, &s_e);

	//printf("N=%d S=%d W=%d E=%d NE=%d NW=%d SE=%d SW=%d my_rank=%d\n", North, South, West, East, n_e, n_w, s_e, s_w, my_rank);

	int counter=0, ret;
	int global_flag=0, local_flag=0;

	MPI_Barrier(new_comm);
	//Start Time counter
	local_start_time = MPI_Wtime();
	
	//while( 1 ){
	while( counter < iterations ) {

		MPI_Isend(&matrix[(total_p_cols) + pix], 1, neighbor_col, West, 0, new_comm, &request_isend_w);
		MPI_Isend(&matrix[(total_p_cols) + p_cols ], 1, neighbor_col, East, 0, new_comm, &request_isend_e);
		MPI_Isend(&matrix[(total_p_cols) + pix], 1, neighbor_row, North, 0, new_comm, &request_isend_n);
		MPI_Isend(&matrix[p_rows*(total_p_cols) + pix], 1, neighbor_row, South, 0, new_comm, &request_isend_s);

		MPI_Isend(&matrix[(total_p_cols) + pix], pix, MPI_UNSIGNED_CHAR, n_w, 0, new_comm, &request_isend_nw);
		MPI_Isend(&matrix[(total_p_cols) + p_cols], pix, MPI_UNSIGNED_CHAR, n_e, 0, new_comm, &request_isend_ne);
		MPI_Isend(&matrix[p_rows*(total_p_cols) + pix], pix, MPI_UNSIGNED_CHAR, s_w, 0, new_comm, &request_isend_sw);
		MPI_Isend(&matrix[p_rows*(total_p_cols) + p_cols], pix, MPI_UNSIGNED_CHAR, s_e, 0, new_comm, &request_isend_se);
		
		MPI_Irecv(&matrix[0], pix, MPI_UNSIGNED_CHAR, n_w, MPI_ANY_TAG, new_comm, &request_irecv_nw);
		MPI_Irecv(&matrix[(total_p_cols)-pix], pix, MPI_UNSIGNED_CHAR, n_e, MPI_ANY_TAG, new_comm, &request_irecv_ne);
		MPI_Irecv(&matrix[(p_rows+1)*(total_p_cols)], pix, MPI_UNSIGNED_CHAR, s_w, MPI_ANY_TAG, new_comm, &request_irecv_sw);
		MPI_Irecv(&matrix[(p_rows+1)*(total_p_cols) + p_cols + pix], pix, MPI_UNSIGNED_CHAR, s_e, MPI_ANY_TAG, new_comm, &request_irecv_se);

		MPI_Irecv(&matrix[(total_p_cols)], 1, neighbor_col, West, MPI_ANY_TAG, new_comm, &request_irecv_w);
		MPI_Irecv(&matrix[(total_p_cols) + p_cols + pix], 1, neighbor_col, East, MPI_ANY_TAG, new_comm, &request_irecv_e);
		MPI_Irecv(&matrix[pix], 1, neighbor_row, North, MPI_ANY_TAG, new_comm, &request_irecv_n);
		MPI_Irecv(&matrix[(p_rows+1)*(total_p_cols) + pix], 1, neighbor_row, South, MPI_ANY_TAG, new_comm, &request_irecv_s);

		//i=2 , j=2 , i < (p_rows+2-2) , j < (p_cols+2 -2)
		

#	pragma omp parallel for num_threads(threads_num) shared(p_rows,p_cols,matrix,h,pix,new_matrix,total_p_cols) private(i,j) collapse(2)
		for( i = 2; i < p_rows; i++) {
			for( j = 2*pix; j < p_cols; j+=1) {
					/*for( p = -S; p <= S; p++) {
						for( q = -S*pix; q <= S*pix; q+=pix ) {
							total_f += (matrix[(i+p)*(total_p_cols) + (j+rgb+q)] * h[p+1][q/pix+1]);
						}
					}*/
				new_matrix[i*(total_p_cols) + j]= (unsigned char) ( (matrix[(i-1)*(total_p_cols) + (j-pix)] * h[S-1][S-1]) +
																	(matrix[(i-1)*(total_p_cols) + j] * h[S-1][S]) +
																	(matrix[(i-1)*(total_p_cols) + (j+pix)] * h[S-1][S+1]) +
																	(matrix[(i)*(total_p_cols) + (j-pix)] * h[S][S-1]) +
																	(matrix[(i)*(total_p_cols) + j] * h[S][S]) +
																	(matrix[(i)*(total_p_cols) + (j+pix)] * h[S][S+1]) +
																	(matrix[(i+1)*(total_p_cols) + (j-pix)] * h[S+1][S-1]) +
																	(matrix[(i+1)*(total_p_cols) + j] * h[S+1][S]) +
																	(matrix[(i+1)*(total_p_cols) + (j+pix)] * h[S+1][S+1])
																  );
			}
	    }
    	MPI_Wait(&request_irecv_n, &status);
   		MPI_Wait(&request_irecv_s, &status);
   		MPI_Wait(&request_irecv_e, &status);
   		MPI_Wait(&request_irecv_w, &status);
   		MPI_Wait(&request_irecv_nw, &status);
   		MPI_Wait(&request_irecv_ne, &status);
   		MPI_Wait(&request_irecv_sw, &status);
   		MPI_Wait(&request_irecv_se, &status);


   		//p_rows+2-1, p_cols+2-1
#	pragma omp parallel for num_threads(threads_num) shared(p_rows,p_cols,matrix,h,pix,new_matrix,total_p_cols) private(i,j) collapse(2)
   		for( i = 1; i < (p_rows+1); i++) {
			for( j = pix; j < (p_cols+pix); j+=1) {
				if(i==1 || (i == p_rows ) || (( j >= pix) && (j < (2*pix))) || ( (j >= p_cols) && (j < (p_cols + pix))) ){
					/*if(i==1 || i == p_rows || j == 1 || j == p_cols){
						for( p = -S; p <= S; p++) {
							for( q = -S*pix; q <= S*pix; q+=pix ) {
								total_f += (matrix[(i+p)*(total_p_cols) + (j+rgb+q)] * h[p+1][q/pix+1]);
							}
						}
					*/
					new_matrix[i*(total_p_cols) + j]= (unsigned char) ( (matrix[(i-1)*(total_p_cols) + (j-pix)] * h[S-1][S-1]) +
																		(matrix[(i-1)*(total_p_cols) + j] * h[S-1][S]) +
																		(matrix[(i-1)*(total_p_cols) + (j+pix)] * h[S-1][S+1]) +
																		(matrix[(i)*(total_p_cols) + (j-pix)] * h[S][S-1]) +
																		(matrix[(i)*(total_p_cols) + j] * h[S][S]) +
																		(matrix[(i)*(total_p_cols) + (j+pix)] * h[S][S+1]) +
																		(matrix[(i+1)*(total_p_cols) + (j-pix)] * h[S+1][S-1]) +
																		(matrix[(i+1)*(total_p_cols) + j] * h[S+1][S]) +
																		(matrix[(i+1)*(total_p_cols) + (j+pix)] * h[S+1][S+1])
																  	  );
				}
			}
		}

   		MPI_Wait(&request_isend_n, &status);
   		MPI_Wait(&request_isend_s, &status);
   		MPI_Wait(&request_isend_e, &status);
   		MPI_Wait(&request_isend_w, &status);
   		MPI_Wait(&request_isend_nw, &status);
   		MPI_Wait(&request_isend_ne, &status);
   		MPI_Wait(&request_isend_sw, &status);
   		MPI_Wait(&request_isend_se, &status);
   	
   		//sigrisi to matrix me ton new_matrix
   		ret = memcmp(matrix, new_matrix, (p_rows+2)*(total_p_cols));

		if(ret == 0) {
			local_flag = 1;
		}

		//ana 5 h 10 o elenxos gt epivarinei tis epanalipseis
		if((counter % 5) == 0)
			MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_MIN, new_comm);

		if(global_flag == 1) break;

		//anathesi
   		temp=matrix;
		matrix=new_matrix;
		new_matrix=temp;

   		counter++;
   	}


   	local_end_time = MPI_Wtime();
	local_elapsed = local_end_time - local_start_time;

	MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, new_comm);
	if(my_rank==0)
    	printf("Total MPI_Wtime is: %1.8f\n", elapsed);


	/*********************************************end of internal code*********************************************/

	

	//to matrix exei to new matrix
#	pragma omp parallel for num_threads(threads_num) shared(p_rows,p_cols,matrix,pix,total_p_cols) private(i,j) collapse(2)
	for(i=0; i< p_rows; i++){
		for (j = 0; j < p_cols; j++){
			read_buffer[i*p_cols+j]=matrix[(i+1)*(total_p_cols) + (j+pix)];
		}
	}

	//open new file for writing
	file_open_error = MPI_File_open(MPI_COMM_WORLD, output_filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh2);
 
	// https://www.hpc.ntnu.no/display/hpc/Reading+from+MPI+Files
  	if (file_open_error != MPI_SUCCESS) {
    	MPI_Abort(MPI_COMM_WORLD, file_open_error);
  	}

	MPI_File_set_view(fh2, my_offset, MPI_CHAR, subarray, "native", MPI_INFO_NULL);

	MPI_File_write_all(fh2, read_buffer, num_of_bytes, MPI_CHAR, &status); 

	MPI_Get_count(&status, MPI_CHAR, &count);
    	//printf("%3d: read %d bytes\n", my_rank, count);


    //free types and pointers and communicator
    MPI_Type_free(&subarray);
	MPI_Type_free(&neighbor_col);
	MPI_Type_free(&neighbor_row);
	MPI_Comm_free(&new_comm);

	MPI_File_close(&fh);
	MPI_File_close(&fh2);

	free(filename);
	free(type);
	free(output_filename);
	free(read_buffer);
	free(matrix);
	free(new_matrix);
	purge_filter(h);

	MPI_Finalize();

	return 0;
}
