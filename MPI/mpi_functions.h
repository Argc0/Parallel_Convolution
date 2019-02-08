#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <mpi.h>

#define S 1
#define D 2
#define MAX_NUM_C 100
#define RGB "rgb"
#define GREY "grey"

/* avoid linker errors
   https://stackoverflow.com/questions/9196801/how-to-define-an-array-of-strings-of-characters-in-header-file
*/
extern const int Blur[S+2][S+2];
extern const int Sharpen[S+2][S+2];
extern const int Edge[S+2][S+2];
extern const int EdgeH[S+2][S+2];
extern const int EdgeV[S+2][S+2];
extern const int GradientH[S+2][S+2];
extern const int GradientV[S+2][S+2];
extern const int SobelH[S+2][S+2];
extern const int SobelV[S+2][S+2];
extern const int Emboss[S+2][S+2];

/* Get arguments */
void get_args( int argc, char** argv, char** filename, int *rows, int *columns, char **type, int *selection, int *iterations);
/* Return chosen filter normalized */
float ** chosen_filter_normalizer(const int filter[S+2][S+2]);
/* Choose filter to be used */
float ** filter_normalizer(int choice);
/* Free chosen filter */
void purge_filter(float ** filter);
/* Initialization of matrix */
void init_matrix(unsigned char* matrix, unsigned char* new_matrix, unsigned char* read_buffer, int p_rows, int p_cols, int pix);
/* Initialization of neighbours */
void init_neighbours(MPI_Comm new_comm, int* dim_size, int my_rank, int ndims, int* North, int* South, int* East, int* West, int* n_w, int* n_e, int* s_w, int* s_e);
/* Set output filename */
void set_output_filename(char** output_filename, char* type, int rows, int columns, int selection);
#endif
