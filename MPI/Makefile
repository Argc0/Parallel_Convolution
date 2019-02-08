OBJS    = mpi_main.o mpi_functions.o filters.o
SOURCE  = mpi_main.c mpi_functions.c filters.c
HEADER  = mpi_functions.h
OUT     = mpi_exe
CC      = mpicc
FLAGS   = -g -Wall -std=c99
# -g option enables debugging mode
# -c flag generates object code for separate files

all: $(OBJS)
	$(CC) $(FLAGS) $(SOURCE) -o $(OUT) -lm

# clean house
clean:
	rm -f $(OBJS) $(OUT)

# do a bit of accounting
count:
	wc $(SOURCE) $(HEADER)
