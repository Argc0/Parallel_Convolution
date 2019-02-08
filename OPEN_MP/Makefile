OBJS    = openmp_main.o open_mp_functions.o filters.o
SOURCE  = openmp_main.c open_mp_functions.c filters.c
HEADER  = open_mp_functions.h
OUT     = openmp_exe
CC      = mpicc
FLAGS   = -g -Wall -fopenmp -std=c99
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
