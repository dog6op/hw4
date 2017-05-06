CC=mpicc
CFLAGS=-O3 

all: mpi_solved1 mpi_solved2 mpi_solved3 mpi_solved4 mpi_solved5 mpi_solved6 jacobi-mpi2D ssort 

mpi_solved1: mpi_solved1.c
	$(CC) $(CFLAGS) mpi_solved1.c -o mpi_solved1
 
mpi_solved2: mpi_solved2.c
	$(CC) $(CFLAGS) mpi_solved2.c -o mpi_solved2

mpi_solved3: mpi_solved3.c
	$(CC) $(CFLAGS) mpi_solved3.c -o mpi_solved3

mpi_solved4: mpi_solved4.c
	$(CC) $(CFLAGS) mpi_solved4.c -o mpi_solved4

mpi_solved5: mpi_solved5.c
	$(CC) $(CFLAGS) mpi_solved5.c -o mpi_solved5

mpi_solved6: mpi_solved6.c
	$(CC) $(CFLAGS) mpi_solved6.c -o mpi_solved6

mpi_solved7: mpi_solved7.c
	$(CC) $(CFLAGS) mpi_solved7.c -o mpi_solved7

jacobi-mpi2D: jacobi-mpi2D.c
	$(CC) $(CFLAGS) jacobi-mpi2D.c -o jacobi-mpi2D

ssort: ssort.c
	$(CC) $(CFLAGS) ssort.c -o ssort

clean:
	rm mpi_solved1 mpi_solved2 mpi_solved3 mpi_solved4 mpi_solved5 mpi_solved6 mpi_solved7 jacobi-mpi2D ssort