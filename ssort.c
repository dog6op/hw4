/* Parallel sample sort
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>


static int compare(const void *a, const void *b)
{
  int *da = (int *)a;
  int *db = (int *)b;

  if (*da > *db)
    return 1;
  else if (*da < *db)
    return -1;
  else
    return 0;
}

int main( int argc, char *argv[])
{
  int rank, size;
  int root = 0;
  int i, N, j, n;
  int *vec, *sample, *all, *splitter, *send, *recv, *sdispls, *rdispls, *newvec;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  /* Number of random numbers per processor (this should be increased
   * for actual tests or could be passed in through the command line */
  N = 100;

  vec = calloc(N, sizeof(int));
  sample = calloc(size - 1, sizeof(int));
  splitter = calloc(size - 1, sizeof(int));
  send = calloc(size, sizeof(int));
  recv = calloc(size, sizeof(int));
  sdispls = calloc(size, sizeof(int));
  rdispls = calloc(size, sizeof(int));
    
  /* seed random number generator differently on every core */
  srand((unsigned int) (rank + 393919));

  /* fill vector with random integers */
  for (i = 0; i < N; ++i) {
    vec[i] = rand();
  }
  printf("rank: %d, first entry: %d\n", rank, vec[0]);

  /* sort locally */
  qsort(vec, N, sizeof(int), compare);

  /* randomly sample s entries from vector or select local splitters,
   * i.e., every N/P-th entry of the sorted vector */
  for(i = 0; i < size - 1; i++){
      sample[i] = vec[(i+1) * N/size];
  }

  /* every processor communicates the selected entries
   * to the root processor; use for instance an MPI_Gather */
  if(rank == root){
        all = calloc(size * (size - 1), sizeof(int));
  }
  MPI_Gather(sample, size - 1, MPI_INT, all, size - 1, MPI_INT, root, MPI_COMM_WORLD);
    
  /* root processor does a sort, determinates splitters that
   * split the data into P buckets of approximately the same size */
  if (rank == root) {
      qsort(all, size * (size - 1), sizeof(int), compare);
      for (i = 0; i < size - 1; i++){
          splitter[i] = all[(i+1)*(size-1)-1];
      }
   }

  /* root process broadcasts splitters */
  MPI_Bcast(splitter, size - 1, MPI_INT, root, MPI_COMM_WORLD);

  /* every processor uses the obtained splitters to decide
   * which integers need to be sent to which other processor (local bins) */
  i = 0;
  for (j = 0; j < N; j++) {
        if (vec[j] < splitter[i]){
            send[i]++;
        }
        else if (i < size - 1) {
            i++;
            send[i]++;
        }
        else {
            send[size - 1]++;
        }
   }
    
  /* send and receive: either you use MPI_AlltoallV, or
   * (and that might be easier), use an MPI_Alltoall to share
   * with every processor how many integers it should expect,
   * and then use MPI_Send and MPI_Recv to exchange the data */
  MPI_Alltoall(send, 1, MPI_INT, recv, 1, MPI_INT, MPI_COMM_WORLD);
  n = 0;
  rdispls[0] = 0;
  rdispls[0] = 0;
  for (i = 0; i < size; i++){
        n += recv[i];
  }
  for(i = 1; i < size; i++){
      sdispls[i] = sdispls[i-1] + send[i-1];
      rdispls[i] = rdispls[i-1] + recv[i-1];
  }
  newvec = calloc(n, sizeof(int));
  MPI_Alltoallv(vec, send, sdispls, MPI_INT, newvec, recv, rdispls, MPI_INT, MPI_COMM_WORLD);

  /* do a local sort */
  qsort(newvec, n, sizeof(int), compare);
    
  /* every processor writes its result to a file */
  FILE* fd = NULL;
  char filename[256];
  snprintf(filename, 256, "output%02d.txt", rank);
  fd = fopen(filename,"w+");
  if(NULL == fd) {
        printf("Error opening file \n");
        return 1;
  }
  fprintf(fd, "rank %d\n", rank);
  for(i = 0; i < n; i++) {
        fprintf(fd, "%d\n", newvec[i]);
  }
  fclose(fd);

  free(vec);
  free(sample);
  free(all);
  free(splitter);
  free(send);
  free(recv);
  free(sdispls);
  free(rdispls);
  free(newvec);
    
  MPI_Finalize();
  return 0;
}
