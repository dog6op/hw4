/* MPI-parallel Jacobi smoothing to solve -u''=f
 * Global vector has N unknowns, each processor works with its
 * part, which has lN = N/p unknowns.
 * Author: Georg Stadler
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "util.h"
#include <string.h>

/* compuate global residual, assuming ghost values are updated */
double compute_residual(double *lu, int lN, double invhsq)
{
    int i;
    double tmp, gres = 0.0, lres = 0.0;
    
    for (i = 1; i <= lN; i++){
        for (j = 1; j <= lN; j++){
            tmp = ((4.0*lu(i*(lN+2)+j)- lu((i-1)*(lN+2)+j) - lu(i*(lN+2)+(j-1)) - lu((i+1)*(lN+2)+j) - lu(i*[lN+2]+(j+1))) * invhsq - 1);
            lres += tmp * tmp;
        }
    }
    /* use allreduce for convenience; a reduce would also be sufficient */
    MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(gres);
}


int main(int argc, char * argv[])
{
    int mpirank, i, j, p, N, lN, iter, max_iters;
    MPI_Status status, status1;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    /* get name of host running MPI process */
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);
    
    sscanf(argv[1], "%d", &N);
    sscanf(argv[2], "%d", &max_iters);
    
    /* compute number of unknowns handled by each process */
    lN = N / sqrt(p);
    if ((N % p != 0) && mpirank == 0 ) {
        printf("N: %d, local N: %d\n", N, lN);
        printf("Exiting. N must be a multiple of p\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    timestamp_type time1, time2;
    get_timestamp(&time1);
    
    /* Allocation of vectors, including left/upper and right/lower ghost points */
    double * lu    = (double *) calloc(sizeof(double), (lN + 2) * (lN + 2));
    double * lunew = (double *) calloc(sizeof(double), (lN + 2) * (lN + 2));
    double * lutemp;
    
    double h = 1.0 / (N + 1);
    double hsq = h * h;
    double invhsq = 1./hsq;
    double gres, gres0, tol = 1e-5;
    
    /* initial residual */
    gres0 = compute_residual(lu, lN, invhsq);
    gres = gres0;
    
    for (iter = 0; iter < max_iters && gres/gres0 > tol; iter++) {
        
        /* Jacobi step for local points */
        for (i = 1; i <= lN; i++){
            for(j = 1; j <= lN; j++){
                lunew[i*(lN+2)+j]  = 0.25 * (hsq + lu((i-1)*(lN+2)+j) + lu(i*(lN+2)+(j-1)) + lu((i+1)*(lN+2)+j) + lu(i*[lN+2]+(j+1)));
            }
        }
        
        
        /* communicate ghost values */
        if (mpirank % sqrt(p) < sqrt(p) - 1) {
            for(i = 1; i <= lN; i++){
                MPI_Send(&(lunew[(i+1)*(lN+2)-2]), 1, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD);
                MPI_Recv(&(lunew[i*(lN+2)]), 1, MPI_DOUBLE, mpirank+1, 123, MPI_COMM_WORLD, &status);
            }
        }
        if (mpirank % sqrt(p) > 0) {
            for(i = 1; i <= lN; i++){
                MPI_Send(&(lunew[i*(lN+2)+1]), 1, MPI_DOUBLE, mpirank-1, 123, MPI_COMM_WORLD);
                MPI_Recv(&(lunew[(i+1)*(lN+2)-1]), 1, MPI_DOUBLE, mpirank-1, 124, MPI_COMM_WORLD, &status1);
            }
        }
        if (mpirank < p - sqrt(p)) {
            for(i = 1; i <= lN; i++){
                MPI_Send(&(lunew[i]), 1, MPI_DOUBLE, mpirank+sqrt(p), 126, MPI_COMM_WORLD);
                MPI_Recv(&(lunew[lN*(lN+2)+i]), 1, MPI_DOUBLE, mpirank+sqrt(p), 125, MPI_COMM_WORLD, &status1);
            }
        }
        if (mpirank >= sqrt(p)) {
            for(i = 1; i <= lN; i++){
                MPI_Send(&(lunew[(lN-1)*(lN+2)+i]), 1, MPI_DOUBLE, mpirank-sqrt(p), 125, MPI_COMM_WORLD);
                MPI_Recv(&(lunew[i]), 1, MPI_DOUBLE, mpirank-sqrt(p), 126, MPI_COMM_WORLD, &status1);
            }
        }
        
        
        /* copy newu to u using pointer flipping */
        lutemp = lu; lu = lunew; lunew = lutemp;
        if (0 == (iter % 10)) {
            gres = compute_residual(lu, lN, invhsq);
            if (0 == mpirank) {
                printf("Iter %d: Residual: %g\n", iter, gres);
            }
        }
    }
    
    /* Clean up */
    free(lu);
    free(lunew);
    
    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    get_timestamp(&time2);
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    if (0 == mpirank) {
        printf("Time elapsed is %f seconds.\n", elapsed);
    }
    MPI_Finalize();
    return 0;
}