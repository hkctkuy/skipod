#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char* argv[]) {
    int matrix_size = atoi(argv[1]);
    int i, j;

    double **matrix = (double**)malloc(matrix_size * sizeof(double*));
    for (i = 0; i < matrix_size; i++) {
        matrix[i] = (double*)malloc(matrix_size * sizeof(double));
        for(j = 0; j < matrix_size ; j++) {
	    matrix[i][j] = i + j;
	}
    }
    double *vector = (double*)malloc(matrix_size * sizeof(double));
    for (i = 0; i < matrix_size; i++) {
        vector[i] = i;
    }
    double *result = (double*)malloc(matrix_size * sizeof(double));

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the number of processes and the rank of the current process
    int num_procs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // The matrix and vector are divided into chunks and distributed among the processes
    int chunk_size = matrix_size / num_procs;
    int start = rank * chunk_size;
    int end = start + chunk_size;

    double start_time = MPI_Wtime();
    for(j = start; j < end; j++) {
        result[j] = 0;
        for(i = 0; i < matrix_size; i++) {
            result[j] += vector[i] * matrix[i][j];   
        }
    }
	MPI_Request request;
	if (rank) {
		// Send result to master
		MPI_Isend(
			result + start,
			chunk_size,
			MPI_DOUBLE,
			0,
			0,
			MPI_COMM_WORLD,
			&request
		);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
		// Collect result by master
		int rank;
		for (rank = 1; rank < num_procs; rank++) {
			MPI_Irecv(
				result + rank * chunk_size,
				chunk_size,
				MPI_DOUBLE,
				rank,
				0,
				MPI_COMM_WORLD,
				&request
			);
		}

		// Show result
        double end_time = MPI_Wtime();
		printf("Time: %f\n", end_time - start_time);
		printf("Result:\n");
		for (i = 0; i < matrix_size; i++) {
			printf("%f\n", result[i]);
		}
    }

    // Finalize MPI
    MPI_Finalize();	

    return 0;
}
