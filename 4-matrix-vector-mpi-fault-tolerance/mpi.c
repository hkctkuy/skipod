#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <setjmp.h>

#include <mpi.h>
#include <mpi-ext.h>
// Rank of process to be killed
int TO_KILL = 1;
MPI_Comm MPI_COMM_CUSTOM;
jmp_buf jbuf;

void mul(
	int N,
	double **matrix,
	double *vector,
	double *result
) {
    // Get the number of processes and the rank of the current process
    int size, rank;
    MPI_Comm_size(MPI_COMM_CUSTOM, &size);
    MPI_Comm_rank(MPI_COMM_CUSTOM, &rank);

    // The matrix and vector are divided into chunks and distributed among the processes
    int chunk_size = N / size;
    int start = rank * chunk_size;
    int end = start + chunk_size;

    int i, j;
	double start_time = MPI_Wtime();
    for(j = start; j < end; j++) {
        result[j] = 0;
        for(i = 0; i < N; i++) {
            result[j] += vector[i] * matrix[i][j];   
        }
    }
	if (rank == TO_KILL) {
		printf("[Process %d]: has been killed\n", rank);
        raise(SIGKILL);
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
			MPI_COMM_CUSTOM,
			&request
		);
	}
	// Wait sending
	MPI_Barrier(MPI_COMM_CUSTOM);
    if (rank == 0) {
		// Collect result by master
		int rank;
		for (rank = 1; rank < size; rank++) {
			MPI_Irecv(
				result + rank * chunk_size,
				chunk_size,
				MPI_DOUBLE,
				rank,
				0,
				MPI_COMM_CUSTOM,
				&request
			);
		}
		// Show result
        double end_time = MPI_Wtime();
		printf("Time: %f\n", end_time - start_time);
		printf("Result:\n");
		for (i = 0; i < N; i++) {
			printf("%f\n", result[i]);
		}
    }

    // Finalize MPI
    MPI_Finalize();
	return;
}

static void errhandler(MPI_Comm *comm, int *err, ...) {
    // Do not kill anyone else 
	TO_KILL = -1;

    int len;
    char errstr[MPI_MAX_ERROR_STRING];

    MPI_Error_string(*err, errstr, &len);
    errstr[len] = 0;

    printf("Captured error: %s\n", errstr);

    MPIX_Comm_shrink(*comm, &MPI_COMM_CUSTOM);
    MPI_Barrier(MPI_COMM_CUSTOM);

    longjmp(jbuf, 0);
}

int main(int argc, char* argv[]) {
    int N = atoi(argv[1]);
    int i, j;

    double **matrix = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) {
        matrix[i] = (double*)malloc(N * sizeof(double));
        for(j = 0; j < N ; j++) {
				matrix[i][j] = i + j;
			}
    }
    double *vector = (double*)malloc(N * sizeof(double));
    for (i = 0; i < N; i++) {
        vector[i] = i;
    }
    double *result = (double*)malloc(N * sizeof(double));

    // Initialize MPI
    MPI_Init(&argc, &argv);
	MPI_COMM_CUSTOM = MPI_COMM_WORLD;
	// Add error handler
	MPI_Errhandler errh;
    MPI_Comm_create_errhandler(errhandler, &errh);
    MPI_Comm_set_errhandler(MPI_COMM_CUSTOM, errh);
	
	MPI_Barrier(MPI_COMM_CUSTOM);

	// Run calculation
	setjmp(jbuf);
	mul(N, matrix, vector, result);

    return 0;
}
