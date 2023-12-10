#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

MPI_Status status;
MPI_Request reqs[2];
MPI_Status st[2];

int main(int argc, char *argv[]) {

	int rank, size;
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int tmp, num = rank;
    int row_sum = num;
    int sum = num;

    printf("[Process %d]: has local num %d\n", rank, num);
    
	MPI_Barrier(MPI_COMM_WORLD);

    // Step 1
    switch (rank % 4) {
        case 0:
            MPI_Send(&row_sum, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            break;
        case 1:
            MPI_Recv(&tmp, 1, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sum += tmp;
            row_sum += tmp;
            break;
        case 2:
            MPI_Recv(&tmp, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            row_sum += tmp;
            break;
        case 3:
            MPI_Send(&row_sum, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
            break;
    }

	// Step 2
    switch (rank %4) {
        case 1:
            MPI_Isend(&row_sum, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &reqs[0]);
            MPI_Irecv(&tmp, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &reqs[1]);
            MPI_Waitall(2, reqs, st);
            row_sum += tmp;
            break;
        case 2:
            MPI_Isend(&row_sum, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &reqs[0]);
            MPI_Irecv(&tmp, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &reqs[1]);
            MPI_Waitall(2, reqs, st);
            sum += tmp;
            row_sum += tmp;
            break;
    }

    // Step 3
    switch (rank % 4) {
        case 1:
            MPI_Send(&row_sum, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
            break;
        case 0:
            MPI_Recv(&tmp, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            row_sum = tmp;
            break;
        case 3:
            MPI_Recv(&tmp, 1, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            row_sum = tmp;
            sum = row_sum;
            break;
        case 2:
            MPI_Send(&row_sum, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            break;
    }

    // Steps 4-6
    if (rank / 4 != 0) {
        MPI_Recv(&tmp, 1, MPI_INT, rank - 4, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        sum += tmp;
        row_sum += tmp;
    }
    if (rank / 4 != 3) {
        MPI_Send(&row_sum, 1, MPI_INT, rank + 4, 0, MPI_COMM_WORLD);
    }

    printf("[Process %d]: has received sum: %d \n", rank, sum);
	MPI_Barrier(MPI_COMM_WORLD);

	return 0;
}
