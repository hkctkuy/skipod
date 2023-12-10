#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int main(int argc, char* argv[]) {
    int i, j;
    int N = atoi(argv[1]);

    double matrix[N*N];
    memset(matrix, 1, N*N);
    
    double vector[N];
    memset(vector, 1, N);
    
    double result[N];
    memset(result, 0, N);

    double start = omp_get_wtime();

    #pragma omp parallel
    {
        double result_private[N];
        memset(result_private, 0, N);
        int i, j;
        #pragma omp for
        for(i = 0; i < N; i++) {
            for(j = 0; j < N; j++) {
                result_private[i] += matrix[i * N + j] * vector[j];
            }
        }
        #pragma omp critical
        {
            for(i = 0; i < N; i++) result[i] += result_private[i];
        }
    }
    double end = omp_get_wtime();

    printf("%f\n", end - start);
    return 0;
}

