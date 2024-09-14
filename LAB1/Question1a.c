#include <stdio.h>
#include <stdlib.h>

// Declare the LAPACK function dgesv for solving linear systems
extern void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

int main() {
    int n = 100;               // Matrix dimension (n x n)
    int lda = n, ldb = n, nrhs = 1, info;   // lda and ldb are leading dimensions, nrhs = number of RHS vectors
    int *ipiv = (int *)malloc(n * sizeof(int));  // Array for pivot indices (size n)
    double *A = (double *)malloc(n * n * sizeof(double));  // Matrix A (size n x n)
    double *x = (double *)malloc(n * sizeof(double));      // Initial solution vector x (size n)
    double *b = (double *)malloc(n * sizeof(double));      // Right-hand side vector b (size n)

    // Initialize matrix A according to the problem specification
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // If i == j (diagonal) or j == n - 1 (last column), A[i][j] = 1
            if (i == j || j == n - 1) {
                A[i * n + j] = 1.0;
            }
            // If j > i, A[i][j] = -10
            else if (j > i) {
                A[i * n + j] = -10.0;
            }
            // Otherwise, A[i][j] = 0 (lower triangular part)
            else {
                A[i * n + j] = 0.0;
            }
        }
        // Initialize vector x such that x_i = i/n
        x[i] = (double)(i + 1) / n;
    }

    // Compute the right-hand side vector b = A * x
    for (int i = 0; i < n; i++) {
        b[i] = 0.0;  // Initialize b[i] to 0
        for (int j = 0; j < n; j++) {
            b[i] += A[i * n + j] * x[j];  // b[i] = Sum of A[i][j] * x[j]
        }
    }

    // Call LAPACK's dgesv to solve the system A * x' = b
    dgesv_(&n, &nrhs, A, &lda, ipiv, b, &ldb, &info);

    // Check if the solution was successful
    if (info == 0) {
        // Print the solution vector x' (now stored in b)
        printf("Solution x':\n");
        for (int i = 0; i < n; i++) {
            printf("x'[%d] = %e\n", i, b[i]);
        }
    } else {
        // Handle errors from dgesv
        printf("dgesv failed with info = %d\n", info);
        // Free memory and return with error code
        free(ipiv);
        free(A);
        free(x);
        free(b);
        return -1;
    }

    // Compare the original vector x with the computed solution x' (stored in b)
    printf("\nComparison between x and x':\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %e, x'[%d] = %e, difference = %e\n", i, x[i], i, b[i], x[i] - b[i]);
    }

    // Free dynamically allocated memory
    free(ipiv);
    free(A);
    free(x);
    free(b);

    return 0;  // Program finished successfully
}
