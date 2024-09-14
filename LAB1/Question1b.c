#include <stdio.h>
#include <stdlib.h>

// Declare the LAPACK function dgesv for solving linear systems
extern void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

int main() {
    int n = 2;  // Matrix dimension for a 2x2 system
    int lda = n, ldb = n, nrhs = 1, info;  // lda = leading dimension of A, nrhs = number of RHS vectors, ldb = leading dimension of b
    int *ipiv = (int *)malloc(n * sizeof(int));  // Array for pivot indices (size n)
    
    // Matrix A is stored in column-major order as required by LAPACK
    // A = [4.5  3.1]
    //     [1.6  1.1]
    double A[4] = {4.5, 1.6, 3.1, 1.1};
    
    // Right-hand side vector b for the system A * x = b
    double b[2] = {19.249, 6.843};  // Original b vector
    
    // Right-hand side vector c for the system A * y = c
    double c[2] = {19.25, 6.84};  // Another right-hand side vector c

    // Solve the system A * x = b using LAPACK's dgesv function
    dgesv_(&n, &nrhs, A, &lda, ipiv, b, &ldb, &info);
    if (info == 0) {
        // Solution was successful, print the solution x
        printf("Solution for A * x = b:\n");
        printf("x1 = %f, x2 = %f\n", b[0], b[1]);  // b now contains the solution vector x
    } else {
        // Error handling: if dgesv fails, print the error code
        printf("dgesv failed with info = %d\n", info);
    }

    // Solve the system A * y = c using the same matrix A (reusing the factorization)
    dgesv_(&n, &nrhs, A, &lda, ipiv, c, &ldb, &info);
    if (info == 0) {
        // Solution was successful, print the solution y
        printf("Solution for A * y = c:\n");
        printf("y1 = %f, y2 = %f\n", c[0], c[1]);  // c now contains the solution vector y
    } else {
        // Error handling: if dgesv fails, print the error code
        printf("dgesv failed with info = %d\n", info);
    }

    // Compare the two solutions (x and y)
    printf("\nComparison between solutions:\n");
    printf("Difference in x1 and y1: %f\n", b[0] - c[0]);  // Difference between x1 and y1
    printf("Difference in x2 and y2: %f\n", b[1] - c[1]);  // Difference between x2 and y2

    // Free dynamically allocated memory
    free(ipiv);

    return 0;  // Program finished successfully
}
