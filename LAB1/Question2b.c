#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern void dgbsv_(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);

int main() {
    int n = 100;
    int kl = 1, ku = 1;  // Number of sub-diagonals and super-diagonals
    int ldab = 2 * kl + ku + 1;  // Leading dimension of banded matrix storage
    int nrhs = 1, info;
    int k = 5;  // Number of right-hand sides to solve
    int *ipiv = (int *)malloc(n * sizeof(int));

    // Allocate memory for banded storage of A
    double *ab = (double *)malloc(ldab * n * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));
    double *x_prev = (double *)malloc(n * sizeof(double));  // Stores the solution of the previous system

    // Initialize banded matrix A in storage
    for (int i = 0; i < ldab * n; i++) ab[i] = 0.0;  // Initialize with zeros

    for (int i = 0; i < n; i++) {
        // Diagonal elements
        ab[kl + ku + i * ldab] = 8.0;
        // Sub-diagonal
        if (i > 0) ab[kl + (ku - 1) + i * ldab] = -1.0;
        // Super-diagonal
        if (i < n - 1) ab[kl + ku + 1 + i * ldab] = 2.0;
    }

    // Initialize right-hand side vector b
    b[0] = 9.0;
    b[1] = 10.0;
    for (int i = 2; i < n - 2; i++) {
        b[i] = 13.0;
    }
    b[n - 2] = 11.0;
    b[n - 1] = 12.0;

    // Save original b for future use
    double *b_orig = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        b_orig[i] = b[i];
    }

    // Perform LU factorization only once
    dgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &n, &info);
    if (info != 0) {
        printf("LU factorization failed with info = %d\n", info);
        return -1;
    }

    // Solve for k right-hand sides
    for (int i = 1; i <= k; i++) {
        printf("\nSolving system %d...\n", i);

        // For b_1, use the original b
        if (i == 1) {
            for (int j = 0; j < n; j++) {
                b[j] = b_orig[j];
            }
        } else {
            // For b_k, use b_k = b_orig + x_(k-1)
            for (int j = 0; j < n; j++) {
                b[j] = b_orig[j] + x_prev[j];
            }
        }

        // Solve the system using the already factorized A
        dgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &n, &info);
        if (info == 0) {
            printf("Solution x_%d:\n", i);
            for (int j = 0; j < n; j++) {
                printf("x[%d] = %f\n", j, b[j]);
                x_prev[j] = b[j];  // Store the solution to use for b_(k+1)
            }
        } else {
            printf("Solving system %d failed with info = %d\n", i, info);
        }
    }

    // Free allocated memory
    free(ipiv);
    free(ab);
    free(b);
    free(b_orig);
    free(x_prev);

    return 0;
}
