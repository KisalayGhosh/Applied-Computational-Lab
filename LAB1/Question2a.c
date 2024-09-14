#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Declare the LAPACK function dgbsv for solving banded linear systems
extern void dgbsv_(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);

int main() {
    int n = 100;  // Matrix dimension (n x n)
    int kl = 1, ku = 1;  // Number of sub-diagonals (kl) and super-diagonals (ku)
    int ldab = 2 * kl + ku + 1;  // Leading dimension of the banded matrix storage (formula: 2*kl + ku + 1)
    int nrhs = 1, info;  // nrhs = number of right-hand side vectors, info = status code
    int *ipiv = (int *)malloc(n * sizeof(int));  // Array for pivot indices (size n)

    // Allocate memory for the banded matrix A in the special banded storage format
    double *ab = (double *)malloc(ldab * n * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));  // Right-hand side vector
    double *x_exact = (double *)malloc(n * sizeof(double));  // Exact solution vector for comparison
    double *error = (double *)malloc(n * sizeof(double));  // Error vector to compute L2 norm

    // Initialize the banded matrix A (stored in banded format)
    // Fill all entries in the banded matrix storage with zeros initially
    for (int i = 0; i < ldab * n; i++) ab[i] = 0.0;

    for (int i = 0; i < n; i++) {
        // Diagonal elements (8.0 on the main diagonal)
        ab[kl + ku + i * ldab] = 8.0;
        // Sub-diagonal (-1.0 on the lower band)
        if (i > 0) ab[kl + (ku - 1) + i * ldab] = -1.0;
        // Super-diagonal (2.0 on the upper band)
        if (i < n - 1) ab[kl + ku + 1 + i * ldab] = 2.0;
    }

    // Initialize the right-hand side vector b based on the problem description
    b[0] = 9.0;
    b[1] = 10.0;
    for (int i = 2; i < n - 2; i++) {
        b[i] = 13.0;  // Middle elements of b
    }
    b[n - 2] = 11.0;  // Second last element
    b[n - 1] = 12.0;  // Last element

    // Initialize the exact solution vector x_exact (all ones)
    for (int i = 0; i < n; i++) {
        x_exact[i] = 1.0;
    }

    // Call LAPACK's dgbsv to solve A * x = b using the banded matrix storage
    dgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &n, &info);

    // Check if the solution was successful
    if (info == 0) {
        // Print the computed solution vector x (now stored in b)
        printf("Solution x:\n");
        for (int i = 0; i < n; i++) {
            printf("x[%d] = %f\n", i, b[i]);
        }
    } else {
        // Error handling: if dgbsv fails, print the error code
        printf("dgbsv failed with info = %d\n", info);
        return -1;
    }

    // Compute the error vector (difference between computed and exact solution) and the L2 norm of the error
    double l2_norm = 0.0;
    for (int i = 0; i < n; i++) {
        error[i] = b[i] - x_exact[i];  // Compute the error for each element
        l2_norm += error[i] * error[i];  // Sum the squares of the errors
    }
    l2_norm = sqrt(l2_norm);  // Take the square root to compute the L2 norm

    // Print the L2 norm of the error
    printf("\nL2 norm of the error: %f\n", l2_norm);

    // Free the dynamically allocated memory
    free(ipiv);
    free(ab);
    free(b);
    free(x_exact);
    free(error);

    return 0;  // Program finished successfully
}
