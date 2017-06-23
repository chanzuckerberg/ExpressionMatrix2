#ifndef CZI_EXPRESSION_MATRIX2_LAPACK_HPP
#define CZI_EXPRESSION_MATRIX2_LAPACK_HPP



// Declaration of some Lapack routines we need.
// This adding library lapack to the link line.
// Requires ubuntu package liblapack-dev.
// If not installed, use "apt install liblapack-dev" to install it.



// LAPACK routines for Householder QR factorization.
// These routines can be used to orthogonalize a set of n vectors, each of length m (n<=m).
// The vectors to be orthogonalized should be stored as the columns of A.
// After these routines are called, the columns of A will contain a set of n
// orthonormal vectors that span the same vector subspace of the initial set of vectors.
// This process is numerically much more stable than the simple Gram-Schmidt
// orthogonalizatiion procedure.
extern "C" void dgeqrf_(
    const int& m,
    const int& n,
    double* a,
    const int& lda,
    double* tau,
    double* work,
    const int& lwork,
    int& info
    );
extern "C" void dorgqr_(
    const int& m,
    const int& n,
    const int& k,
    double* a,
    const int& lda,
    double* tau,
    double* work,
    const int& lwork,
    int& info
    );


#endif
