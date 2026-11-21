#ifndef EXVECTRDSP_CAL_LINALG_HPP
#define EXVECTRDSP_CAL_LINALG_HPP

/**
 * Small dense linear algebra in double precision for the calibration fits (n <= 10).
 * Row-major arrays, no allocation. Pure code: <math.h> only.
 */

#include <stddef.h>

namespace VCTR {
namespace DSP {
namespace Calib {

constexpr size_t kLinalgMaxN = 10;

/**
 * Solves A x = b for an n x n system (Gaussian elimination with partial pivoting). A and b are overwritten; the
 * solution is returned in b. Returns false if n is out of range, an entry is not finite, or a pivot is below
 * 1e-14 times the largest absolute entry of A (singular for practical purposes).
 */
bool solveLinear(double *A, double *b, size_t n);

/**
 * Eigenvalues of a symmetric n x n matrix (cyclic Jacobi), ascending in eig. Only the upper triangle is read.
 * Returns false if n is out of range or an entry is not finite.
 */
bool symmetricEigenvalues(const double *A, size_t n, double *eig);

/// Inverse of a 3x3 matrix (row-major). Returns false if |det| < 1e-12 or not finite.
bool invert3(const double *M, double *Minv);

/// C = A * B for 3x3 row-major matrices (C may not alias A or B).
void mul3(const double *A, const double *B, double *C);

} // namespace Calib
} // namespace DSP
} // namespace VCTR

#endif
