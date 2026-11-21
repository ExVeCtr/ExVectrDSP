#include "ExVectrDSP/calibration/cal_linalg.hpp"
#include "ExVectrDSP/calibration/calibration_result.hpp"

#include <math.h>

namespace VCTR {
namespace DSP {
namespace Calib {

const char *toString(CalResult r) {
  switch (r) {
  case CalResult::None: return "NONE";
  case CalResult::Pass: return "PASS";
  case CalResult::Warn: return "WARN";
  case CalResult::FailResidual: return "FAIL_RESIDUAL";
  case CalResult::FailCoverage: return "FAIL_COVERAGE";
  case CalResult::FailRange: return "FAIL_RANGE";
  case CalResult::Motion: return "MOTION";
  case CalResult::Timeout: return "TIMEOUT";
  case CalResult::Numeric: return "NUMERIC";
  }
  return "?";
}

bool solveLinear(double *A, double *b, size_t n) {
  if (n == 0 || n > kLinalgMaxN)
    return false;
  double maxAbs = 0;
  for (size_t i = 0; i < n * n; i++) {
    if (!isfinite(A[i]))
      return false;
    if (fabs(A[i]) > maxAbs)
      maxAbs = fabs(A[i]);
  }
  for (size_t i = 0; i < n; i++)
    if (!isfinite(b[i]))
      return false;
  if (maxAbs <= 0)
    return false;
  const double tol = 1e-14 * maxAbs;

  for (size_t col = 0; col < n; col++) {
    size_t piv = col;
    for (size_t r = col + 1; r < n; r++)
      if (fabs(A[r * n + col]) > fabs(A[piv * n + col]))
        piv = r;
    if (fabs(A[piv * n + col]) < tol)
      return false;
    if (piv != col) {
      for (size_t c = 0; c < n; c++) {
        const double t = A[col * n + c];
        A[col * n + c] = A[piv * n + c];
        A[piv * n + c] = t;
      }
      const double t = b[col];
      b[col] = b[piv];
      b[piv] = t;
    }
    for (size_t r = col + 1; r < n; r++) {
      const double f = A[r * n + col] / A[col * n + col];
      if (f == 0)
        continue;
      for (size_t c = col; c < n; c++)
        A[r * n + c] -= f * A[col * n + c];
      b[r] -= f * b[col];
    }
  }
  for (size_t i = n; i-- > 0;) {
    double s = b[i];
    for (size_t c = i + 1; c < n; c++)
      s -= A[i * n + c] * b[c];
    b[i] = s / A[i * n + i];
    if (!isfinite(b[i]))
      return false;
  }
  return true;
}

bool symmetricEigenvalues(const double *Ain, size_t n, double *eig) {
  if (n == 0 || n > kLinalgMaxN)
    return false;
  double a[kLinalgMaxN * kLinalgMaxN];
  for (size_t r = 0; r < n; r++)
    for (size_t c = 0; c < n; c++) {
      const double v = (c >= r) ? Ain[r * n + c] : Ain[c * n + r];
      if (!isfinite(v))
        return false;
      a[r * n + c] = v;
    }

  for (int sweep = 0; sweep < 60; sweep++) {
    double off = 0, diag = 0;
    for (size_t r = 0; r < n; r++) {
      diag += a[r * n + r] * a[r * n + r];
      for (size_t c = r + 1; c < n; c++)
        off += a[r * n + c] * a[r * n + c];
    }
    if (off <= 1e-30 * (diag + 1e-300))
      break;
    for (size_t p = 0; p < n; p++) {
      for (size_t q = p + 1; q < n; q++) {
        const double apq = a[p * n + q];
        if (fabs(apq) < 1e-300)
          continue;
        const double app = a[p * n + p], aqq = a[q * n + q];
        const double theta = (aqq - app) / (2 * apq);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1));
        const double c = 1 / sqrt(t * t + 1), s = t * c;
        for (size_t k = 0; k < n; k++) {
          const double akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (size_t k = 0; k < n; k++) {
          const double apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
      }
    }
  }
  for (size_t i = 0; i < n; i++)
    eig[i] = a[i * n + i];
  for (size_t i = 1; i < n; i++) { // insertion sort, ascending
    const double v = eig[i];
    size_t j = i;
    while (j > 0 && eig[j - 1] > v) {
      eig[j] = eig[j - 1];
      j--;
    }
    eig[j] = v;
  }
  return true;
}

bool invert3(const double *M, double *Mi) {
  const double det = M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) +
                     M[2] * (M[3] * M[7] - M[4] * M[6]);
  if (!isfinite(det) || fabs(det) < 1e-12)
    return false;
  const double id = 1 / det;
  Mi[0] = (M[4] * M[8] - M[5] * M[7]) * id;
  Mi[1] = (M[2] * M[7] - M[1] * M[8]) * id;
  Mi[2] = (M[1] * M[5] - M[2] * M[4]) * id;
  Mi[3] = (M[5] * M[6] - M[3] * M[8]) * id;
  Mi[4] = (M[0] * M[8] - M[2] * M[6]) * id;
  Mi[5] = (M[2] * M[3] - M[0] * M[5]) * id;
  Mi[6] = (M[3] * M[7] - M[4] * M[6]) * id;
  Mi[7] = (M[1] * M[6] - M[0] * M[7]) * id;
  Mi[8] = (M[0] * M[4] - M[1] * M[3]) * id;
  return true;
}

void mul3(const double *A, const double *B, double *C) {
  for (int r = 0; r < 3; r++)
    for (int c = 0; c < 3; c++)
      C[r * 3 + c] = A[r * 3 + 0] * B[0 * 3 + c] + A[r * 3 + 1] * B[1 * 3 + c] + A[r * 3 + 2] * B[2 * 3 + c];
}

} // namespace Calib
} // namespace DSP
} // namespace VCTR
