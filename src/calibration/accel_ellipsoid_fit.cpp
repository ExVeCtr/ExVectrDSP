#include "ExVectrDSP/calibration/accel_ellipsoid_fit.hpp"
#include "ExVectrDSP/calibration/cal_linalg.hpp"

#include <math.h>

namespace VCTR {
namespace DSP {
namespace Calib {

namespace {

// Parameter vector (normalised units, x / g): p = [b0, b1, b2, a00, a11, a22, a01, a02, a12].
constexpr size_t kMaxParams = 9;

void unpack(const double *p, double *b, double *A) {
  b[0] = p[0];
  b[1] = p[1];
  b[2] = p[2];
  A[0] = p[3];
  A[4] = p[4];
  A[8] = p[5];
  A[1] = A[3] = p[6];
  A[2] = A[6] = p[7];
  A[5] = A[7] = p[8];
}

/// Residual r = |A y| - 1 and its gradient over the first np parameters (normalised units).
double residualAndGradient(const double *p, size_t np, const double *xg, double *grad) {
  double b[3], A[9];
  unpack(p, b, A);
  const double y[3] = {xg[0] - b[0], xg[1] - b[1], xg[2] - b[2]};
  const double u[3] = {A[0] * y[0] + A[1] * y[1] + A[2] * y[2], A[3] * y[0] + A[4] * y[1] + A[5] * y[2],
                       A[6] * y[0] + A[7] * y[1] + A[8] * y[2]};
  const double nrm = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  if (grad != nullptr) {
    const double inv = nrm > 1e-12 ? 1 / nrm : 0;
    // d|u|/db = -(A^T u)/|u| = -(A u)/|u| (A symmetric)
    grad[0] = -(A[0] * u[0] + A[3] * u[1] + A[6] * u[2]) * inv;
    grad[1] = -(A[1] * u[0] + A[4] * u[1] + A[7] * u[2]) * inv;
    grad[2] = -(A[2] * u[0] + A[5] * u[1] + A[8] * u[2]) * inv;
    grad[3] = u[0] * y[0] * inv;
    grad[4] = u[1] * y[1] * inv;
    grad[5] = u[2] * y[2] * inv;
    if (np > 6) {
      grad[6] = (u[0] * y[1] + u[1] * y[0]) * inv;
      grad[7] = (u[0] * y[2] + u[2] * y[0]) * inv;
      grad[8] = (u[1] * y[2] + u[2] * y[1]) * inv;
    }
  }
  return nrm - 1.0;
}

double cost(const double *p, size_t np, const double (*poses)[3], size_t n, double g) {
  double s = 0;
  for (size_t i = 0; i < n; i++) {
    const double xg[3] = {poses[i][0] / g, poses[i][1] / g, poses[i][2] / g};
    const double r = residualAndGradient(p, np, xg, nullptr);
    s += r * r;
  }
  return s;
}

void normalEquations(const double *p, size_t np, const double (*poses)[3], size_t n, double g, double *JtJ,
                     double *Jtr) {
  for (size_t i = 0; i < np * np; i++)
    JtJ[i] = 0;
  for (size_t i = 0; i < np; i++)
    Jtr[i] = 0;
  for (size_t i = 0; i < n; i++) {
    const double xg[3] = {poses[i][0] / g, poses[i][1] / g, poses[i][2] / g};
    double grad[kMaxParams];
    const double r = residualAndGradient(p, np, xg, grad);
    for (size_t a = 0; a < np; a++) {
      Jtr[a] += grad[a] * r;
      for (size_t c = a; c < np; c++)
        JtJ[a * np + c] += grad[a] * grad[c];
    }
  }
  for (size_t a = 0; a < np; a++)
    for (size_t c = 0; c < a; c++)
      JtJ[a * np + c] = JtJ[c * np + a];
}

} // namespace

const char *accelParamName(int idx) {
  static const char *const names[9] = {"offset X", "offset Y", "offset Z", "scale XX", "scale YY",
                                       "scale ZZ",  "skew XY",  "skew XZ",  "skew YZ"};
  return (idx >= 0 && idx < 9) ? names[idx] : "?";
}

double accelParamLimit(int idx, const AccelFitConfig &cfg) { return idx < 3 ? cfg.maxOffsetSigma : cfg.maxScaleSigma; }

const char *accelParamAdvice(int idx) {
  static const char *const adv[9] = {"gravity along X (+X and -X)",   "gravity along Y (+Y and -Y)",
                                     "gravity along Z (+Z and -Z)",   "gravity along X (+X and -X)",
                                     "gravity along Y (+Y and -Y)",   "gravity along Z (+Z and -Z)",
                                     "gravity at 45 deg between X and Y", "gravity at 45 deg between X and Z",
                                     "gravity at 45 deg between Y and Z"};
  return (idx >= 0 && idx < 9) ? adv[idx] : "";
}

size_t accelWeakParams(const AccelFitResult &r, const AccelFitConfig &cfg, AccelWeakParam *out, size_t max) {
  if (!r.sigmaValid || max == 0)
    return 0;
  const size_t np = r.diagonal ? 6 : 9;
  AccelWeakParam all[9];
  for (size_t i = 0; i < np; i++) {
    all[i].idx = int(i);
    all[i].sigma = r.sigma[i];
    all[i].limit = accelParamLimit(int(i), cfg);
  }
  // insertion sort by sigma / limit, descending
  for (size_t i = 1; i < np; i++)
    for (size_t k = i; k > 0 && all[k].sigma / all[k].limit > all[k - 1].sigma / all[k - 1].limit; k--) {
      const AccelWeakParam t = all[k];
      all[k] = all[k - 1];
      all[k - 1] = t;
    }
  size_t cnt = 0;
  for (size_t i = 0; i < np && cnt < max; i++)
    if (all[i].sigma >= all[i].limit)
      out[cnt++] = all[i];
  if (cnt == 0)
    out[cnt++] = all[0];
  return cnt;
}

double accelCalibratedNorm(const AccelFitResult &r, const double x[3]) {
  const double y[3] = {x[0] - r.offset[0], x[1] - r.offset[1], x[2] - r.offset[2]};
  const double u[3] = {r.A[0] * y[0] + r.A[1] * y[1] + r.A[2] * y[2], r.A[3] * y[0] + r.A[4] * y[1] + r.A[5] * y[2],
                       r.A[6] * y[0] + r.A[7] * y[1] + r.A[8] * y[2]};
  return sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
}

AccelFitResult fitAccel(const double (*poses)[3], size_t n, const AccelFitConfig &cfg) {
  AccelFitResult res;
  res.nUsed = n;
  const double g = cfg.gRef;

  if (n >= cfg.minPosesFull)
    res.checks |= kAccelCheckPosesFull;

  for (size_t i = 0; i < n; i++) {
    if (!isfinite(poses[i][0]) || !isfinite(poses[i][1]) || !isfinite(poses[i][2])) {
      res.result = CalResult::Numeric;
      return res;
    }
  }
  if (n < cfg.minPosesDiag || !(g > 0)) {
    res.result = CalResult::FailCoverage;
    return res;
  }

  res.diagonal = n < cfg.minPosesFull;
  const size_t np = res.diagonal ? 6 : 9;

  double p[kMaxParams] = {0, 0, 0, 1, 1, 1, 0, 0, 0};
  double JtJ[kMaxParams * kMaxParams], Jtr[kMaxParams];
  double c = cost(p, np, poses, n, g);
  double lambda = 1e-3;
  bool solved = isfinite(c);
  int it = 0;
  bool converged = false;
  for (; solved && !converged && it < cfg.maxIterations; it++) {
    normalEquations(p, np, poses, n, g, JtJ, Jtr);
    bool improved = false;
    for (int tries = 0; tries < 12 && !improved; tries++) {
      double M[kMaxParams * kMaxParams], d[kMaxParams];
      for (size_t a = 0; a < np * np; a++)
        M[a] = JtJ[a];
      for (size_t a = 0; a < np; a++) {
        M[a * np + a] += lambda * (JtJ[a * np + a] > 1e-12 ? JtJ[a * np + a] : 1e-12);
        d[a] = -Jtr[a];
      }
      if (!solveLinear(M, d, np)) {
        lambda *= 10;
        continue;
      }
      double pn[kMaxParams];
      for (size_t a = 0; a < kMaxParams; a++)
        pn[a] = p[a];
      double step = 0;
      for (size_t a = 0; a < np; a++) {
        pn[a] += d[a];
        step += d[a] * d[a];
      }
      const double cn = cost(pn, np, poses, n, g);
      if (isfinite(cn) && cn <= c) {
        const double rel = (c - cn) / (c + 1e-30);
        for (size_t a = 0; a < kMaxParams; a++)
          p[a] = pn[a];
        c = cn;
        lambda = lambda / 10 > 1e-12 ? lambda / 10 : 1e-12;
        improved = true;
        if (step < 1e-24 || rel < 1e-15)
          converged = true;
      } else {
        lambda *= 10;
      }
    }
    if (!improved)
      break; // no descent step any more: at the minimum (or stuck)
  }
  res.iterations = it;

  // Condition number of the normal matrix at the solution.
  normalEquations(p, np, poses, n, g, JtJ, Jtr);
  double eig[kMaxParams];
  if (symmetricEigenvalues(JtJ, np, eig) && eig[0] > 0 && isfinite(eig[np - 1]))
    res.condition = eig[np - 1] / eig[0];
  else
    res.condition = INFINITY;

  for (size_t a = 0; a < kMaxParams; a++)
    if (!isfinite(p[a]))
      solved = false;

  double b[3], A[9];
  unpack(p, b, A);
  for (int k = 0; k < 3; k++)
    res.offset[k] = b[k] * g;
  for (int k = 0; k < 9; k++)
    res.A[k] = A[k];

  if (solved)
    res.checks |= kAccelCheckSolved;
  if (res.condition < cfg.maxCondition)
    res.checks |= kAccelCheckCondition;

  // Residuals in m/s^2.
  double sum2 = 0, rmax = 0;
  for (size_t i = 0; i < n; i++) {
    const double r = accelCalibratedNorm(res, poses[i]) - g;
    sum2 += r * r;
    if (fabs(r) > rmax || !isfinite(r))
      rmax = isfinite(r) ? fabs(r) : INFINITY;
  }
  res.residualRms = sqrt(sum2 / double(n));
  res.residualMax = rmax;
  if (res.residualRms < cfg.maxResidualRms)
    res.checks |= kAccelCheckResidualRms;
  if (res.residualMax < cfg.maxResidualMax)
    res.checks |= kAccelCheckResidualMax;

  // Parameter standard errors: cov = sigma^2 (J^T J)^-1 in normalised units, offsets back to m/s^2.
  {
    const double sig = (res.residualRms > cfg.noiseFloor ? res.residualRms : cfg.noiseFloor) / g;
    bool ok = solved;
    double col[kMaxParams], M[kMaxParams * kMaxParams], var[kMaxParams];
    for (size_t k = 0; k < np && ok; k++) {
      for (size_t a = 0; a < np * np; a++)
        M[a] = JtJ[a];
      for (size_t a = 0; a < np; a++)
        col[a] = (a == k) ? 1.0 : 0.0;
      if (!solveLinear(M, col, np) || !(col[k] > 0) || !isfinite(col[k]))
        ok = false;
      else
        var[k] = col[k];
    }
    res.sigmaValid = ok;
    if (ok) {
      bool well = true;
      for (size_t k = 0; k < np; k++) {
        const double sd = sig * sqrt(var[k]);
        res.sigma[k] = k < 3 ? sd * g : sd;
        if (!(res.sigma[k] < accelParamLimit(int(k), cfg)))
          well = false;
      }
      res.wellDetermined = well;
    }
    if (res.wellDetermined)
      res.checks |= kAccelCheckCoverage;
  }

  bool eigOk = symmetricEigenvalues(res.A, 3, res.scaleEig);
  if (eigOk && res.scaleEig[0] >= cfg.scaleMin && res.scaleEig[2] <= cfg.scaleMax)
    res.checks |= kAccelCheckScale;
  if (fabs(res.A[1]) < cfg.maxOffDiag && fabs(res.A[2]) < cfg.maxOffDiag && fabs(res.A[5]) < cfg.maxOffDiag)
    res.checks |= kAccelCheckOffDiag;
  bool offOk = true, offNominal = true;
  for (int k = 0; k < 3; k++) {
    if (!(fabs(res.offset[k]) < cfg.maxOffset))
      offOk = false;
    if (!(fabs(res.offset[k]) < cfg.warnOffset))
      offNominal = false;
  }
  if (offOk)
    res.checks |= kAccelCheckOffset;
  if (offNominal)
    res.checks |= kAccelCheckOffsetNominal;

  if (!solved || !res.sigmaValid)
    res.result = CalResult::Numeric; // not solvable, or singular normal matrix (a parameter is not constrained at all)
  else if (!(res.checks & kAccelCheckResidualRms) || !(res.checks & kAccelCheckResidualMax))
    res.result = CalResult::FailResidual; // before the distribution: a bad pose also inflates the standard errors
  else if (!res.wellDetermined)
    res.result = CalResult::FailCoverage;
  else if (!(res.checks & kAccelCheckCondition))
    res.result = CalResult::Numeric;
  else if (!(res.checks & kAccelCheckScale) || !(res.checks & kAccelCheckOffDiag) || !(res.checks & kAccelCheckOffset))
    res.result = CalResult::FailRange;
  else if (res.diagonal || !offNominal)
    res.result = CalResult::Warn;
  else
    res.result = CalResult::Pass;
  return res;
}

} // namespace Calib
} // namespace DSP
} // namespace VCTR
