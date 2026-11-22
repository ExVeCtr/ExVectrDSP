#include "ExVectrDSP/calibration/mag_ellipsoid_fit.hpp"
#include "ExVectrDSP/calibration/cal_linalg.hpp"

#include <math.h>
#include <stdio.h>

namespace VCTR {
namespace DSP {
namespace Calib {

namespace {

constexpr size_t kNp = 9; // b0 b1 b2 w00 w11 w22 w01 w02 w12

void unpack(const double *p, double *b, double *W) {
  b[0] = p[0];
  b[1] = p[1];
  b[2] = p[2];
  W[0] = p[3];
  W[4] = p[4];
  W[8] = p[5];
  W[1] = W[3] = p[6];
  W[2] = W[6] = p[7];
  W[5] = W[7] = p[8];
}

void pack(const double *b, const double *W, double *p) {
  p[0] = b[0];
  p[1] = b[1];
  p[2] = b[2];
  p[3] = W[0];
  p[4] = W[4];
  p[5] = W[8];
  p[6] = W[1];
  p[7] = W[2];
  p[8] = W[5];
}

/// Residual r = |W (x - b)| - 1 and its gradient (normalised units).
double residualAndGradient(const double *p, const double *xn, double *grad) {
  double b[3], W[9];
  unpack(p, b, W);
  const double y[3] = {xn[0] - b[0], xn[1] - b[1], xn[2] - b[2]};
  const double u[3] = {W[0] * y[0] + W[1] * y[1] + W[2] * y[2], W[3] * y[0] + W[4] * y[1] + W[5] * y[2],
                       W[6] * y[0] + W[7] * y[1] + W[8] * y[2]};
  const double nrm = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  if (grad != nullptr) {
    const double inv = nrm > 1e-12 ? 1 / nrm : 0;
    grad[0] = -(W[0] * u[0] + W[3] * u[1] + W[6] * u[2]) * inv;
    grad[1] = -(W[1] * u[0] + W[4] * u[1] + W[7] * u[2]) * inv;
    grad[2] = -(W[2] * u[0] + W[5] * u[1] + W[8] * u[2]) * inv;
    grad[3] = u[0] * y[0] * inv;
    grad[4] = u[1] * y[1] * inv;
    grad[5] = u[2] * y[2] * inv;
    grad[6] = (u[0] * y[1] + u[1] * y[0]) * inv;
    grad[7] = (u[0] * y[2] + u[2] * y[0]) * inv;
    grad[8] = (u[1] * y[2] + u[2] * y[1]) * inv;
  }
  return nrm - 1.0;
}

/// Symmetric positive-definite 3x3 square root (Denman-Beavers iteration). False if not converging.
bool sqrtSym3(const double *M, double *S) {
  double Y[9], Z[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  for (int i = 0; i < 9; i++)
    Y[i] = M[i];
  for (int it = 0; it < 60; it++) {
    double Yi[9], Zi[9];
    if (!invert3(Y, Yi) || !invert3(Z, Zi))
      return false;
    double delta = 0;
    for (int i = 0; i < 9; i++) {
      const double yn = 0.5 * (Y[i] + Zi[i]);
      const double zn = 0.5 * (Z[i] + Yi[i]);
      delta += fabs(yn - Y[i]);
      Y[i] = yn;
      Z[i] = zn;
    }
    if (delta < 1e-14)
      break;
  }
  for (int i = 0; i < 9; i++) {
    if (!isfinite(Y[i]))
      return false;
    S[i] = Y[i];
  }
  // symmetrise
  S[1] = S[3] = 0.5 * (Y[1] + Y[3]);
  S[2] = S[6] = 0.5 * (Y[2] + Y[6]);
  S[5] = S[7] = 0.5 * (Y[5] + Y[7]);
  return true;
}

} // namespace

const char *magFailReason(const MagEllipsoidResult &r) {
  const uint16_t c = r.checks;
  if (!(c & kMagCheckSamples))
    return "too few stored samples, keep tumbling";
  if (!(c & kMagCheckCoverage))
    return "directions not covered, tumble through all orientations";
  if (!(c & kMagCheckSolved))
    return "fit not solvable (degenerate data)";
  if (!(c & kMagCheckRejected))
    return "too many disturbed samples (metal or current nearby), move away";
  if (!(c & kMagCheckResidualRms) || !(c & kMagCheckResidualMax))
    return "field length varies too much after the fit (disturbance or distortion beyond the model)";
  if (!(c & kMagCheckWellDetermined))
    return "parameters not well determined, tumble through more varied orientations";
  if (!(c & kMagCheckSoftIron))
    return "soft-iron distortion too large (axis ratio or eigenvalues out of range)";
  if (!(c & kMagCheckOffset))
    return "hard-iron offset too large";
  if (!(c & kMagCheckResidualNominal)) {
    static char buf[64]; // single-threaded operator log
    snprintf(buf, sizeof(buf), "field noise or disturbance during the tumble (rms %.1f %%)", r.residualRms * 100.0);
    return buf;
  }
  if (!(c & kMagCheckRejectedNominal))
    return "some disturbed samples were rejected";
  return "all checks passed";
}

double magCalibrate(const MagEllipsoidResult &r, const double x[3], double out[3]) {
  const double y[3] = {x[0] - r.offset[0], x[1] - r.offset[1], x[2] - r.offset[2]};
  out[0] = r.A[0] * y[0] + r.A[1] * y[1] + r.A[2] * y[2];
  out[1] = r.A[3] * y[0] + r.A[4] * y[1] + r.A[5] * y[2];
  out[2] = r.A[6] * y[0] + r.A[7] * y[1] + r.A[8] * y[2];
  return sqrt(out[0] * out[0] + out[1] * out[1] + out[2] * out[2]);
}

MagEllipsoidFit::MagEllipsoidFit(float (*buffer)[3], size_t capacity, const MagEllipsoidConfig &cfg)
    : buf_(buffer), cap_(buffer != nullptr ? (capacity < kMaxStored ? capacity : kMaxStored) : 0), cfg_(cfg) {
  reset();
}

void MagEllipsoidFit::reset() {
  for (int i = 0; i < 81; i++)
    eN_[i] = 0;
  for (int i = 0; i < 9; i++)
    eY_[i] = 0;
  for (int i = 0; i < 16; i++)
    sN_[i] = 0;
  for (int i = 0; i < 4; i++)
    sY_[i] = 0;
  shift_[0] = shift_[1] = shift_[2] = 0;
  haveShift_ = false;
  nAcc_ = 0;
  nStored_ = 0;
  centre_[0] = centre_[1] = centre_[2] = 0;
  haveLast_ = false;
  coverage_.reset();
}

void MagEllipsoidFit::toNormalised(const float *x, double xn[3]) const {
  for (int k = 0; k < 3; k++)
    xn[k] = (double(x[k]) - shift_[k]) / cfg_.normScale;
}

void MagEllipsoidFit::addSample(const double m[3]) {
  if (!isfinite(m[0]) || !isfinite(m[1]) || !isfinite(m[2]))
    return;
  if (!haveShift_) {
    for (int k = 0; k < 3; k++)
      shift_[k] = m[k];
    haveShift_ = true;
  }
  const double x = (m[0] - shift_[0]) / cfg_.normScale, y = (m[1] - shift_[1]) / cfg_.normScale,
               z = (m[2] - shift_[2]) / cfg_.normScale;
  const double phi[9] = {y * y - x * x, z * z - x * x, x * y, x * z, y * z, x, y, z, 1.0};
  const double target = -x * x;
  for (int r = 0; r < 9; r++) {
    for (int c = 0; c < 9; c++)
      eN_[r * 9 + c] += phi[r] * phi[c];
    eY_[r] += phi[r] * target;
  }
  const double m4[4] = {x, y, z, 1.0};
  const double sq = x * x + y * y + z * z;
  for (int r = 0; r < 4; r++) {
    for (int c = 0; c < 4; c++)
      sN_[r * 4 + c] += m4[r] * m4[c];
    sY_[r] += m4[r] * sq;
  }
  nAcc_++;

  // Direction-thinned storage about the current centre.
  double d[3] = {m[0] - centre_[0], m[1] - centre_[1], m[2] - centre_[2]};
  const double n = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  if (!(n > 1e-12) || nStored_ >= cap_)
    return;
  d[0] /= n;
  d[1] /= n;
  d[2] /= n;
  if (haveLast_) {
    const double dot = d[0] * lastStoredDir_[0] + d[1] * lastStoredDir_[1] + d[2] * lastStoredDir_[2];
    if (dot > cos(cfg_.storeAngleRad))
      return;
  }
  buf_[nStored_][0] = float(m[0]);
  buf_[nStored_][1] = float(m[1]);
  buf_[nStored_][2] = float(m[2]);
  nStored_++;
  lastStoredDir_[0] = d[0];
  lastStoredDir_[1] = d[1];
  lastStoredDir_[2] = d[2];
  haveLast_ = true;
  coverage_.add(d);
}

bool MagEllipsoidFit::solveSphereNormalised(Params &p) const {
  if (nAcc_ < 4)
    return false;
  double A[16], b[4];
  for (int i = 0; i < 16; i++)
    A[i] = sN_[i];
  for (int i = 0; i < 4; i++)
    b[i] = sY_[i];
  if (!solveLinear(A, b, 4))
    return false;
  const double c[3] = {0.5 * b[0], 0.5 * b[1], 0.5 * b[2]};
  const double r2 = b[3] + c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
  if (!(r2 > 0) || !isfinite(r2))
    return false;
  const double invR = 1.0 / sqrt(r2);
  const double W[9] = {invR, 0, 0, 0, invR, 0, 0, 0, invR};
  pack(c, W, p.v);
  return true;
}

bool MagEllipsoidFit::solveAlgebraic(Params &p) const {
  if (nAcc_ < 10)
    return false;
  double N[81], th[9];
  for (int i = 0; i < 81; i++)
    N[i] = eN_[i];
  for (int i = 0; i < 9; i++)
    th[i] = eY_[i];
  if (!solveLinear(N, th, 9))
    return false;
  const double a = 1.0 - th[0] - th[1];
  const double Q[9] = {a, th[2] / 2, th[3] / 2, th[2] / 2, th[0], th[4] / 2, th[3] / 2, th[4] / 2, th[1]};
  const double lin[3] = {th[5], th[6], th[7]};
  const double c0 = th[8];
  double eig[3];
  if (!symmetricEigenvalues(Q, 3, eig) || !(eig[0] > 1e-9))
    return false; // not an ellipsoid
  double Qi[9];
  if (!invert3(Q, Qi))
    return false;
  double b[3];
  for (int r = 0; r < 3; r++)
    b[r] = -0.5 * (Qi[3 * r] * lin[0] + Qi[3 * r + 1] * lin[1] + Qi[3 * r + 2] * lin[2]);
  double Qb[3];
  for (int r = 0; r < 3; r++)
    Qb[r] = Q[3 * r] * b[0] + Q[3 * r + 1] * b[1] + Q[3 * r + 2] * b[2];
  const double k = b[0] * Qb[0] + b[1] * Qb[1] + b[2] * Qb[2] - c0;
  if (!(k > 1e-12) || !isfinite(k))
    return false;
  double M[9], W[9];
  for (int i = 0; i < 9; i++)
    M[i] = Q[i] / k;
  if (!sqrtSym3(M, W))
    return false;
  pack(b, W, p.v);
  for (int i = 0; i < 9; i++)
    if (!isfinite(p.v[i]))
      return false;
  return true;
}

void MagEllipsoidFit::recomputeCoverage() {
  coverage_.reset();
  for (size_t i = 0; i < nStored_; i++) {
    const double d[3] = {buf_[i][0] - centre_[0], buf_[i][1] - centre_[1], buf_[i][2] - centre_[2]};
    coverage_.add(d);
  }
}

bool MagEllipsoidFit::updateEstimate() {
  Params p;
  bool ok = cfg_.softIron && solveAlgebraic(p);
  if (!ok)
    ok = solveSphereNormalised(p);
  if (!ok)
    return false;
  for (int k = 0; k < 3; k++)
    centre_[k] = shift_[k] + cfg_.normScale * p.v[k];
  recomputeCoverage();
  return true;
}

void MagEllipsoidFit::refine(Params &p, const uint8_t *mask, int &iterations) const {
  auto cost = [&](const double *pp) {
    double s = 0;
    for (size_t i = 0; i < nStored_; i++) {
      if (!mask[i])
        continue;
      double xn[3];
      toNormalised(buf_[i], xn);
      const double r = residualAndGradient(pp, xn, nullptr);
      s += r * r;
    }
    return s;
  };
  double lambda = 1e-3;
  double c = cost(p.v);
  if (!isfinite(c))
    return;
  bool converged = false;
  int it = 0;
  for (; !converged && it < cfg_.maxIterations; it++) {
    double JtJ[kNp * kNp], Jtr[kNp];
    for (size_t i = 0; i < kNp * kNp; i++)
      JtJ[i] = 0;
    for (size_t i = 0; i < kNp; i++)
      Jtr[i] = 0;
    for (size_t i = 0; i < nStored_; i++) {
      if (!mask[i])
        continue;
      double xn[3], g[kNp];
      toNormalised(buf_[i], xn);
      const double r = residualAndGradient(p.v, xn, g);
      for (size_t a = 0; a < kNp; a++) {
        Jtr[a] += g[a] * r;
        for (size_t b = a; b < kNp; b++)
          JtJ[a * kNp + b] += g[a] * g[b];
      }
    }
    for (size_t a = 0; a < kNp; a++)
      for (size_t b = 0; b < a; b++)
        JtJ[a * kNp + b] = JtJ[b * kNp + a];
    bool improved = false;
    for (int tries = 0; tries < 12 && !improved; tries++) {
      double M[kNp * kNp], d[kNp];
      for (size_t a = 0; a < kNp * kNp; a++)
        M[a] = JtJ[a];
      for (size_t a = 0; a < kNp; a++) {
        M[a * kNp + a] += lambda * (JtJ[a * kNp + a] > 1e-12 ? JtJ[a * kNp + a] : 1e-12);
        d[a] = -Jtr[a];
      }
      if (!solveLinear(M, d, kNp)) {
        lambda *= 10;
        continue;
      }
      double pn[kNp], step = 0;
      for (size_t a = 0; a < kNp; a++) {
        pn[a] = p.v[a] + d[a];
        step += d[a] * d[a];
      }
      const double cn = cost(pn);
      if (isfinite(cn) && cn <= c) {
        const double rel = (c - cn) / (c + 1e-30);
        for (size_t a = 0; a < kNp; a++)
          p.v[a] = pn[a];
        c = cn;
        lambda = lambda / 10 > 1e-12 ? lambda / 10 : 1e-12;
        improved = true;
        if (step < 1e-24 || rel < 1e-12)
          converged = true;
      } else {
        lambda *= 10;
      }
    }
    if (!improved)
      break;
  }
  iterations += it;
}

MagEllipsoidResult MagEllipsoidFit::evaluate() {
  MagEllipsoidResult res;
  res.softIron = cfg_.softIron;
  res.nAccumulated = nAcc_;
  res.nStored = nStored_;
  if (nStored_ >= cfg_.minStored)
    res.checks |= kMagCheckSamples;

  Params p;
  bool init = cfg_.softIron && solveAlgebraic(p);
  if (!init)
    init = solveSphereNormalised(p);
  if (!init) {
    recomputeCoverage();
    res.coverage = coverage_.fractionAtLeast(cfg_.minBinCount);
    if (res.coverage >= cfg_.minCoverage)
      res.checks |= kMagCheckCoverage;
    res.result = (nAcc_ < 10) ? CalResult::FailCoverage : CalResult::Numeric;
    return res;
  }

  for (size_t i = 0; i < nStored_; i++)
    mask_[i] = 1;

  int iterations = 0;
  if (cfg_.softIron && nStored_ >= 12) {
    refine(p, mask_, iterations);
    for (int pass = 0; pass < cfg_.outlierPasses; pass++) {
      // robust sigma from the median absolute residual over ALL stored samples
      for (size_t i = 0; i < nStored_; i++) {
        double xn[3];
        toNormalised(buf_[i], xn);
        absRes_[i] = float(fabs(residualAndGradient(p.v, xn, nullptr)));
      }
      float tmp[kMaxStored];
      for (size_t i = 0; i < nStored_; i++)
        tmp[i] = absRes_[i];
      for (size_t gap = nStored_ / 2; gap > 0; gap /= 2) // shell sort
        for (size_t i = gap; i < nStored_; i++) {
          const float v = tmp[i];
          size_t j = i;
          for (; j >= gap && tmp[j - gap] > v; j -= gap)
            tmp[j] = tmp[j - gap];
          tmp[j] = v;
        }
      const double med = nStored_ % 2 ? tmp[nStored_ / 2] : 0.5 * (tmp[nStored_ / 2 - 1] + tmp[nStored_ / 2]);
      double sigma = 1.4826 * med;
      if (sigma < cfg_.noiseFloor)
        sigma = cfg_.noiseFloor;
      size_t nIn = 0;
      for (size_t i = 0; i < nStored_; i++) {
        mask_[i] = absRes_[i] <= cfg_.outlierSigmas * sigma ? 1 : 0;
        nIn += mask_[i];
      }
      if (nIn < 30) { // not enough inliers: keep everything
        for (size_t i = 0; i < nStored_; i++)
          mask_[i] = 1;
        break;
      }
      refine(p, mask_, iterations);
    }
  }
  res.iterations = iterations;

  // Physical parameters.
  double bn[3], W[9];
  unpack(p.v, bn, W);
  bool finiteOk = true;
  for (int i = 0; i < 9; i++)
    if (!isfinite(p.v[i]))
      finiteOk = false;
  const double aScale = cfg_.rRef / cfg_.normScale;
  for (int k = 0; k < 3; k++)
    res.offset[k] = shift_[k] + cfg_.normScale * bn[k];
  for (int k = 0; k < 9; k++)
    res.A[k] = aScale * W[k];
  for (int k = 0; k < 3; k++)
    centre_[k] = res.offset[k];
  recomputeCoverage();
  res.coverage = coverage_.fractionAtLeast(cfg_.minBinCount);
  if (res.coverage >= cfg_.minCoverage)
    res.checks |= kMagCheckCoverage;

  if (!finiteOk || !symmetricEigenvalues(res.A, 3, res.eig)) {
    res.result = CalResult::Numeric;
    return res;
  }
  res.checks |= kMagCheckSolved;
  res.scale = (res.eig[0] + res.eig[1] + res.eig[2]) / 3.0;
  res.radius = res.scale > 0 ? cfg_.rRef / res.scale : 0;
  res.axisRatio = res.eig[0] > 0 ? res.eig[2] / res.eig[0] : INFINITY;

  // Residuals over the inliers.
  double sum2 = 0, rmax = 0, sumRaw = 0, sumRaw2 = 0;
  size_t nIn = 0;
  for (size_t i = 0; i < nStored_; i++) {
    const double x[3] = {buf_[i][0], buf_[i][1], buf_[i][2]};
    const double raw = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    sumRaw += raw;
    sumRaw2 += raw * raw;
    if (!mask_[i])
      continue;
    double cal[3];
    const double r = magCalibrate(res, x, cal) / cfg_.rRef - 1.0;
    sum2 += r * r;
    if (!isfinite(r) || fabs(r) > rmax)
      rmax = isfinite(r) ? fabs(r) : INFINITY;
    nIn++;
  }
  res.nInliers = nIn;
  res.nRejected = nStored_ - nIn;
  if (nIn > 0) {
    res.residualRms = sqrt(sum2 / double(nIn));
    res.residualMax = rmax;
    const double mean = sumRaw / double(nStored_);
    const double var = sumRaw2 / double(nStored_) - mean * mean;
    res.rawSpread = mean > 0 ? sqrt(var > 0 ? var : 0) / mean : 0;
    if (res.residualRms < cfg_.maxResidualRms)
      res.checks |= kMagCheckResidualRms;
    if (res.residualMax < cfg_.maxResidualMax)
      res.checks |= kMagCheckResidualMax;
    if (res.residualRms < cfg_.warnResidualRms)
      res.checks |= kMagCheckResidualNominal;
  }
  const double rejFrac = nStored_ > 0 ? double(res.nRejected) / double(nStored_) : 0;
  if (rejFrac <= cfg_.maxRejectedFraction)
    res.checks |= kMagCheckRejected;
  if (rejFrac <= cfg_.warnRejectedFraction)
    res.checks |= kMagCheckRejectedNominal;

  // Standard errors from sigma^2 (J'J)^-1 on the inliers (soft-iron mode only: the sphere mode has no matrix fit).
  if (cfg_.softIron) {
    double JtJ[kNp * kNp];
    for (size_t i = 0; i < kNp * kNp; i++)
      JtJ[i] = 0;
    for (size_t i = 0; i < nStored_; i++) {
      if (!mask_[i])
        continue;
      double xn[3], g[kNp];
      toNormalised(buf_[i], xn);
      residualAndGradient(p.v, xn, g);
      for (size_t a = 0; a < kNp; a++)
        for (size_t b = 0; b < kNp; b++)
          JtJ[a * kNp + b] += g[a] * g[b];
    }
    const double sig = res.residualRms > cfg_.noiseFloor ? res.residualRms : cfg_.noiseFloor;
    const double meanW = (W[0] + W[4] + W[8]) / 3.0;
    bool ok = nIn >= kNp && meanW > 0;
    double sd[kNp];
    for (size_t k = 0; k < kNp && ok; k++) {
      double M[kNp * kNp], col[kNp];
      for (size_t a = 0; a < kNp * kNp; a++)
        M[a] = JtJ[a];
      for (size_t a = 0; a < kNp; a++)
        col[a] = a == k ? 1.0 : 0.0;
      if (!solveLinear(M, col, kNp) || !(col[k] > 0) || !isfinite(col[k]))
        ok = false;
      else
        sd[k] = sig * sqrt(col[k]) / meanW; // relative to the radius (offset) or the mean eigenvalue (matrix)
    }
    res.sigmaValid = ok;
    if (ok) {
      double so = 0, sm = 0;
      for (size_t k = 0; k < 3; k++)
        if (sd[k] > so)
          so = sd[k];
      for (size_t k = 3; k < kNp; k++)
        if (sd[k] > sm)
          sm = sd[k];
      res.sigmaOffset = so;
      res.sigmaMatrix = sm;
      if (so < cfg_.maxOffsetSigma && sm < cfg_.maxMatrixSigma)
        res.checks |= kMagCheckWellDetermined;
    }
  } else {
    res.checks |= kMagCheckWellDetermined;
  }

  if (res.eig[0] > 0 && res.eig[0] >= cfg_.eigMinRel * res.scale && res.eig[2] <= cfg_.eigMaxRel * res.scale &&
      res.axisRatio < cfg_.maxAxisRatio)
    res.checks |= kMagCheckSoftIron;
  const double off = sqrt(res.offset[0] * res.offset[0] + res.offset[1] * res.offset[1] +
                          res.offset[2] * res.offset[2]);
  if (res.radius > 0 && off < cfg_.maxOffsetRatio * res.radius)
    res.checks |= kMagCheckOffset;

  if (!(res.checks & kMagCheckSamples) || !(res.checks & kMagCheckCoverage))
    res.result = CalResult::FailCoverage;
  else if (cfg_.softIron && !res.sigmaValid)
    res.result = CalResult::Numeric;
  else if (!(res.checks & kMagCheckResidualRms) || !(res.checks & kMagCheckResidualMax) ||
           !(res.checks & kMagCheckRejected))
    res.result = CalResult::FailResidual;
  else if (!(res.checks & kMagCheckWellDetermined))
    res.result = CalResult::FailCoverage;
  else if (!(res.checks & kMagCheckSoftIron) || !(res.checks & kMagCheckOffset))
    res.result = CalResult::FailRange;
  else if (!(res.checks & kMagCheckRejectedNominal) || !(res.checks & kMagCheckResidualNominal))
    res.result = CalResult::Warn;
  else
    res.result = CalResult::Pass;
  return res;
}

} // namespace Calib
} // namespace DSP
} // namespace VCTR
