#include "ExVectrDSP/calibration/mag_sphere_fit.hpp"
#include "ExVectrDSP/calibration/cal_linalg.hpp"

#include <math.h>

namespace VCTR {
namespace DSP {
namespace Calib {

MagSphereFit::MagSphereFit(float (*buffer)[3], size_t capacity, const MagFitConfig &cfg)
    : buf_(buffer), cap_(buffer != nullptr ? capacity : 0), cfg_(cfg) {
  reset();
}

void MagSphereFit::reset() {
  for (int i = 0; i < 16; i++)
    xTx_[i] = 0;
  for (int i = 0; i < 4; i++)
    xTy_[i] = 0;
  nAcc_ = 0;
  nStored_ = 0;
  centre_[0] = centre_[1] = centre_[2] = 0;
  radius_ = 0;
  haveEstimate_ = false;
  haveLast_ = false;
  coverage_.reset();
}

void MagSphereFit::addSample(const double m[3]) {
  if (!isfinite(m[0]) || !isfinite(m[1]) || !isfinite(m[2]))
    return;
  const double m4[4] = {m[0], m[1], m[2], 1.0};
  const double sq = m[0] * m[0] + m[1] * m[1] + m[2] * m[2];
  for (int r = 0; r < 4; r++) {
    for (int c = 0; c < 4; c++)
      xTx_[r * 4 + c] += m4[r] * m4[c];
    xTy_[r] += m4[r] * sq;
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

bool MagSphereFit::solveSphere(const double *xTx, const double *xTy, double centre[3], double &radius) {
  double A[16], b[4];
  for (int i = 0; i < 16; i++)
    A[i] = xTx[i];
  for (int i = 0; i < 4; i++)
    b[i] = xTy[i];
  if (!solveLinear(A, b, 4))
    return false;
  centre[0] = 0.5 * b[0];
  centre[1] = 0.5 * b[1];
  centre[2] = 0.5 * b[2];
  const double r2 = b[3] + centre[0] * centre[0] + centre[1] * centre[1] + centre[2] * centre[2];
  if (!(r2 > 0) || !isfinite(r2))
    return false;
  radius = sqrt(r2);
  return true;
}

void MagSphereFit::recomputeCoverage() {
  coverage_.reset();
  for (size_t i = 0; i < nStored_; i++) {
    const double d[3] = {buf_[i][0] - centre_[0], buf_[i][1] - centre_[1], buf_[i][2] - centre_[2]};
    coverage_.add(d);
  }
}

bool MagSphereFit::updateEstimate() {
  if (nAcc_ < 4)
    return false;
  double c[3], r;
  if (!solveSphere(xTx_, xTy_, c, r))
    return false;
  centre_[0] = c[0];
  centre_[1] = c[1];
  centre_[2] = c[2];
  radius_ = r;
  haveEstimate_ = true;
  recomputeCoverage();
  return true;
}

MagFitResult MagSphereFit::evaluate() {
  MagFitResult res;
  res.nAccumulated = nAcc_;
  const bool solved = updateEstimate();
  res.nStored = nStored_;
  res.coverage = coverage_.fractionAtLeast(cfg_.minBinCount);
  if (nStored_ >= cfg_.minStored)
    res.checks |= kMagCheckSamples;
  if (res.coverage >= cfg_.minCoverage)
    res.checks |= kMagCheckCoverage;

  if (!solved) {
    res.result = (nAcc_ < 4) ? CalResult::FailCoverage : CalResult::Numeric;
    return res;
  }
  res.checks |= kMagCheckSolved;
  for (int k = 0; k < 3; k++)
    res.offset[k] = centre_[k];
  res.radius = radius_;
  res.scale = cfg_.rRef / radius_;

  double sum2 = 0, rmax = 0, sumRaw = 0, sumRaw2 = 0;
  for (size_t i = 0; i < nStored_; i++) {
    const double d[3] = {buf_[i][0] - centre_[0], buf_[i][1] - centre_[1], buf_[i][2] - centre_[2]};
    const double nrm = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    const double r = nrm * res.scale / cfg_.rRef - 1.0;
    sum2 += r * r;
    if (fabs(r) > rmax)
      rmax = fabs(r);
    const double raw = sqrt(double(buf_[i][0]) * buf_[i][0] + double(buf_[i][1]) * buf_[i][1] +
                            double(buf_[i][2]) * buf_[i][2]);
    sumRaw += raw;
    sumRaw2 += raw * raw;
  }
  if (nStored_ > 0) {
    res.residualRms = sqrt(sum2 / double(nStored_));
    res.residualMax = rmax;
    const double mean = sumRaw / double(nStored_);
    const double var = sumRaw2 / double(nStored_) - mean * mean;
    res.rawSpread = mean > 0 ? sqrt(var > 0 ? var : 0) / mean : 0;
  }
  if (nStored_ > 0 && res.residualRms < cfg_.maxResidualRms)
    res.checks |= kMagCheckResidualRms;
  if (nStored_ > 0 && res.residualMax < cfg_.maxResidualMax)
    res.checks |= kMagCheckResidualMax;
  const double off = sqrt(centre_[0] * centre_[0] + centre_[1] * centre_[1] + centre_[2] * centre_[2]);
  if (off < cfg_.maxOffsetRatio * radius_)
    res.checks |= kMagCheckOffset;

  if (!(res.checks & kMagCheckSamples) || !(res.checks & kMagCheckCoverage))
    res.result = CalResult::FailCoverage;
  else if (!(res.checks & kMagCheckResidualRms) || !(res.checks & kMagCheckResidualMax))
    res.result = CalResult::FailResidual;
  else if (!(res.checks & kMagCheckOffset))
    res.result = CalResult::FailRange;
  else
    res.result = CalResult::Pass;
  return res;
}

} // namespace Calib
} // namespace DSP
} // namespace VCTR
