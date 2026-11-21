#include "ExVectrDSP/calibration/stillness_detector.hpp"

#include <math.h>

namespace VCTR {
namespace DSP {
namespace Calib {

namespace {
double stdOf(double sum, double sum2, uint32_t n) {
  if (n < 2)
    return 0;
  const double mean = sum / n;
  const double var = sum2 / n - mean * mean;
  return var > 0 ? sqrt(var) : 0;
}
} // namespace

void MotionStats::addGyro(const double w[3], const double bias[3]) {
  double n2 = 0;
  for (int k = 0; k < 3; k++) {
    gSum[k] += w[k];
    gSum2[k] += w[k] * w[k];
    if (w[k] < gMin[k])
      gMin[k] = w[k];
    if (w[k] > gMax[k])
      gMax[k] = w[k];
    const double c = w[k] - bias[k];
    n2 += c * c;
  }
  const double n = sqrt(n2);
  if (n > gNormMax)
    gNormMax = n;
  nGyro++;
}

void MotionStats::addAccel(const double a[3]) {
  double n2 = 0;
  for (int k = 0; k < 3; k++) {
    aSum[k] += a[k];
    aSum2[k] += a[k] * a[k];
    n2 += a[k] * a[k];
  }
  const double n = sqrt(n2);
  aNormSum += n;
  aNormSum2 += n * n;
  nAccel++;
}

void MotionStats::merge(const MotionStats &o) {
  nGyro += o.nGyro;
  nAccel += o.nAccel;
  for (int k = 0; k < 3; k++) {
    gSum[k] += o.gSum[k];
    gSum2[k] += o.gSum2[k];
    if (o.gMin[k] < gMin[k])
      gMin[k] = o.gMin[k];
    if (o.gMax[k] > gMax[k])
      gMax[k] = o.gMax[k];
    aSum[k] += o.aSum[k];
    aSum2[k] += o.aSum2[k];
  }
  if (o.gNormMax > gNormMax)
    gNormMax = o.gNormMax;
  aNormSum += o.aNormSum;
  aNormSum2 += o.aNormSum2;
}

double MotionStats::gyroStd(int k) const { return stdOf(gSum[k], gSum2[k], nGyro); }
double MotionStats::accelStd(int k) const { return stdOf(aSum[k], aSum2[k], nAccel); }
double MotionStats::accelNormStd() const { return stdOf(aNormSum, aNormSum2, nAccel); }

StillnessDetector::StillnessDetector(const StillnessConfig &cfg) : cfg_(cfg) {
  if (cfg_.windowBlocks == 0)
    cfg_.windowBlocks = 1;
  if (cfg_.windowBlocks > kMaxBlocks)
    cfg_.windowBlocks = kMaxBlocks;
  if (cfg_.blockNs <= 0)
    cfg_.blockNs = 100000000LL;
  reset();
}

void StillnessDetector::reset() {
  head_ = 0;
  count_ = 0;
  cur_.clear();
  started_ = false;
  still_ = false;
  winGyroNormMax_ = 0;
  winAccelStdMax_ = 0;
}

void StillnessDetector::setGyroBias(const double bias[3]) {
  for (int k = 0; k < 3; k++)
    bias_[k] = bias[k];
}

void StillnessDetector::advance(int64_t tNs) {
  if (!started_) {
    started_ = true;
    blockStart_ = tNs;
    return;
  }
  if (tNs < blockStart_) { // time went backwards: start over
    reset();
    started_ = true;
    blockStart_ = tNs;
    return;
  }
  if (tNs - blockStart_ < cfg_.blockNs)
    return;
  const int64_t elapsed = (tNs - blockStart_) / cfg_.blockNs;
  if (elapsed > 2) { // gap without samples: the window is not continuous
    head_ = 0;
    count_ = 0;
    cur_.clear();
    still_ = false;
    blockStart_ = tNs;
    return;
  }
  closeBlock();
  blockStart_ += cfg_.blockNs;
  if (tNs - blockStart_ >= cfg_.blockNs) { // one empty block in between
    closeBlock();
    blockStart_ += cfg_.blockNs;
  }
}

void StillnessDetector::closeBlock() {
  blocks_[head_] = cur_;
  head_ = uint8_t((head_ + 1) % cfg_.windowBlocks);
  if (count_ < cfg_.windowBlocks)
    count_++;
  cur_.clear();
  evaluate();
}

void StillnessDetector::evaluate() {
  if (count_ < cfg_.windowBlocks) {
    still_ = false;
    return;
  }
  MotionStats w;
  bool emptyBlock = false;
  for (uint8_t i = 0; i < cfg_.windowBlocks; i++) {
    if (blocks_[i].nAccel == 0 && blocks_[i].nGyro == 0)
      emptyBlock = true;
    w.merge(blocks_[i]);
  }
  winGyroNormMax_ = w.gNormMax;
  winAccelStdMax_ = 0;
  bool ok = !emptyBlock && w.nAccel >= 2;
  if (cfg_.maxGyroNorm > 0 && !(w.nGyro > 0 && w.gNormMax < cfg_.maxGyroNorm))
    ok = false;
  for (int k = 0; k < 3; k++) {
    const double as = w.accelStd(k);
    if (as > winAccelStdMax_)
      winAccelStdMax_ = as;
    if (cfg_.maxAccelStd > 0 && !(as < cfg_.maxAccelStd))
      ok = false;
    if (cfg_.maxGyroStd > 0 && !(w.gyroStd(k) < cfg_.maxGyroStd))
      ok = false;
    if (cfg_.maxGyroP2P > 0 && !(w.gyroP2P(k) < cfg_.maxGyroP2P))
      ok = false;
  }
  if (cfg_.maxAccelNormStd > 0 && !(w.accelNormStd() < cfg_.maxAccelNormStd))
    ok = false;
  still_ = ok;
}

void StillnessDetector::addGyro(int64_t tNs, const double w[3]) {
  if (!isfinite(w[0]) || !isfinite(w[1]) || !isfinite(w[2]))
    return;
  advance(tNs);
  cur_.addGyro(w, bias_);
}

void StillnessDetector::addAccel(int64_t tNs, const double a[3]) {
  if (!isfinite(a[0]) || !isfinite(a[1]) || !isfinite(a[2]))
    return;
  advance(tNs);
  cur_.addAccel(a);
}

} // namespace Calib
} // namespace DSP
} // namespace VCTR
