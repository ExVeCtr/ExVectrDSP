#include "ExVectrDSP/calibration/gyro_bias_collector.hpp"

#include <math.h>
#include <stdio.h>

namespace VCTR {
namespace DSP {
namespace Calib {

namespace {
const double kZero[3] = {0, 0, 0};
}

void GyroBiasCollector::start(int64_t nowNs, int64_t deadlineNs, int64_t windowNs) {
  cfg_.windowNs = windowNs > 0 ? windowNs : baseWindowNs_;
  state_ = State::Collecting;
  deadline_ = deadlineNs;
  restarts_ = 0;
  bestRun_ = rejected_ = 0;
  worstViol_ = GyroViolation();
  res_ = GyroBiasResult();
  diag_ = GyroBlockDiagnostics();
  last_ = GyroLastBlock();
  blockStarted_ = false;
  restartWindow();
  (void)nowNs;
}

void GyroBiasCollector::restartWindow() {
  block_.clear();
  window_.clear();
  for (int k = 0; k < 3; k++)
    blockMagSum_[k] = winMagSum_[k] = 0;
  blockMagN_ = winMagN_ = 0;
  blockTMin_ = winTMin_ = 1e300;
  blockTMax_ = winTMax_ = -1e300;
  blockTSum_ = winTSum_ = 0;
  blockTN_ = winTN_ = 0;
  windowCollectedNs_ = 0;
}

double GyroBiasCollector::progress() const {
  if (state_ == State::Done)
    return 1.0;
  if (state_ != State::Collecting || cfg_.windowNs <= 0)
    return 0;
  const double p = double(windowCollectedNs_) / double(cfg_.windowNs);
  return p > 1 ? 1 : p;
}

void GyroBiasCollector::advance(int64_t tNs) {
  if (state_ != State::Collecting)
    return;
  if (!blockStarted_) {
    blockStarted_ = true;
    blockStart_ = tNs;
    return;
  }
  if (tNs < blockStart_) {
    blockStart_ = tNs;
    return;
  }
  if (tNs - blockStart_ < cfg_.blockNs)
    return;
  if (tNs - blockStart_ > 2 * cfg_.blockNs) { // data gap: the window is not continuous
    restarts_++;
    restartWindow();
    blockStart_ = tNs;
    return;
  }
  closeBlock();
  blockStart_ += cfg_.blockNs;
}

void GyroBiasCollector::closeBlock() {
  bool ok = block_.nGyro >= 2 && block_.nAccel >= 2;
  if (block_.nGyro >= 2 && block_.nAccel >= 2) { // diagnostics: the worst values seen in any closed block
    diag_.blocks++;
    double mgs = 0, mgp = 0, mas = 0;
    for (int k = 0; k < 3; k++) {
      if (block_.gyroStd(k) > mgs) mgs = block_.gyroStd(k);
      if (block_.gyroP2P(k) > mgp) mgp = block_.gyroP2P(k);
      if (block_.accelStd(k) > mas) mas = block_.accelStd(k);
    }
    if (mgs < diag_.minGyroStd) diag_.minGyroStd = mgs;
    if (mgp < diag_.minGyroP2P) diag_.minGyroP2P = mgp;
    if (mas < diag_.minAccelAxisStd) diag_.minAccelAxisStd = mas;
    if (block_.accelNormStd() < diag_.minAccelNormStd) diag_.minAccelNormStd = block_.accelNormStd();
    for (int k = 0; k < 3; k++) {
      if (block_.gyroStd(k) > diag_.gyroStd)
        diag_.gyroStd = block_.gyroStd(k);
      if (block_.gyroP2P(k) > diag_.gyroP2P)
        diag_.gyroP2P = block_.gyroP2P(k);
      if (block_.accelStd(k) > diag_.accelAxisStd)
        diag_.accelAxisStd = block_.accelStd(k);
    }
    if (block_.accelNormStd() > diag_.accelNormStd)
      diag_.accelNormStd = block_.accelNormStd();
  } else {
    diag_.shortBlocks++;
  }
  // Limit check. Every limit is evaluated (not first-fail) so the worst violation by ratio is recorded.
  GyroViolation blockViol;
  const auto check = [&](GyroViolation::Kind kind, int axis, double value, double limit) {
    if (value < limit)
      return;
    ok = false;
    const double ratio = limit > 0 ? value / limit : 1e300;
    if (ratio > blockViol.ratio || !(ratio == ratio)) {
      blockViol.kind = kind;
      blockViol.axis = axis;
      blockViol.value = value;
      blockViol.limit = limit;
      blockViol.ratio = ratio;
    }
  };
  if (block_.nGyro < 2 || block_.nAccel < 2) {
    ok = false;
    blockViol.kind = GyroViolation::ShortBlock;
    blockViol.value = double(block_.nGyro < block_.nAccel ? block_.nGyro : block_.nAccel);
    blockViol.limit = 2;
    blockViol.ratio = 1;
  } else {
    for (int k = 0; k < 3; k++) {
      check(GyroViolation::GyroStd, k, block_.gyroStd(k), cfg_.maxBlockStd);
      check(GyroViolation::GyroP2P, k, block_.gyroP2P(k), cfg_.maxBlockP2P);
      if (cfg_.maxAccelAxisStd > 0)
        check(GyroViolation::AccelAxisStd, k, block_.accelStd(k), cfg_.maxAccelAxisStd);
    }
    check(GyroViolation::AccelNormStd, -1, block_.accelNormStd(), cfg_.maxAccelNormStd);
  }

  if (block_.nGyro >= 2 && block_.nAccel >= 2) {
    last_.valid = true;
    last_.accepted = ok;
    last_.gyroStd = last_.gyroP2P = last_.accelAxisStd = 0;
    for (int k = 0; k < 3; k++) {
      if (block_.gyroStd(k) > last_.gyroStd) last_.gyroStd = block_.gyroStd(k);
      if (block_.gyroP2P(k) > last_.gyroP2P) last_.gyroP2P = block_.gyroP2P(k);
      if (block_.accelStd(k) > last_.accelAxisStd) last_.accelAxisStd = block_.accelStd(k);
    }
    last_.accelNormStd = block_.accelNormStd();
  }

  if (!ok) {
    restarts_++;
    rejected_++;
    if (blockViol.ratio > worstViol_.ratio)
      worstViol_ = blockViol;
    restartWindow();
    return;
  }
  window_.merge(block_);
  for (int k = 0; k < 3; k++)
    winMagSum_[k] += blockMagSum_[k];
  winMagN_ += blockMagN_;
  if (blockTN_ > 0) {
    if (blockTMin_ < winTMin_)
      winTMin_ = blockTMin_;
    if (blockTMax_ > winTMax_)
      winTMax_ = blockTMax_;
    winTSum_ += blockTSum_;
    winTN_ += blockTN_;
  }
  windowCollectedNs_ += cfg_.blockNs;
  if (cfg_.blockNs > 0 && uint32_t(windowCollectedNs_ / cfg_.blockNs) > bestRun_)
    bestRun_ = uint32_t(windowCollectedNs_ / cfg_.blockNs);
  block_.clear();
  for (int k = 0; k < 3; k++)
    blockMagSum_[k] = 0;
  blockMagN_ = 0;
  blockTMin_ = 1e300;
  blockTMax_ = -1e300;
  blockTSum_ = 0;
  blockTN_ = 0;
}

void GyroBiasCollector::addGyro(int64_t tNs, const double w[3]) {
  if (state_ != State::Collecting || !isfinite(w[0]) || !isfinite(w[1]) || !isfinite(w[2]))
    return;
  advance(tNs);
  block_.addGyro(w, kZero);
}

void GyroBiasCollector::addAccel(int64_t tNs, const double a[3]) {
  if (state_ != State::Collecting || !isfinite(a[0]) || !isfinite(a[1]) || !isfinite(a[2]))
    return;
  advance(tNs);
  block_.addAccel(a);
}

void GyroBiasCollector::addMag(int64_t tNs, const double m[3]) {
  if (state_ != State::Collecting || !isfinite(m[0]) || !isfinite(m[1]) || !isfinite(m[2]))
    return;
  (void)tNs;
  for (int k = 0; k < 3; k++)
    blockMagSum_[k] += m[k];
  blockMagN_++;
}

void GyroBiasCollector::addTemperature(int64_t tNs, double temp) {
  if (state_ != State::Collecting || !isfinite(temp))
    return;
  (void)tNs;
  if (temp < blockTMin_)
    blockTMin_ = temp;
  if (temp > blockTMax_)
    blockTMax_ = temp;
  blockTSum_ += temp;
  blockTN_++;
}

void GyroBiasCollector::finish(int64_t nowNs, bool complete) {
  (void)nowNs;
  state_ = State::Done;
  res_ = GyroBiasResult();
  res_.restarts = restarts_;
  res_.bestRun = bestRun_;
  res_.rejectedBlocks = rejected_;
  res_.worstViolation = worstViol_;
  res_.worst = diag_;
  if (!complete) {
    res_.result = restarts_ > 0 ? CalResult::Motion : CalResult::Timeout;
    return;
  }
  res_.checks |= kGyroCheckComplete;
  res_.nGyro = window_.nGyro;
  res_.nAccel = window_.nAccel;
  res_.nMag = winMagN_;
  bool biasOk = true, noiseOk = true;
  for (int k = 0; k < 3; k++) {
    res_.bias[k] = window_.gyroMean(k);
    res_.noiseStd[k] = window_.gyroStd(k);
    res_.accelMean[k] = window_.accelMean(k);
    res_.magMean[k] = winMagN_ > 0 ? winMagSum_[k] / winMagN_ : 0;
    if (!(fabs(res_.bias[k]) < cfg_.maxBias))
      biasOk = false;
    if (!(res_.noiseStd[k] < cfg_.maxNoiseStd))
      noiseOk = false;
  }
  res_.haveTemp = winTN_ > 0;
  if (res_.haveTemp) {
    res_.tempDelta = winTMax_ - winTMin_;
    res_.tempMean = winTSum_ / winTN_;
  }
  const bool tempOk = !res_.haveTemp || res_.tempDelta < cfg_.maxTempDelta;
  if (biasOk)
    res_.checks |= kGyroCheckBias;
  if (noiseOk)
    res_.checks |= kGyroCheckNoise;
  if (tempOk)
    res_.checks |= kGyroCheckTemp;
  if (!biasOk)
    res_.result = CalResult::FailRange;
  else if (!noiseOk || !tempOk)
    res_.result = CalResult::FailResidual;
  else
    res_.result = CalResult::Pass;
}

int describeGyroRun(const GyroBiasResult &r, const GyroBiasConfig &cfg, char *buf, size_t n) {
  const uint32_t need = cfg.blockNs > 0 ? uint32_t(cfg.windowNs / cfg.blockNs) : 0;
  const GyroViolation &v = r.worstViolation;
  char head[64];
  snprintf(head, sizeof head, "best %lu of %lu blocks, %lu rejected", (unsigned long)r.bestRun, (unsigned long)need,
           (unsigned long)r.rejectedBlocks);
  if (v.kind == GyroViolation::None)
    return snprintf(buf, n, "%s", head);
  static const char *const axisName = "XYZ";
  char axis[4] = {0, 0, 0, 0};
  if (v.axis >= 0 && v.axis < 3) {
    axis[0] = ' ';
    axis[1] = axisName[v.axis];
  }
  const double dg = 180.0 / 3.14159265358979;
  switch (v.kind) {
  case GyroViolation::GyroStd:
    return snprintf(buf, n, "%s, worst: gyro std%s %.2f deg/s > limit %.2f deg/s", head, axis, v.value * dg, v.limit * dg);
  case GyroViolation::GyroP2P:
    return snprintf(buf, n, "%s, worst: gyro p2p%s %.2f deg/s > limit %.2f deg/s", head, axis, v.value * dg, v.limit * dg);
  case GyroViolation::AccelAxisStd:
    return snprintf(buf, n, "%s, worst: accel std%s %.3f m/s^2 > limit %.3f m/s^2", head, axis, v.value, v.limit);
  case GyroViolation::AccelNormStd:
    return snprintf(buf, n, "%s, worst: |a| std %.3f m/s^2 > limit %.3f m/s^2", head, v.value, v.limit);
  default:
    return snprintf(buf, n, "%s, worst: block with only %.0f gyro/accel samples (need 2)", head, v.value);
  }
}

bool GyroBiasCollector::poll(int64_t nowNs) {
  if (state_ != State::Collecting)
    return false;
  if (windowCollectedNs_ >= cfg_.windowNs) {
    finish(nowNs, true);
    return true;
  }
  if (nowNs >= deadline_) {
    finish(nowNs, false);
    return true;
  }
  return false;
}

} // namespace Calib
} // namespace DSP
} // namespace VCTR
