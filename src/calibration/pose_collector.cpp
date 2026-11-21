#include "ExVectrDSP/calibration/pose_collector.hpp"

#include <math.h>

namespace VCTR {
namespace DSP {
namespace Calib {

const char *toString(PoseEvent e) {
  switch (e) {
  case PoseEvent::None: return "none";
  case PoseEvent::Averaging: return "averaging";
  case PoseEvent::Captured: return "captured";
  case PoseEvent::Duplicate: return "duplicate";
  case PoseEvent::Aborted: return "aborted";
  case PoseEvent::Full: return "full";
  }
  return "?";
}

void PoseCollector::reset() {
  state_ = State::WaitStill;
  avgN_ = 0;
  n_ = 0;
  lastNearest_ = 0;
  for (size_t i = 0; i < kMaxPoses; i++)
    for (int k = 0; k < 3; k++)
      mean_[i][k] = std_[i][k] = 0;
}

PoseEvent PoseCollector::addAccel(int64_t tNs, const double a[3], bool still) {
  if (!isfinite(a[0]) || !isfinite(a[1]) || !isfinite(a[2]))
    return PoseEvent::None;

  switch (state_) {
  case State::WaitMotion:
    if (!still)
      state_ = State::WaitStill;
    return PoseEvent::None;

  case State::WaitStill:
    if (!still)
      return PoseEvent::None;
    state_ = State::Averaging;
    avgStart_ = tNs;
    avgN_ = 0;
    for (int k = 0; k < 3; k++)
      avgSum_[k] = avgSum2_[k] = 0;
    // fall through: this sample is the first of the average
    for (int k = 0; k < 3; k++) {
      avgSum_[k] += a[k];
      avgSum2_[k] += a[k] * a[k];
    }
    avgN_++;
    return PoseEvent::Averaging;

  case State::Averaging:
    break;
  }

  if (!still) {
    state_ = State::WaitStill;
    return PoseEvent::Aborted;
  }
  for (int k = 0; k < 3; k++) {
    avgSum_[k] += a[k];
    avgSum2_[k] += a[k] * a[k];
  }
  avgN_++;
  if (tNs - avgStart_ < cfg_.averageNs)
    return PoseEvent::None;

  // Average complete.
  state_ = State::WaitMotion;
  double m[3], s[3];
  for (int k = 0; k < 3; k++) {
    m[k] = avgSum_[k] / avgN_;
    const double var = avgSum2_[k] / avgN_ - m[k] * m[k];
    s[k] = var > 0 ? sqrt(var) : 0;
    lastAvg_[k] = m[k];
  }
  const double mn = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
  double nearest = INFINITY;
  for (size_t i = 0; i < n_; i++) {
    const double pn = sqrt(mean_[i][0] * mean_[i][0] + mean_[i][1] * mean_[i][1] + mean_[i][2] * mean_[i][2]);
    double c = (m[0] * mean_[i][0] + m[1] * mean_[i][1] + m[2] * mean_[i][2]) / (mn * pn);
    c = c > 1 ? 1 : (c < -1 ? -1 : c);
    const double ang = acos(c);
    if (ang < nearest)
      nearest = ang;
  }
  lastNearest_ = nearest;
  if (!(mn > 1e-9))
    return PoseEvent::Duplicate;
  if (n_ > 0 && !(nearest > cfg_.newPoseAngleRad))
    return PoseEvent::Duplicate;
  if (n_ >= kMaxPoses)
    return PoseEvent::Full;
  for (int k = 0; k < 3; k++) {
    mean_[n_][k] = m[k];
    std_[n_][k] = s[k];
  }
  n_++;
  return PoseEvent::Captured;
}

} // namespace Calib
} // namespace DSP
} // namespace VCTR
