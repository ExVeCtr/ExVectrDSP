#ifndef EXVECTRDSP_POSE_COLLECTOR_HPP
#define EXVECTRDSP_POSE_COLLECTOR_HPP

/**
 * Automatic static pose capture for the accelerometer fit (workspace docs/calibration-design.md section 2.4).
 *
 * Fed with accelerometer samples and the still flag of a StillnessDetector (which already requires 1.5 s of
 * stillness). When still, the samples are averaged for averageNs (2 s); if the vehicle moves during the average the
 * pose is dropped (Aborted). A completed average counts as a NEW pose if its direction differs by more than
 * newPoseAngle (20 degrees) from every pose captured so far; otherwise it is a Duplicate and not stored. After any
 * completed average the collector waits for motion (still flag false) before it averages again, so a vehicle left
 * lying produces one event, not a stream.
 *
 * Pure code: no Arduino, no ExVectrCore. Host-tested (flight-software test/test_cal_stillness_pose).
 */

#include <stddef.h>
#include <stdint.h>

namespace VCTR {
namespace DSP {
namespace Calib {

struct PoseConfig {
  int64_t averageNs = 2000000000LL;                       ///< 2 s
  double newPoseAngleRad = 20.0 * 3.14159265358979 / 180.0; ///< 20 degrees
};

enum class PoseEvent : uint8_t {
  None = 0,
  Averaging, ///< stillness reached, averaging started
  Captured,  ///< new pose stored
  Duplicate, ///< average done, but within newPoseAngle of a stored pose: not stored
  Aborted,   ///< motion during the average: dropped
  Full,      ///< average done, but kMaxPoses are stored already
};

const char *toString(PoseEvent e);

class PoseCollector {
public:
  static constexpr size_t kMaxPoses = 24;

  explicit PoseCollector(const PoseConfig &cfg = PoseConfig()) : cfg_(cfg) { reset(); }

  void reset();

  /// One accelerometer sample (raw sensor frame, m/s^2) with the current still flag.
  PoseEvent addAccel(int64_t tNs, const double a[3], bool still);

  size_t count() const { return n_; }
  const double (*means() const)[3] { return mean_; }
  const double *mean(size_t i) const { return mean_[i]; }
  const double *stdDev(size_t i) const { return std_[i]; }
  bool averaging() const { return state_ == State::Averaging; }
  /// Angle (rad) of the last completed average to the nearest stored pose (diagnostics; inf if none stored).
  double lastNearestAngle() const { return lastNearest_; }
  const double *lastAverage() const { return lastAvg_; }

private:
  enum class State : uint8_t { WaitStill, Averaging, WaitMotion };

  PoseConfig cfg_;
  State state_ = State::WaitStill;
  int64_t avgStart_ = 0;
  uint32_t avgN_ = 0;
  double avgSum_[3] = {0, 0, 0}, avgSum2_[3] = {0, 0, 0};
  size_t n_ = 0;
  double mean_[kMaxPoses][3];
  double std_[kMaxPoses][3];
  double lastAvg_[3] = {0, 0, 0};
  double lastNearest_ = 0;
};

} // namespace Calib
} // namespace DSP
} // namespace VCTR

#endif
