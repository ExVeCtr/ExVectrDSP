#ifndef EXVECTRDSP_STILLNESS_DETECTOR_HPP
#define EXVECTRDSP_STILLNESS_DETECTOR_HPP

/**
 * Stillness detection from gyro and accelerometer samples (workspace docs/calibration-design.md sections 2.3, 2.4).
 *
 * Samples are grouped into time blocks (blockNs, by sample timestamp). The detector reports still() when the last
 * windowBlocks complete blocks are present without gaps and the statistics over that window are inside every
 * enabled limit (a limit <= 0 is disabled):
 *   - maxGyroNorm: largest |w - gyroBias| of any sample, rad/s (calibrated rate; the bias is set by the caller)
 *   - maxGyroStd: per-axis standard deviation of the gyro, rad/s
 *   - maxGyroP2P: per-axis peak-to-peak of the gyro, rad/s
 *   - maxAccelStd: per-axis standard deviation of the accelerometer, m/s^2
 *   - maxAccelNormStd: standard deviation of |a|, m/s^2
 * Defaults: the accelerometer pose capture rule (still for 1.5 s: |w_cal| < 0.02 rad/s, accel per-axis std
 * < 0.04 m/s^2).
 *
 * A gap of more than two blocks without samples empties the window. Time is int64 ns, injected.
 * Pure code: no Arduino, no ExVectrCore. Host-tested (flight-software test/test_cal_stillness_pose).
 */

#include <stddef.h>
#include <stdint.h>

namespace VCTR {
namespace DSP {
namespace Calib {

struct StillnessConfig {
  int64_t blockNs = 100000000LL; ///< 0.1 s
  uint8_t windowBlocks = 15;     ///< 1.5 s
  double maxGyroNorm = 0.02;     ///< rad/s
  double maxGyroStd = 0;         ///< disabled
  double maxGyroP2P = 0;         ///< disabled
  double maxAccelStd = 0.04;     ///< m/s^2
  double maxAccelNormStd = 0;    ///< disabled
};

/// Sums of one block or window.
struct MotionStats {
  uint32_t nGyro = 0, nAccel = 0;
  double gSum[3] = {0, 0, 0}, gSum2[3] = {0, 0, 0};
  double gMin[3] = {1e300, 1e300, 1e300}, gMax[3] = {-1e300, -1e300, -1e300};
  double gNormMax = 0;
  double aSum[3] = {0, 0, 0}, aSum2[3] = {0, 0, 0};
  double aNormSum = 0, aNormSum2 = 0;

  void clear() { *this = MotionStats(); }
  void addGyro(const double w[3], const double bias[3]);
  void addAccel(const double a[3]);
  void merge(const MotionStats &o);
  double gyroStd(int k) const;
  double gyroP2P(int k) const { return nGyro > 0 ? gMax[k] - gMin[k] : 0; }
  double accelStd(int k) const;
  double accelNormStd() const;
  double gyroMean(int k) const { return nGyro > 0 ? gSum[k] / nGyro : 0; }
  double accelMean(int k) const { return nAccel > 0 ? aSum[k] / nAccel : 0; }
};

class StillnessDetector {
public:
  static constexpr uint8_t kMaxBlocks = 16;

  explicit StillnessDetector(const StillnessConfig &cfg = StillnessConfig());

  void reset();
  void setGyroBias(const double bias[3]);

  void addGyro(int64_t tNs, const double w[3]);
  void addAccel(int64_t tNs, const double a[3]);

  bool still() const { return still_; }
  /// Largest per-axis gyro / accel std and largest |w_cal| over the current window (diagnostics).
  double windowGyroNormMax() const { return winGyroNormMax_; }
  double windowAccelStdMax() const { return winAccelStdMax_; }

private:
  void advance(int64_t tNs);
  void closeBlock();
  void evaluate();

  StillnessConfig cfg_;
  double bias_[3] = {0, 0, 0};
  MotionStats blocks_[kMaxBlocks];
  uint8_t head_ = 0, count_ = 0;
  MotionStats cur_;
  int64_t blockStart_ = 0;
  bool started_ = false;
  bool still_ = false;
  double winGyroNormMax_ = 0, winAccelStdMax_ = 0;
};

} // namespace Calib
} // namespace DSP
} // namespace VCTR

#endif
