#ifndef EXVECTRDSP_GYRO_BIAS_COLLECTOR_HPP
#define EXVECTRDSP_GYRO_BIAS_COLLECTOR_HPP

/**
 * Stationary gyro bias measurement (workspace docs/calibration-design.md section 2.3).
 *
 * The bias is the mean of the raw gyro over a still window of windowNs. Samples are checked in blocks of blockNs
 * (1 s): a block with any gyro axis std > maxBlockStd, or peak-to-peak > maxBlockP2P, or |a| std > maxAccelNormStd,
 * or any accel axis std > maxAccelAxisStd restarts the window (the samples so far are dropped, restarts are
 * counted), so the run is accepted as soon as windowNs/blockNs CONSECUTIVE good blocks occur anywhere before the
 * deadline given to start() (isolated bad blocks cost only the blocks collected before them). Otherwise the result is
 * Motion (if a restart happened) or Timeout. The worst violation is kept in the result (worstViolation).
 *
 * Accept: |bias| < maxBias per axis (FailRange), gyro noise std per axis < maxNoiseStd (FailResidual), IMU
 * temperature change over the window < maxTempDelta (FailResidual; only if temperatures were fed).
 *
 * The mean raw accelerometer and magnetometer over the same window are kept (zero motion), so the caller can check
 * the active accel/mag calibrations against gravity and the field radius.
 *
 * Pure code: no Arduino, no ExVectrCore. Host-tested (flight-software test/test_cal_gyro).
 */

#include <stddef.h>
#include <stdint.h>

#include "ExVectrDSP/calibration/calibration_result.hpp"
#include "ExVectrDSP/calibration/stillness_detector.hpp"

namespace VCTR {
namespace DSP {
namespace Calib {

/**
 * Limits of the 1 s still blocks (owner decision 2026-11-22: the original 0.25 deg/s / 1 deg/s / 0.05 m/s^2 limits
 * rejected every block on the owner's bench, so they were loosened a lot). The block check only has to tell "resting
 * on the ground, hand-touched or on a vibrating table" from "being carried or rotated": the residual gyro bias error
 * of a slightly noisy window is tracked by the attitude estimator anyway, and a bias beyond maxBias still fails.
 * The accelerometer limits stay as the guard against a slowly rotating or handled ship (a steady rotation has a
 * low gyro std but turns the gravity vector): 10 deg/s gives about 0.5 m/s^2 std, 5 deg/s about 0.25.
 */
/// Still-window length per run kind (named, not magic). Every judging limit below is the same for both; only the
/// number of consecutive good 1 s blocks that must occur differs (a failed block restarts the window, never relaxed).
constexpr int64_t kGyroWindowStartupNs = 5000000000LL;   ///< boot zero, ZeroEstimators, manual runs: 5 good blocks
constexpr int64_t kGyroWindowCountdownNs = 3000000000LL; ///< the go countdown (Initialisation phase of 5 s): 3 good blocks

struct GyroBiasConfig {
  int64_t windowNs = kGyroWindowStartupNs;  ///< still data (= windowNs / blockNs consecutive good blocks anywhere before the deadline)
  int64_t blockNs = 1000000000LL;   ///< 1 s motion check blocks
  double maxBlockStd = 1.5 * 3.14159265358979 / 180.0;    ///< rad/s per gyro axis
  double maxBlockP2P = 8.0 * 3.14159265358979 / 180.0;    ///< rad/s per gyro axis (white noise of std s has p2p about 7.5 s per 1000 samples: 5 deg/s failed the owner-like 0.7 deg/s)
  double maxAccelNormStd = 0.3;     ///< m/s^2
  double maxAccelAxisStd = 0.3;     ///< m/s^2 per axis
  double maxBias = 5.0 * 3.14159265358979 / 180.0;        ///< rad/s per axis (MPU9250 ZRO +-5 deg/s)
  double maxNoiseStd = 1.5 * 3.14159265358979 / 180.0;    ///< rad/s per axis (same as the block limit; the window std
                                                           ///< of accepted blocks is below it up to the bias drift)
  double maxTempDelta = 0.5;        ///< K
};

/// Pass bits of GyroBiasResult::checks.
enum GyroCheck : uint16_t {
  kGyroCheckComplete = 1u << 0, ///< a full still window was collected before the deadline
  kGyroCheckBias = 1u << 1,
  kGyroCheckNoise = 1u << 2,
  kGyroCheckTemp = 1u << 3,
};

/// Worst values seen in any closed 1 s block of a run (diagnostics for tuning the limits; rad/s and m/s^2).
struct GyroBlockDiagnostics {
  double gyroStd = 0;       ///< largest per-axis gyro std
  double gyroP2P = 0;       ///< largest per-axis gyro peak-to-peak
  double accelAxisStd = 0;  ///< largest per-axis accel std
  double accelNormStd = 0;  ///< largest std of |a|
  uint32_t shortBlocks = 0; ///< blocks with fewer than 2 gyro or accel samples
  uint32_t blocks = 0;      ///< closed blocks with enough samples
  // The quietest block seen (the sensor's own noise floor when the vehicle is really still).
  double minGyroStd = 1e300;
  double minGyroP2P = 1e300;
  double minAccelAxisStd = 1e300;
  double minAccelNormStd = 1e300;
};

/// Which block limit was violated (the worst violation of a run, for the failure log and for tuning).
struct GyroViolation {
  enum Kind : uint8_t { None, GyroStd, GyroP2P, AccelAxisStd, AccelNormStd, ShortBlock };
  Kind kind = None;
  int axis = -1;        ///< 0..2 (x,y,z) for the per-axis limits, else -1
  double value = 0;     ///< measured (rad/s, m/s^2; sample count for ShortBlock)
  double limit = 0;     ///< the limit it exceeded
  double ratio = 0;     ///< value / limit, the key for "worst"
};

/// The last closed 1 s block (live feedback while a run collects).
struct GyroLastBlock {
  bool valid = false;
  bool accepted = false;
  double gyroStd = 0, gyroP2P = 0, accelAxisStd = 0, accelNormStd = 0; ///< worst axis / of |a|
};

struct GyroBiasResult {
  CalResult result = CalResult::None;
  uint16_t checks = 0;
  double bias[3] = {0, 0, 0};   ///< raw sensor frame, rad/s
  double noiseStd[3] = {0, 0, 0};
  double tempDelta = 0;         ///< K (0 if no temperature was fed)
  double tempMean = 0;          ///< K or degC as fed
  bool haveTemp = false;
  double accelMean[3] = {0, 0, 0}; ///< raw accel over the window, m/s^2
  uint32_t nAccel = 0;
  double magMean[3] = {0, 0, 0};   ///< raw mag over the window
  uint32_t nMag = 0;
  uint32_t restarts = 0;
  uint32_t nGyro = 0;
  uint32_t bestRun = 0;         ///< most consecutive good blocks seen (needed: windowNs / blockNs)
  uint32_t rejectedBlocks = 0;
  GyroViolation worstViolation; ///< the rejected block that exceeded its limit by the largest factor
  GyroBlockDiagnostics worst; ///< over all blocks of the run, accepted or not
};

/// One-line description of a run's blocks for the log, e.g. "best 3 of 5 blocks, 4 rejected, worst: gyro std X 2.10
/// deg/s > limit 1.50 deg/s". Returns the snprintf result.
int describeGyroRun(const GyroBiasResult &r, const GyroBiasConfig &cfg, char *buf, size_t n);

class GyroBiasCollector {
public:
  enum class State : uint8_t { Idle, Collecting, Done };

  explicit GyroBiasCollector(const GyroBiasConfig &cfg = GyroBiasConfig()) : cfg_(cfg), baseWindowNs_(cfg.windowNs) {}

  /// Starts a new measurement; it must complete by deadlineNs. `windowNs` > 0 overrides the still-window length for
  /// this run only (kGyroWindowCountdownNs for the go countdown); 0 = the configured window. All limits unchanged.
  void start(int64_t nowNs, int64_t deadlineNs, int64_t windowNs = 0);
  void cancel() { state_ = State::Idle; }

  void addGyro(int64_t tNs, const double w[3]);
  void addAccel(int64_t tNs, const double a[3]);
  void addMag(int64_t tNs, const double m[3]);
  void addTemperature(int64_t tNs, double temp);

  /// Advances time (checks the deadline). Returns true once, when the state becomes Done.
  bool poll(int64_t nowNs);

  State state() const { return state_; }
  bool collecting() const { return state_ == State::Collecting; }
  const GyroBiasResult &result() const { return res_; }
  const GyroLastBlock &lastBlock() const { return last_; }
  /// Fraction (0..1) of the window collected so far.
  double progress() const;
  uint32_t restarts() const { return restarts_; }
  const GyroBiasConfig &config() const { return cfg_; }

private:
  void advance(int64_t tNs);
  void closeBlock();
  void restartWindow();
  void finish(int64_t nowNs, bool complete);

  GyroBiasConfig cfg_;
  int64_t baseWindowNs_ = kGyroWindowStartupNs;
  State state_ = State::Idle;
  int64_t deadline_ = 0;
  int64_t blockStart_ = 0;
  bool blockStarted_ = false;
  MotionStats block_, window_;
  double blockMagSum_[3] = {0, 0, 0}, winMagSum_[3] = {0, 0, 0};
  uint32_t blockMagN_ = 0, winMagN_ = 0;
  double blockTMin_ = 1e300, blockTMax_ = -1e300, blockTSum_ = 0;
  uint32_t blockTN_ = 0;
  double winTMin_ = 1e300, winTMax_ = -1e300, winTSum_ = 0;
  uint32_t winTN_ = 0;
  int64_t windowCollectedNs_ = 0;
  uint32_t restarts_ = 0;
  uint32_t bestRun_ = 0, rejected_ = 0;
  GyroViolation worstViol_;
  GyroBiasResult res_;
  GyroBlockDiagnostics diag_;
  GyroLastBlock last_;
};

} // namespace Calib
} // namespace DSP
} // namespace VCTR

#endif
