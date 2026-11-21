#ifndef EXVECTRDSP_MAG_SPHERE_FIT_HPP
#define EXVECTRDSP_MAG_SPHERE_FIT_HPP

/**
 * Magnetometer hard-iron (bias) calibration by a least-squares SPHERE fit (owner decision 2026-11-21: bias only, no
 * soft-iron ellipsoid; workspace docs/calibration-design.md section 6.1).
 *
 * Algorithm (the vehicle's former MagnetometerCalibrationTask, moved here unchanged in its maths but accumulated in
 * double): for every sample m, |m|^2 = 2 c . m + (R^2 - |c|^2) is linear in beta = [2c, R^2 - |c|^2]; the 4x4 normal
 * matrix sum(m4 m4^T) and sum(m4 |m|^2), m4 = [m, 1], are accumulated online (no sample limit), and solved for the
 * centre c (hard-iron offset, raw units) and the radius R. The calibration is
 *
 *   b = c,   M_cal = s * I with s = rRef / R   (|B_cal| = rRef, 0.05 driver units, keeps the attitude EKF gate)
 *
 * Quality (added 2026-11-21): a direction-thinned copy of the samples (a sample is stored only if its direction about
 * the current centre differs by more than storeAngle from the last stored one) is kept in a caller-provided buffer.
 * From it: the residual of |s (m - c)| / rRef - 1 (RMS and max), and the direction coverage over the 64 equal-area
 * bins of SphereCoverage. Checks: coverage, stored sample count, residual RMS and max, hard-iron offset |c| < ratio * R.
 *
 * Pure code: no Arduino, no ExVectrCore. Host-tested (flight-software test/test_cal_mag_fit).
 */

#include <stddef.h>
#include <stdint.h>

#include "ExVectrDSP/calibration/calibration_result.hpp"
#include "ExVectrDSP/calibration/sphere_coverage.hpp"

namespace VCTR {
namespace DSP {
namespace Calib {

struct MagFitConfig {
  double rRef = 0.05;           ///< calibrated field radius, driver units (decision D2)
  double minCoverage = 0.75;    ///< fraction of the 64 bins with at least minBinCount stored samples
  uint8_t minBinCount = 2;
  size_t minStored = 300;       ///< stored (direction-thinned) samples
  double maxResidualRms = 0.03; ///< of |B_cal| / rRef - 1 (sphere fit: soft iron and axis mismatch remain)
  double maxResidualMax = 0.10;
  double maxOffsetRatio = 0.8;  ///< |c| < ratio * R
  double storeAngleRad = 5.0 * 3.14159265358979 / 180.0;
};

/// Pass bits of MagFitResult::checks.
enum MagCheck : uint16_t {
  kMagCheckSolved = 1u << 0,
  kMagCheckSamples = 1u << 1,
  kMagCheckCoverage = 1u << 2,
  kMagCheckResidualRms = 1u << 3,
  kMagCheckResidualMax = 1u << 4,
  kMagCheckOffset = 1u << 5,
};

struct MagFitResult {
  CalResult result = CalResult::None;
  uint16_t checks = 0;
  double offset[3] = {0, 0, 0}; ///< hard-iron centre c, raw units
  double radius = 0;            ///< fitted raw field radius R
  double scale = 1;             ///< s = rRef / R
  size_t nAccumulated = 0;      ///< samples in the normal matrix
  size_t nStored = 0;           ///< direction-thinned samples
  double coverage = 0;          ///< fraction of bins with >= minBinCount stored samples
  double residualRms = 0;       ///< of |B_cal|/rRef - 1
  double residualMax = 0;
  double rawSpread = 0;         ///< std/mean of the raw |m| over the stored samples (before calibration)
};

class MagSphereFit {
public:
  /// buffer: storage for up to capacity stored samples (the vehicle puts it in DMAMEM).
  MagSphereFit(float (*buffer)[3], size_t capacity, const MagFitConfig &cfg = MagFitConfig());

  void reset();

  /// Adds one raw sample (raw units). Non-finite samples are ignored.
  void addSample(const double m[3]);

  /// Solves the normal equations with the samples so far and updates the centre used for thinning and coverage.
  /// Returns false if not solvable yet (fewer than 4 samples, singular, or R^2 <= 0).
  bool updateEstimate();

  /// Final fit and checks over everything collected (also usable for a live candidate).
  MagFitResult evaluate();

  size_t storedCount() const { return nStored_; }
  size_t accumulatedCount() const { return nAcc_; }
  /// Coverage of the stored samples about the current centre (recomputed by updateEstimate / evaluate).
  const SphereCoverage &coverage() const { return coverage_; }
  const MagFitConfig &config() const { return cfg_; }

  /// Least-squares sphere through n points: centre and radius. Returns false if not solvable.
  static bool solveSphere(const double *xTx /*4x4*/, const double *xTy /*4*/, double centre[3], double &radius);

private:
  void recomputeCoverage();

  float (*buf_)[3];
  size_t cap_;
  MagFitConfig cfg_;
  double xTx_[16];
  double xTy_[4];
  size_t nAcc_ = 0;
  size_t nStored_ = 0;
  double centre_[3] = {0, 0, 0};
  double radius_ = 0;
  bool haveEstimate_ = false;
  double lastStoredDir_[3] = {0, 0, 0};
  bool haveLast_ = false;
  SphereCoverage coverage_;
};

} // namespace Calib
} // namespace DSP
} // namespace VCTR

#endif
