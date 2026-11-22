#ifndef EXVECTRDSP_MAG_ELLIPSOID_FIT_HPP
#define EXVECTRDSP_MAG_ELLIPSOID_FIT_HPP

/**
 * Magnetometer hard-iron + soft-iron calibration by an ELLIPSOID fit (owner decision 2026-11-21, supersedes the
 * "bias only" decision 6.1 of the workspace docs/calibration-design.md; method: section 2.5).
 *
 * Model: m_cal = A (x - b), offset b (hard iron, raw driver units) and a SYMMETRIC positive-definite 3x3 matrix A
 * (soft iron: axis scale and non-orthogonality), scaled so that the calibrated field has the constant length rRef
 * (0.05 driver units, keeps the attitude filter gate 0.01..0.1). The rotation part of A is not observable from field
 * lengths and stays in the mounting matrix T (vehicle side). The attitude filter uses the field direction only, so a
 * soft-iron distortion that is not removed makes tipping the ship look like a rotation.
 *
 * Algorithm:
 *  1. Online accumulation (all samples, no sample limit, double precision): the data are normalised first (shifted by
 *     the first sample, divided by normScale) and the normal matrix of the algebraic quadric fit is accumulated.
 *     Quadric x'Qx + p'x + c = 0 with the linear constraint trace(Q) = 1 (a sphere-like normalisation that cannot be
 *     degenerate for ellipsoid data and is shift invariant), solved as a 9x9 linear least-squares problem. Q, p give
 *     the centre b = -Q^-1 p / 2 and, with the matrix square root, A.
 *  2. Gauss-Newton (Levenberg-Marquardt damped) refinement of the 9 parameters (b, A) minimising
 *     sum (|A (x_i - b)| / rRef - 1)^2 on the direction-thinned stored samples (the same buffer as MagSphereFit),
 *     started from step 1 (or from the sphere fit if the quadric is not an ellipsoid).
 *  3. Outlier rejection: outlierPasses (2) passes, each marks samples with |residual| > outlierSigmas (3) x sigma as
 *     disturbed (sigma = 1.4826 x median |residual|, at least noiseFloor) and refits on the inliers. The number of
 *     rejected samples is reported; too many rejected is a failure (maxRejectedFraction), some is a warning.
 *  4. Acceptance from the fit RESULT (not from where the ship was pointed): residual rms and max on the INLIERS,
 *     direction coverage (SphereCoverage, 64 bins), stored sample count, parameter standard errors from the
 *     Gauss-Newton covariance sigma^2 (J'J)^-1 (like the accelerometer fit), and soft-iron plausibility
 *     (eigenvalues of A relative to their mean within [kMagEigMinRel, kMagEigMaxRel], axis ratio max/min below
 *     kMagAxisRatioMax) and hard-iron |b| < maxOffsetRatio R.
 *
 * softIron = false gives the former sphere behaviour (A = s I, s = rRef / R, algebraic sphere fit, no refinement and
 * no outlier rejection) with the same acceptance limits, to compare both on the same data.
 *
 * Pure code: no Arduino, no ExVectrCore. Host-tested (flight-software test/test_cal_mag_ellipsoid).
 */

#include <stddef.h>
#include <stdint.h>

#include "ExVectrDSP/calibration/calibration_result.hpp"
#include "ExVectrDSP/calibration/sphere_coverage.hpp"

namespace VCTR {
namespace DSP {
namespace Calib {

/// Soft-iron plausibility, shared with the vehicle's stored-parameter check (Arm/MissionBegin).
constexpr double kMagAxisRatioMax = 1.5; ///< max/min eigenvalue of A
constexpr double kMagEigMinRel = 0.6;    ///< every eigenvalue of A >= this x their mean
constexpr double kMagEigMaxRel = 1.6;    ///< ... and <= this x their mean

struct MagEllipsoidConfig {
  double rRef = 0.05;            ///< calibrated field length, driver units (decision D2)
  bool softIron = true;          ///< false: sphere (A = s I) for comparison with the old behaviour
  double normScale = 0.3;        ///< raw units divided by this before accumulation (typical field radius)
  double minCoverage = 0.75;     ///< fraction of the 64 bins with at least minBinCount stored samples
  uint8_t minBinCount = 2;
  size_t minStored = 300;        ///< stored (direction-thinned) samples
  double maxResidualRms = 0.035; ///< of |B_cal| / rRef - 1, inliers (3.5 %): above this FailResidual
  double warnResidualRms = 0.015; ///< rms between this (1.5 %) and maxResidualRms: Warn (field noise or disturbance)
  double maxResidualMax = 0.10;  ///< inliers (10 %)
  double maxOffsetRatio = 0.8;   ///< |b| < ratio * R
  double maxAxisRatio = kMagAxisRatioMax;
  double eigMinRel = kMagEigMinRel;
  double eigMaxRel = kMagEigMaxRel;
  double maxOffsetSigma = 0.01;  ///< standard error of each offset component, as a fraction of the field radius R
  double maxMatrixSigma = 0.02;  ///< standard error of each entry of A, as a fraction of the mean eigenvalue of A
  double noiseFloor = 0.002;     ///< lower bound of the residual sigma (relative to rRef)
  double outlierSigmas = 3.0;
  int outlierPasses = 2;
  double warnRejectedFraction = 0.05; ///< more rejected than this: Warn
  double maxRejectedFraction = 0.15;  ///< more rejected than this: FailResidual (the data are mostly disturbed)
  int maxIterations = 40;
  double storeAngleRad = 5.0 * 3.14159265358979 / 180.0;
};

/// Pass bits of MagEllipsoidResult::checks.
enum MagCheck : uint16_t {
  kMagCheckSolved = 1u << 0,
  kMagCheckSamples = 1u << 1,
  kMagCheckCoverage = 1u << 2,
  kMagCheckResidualRms = 1u << 3,
  kMagCheckResidualMax = 1u << 4,
  kMagCheckOffset = 1u << 5,
  kMagCheckWellDetermined = 1u << 6, ///< standard errors within limits (verdict FAIL_COVERAGE when not)
  kMagCheckSoftIron = 1u << 7,       ///< eigenvalues of A and axis ratio plausible (verdict FAIL_RANGE when not)
  kMagCheckRejected = 1u << 8,       ///< rejected fraction <= maxRejectedFraction
  kMagCheckRejectedNominal = 1u << 9, ///< rejected fraction <= warnRejectedFraction (warning only)
  kMagCheckResidualNominal = 1u << 10, ///< residual rms < warnResidualRms (warning only)
};

struct MagEllipsoidResult {
  CalResult result = CalResult::None;
  uint16_t checks = 0;
  bool softIron = true;                      ///< the mode that produced this result
  double offset[3] = {0, 0, 0};              ///< hard-iron b, raw driver units
  double A[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; ///< symmetric, row-major, sensor frame: m_cal = A (x - b), |m_cal| = rRef
  double eig[3] = {1, 1, 1};                 ///< eigenvalues of A, ascending
  double scale = 1;                          ///< mean eigenvalue of A (= s = rRef / R for a sphere)
  double radius = 0;                         ///< mean raw field radius rRef / scale
  double axisRatio = 1;                      ///< max / min eigenvalue of A
  size_t nAccumulated = 0;                   ///< samples in the normal matrix
  size_t nStored = 0;                        ///< direction-thinned samples
  size_t nInliers = 0;                       ///< stored samples used for the final fit
  size_t nRejected = 0;                      ///< stored samples rejected as disturbed
  double coverage = 0;                       ///< fraction of bins with >= minBinCount stored samples
  double residualRms = 0;                    ///< of |B_cal|/rRef - 1 over the inliers
  double residualMax = 0;                    ///< over the inliers
  double rawSpread = 0;                      ///< std/mean of the raw |x| over the stored samples (before calibration)
  double sigmaOffset = 0;                    ///< worst offset standard error, fraction of R (valid if sigmaValid)
  double sigmaMatrix = 0;                    ///< worst entry standard error of A, fraction of the mean eigenvalue
  bool sigmaValid = false;                   ///< the Jacobian was invertible
  int iterations = 0;
};

/// Short text naming the first failed check (for the operator log), or "all checks passed".
const char *magFailReason(const MagEllipsoidResult &r);

class MagEllipsoidFit {
public:
  static constexpr size_t kMaxStored = 600;

  /// buffer: storage for up to capacity (<= kMaxStored) stored samples (the vehicle puts it in DMAMEM).
  MagEllipsoidFit(float (*buffer)[3], size_t capacity, const MagEllipsoidConfig &cfg = MagEllipsoidConfig());

  void reset();

  /// Adds one raw sample (raw units). Non-finite samples are ignored.
  void addSample(const double m[3]);

  /// Refreshes the centre used for thinning and coverage from the samples so far (cheap algebraic solve).
  /// Returns false if not solvable yet.
  bool updateEstimate();

  /// Final fit and checks over everything collected.
  MagEllipsoidResult evaluate();

  size_t storedCount() const { return nStored_; }
  size_t accumulatedCount() const { return nAcc_; }
  const SphereCoverage &coverage() const { return coverage_; }
  const MagEllipsoidConfig &config() const { return cfg_; }

private:
  struct Params { ///< normalised units: b (3) and W symmetric (6), |W (x_n - b)| = 1
    double v[9];
  };
  bool solveAlgebraic(Params &p) const;
  bool solveSphereNormalised(Params &p) const;
  void recomputeCoverage();
  void toNormalised(const float *x, double xn[3]) const;
  void refine(Params &p, const uint8_t *mask, int &iterations) const;

  float (*buf_)[3];
  size_t cap_;
  MagEllipsoidConfig cfg_;
  // Algebraic ellipsoid normal matrix (9x9, regressors [y2-x2, z2-x2, xy, xz, yz, x, y, z, 1], target -x2).
  double eN_[81];
  double eY_[9];
  // Sphere normal matrix (normalised data).
  double sN_[16];
  double sY_[4];
  double shift_[3] = {0, 0, 0};
  bool haveShift_ = false;
  size_t nAcc_ = 0;
  size_t nStored_ = 0;
  double centre_[3] = {0, 0, 0}; ///< raw units, for thinning and coverage
  double lastStoredDir_[3] = {0, 0, 0};
  bool haveLast_ = false;
  SphereCoverage coverage_;
  uint8_t mask_[kMaxStored];
  float absRes_[kMaxStored];
};

/// Calibrated field of a sample: A (x - b), writes out[3]; returns its length.
double magCalibrate(const MagEllipsoidResult &r, const double x[3], double out[3]);

} // namespace Calib
} // namespace DSP
} // namespace VCTR

#endif
