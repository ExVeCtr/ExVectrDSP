#ifndef EXVECTRDSP_ACCEL_ELLIPSOID_FIT_HPP
#define EXVECTRDSP_ACCEL_ELLIPSOID_FIT_HPP

/**
 * Accelerometer multi-pose magnitude fit (workspace docs/calibration-design.md section 2.4, owner decision D3 = A).
 *
 * Model: a_cal = A * (x - b), with offset b (3) and a SYMMETRIC 3x3 matrix A (scale and non-orthogonality, 6 values):
 * 9 parameters. Fitted by Levenberg-Marquardt (Gauss-Newton with damping) to static pose means x_i, minimising
 * sum_i (|A (x_i - b)| - g)^2, starting from b = 0, A = I. The poses only need to be different and still: no
 * orientation has to be known or aligned. The rotation part of the accelerometer matrix is not observable from
 * magnitudes and stays in the mounting matrix (vehicle side).
 *
 * With fewer than minPosesFull (9) but at least minPosesDiag (6) poses, a diagonal fit (b + diag(A), 6 parameters)
 * is used and the result is at best Warn.
 *
 * Internally the fit runs in units of g (x / g), so the condition number of the normal matrix compares offset and
 * scale columns of equal size.
 *
 * Pure code: no Arduino, no ExVectrCore. Host-tested (flight-software test/test_cal_accel_fit).
 */

#include <stddef.h>
#include <stdint.h>

#include "ExVectrDSP/calibration/calibration_result.hpp"

namespace VCTR {
namespace DSP {
namespace Calib {

/// Accept limit of the accelerometer offset |b_k| in m/s^2 (owner decision 2026-11-21: 2.5, was 1.5). The vehicle's
/// stored-parameter check (Arm/MissionBegin) uses this same constant.
constexpr double kAccelOffsetLimit = 2.5;

struct AccelFitConfig {
  double gRef = 9.81;              ///< gravity reference, m/s^2 (CAL_G_REF)
  size_t minPosesFull = 9;         ///< poses for the full 9-parameter fit
  size_t minPosesDiag = 6;         ///< poses for the diagonal fallback
  double coverageAngleRad = 25.0 * 3.14159265358979 / 180.0; ///< each sensor axis within this of +g and -g once
  double maxResidualRms = 0.03;    ///< m/s^2
  double maxResidualMax = 0.06;    ///< m/s^2
  double scaleMin = 0.95;          ///< eigenvalues of A
  double scaleMax = 1.05;
  double maxOffDiag = 0.02;        ///< |A_jk|, j != k
  double maxOffset = kAccelOffsetLimit; ///< |b_k|, m/s^2
  double warnOffset = 0.8;         ///< |b_k| above this: Warn (beyond the MPU9250 datasheet zero-g tolerance)
  double maxCondition = 1e4;       ///< condition number of the (normalised) normal matrix
  int maxIterations = 50;
};

/// Pass bits of AccelFitResult::checks (bit set = check passed).
enum AccelCheck : uint16_t {
  kAccelCheckPosesFull = 1u << 0, ///< n >= minPosesFull (full fit, not the diagonal fallback)
  kAccelCheckCoverage = 1u << 1,  ///< every sensor axis within coverageAngle of +g and of -g at least once
  kAccelCheckSolved = 1u << 2,    ///< the iteration produced finite parameters
  kAccelCheckCondition = 1u << 3, ///< condition number < maxCondition
  kAccelCheckResidualRms = 1u << 4,
  kAccelCheckResidualMax = 1u << 5,
  kAccelCheckScale = 1u << 6,     ///< eigenvalues of A in [scaleMin, scaleMax]
  kAccelCheckOffDiag = 1u << 7,
  kAccelCheckOffset = 1u << 8,    ///< |b_k| < maxOffset
  kAccelCheckOffsetNominal = 1u << 9, ///< |b_k| < warnOffset (warning only)
};

struct AccelFitResult {
  CalResult result = CalResult::None;
  uint16_t checks = 0;     ///< AccelCheck pass bits
  bool diagonal = false;   ///< diagonal fallback used
  size_t nUsed = 0;        ///< poses used
  double offset[3] = {0, 0, 0}; ///< b, raw sensor frame, m/s^2
  double A[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; ///< symmetric, row-major, sensor frame
  double scaleEig[3] = {1, 1, 1}; ///< eigenvalues of A, ascending
  double residualRms = 0;  ///< of |a_cal| - g over the poses, m/s^2
  double residualMax = 0;  ///< max | |a_cal| - g |, m/s^2
  double condition = 0;    ///< of the normalised normal matrix (inf if singular)
  int iterations = 0;
};

/**
 * Fits b and A to n pose means (raw sensor frame, m/s^2). n may be any number; the verdict covers too few poses.
 * Verdict order: too few poses / coverage -> FailCoverage; not solvable or condition -> Numeric; residuals ->
 * FailResidual; scale, off-diagonal, offset -> FailRange; else Warn (diagonal fallback or offset above warnOffset)
 * or Pass.
 */
AccelFitResult fitAccel(const double (*poses)[3], size_t n, const AccelFitConfig &cfg = AccelFitConfig());

/// |A (x - b)| with the parameters of r.
double accelCalibratedNorm(const AccelFitResult &r, const double x[3]);

/// Bit k set = direction k covered by at least one pose within coverageAngleRad of it. Directions in the sensor frame
/// as the fit sees it (raw pose means): bit 0 +X, 1 -X, 2 +Y, 3 -Y, 4 +Z, 5 -Z (the axis pointing along +g
/// reading, i.e. the axis the specific force points to). 0x3F = all covered.
uint8_t accelCoverageMask(const double (*poses)[3], size_t n, double coverageAngleRad);

/// True if every sensor axis lies within coverageAngleRad of +g and of -g for at least one pose direction
/// (accelCoverageMask == 0x3F).
bool accelPoseCoverage(const double (*poses)[3], size_t n, double coverageAngleRad);

/// Number of covered directions (0..6) in a mask.
int accelCoverageCount(uint8_t mask);

/// Writes the missing directions of a mask as "+X, -Z" (bit order) into buf (empty string if none missing).
void accelMissingText(uint8_t mask, char *buf, size_t bufLen);

/// True when the operator can stop and fit: at least minPoses poses and all six directions covered.
bool accelCoverageComplete(size_t nPoses, uint8_t mask, size_t minPoses);

} // namespace Calib
} // namespace DSP
} // namespace VCTR

#endif
