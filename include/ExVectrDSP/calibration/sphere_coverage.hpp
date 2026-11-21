#ifndef EXVECTRDSP_SPHERE_COVERAGE_HPP
#define EXVECTRDSP_SPHERE_COVERAGE_HPP

/**
 * Direction coverage on the unit sphere: 64 equal-area bins on a Fibonacci lattice
 * (workspace docs/calibration-design.md section 3.5, shared definition with the ground station):
 *
 *   bin centre i (0..63): z = 1 - (2i+1)/64, r = sqrt(1 - z^2), phi = i * pi * (3 - sqrt(5)),
 *   c_i = (r cos phi, r sin phi, z); a direction u belongs to the bin with the largest u . c_i.
 *
 * Pure code: no Arduino, no ExVectrCore. Host-tested (flight-software test/test_cal_coverage).
 */

#include <stddef.h>
#include <stdint.h>

namespace VCTR {
namespace DSP {
namespace Calib {

class SphereCoverage {
public:
  static constexpr size_t kBins = 64;

  /// Centre of bin i (unit vector).
  static void binCentre(size_t i, double c[3]);

  /// Bin of direction v (need not be normalised). Returns -1 for a zero or non-finite vector.
  static int binOf(const double v[3]);

  void reset();

  /// Adds one direction; counts saturate at 255. Returns the bin, or -1 if v was rejected.
  int add(const double v[3]);

  uint8_t count(size_t i) const { return i < kBins ? counts_[i] : 0; }

  /// Fraction (0..1) of bins with at least k samples.
  double fractionAtLeast(uint8_t k) const;

private:
  uint8_t counts_[kBins] = {0};
};

} // namespace Calib
} // namespace DSP
} // namespace VCTR

#endif
