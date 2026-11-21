#include "ExVectrDSP/calibration/sphere_coverage.hpp"

#include <math.h>

namespace VCTR {
namespace DSP {
namespace Calib {

namespace {
// Bin centres, computed once on first use (64 x 3 floats).
float gCentres[SphereCoverage::kBins][3];
bool gCentresReady = false;

void ensureCentres() {
  if (gCentresReady)
    return;
  for (size_t i = 0; i < SphereCoverage::kBins; i++) {
    double c[3];
    SphereCoverage::binCentre(i, c);
    gCentres[i][0] = float(c[0]);
    gCentres[i][1] = float(c[1]);
    gCentres[i][2] = float(c[2]);
  }
  gCentresReady = true;
}
} // namespace

void SphereCoverage::binCentre(size_t i, double c[3]) {
  const double n = double(kBins);
  const double z = 1.0 - (2.0 * double(i) + 1.0) / n;
  const double r = sqrt(1.0 - z * z);
  const double phi = double(i) * 3.14159265358979323846 * (3.0 - sqrt(5.0));
  c[0] = r * cos(phi);
  c[1] = r * sin(phi);
  c[2] = z;
}

int SphereCoverage::binOf(const double v[3]) {
  if (!isfinite(v[0]) || !isfinite(v[1]) || !isfinite(v[2]))
    return -1;
  if (v[0] == 0 && v[1] == 0 && v[2] == 0)
    return -1;
  ensureCentres();
  int best = 0;
  double bestDot = -INFINITY;
  for (size_t i = 0; i < kBins; i++) {
    const double d = v[0] * gCentres[i][0] + v[1] * gCentres[i][1] + v[2] * gCentres[i][2];
    if (d > bestDot) {
      bestDot = d;
      best = int(i);
    }
  }
  return best;
}

void SphereCoverage::reset() {
  for (size_t i = 0; i < kBins; i++)
    counts_[i] = 0;
}

int SphereCoverage::add(const double v[3]) {
  const int b = binOf(v);
  if (b >= 0 && counts_[b] < 255)
    counts_[b]++;
  return b;
}

double SphereCoverage::fractionAtLeast(uint8_t k) const {
  size_t n = 0;
  for (size_t i = 0; i < kBins; i++)
    if (counts_[i] >= k)
      n++;
  return double(n) / double(kBins);
}

} // namespace Calib
} // namespace DSP
} // namespace VCTR
