#ifndef EXVECTRDSP_CALIBRATION_RESULT_HPP
#define EXVECTRDSP_CALIBRATION_RESULT_HPP

/**
 * Common result codes of the sensor calibration fits (ExVectrDSP calibration module).
 *
 * The values follow the planned `CalibResult` packet enum (workspace docs/calibration-design.md section 3.5), so a
 * later telemetry packet can carry them unchanged. Pure code: no Arduino, no ExVectrCore.
 */

#include <stdint.h>

namespace VCTR {
namespace DSP {
namespace Calib {

enum class CalResult : uint8_t {
  None = 0,      ///< no result yet
  Pass,          ///< all checks passed
  Warn,          ///< usable, but a warning check failed (for example the accel diagonal fallback)
  FailResidual,  ///< fit residuals (or sample noise) above the limit
  FailCoverage,  ///< too few poses/samples or directions not covered
  FailRange,     ///< parameters outside the plausible range
  Motion,        ///< the sensor was not still long enough
  Timeout,       ///< the session ran out of time
  Numeric,       ///< fit not solvable (degenerate input, singular or ill-conditioned)
};

/// True for Pass and Warn: a candidate with this result may be applied (owner decision D6).
inline bool acceptable(CalResult r) { return r == CalResult::Pass || r == CalResult::Warn; }

/// Short upper-case text for logs, for example "PASS".
const char *toString(CalResult r);

} // namespace Calib
} // namespace DSP
} // namespace VCTR

#endif
