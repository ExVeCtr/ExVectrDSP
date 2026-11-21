#ifndef EXVECTRDSP_BODY_SIM_MODEL_HPP
#define EXVECTRDSP_BODY_SIM_MODEL_HPP

/**
 * Pure physics model of the Starship body simulator (feature/sim-aero, 2026-11-21).
 *
 * No scheduler, no topics, no wall clock, no random numbers: the caller drives it with a fixed step, so identical
 * inputs give bit-identical runs (host-tested in the flight-software repo, test/test_sim_aero). `BodySimulator`
 * (body_simulator.hpp) is the scheduler task that wraps it on the vehicle.
 *
 * Frames and units: body frame +Z nose, belly -X, flaps TL/BL on +Y and TR/BR on -Y, origin at the centre of
 * mass; world NWU; q = body -> world [w,x,y,z]; attitude state [omega(3) body rad/s, q(4)], position state
 * [V(3), P(3)] world m/s, m. SI throughout. Flap angle: 0 = fully extended along +-Y, 90 deg = folded toward +X.
 *
 * Aerodynamic model (docs/libs/ExVectrDSP.md, "Body simulator model"):
 * - Body crossflow: strip theory over the real length (tube + parabolic nose), each strip
 *   dF = 1/2 rho eta Cdc d(z) dz |u_perp| u_perp with the local air velocity u = u_cg - omega x r, so the
 *   distributed pitch/yaw damping follows from the geometry.
 * - Body potential (slender-body, Jorgensen) normal force 1/2 rho V^2 A_base sin(2 alpha) cos(alpha/2) at the
 *   leading end (nose centre of pressure when the nose leads, the base face when the tail leads): this gives the
 *   angle-of-attack dependent centre of pressure and the belly-down pitch stiffness.
 * - Body axial force 1/2 rho Cd_axial A_base u_z |u_z| (separate coefficients for tail-first and nose-first flow).
 * - Each flap: flat plate normal force 1/2 rho |u|^2 A C_N(alpha) along the plate normal, C_N blended from the
 *   attached (C_N_alpha sin a cos a) to the stalled (C_N90 sin a) branch, centre of pressure moving along the
 *   chord from 1/4 (small alpha) to 1/2 (alpha = 90 deg): x_cp/c = 0.5 -/+ 0.25 cos(alpha). The flap geometry in
 *   the aerodynamics uses the effective angle a_eff = smoothmax(a, 20 deg) (an extended flap still sees a small
 *   effective angle, BodySimParams::flapEffectiveMinAngleRad); the servo angle stays the true angle.
 * - One global calibration factor scales every aerodynamic force (see BodySimParams::aeroCalibration).
 */

#include <stdint.h>

#include "ExVectrControl/physics_bodysim.hpp"
#include "ExVectrMath.hpp"

namespace VCTR {
namespace DSP {

/// One flap pair (left and right flap identical, mirrored about the body X-Z plane).
struct BodySimFlapGeometry {
  float width_m;              ///< Radial extent from the hinge line to the outer edge (the chord for belly-down flow), m.
  float height_m;             ///< Extent along the body axis, m.
  float chamfer_m;            ///< Length of the 45 deg cut at the outer nose-side corner, m (0 = rectangle).
  float bottomFromShipBottom_m; ///< Lower (tail-side) edge above the bottom of the ship, m.
  float effectiveness;        ///< Multiplies this pair's aerodynamic force (body shielding, gaps). 1 = free plate.
};

/**
 * All parameters of the model with units and physical meaning. Defaults: owner's ship facts (2026-11-21) and
 * textbook coefficients, with `aeroCalibration` chosen so that the belly-down steady fall at neutral flaps
 * (top 50 deg, bottom 40 deg) is 6.4 m/s (real flight, 21..25 km/h, owner 2026-11-21).
 */
struct BodySimParams {
  // ---- Mass properties (measured, double swing; owner 2026-11-21) ----
  float mass_kg = 1.1f;   ///< Total mass incl. battery and motors (= VEHICLE_MASS_KG).
  float inertiaXX = 0.044f; ///< kg m^2 about the CoM, body X (pitch/yaw axes are X and Y).
  float inertiaYY = 0.044f; ///< kg m^2, body Y.
  float inertiaZZ = 0.004f; ///< kg m^2, body Z (long axis).
  float inertiaXY = 0.0f, inertiaXZ = 0.0f, inertiaYZ = 0.0f; ///< Products of inertia, kg m^2 (tensor entries).

  // ---- Geometry (owner 2026-11-21). Origin of the body frame = CoM on the axis. ----
  float bodyRadius_m = 0.065f;     ///< Tube radius (diameter 13 cm).
  float tubeLength_m = 0.54f;      ///< Cylindrical tube, bottom face to the nose base.
  float noseLength_m = 0.20f;      ///< Parabolic nose (r = R (2 s - s^2), s = distance from the tip / length).
  float cgFromBottom_m = 0.37f;    ///< CoM above the bottom of the ship ("basically the middle" of 0.74 m).
  float flapHingeRadius_m = 0.065f; ///< Hinge line distance from the axis (on the tube surface).
  BodySimFlapGeometry flapTop = {0.07f, 0.14f, 0.08f, 0.53f, 1.0f};    ///< TL, TR.
  BodySimFlapGeometry flapBottom = {0.08f, 0.21f, 0.08f, 0.01f, 1.0f}; ///< BL, BR: spans 0.01..0.22 m above the bottom.
  /**
   * Inward cant of the TOP flaps (owner 2026-11-21, 10 deg): each top flap is mounted tilted about the
   * horizontal line through its LOWER hinge point so that its upper end moves toward the body axis (the top
   * connection sits close to the nose cone). The lower hinge point stays at (r = flapHingeRadius, z = lower edge);
   * the hinge line, plate normal, centroid, centre of pressure and area direction are rotated with it; the flap
   * angle convention is unchanged (rotation about the tilted hinge line). Bottom flaps are not canted.
   */
  float topFlapCantRad = 0.17453f; ///< 10 deg.
  /**
   * Fold direction of the flaps. false (default, CONFIRMED by the owner 2026-11-21: folded flaps point toward +X,
   * the back): an increasing angle folds the flap toward the back (+X, downstream in the belly-down fall); true:
   * toward the belly (-X, upstream), kept for comparison.
   * Only the sign of the flaps' sideways (Y) force changes, so only the azimuth (tau_x) response to the
   * controller's u0 mixer flips: with false it is reversed against ControlAttitudeFlaps at every angle of attack
   * tested, with true it has the controller's sign; pitch, roll, drag and terminal speed are identical
   * (test_sim_aero "control axis signs", ExVectrDSP.md Known issues).
   */
  bool flapFoldTowardBelly = false;

  // ---- Air ----
  float airDensity_kgm3 = 1.2f; ///< rho.

  // ---- Body aerodynamics ----
  float cdCrossflow = 1.2f;     ///< 2-D circular cylinder crossflow drag coefficient (subcritical Re ~ 7e4).
  float crossflowEta = 0.65f;   ///< Finite-length factor on the crossflow drag (L/D ~ 5.7).
  float potentialFactor = 1.0f; ///< Multiplies the slender-body potential normal force (0 disables it).
  float cdAxialTailFirst = 0.9f; ///< Axial drag coefficient on the base area, flow from the tail (falling upright).
  float cdAxialNoseFirst = 0.35f; ///< Axial drag coefficient on the base area, flow from the nose.
  int bodyStrips = 24;          ///< Number of strips over the body length (crossflow and damping), 4..64.

  // ---- Flap aerodynamics (flat plate) ----
  float flapCnAlpha = 3.0f;           ///< Attached normal-force slope, 1/rad (Helmbold, plate AR ~ 2.5).
  float flapCn90 = 1.17f;             ///< Normal-force coefficient broadside (alpha = 90 deg).
  float flapStallAngle_rad = 0.5236f; ///< 30 deg: centre of the attached -> stalled blend (low-AR plates stall late).
  float flapStallWidth_rad = 0.0873f; ///< 5 deg: width of the blend (logistic).

  /**
   * Effective aerodynamic flap angle (owner 2026-11-21). The flow around the body means a fully extended flap
   * (fold angle 0) still sees a small effective angle instead of producing pure drag. The aerodynamic calculation
   * (plate normal = force direction, C_N argument, centre of pressure) uses
   *   a_eff = smoothmax(a, a_min) = 0.5 (a + a_min + sqrt((a - a_min)^2 + k^2)),  k = flapEffectiveSmoothRad,
   * which is C1 (C-infinity) and converges to a_min (+ k^2 / (4 a_min)) for a -> 0 and to a far above a_min. The
   * servo angle (flapAngle(), telemetry) stays the true angle. flapEffectiveMinAngleRad <= 0 disables it.
   */
  float flapEffectiveMinAngleRad = 0.349066f; ///< a_min, 20 deg.
  float flapEffectiveSmoothRad = 0.087266f;   ///< k, 5 deg: transition width of the smooth max / min.
  /// Same treatment at the retracted end: a_eff = smoothmin(a_eff, 90 deg - margin). 0 = off (a folded flap is
  /// edge-on to the belly-down flow; kept as an option).
  float flapEffectiveRetractMarginRad = 0.0f;

  // ---- Candidate flap/body interaction effects (azimuth sign investigation, 2026-11-21; all off by default) ----
  /**
   * (a) Sidewash at the flap: the local crossflow is turned toward the body axis (gain > 0) or away from it
   * (gain < 0) by delta = gain (R / r_c)^2 sin(alpha_body), r_c = distance of the flap's reference point from the
   * axis, alpha_body = angle between the flow and the body axis (90 deg belly-down). Implemented as an extra
   * velocity delta |u_perp| along -r_hat (radial unit vector of the flap point).
   */
  float flapSidewashGain = 0.0f;
  /**
   * (b) Body-crossflow induced flow: adds gain x the 2-D potential-flow perturbation of a circular cylinder of
   * radius R in the body crossflow, evaluated at the flap's reference point: (u - i v)_p = -V R^2 / z^2 with
   * V = u_x + i u_y (local crossflow) and z = x + i y. 1 = full potential flow, 0 = off.
   */
  float flapBodyCrossflowGain = 0.0f;
  /**
   * (c) Flap/body interference: the effective flap area is multiplied by 1 - factor sin^2(a_eff), i.e. a flap
   * folding against the body side (into its boundary layer / wake) loses area. 0 = off, 1 = no force when folded.
   * ON by default since 2026-11-21 (owner-approved, 0.7): it is the only candidate effect that gives the flap
   * controller's expected yaw (azimuth) sign at the neutral angles 50/40 deg over AoA -20..+20 deg, and the real
   * ship's heading loop works. Fitted to the sign, not measured: the owner plans a bench pivot test in an airflow.
   */
  float flapFoldShieldFactor = 0.7f;

  /**
   * Global calibration factor on every aerodynamic force and moment. Without it (1.0) the model falls belly-down
   * at neutral flaps at about 13.5 m/s; the real ship falls at 6.4 m/s (21..25 km/h, flight position filter), which
   * needs a factor of 4.64 with the fold shield 0.7 (4.00 without it; an effective drag area of 0.43 m^2 instead of
   * 0.097 m^2). That is NOT a plausible
   * coefficient for the documented 1.1 kg ship of 13 x 74 cm: it stands for an unexplained difference (mass, the
   * speed measurement, or missing drag sources) to be resolved with the owner.
   * Because it scales all aerodynamic forces equally, it also scales stiffness, damping and flap authority.
   */
  float aeroCalibration = 4.64f;

  // ---- Actuators ----
  /// True: the model applies the servo rate limits itself (host harness). False on the vehicle, where the flap
  /// driver (ServoAngleSimulator) and the TVC driver already publish rate-limited "actual" values.
  bool applyActuatorRateLimits = true;
  float flapServoRate_radps = 5.23599f;   ///< 300 deg/s (9 g servo, 60 deg per 0.2 s, owner's choice).
  float flapServoLag_s = 0.0f;            ///< Optional first-order lag of the flap angle, s (0 = none).
  float flapServoDeadTime_s = 0.0f;       ///< Optional dead time of the flap angle, s (0 = none, max 0.25 s).
  float flapMinAngle_rad = 0.0f;          ///< Mechanical limits of the modelled flap angle.
  float flapMaxAngle_rad = 1.5708f;       ///< 90 deg.
  float tvcServoRate_radps = 5.23599f;    ///< 300 deg/s vane servo rate (same 9 g servos).
  float tvcServoToThrustGain = 15.0f;     ///< Servo angle per thrust angle (StarshipTVC tvcFactor_ default 15).
  float tvcThrustLimit_N = 14.0f;         ///< |F| limit (TVC_THRUST_LIMIT_N).
  float tvcAngleLimit_rad = 0.05236f;     ///< Thrust angle limit from body +Z (45 deg / 15 = 3 deg).
  float tvcFromCg_m = -0.35f;             ///< Thrust application point on the axis, z (body frame, m).

  // ---- Wind (steady, no gusts; off by default) ----
  float windSpeed_mps = 0.0f;     ///< Horizontal wind speed, m/s (>= 0).
  float windDirection_rad = 0.0f; ///< Direction the air moves TOWARD, counterclockwise from North toward West (NWU).

  // ---- Gusts (simulation only, owner-approved exception 2026-11-21; deterministic, off by default) ----
  /**
   * Gust model: the gust is a deterministic sum of 8 sinusoids with equal amplitude strength/8, frequencies
   * (0.5 .. 2) / gustPeriod_s and phases drawn from a xorshift32 generator seeded with gustSeed (no wall clock, no
   * global random state), evaluated at the simulated time. It is added to the mean wind in the world frame, along
   * the mean wind direction (a head or tail gust); a second independent process of the same kind, scaled by
   * gustSpread, acts horizontally across the mean direction. |gust| along the mean never exceeds gustStrength_mps
   * (hard bound; realised peaks are about 0.6..0.8 of it), RMS = strength / 4, correlation time about the period.
   * Same seed and parameters give a bit-identical gust series.
   */
  float gustStrength_mps = 0.0f; ///< Peak (bound) of the extra speed along the mean direction, m/s (0..10; 0 = no gusts).
  float gustPeriod_s = 5.0f;     ///< Gust timescale, s (0.5..30): the centre period of the sinusoids.
  uint32_t gustSeed = 1;         ///< Seed of the phases and frequency ratios.
  float gustSpread = 0.3f;       ///< Amplitude of the crosswind gust process relative to the along-wind one (0..1).

  // ---- Integration ----
  float step_s = 0.002f;     ///< Fixed integration step, s.
  float maxAdvance_s = 0.05f; ///< advance() integrates at most this much per call (bounded catch-up), s.
};

/// Output of one force evaluation (body frame), kept for logging and tests.
struct BodySimLoads {
  Math::Vector<float, 3> aeroForce;   ///< Sum of aerodynamic forces, N.
  Math::Vector<float, 3> aeroTorque;  ///< Sum of aerodynamic moments about the CoM, N m.
  Math::Vector<float, 3> flapForce[4]; ///< Per flap TL, TR, BL, BR, N.
  Math::Vector<float, 3> airVelocity; ///< Air velocity relative to the CoM, body frame, m/s.
};

class BodySimModel {
 public:
  explicit BodySimModel(const BodySimParams &params = BodySimParams());

  /// Parameters (change, then call applyParams() for mass, inertia or geometry changes).
  BodySimParams &params() { return params_; }
  const BodySimParams &params() const { return params_; }
  void applyParams();

  /// Upright (q = identity) at the origin, at rest; actuators keep their state.
  void resetToStartPose();
  void setAttitudeState(const Math::Vector<float, 7> &s) { body_.setAttitudeState(s); }
  void setPositionState(const Math::Vector<float, 6> &s) { body_.setPositionState(s); }

  /**
   * Thrust input: force vector in the body frame (its magnitude is the thrust, N) and the roll torque about body
   * Z (N m). Both act only while `enabled`. The angle and magnitude limits are applied here.
   */
  void setThrust(const Math::Vector<float, 3> &forceBody, float rollTorque, bool enabled);
  /// Flap angle inputs TL, TR, BL, BR, rad (commanded angles in the harness, "true" angles on the vehicle).
  void setFlapInput(float tl, float tr, float bl, float br);
  void setWind(float speed_mps, float direction_rad);
  /// Gust settings (clamped: strength 0..10 m/s, period 0.5..30 s); non-finite values are ignored. Rebuilds the
  /// gust table; the gust series depends only on these values and the simulated time.
  void setGust(float strength_mps, float period_s, uint32_t seed);
  /// Wind and gusts off, defaults for the period and seed.
  void clearWind();
  /// Gust speed along (index 0) and across (index 1) the mean wind direction at the current simulated time, m/s.
  void gustAt(double t_s, float &along, float &across) const;

  /// One fixed step of params().step_s.
  void step();
  /**
   * Integrates `elapsed_s` of real time in fixed steps (the remainder is carried to the next call). Elapsed time
   * above maxAdvance_s is dropped (bounded catch-up, no huge step after a pause). Returns the number of steps.
   */
  int advance(float elapsed_s);
  /// Drops any carried remainder (call when the simulation is (re)started after a pause).
  void clearTimeCarry() { carry_s_ = 0.0; }

  const Math::Vector<float, 7> &attitudeState() { return body_.getAttitudeState(); }
  const Math::Vector<float, 6> &positionState() { return body_.getPositionState(); }
  /// Modelled flap angles TL, TR, BL, BR after the servo model, rad.
  float flapAngle(int i) const { return flapAngle_[i]; }
  /// Modelled thrust angles [about body X, about body Y] after the TVC slew model, rad.
  float tvcAngleX() const { return tvcAngleX_; }
  float tvcAngleY() const { return tvcAngleY_; }
  const BodySimLoads &loads() const { return loads_; }
  /// Simulated time since construction, s (fixed steps x step_s).
  double time() const { return time_s_; }
  /// World wind vector (NWU, air velocity), m/s.
  Math::Vector<float, 3> windWorld() const;

  /// Flat-plate normal-force coefficient at angle of attack alpha (0..pi/2), exposed for tests.
  float flapCn(float alpha) const;
  /// Effective aerodynamic flap angle for the true (servo) angle `a`, rad (see flapEffectiveMinAngleRad).
  float flapAeroAngle(float a) const;
  /// Plate geometry of flap `index` (TL, TR, BL, BR) at flap angle `angle` in the body frame: area centroid and unit
  /// normal (cant included), exposed for tests.
  void flapPlate(int index, float angle, Math::Vector<float, 3> &centroid, Math::Vector<float, 3> &normal) const;
  /// Computes the aerodynamic loads for the current state and flap angles (does not integrate).
  void computeAero();

 private:
  BodySimParams params_;
  CTRL::PhysicsBodySim body_;

  // Thrust
  Math::Vector<float, 3> thrustCmd_;
  float thrustMag_ = 0.0f;
  float rollTorque_ = 0.0f;
  bool thrustEnabled_ = false;
  float tvcAngleX_ = 0.0f, tvcAngleY_ = 0.0f;

  // Flaps
  float flapInput_[4] = {1.5708f, 1.5708f, 1.5708f, 1.5708f};
  float flapRateLimited_[4] = {1.5708f, 1.5708f, 1.5708f, 1.5708f};
  float flapAngle_[4] = {1.5708f, 1.5708f, 1.5708f, 1.5708f};
  static constexpr int kDelayLen = 128;
  float flapDelay_[4][kDelayLen] = {};
  int delayIdx_ = 0;
  bool delayPrimed_ = false;

  // Derived geometry (applyParams)
  float zBottom_ = 0.0f, zNoseBase_ = 0.0f, zTip_ = 0.0f, zNoseCp_ = 0.0f;
  float baseArea_ = 0.0f;
  float flapArea_[2] = {0, 0};
  float flapSc_[2] = {0, 0};   ///< Area centroid distance from the lower edge along the flap (hinge line), m.
  float flapZLow_[2] = {0, 0}; ///< Lower hinge point z (body frame), m.
  float flapCant_[2] = {0, 0}; ///< Inward cant of the pair, rad.
  /// Orthonormal plate frame of flap `index` at `angle`: normal n, radial d (hinge -> outer edge), axial e (along
  /// the hinge line, lower -> upper end). Cant about the horizontal X line through the lower hinge point.
  void flapFrame(int index, float angle, Math::Vector<float, 3> &n, Math::Vector<float, 3> &d,
                 Math::Vector<float, 3> &e) const;

  static constexpr int kGustTerms = 8;
  float gustOmega_[2][kGustTerms] = {};
  float gustPhase_[2][kGustTerms] = {};
  void buildGust();

  BodySimLoads loads_;
  double carry_s_ = 0.0;
  double time_s_ = 0.0;

  void updateActuators(float dt);
  float bodyDiameterAt(float z) const;
  void addFlap(int index, float angle, const Math::Vector<float, 3> &u0, const Math::Vector<float, 3> &omega);
  void groundContact(const Math::Vector<float, 6> &posBefore, const Math::Vector<float, 7> &attBefore);
};

} // namespace DSP
} // namespace VCTR

#endif
