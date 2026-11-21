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
 *   chord from 1/4 (small alpha) to 1/2 (alpha = 90 deg): x_cp/c = 0.5 -/+ 0.25 cos(alpha).
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
 * (top 50 deg, bottom 40 deg) is 7.7 m/s (flight 4).
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
  BodySimFlapGeometry flapBottom = {0.08f, 0.21f, 0.08f, 0.10f, 1.0f}; ///< BL, BR.

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
   * Global calibration factor on every aerodynamic force and moment. Without it (1.0) the model falls belly-down
   * at neutral flaps at 13.5 m/s; flight 4 fell at 7.7 m/s, which needs 3.08 (a drag area of 0.30 m^2 instead of
   * 0.097 m^2). That is NOT a plausible coefficient for the documented 1.1 kg ship of 13 x 74 cm: it stands for an
   * unexplained difference (mass, the 7.7 m/s figure, or missing drag sources) to be resolved with the owner.
   * Because it scales all aerodynamic forces equally, it also scales stiffness, damping and flap authority.
   */
  float aeroCalibration = 3.08f;

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
  float flapZc_[2] = {0, 0};

  BodySimLoads loads_;
  double carry_s_ = 0.0;
  double time_s_ = 0.0;

  void updateActuators(float dt);
  float bodyDiameterAt(float z) const;
  void addFlap(int index, float angle, const Math::Vector<float, 3> &u0, const Math::Vector<float, 3> &omega);
  void groundContact();
};

} // namespace DSP
} // namespace VCTR

#endif
