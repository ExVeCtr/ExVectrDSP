#ifndef EXVECTRDSP_SIMULATOR_HPP
#define EXVECTRDSP_SIMULATOR_HPP

#include "ExVectrControl/starship/control_attitude_flaps.hpp"
#include "ExVectrCore/task_types.hpp"
#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"
#include "ExVectrDSP/body_sim_model.hpp"
#include "ExVectrMath.hpp"

namespace VCTR {

namespace DSP {

/**
 * @brief Simulation task: replaces the estimators in simulation mode. It feeds the TVC "actual thrust" and the true
 * flap angles into the pure physics model `BodySimModel` (body_sim_model.hpp: geometry, aerodynamics, rigid body,
 * all parameters in `BodySimParams`) and publishes the same topic types as the estimators.
 *
 * Time: the model integrates in fixed steps (`BodySimParams::step_s`, 2 ms) of the wall-clock time elapsed since
 * the previous run, at most `maxAdvance_s` (50 ms) per run. A gap of more than 100 ms (task paused, simulation
 * restarted) is not integrated at all: the simulation simply continues from where it stopped.
 */
class BodySimulator : public Core::Task_Periodic {
protected:
  Core::Topic<Core::Timestamped<Math::Vector<float, 7>>> attitudeTopic_;
  Core::Topic<Core::Timestamped<Math::Vector<float, 6>>> positionTopic_;
  Core::Topic<Core::Timestamped<Math::Vector<float, 2>>> groundTopic_;

  // Where TVC input is published. In form: [X, Y, Z, T], where X, Y, Z show the
  // thrust vector in body frame (magnitude of vector is thrust magnitude) and T
  // is the roll torque about body Z in N m.
  Core::Simple_Subscriber<Math::Vector<float, 4>> tvcSubr_;
  Core::Simple_Subscriber<CTRL::ControlAttitudeFlapSetting> flapSettingsSubr_;

  BodySimModel model_;

  Math::Vector<float, 3> tvcForce_;  ///< last thrust vector received, body frame, N
  float tvcRollTorque_ = 0;          ///< last roll torque received, N m

  /// Timestamp of the published state (time of the last run).
  int64_t stateTimestamp_ = 0;

  bool zeroingMode_ = false; ///< While true the body is held at the start pose every run.
  bool tvcEnabled_ = true;   ///< Thrust and roll torque act only while true.
  bool flapsEnabled_ = false; ///< Informational: the flap driver publishes true angles only while enabled.

  Math::Vector<float, 7> attitudeState_; ///< [omega(3) body rad/s, q(4) body -> world]
  Math::Vector<float, 6> positionState_; ///< [V(3), P(3)] world NWU, m/s and m

public:
  /// Gap above which the elapsed wall-clock time is not integrated (pause, restart).
  static constexpr int64_t kMaxGapNs = 100 * Core::MILLISECONDS;

  /**
   * @param mass kg, @param inertiaTensor diagonal of the inertia tensor about the CoM, kg m^2,
   * @param tvcThrustLimit_N, tvcAngleLimit_Rad thrust magnitude and cone limits.
   * All other model parameters are the defaults of BodySimParams (see params()).
   */
  BodySimulator(float mass, Math::Vector<float, 3> inertiaTensor,
                float tvcThrustLimit_N, float tvcAngleLimit_Rad,
                int64_t simulationInterval = 0.01 * Core::SECONDS);

  /**
   * @brief Subscribes to the TVC actual thrust topic [X, Y, Z, T] (body N, roll torque N m).
   * @param tvcPosition Thrust application point in the body frame; only its z (on the axis) is used.
   */
  void setTVCInputTopic(Core::Topic<Math::Vector<float, 4>> &tvcTopic,
                        const Math::Vector<float, 3> &tvcPosition,
                        float tvcAngleLimit, float tvcForceLimit);

  void setFlapSettingInputTopic(
      Core::Topic<CTRL::ControlAttitudeFlapSetting> &flapTopic) {
    flapSettingsSubr_.subscribe(flapTopic);
  }

  /// One simulation run: take the inputs, integrate the elapsed time, apply zeroing mode.
  void update();

  /**
   * @brief While enabled the body is held at the start pose every tick (resetToStartPose()).
   */
  void enableZeroingMode(bool enable) { zeroingMode_ = enable; }

  /**
   * @brief Puts the simulated body back to the start pose: upright (q = identity), at the reference origin, all
   * velocities zero. The same pose is held by zeroing mode, so a reset and the end of a zeroing window leave the
   * body in the same state.
   */
  void resetToStartPose();

  void enableTVC(bool enable) { tvcEnabled_ = enable; }
  void enableFlaps(bool enable) { flapsEnabled_ = enable; }

  /**
   * @brief Steady wind of the simulated world (simulation only; off by default).
   * @param speed_mps horizontal speed, m/s (clamped to >= 0); @param direction_rad direction the air moves
   * TOWARD, counterclockwise from North toward West (NWU). Non-finite values leave the setting unchanged.
   */
  void setWind(float speed_mps, float direction_rad) { model_.setWind(speed_mps, direction_rad); }
  float windSpeed() const { return model_.params().windSpeed_mps; }
  float windDirection() const { return model_.params().windDirection_rad; }

  /// Model parameters (change only while the simulation is paused; call model().applyParams() after mass,
  /// inertia or geometry changes).
  BodySimModel &model() { return model_; }

  const VCTR::Math::Vector<float, 7> &getAttitudeState() { return attitudeState_; }
  const VCTR::Math::Vector<float, 6> &getPositionState() { return positionState_; }

  /// @return attitude [W, Q] (body rates, q body -> world).
  Core::Topic<Core::Timestamped<Math::Vector<float, 7>>> &getAttitudeEstTopic() { return attitudeTopic_; }
  /// @return [V, P] in the world frame (NWU).
  Core::Topic<Core::Timestamped<Math::Vector<float, 6>>> &getStateEstTopic() { return positionTopic_; }
  /// @return ground [V, P]: vertical velocity m/s, height above the ground m.
  Core::Topic<Core::Timestamped<Math::Vector<float, 2>>> &getGroundTopic() { return groundTopic_; }

  void setAttitudeState(const Math::Vector<float, 7> &attitudeState) { model_.setAttitudeState(attitudeState); }
  void setPositionState(const Math::Vector<float, 6> &positionState) { model_.setPositionState(positionState); }

private:
  void taskInit() override;
  void taskThread() override;
};

} // namespace DSP

} // namespace VCTR

#endif
