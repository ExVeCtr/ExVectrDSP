#include "ExVectrDSP/body_simulator.hpp"

#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrMath.hpp"

namespace VCTR {

namespace DSP {

namespace {
BodySimParams vehicleParams(float mass, const Math::Vector<float, 3> &inertia, float thrustLimit, float angleLimit) {
  BodySimParams p;
  p.mass_kg = mass;
  p.inertiaXX = inertia(0);
  p.inertiaYY = inertia(1);
  p.inertiaZZ = inertia(2);
  p.tvcThrustLimit_N = thrustLimit;
  p.tvcAngleLimit_rad = angleLimit;
  // On the vehicle the flap driver and the TVC driver already publish rate-limited "actual" values.
  p.applyActuatorRateLimits = false;
  return p;
}
} // namespace

BodySimulator::BodySimulator(float mass, Math::Vector<float, 3> inertiaTensor,
                             float tvcThrustLimit_N, float tvcAngleLimit_Rad,
                             int64_t simulationInterval)
    : Core::Task_Periodic("Simlation task", simulationInterval),
      model_(vehicleParams(mass, inertiaTensor, tvcThrustLimit_N, tvcAngleLimit_Rad)) {
  Core::getSystemScheduler().addTask(*this);
  tvcForce_ = {0, 0, 0}; // no thrust until the first TVC message (was +m g before: a phantom hover force)
  attitudeState_ = {0, 0, 0, 1, 0, 0, 0};
  positionState_ = {0, 0, 0, 0, 0, 0};
}

void BodySimulator::setTVCInputTopic(Core::Topic<Math::Vector<float, 4>> &tvcTopic,
                                     const Math::Vector<float, 3> &tvcPosition, float tvcAngleLimit,
                                     float tvcForceLimit) {
  model_.params().tvcAngleLimit_rad = tvcAngleLimit;
  model_.params().tvcThrustLimit_N = tvcForceLimit;
  model_.params().tvcFromCg_m = tvcPosition(2);
  tvcSubr_.subscribe(tvcTopic);
}

void BodySimulator::resetToStartPose() {
  model_.resetToStartPose();
  model_.clearTimeCarry();
  attitudeState_ = model_.attitudeState();
  positionState_ = model_.positionState();
}

void BodySimulator::update() {
  const int64_t now = Core::NowNs();
  const int64_t elapsedNs = now - stateTimestamp_;
  stateTimestamp_ = now;

  if (tvcSubr_.isDataNew()) {
    const auto in = tvcSubr_.getItem();
    tvcForce_ = {in(0), in(1), in(2)};
    tvcRollTorque_ = in(3);
  }
  model_.setThrust(tvcForce_, tvcRollTorque_, tvcEnabled_);

  if (flapSettingsSubr_.isDataNew()) {
    const auto f = flapSettingsSubr_.getItem();
    model_.setFlapInput(f.flapTLAngle_Rad, f.flapTRAngle_Rad, f.flapBLAngle_Rad, f.flapBRAngle_Rad);
  }

  if (elapsedNs > 0 && elapsedNs <= kMaxGapNs) {
    model_.advance(float(double(elapsedNs) / double(Core::SECONDS)));
  } else {
    model_.clearTimeCarry(); // first run after a pause or restart: continue without a catch-up step
  }

  if (zeroingMode_) {
    model_.resetToStartPose();
  }

  attitudeState_ = model_.attitudeState();
  positionState_ = model_.positionState();
}

void BodySimulator::taskInit() { stateTimestamp_ = Core::NowNs(); }

void BodySimulator::taskThread() {
  update();

  Core::Timestamped<Math::Vector<float, 7>> attState;
  attState.timestamp = stateTimestamp_;
  attState.data = attitudeState_;
  Core::Timestamped<Math::Vector<float, 6>> posState;
  posState.timestamp = stateTimestamp_;
  posState.data = positionState_;
  Core::Timestamped<Math::Vector<float, 2>> groundState;
  groundState.timestamp = stateTimestamp_;
  groundState.data = {positionState_(2), positionState_(5)}; // [vertical velocity, height above the ground]
  attitudeTopic_.publish(attState);
  positionTopic_.publish(posState);
  groundTopic_.publish(groundState);
}

} // namespace DSP

} // namespace VCTR
