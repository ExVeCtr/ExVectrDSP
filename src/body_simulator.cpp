#include "ExVectrDSP/body_simulator.hpp"

#include "ExVectrCore/print.hpp"
#include "ExVectrCore/random.h"
#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrDSP/value_covariance.hpp"
#include "ExVectrMath.hpp"

namespace VCTR {

namespace DSP {

BodySimulator::BodySimulator(float mass, Math::Vector<float, 3> inertiaTensor,
                             float tvcThrustLimit_N, float tvcAngleLimit_Rad,
                             int64_t simulationInterval)
    : Core::Task_Periodic("Simlation task", simulationInterval),
      bodySim_(mass, inertiaTensor) {
  Core::getSystemScheduler().addTask(*this);

  mass_ = mass;
  inertiaTensor_ = Math::Matrix<float, 3, 3>::eye();
  inertiaTensor_(0, 0) = inertiaTensor(0);
  inertiaTensor_(1, 1) = inertiaTensor(1);
  inertiaTensor_(2, 2) = inertiaTensor(2);

  bodySim_.setInertiaTensor(inertiaTensor_);

  attitudeState_ = {0, 0, 0, 1, 0, 0, 0};
  positionState_ = {0, 0, 0, 0, 0, 0};

  totalBodyForce = Math::GRAVITY_3F * mass_;
  totalBodyTorque = 0;

  tvcPosition_ = {0, 0, -0.35};

  tvcAngleLimit_ = tvcAngleLimit_Rad;
  tvcForceLimit_ = tvcThrustLimit_N;
}

void BodySimulator::update() {
  using namespace VCTR::Math;

  float dTime = float(Core::NOW() - stateTimestamp_) / Core::SECONDS;
  stateTimestamp_ = Core::NOW();

  // gather force and torque from the thrust vector control input
  if (tvcSubr_.isDataNew()) {
    auto tvcInput = tvcSubr_.getItem();
    Math::Vector<float, 3> tvcForce = {tvcInput(0), tvcInput(1), tvcInput(2)};

    // output TVC force to LOG_MSG
    LOG_MSG("TVC Force: %.2f %.2f %.2f\n", tvcForce(0), tvcForce(1),
            tvcForce(2));

    // Limit the TVC angle to the specified limit. We do this by checking if the
    // vector is outside the limit and if so, we rotate the vector back to the
    // limit.
    auto tvcAngle = tvcForce.getAngleTo(Math::Vector<float, 3>{0, 0, 1});
    auto tvcRotationAxis =
        (tvcForce.cross(Vector<float, 3>{0, 0, 1})).normalize();
    if (tvcAngle > tvcAngleLimit_) {
      auto tvcRotation =
          Math::Quat_F(tvcRotationAxis, tvcAngle - tvcAngleLimit_);
      tvcForce = tvcRotation.conjugate().rotate(tvcForce);
    }

    // Limit the TVC force to the specified limit. We do this by checking if the
    // vector is outside the limit and if so, we scale the vector back to the
    // limit.
    auto tvcMagnitude = tvcForce.magnitude();
    if (tvcMagnitude > tvcForceLimit_) {
      tvcForce = tvcForce.normalize() * tvcForceLimit_;
    }

    tvcInput(0) = tvcForce(0);
    tvcInput(1) = tvcForce(1);
    tvcInput(2) = tvcForce(2);

    totalBodyForce = 0;
    totalBodyTorque = 0;

    tvcZTorque_ = tvcInput(3);

    if (tvcEnabled_) {
      totalBodyForce = totalBodyForce + tvcForce;
    }
  }

  if (Core::NOW() - lastWindUpdateTime_ > 1 * Core::SECONDS) {
    lastWindUpdateTime_ = Core::NOW();
    const float maxWindSpeed = 4.0f;
    const float minWindSpeed = 0.0f;
    const Math::Vector<float, 3> windDirection = {0, 1, 0};
    const Math::Vector<float, 3> windBias = {-1, 0, 0};
    windVelocity = windBias +
                   windDirection.normalize() *
                       Core::randGen(maxWindSpeed - minWindSpeed, Core::NOW()) +
                   minWindSpeed;
  }
  // windVelocity = 0;

  auto attitude = bodySim_.getAttitude();
  auto velocity = bodySim_.getVelocity();
  auto rotationVel = bodySim_.getAngularVelocity();
  auto velocityBody = attitude.conjugate().rotate(velocity);

  bodySim_.clearForces();
  bodySim_.addForce(-Math::GRAVITY_3F, 0, false, true);
  bodySim_.addForce(totalBodyForce, tvcPosition_);
  bodySim_.addTorque({0, 0, tvcZTorque_}, 0);

  if (true) {
    auto velMag = velocity.magnitude();

    auto airFlowVector = -velocityBody - windVelocity;

    const float bodyCylZAxisOffset = 0.1;
    auto afCylTop =
        airFlowVector -
        getLocalAirFlowVector(rotationVel, {0, 0, bodyCylZAxisOffset});
    auto afCylBottom =
        airFlowVector -
        getLocalAirFlowVector(rotationVel, {0, 0, -bodyCylZAxisOffset});

    auto cylBodyForceTop = calcCylForce(afCylTop) / 2;
    auto cylBodyForceBottom = calcCylForce(afCylBottom) / 2;
    auto cylTopForce =
        calcTopBottomForce(airFlowVector, Math::Vector_F{0, 0, 1});
    auto cylBottomForce =
        calcTopBottomForce(airFlowVector, Math::Vector_F{0, 0, -1});

    LOG_MSG("VelBody: %.2f %.2f %.2f |%.2f|\n", airFlowVector(0),
            airFlowVector(1), airFlowVector(2), velMag);

    auto velToPosTensor = Matrix_F<3, 3>{
        0, 0, 0,             //
        0, 0, 0,             //
        0, 0, -1.0f / 200.0f //
    };

    bodySim_.addForce(cylBodyForceTop,
                      velToPosTensor * airFlowVector +
                          Math::Vector_F{0, 0, bodyCylZAxisOffset});
    bodySim_.addForce(cylBodyForceBottom,
                      velToPosTensor * airFlowVector +
                          Math::Vector_F{0, 0, -bodyCylZAxisOffset});

    bodySim_.addForce(cylTopForce, {0, 0, 0.35});
    bodySim_.addForce(cylBottomForce, {0, 0, -0.35});

    const float mountY = 0.08;
    const float mountZ = 0.26;
    const float copY = 0.03;

    auto flapSettings = flapSettingsSubr_.getItem();

    auto tlRot = Quat<float>(Vector_F{0, 0, 1}, -flapSettings.flapTLAngle_Rad);
    auto trRot = Quat<float>(Vector_F{0, 0, 1}, flapSettings.flapTRAngle_Rad);
    auto blRot = Quat<float>(Vector_F{0, 0, 1}, -flapSettings.flapBLAngle_Rad);
    auto brRot = Quat<float>(Vector_F{0, 0, 1}, flapSettings.flapBRAngle_Rad);

    auto tlPos = tlRot.rotate(Vector_F{0.00, copY, 0.0}) +
                 Vector_F{0.00, mountY, mountZ};
    auto trPos = trRot.rotate(Vector_F{0.00, -copY, 0.0}) +
                 Vector_F{0.00, -mountY, mountZ};
    auto blPos = blRot.rotate(Vector_F{0.00, copY, 0.0}) +
                 Vector_F{0.00, mountY, -mountZ};
    auto brPos = brRot.rotate(Vector_F{0.00, -copY, 0.0}) +
                 Vector_F{0.00, -mountY, -mountZ};

    auto tlaf = airFlowVector - getLocalAirFlowVector(rotationVel, tlPos);
    auto traf = airFlowVector - getLocalAirFlowVector(rotationVel, trPos);
    auto blaf = airFlowVector - getLocalAirFlowVector(rotationVel, blPos);
    auto braf = airFlowVector - getLocalAirFlowVector(rotationVel, brPos);

    calcFlapForces(tlaf, traf, blaf, braf);
    auto tlFF = flapForces_[0];
    auto trFF = flapForces_[1];
    auto blFF = flapForces_[2];
    auto brFF = flapForces_[3];

    // LOG_MSG("TL Pos: %.2f %.2f %.2f\n", tlPos(0), tlPos(1), tlPos(2));

    bodySim_.addForce(tlFF, tlPos);
    bodySim_.addForce(trFF, trPos);
    bodySim_.addForce(blFF, blPos);
    bodySim_.addForce(brFF, brPos);

    // bodySim_.addTorque({0, 0, 0.1}, 0, false);*/
  }

  // if (flapsEnabled_)
  // bodySim_.addForce(dragCoeff * finVel, {0, 0, 0.4}); // Add the drag force
  // in reference frame

  // LOG_MSG("Flaps enabled: %d\n", flapsEnabled_); // Print the flaps enabled
  // to the console

  // auto totalForce = bodySim_.getNavForceSum();
  // LOG_MSG("Force: %.2f %.2f %.2f |%.1f|\n", totalForce(0), totalForce(1),
  // totalForce(2), totalForce.magnitude());

  bodySim_.simulateForTime(
      dTime * Core::SECONDS); // Simulate the body simulator for the given time

  positionState_ = bodySim_.getPositionState(); // Get the position state from
                                                // the body simulator
  attitudeState_ = bodySim_.getAttitudeState(); // Get the attitude state from
                                                // the body simulator

  // LOG_MSG("Dtime: %.2f\n", dTime);
  if (zeroingMode_) {
    positionState_ = {0, 0, 0, 1, 0, 0};
    attitudeState_ = {0, 0, 0, 0, 0, 0, 1};
    totalBodyForce = {0, 0, 0};
    totalBodyTorque = {0, 0, 0};
    attitudeState_ = attitudeState_.normalize(3, 7);

    bodySim_.setPositionState(positionState_);
    bodySim_.setAttitudeState(attitudeState_);
  } else {
    // positionState_ = {0, 0, 0, 1, 0, 1};
    // bodySim_.setPositionState(positionState_);
  }

  if (positionState_(5) < 0) {
    if (positionState_(2) < 0)
      positionState_(2) = 0;
    positionState_(0) = 0;
    positionState_(1) = 0;
    positionState_(5) = 0;

    attitudeState_(0) = 0;
    attitudeState_(1) = 0;
    attitudeState_(2) = 0;

    bodySim_.setPositionState(positionState_);
    bodySim_.setAttitudeState(attitudeState_);
  }
}

const VCTR::Math::Vector<float, 7> &BodySimulator::getAttitudeState() {
  return attitudeState_;
}

const VCTR::Math::Vector<float, 6> &BodySimulator::getPositionState() {
  return positionState_;
}

Core::Topic<Core::Timestamped<Math::Vector<float, 7>>> &
BodySimulator::getAttitudeEstTopic() {
  return attitudeTopic_;
}

Core::Topic<Core::Timestamped<Math::Vector<float, 6>>> &
BodySimulator::getStateEstTopic() {
  return positionTopic_;
}

Core::Topic<Core::Timestamped<Math::Vector<float, 2>>> &
BodySimulator::getGroundTopic() {
  return groundTopic_;
}

void BodySimulator::taskInit() {
  stateTimestamp_ =
      Core::NOW(); // Initialize the state timestamp to the current time
}

void BodySimulator::taskThread() {
  // Update the state estimation using any new sensor information
  update();

  // Publish the state estimation to the topic
  Core::Timestamped<Math::Vector<float, 7>> attState;
  attState.timestamp = stateTimestamp_;
  attState.data = attitudeState_;
  Core::Timestamped<Math::Vector<float, 6>> posState;
  posState.timestamp = stateTimestamp_;
  posState.data = positionState_;
  Core::Timestamped<Math::Vector<float, 2>> groundState;
  groundState.timestamp = stateTimestamp_;
  groundState.data = {
      positionState_(2),
      positionState_(5)}; // Formed as: [V, P], where V is the vertical velocity
                          // in m/s and P is the ground distance in meters
  attitudeTopic_.publish(attState);
  positionTopic_.publish(posState);
  groundTopic_.publish(groundState);
}

void BodySimulator::calcFlapForces(const Math::Vector_F &tlVelocity,
                                   const Math::Vector_F &trVelocity,
                                   const Math::Vector_F &blVelocity,
                                   const Math::Vector_F &brVelocity) {
  auto flapSettings = flapSettingsSubr_.getItem();

  //   LOG_MSG("Flap Settings: %.2f %.2f %.2f %.2f\n",
  //           flapSettings.flapTLAngle_Rad / DEGREES,
  //           flapSettings.flapTRAngle_Rad / DEGREES,
  //           flapSettings.flapBLAngle_Rad / DEGREES,
  //           flapSettings.flapBRAngle_Rad / DEGREES);

  // Limits the input value to the given range
  auto limitAngle = [](float &angle, float min, float max) {
    if (angle < min)
      angle = min;
    else if (angle > max)
      angle = max;
  };

  const auto minFlapAngle = 10 * DEGREES;
  const auto maxFlapAngle = 90 * DEGREES;

  limitAngle(flapSettings.flapTLAngle_Rad, minFlapAngle, maxFlapAngle);
  limitAngle(flapSettings.flapTRAngle_Rad, minFlapAngle, maxFlapAngle);
  limitAngle(flapSettings.flapBLAngle_Rad, minFlapAngle, maxFlapAngle);
  limitAngle(flapSettings.flapBRAngle_Rad, minFlapAngle, maxFlapAngle);

  // Lets first calculate the rotations
  auto tlRot =
      Math::Quat<float>(Math::Vector_F{0, 0, 1}, -flapSettings.flapTLAngle_Rad);
  auto trRot =
      Math::Quat<float>(Math::Vector_F{0, 0, 1}, flapSettings.flapTRAngle_Rad);
  auto blRot =
      Math::Quat<float>(Math::Vector_F{0, 0, 1}, -flapSettings.flapBLAngle_Rad);
  auto brRot =
      Math::Quat<float>(Math::Vector_F{0, 0, 1}, flapSettings.flapBRAngle_Rad);

  // Now we calculate the surface normals
  auto tlSN = tlRot.rotate(Math::Vector_F{1, 0, 0});
  auto trSN = trRot.rotate(Math::Vector_F{1, 0, 0});
  auto blSN = blRot.rotate(Math::Vector_F{1, 0, 0});
  auto brSN = brRot.rotate(Math::Vector_F{1, 0, 0});

  // Now we calculate the respective force magnitudes
  const float tfArea = 7.49e-3;
  const float bfArea = 9.12e-3;
  const float lambda =
      1.44; // Air density multiplied by flat surface drag coefficient
  const float adjustmentFactor = 1;

  auto tlForce =
      calcFlapForce(tlSN, tlVelocity, tfArea * lambda * adjustmentFactor);
  auto trForce =
      calcFlapForce(trSN, trVelocity, tfArea * lambda * adjustmentFactor);
  auto blForce =
      calcFlapForce(blSN, blVelocity, bfArea * lambda * adjustmentFactor);
  auto brForce =
      calcFlapForce(brSN, brVelocity, bfArea * lambda * adjustmentFactor);

  flapForces_[0] = tlForce;
  flapForces_[1] = trForce;
  flapForces_[2] = blForce;
  flapForces_[3] = brForce;
}

Math::Vector_F BodySimulator::calcFlapForce(const Math::Vector_F &surfaceNormal,
                                            const Math::Vector_F &airFlowVector,
                                            float basis) {
  float angle = airFlowVector.getAngleTo(surfaceNormal);
  // LOG_MSG("Flap angle: %.2f\n", angle/DEGREES);

  // if (angle < 90*DEGREES)
  //     return 0;
  auto velSurface = airFlowVector.getProjectionOn(surfaceNormal);
  auto velSurfaceMag = velSurface.magnitude();

  return velSurface * (velSurfaceMag * 0.5 * basis);
}

Math::Vector_F
BodySimulator::calcCylForce(const Math::Vector_F &airFlowVector) {
  const float cylSideArea = 0.12 * 0.5; // Side surface assuming flat
  const float cylCoeff = 1.1;

  auto velProj = airFlowVector.getProjectionOn(Math::Vector_F{0, 0, 1}, true);

  return velProj * (0.5 * 1.2 * cylCoeff * cylSideArea * velProj.magnitude());
}

Math::Vector_F
BodySimulator::calcTopBottomForce(const Math::Vector_F &airFlowVector,
                                  const Math::Vector_F &normal) {
  const float cylFlatArea = 2 * M_PI * (0.12 / 2) * (0.12 / 2);
  const float flatCoeff = 1.17;

  auto angle = airFlowVector.getAngleTo(normal);
  if (angle > 90 * DEGREES) {
    return 0;
  }

  auto velProj = airFlowVector.getProjectionOn(normal, false);
  auto velProjDrag = airFlowVector.getProjectionOn(normal, true);

  return velProj * (0.5 * 1.2 * flatCoeff * cylFlatArea * velProj.magnitude());
}

Math::Vector<float, 3>
BodySimulator::calcForceFromTVC(const Math::Vector<float, 4> &tvcInput) {
  return {tvcInput(0), tvcInput(1), tvcInput(2)}; // Force vector in body frame
}

Math::Vector<float, 3>
BodySimulator::calcTorqueFromTVC(const Math::Vector<float, 4> &tvcInput,
                                 const Math::Vector<float, 3> &tvcPosition) {
  Math::Vector<float, 3> forceVector = {tvcInput(0), tvcInput(1), tvcInput(2)};
  Math::Vector<float, 3> torqueVector = forceVector.cross(tvcPosition);
  Math::Vector<float, 3> rollTorque = {0, 0, tvcInput(3)};

  return rollTorque; // Add the roll torque to the torque vector
}

Math::Vector<float, 3>
BodySimulator::getLocalAirFlowVector(const Math::Vector<float, 3> &rotationVel,
                                     const Math::Vector<float, 3> &position) {
  // Calculate the additional air flow vector in local body frame due to
  // rotation
  return rotationVel.cross(position);
}

} // namespace DSP

} // namespace VCTR