#include "ExVectrDSP/body_sim_model.hpp"

#include <math.h>

namespace VCTR {
namespace DSP {

namespace {
constexpr float kGravity = 9.807f; // = Math::GRAVITY
constexpr float kPi = 3.14159265f;

inline float clampf(float v, float lo, float hi) { return v < lo ? lo : (v > hi ? hi : v); }

/// Moves `cur` toward `target` by at most `maxStep`.
inline float slew(float cur, float target, float maxStep) {
  const float d = target - cur;
  if (d > maxStep)
    return cur + maxStep;
  if (d < -maxStep)
    return cur - maxStep;
  return target;
}
} // namespace

BodySimModel::BodySimModel(const BodySimParams &params)
    : params_(params), body_(params.mass_kg, Math::Matrix<float, 3, 3>::eye()) {
  thrustCmd_ = {0, 0, 0};
  loads_.aeroForce = {0, 0, 0};
  loads_.aeroTorque = {0, 0, 0};
  loads_.airVelocity = {0, 0, 0};
  for (int i = 0; i < 4; i++)
    loads_.flapForce[i] = {0, 0, 0};
  applyParams();
  resetToStartPose();
}

void BodySimModel::applyParams() {
  const BodySimParams &p = params_;
  body_.setMass(p.mass_kg);
  Math::Matrix<float, 3, 3> I = Math::Matrix<float, 3, 3>::eye();
  I(0, 0) = p.inertiaXX;
  I(1, 1) = p.inertiaYY;
  I(2, 2) = p.inertiaZZ;
  I(0, 1) = I(1, 0) = p.inertiaXY;
  I(0, 2) = I(2, 0) = p.inertiaXZ;
  I(1, 2) = I(2, 1) = p.inertiaYZ;
  body_.setInertiaTensor(I);

  zBottom_ = -p.cgFromBottom_m;
  zNoseBase_ = p.tubeLength_m - p.cgFromBottom_m;
  zTip_ = zNoseBase_ + p.noseLength_m;
  // Slender-body potential normal force of a nose acts at x = L (1 - V_nose / (A_base L)) from the tip; for the
  // parabolic nose r = R (2s - s^2): V_nose = 8/15 A_base L.
  zNoseCp_ = zTip_ - p.noseLength_m * (1.0f - 8.0f / 15.0f);
  baseArea_ = kPi * p.bodyRadius_m * p.bodyRadius_m;

  const BodySimFlapGeometry *g[2] = {&p.flapTop, &p.flapBottom};
  for (int k = 0; k < 2; k++) {
    const float w = g[k]->width_m, h = g[k]->height_m, b = g[k]->bottomFromShipBottom_m;
    float c = g[k]->chamfer_m / sqrtf(2.0f); // legs of the 45 deg cut
    c = c < w ? c : w;
    c = c < h ? c : h;
    const float aRect = w * h, aTri = 0.5f * c * c;
    flapArea_[k] = aRect - aTri;
    // Centroid along the axis: rectangle minus the corner triangle at the nose-side outer corner.
    const float zc = (aRect * (b + 0.5f * h) - aTri * (b + h - c / 3.0f)) / flapArea_[k];
    flapZc_[k] = zc - p.cgFromBottom_m;
  }
}

void BodySimModel::resetToStartPose() {
  body_.setAttitudeState({0, 0, 0, 1, 0, 0, 0});
  body_.setPositionState({0, 0, 0, 0, 0, 0});
}

void BodySimModel::setThrust(const Math::Vector<float, 3> &forceBody, float rollTorque, bool enabled) {
  Math::Vector<float, 3> f = forceBody;
  if (!(isfinite(f(0)) && isfinite(f(1)) && isfinite(f(2))))
    f = {0, 0, 0};
  // Magnitude limit.
  float mag = f.magnitude();
  if (mag > params_.tvcThrustLimit_N) {
    f = f * (params_.tvcThrustLimit_N / mag);
    mag = params_.tvcThrustLimit_N;
  }
  // Cone limit about body +Z: keep the azimuth, put the vector on the cone (the old clamp rotated away from +Z).
  const float lat = sqrtf(f(0) * f(0) + f(1) * f(1));
  if (mag > 1e-6f && atan2f(lat, f(2)) > params_.tvcAngleLimit_rad && lat > 1e-9f) {
    const float s = sinf(params_.tvcAngleLimit_rad) * mag / lat;
    f = {f(0) * s, f(1) * s, cosf(params_.tvcAngleLimit_rad) * mag};
  }
  thrustCmd_ = f;
  thrustMag_ = mag;
  rollTorque_ = isfinite(rollTorque) ? rollTorque : 0.0f;
  thrustEnabled_ = enabled;
}

void BodySimModel::setFlapInput(float tl, float tr, float bl, float br) {
  const float in[4] = {tl, tr, bl, br};
  for (int i = 0; i < 4; i++)
    if (isfinite(in[i]))
      flapInput_[i] = clampf(in[i], params_.flapMinAngle_rad, params_.flapMaxAngle_rad);
}

void BodySimModel::setWind(float speed_mps, float direction_rad) {
  if (isfinite(speed_mps))
    params_.windSpeed_mps = speed_mps < 0.0f ? 0.0f : speed_mps;
  if (isfinite(direction_rad))
    params_.windDirection_rad = direction_rad;
}

Math::Vector<float, 3> BodySimModel::windWorld() const {
  return {params_.windSpeed_mps * cosf(params_.windDirection_rad),
          params_.windSpeed_mps * sinf(params_.windDirection_rad), 0.0f};
}

float BodySimModel::flapCn(float alpha) const {
  const float sa = sinf(alpha), ca = cosf(alpha);
  const float sigma = 1.0f / (1.0f + expf(-(alpha - params_.flapStallAngle_rad) / params_.flapStallWidth_rad));
  return (1.0f - sigma) * params_.flapCnAlpha * sa * ca + sigma * params_.flapCn90 * sa;
}

float BodySimModel::bodyDiameterAt(float z) const {
  const float d = 2.0f * params_.bodyRadius_m;
  if (z <= zNoseBase_)
    return d;
  float s = (zTip_ - z) / params_.noseLength_m;
  s = clampf(s, 0.0f, 1.0f);
  return d * (2.0f * s - s * s);
}

void BodySimModel::updateActuators(float dt) {
  const BodySimParams &p = params_;

  // Flaps: dead time -> rate limit -> first-order lag.
  int delaySteps = int(p.flapServoDeadTime_s / p.step_s + 0.5f);
  delaySteps = delaySteps < 0 ? 0 : (delaySteps > kDelayLen - 1 ? kDelayLen - 1 : delaySteps);
  if (!delayPrimed_) {
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < kDelayLen; j++)
        flapDelay_[i][j] = flapInput_[i];
    delayPrimed_ = true;
  }
  for (int i = 0; i < 4; i++)
    flapDelay_[i][delayIdx_] = flapInput_[i];
  const int readIdx = (delayIdx_ - delaySteps + kDelayLen) % kDelayLen;
  delayIdx_ = (delayIdx_ + 1) % kDelayLen;

  for (int i = 0; i < 4; i++) {
    const float target = flapDelay_[i][readIdx];
    flapRateLimited_[i] = p.applyActuatorRateLimits ? slew(flapRateLimited_[i], target, p.flapServoRate_radps * dt)
                                                    : target;
    if (p.flapServoLag_s > 0.0f)
      flapAngle_[i] += (flapRateLimited_[i] - flapAngle_[i]) * (dt / (p.flapServoLag_s + dt));
    else
      flapAngle_[i] = flapRateLimited_[i];
  }

  // TVC: slew of the thrust angles = servo rate / servo-to-thrust gain (servo space, review F6).
  const float tgtX = atan2f(thrustCmd_(1), thrustCmd_(2));
  const float tgtY = atan2f(thrustCmd_(0), thrustCmd_(2));
  if (p.applyActuatorRateLimits && thrustMag_ > 1e-6f) {
    const float maxStep = p.tvcServoRate_radps / p.tvcServoToThrustGain * dt;
    tvcAngleX_ = slew(tvcAngleX_, tgtX, maxStep);
    tvcAngleY_ = slew(tvcAngleY_, tgtY, maxStep);
  } else if (thrustMag_ > 1e-6f) {
    tvcAngleX_ = tgtX;
    tvcAngleY_ = tgtY;
  }
}

void BodySimModel::addFlap(int index, float angle, const Math::Vector<float, 3> &u0,
                           const Math::Vector<float, 3> &omega) {
  const BodySimParams &p = params_;
  const int pair = index < 2 ? 0 : 1;
  const BodySimFlapGeometry &g = pair == 0 ? p.flapTop : p.flapBottom;
  const float side = (index == 0 || index == 2) ? 1.0f : -1.0f; // TL, BL on +Y
  const float ca = cosf(angle), sa = sinf(angle);
  // TL/BL rotate by -angle about Z, TR/BR by +angle: normal n = Rz(-+a) x, radial d = Rz(-+a) (+-y).
  const Math::Vector<float, 3> n = {ca, -side * sa, 0.0f};
  const Math::Vector<float, 3> d = {sa, side * ca, 0.0f};
  const Math::Vector<float, 3> hinge = {0.0f, side * p.flapHingeRadius_m, flapZc_[pair]};

  const Math::Vector<float, 3> pc = hinge + d * (0.5f * g.width_m);
  const Math::Vector<float, 3> u = u0 - omega.cross(pc);
  const float speed2 = u.dot(u);
  if (speed2 < 1e-6f) {
    loads_.flapForce[index] = {0, 0, 0};
    return;
  }
  const float speed = sqrtf(speed2);
  const float un = u.dot(n), ud = u.dot(d);
  const float alpha = asinf(clampf(fabsf(un) / speed, 0.0f, 1.0f));
  const float fmag = 0.5f * p.airDensity_kgm3 * p.aeroCalibration * g.effectiveness * flapArea_[pair] * speed2 *
                     flapCn(alpha);
  const Math::Vector<float, 3> F = n * (un >= 0.0f ? fmag : -fmag);
  // Chordwise centre of pressure from the leading edge: 1/4 chord at small alpha, 1/2 at 90 deg.
  const float xLead = g.width_m * (0.5f - 0.25f * cosf(alpha));
  const float xcp = ud >= 0.0f ? xLead : g.width_m - xLead; // flow hinge -> tip: hinge edge leads
  const Math::Vector<float, 3> pcp = hinge + d * xcp;

  loads_.flapForce[index] = F;
  loads_.aeroForce = loads_.aeroForce + F;
  loads_.aeroTorque = loads_.aeroTorque + pcp.cross(F);
}

void BodySimModel::computeAero() {
  const BodySimParams &p = params_;
  const Math::Quat_F q = body_.getAttitude();
  const Math::Vector<float, 3> v = body_.getVelocity();
  const Math::Vector<float, 3> omega = body_.getAngularVelocity();
  const Math::Vector<float, 3> u0 = q.conjugate().rotate(windWorld() - v); // air velocity at the CoM, body
  const float k = 0.5f * p.airDensity_kgm3 * p.aeroCalibration;

  loads_.airVelocity = u0;
  loads_.aeroForce = {0, 0, 0};
  loads_.aeroTorque = {0, 0, 0};

  // Body crossflow strips (with the rotation-induced flow: distributed damping).
  int n = p.bodyStrips < 4 ? 4 : (p.bodyStrips > 64 ? 64 : p.bodyStrips);
  const float len = zTip_ - zBottom_, dz = len / float(n);
  float fx = 0, fy = 0, tx = 0, ty = 0;
  for (int i = 0; i < n; i++) {
    const float z = zBottom_ + (float(i) + 0.5f) * dz;
    // u = u0 - omega x r, r = (0, 0, z): omega x r = (wy z, -wx z, 0)
    const float ux = u0(0) - omega(1) * z;
    const float uy = u0(1) + omega(0) * z;
    const float up = sqrtf(ux * ux + uy * uy);
    const float c = k * p.crossflowEta * p.cdCrossflow * bodyDiameterAt(z) * dz * up;
    const float sfx = c * ux, sfy = c * uy;
    fx += sfx;
    fy += sfy;
    tx += -z * sfy; // r x F = (0,0,z) x (Fx,Fy,0) = (-z Fy, z Fx, 0)
    ty += z * sfx;
  }
  loads_.aeroForce = loads_.aeroForce + Math::Vector<float, 3>{fx, fy, 0.0f};
  loads_.aeroTorque = loads_.aeroTorque + Math::Vector<float, 3>{tx, ty, 0.0f};

  // Slender-body potential normal force at the leading end, and axial force.
  const float uPerp = sqrtf(u0(0) * u0(0) + u0(1) * u0(1));
  const float speed = sqrtf(uPerp * uPerp + u0(2) * u0(2));
  if (speed > 1e-3f && uPerp > 1e-6f && p.potentialFactor != 0.0f) {
    const float alpha = atan2f(uPerp, fabsf(u0(2))); // 0 = axial flow, 90 deg = broadside
    const float fp = k * p.potentialFactor * baseArea_ * speed * speed * sinf(2.0f * alpha) * cosf(0.5f * alpha);
    const float zLead = u0(2) < 0.0f ? zNoseCp_ : zBottom_; // air from the nose (u_z < 0): nose leads
    const float px = fp * u0(0) / uPerp, py = fp * u0(1) / uPerp;
    loads_.aeroForce = loads_.aeroForce + Math::Vector<float, 3>{px, py, 0.0f};
    loads_.aeroTorque = loads_.aeroTorque + Math::Vector<float, 3>{-zLead * py, zLead * px, 0.0f};
  }
  const float cdAx = u0(2) > 0.0f ? p.cdAxialTailFirst : p.cdAxialNoseFirst;
  loads_.aeroForce(2) += k * cdAx * baseArea_ * u0(2) * fabsf(u0(2)); // on the axis: no moment

  for (int i = 0; i < 4; i++)
    addFlap(i, flapAngle_[i], u0, omega);
}

void BodySimModel::groundContact(const Math::Vector<float, 6> &posBefore, const Math::Vector<float, 7> &attBefore) {
  Math::Vector<float, 6> &pos = body_.getPositionState();
  if (pos(5) >= 0.0f)
    return;
  // Ground (height 0 of the CoM): stop, no sliding, no rotation, no tip-over (simplest contact model). A body that
  // was already resting keeps its previous horizontal position and attitude, so wind or residual torques cannot
  // move it along the ground; a body arriving from the air stops where it hits.
  const bool wasResting = posBefore(5) <= 0.0f;
  pos(5) = 0.0f;
  pos(0) = 0.0f;
  pos(1) = 0.0f;
  if (pos(2) < 0.0f)
    pos(2) = 0.0f;
  Math::Vector<float, 7> &att = body_.getAttitudeState();
  if (wasResting) {
    pos(3) = posBefore(3);
    pos(4) = posBefore(4);
    for (int i = 3; i < 7; i++)
      att(i) = attBefore(i);
  }
  att(0) = 0.0f;
  att(1) = 0.0f;
  att(2) = 0.0f;
}

void BodySimModel::step() {
  const float dt = params_.step_s;
  updateActuators(dt);
  computeAero();

  body_.clearForces();
  body_.addForce({0.0f, 0.0f, -kGravity}, {0, 0, 0}, false, true); // gravity, world frame, as acceleration
  if (thrustEnabled_ && thrustMag_ > 1e-6f) {
    Math::Vector<float, 3> f = thrustCmd_;
    if (params_.applyActuatorRateLimits) // rebuild the vector from the slewed angles (as StarshipTVC does)
      f = {sinf(tvcAngleY_) * thrustMag_, sinf(tvcAngleX_) * thrustMag_,
           cosf(tvcAngleX_) * cosf(tvcAngleY_) * thrustMag_};
    body_.addForce(f, {0.0f, 0.0f, params_.tvcFromCg_m});
    body_.addTorque({0.0f, 0.0f, rollTorque_}); // gated with the thrust (review F13)
  }
  body_.addForce(loads_.aeroForce);
  body_.addTorque(loads_.aeroTorque);
  const Math::Vector<float, 6> posBefore = body_.getPositionState();
  const Math::Vector<float, 7> attBefore = body_.getAttitudeState();
  body_.step(dt);
  groundContact(posBefore, attBefore);
  time_s_ += double(dt);
}

int BodySimModel::advance(float elapsed_s) {
  if (!(elapsed_s > 0.0f) || !isfinite(elapsed_s))
    return 0;
  if (elapsed_s > params_.maxAdvance_s)
    elapsed_s = params_.maxAdvance_s;
  carry_s_ += double(elapsed_s);
  int steps = 0;
  const double h = double(params_.step_s);
  while (carry_s_ >= h - 1e-6) { // 1 us tolerance: float step and float elapsed
    step();
    carry_s_ -= h;
    steps++;
  }
  if (carry_s_ < 0.0)
    carry_s_ = 0.0;
  return steps;
}

} // namespace DSP
} // namespace VCTR
