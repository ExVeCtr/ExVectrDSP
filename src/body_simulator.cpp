#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/print.hpp"

#include "ExVectrMath.hpp"

#include "ExVectrDSP/value_covariance.hpp"

#include "ExVectrDSP/body_simulator.hpp"

namespace VCTR
{

    namespace DSP
    {

        BodySimulator::BodySimulator(float mass, Math::Vector<float, 3> inertiaTensor, float tvcThrustLimit_N, float tvcAngleLimit_Rad, int64_t simulationInterval) :
            Core::Task_Periodic("Simlation task", simulationInterval),
            bodySim_(mass, inertiaTensor)
        {

            Core::getSystemScheduler().addTask(*this);

            mass_ = mass;
            inertiaTensor_ = Math::Matrix<float, 3, 3>::eye();
            inertiaTensor_(0, 0) = inertiaTensor(0);
            inertiaTensor_(1, 1) = inertiaTensor(1);
            inertiaTensor_(2, 2) = inertiaTensor(2);

            bodySim_.setInertiaTensor(inertiaTensor_); // Set the inertia tensor of the body simulator

            attitudeState_ = {0, 0, 0, 1, 0, 0, 0};
            positionState_ = {0, 0, 0, 0, 0, 0};

            totalBodyForce = Math::GRAVITY_3F * mass_;
            totalBodyTorque = 0;

            tvcPosition_ = {0, 0, -0.35}; // Position of the thrust vector control in body frame

            tvcAngleLimit_ = tvcAngleLimit_Rad; // Set the thrust vector control angle limit
            tvcForceLimit_ = tvcThrustLimit_N; // Set the thrust vector control force limit


        }

        void BodySimulator::update() {

            float dTime = float(Core::NOW() - stateTimestamp_) / Core::SECONDS;
            stateTimestamp_ = Core::NOW();

            // gather force and torque from the thrust vector control input
            if (tvcSubr_.isDataNew()) {

                auto tvcInput = tvcSubr_.getItem();
                Math::Vector<float, 3> tvcVector = {tvcInput(0), tvcInput(1), tvcInput(2)}; //Force vector

                //Limit the TVC angle to the specified limit. We do this by checking if the vector is outside the limit and if so, we rotate the vector back to the limit.
                auto tvcAngle = tvcVector.getAngleTo(Math::Vector<float, 3>({0, 0, 1}));
                auto tvcRotationAxis = (tvcVector.cross(Math::Vector<float, 3>({0, 0, 1}))).normalize(); //Rotation axis is the cross product of the vector and the Z-Axis.
                if (tvcAngle > tvcAngleLimit_) {
                    auto tvcRotation = Math::Quat_F(tvcRotationAxis, tvcAngle - tvcAngleLimit_).conjugate();
                    tvcVector = tvcRotation.rotate(tvcVector); //Rotate the vector back to the limit.
                }

                //Limit the TVC force to the specified limit. We do this by checking if the vector is outside the limit and if so, we scale the vector back to the limit.
                auto tvcMagnitude = tvcVector.magnitude();
                if (tvcMagnitude > tvcForceLimit_) {
                    tvcVector = tvcVector.normalize() * tvcForceLimit_; //Scale the vector back to the limit.
                }

                tvcInput(0) = tvcVector(0);
                tvcInput(1) = tvcVector(1);
                tvcInput(2) = tvcVector(2); //Update the input vector with the limited vector.

                //LOG_MSG("TVC Force: %.2f %.2f %.2f |%.1f|\n", tvcVector(0), tvcVector(1), tvcVector(2), tvcVector.magnitude()); // Print the force vector to the console

                totalBodyForce = 0;
                totalBodyTorque = 0;

                tvcZTorque_ = tvcInput(3)*1.5; // Get the torque around the Z axis from the thrust vector control input in body frame

                auto tvcForce = calcForceFromTVC(tvcInput);
                auto tvcTorque = -calcTorqueFromTVC(tvcInput, tvcPosition_);

                if (tvcEnabled_) {
                    totalBodyForce = totalBodyForce + tvcForce; // Add the force vector to the total force vector
                    totalBodyTorque = totalBodyTorque + tvcTorque; // Add the torque vector to the total torque vector
                }

                tvcAngle = tvcVector.getAngleTo(Math::Vector<float, 3>({0, 0, 1})); // Get the angle of the force vector to the Z-axis

                //LOG_MSG("TVC Force: %.2f %.2f %.2f |%.1f|, Angle: %.2f\n", totalBodyForce(0), totalBodyForce(1), totalBodyForce(2), totalBodyForce.magnitude(), tvcAngle/DEGREES); // Print the force vector and torque to the console

            }

            auto attitude = bodySim_.getAttitude(); // Get the attitude from the body simulator  
            auto angularVelocity = bodySim_.getAngularVelocity(); // Get the angular velocity from the body simulator
            auto velocity = bodySim_.getVelocity(); // Get the velocity from the body simulator
            auto velocityBody = attitude.rotate(velocity);

            //LOG_MSG("AngVel: %.2f %.2f %.2f |%.1f|\n", angularVelocity(0)/DEGREES, angularVelocity(1)/DEGREES, angularVelocity(2)/DEGREES, angularVelocity.magnitude()/DEGREES); // Print the angular velocity to the console

            bodySim_.clearForces();
            bodySim_.addForce(-Math::GRAVITY_3F, 0, false, true);
            bodySim_.addForce(totalBodyForce, tvcPosition_); // Add the tvc force
            bodySim_.addTorque({0, 0, tvcZTorque_}, 0); // Add the tvc twist torque in body frame
            
            float angVelDragCoeff = 0.01; // Coefficient for the angular velocity drag force
            bodySim_.setAngularVelocity(angularVelocity * (1 - angVelDragCoeff));

            if (velocity(2) < -1 || true) { // We are moving down. Start simulating the pseudo drag force for the belly flop.

                const float lambda = 0.0436; // drag pseudo coefficient
                auto velMag = velocity.magnitude(); // Get the magnitude of the velocity vector
                auto dragForce = (-velocity) * (lambda * velMag); // Calculate the drag force in body frame

                dragForce = attitude.rotate(dragForce); // Rotate the drag force to the body frame

                //LOG_MSG("Drag force: %.2f %.2f %.2f |%.1f|\n", dragForce(0), dragForce(1), dragForce(2), dragForce.magnitude()); // Print the drag force to the console

                bodySim_.addForce(dragForce, {0.05, 0, 0}, true);

                /*auto cylBodyForce = calcCylForce(velocityBody);
                auto cylTopForce = calcTopBottomForce(velocityBody, Math::Vector_F({0, 0, 1}));
                auto cylBottomForce = calcTopBottomForce(velocityBody, Math::Vector_F({0, 0, -1}));


                LOG_MSG("Vel: %.2f %.2f %.2f |%.1f|\n", velocityBody(0), velocityBody(1), velocityBody(2), velocityBody.magnitude());
                LOG_MSG("Force: %.2f %.2f %.2f |%.1f|\n", cylBottomForce(0), cylBottomForce(1), cylBottomForce(2), cylBottomForce.magnitude());

                cylBodyForce = attitude.rotate(cylBodyForce);
                cylTopForce = attitude.rotate(cylTopForce);
                cylBottomForce = attitude.rotate(cylBottomForce);

                bodySim_.addForce(cylBodyForce, {0, 0, 0}, true);
                bodySim_.addForce(cylTopForce, {0, 0, 0.35}, false);
                bodySim_.addForce(cylBottomForce, {0, 0, -0.35}, true);*/

                /*calcFlapForces(velocityBody);

                auto tlFF = attitude.rotate(flapForces_[0]);
                auto trFF = attitude.rotate(flapForces_[1]);
                auto blFF = attitude.rotate(flapForces_[2]);
                auto brFF = attitude.rotate(flapForces_[3]);

                bodySim_.addForce(tlFF, {0.05, 0.02, 0.35});
                bodySim_.addForce(trFF, {-0.05, 0.02, 0.35});
                bodySim_.addForce(blFF, {0.05, 0.02, -0.35});
                bodySim_.addForce(brFF, {-0.05, 0.02, -0.35});

                //bodySim_.addTorque({0, 0, 0.1}, 0, false);*/

            }

            //if (flapsEnabled_)
                //bodySim_.addForce(dragCoeff * finVel, {0, 0, 0.4}); // Add the drag force in reference frame

            //LOG_MSG("Flaps enabled: %d\n", flapsEnabled_); // Print the flaps enabled to the console

            //auto totalForce = bodySim_.getNavForceSum();
            //LOG_MSG("Force: %.2f %.2f %.2f |%.1f|\n", totalForce(0), totalForce(1), totalForce(2), totalForce.magnitude());

            bodySim_.simulateForTime(dTime * Core::SECONDS); // Simulate the body simulator for the given time

            positionState_ = bodySim_.getPositionState(); // Get the position state from the body simulator
            attitudeState_ = bodySim_.getAttitudeState(); // Get the attitude state from the body simulator

            
            //LOG_MSG("Dtime: %.2f\n", dTime);
            if (zeroingMode_) {
                //LOG_MSG("Zeroing mode enabled\n");
                positionState_ = {0, 0, 0, 1, 0, 0}; // Set the position and velocity to 0
                attitudeState_ = {0, 0, 0, 1, 0, 0, 0}; // Set the attitude to the reference frame
                totalBodyForce = {0, 0, 0}; // Set the force to 0
                totalBodyTorque = {0, 0, 0}; // Set the torque to 0
                //zeroingMode_ = false; // Disable the zeroing mod
                attitudeState_ = attitudeState_.normalize();

                bodySim_.setPositionState(positionState_); // Set the position state to the body simulator
                bodySim_.setAttitudeState(attitudeState_); // Set the attitude state to the body simulator

            }

            if (positionState_(5) < 0) {

                if (positionState_(2) < 0) // Only limit the velocity to 0 if the position is below 0. We want to be able to go back up.
                    positionState_(2) = 0;
                positionState_(0) = 0;
                positionState_(1) = 0;
                positionState_(5) = 0;

                attitudeState_(0) = 0; // Set the X velocity to 0
                attitudeState_(1) = 0; // Set the Y velocity to 0
                attitudeState_(2) = 0; // Set the Z velocity to 0

                bodySim_.setPositionState(positionState_); // Set the position state to the body simulator
                bodySim_.setAttitudeState(attitudeState_); // Set the attitude state to the body simulator

            }


        }

        const VCTR::Math::Vector<float, 7> &BodySimulator::getAttitudeState() {
            return attitudeState_;
        }

        const VCTR::Math::Vector<float, 6> &BodySimulator::getPositionState() {
            return positionState_;
        }


        Core::Topic<Core::Timestamped<Math::Vector<float, 7>>>& BodySimulator::getAttitudeEstTopic() {
            return attitudeTopic_;
        }

        Core::Topic<Core::Timestamped<Math::Vector<float, 6>>>& BodySimulator::getStateEstTopic() {
            return positionTopic_;
        }


        void BodySimulator::taskInit() {
            stateTimestamp_ = Core::NOW(); // Initialize the state timestamp to the current time
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
            attitudeTopic_.publish(attState);
            positionTopic_.publish(posState);
        }

        void BodySimulator::calcFlapForces(const Math::Vector_F& velocity) {

            auto flapSettings = flapSettingsSubr_.getItem();

            // Lets setup some values we will need 


            // Lets first calculate the rotations
            auto tlRot = Math::Quat<float>(Math::Vector_F({0, 0, 1}), -flapSettings.flapTLAngle_Rad);
            auto trRot = Math::Quat<float>(Math::Vector_F({0, 0, 1}), flapSettings.flapTRAngle_Rad);
            auto blRot = Math::Quat<float>(Math::Vector_F({0, 0, 1}), -flapSettings.flapBLAngle_Rad);
            auto brRot = Math::Quat<float>(Math::Vector_F({0, 0, 1}), flapSettings.flapBRAngle_Rad);

            //Now we calculate the surface normals
            auto tlSN = tlRot.rotate(Math::Vector_F({1, 0, 0}));
            auto trSN = trRot.rotate(Math::Vector_F({1, 0, 0}));
            auto blSN = blRot.rotate(Math::Vector_F({1, 0, 0}));
            auto brSN = brRot.rotate(Math::Vector_F({1, 0, 0}));

            //Now we calculate the respective force magnitudes
            const float tfArea = 2.49e-3;
            const float bfArea = 9.12e-3;
            const float lambda = 1.44; //Air density multiplied by flat surface drag coefficient

            auto tlForce = calcFlapForce(tlSN, velocity, tfArea*lambda);
            auto trForce = calcFlapForce(trSN, velocity, tfArea*lambda);
            auto blForce = calcFlapForce(blSN, velocity, bfArea*lambda);
            auto brForce = calcFlapForce(brSN, velocity, bfArea*lambda);

            flapForces_[0] = tlForce;
            flapForces_[1] = trForce;
            flapForces_[0] = blForce;
            flapForces_[1] = brForce;

        }

        Math::Vector_F BodySimulator::calcFlapForce(const Math::Vector_F& surfaceNormal, const Math::Vector_F& velocity, float basis) {

            float angle = velocity.getAngleTo(surfaceNormal);

            //if (angle > 90*DEGREES)
            //    return 0;
            auto velSurface = velocity.getProjectionOn(surfaceNormal);
            auto velSurfaceMag = velSurface.magnitude();

            return velSurface * (velSurfaceMag * 0.5 * basis);

        }

        Math::Vector_F BodySimulator::calcCylForce(const Math::Vector_F& velocityBody) {

            const float cylSideArea = 0.12*0.5; //Side surface assuming flat
            const float cylCoeff = 1.17;

            auto velProj = velocityBody.getProjectionOn(Math::Vector_F({0, 0, 1}), true);

            return velProj * (0.5 * 1.2 * cylCoeff * cylSideArea * velProj.magnitude());

        }

        Math::Vector_F BodySimulator::calcTopBottomForce(const Math::Vector_F& velocityBody, const Math::Vector_F& normal) {

            const float cylFlatArea = 2 * M_PI * (0.12/2) * (0.12/2); //Top bottom surfaces
            const float flatCoeff = 1.17;

            auto angle = velocityBody.getAngleTo(normal);
            if (angle > 90*DEGREES)
                return 0;

            auto velProj = velocityBody.getProjectionOn(normal, false);// * 0.7 + velocityBody*0.3;

            return velProj * (0.5 * 1.2 * flatCoeff * cylFlatArea * velProj.magnitude());

        }

        Math::Vector<float, 3> BodySimulator::calcForceFromTVC(const Math::Vector<float, 4> &tvcInput) {

            return {tvcInput(0), tvcInput(1), tvcInput(2)}; //Force vector in body frame

        }

        Math::Vector<float, 3> BodySimulator::calcTorqueFromTVC(const Math::Vector<float, 4> &tvcInput, const Math::Vector<float, 3>& tvcPosition) {

            Math::Vector<float, 3> forceVector = {tvcInput(0), tvcInput(1), tvcInput(2)}; //Force vector in body frame
            Math::Vector<float, 3> torqueVector = forceVector.cross(tvcPosition); //Torque vector in body frame
            Math::Vector<float, 3> rollTorque = {0, 0, tvcInput(3)}; //Roll torque in body frame. 

            return torqueVector - rollTorque; //Add the roll torque to the torque vector

        }


    }

}