#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/print.hpp"

#include "ExVectrMath.hpp"

#include "ExVectrDSP/value_covariance.hpp"

#include "ExVectrDSP/body_simulator.hpp"

namespace VCTR
{

    namespace DSP
    {

        BodySimulator::BodySimulator(float mass, Math::Vector<float, 3> inertiaTensor, int64_t simulationInterval) :
            Core::Task_Periodic("Simlation task", simulationInterval)
        {

            Core::getSystemScheduler().addTask(*this);

            mass_ = mass;
            inertiaTensor_ = Math::Matrix<float, 3, 3>::eye();
            inertiaTensor_(0, 0) = inertiaTensor(0);
            inertiaTensor_(1, 1) = inertiaTensor(1);
            inertiaTensor_(2, 2) = inertiaTensor(2);

            attitudeState_ = {0, 0, 0, 1, 0, 0, 0};
            positionState_ = {0, 0, 0, 0, 0, 0};

            totalBodyForce = Math::GRAVITY_3F * mass_;
            totalBodyTorque = 0;

            tvcPosition_ = {0, 0, -0.35}; // Position of the thrust vector control in body frame

        }

        void BodySimulator::update() {

            float dTime = float(Core::NOW() - stateTimestamp_) / Core::SECONDS;
            stateTimestamp_ = Core::NOW();

            float hdtsq = 0.5 * dTime * dTime;

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
                    tvcVector = tvcVector * (tvcForceLimit_ / tvcMagnitude); //Scale the vector back to the limit.
                }

                tvcInput(0) = tvcVector(0);
                tvcInput(1) = tvcVector(1);
                tvcInput(2) = tvcVector(2); //Update the input vector with the limited vector.

                //LOG_MSG("TVC Force: %.2f %.2f %.2f\n", tvcInput(0), tvcInput(1), tvcInput(2)); // Print the force vector to the console

                auto tvcForce = calcForceFromTVC(tvcInput);
                auto tvcTorque = -calcTorqueFromTVC(tvcInput, tvcPosition_);

                totalBodyForce = tvcForce;
                totalBodyTorque = tvcTorque;

            }


            
            //LOG_MSG("Dtime: %.2f\n", dTime);
            if (zeroingMode_) {
                //LOG_MSG("Zeroing mode enabled\n");
                positionState_ = {0, 0, 0, 0, 0, 0}; // Set the position and velocity to 0
                attitudeState_ = {0, 0, 0, 1, 0, 0, 0}; // Set the attitude to the reference frame
                totalBodyForce = {0, 0, 0}; // Set the force to 0
                totalBodyTorque = {0, 0, 0}; // Set the torque to 0
                //zeroingMode_ = false; // Disable the zeroing mode

                attitudeState_ = attitudeState_.normalize();
            }

            //totalBodyForce = {0, 0, 9.81};

            //positionState_(2) = 1;

            // Calculate the position state transition and input models
            Math::Matrix<float, 6, 6> F = {
                1, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0,
                dTime, 0, 0, 1, 0, 0,
                0, dTime, 0, 0, 1, 0,
                0, 0, dTime, 0, 0, 1
            };

            Math::Matrix<float, 6, 3> B = {
                dTime / mass_, 0, 0,
                0, dTime / mass_, 0,
                0, 0, dTime / mass_,
                hdtsq, 0, 0,
                0, hdtsq, 0,
                0, 0, hdtsq
            };

            // Calculate total force in reference frame from body frame
            Math::Quat_F attQuat = attitudeState_.block<4, 1>(3, 0); // Quaternion in body frame
            //Math::Vector<float, 3> additionalForce = {0, 1, 0}; // Additional force in body frame
            Math::Vector<float, 3> totalBodyForceRef = attQuat.conjugate().rotate(totalBodyForce) - Math::GRAVITY_3F * mass_;// + additionalForce; // Rotate the force vector to reference frame
            //LOG_MSG("Total body force: %.2f %.2f %.2f\n", totalBodyForceRef(0), totalBodyForceRef(1), totalBodyForceRef(2)); // Print the force vector to the console
            // Calculate the new position state
            positionState_ = F * positionState_ + B * totalBodyForceRef;
            //LOG_MSG("Accel input: %.2f %.2f %.2f %.2f %.2f %.2f\n", stateAccelInput(0), stateAccelInput(1), stateAccelInput(2), stateAccelInput(3), stateAccelInput(4), stateAccelInput(5)); // Print the acceleration input to the console
            //LOG_MSG("Position state: %.2f %.2f %.2f %.2f %.2f %.2f\n", positionState_(0), positionState_(1), positionState_(2), positionState_(3), positionState_(4), positionState_(5)); // Print the position state to the console
            // Do a ground collision check. If the Z-Position is below 0, set it to 0 and set the velocity z axis to 0.
            if (positionState_(5) < 0) {

                if (positionState_(2) < 0) // Only limit the velocity to 0 if the position is below 0. We want to be able to go back up.
                    positionState_(2) = 0;
                positionState_(0) = 0;
                positionState_(1) = 0;
                positionState_(5) = 0;

                attitudeState_(0) = 0; // Set the X velocity to 0
                attitudeState_(1) = 0; // Set the Y velocity to 0
                attitudeState_(2) = 0; // Set the Z velocity to 0

            }

            // ##### Testing purposes only #####
            //positionState_ = {0, 0, 0, 0, 0, 0.5}; // Set the position and velocity to 0
            //positionState_(0) = 0; // Set the X position to 0
            //positionState_(1) = 0; // Set the Y position to 0
            //positionState_(2) = 0; // Set the Z position to 0
            //positionState_(3) = 0; // Set the Z position to 0
            //positionState_(4) = 0; // Set the Z velocity to 0
            //positionState_(5) = 0.5; // Set the Z velocity to 0
            // Retrieve the angular velocity and quaternion from the attitude state
            Math::Vector<float, 3> angVel = attitudeState_.block<3, 1>(0, 0); // Angular velocity in body frame
            /*if (rand() % 100 == 0) {
                angVel = angVel + Math::Vector<float, 3>({0, 0.01, 0}); // Add some noise to the angular velocity
                LOG_MSG("TORQUE KICK!\n");
            }*/
            //Math::Quat_F attQuat = attitudeState_.block<4, 1>(3, 0); // Quaternion in body frame

            // Update the attitude and angular velocity
            attQuat = attQuat * Math::Quat_F(angVel.normalize(), angVel.magnitude() * dTime); // Update quaternion using angular velocity
            angVel(0) += totalBodyTorque(0) / inertiaTensor_(0, 0) * dTime; // Update angular velocity using torque and inertia tensor
            angVel(1) += totalBodyTorque(1) / inertiaTensor_(1, 1) * dTime; // Update angular velocity using torque and inertia tensor
            angVel(2) += totalBodyTorque(2) / inertiaTensor_(2, 2) * dTime; // Update angular velocity using torque and inertia tensor
            attQuat.normalize(); // Normalize the quaternion

            //Update the attitude state
            attitudeState_.block(angVel, 0, 0, 0, 0);
            attitudeState_.block(attQuat, 3, 0, 0, 0); // Update the quaternion in the state vector

            // ###### Testing purposes only #####
            //attitudeState_ = {0, 0, 0, 1, 0, 0, 0};


            //LOG_MSG("Attitude: %.2f %.2f %.2f %.2f\n", attitudeState_(3), attitudeState_(4), attitudeState_(5), attitudeState_(6));

        }

        VCTR::Math::Vector<float, 7> &BodySimulator::getAttitudeState() {
            return attitudeState_;
        }

        VCTR::Math::Vector<float, 6> &BodySimulator::getPositionState() {
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