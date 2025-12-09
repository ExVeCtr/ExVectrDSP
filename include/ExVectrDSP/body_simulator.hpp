#ifndef EXVECTRDSP_SIMULATOR_HPP
#define EXVECTRDSP_SIMULATOR_HPP

#include "ExVectrMath.hpp"

#include "ExVectrCore/task_types.hpp"

#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"

#include "ExVectrControl/physics_bodysim.hpp"
#include "ExVectrControl/starship/control_attitude_flaps.hpp"

#include "value_covariance.hpp"

namespace VCTR
{

    namespace DSP
    {

        /**
         * @brief A class implementing a quaternion based EKF to estimate an rotation from world to body or aka attitude from IMU data.
         * @note The internal state vector is formed as: [B, Q], where B is the gyro bias in sensor frame and Q is a unit quaternion rotation from the reference frame to body frame.
         */
        class BodySimulator : public Core::Task_Periodic
        {
        protected:
            Core::Topic<Core::Timestamped<Math::Vector<float, 7>>> attitudeTopic_;
            Core::Topic<Core::Timestamped<Math::Vector<float, 6>>> positionTopic_;
            Core::Topic<Core::Timestamped<Math::Vector<float, 2>>> groundTopic_;

            // Where TVC input is published. In form: [X, Y, Z, T], where X, Y, Z show the thrust vector in body frame (magnitude of vector is thrust magnitude) and T is the roll torque (Z-Axis) angle in radians.
            Core::Simple_Subscriber<Math::Vector<float, 4>> tvcSubr_;
            Math::Vector<float, 3> tvcPosition_;
            float tvcAngleLimit_ = 8 * 3.14 / 180; // 30 degrees in radians
            float tvcForceLimit_ = 12;             // 1000 N

            Core::Simple_Subscriber<CTRL::ControlAttitudeFlapSetting> flapSettingsSubr_;
            Math::Vector<float, 3> flapForces_[4]; // tl, tr, bl, br

            /// @brief Timestamp of state estimation.
            int64_t stateTimestamp_;

            bool zeroingMode_ = false; // When enabled, then the filter attempts to zero the position and velocity, the reference position is set to the current position and finally the barometer and accelerometer are calibrated assuming no motion.

            Math::Vector<float, 7> attitudeState_; // Formed as: [V, Q], where V is the angular velocity in body frame and Q is a unit quaternion rotation from the reference frame to body frame.
            Math::Vector<float, 6> positionState_; // Formed as: [V, P], where V is the current velocity, P is the current position in reference frame

            CTRL::PhysicsBodySim bodySim_; // The body simulator object that simulates the forces and torques on the vehicle.

            Math::Vector<float, 3> totalBodyForce;  // Total force in body frame
            Math::Vector<float, 3> totalBodyTorque; // Total torque in body frame

            float tvcZTorque_ = 0; // The torque around the Z axis from the thrust vector control input in body frame

            float mass_ = 1.0f;                                                          // Mass of the vehicle in kg
            Math::Matrix<float, 3, 3> inertiaTensor_ = Math::Matrix<float, 3, 3>::eye(); // Inertia matrix of the vehicle in kg*m^2

            bool tvcEnabled_ = true;    // If the actuators are enabled or not. This is set to true when the actuators are enabled.
            bool flapsEnabled_ = false; // If the actuators are active or not. This is set to true when the actuators are active.

        public:
            /**
             * @brief Standard constructor. Sets state to 0 and covariance to 1000 as starting values.
             */
            BodySimulator(float mass, Math::Vector<float, 3> inertiaTensor, float tvcThrustLimit_N, float tvcAngleLimit_Rad, int64_t simulationInterval = 0.01 * Core::SECONDS);

            /**
             * @brief Subsribes to a topic to which the thrust vector control output is published. In form: [X, Y, Z, T], where X, Y, Z show the thrust vector in body frame (magnitude of vector is thrust magnitude) and T is the roll torque (Z-Axis) angle in radians.
             * @param tvcTopic The topic to which the thrust vector control output is published.
             * @param tvcPosition The position of the thrust vector control in body frame.
             */
            void setTVCInputTopic(Core::Topic<Math::Vector<float, 4>> &tvcTopic, const Math::Vector<float, 3> &tvcPosition, float tvcAngleLimit, float tvcForceLimit)
            {
                tvcAngleLimit_ = tvcAngleLimit;
                tvcForceLimit_ = tvcForceLimit;
                tvcSubr_.subscribe(tvcTopic);
                tvcPosition_ = tvcPosition;
            }

            void setFlapSettingInputTopic(Core::Topic<CTRL::ControlAttitudeFlapSetting> &flapTopic)
            {
                flapSettingsSubr_.subscribe(flapTopic);
            }

            /**
             * @brief Updates the current attitude estimation using any new sensor information.
             */
            void update();

            /**
             * @brief WHen enabled then the position and velocity are set to 0
             */
            void enableZeroingMode(bool enable) { zeroingMode_ = enable; }

            void enableTVC(bool enable) { tvcEnabled_ = enable; } // Enable or disable the thrust vector control

            void enableFlaps(bool enable) { flapsEnabled_ = enable; } // Enable or disable the thrust vector control

            /**
             * @note Current state estimation. Formed as: [W, Q], where W is the angular velocity in body frame and Q is a unit quaternion rotation from the reference frame to body frame.
             * @returns estimation for state.
             */
            const VCTR::Math::Vector<float, 7> &getAttitudeState();

            /**
             * @note Current state estimation. Formed as: [V, P], where V is the current velocity, P is the current position in reference frame (X-North, Y-West, Z-Up).
             * @returns estimation for state.
             */
            const VCTR::Math::Vector<float, 6> &getPositionState();

            /**
             * @note Current ground distance estimation. Formed as: [V, P], where V is the vertical velocity in m/s and P is the ground distance in meters.
             */
            const VCTR::Math::Vector<float, 2> &getGroundState();

            /**
             * @return The attitude state estimation from reference to body frame in form: [W, Q].
             */
            Core::Topic<Core::Timestamped<Math::Vector<float, 7>>> &getAttitudeEstTopic();

            /**
             * @return The state estimation for velocity and position in form: [V, P].
             */
            Core::Topic<Core::Timestamped<Math::Vector<float, 6>>> &getStateEstTopic();

            /**
             * @return The ground distance estimation in form: [V, P], where V is the vertical velocity in m/s and P is the ground distance in meters.
             */
            Core::Topic<Core::Timestamped<Math::Vector<float, 2>>> &getGroundTopic();

            void setAttitudeState(const Math::Vector<float, 7> &attitudeState) { bodySim_.setAttitudeState(attitudeState); } // Set the attitude state in reference frame
            void setPositionState(const Math::Vector<float, 6> &positionState) { bodySim_.setPositionState(positionState); } // Set the position state in reference frame

            /**
             * @return The attitude state estimation cov in form: [W, Q].
             */
            // Core::Topic<Core::Timestamped<Math::Matrix<float, 7, 7>>>& getAttitudeCovTopic();

        private:
            void taskInit() override;

            void taskThread() override;

            void calcFlapForces(const Math::Vector_F &velocity);
            Math::Vector_F calcFlapForce(const Math::Vector_F &surfaceNormal, const Math::Vector_F &velocityBody, float basis);

            Math::Vector_F calcCylForce(const Math::Vector_F &velocityBody);
            Math::Vector_F calcTopBottomForce(const Math::Vector_F &velocityBody, const Math::Vector_F &normal);

            Math::Vector<float, 3> calcForceFromTVC(const Math::Vector<float, 4> &tvcInput);

            /// @brief Calculates the torque from the thrust vector control input. The torque is calculated as a cross product of the thrust vector and the distance from the center of gravity to the thrust vector.
            /// @param tvcInput The thrust vector control input in form: [X, Y, Z, T], where X, Y, Z show the thrust vector in body frame (magnitude of vector is thrust magnitude) and T is the roll torque (Z-Axis) angle in radians.
            /// @param tvcPosition The position of the thrust vector control in body frame.
            Math::Vector<float, 3> calcTorqueFromTVC(const Math::Vector<float, 4> &tvcInput, const Math::Vector<float, 3> &tvcPosition);
        };

    }

}

#endif