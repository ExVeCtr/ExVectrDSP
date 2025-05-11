#ifndef EXVECTRDSP_SIMULATOR_HPP
#define EXVECTRDSP_SIMULATOR_HPP

#include "ExVectrCore/task_types.hpp"

#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"

#include "ExVectrMath.hpp"

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

            //Where TVC input is published. In form: [X, Y, Z, T], where X, Y, Z show the thrust vector in body frame (magnitude of vector is thrust magnitude) and T is the roll torque (Z-Axis) angle in radians.
            Core::Simple_Subscriber<Math::Vector<float, 4>> tvcSubr_;
            Math::Vector<float, 3> tvcPosition_;
            float tvcAngleLimit_ = 10 * 3.14/180; // 30 degrees in radians
            float tvcForceLimit_ = 12; // 1000 N

            /// @brief Timestamp of state estimation.
            int64_t stateTimestamp_;

            bool zeroingMode_ = false; // When enabled, then the filter attempts to zero the position and velocity, the reference position is set to the current position and finally the barometer and accelerometer are calibrated assuming no motion.
            
            Math::Vector<float, 7> attitudeState_; // Formed as: [V, Q], where V is the angular velocity in body frame and Q is a unit quaternion rotation from the reference frame to body frame.
            Math::Vector<float, 6> positionState_; // Formed as: [V, P], where V is the current velocity, P is the current position in reference frame

            Math::Vector<float, 3> totalBodyForce; // Total force in body frame
            Math::Vector<float, 3> totalBodyTorque; // Total torque in body frame

            float mass_ = 1.0f; // Mass of the vehicle in kg
            Math::Matrix<float, 3, 3> inertiaTensor_ = Math::Matrix<float, 3, 3>::eye(); // Inertia matrix of the vehicle in kg*m^2

        public:
            /**
             * @brief Standard constructor. Sets state to 0 and covariance to 1000 as starting values.
             */
            BodySimulator(float mass, Math::Vector<float, 3> inertiaTensor, int64_t simulationInterval = 0.01*Core::SECONDS);

            /**
             * @brief Subsribes to a topic to which the thrust vector control output is published. In form: [X, Y, Z, T], where X, Y, Z show the thrust vector in body frame (magnitude of vector is thrust magnitude) and T is the roll torque (Z-Axis) angle in radians.
             * @param tvcTopic The topic to which the thrust vector control output is published.
             * @param tvcPosition The position of the thrust vector control in body frame.
             */
            void setTVCInputTopic(Core::Topic<Math::Vector<float, 4>> &tvcTopic, const Math::Vector<float, 3>& tvcPosition, float tvcAngleLimit, float tvcForceLimit) {
                tvcAngleLimit_ = tvcAngleLimit;
                tvcForceLimit_ = tvcForceLimit;
                tvcSubr_.subscribe(tvcTopic);
                tvcPosition_ = tvcPosition;
            }

            /**
             * @brief Updates the current attitude estimation using any new sensor information.
             */
            void update();

            /**
             * @brief WHen enabled then the position and velocity are set to 0
             */
            void enableZeroingMode(bool enable) { zeroingMode_ = enable; }


            /**
             * @note Current state estimation. Formed as: [W, Q], where W is the angular velocity in body frame and Q is a unit quaternion rotation from the reference frame to body frame.
             * @returns estimation for state.
             */
            VCTR::Math::Vector<float, 7> &getAttitudeState();

            /**
             * @note Current state estimation. Formed as: [V, P], where V is the current velocity, P is the current position in reference frame (X-North, Y-West, Z-Up).
             * @returns estimation for state.
             */
            VCTR::Math::Vector<float, 6> &getPositionState();


            /**
             * @return The attitude state estimation from reference to body frame in form: [W, Q].
             */
            Core::Topic<Core::Timestamped<Math::Vector<float, 7>>>& getAttitudeEstTopic();

            /**
             * @return The state estimation for velocity and position in form: [V, P].
             */
            Core::Topic<Core::Timestamped<Math::Vector<float, 6>>>& getStateEstTopic();

            /**
             * @return The attitude state estimation cov in form: [W, Q].
             */
            //Core::Topic<Core::Timestamped<Math::Matrix<float, 7, 7>>>& getAttitudeCovTopic();


        private:

            void taskInit() override;

            void taskThread() override;

            Math::Vector<float, 3> calcForceFromTVC(const Math::Vector<float, 4> &tvcInput);

            /// @brief Calculates the torque from the thrust vector control input. The torque is calculated as a cross product of the thrust vector and the distance from the center of gravity to the thrust vector.
            /// @param tvcInput The thrust vector control input in form: [X, Y, Z, T], where X, Y, Z show the thrust vector in body frame (magnitude of vector is thrust magnitude) and T is the roll torque (Z-Axis) angle in radians.
            /// @param tvcPosition The position of the thrust vector control in body frame.
            Math::Vector<float, 3> calcTorqueFromTVC(const Math::Vector<float, 4> &tvcInput, const Math::Vector<float, 3>& tvcPosition);

        };

    }
    
}

#endif