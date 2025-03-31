#ifndef ExVectrDSP_IMUATTITUDEEKF_H
#define ExVectrDSP_IMUATTITUDEEKF_H

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
        class IMUAttitudeEKF
        {
        protected:
            /// @brief Timestamp of state estimation.
            //int64_t stateTimestamp_;
            /// @brief Current state estimation. Formed as: [B, Q], where B is the gyro bias in sensor frame and Q is a unit quaternion rotation from the reference frame to body frame.
            VCTR::Math::Vector<float, 7> x_;
            /// @brief Current state estimation covariance.
            VCTR::Math::Matrix<float, 7, 7> p_;
            /// @brief Timestamp of last gyro update.
            int64_t gyroUpdateTimestamp_ = 0;

            ///@brief Rotation transform from gyro sensor frame to body frame.
            Math::Matrix<float, 3, 3> gyroRot_;
            ///@brief Rotation transform from acc sensor frame to body frame.
            Math::Matrix<float, 3, 3> accRot_;
            ///@brief Rotation transform from mag sensor frame to body frame.
            Math::Matrix<float, 3, 3> magTransform_; Math::Matrix<float, 3, 1> magBias_;

            Core::Simple_Subscriber<Core::Timestamped<ValueCov<float, 3>>> gyroSubr_;
            Core::Simple_Subscriber<Core::Timestamped<ValueCov<float, 3>>> accSubr_;
            Core::Simple_Subscriber<Core::Timestamped<ValueCov<float, 3>>> magSubr_;

            bool gyroInitialised_ = false;
            bool accInitialised_ = false;
            bool magInitialised_ = false;

            Core::Timestamped<ValueCov<float, 3>> lastGyroData_;
            Core::Timestamped<ValueCov<float, 3>> lastAccData_;
            Core::Timestamped<ValueCov<float, 3>> lastMagData_;
            
            /// @brief The topic to which the bias estimations are published. The bias is in body frame.
            Core::Topic<Core::Timestamped<ValueCov<float, 3>>> biasTopic_;

            /// @brief The topic to which the attitude estimations are published. Formed as: [W, Q], where W is the angular velocity in body frame and Q is a unit quaternion rotation from the reference frame to body frame.
            Core::Topic<Core::Timestamped<ValueCov<float, 7>>> attitudeTopic_;

        public:
            /**
             * @brief Standard constructor. Sets state to 0 and covariance to 1000 as starting values.
             */
            IMUAttitudeEKF();

            /**
             * @brief Predicts system state upto given time.
             * @param time Defaults to current time.
             */
            void predict(int64_t time = VCTR::Core::NOW());

            /**
             * @brief Sets the gyro input topic for the IMUAttitudeEKF.
             * @note The gyro input should be in [rad/s] in body frame.
             * @param gyroTopic Topic for gyroscope input.
             * @param gyroRot Rotation transform from gyro sensor frame to body frame.
             */
            void setGyroInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &gyroTopic, const VCTR::Math::Matrix<float, 3, 3> &gyroRot = 1);

            /**
             * @brief Sets the acc input topic for the IMUAttitudeEKF.
             * @note The acc input should be in [m/s^2] in body frame.
             * @param accTopic Topic for accelerometer input.
             * @param accRot Rotation transform from acc sensor frame to body frame.
             */
            void setAccInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &accTopic, const VCTR::Math::Matrix<float, 3, 3> &accRot = 1);

            /**
             * @brief Sets the mag input topic for the IMUAttitudeEKF.
             * @note The mag input should be in [tesla] in body frame. Sensor calibration model: mag = magTransform * (worldMag - magBias)
             * @param magTopic Topic for magnetometer input.
             * @param magTransform Rotation transform from mag sensor frame to body frame and also scaling.
             * @param magBias Bias of magnetometer.
             */
            void setMagInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &magTopic, const VCTR::Math::Matrix<float, 3, 3> &magTransform = 1, const VCTR::Math::Matrix<float, 3, 1> &magBias = 0);

            /**
             * @brief Updates the current attitude estimation using any new sensor information.
             */
            void update();

            /**
             * @brief Updates the current attitude estimation using an obvervation of the angular velocity in body frame.
             * @param gyroData 
             */
            void updateGyro(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>>& gyroData);

            /**
             * @brief Updates the current attitude estimation using an obvervation of the acceleration in body frame.
             * @param accData
             */
            void updateAcc(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>>& accData);

            /**
             * @brief Updates the current attitude estimation using an obvervation of the magnetic field in body frame.
             * @param magData
             */
            void updateMag(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>>& magData);

            /**
             * @brief Initialises the state estimation using an obvervation of the acceleration in body frame.
             * @param accData
             */
            void initialiseAcc(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>>& accData);

            /**
             * @brief Initialises the state estimation using an obvervation of the magnetic field in body frame.
             * @param magData
             */
            void initialiseMag(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>>& magData);

            /**
             * Sets the expected process noise covariance for the model prediction.
             * @param Noise Covariance matrix.
             */
            void setProcessNoise(const VCTR::Math::Matrix<float, 7, 7> &noise);

            /**
             * @note Current state estimation. Formed as: [W, Q], where W is the angular velocity in body frame and Q is a unit quaternion rotation from the reference frame to body frame.
             * @returns estimation for state.
             */
            const VCTR::Math::Vector<float, 7> &getState();

            /**
             * @note Current state estimation covariance. Formed as: [W, Q], where W is the angular velocity in body frame and Q is a unit quaternion rotation from the reference frame to body frame.
             * @returns estimation for covariance.
             */
            const VCTR::Math::Matrix<float, 7, 7> &getCovariance();

            /**
             * @brief Sets the current estimation for the state and covariance.
             * @note The internal state vector is formed as: [V, Q], where V is the angular velocity vector and Q is a unit quaternion rotation from the reference frame to body frame.
             * @param state The state vector and covariance to be used.
             */
            void setState(const ValueCov<float, 7> &state);

            /**
             * @return The gyro bias vector in body frame topic.
             */
            Core::Topic<Core::Timestamped<ValueCov<float, 3>>>& getBiasEstTopic();

            /**
             * @return The attitude state estimation from reference to body frame in form: [W, Q].
             */
            Core::Topic<Core::Timestamped<ValueCov<float, 7>>>& getAttitudeEstTopic();

        };

        /**
         * @brief This combines the IMUAttitudeEKF with a periodic task to update the state automatically
         * @note The internal state vector is formed as: [V, Q], where V is the angular velocity vector and Q is a unit quaternion rotation from the reference frame to body frame.
         */
        class IMUAttitudeEKFTask : public IMUAttitudeEKF, public Core::Task_Periodic {
        public:

            IMUAttitudeEKFTask(int64_t period, Core::Scheduler& scheduler = Core::getSystemScheduler());

            void taskCheck() override;

            void taskInit() override;

            void taskThread() override;

        };

    }
    
}

#endif