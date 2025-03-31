#ifndef ExVectrDSP_IMUGPSPOSITIONKF_H
#define ExVectrDSP_IMUGPSPOSITIONKF_H

#include "ExVectrCore/task_types.hpp"

#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"

#include "ExVectrMath.hpp"

#include "value_covariance.hpp"
#include "gnss_data.hpp"

namespace VCTR
{

    namespace Data
    {

        /**
         * @brief A class implementing a standard linear kalman filter for the position using IMU and GPS accel and position sensors.
         * @note The internal state vector is formed as: [V, P], where V is the current velocity, P is the current position in reference frame (X-North, Y-West, Z-Up).
         */
        class IMUGPSPositionKalman
        {
        protected:

            /// @brief Current state estimation. Formed as: [V, P], where V is the current velocity, P is the current position in reference frame (X-North, Y-West, Z-Up).
            VCTR::Math::Vector<float, 6> x_;
            /// @brief Current state estimation covariance.
            VCTR::Math::Matrix<float, 6, 6> p_;

            ///@brief Rotation transform from acc sensor frame to body frame.
            Math::Matrix<float, 3, 3> accRot_;

            Core::Simple_Subscriber<Core::Timestamped<ValueCov<float, 7>>> attSubr_;
            //Core::Simple_Subscriber<Core::Timestamped<ValueCov<float, 3>>> accSubr_;
            Core::Buffer_Subscriber<Core::Timestamped<ValueCov<float, 3>>, 5> accSubr_;
            Core::Simple_Subscriber<Core::Timestamped<ValueCov<float, 1>>> baroSubr_;
            Core::Simple_Subscriber<Core::Timestamped<Data::GNSSData>> gnssSubr_;

            Core::Timestamped<ValueCov<float, 3>> lastAccData_;
            Core::Timestamped<ValueCov<float, 1>> lastBaroData_;

            float seaLevelPressure_ = 101325;

            bool attInitialised_ = false;
            bool accInitialised_ = false;
            bool baroInitialised_ = false;
            bool gnssInitialised_ = false;

            Math::Vector<double, 3> gnssRefPos_ = 0;

            Math::Quat_F att_;

            Core::Topic<Core::Timestamped<ValueCov<float, 6>>> stateEstTopic_;


        public:
            /**
             * @brief Standard constructor. Sets state to 0 and covariance to 1000 as starting values.
             */
            IMUGPSPositionKalman();

            /**
             * @brief Predicts system state upto given time.
             * @param time Defaults to current time.
             */
            void predict(int64_t time = VCTR::Core::NOW());

            /**
             * @brief Sets the attitude input topic.
             * @note The attitude input must be in form: [W, Q], where W is the angular velocity [rad/s] in body frame and Q is a unit quaternion rotation from the reference frame to body frame.
             * @param attitudeTopic Topic for attitude estimate input.
             */
            void setAttitudeInput(Core::Topic<Core::Timestamped<ValueCov<float, 7>>> &attitudeTopic);

            /**
             * @brief Sets the acc input topic.
             * @note The acc input should be raw accelerometer data in [m/s^2] in body frame.
             * @param accTopic Topic for accelerometer input.
             * @param accRot Rotation transform from acc sensor frame to body frame.
             */
            void setAccInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &accTopic, const VCTR::Math::Matrix<float, 3, 3> &accRot = 1);

            /**
             * @brief Sets the barometer input topic.
             * @note The baro input should be in [pa]. The Covariance in meters.
             * @param baroTopic Topic for barometer input.
             * @param seaLevelPressure Reference sea level pressure. This is estimated when using a GPS input.
             */
            void setBaroInput(Core::Topic<Core::Timestamped<ValueCov<float, 1>>> &baroTopic, float seaLevelPressure = 101325);

            /**
             * @brief Sets the GNSS input topic.
             * @note The GNSS input should be in Form: [V, P], where V is the velocity vector (NED) in m/s and P is the position in Form: [Lat, Lon, Alt].
             *       The X, Y values correct the position and velocity. Z is used to correct the barometer.
             * @param gnssTopic Topic for barometer input.
             * @param refPos Reference position for the GNSS data in form [Latitude, Longitude, Altitude]. If not given, the first GNSS data is used as the reference position.
             */
            void setGNSSInput(Core::Topic<Core::Timestamped<Data::GNSSData>> &gnssTopic, Math::Vector<double, 3> refPos = 0);

            /**
             * @brief Sets the reference position for the GNSS data. This is used to calculate the relative position from the GNSS data.
             * @param refPos Reference position for the GNSS data in form [Latitude, Longitude, Altitude].
             */
            void setPositionReference(const Math::Vector<double, 3>& refPos) { gnssRefPos_ = refPos; }

            /**
             * @brief Updates the current attitude estimation using any new sensor information.
             */
            void update();

            /**
             * @brief Updates the current attitude estimation.
             * @param attitudeData 
             */
            void updateAttitude(const Core::Timestamped<VCTR::Data::ValueCov<float, 7>> & attitudeData);

            /**
             * @brief Updates the current attitude estimation using an obvervation of the acceleration in body frame.
             * @param accData
             */
            void updateAcc(const VCTR::Core::Timestamped<VCTR::Data::ValueCov<float, 3U>>& accData);

            /**
             * @brief Updates the current altitude estimation (Z-Axis) using the air pressure measured.
             * @param baroData
             */
            void updateBaro(const VCTR::Core::Timestamped<VCTR::Data::ValueCov<float, 1U>>& baroData); 

            /**
             * @brief Updates the current position estimation using the GNSS data.
             * @param gnssData
             */
            void updateGNSS(const VCTR::Core::Timestamped<VCTR::Data::GNSSData>& gnssData);

            /**
             * @brief Initialises the state estimation using an obvervation of the acceleration in body frame.
             * @param accData
             */
            void initialiseAcc(const VCTR::Core::Timestamped<VCTR::Data::ValueCov<float, 3U>>& accData);

            /**
             * @brief Initialises the state estimation using an obvervation of the air pressure.
             * @param baroData
             */
            void initialiseBaro(const VCTR::Core::Timestamped<VCTR::Data::ValueCov<float, 1U>>& baroData);

            /**
             * @brief Initialises the state estimation using an obvervation of the GNSS data. If otherwise not given, the first GNSS data is used as the reference position (Zero point).
             * @param gnssData
             */
            void initialiseGNSS(const VCTR::Core::Timestamped<VCTR::Data::GNSSData>& gnssData);

            /**
             * Sets the expected process noise covariance for the model prediction.
             * @param Noise Covariance matrix.
             */
            void setProcessNoise(const VCTR::Math::Matrix<float, 6, 6> &noise);

            /**
             * @note The internal state vector is formed as: [V, Q], where V is the angular velocity vector and Q is a unit quaternion rotation from the reference frame to body frame.
             * @returns estimation for state.
             */
            const VCTR::Math::Vector<float, 6> &getState();

            /**
             * @note The internal state vector is formed as: [V, Q], where V is the angular velocity vector and Q is a unit quaternion rotation from the reference frame to body frame.
             * @returns estimation for covariance.
             */
            const VCTR::Math::Matrix<float, 6, 6> &getCovariance();

            /**
             * @brief Sets the current estimation for the state and covariance.
             * @note The internal state vector is formed as: [V, Q], where V is the angular velocity vector and Q is a unit quaternion rotation from the reference frame to body frame.
             * @param state The state vector and covariance to be used.
             */
            void setState(const ValueCov<float, 6> &state);

            /**
             * @return The state estimation for velocity and position in form: [V, P].
             */
            Core::Topic<Core::Timestamped<ValueCov<float, 6>>>& getStateEstTopic();

            /// @returns altutide calculated from the given current air pressure and reference sea level pressure.
            float calcAltitudeFromPressure(float pressure, float seaLevelPressure = 101325);

            /// @returns sea level pressure calculated from the given current air pressure and altitude.
            float calcSealevelPressFromAltitude(float pressure, float altitude);

            /// @returns relative position from the given GNSS data and reference position in a NED frame.
            Math::Vector_F calcRelPosFromGNSS(const Math::Vector<double, 3> &gnssData, const Math::Vector<double, 3> &refPos);

        };

        /**
         * @brief This combines the IMUAttitudeEKF with a periodic task to update the state automatically
         * @note The internal state vector is formed as: [V, Q], where V is the angular velocity vector and Q is a unit quaternion rotation from the reference frame to body frame.
         */
        class IMUGPSPositionKalmanTask : public IMUGPSPositionKalman, public Core::Task_Periodic {
        public:

            IMUGPSPositionKalmanTask(int64_t period, Core::Scheduler& scheduler = Core::getSystemScheduler());

            void taskCheck() override;

            void taskInit() override;

            void taskThread() override;

        };

    }
    
}

#endif