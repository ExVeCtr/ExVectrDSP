#ifndef ExVectrDSP_MAGNETOMETERCALIBRATOR_H
#define ExVectrDSP_MAGNETOMETERCALIBRATOR_H

#include "ExVectrCore/task_types.hpp"

#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"

#include "ExVectrMath.hpp"

#include "value_covariance.hpp"

namespace VCTR
{

    namespace Data
    {

        /**
         * @brief A class implementing a simple min max algorithm for magnetometer calibration. The Model for calibration is this: mag = (magraw - bias) * transform. With Magraw values directly from the sensor and mag values are corrected for bias and scale.
         * @note The internal state vector is formed as: [B, Q], where B is the gyro bias in sensor frame and Q is a unit quaternion rotation from the reference frame to body frame.
         */
        class MagnetometerCalibrator
        {
        protected:

            Core::Simple_Subscriber<Core::Timestamped<ValueCov<float, 3>>> magSubr_;

            Math::Vector_F magMax_ = -1000;
            Math::Vector_F magMin_ = 1000;

            Math::Vector_F magBias_ = 0;
            Math::Matrix<float, 3, 3> magTransform_ = Math::Matrix<float, 3, 3>::eye();


        public:
            /**
             * @brief Standard constructor. Sets state to 0 and covariance to 1000 as starting values.
             */
            MagnetometerCalibrator();

            /**
             * @brief Sets the mag input topic for the calibration.
             * @note The mag input should be in [tesla] in sensor frame.
             * @param magTopic Topic for magnetometer input.
             */
            void setMagInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &magTopic);

            /**
             * @brief Updates the current attitude estimation using any new sensor information.
             */
            void update();

            /**
             * @brief Gets the bias estimation for the magnetometer. (Bias is without the transform)
             * @returns Bias estimation.
             */
            const Math::Vector_F& getBias() { return magBias_; };

            /**
             * @brief Gets the transform estimation for the magnetometer.
             * @returns Transform estimation.
             */
            const Math::Matrix<float, 3, 3>& getTransform() { return magTransform_; };

        };

        /**
         * @brief This combines the IMUAttitudeEKF with a periodic task to update the state automatically
         * @note The internal state vector is formed as: [V, Q], where V is the angular velocity vector and Q is a unit quaternion rotation from the reference frame to body frame.
         */
        class MagnetometerCalibratorTask : public MagnetometerCalibrator, public Core::Task_Periodic {
        public:

            MagnetometerCalibratorTask(int64_t period, Core::Scheduler& scheduler = Core::getSystemScheduler());

            void taskCheck() override;

            void taskInit() override;

            void taskThread() override;

        };

    }
    
}

#endif