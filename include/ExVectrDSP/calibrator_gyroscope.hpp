#ifndef ExVectrDSP_CALIBRATORGYROSCOPE_HPP
#define ExVectrDSP_CALIBRATORGYROSCOPE_HPP

#include "ExVectrCore/task_types.hpp"

#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"

#include "ExVectrData/memory_interface.hpp"
#include "ExVectrData/memory_manager.hpp"

#include "ExVectrMath.hpp"

#include "value_covariance.hpp"
#include "topic_coord_transform.hpp"
#include "calibrator_interface.hpp"


namespace VCTR
{

    namespace DSP
    {

        /**
         * @brief A class implementing a simple averaging algorithm for gyroscope calibration. The Model for calibration is this: sensor = (raw - bias) * transform. With raw values directly from the sensor and sensor values are corrected for bias, scale and rotation.
         */
        class Calibrator_Gyroscope : public Calibrator_Interface<float, 3>
        {
        protected:

            Core::Callback_Subscriber<Core::Timestamped<ValueCov<float, 3>>, Calibrator_Gyroscope> gyroSubr_;

            Math::Vector_F valueSum_ = 0;
            uint32_t valueCounter_ = 0;

            int64_t calibStartTime_ = 0;
            int64_t calibDuration_ = 0;

        public:

            Calibrator_Gyroscope(int64_t calibDuration = 10*Core::SECONDS);

            /**
             * * @brief Starts the calibration process. Can be used to retrigger the calibration process if needed.
             */
            void initCalibration() override;

            /**
             * @brief Sets the gyro input topic for the calibration.
             * @note The gyro input should be in [tesla] in sensor frame.
             * @param gyroTopic Topic for gyroscope input.
             */
            void setGyroInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &gyroTopic, const Math::Matrix<float, 3, 3>& gyroRotation);

            /**
             * @brief Updates the current attitude estimation using any new sensor information.
             */
            void update(Core::Timestamped<ValueCov<float, 3>> const &gyroData);

        };

    }
    
}

#endif