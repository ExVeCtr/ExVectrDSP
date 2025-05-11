#ifndef ExVectrDSP_CALIBRATORMAGNETOMETER_HPP
#define ExVectrDSP_CALIBRATORMAGNETOMETER_HPP

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
         * @brief A class implementing a simple min max algorithm for magnetometer calibration. The Model for calibration is this: mag = (magraw - bias) * transform. With Magraw values directly from the sensor and mag values are corrected for bias, scale and rotation.
         */
        class Calibrator_Magnetometer : public Calibrator_Interface<float, 3>
        {
        protected:

            Core::Callback_Subscriber<Core::Timestamped<ValueCov<float, 3>>, Calibrator_Magnetometer> magSubr_;

            Math::Vector_F magMax_ = -1000;
            Math::Vector_F magMin_ = 1000;

            Math::Vector_F magBias_ = 0;
            Math::Matrix<float, 3, 3> magScale_ = Math::Matrix<float, 3, 3>::eye();


        public:

            void initCalibration() override;

            /**
             * @brief Sets the mag input topic for the calibration.
             * @note The mag input should be in [tesla] in sensor frame.
             * @param magTopic Topic for magnetometer input.
             */
            void setMagInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &magTopic, const Math::Matrix<float, 3, 3>& magRotation);

            /**
             * @brief Updates the current attitude estimation using any new sensor information.
             */
            void update(Core::Timestamped<ValueCov<float, 3>> const &magData);

        };

    }
    
}

#endif