#ifndef EXVECTRDSP_CALIBRATORINTERFACE_HPP
#define EXVECTRDSP_CALIBRATORINTERFACE_HPP


#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"

#include "ExVectrMath/matrix_base.hpp"

#include "ExVectrDSP/value_covariance.hpp"
#include "ExVectrDSP/topic_coord_transform.hpp"


namespace VCTR
{

    namespace DSP
    {

        /**
         * * @brief Interface class for interfacing between sensor calibrators and the rest of the system. Often interface between sensor calibration, topic transfrom and EEPROM storage.
         * * @tparam TYPE Type of the data. (float, double, etc.)
         * * @tparam DIM Dimension of the data. (1, 3, 6, etc.) Generally 3 for 3D sensors (IMU)
         * * @note The output of sensor calibration uses the ValueCov struct, but the value is the bias and the covariance is the scale/rotation matrix. Value = (Raw - bias) * scale/rotation. Raw is directly from sensor and value is calibrated sensor value. 
         */
        template<typename TYPE, size_t DIM = 3>
        class Calibrator_Interface
        {
        protected:

            /// @brief The topic to which calibration data is published.
            Core::Timestamped<ValueCov<TYPE, DIM>> calibrationData_;

            /// @brief This matrix transforms the sensor data from the sensor frame to the body frame. It is used to transform the sensor data to the body frame.
            Math::Matrix<TYPE, DIM, DIM> sensorToBodyMat_;

            /// @brief Flag to indicate if the sensor is calibrated or not. This is set to true when the calibration process is complete.
            bool isCalibrated_ = false; 
            

        public:

            virtual ~Calibrator_Interface() = default;

            /**
             * @brief Starts the calibration process. Can be used to retrigger the calibration process if needed.
             */
            virtual void initCalibration() = 0;

            /**
             * @returns true if the calibration process is complete, false otherwise. The calibration might continue, but at this point the calibration is good.
             */
            bool getIsCalibrated() const {
                return isCalibrated_;
            }

            /**
             * @brief Sets the transformation matrix from the sensor frame to the body frame.
             */
            void setSensorToBodyMat(const Math::Matrix<TYPE, DIM, DIM>& sensorToBodyMat) {
                sensorToBodyMat_ = sensorToBodyMat;
            }

            /**
             * @brief Sets the calibration data. Some calibration systems possibly do not support being set this way, so this is optional.
             */
            virtual void setCalibrationData(const Core::Timestamped<ValueCov<TYPE, DIM>>& calibData) {
                isCalibrated_ = true;
                calibrationData_ = calibData;
            }

            /**
             * @brief This is where calibration data will be published to.
             */
            Core::Timestamped<ValueCov<TYPE, DIM>>& getCalibrationData() {
                return calibrationData_;
            }



        };

    }

}











#endif