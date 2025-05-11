#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/print.hpp"

#include "ExVectrMath.hpp"

#include "ExVectrDSP/value_covariance.hpp"

#include "ExVectrDSP/calibrator_gyroscope.hpp"


namespace VCTR
{

    namespace DSP
    {

        /*MagnetometerCalibrator::MagnetometerCalibrator() 
        {
            //memoryManager_.setMemory(memory, memoryKey);

            magMax_ = -1000;
            magMin_ = 1000;
            magBias_ = 0;
            magTransform_ = Math::Matrix<float, 3, 3>::eye();
        }*/

        Calibrator_Gyroscope::Calibrator_Gyroscope(int64_t calibDuration) : 
            calibDuration_(calibDuration) 
        {
        }

        void Calibrator_Gyroscope::initCalibration() {
            valueCounter_ = 0;
            valueSum_ = 0;
            calibStartTime_ = Core::NOW();
            isCalibrated_ = false;
        }

        void Calibrator_Gyroscope::setGyroInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &magTopic, const Math::Matrix<float, 3, 3>& magRotation)
        {
            gyroSubr_.subscribe(magTopic);
            gyroSubr_.setCallback(this, &Calibrator_Gyroscope::update);

            sensorToBodyMat_ = magRotation;
        }

        void Calibrator_Gyroscope::update(Core::Timestamped<ValueCov<float, 3>> const &gyroData)
        {

            if (isCalibrated_) {
                return;
            } 
                
            if (calibStartTime_ == 0) {
                calibStartTime_ = Core::NOW();
            }   

            valueSum_ = valueSum_ + gyroData.data.val;
            valueCounter_++;

            if (Core::NOW() - calibStartTime_ > calibDuration_) {
                isCalibrated_ = true;
                calibrationData_.data.val = valueSum_ / valueCounter_;
                calibrationData_.data.cov = sensorToBodyMat_;
                calibrationData_.timestamp = Core::NOW();
            }

            //calibrationData_.timestamp = gyroData.timestamp;
            //calibrationData_.data.val = magBias_;
            //calibrationData_.data.cov = sensorToBodyMat_;

        }

        
        // ############################### MagnetometerCalibratorTask ###############################

        /*MagnetometerCalibratorTask::MagnetometerCalibratorTask(int64_t period, Core::Scheduler &scheduler) : Task_Periodic("MagnetometerCalibratorTask", period)
        {
            scheduler.addTask(*this);
        }

        void MagnetometerCalibratorTask::taskCheck()
        {
            if (magSubr_.isDataNew())
            {
                setPaused(false);
            }
        }

        void MagnetometerCalibratorTask::taskInit()
        {
        }

        void MagnetometerCalibratorTask::taskThread()
        {
            update();
            setPaused(true);
        }*/

    }

}