#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/print.hpp"

#include "ExVectrMath.hpp"

#include "ExVectrDSP/value_covariance.hpp"

#include "ExVectrDSP/calibrator_magnetometer.hpp"


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

        void Calibrator_Magnetometer::initCalibration() {
            magMax_ = -1000;
            magMin_ = 1000;
            magBias_ = 0;
            magScale_ = Math::Matrix<float, 3, 3>::eye();
            isCalibrated_ = false;
        }

        void Calibrator_Magnetometer::setMagInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &magTopic, const Math::Matrix<float, 3, 3>& magRotation)
        {
            magSubr_.subscribe(magTopic);
            magSubr_.setCallback(this, &Calibrator_Magnetometer::update);

            sensorToBodyMat_ = magRotation;
        }

        void Calibrator_Magnetometer::update(Core::Timestamped<ValueCov<float, 3>> const &magData)
        {
                
            auto magVal = magData.data.val;

            for (size_t i = 0; i < 3; i++)
            {
                if (magVal(i) < magMin_(i)) magMin_(i) = magVal(i);
                if (magVal(i) > magMax_(i)) magMax_(i) = magVal(i);
            }

            magBias_ = (magMax_ + magMin_) * 0.5;
            magScale_(0, 0) = 0.05/(magMax_(0) - magMin_(0))/2;
            magScale_(1, 1) = 0.05/(magMax_(1) - magMin_(1))/2;
            magScale_(2, 2) = 0.05/(magMax_(2) - magMin_(2))/2;
            //magScale_ = magScale_ * 0.05; //Scale so we reach 50uT which is assumed to be the earths magnetic field strength.

            calibrationData_.timestamp = magData.timestamp;
            calibrationData_.data.val = magBias_;
            calibrationData_.data.cov = sensorToBodyMat_ * magScale_;

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