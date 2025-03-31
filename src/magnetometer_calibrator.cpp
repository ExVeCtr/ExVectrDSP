#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/print.hpp"

#include "ExVectrMath.hpp"

#include "ExVectrDSP/value_covariance.hpp"

#include "ExVectrDSP/magnetometer_calibrator.hpp"

namespace VCTR
{

    namespace Data
    {

        MagnetometerCalibrator::MagnetometerCalibrator()
        {}

        void MagnetometerCalibrator::setMagInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &magTopic)
        {
            magSubr_.subscribe(magTopic);
        }

        void MagnetometerCalibrator::update()
        {

            if (magSubr_.isDataNew()) // We can only use the mag if the acc is already initialised
            {
                
                auto magVal = magSubr_.getItem().data.val;

                for (size_t i = 0; i < 3; i++)
                {
                    if (magVal(i) < magMin_(i)) magMin_(i) = magVal(i);
                    if (magVal(i) > magMax_(i)) magMax_(i) = magVal(i);
                }

                magBias_ = (magMax_ + magMin_) * 0.5;
                magTransform_(0, 0) = 1/(magMax_[0] - magMin_[0])/2;
                magTransform_(1, 1) = 1/(magMax_[1] - magMin_[1])/2;
                magTransform_(2, 2) = 1/(magMax_[2] - magMin_[2])/2;
                magTransform_ = magTransform_ * 0.05; //Scale so we reach 50uT which is assumed to be the earths magnetic field strength.

            }

        }

        
        // ############################### MagnetometerCalibratorTask ###############################

        MagnetometerCalibratorTask::MagnetometerCalibratorTask(int64_t period, Core::Scheduler &scheduler) : Task_Periodic("MagnetometerCalibratorTask", period)
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
        }

    }

}