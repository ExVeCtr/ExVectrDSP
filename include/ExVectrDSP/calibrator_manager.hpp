#ifndef EXVECTRDSP_CALIBRATORMANAGER_HPP
#define EXVECTRDSP_CALIBRATORMANAGER_HPP


#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"
#include "ExVectrCore/task_types.hpp"

#include "ExVectrMath/matrix_base.hpp"

#include "ExVectrData/memory_interface.hpp"
#include "ExVectrData/memory_manager.hpp"

#include "ExVectrDSP/value_covariance.hpp"
#include "ExVectrDSP/topic_coord_transform.hpp"
#include "ExVectrDSP/calibrator_interface.hpp"


namespace VCTR
{

    namespace DSP
    {

        /**
         * * @brief This class takes care of managing calibration data for a given sensor. It will load and save calibration data to memory and update the given topic transform. If there is no calibration data in memory, it will be created and written once the calibration is complete.
         * * @tparam TYPE Type of the data. (float, double, etc.)
         * * @tparam DIM Dimension of the data. (1, 3, 6, etc.) Generally 3 for 3D sensors (IMU)
         * * @note The input of sensor calibration uses the ValueCov struct, but the value is the bias and the covariance is the scale/rotation matrix. Value = (Raw - bias) * scale/rotation. Raw is directly from sensor and value is calibrated sensor value. 
         */
        template<typename TYPE, size_t DIM = 3>
        class Calibrator_Manager : public Core::Task_Periodic
        {
        protected:

            Data::Memory_Manager memoryManager_;
            uint32_t memoryKey_ = 0;

            TopicCoordTransform<TYPE, DIM>& topicCoord_;
            Calibrator_Interface<TYPE, DIM>* calibInterface_ = nullptr;

            int64_t lastCalibTimestamp_ = 0;

            bool memoryHasCalibData_ = false;

            bool saveCalibrationData_ = false;
            bool loadCalibrationData_ = false;

            Core::Timestamped<ValueCov<TYPE, DIM>> calibData_; 

        public:

            /// @brief Constructor for the Calibrator_Manager class.
            /// @param topicCoord The topic to that uses the calibration data.
            /// @param calibInterface The calibrator interface the produced the calibration data.
            /// @param memory The memory to use to load and save calibration data.
            /// @param memoryKey The key to use to identify the calibration data in memory. MUST BE DIFFERENT FOR EACH CALIBRATOR.
            Calibrator_Manager(Calibrator_Interface<TYPE, DIM>& calibInterface, TopicCoordTransform<TYPE, DIM>& topicCoord, Data::Memory_Interface& memory, uint32_t memoryKey) : 
                Calibrator_Manager(topicCoord, memory, memoryKey)
            {
                calibInterface_ = &calibInterface;
            }

            Calibrator_Manager(TopicCoordTransform<TYPE, DIM>& topicCoord, Data::Memory_Interface& memory, uint32_t memoryKey) : 
                Core::Task_Periodic("Calibrator Manager", 500*Core::MILLISECONDS),
                memoryManager_(memory),
                memoryKey_(memoryKey),
                topicCoord_(topicCoord)
            {
                calibInterface_ = nullptr;
            }


            void taskInit() override {

                if (memoryManager_.getMemoryVersion() != memoryManager_.getManagerVersion()) {
                    setInitialised(false);
                    setRelease(Core::NOW() + 1*Core::SECONDS);
                    return;
                }

                memoryHasCalibData_ = memoryManager_.readItem(calibData_.data, memoryKey_);
                if (memoryHasCalibData_) {

                    topicCoord_.setTransform(calibData_.data.cov, calibData_.data.val);

                } else {
                    
                    calibData_.data.cov = 1;
                    calibData_.data.val = 0;
                    memoryManager_.allocateItem(calibData_.data, memoryKey_);

                }

                setPaused(true);

            }

            void saveCalibrationData() {
                saveCalibrationData_ = true;
            }

            void loadCalibrationData() {
                loadCalibrationData_ = true;
            }

            void applyCalibrationData() {
                topicCoord_.setTransform(calibData_.data.cov, calibData_.data.val);
            }

            void taskCheck() override {
                if (calibInterface_ == nullptr) return; //No calibrator interface set.
                if (lastCalibTimestamp_ != calibInterface_->getCalibrationData().timestamp) 
                    setPaused(false);
            }

            void taskThread() override {

                bool newData = false;
                if (calibInterface_ != nullptr && lastCalibTimestamp_ != calibInterface_->getCalibrationData().timestamp) {
                    newData = true;

                    calibData_ = calibInterface_->getCalibrationData();
                    //topicCoord_.setTransform(calibData_.data.cov, calibData_.data.val);

                    lastCalibTimestamp_ = calibInterface_->getCalibrationData().timestamp;

                }


                if (saveCalibrationData_ && newData) {
                    saveCalibrationData_ = false;
                    if (memoryManager_.writeItem(calibData_.data, memoryKey_)) {
                        LOG_MSG("Calibrator Manager: Calibration data saved to internal memory.\n");
                    } else {
                        LOG_MSG("Calibrator Manager: Failed to save calibration data to memory.\n");
                    }
                }
                if (loadCalibrationData_ && !saveCalibrationData_) { //No need to load the calib data again if we just saved it.
                    loadCalibrationData_ = false;
                    if (memoryManager_.readItem(calibData_.data, memoryKey_)) {
                        LOG_MSG("Calibrator Manager: Calibration data loaded from memory. Updating topicCoordTransform\n");
                        topicCoord_.setTransform(calibData_.data.cov, calibData_.data.val);
                    } else {
                        LOG_MSG("Calibrator Manager: Failed to load calibration data from memory.\n");
                    }
                }

            }




        };

    }

}











#endif