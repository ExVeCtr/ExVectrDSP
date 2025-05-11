#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/print.hpp"

#include "ExVectrMath.hpp"

#include "ExVectrSensor/gnss.hpp"

#include "ExVectrDSP/value_covariance.hpp"

#include "ExVectrDSP/imu_gps_position_kf.hpp"

namespace VCTR
{

    namespace DSP
    {

        IMUGPSPositionKalman::IMUGPSPositionKalman()
        {
            setState(ValueCov<float, 6>(0, 1));
            // q_ = 10;
            //x_(3) = 1;
            // gy = VCTR::Core::NOW();
        }

        void IMUGPSPositionKalman::setAttitudeInput(Core::Topic<Core::Timestamped<Math::Vector<float, 7>>> &attitudeTopic)
        {
            attSubr_.subscribe(attitudeTopic);
        }

        void IMUGPSPositionKalman::setAccInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &accTopic)
        {
            accSubr_.subscribe(accTopic);
        }

        void IMUGPSPositionKalman::setBaroInput(Core::Topic<Core::Timestamped<ValueCov<float, 1>>> &baroTopic, float seaLevelPressure)
        {
            baroSubr_.subscribe(baroTopic);
            seaLevelPressure_ = seaLevelPressure;
        }

        void IMUGPSPositionKalman::setGNSSInput(Core::Topic<Core::Timestamped<SNSR::GNSSData>> &gnssTopic, Math::Vector<double, 3> refPos)
        {
            gnssSubr_.subscribe(gnssTopic);
            gnssRefPos_ = refPos;
            if (refPos.magnitude() != 0) gnssInitialised_ = true;
        }

        void IMUGPSPositionKalman::update()
        {

            bool update = false;

            // Update state.
            if (attSubr_.isDataNew())
            {
                if (!attInitialised_)
                {
                    attInitialised_ = true;
                }
                att_.block(attSubr_.getItem().data, 0, 0, 3, 0, 4, 1);
                //update = true;
            }

            if (accSubr_.size() > 0 && attInitialised_) // We can only initialise the acc if we are not moving. We can assume this to be the case if the gyro is not moving.
            {

                Core::ListBuffer<float, 5> accDataX;
                Core::ListBuffer<float, 5> accDataY;
                Core::ListBuffer<float, 5> accDataZ;

                for (size_t i = 0; i < accSubr_.size(); i++)
                {
                    accDataX.placeBack(accSubr_[i].data.val(0));
                    accDataY.placeBack(accSubr_[i].data.val(1));
                    accDataZ.placeBack(accSubr_[i].data.val(2));
                }
                
                auto accDataAvg = accSubr_[0];
                accDataAvg.data.val(0) = accDataX.getAverage();
                accDataAvg.data.val(1) = accDataY.getAverage();
                accDataAvg.data.val(2) = accDataZ.getAverage();

                if (!accInitialised_)
                {
                    accInitialised_ = true;
                    initialiseAcc(accDataAvg);
                }
                else {
                    updateAcc(accDataAvg);
                    update = true;
                }

                accSubr_.clear();

            }

            if (gnssSubr_.isDataNew() && baroInitialised_)
            {
                if (!gnssInitialised_ && gnssSubr_.getItem().data.positionValid && gnssSubr_.getItem().data.positionCov(0) < 20)
                {

                    gnssInitialised_ = true;
                    initialiseGNSS(gnssSubr_.getItem());
                }
                if (gnssInitialised_) {
                    updateGNSS(gnssSubr_.getItem());
                    update = true;
                } else {
                    //Set the position and velocity x and y to zero
                    x_(0) = 0;
                    x_(1) = 0;
                    x_(3) = 0;
                    x_(4) = 0;
                }
            }

            if (baroSubr_.isDataNew())
            {
                if (!baroInitialised_)
                {
                    baroInitialised_ = true;
                    initialiseBaro(baroSubr_.getItem());
                }
                updateBaro(baroSubr_.getItem());
                update = true;
            }

            if (zeroingMode_) {

                x_ = x_ * 0.95;

                updateAccBias(lastGNSSData_, true); //Assume zero motion.

            }

            if (update) { //publish the new data.

                //Symmetrise the cov matrix for stability
                p_ = (p_ + p_.transpose()) / 2;

                Core::Timestamped<Math::Vector<float, 6>> estimationData;
                estimationData.timestamp = Core::NOW();
                estimationData.data.block(x_);

                Core::Timestamped<Math::Matrix<float, 6, 6>> estimationCov;
                estimationCov.timestamp = Core::NOW();
                estimationCov.data.block(p_);
                //estimationData.data.cov.block(p_);
                
                stateEstTopic_.publish(estimationData);
                stateCovTopic_.publish(estimationCov);

                //Core::printM("%.3f\n", x_(5));

            }

        }

        void IMUGPSPositionKalman::predict(int64_t time)
        {
        
            auto accel = lastAccData_;
            accel.timestamp = time;

            updateAcc(accel);

        }


        void IMUGPSPositionKalman::updateAcc(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &accData)
        {

            // Transform acceleration from sensor to body to reference.
            auto sensorRotation = att_.to3x3RotMat();
            auto acc = sensorRotation * (accData.data.val - accelBias_) - Math::GRAVITY_3F;
            auto accCov = sensorRotation * accData.data.cov * sensorRotation.transpose() * 1000;

            //accMedianX_.placeBack(acc(0), true);
            //accMedianY_.placeBack(acc(1), true);
            //accMedianZ_.placeBack(acc(2), true);
            //acc(0) = accMedianX_.getMedian();
            //acc(1) = accMedianY_.getMedian();
            //acc(2) = accMedianZ_.getMedian();

            //acc(2) = 0; // We are only interested in the horizontal acceleration.

            float dt = double(accData.timestamp - lastAccData_.timestamp)/Core::SECONDS;
            float hdtq = 0.5 * dt * dt;

            auto F = VCTR::Math::Matrix<float, 6, 6>({//State transition model for prediction
                1, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0,
                dt, 0, 0, 1, 0, 0,
                0, dt, 0, 0, 1, 0, 
                0, 0, dt, 0, 0, 1
            });

            auto B = VCTR::Math::Matrix<float, 6, 3>({//Control input model.
                dt, 0, 0, 
                0, dt, 0, 
                0, 0, dt,
                hdtq, 0, 0,
                0, hdtq, 0,
                0, 0, hdtq
            });

            B = VCTR::Math::Matrix<float, 6, 3>({//Control input model.
                dt, 0, 0, 
                0, dt, 0, 
                0, 0, dt,
                hdtq, 0, 0,
                0, hdtq, 0,
                0, 0, hdtq
            });

            //accData.data.val.printTo(Core::printM);

            x_ = F * x_ + B * acc; //Prediction using accelerometer
            p_ = F * p_ * F.transpose() + B * accCov * B.transpose();


            lastAccData_ = accData;

        }

        void IMUGPSPositionKalman::updateBaro(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 1U>>& baroData)
        {

            float dt = double(baroData.timestamp - lastBaroData_.timestamp)/Core::SECONDS;
            if (dt < 0.001) return;

            float baroAlt = calcAltitudeFromPressure(baroData.data.val(0), seaLevelPressure_);
            float baroAltLast = calcAltitudeFromPressure(lastBaroData_.data.val(0), seaLevelPressure_);
            auto baroVel = (baroAlt - baroAltLast) / dt;
            
            auto baroCov = baroData.data.cov(0) * 1000;
            auto baroVelCov = baroCov * 2 / dt * 10;

            if (zeroingMode_)
                gnssRefPos_(2) = gnssRefPos_(2) * 0.9 + baroAlt * 0.1; //Update the reference position with the barometer data.

            //Core::printM("%.3f\n", baroVel);

            baroAlt = baroAlt - gnssRefPos_(2); //Relative altitude.

            LOG_MSG("Baro Alt: %.2f, Vel: %.2f\n", baroAlt, baroVel);

            auto baroState = VCTR::Math::Matrix<float, 2, 1>({
                baroVel,
                baroAlt
            });

            auto baroStateCov = VCTR::Math::Matrix<float, 2, 2>({
                baroVelCov, 0,
                0, baroCov
            });

            auto H = VCTR::Math::Matrix<float, 2, 6>({
                0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 1
            });

            auto y = baroState - H * x_; //Innovation
            auto S = H * p_ * H.transpose() + baroStateCov; //Innovation covariance
            auto K = p_ * H.transpose() * S.inverse(); //Kalman gain
            
            x_ = x_ + K * y; //Update state
            p_ = (VCTR::Math::Matrix<float, 6, 6>(1) - K * H) * p_; //Update covariance

            //Core::printM("%.3f %.3f %.3f %.3f\n", baroAlt, x_(5), baroVel, x_(2));

            lastBaroData_ = baroData;

        }

        void IMUGPSPositionKalman::updateGNSS(const VCTR::Core::Timestamped<VCTR::SNSR::GNSSData>& gnssData)
        {

            if (gnssData.data.positionCov.magnitude() > 20 || gnssData.data.velocityCov.magnitude() > 5) return;

            if (zeroingMode_) {

                auto position = gnssData.data.position;
                gnssRefPos_(0) = gnssRefPos_(0) * 0.9 + position(0) * 0.1;
                gnssRefPos_(1) = gnssRefPos_(1) * 0.9 + position(1) * 0.1;

            }

            auto relPos = calcRelPosFromGNSS(gnssData.data.position, gnssRefPos_);
            Math::Vector<float, 4> gnssState;
            gnssState(0) = gnssData.data.velocity(0);
            gnssState(1) = gnssData.data.velocity(1);
            gnssState(2) = relPos(0);
            gnssState(3) = relPos(1);

            auto gnssCov = VCTR::Math::Matrix<float, 4, 4>();
            gnssCov(0, 0) = gnssData.data.velocityCov(0) * 200;
            gnssCov(1, 1) = gnssData.data.velocityCov(1) * 200;
            gnssCov(2, 2) = gnssData.data.positionCov(0) * 200;
            gnssCov(3, 3) = gnssData.data.positionCov(1) * 200;

            auto H = VCTR::Math::Matrix<float, 4, 6>({
                1, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 1, 0
            });

            auto y = gnssState - H * x_; //Innovation
            auto S = H * p_ * H.transpose() + gnssCov; //Innovation covariance
            auto K = p_ * H.transpose() * S.inverse(); //Kalman gain

            x_ = x_ + K * y; //Update state
            //x_(1) = gnssState(1);
            //x_(4) = gnssState(3);
            p_ = (VCTR::Math::Matrix<float, 6, 6>(1) - K * H) * p_; //Update covariance

            //updateAccBias(gnssData, false); //Update the accelerometer bias estimate using the GNSS data.

            lastGNSSData_ = gnssData;

            //Update the sea level pressure estimate.
            //float factor = 0.0001;
            //if (zeroingMode_) factor = 0.1;
            //auto sealevelPres = calcSealevelPressFromAltitude(lastBaroData_.data.val(0), gnssData.data.position(2));
            //seaLevelPressure_ = seaLevelPressure_ * (1.0f - factor) + sealevelPres * factor;


        }

        void IMUGPSPositionKalman::updateAccBias(const VCTR::Core::Timestamped<VCTR::SNSR::GNSSData>& gnssData, float assumeZero) {

            float dTime = double(gnssData.timestamp - lastGNSSData_.timestamp) / Core::SECONDS;

            if (!assumeZero && dTime < 0.001) return;

            //Calulate the acceleration in reference frame using gnss data
            auto accelGNSSWorld = (gnssData.data.velocity - lastGNSSData_.data.velocity) / dTime + Math::GRAVITY_3F; //Acceleration in world frame.
            if (assumeZero)
                accelGNSSWorld = Math::GRAVITY_3F;
            auto accelGNSSBody = att_.rotate(accelGNSSWorld); //Transform to body frame 

            float factor = 0.1;
            //if (zeroingMode_) factor = 0.1;
            accelBias_ = accelBias_ * (1.0f - factor) - (accelGNSSBody - lastAccData_.data.val) * factor; //Update the bias estimate

            //LOG_MSG("Acc Bias: %.3f %.3f %.3f\n", accelBias_(0), accelBias_(1), accelBias_(2));

        }

        void IMUGPSPositionKalman::initialiseAcc(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &accData)
        {

            //Math::Quat<float> quat = x_.block<4, 1>(3, 0);

            // Transform acceleration from sensor to body to reference.
            //auto acc = quat.rotate(accRot_ * accData.data.val);
            //auto accCov = accRot_ * accData.data.cov * accRot_.transpose() * 1000;

            lastAccData_ = accData;

        }

        void IMUGPSPositionKalman::initialiseBaro(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 1U>>& baroData)
        {
            auto altitude = calcAltitudeFromPressure(baroData.data.val(0), seaLevelPressure_);
            gnssRefPos_(2) = altitude;
            x_(5) = altitude - gnssRefPos_(2);
            lastBaroData_ = baroData;
        }

        void IMUGPSPositionKalman::initialiseGNSS(const VCTR::Core::Timestamped<VCTR::SNSR::GNSSData>& gnssData)
        {
            gnssRefPos_ = gnssData.data.position;
            seaLevelPressure_ = calcSealevelPressFromAltitude(lastBaroData_.data.val(0), x_(5) + gnssRefPos_(2));
        }

        void IMUGPSPositionKalman::setPositionReference() {

            if (gnssInitialised_ && baroInitialised_) {
                seaLevelPressure_ = calcSealevelPressFromAltitude(lastBaroData_.data.val(0), lastGNSSData_.data.position(2));
            }

            if (gnssInitialised_) {
                gnssRefPos_(0) = lastGNSSData_.data.position(0);
                gnssRefPos_(1) = lastGNSSData_.data.position(1);
            } else {
                gnssRefPos_(0) = 0;
                gnssRefPos_(1) = 0;
            }

            if (baroInitialised_) {
                gnssRefPos_(2) = calcAltitudeFromPressure(lastBaroData_.data.val(0), seaLevelPressure_);
            } else {
                gnssRefPos_(2) = 0;
            }

            x_ = 0;

        }

        void IMUGPSPositionKalman::setProcessNoise(const VCTR::Math::Matrix<float, 6, 6> &noise)
        {
            // q_ = noise;
        }

        const VCTR::Math::Vector<float, 6> &IMUGPSPositionKalman::getState()
        {
            return x_;
        }

        const VCTR::Math::Matrix<float, 6, 6> &IMUGPSPositionKalman::getCovariance()
        {
            return p_;
        }

        void IMUGPSPositionKalman::setState(const ValueCov<float, 6> &state)
        {
            x_ = state.val;
            p_ = state.cov;
        }

        Core::Topic<Core::Timestamped<Math::Vector<float, 6>>>& IMUGPSPositionKalman::getStateEstTopic() 
        {
            return stateEstTopic_;
        }

        Core::Topic<Core::Timestamped<Math::Matrix<float, 6, 6>>>& IMUGPSPositionKalman::getStateCovTopic()
        {
            return stateCovTopic_;
        }

        float IMUGPSPositionKalman::calcAltitudeFromPressure(float pressure, float seaLevelPressure)
        {
            return 44330 * (1 - pow(pressure / seaLevelPressure, 0.1903));
        }

        float IMUGPSPositionKalman::calcSealevelPressFromAltitude(float pressure, float altitude)
        {
            return pressure / pow(1 - altitude / 44330, 5.255);
        }

        Math::Vector_F IMUGPSPositionKalman::calcRelPosFromGNSS(const Math::Vector<double, 3> &gnssData, const Math::Vector<double, 3> &refPos) 
        {
            
            //Using data from the WGS84 ellipsoid.
            const double a = 6371000.0; // Semi major axis
            const double b = 6356752.3142; // Semi minor axis

            float north = (gnssData(0) - refPos(0)) * b;
            float west = -(gnssData(1) - refPos(1)) * a * cos(gnssData(0));
            float up = gnssData(2) - refPos(2);

            return Math::Vector_F({north, west, up});

        }


        // ############################### IMUGPSPositionKalmanTask ###############################

        IMUGPSPositionKalmanTask::IMUGPSPositionKalmanTask(int64_t period, Core::Scheduler &scheduler) : Task_Periodic("IMUGPSPositionKalmanTask", period)
        {
            scheduler.addTask(*this);
        }

        void IMUGPSPositionKalmanTask::taskCheck()
        {
            if (attSubr_.isDataNew() || accSubr_.size() > 0)
            {
                setPaused(false);
            }
        }

        void IMUGPSPositionKalmanTask::taskInit()
        {
        }

        void IMUGPSPositionKalmanTask::taskThread()
        {

            predict();

            update();

            setPaused(true);
        }

    }

}