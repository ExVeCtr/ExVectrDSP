#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/print.hpp"

#include "ExVectrMath.hpp"

#include "ExVectrDSP/value_covariance.hpp"

#include "ExVectrDSP/imu_gps_position_kf.hpp"

namespace VCTR
{

    namespace Data
    {

        IMUGPSPositionKalman::IMUGPSPositionKalman()
        {
            setState(ValueCov<float, 6>(0, 1));
            // q_ = 10;
            //x_(3) = 1;
            // gy = VCTR::Core::NOW();
        }

        void IMUGPSPositionKalman::setAttitudeInput(Core::Topic<Core::Timestamped<ValueCov<float, 7>>> &attitudeTopic)
        {
            attSubr_.subscribe(attitudeTopic);
        }

        void IMUGPSPositionKalman::setAccInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &accTopic, const VCTR::Math::Matrix<float, 3, 3> &accRot)
        {
            accSubr_.subscribe(accTopic);
            accRot_ = accRot;
        }

        void IMUGPSPositionKalman::setBaroInput(Core::Topic<Core::Timestamped<ValueCov<float, 1>>> &baroTopic, float seaLevelPressure)
        {
            baroSubr_.subscribe(baroTopic);
            seaLevelPressure_ = seaLevelPressure;
        }

        void IMUGPSPositionKalman::setGNSSInput(Core::Topic<Core::Timestamped<Data::GNSSData>> &gnssTopic, Math::Vector<double, 3> refPos)
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
                att_.block(attSubr_.getItem().data.val, 0, 0, 3, 0, 4, 1);
                //update = true;
            }

            if (accSubr_.size() > 0 && attInitialised_) // We can only initialise the acc if we are not moving. We can assume this to be the case if the gyro is not moving.
            {

                Core::ListBuffer<Math::Vector_F, 5> accData;

                for (size_t i = 0; i < accSubr_.size(); i++)
                {
                    accData.placeBack(accSubr_[i].data.val);
                }
                
                auto accDataAvg = accSubr_[0];
                accDataAvg.data.val = accData.getAverage();

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

            if (update) { //publish the new data.

                Core::Timestamped<ValueCov<float, 6>> estimationData;
                estimationData.timestamp = Core::NOW();
                estimationData.data.val.block(x_);
                estimationData.data.cov.block(p_);
                
                stateEstTopic_.publish(estimationData);

                //Core::printM("%.3f\n", x_(5));

            }

        }

        void IMUGPSPositionKalman::predict(int64_t time)
        {
        
            auto accel = lastAccData_;
            accel.timestamp = time;

            updateAcc(accel);

        }


        void IMUGPSPositionKalman::updateAcc(const VCTR::Core::Timestamped<VCTR::Data::ValueCov<float, 3U>> &accData)
        {

            // Transform acceleration from sensor to body to reference.
            auto sensorRotation = att_.to3x3RotMat() * accRot_;
            auto acc = sensorRotation * accData.data.val;
            auto accCov = sensorRotation * accData.data.cov * sensorRotation.transpose() * 1000;

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

            //accData.data.val.printTo(Core::printM);

            x_ = F * x_ + B * (acc - Math::GRAVITY_3F); //Prediction using accelerometer
            p_ = F * p_ * F.transpose() + B * accCov * B.transpose();


            lastAccData_ = accData;

        }

        void IMUGPSPositionKalman::updateBaro(const VCTR::Core::Timestamped<VCTR::Data::ValueCov<float, 1U>>& baroData)
        {

            float dt = double(baroData.timestamp - lastBaroData_.timestamp)/Core::SECONDS;
            if (dt < 0.001) return;

            float baroAlt = calcAltitudeFromPressure(baroData.data.val(0), seaLevelPressure_) - gnssRefPos_(2);
            float baroAltLast = calcAltitudeFromPressure(lastBaroData_.data.val(0), seaLevelPressure_) - gnssRefPos_(2);
            auto baroVel = (baroAlt - baroAltLast) / dt;
            auto baroVelCov = baroData.data.cov(0) * 2 / dt;

            //Core::printM("%.3f\n", baroVel);

            auto baroState = VCTR::Math::Matrix<float, 2, 1>({
                baroVel,
                baroAlt
            });

            auto baroStateCov = VCTR::Math::Matrix<float, 2, 2>({
                baroVelCov, 0,
                0, baroVelCov
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

            //Core::printM("%.3f\n", x_(5));

            lastBaroData_ = baroData;

        }

        void IMUGPSPositionKalman::updateGNSS(const VCTR::Core::Timestamped<VCTR::Data::GNSSData>& gnssData)
        {

            //if (gnssData.data.positionCov.magnitude() > 20 || gnssData.data.velocityCov.magnitude() > 5) return;

            auto relPos = calcRelPosFromGNSS(gnssData.data.position, gnssRefPos_);
            Math::Vector<float, 4> gnssState;
            gnssState(0) = gnssData.data.velocity(0);
            gnssState(1) = gnssData.data.velocity(1);
            gnssState(2) = relPos(0);
            gnssState(3) = relPos(1);

            auto gnssCov = VCTR::Math::Matrix<float, 4, 4>();
            gnssCov(0, 0) = gnssData.data.velocityCov(0);
            gnssCov(1, 1) = gnssData.data.velocityCov(1);
            gnssCov(2, 2) = gnssData.data.positionCov(0);
            gnssCov(3, 3) = gnssData.data.positionCov(1);
            gnssCov = gnssCov * 1000;

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
            p_ = (VCTR::Math::Matrix<float, 6, 6>(1) - K * H) * p_; //Update covariance

            //Update the sea level pressure estimate.
            auto sealevelPres = calcSealevelPressFromAltitude(lastBaroData_.data.val(0), gnssData.data.position(2));
            seaLevelPressure_ = seaLevelPressure_ * 0.9 + sealevelPres * 0.1;

        }

        void IMUGPSPositionKalman::initialiseAcc(const VCTR::Core::Timestamped<VCTR::Data::ValueCov<float, 3U>> &accData)
        {

            //Math::Quat<float> quat = x_.block<4, 1>(3, 0);

            // Transform acceleration from sensor to body to reference.
            //auto acc = quat.rotate(accRot_ * accData.data.val);
            //auto accCov = accRot_ * accData.data.cov * accRot_.transpose() * 1000;

            lastAccData_ = accData;

        }

        void IMUGPSPositionKalman::initialiseBaro(const VCTR::Core::Timestamped<VCTR::Data::ValueCov<float, 1U>>& baroData)
        {
            x_(5) = calcAltitudeFromPressure(baroData.data.val(0), seaLevelPressure_) - gnssRefPos_(2);
            lastBaroData_ = baroData;
        }

        void IMUGPSPositionKalman::initialiseGNSS(const VCTR::Core::Timestamped<VCTR::Data::GNSSData>& gnssData)
        {
            gnssRefPos_ = gnssData.data.position;
            seaLevelPressure_ = calcSealevelPressFromAltitude(lastBaroData_.data.val(0), x_(5) + gnssRefPos_(2));
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

        Core::Topic<Core::Timestamped<ValueCov<float, 6>>>& IMUGPSPositionKalman::getStateEstTopic() 
        {
            return stateEstTopic_;
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