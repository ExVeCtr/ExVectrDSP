#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/print.hpp"

#include "ExVectrMath.hpp"

#include "ExVectrDSP/value_covariance.hpp"

#include "ExVectrDSP/imu_attitude_ekf.hpp"

namespace VCTR
{

    namespace DSP
    {

        IMUAttitudeEKF::IMUAttitudeEKF()
        {
            setState(ValueCov<float, 7>(0, 1));
            // q_ = 10;
            x_(3) = 1;
            // gy = VCTR::Core::NOW();
            gyroUpdateTimestamp_ = VCTR::Core::NOW();
        }

        void IMUAttitudeEKF::setGyroInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &gyroTopic, const VCTR::Math::Matrix<float, 3, 3> &gyroRot)
        {
            gyroSubr_.subscribe(gyroTopic);
            gyroRot_ = gyroRot;
        }

        void IMUAttitudeEKF::setAccInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &accTopic, const VCTR::Math::Matrix<float, 3, 3> &accRot)
        {
            accSubr_.subscribe(accTopic);
            accRot_ = accRot;
        }

        void IMUAttitudeEKF::setMagInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &magTopic, const VCTR::Math::Matrix<float, 3, 3> &magTransform, const VCTR::Math::Matrix<float, 3, 1> &magBias)
        {
            magSubr_.subscribe(magTopic);
            magTransform_ = magTransform;
            magBias_ = magBias;
        }

        void IMUAttitudeEKF::update()
        {

            bool update = false;

            // Update state.
            if (gyroSubr_.isDataNew())
            {
                if (!gyroInitialised_)
                {
                    gyroInitialised_ = true;
                }
                updateGyro(gyroSubr_.getItem());
                update = true;
            }

            if (accSubr_.isDataNew() && (accInitialised_ || (lastGyroData_.data.val.magnitude() < 0.1 && gyroInitialised_))) // We can only initialise the acc if we are not moving. We can assume this to be the case if the gyro is not moving.
            {
                if (!accInitialised_)
                {
                    accInitialised_ = true;
                    initialiseAcc(accSubr_.getItem());
                }
                else {
                    updateAcc(accSubr_.getItem());
                    update = true;
                }
            }

            if (magSubr_.isDataNew() && accInitialised_) // We can only use the mag if the acc is already initialised
            {
                if (!magInitialised_)
                {
                    magInitialised_ = true;
                    initialiseMag(magSubr_.getItem());
                }
                else {
                    updateMag(magSubr_.getItem());
                    update = true;
                }
            }

            if (update) { //Fix covarience and publish the new data.
                
                //Symmetrise the cov matrix for stability
                auto P = p_.block<4, 4>(3, 3);
                P = (P + P.transpose()) / 2;
                p_.block(P, 3, 3);

                Core::Timestamped<ValueCov<float, 7>> attitudeData;
                attitudeData.timestamp = Core::NOW();
                attitudeData.data.val.block(x_, 3, 0, 3, 0);
                attitudeData.data.cov.block(p_, 3, 3, 3, 3); //Copy the lower right part from p for attitude quat into cov mat.
                attitudeData.data.val.block(lastGyroData_.data.val - x_.block<3, 1>(0, 0), 0, 0, 0, 0, 3, 1);
                attitudeData.data.cov.block(lastGyroData_.data.cov * 2, 0, 0, 0, 0, 3, 3);
                attitudeTopic_.publish(attitudeData);

                Core::Timestamped<ValueCov<float, 3>> biasData;
                biasData.timestamp = attitudeData.timestamp;
                biasData.data.val.block(x_, 0, 0, 0, 0, 3, 1);
                biasData.data.cov.block(p_, 0, 0, 0, 0, 3, 3); //Copy the lower right part from p for attitude quat into cov mat.
                biasTopic_.publish(biasData);

            }

        }

        void IMUAttitudeEKF::predict(int64_t time)
        {
        
            auto gyro = lastGyroData_;
            gyro.timestamp = time;

            updateGyro(gyro);

        }

        void IMUAttitudeEKF::updateGyro(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &gyroData)
        {

            float dTime = static_cast<float>(gyroData.timestamp - gyroUpdateTimestamp_) / VCTR::Core::SECONDS;
            gyroUpdateTimestamp_ = gyroData.timestamp;

            float dTimeHalf = dTime / 2;

            auto quat = x_.block<4, 1>(3, 0);
            auto bias = x_.block<3, 1>(0, 0);

            // Transform gyro from sensor to body
            auto gyro = gyroRot_ * (gyroData.data.val - bias);
            auto gyroCov = gyroRot_ * gyroData.data.cov * gyroRot_.transpose() * 1000;

            // Non linear process model
            x_[3][0] = quat[0][0] - dTimeHalf * (gyro[0][0] * quat[1][0] + gyro[1][0] * quat[2][0] + gyro[2][0] * quat[3][0]);
            x_[4][0] = quat[1][0] + dTimeHalf * (gyro[0][0] * quat[0][0] - gyro[1][0] * quat[3][0] + gyro[2][0] * quat[2][0]);
            x_[5][0] = quat[2][0] + dTimeHalf * (gyro[0][0] * quat[3][0] + gyro[1][0] * quat[0][0] - gyro[2][0] * quat[1][0]);
            x_[6][0] = quat[3][0] - dTimeHalf * (gyro[0][0] * quat[2][0] - gyro[1][0] * quat[1][0] - gyro[2][0] * quat[0][0]);

            float norm = sqrtf(x_[3][0] * x_[3][0] + x_[4][0] * x_[4][0] + x_[5][0] * x_[5][0] + x_[6][0] * x_[6][0]);
            if (x_[0][0] < 0) // Normalize so the w component is always positive
            {
                norm = -norm;
            }
            x_[3][0] /= norm;
            x_[4][0] /= norm;
            x_[5][0] /= norm;
            x_[6][0] /= norm;

            // Jacobian of process model
            auto F = VCTR::Math::Matrix<float, 4, 4>::eye();
            F[0][1] = -dTimeHalf * gyro(0);
            F[0][2] = -dTimeHalf * gyro(1);
            F[0][3] = -dTimeHalf * gyro(2);

            F[1][0] = dTimeHalf * gyro(0);
            F[1][2] = dTimeHalf * gyro(2);
            F[1][3] = -dTimeHalf * gyro(1);

            F[2][0] = dTimeHalf * gyro(1);
            F[2][1] = -dTimeHalf * gyro(2);
            F[2][3] = dTimeHalf * gyro(0);

            F[3][0] = dTimeHalf * gyro(2);
            F[3][1] = dTimeHalf * gyro(1);
            F[3][2] = -dTimeHalf * gyro(0);

            // Covariance of process model
            auto W = VCTR::Math::Matrix<float, 4, 3>::eye();
            W[0][0] = -dTimeHalf * quat[1][0];
            W[0][1] = -dTimeHalf * quat[2][0];
            W[0][2] = -dTimeHalf * quat[3][0];

            W[1][0] = dTimeHalf * quat[0][0];
            W[1][1] = -dTimeHalf * quat[3][0];
            W[1][2] = dTimeHalf * quat[2][0];

            W[2][0] = dTimeHalf * quat[3][0];
            W[2][1] = dTimeHalf * quat[0][0];
            W[2][2] = -dTimeHalf * quat[1][0];

            W[3][0] = -dTimeHalf * quat[2][0];
            W[3][1] = dTimeHalf * quat[1][0];
            W[3][2] = dTimeHalf * quat[0][0];

            auto P = F * p_.block<4, 4>(3, 3) * F.transpose() + W * gyroCov * W.transpose();
            p_.block(P, 3, 3);

            lastGyroData_ = gyroData;

        }

        void IMUAttitudeEKF::updateAcc(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &accData)
        {

            Math::Quat<float> quat = x_.block<4, 1>(3, 0);

            // Transform acc from sensor to body and normalize to get observed gravity vector
            auto acc = (accRot_ * accData.data.val).normalize();
            auto accCov = accRot_ * accData.data.cov * accRot_.transpose() * 1000;

            // Gravity reference vector
            // auto g = VCTR::Math::Matrix<float, 3, 1>({0, 0, 1});

            // Gravity measurement model
            auto h = Math::Matrix<float, 3, 1>({2 * (quat[1][0] * quat[3][0] - quat[0][0] * quat[2][0]),
                                                2 * (quat[0][0] * quat[1][0] + quat[2][0] * quat[3][0]),
                                                2 * (0.5f - quat[1][0] * quat[1][0] - quat[2][0] * quat[2][0])});

            // Jacobian of measurement model
            auto H = VCTR::Math::Matrix<float, 3, 4>({-2 * quat[2][0], 2 * quat[3][0], -2 * quat[0][0], 2 * quat[1][0],
                                                      2 * quat[1][0], 2 * quat[0][0], 2 * quat[3][0], 2 * quat[2][0],
                                                      0, -4 * quat[1][0], -4 * quat[2][0], 0});

            // Residual
            auto y = acc - h;

            // Innovation covariance
            auto S = H * p_.block<4, 4>(3, 3) * H.transpose() + accCov;

            // Kalman gain
            auto K = p_.block<4, 4>(3, 3) * H.transpose() * S.inverse();

            // Update state and covariance
            auto x = quat + K * y;
            auto P = (VCTR::Math::Matrix<float, 4, 4>::eye() - K * H) * p_.block<4, 4>(3, 3);

            // Normalize quaternion and update matricies
            float norm = sqrtf(x[0][0] * x[0][0] + x[1][0] * x[1][0] + x[2][0] * x[2][0] + x[3][0] * x[3][0]);
            if (x[0][0] < 0) // Normalize so the w component is always positive
            {
                norm = -norm;
            }
            x_[3][0] = x[0][0] / norm;
            x_[4][0] = x[1][0] / norm;
            x_[5][0] = x[2][0] / norm;
            x_[6][0] = x[3][0] / norm;
            
            p_.block(P, 3, 3);

            // Update gyroBias using a simple low pass filter
            //auto gyroBias = x_.block<3, 1>(0, 0);

            //auto gyroModelRot = -lastAccData_.data.val.cross(acc).normalize() * (lastAccData_.data.val.getAngleTo(Math::Vector_F(acc)))  / (float(accData.timestamp - lastAccData_.timestamp)/Core::SECONDS);
            //auto toPrint = lastAccData_.data.val * Math::Vector_F(acc);
            //auto measuredGyroBias = lastGyroData_.data.val - gyroModelRot;
            //measuredGyroBias = Math::Quat_F(x).conjugate().rotate(measuredGyroBias);
            //measuredGyroBias(2) = 0;
            //measuredGyroBias = Math::Quat_F(x).rotate(measuredGyroBias);
            //gyroBias = gyroBias * 0.9999f + measuredGyroBias * 0.0001f;

            //x_.block(gyroBias, 0, 0);

            lastAccData_.data.val = acc;
            lastAccData_.data.cov = accCov;
            lastAccData_.timestamp = accData.timestamp;

        }

        void IMUAttitudeEKF::updateMag(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &magData)
        {

            VCTR::Math::Quat<float> quat = x_.block<4, 1>(3, 0);
            auto quatMat = quat.operator VCTR::Math::Matrix<float, 3U, 3U>(); // Convert to matrix for rotation
            //VCTR::Math::Quat<float> quatConj = quat.conjugate();

            // Transform mag from sensor to reference frame and project onto horizontal plane while normalizing, then rotate to body frame
            auto mag = quatMat * magTransform_ * (magData.data.val - magBias_);                  // Rotate to reference frame
            auto magCov = magTransform_ * magData.data.cov * magTransform_.transpose() * 1000; // Rotate to body frame
            // magCov = quatMat * magCov * quatMat.transpose(); //Rotate to reference frame

            //magCov.printTo(Core::printM);

            //Core::printM("Mag: %f \n", mag.magnitude());

            // Check if mag measurement is within normal earth bounds. Leave if not
            if (mag.magnitude() > 0.5 || mag.magnitude() < 0.01)
            {
                return;
            }

            mag[2][0] = 0;
            magCov[2][2] = 0;
            mag = mag.normalize();  // Normalize to get observed magnetic field vector
            //mag = quat.rotate(mag); // Rotate to body frame




            /// ############ below is a simpler way to fuse the mag data using complementary type like filter ############

            //auto magRef = VCTR::Math::Matrix<float, 3, 1>({1, 0, 0}); // Reference magnetic field vector

            // Calculate the rotation vector between the reference and observed magnetic field vectors
            auto rotAngle = atan2f(mag[1][0], mag[0][0]);

            // Convert rotation vector to quaternion but only a little bit
            auto rotQuat = VCTR::Math::Quat<float>({0, 0, 1}, -rotAngle * 0.005f);

            // Fuse the rotation quaternion with the current state quaternion
            auto newQuat = quat * rotQuat;

            // Update state matricies
            x_.block(newQuat, 3, 0);


            /// ############ below is the standard EKF but not currently use as it has issues ############
            
            /*
            // Magnetic Measurement model
            auto h = VCTR::Math::Matrix<float, 3, 1>({2 * (0.5f - quat[2][0] * quat[2][0] - quat[3][0] * quat[3][0]),
                                                      2 * (quat[1][0] * quat[2][0] - quat[0][0] * quat[3][0]),
                                                      2 * (quat[0][0] * quat[2][0] + quat[1][0] * quat[3][0])});

            // Jacobian of measurement model
            auto H = VCTR::Math::Matrix<float, 3, 4>({0, 0, -4 * quat[2][0], -4 * quat[3][0],
                                                      -2 * quat[3][0], 2 * quat[2][0], 2 * quat[1][0], -2 * quat[0][0],
                                                      2 * quat[2][0], 2 * quat[3][0], 2 * quat[0][0], 2 * quat[1][0]});

            // Residual
            auto y = mag - h;

            // Innovation covariance
            auto S = H * p_.block<4, 4>(3, 3) * H.transpose() + magCov;

            // Kalman gain
            auto K = p_.block<4, 4>(3, 3) * H.transpose() * S.inverse();

            // Update state and covariance
            auto x = quat + K * y;
            auto P = (VCTR::Math::Matrix<float, 4, 4>::eye() - K * H) * p_.block<4, 4>(3, 3);

            // Normalize quaternion and update matricies
            float norm = sqrt(x[0][0] * x[0][0] + x[1][0] * x[1][0] + x[2][0] * x[2][0] + x[3][0] * x[3][0]);
            if (x[0][0] < 0) // Normalize so the w component is always positive
            {
                norm = -norm;
            }
            x_[3][0] = x[0][0] / norm;
            x_[4][0] = x[1][0] / norm;
            x_[5][0] = x[2][0] / norm;
            x_[6][0] = x[3][0] / norm;

            p_.block(P, 3, 3);

            // Update gyroBias using a simple low pass filter
            auto gyroBias = x_.block<3, 1>(0, 0);
            auto gyroModelRot = lastMagData_.data.val.cross(mag) / (float(magData.timestamp - lastMagData_.timestamp)/Core::SECONDS);
            gyroBias = gyroBias * 0.999 + (lastGyroData_.data.val - gyroModelRot) * 0.001;

            x_.block(gyroBias, 0, 0);

            lastMagData_.data.val = mag;
            lastMagData_.data.cov = magCov;
            lastMagData_.timestamp = magData.timestamp;*/

        }

        void IMUAttitudeEKF::initialiseAcc(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &accData)
        {

            Math::Quat<float> quat = x_.block<4, 1>(3, 0);

            // Transform acc from sensor to body and normalize to get observed gravity vector
            auto acc = (accRot_ * accData.data.val).normalize();
            auto accCov = accRot_ * accData.data.cov * accRot_.transpose() * 100;

            // Gravity reference
            auto g = Math::Matrix<float, 3, 1>({0, 0, 1});

            // Gravity measurement in reference frame
            auto h = quat.conjugate().rotate(g);

            // Get the rotation vector from measured gravity to reference gravity
            auto rotVec = h.cross(acc);
            auto rotAngle = h.getAngleTo(Math::Vector<float, 3>(acc));

            // Rotate the current quaternion by the rotation vector
            auto q = quat * Math::Quat<float>(rotVec.normalize(), rotAngle);

            // Update state
            x_.block(q, 3, 0);

            lastAccData_.data.val = acc;
            lastAccData_.data.cov = accCov;
            lastAccData_.timestamp = accData.timestamp;

        }

        void IMUAttitudeEKF::initialiseMag(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &magData)
        {

            Math::Quat<float> quat = x_.block<4, 1>(3, 0);

            // Transform acc from sensor to body and normalize to get observed gravity vector
            auto mag = (magTransform_ * magData.data.val).normalize();
            //auto accCov = accRot_ * accData.data.cov * accRot_.transpose() * 100;

            //Project mag to horizontal plane
            mag = quat.conjugate().rotate(mag); //Rotate to reference frame
            mag(2) = 0;
            mag = quat.rotate(mag.normalize()); //Normalize and rotate back to body frame

            // Rotation angle in horizontal plane
            auto rotAngle = atan2f(mag[1][0], mag[0][0]);

            // Convert rotation vector to quaternion but only a little bit
            auto rotQuat = VCTR::Math::Quat<float>({0, 0, 1}, -rotAngle);

            // Fuse the rotation quaternion with the current state quaternion
            auto newQuat = quat * rotQuat;

            // Update state matricies
            x_.block(newQuat, 3, 0);

            lastMagData_.data.val = mag;
            //lastMagData_.data.cov = magCov;
            lastMagData_.timestamp = magData.timestamp;

        }

        void IMUAttitudeEKF::setProcessNoise(const VCTR::Math::Matrix<float, 7, 7> &noise)
        {
            // q_ = noise;
        }

        const VCTR::Math::Vector<float, 7> &IMUAttitudeEKF::getState()
        {
            return x_;
        }

        const VCTR::Math::Matrix<float, 7, 7> &IMUAttitudeEKF::getCovariance()
        {
            return p_;
        }

        void IMUAttitudeEKF::setState(const ValueCov<float, 7> &state)
        {
            x_ = state.val;
            p_ = state.cov;
        }

        Core::Topic<Core::Timestamped<ValueCov<float, 3>>>& IMUAttitudeEKF::getBiasEstTopic() 
        {
            return biasTopic_;
        }

        Core::Topic<Core::Timestamped<ValueCov<float, 7>>>& IMUAttitudeEKF::getAttitudeEstTopic() 
        {
            return attitudeTopic_;
        }


        // ############################### IMUAttitudeEKFTask ###############################

        IMUAttitudeEKFTask::IMUAttitudeEKFTask(int64_t period, Core::Scheduler &scheduler) : Task_Periodic("IMUAttitudeEKFTask", period)
        {
            scheduler.addTask(*this);
        }

        void IMUAttitudeEKFTask::taskCheck()
        {
            if (gyroSubr_.isDataNew() || accSubr_.isDataNew() || magSubr_.isDataNew())
            {
                setPaused(false);
            }
        }

        void IMUAttitudeEKFTask::taskInit()
        {
        }

        void IMUAttitudeEKFTask::taskThread()
        {

            predict();

            update();

            setPaused(true);
        }

    }

}