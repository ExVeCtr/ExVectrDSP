#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/print.hpp"

#include "ExVectrMath.hpp"

#include "ExVectrDSP/value_covariance.hpp"

#include "ExVectrDSP/imu_attitude_cf.hpp"

namespace VCTR
{

    namespace DSP
    {

        IMUAttitudeCF::IMUAttitudeCF()
        {
            setState(ValueCov<float, 7>(0, 1));
            // q_ = 10;
            x_(3) = 1;
            // gy = VCTR::Core::NOW();
            gyroUpdateTimestamp_ = VCTR::Core::NOW();
        }

        void IMUAttitudeCF::setGyroInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &gyroTopic)
        {
            gyroSubr_.subscribe(gyroTopic);
        }

        void IMUAttitudeCF::setAccInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &accTopic)
        {
            accSubr_.subscribe(accTopic);
        }

        void IMUAttitudeCF::setMagInput(Core::Topic<Core::Timestamped<ValueCov<float, 3>>> &magTopic)
        {
            magSubr_.subscribe(magTopic);
        }

        void IMUAttitudeCF::update()
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

            predict(VCTR::Core::NOW());

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

            if (update) { //Fix covariance and publish the new data.
                
                //Symmetrise the cov matrix for stability
                auto P = p_.block<4, 4>(3, 3);
                P = (P + P.transpose()) / 2;
                p_.block(P, 3, 3);

                Core::Timestamped<Math::Vector<float, 7>> attitudeData;
                attitudeData.timestamp = Core::NOW();
                attitudeData.data.block(x_, 3, 0, 3, 0);
                //attitudeData.data.cov.block(p_, 3, 3, 3, 3); //Copy the lower right part from p for attitude quat into cov mat.
                attitudeData.data.block(lastGyroData_.data.val - x_.block<3, 1>(0, 0), 0, 0, 0, 0, 3, 1);
                //attitudeData.data.cov.block(lastGyroData_.data.cov * 2, 0, 0, 0, 0, 3, 3);
                attitudeTopic_.publish(attitudeData);

                Core::Timestamped<Math::Vector<float, 3>> biasData;
                biasData.timestamp = attitudeData.timestamp;
                biasData.data.block(x_, 0, 0, 0, 0, 3, 1);
                //biasData.data.cov.block(p_, 0, 0, 0, 0, 3, 3); //Copy the lower right part from p for attitude quat into cov mat.
                biasTopic_.publish(biasData);

            }

        }

        void IMUAttitudeCF::predict(int64_t time)
        {
        
            auto gyro = lastGyroData_;
            gyro.timestamp = time;

            updateGyro(gyro);

        }

        void IMUAttitudeCF::updateGyro(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &gyroData)
        {

            float dTime = static_cast<float>(gyroData.timestamp - lastGyroData_.timestamp) / VCTR::Core::SECONDS;
            //lastGyroData_.timestamp = gyroData.timestamp;

            float dTimeHalf = dTime / 2;

            auto quat = x_.block<4, 1>(3, 0);
            auto bias = x_.block<3, 1>(0, 0);

            // Transform gyro from sensor to body
            auto gyro = gyroData.data.val - bias;
            auto gyroCov = gyroData.data.cov * 1000;

            // Non linear process model
            //x_[3][0] = quat[0][0] - dTimeHalf * (gyro[0][0] * quat[1][0] + gyro[1][0] * quat[2][0] + gyro[2][0] * quat[3][0]);
            //x_[4][0] = quat[1][0] + dTimeHalf * (gyro[0][0] * quat[0][0] - gyro[1][0] * quat[3][0] + gyro[2][0] * quat[2][0]);
            ///x_[5][0] = quat[2][0] + dTimeHalf * (gyro[0][0] * quat[3][0] + gyro[1][0] * quat[0][0] - gyro[2][0] * quat[1][0]);
            //x_[6][0] = quat[3][0] - dTimeHalf * (gyro[0][0] * quat[2][0] - gyro[1][0] * quat[1][0] - gyro[2][0] * quat[0][0]);

            //Simple rotation via quaternion multiplication and angle axis rotation
            auto axis = gyro.normalize();
            //axis.printTo(Core::printM);
            auto angle = gyro.magnitude() * dTime;
            //Core::printM("Angle: %.2f\n", angle * 180 / 3.14);
            auto rotQuat = Math::Quat<float>(axis, angle);
            //rotQuat.printTo(Core::printM);
            auto newQuat = rotQuat * Math::Quat_F(quat);
            x_[3][0] = newQuat[0][0];
            x_[4][0] = newQuat[1][0];
            x_[5][0] = newQuat[2][0];
            x_[6][0] = newQuat[3][0];

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

        void IMUAttitudeCF::updateAcc(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &accData)
        {

            Math::Quat<float> quat = x_.block<4, 1>(3, 0);

            // Transform acc from sensor to body and normalize to get observed gravity vector
            auto acc = (accData.data.val).normalize();
            auto accErrorIndex = (quat.conjugate().rotate(accData.data.val) - Math::GRAVITY_3F).magnitude();
            //LOG_MSG("Acc Error Index: %.2f\n", accErrorIndex);
            auto accCov = accData.data.cov * (5000 * (accErrorIndex/10 + 1)); //Increase the covariance if the accel is far off from expected gravity vector.

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
            auto gyroBias = x_.block<3, 1>(0, 0);

            auto gyroModelRot = -lastAccData_.data.val.cross(acc).normalize() * (lastAccData_.data.val.getAngleTo(Math::Vector_F(acc)))  / (float(accData.timestamp - lastAccData_.timestamp)/Core::SECONDS);
            //auto toPrint = lastAccData_.data.val * Math::Vector_F(acc);
            //auto measuredGyroBias = lastGyroData_.data.val - gyroModelRot;
            //gyroModelRot = Math::Quat_F(x).conjugate().rotate(gyroModelRot);
            //gyroModelRot(2) = 0;
            //gyroModelRot = Math::Quat_F(x).rotate(gyroModelRot);
            auto measuredGyroBias = lastGyroData_.data.val - gyroModelRot;
            gyroBias = gyroBias * 0.99995f + measuredGyroBias * 0.00005f;

            //x_.block(gyroBias, 0, 0);

            lastAccData_.data.val = acc;
            lastAccData_.data.cov = accCov;
            lastAccData_.timestamp = accData.timestamp;

        }

        void IMUAttitudeCF::updateMag(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &magData)
        {

            VCTR::Math::Quat<float> quat = x_.block<4, 1>(3, 0);
            auto quatMat = quat.operator VCTR::Math::Matrix<float, 3U, 3U>(); // Convert to matrix for rotation
            //VCTR::Math::Quat<float> quatConj = quat.conjugate();

            // Transform mag from sensor to reference frame and project onto horizontal plane while normalizing, then rotate to body frame
            auto mag = quatMat * magData.data.val;                  // Rotate to reference frame
            auto magCov = magData.data.cov * 1000; // Rotate to body frame
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

        void IMUAttitudeCF::initialiseAcc(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &accData)
        {

            Math::Quat<float> quat = x_.block<4, 1>(3, 0);

            // Transform acc from sensor to body and normalize to get observed gravity vector
            auto acc = (accData.data.val).normalize();
            auto accCov = accData.data.cov * 100;

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

        void IMUAttitudeCF::initialiseMag(const VCTR::Core::Timestamped<VCTR::DSP::ValueCov<float, 3U>> &magData)
        {

            Math::Quat<float> quat = x_.block<4, 1>(3, 0);

            // Transform acc from sensor to body and normalize to get observed gravity vector
            auto mag = (magData.data.val).normalize();
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

        void IMUAttitudeCF::setProcessNoise(const VCTR::Math::Matrix<float, 7, 7> &noise)
        {
            // q_ = noise;
        }

        const VCTR::Math::Vector<float, 7> &IMUAttitudeCF::getState()
        {
            return x_;
        }

        const VCTR::Math::Matrix<float, 7, 7> &IMUAttitudeCF::getCovariance()
        {
            return p_;
        }

        void IMUAttitudeCF::setState(const ValueCov<float, 7> &state)
        {
            x_ = state.val;
            p_ = state.cov;
        }

        Core::Topic<Core::Timestamped<Math::Vector<float, 3>>>& IMUAttitudeCF::getBiasEstTopic() 
        {
            return biasTopic_;
        }

        Core::Topic<Core::Timestamped<Math::Vector<float, 7>>>& IMUAttitudeCF::getAttitudeEstTopic() 
        {
            return attitudeTopic_;
        }


        // ############################### IMUAttitudeCFTask ###############################

        IMUAttitudeCFTask::IMUAttitudeCFTask(int64_t period, Core::Scheduler &scheduler) : Task_Periodic("IMUAttitudeCFTask", period)
        {
            scheduler.addTask(*this);
        }

        void IMUAttitudeCFTask::taskCheck()
        {
            if (gyroSubr_.isDataNew() || accSubr_.isDataNew() || magSubr_.isDataNew())
            {
                setPaused(false);
            }
        }

        void IMUAttitudeCFTask::taskInit()
        {
        }

        void IMUAttitudeCFTask::taskThread()
        {

            //if (!gyroSubr_.isDataNew())
            //{
            //    //predict(); //Predict up to this point to keep the state up to date for corrections
            //}

            update();

            setPaused(true);
        }

    }

}