#ifndef ExVectrDSP_KALMANFILTER_H
#define ExVectrDSP_KALMANFILTER_H

#include "ExVectrMath/matrix_base.hpp"
#include "ExVectrMath/matrix_vector.hpp"

#include "value_covariance.hpp"

namespace VCTR
{

    namespace DSP
    {

        /**
         * A class implementing a kalman filter.
         * @tparam LEN Length of the state vector. E.g. for simple position in xyz this is 3, with position, velocity and accel is 9.
         * @tparam TYPE The type of values used. Usually float and defaults to this.
        */
        template<size_t LEN, typename TYPE = float>
        class KalmanFilter {    
        private:

            /// @brief Current state estimation.
            VCTR::Math::Vector<TYPE, LEN> x_;
            /// @brief Current state estimation covariance.
            VCTR::Math::Matrix<TYPE, LEN, LEN> p_;

            /// @brief Matrix modelling the system transition (aka prediction).
            VCTR::Math::Matrix<TYPE, LEN, LEN> f_;
            /// @brief Matrix modelling the system prediction covariance aka process noise. How much we trust the prediction.
            VCTR::Math::Matrix<TYPE, LEN, LEN> q_;

        public:

            /**
             * @param transition Transition matrix. Transforms state into new state prediction.
             * @param noise Noise of the state prediction. 
            */
            KalmanFilter(const VCTR::Math::Matrix<TYPE, LEN, LEN>& transition = 1, const VCTR::Math::Matrix<TYPE, LEN, LEN>& noise = 1);

            /**
             * Predicts the system state using the given input and input model.
             * @param inputModel Models the input to the system. E.g. force to acceleration on system position. Defaults to eye matrix.
             * @param input Input to the system. E.g. force.
            */
            template<size_t CCOLS = LEN> 
            void predict(const VCTR::Math::Matrix<TYPE, LEN, CCOLS>& inputModel = 1, const VCTR::Math::Vector<TYPE, CCOLS>& input = 0);

            /**
             * Updates the system estimation/prediction using the given observation model and input/observation/sensor.
             * @param observationModel Models the system state vector to the expected observation/sensor output.
             * @param input Observation/sensor used to correct system. Is a ValueCov object so expects the covariance. Easy use by -> ValueCov{input, covariance}.
            */
            template<size_t UCOLS> 
            void update(const VCTR::Math::Matrix<TYPE, LEN, UCOLS>& observationModel, const ValueCov<TYPE, UCOLS>& input);

            /**
             * Sets the transition matrix for prediction-
             * @param transition Transition matrix. Transforms state into new state prediction.
            */
            void setTransitionMatrix(const VCTR::Math::Matrix<TYPE, LEN, LEN>& transition);

            /**
             * Sets the transition matrix for prediction-
             * @param transition Transition matrix. Transforms state into new state prediction.
            */
            void setProcessNoise(const VCTR::Math::Matrix<TYPE, LEN, LEN>& noise);

            /**
             * @returns estimation for state.
            */
            const VCTR::Math::Vector<TYPE, LEN>& getState();

            /**
             * @returns estimation for covariance.
            */
            const VCTR::Math::Matrix<TYPE, LEN, LEN>& getCovariance();

            /**
             * Sets the current estimation for the state and covariance.
             * @param state The state vector and covariance to be used.
            */
            void setState(const ValueCov<TYPE, LEN>& state);

        };

        template<size_t LEN, typename TYPE>
        KalmanFilter<LEN, TYPE>::KalmanFilter(const VCTR::Math::Matrix<TYPE, LEN, LEN>& transition, const VCTR::Math::Matrix<TYPE, LEN, LEN>& noise) {
            setTransitionMatrix(transition);
            setProcessNoise(noise);
            setState(ValueCov<TYPE, LEN>{0, 1000});
        }
        
        template<size_t LEN, typename TYPE>
        template<size_t CCOLS> 
        void KalmanFilter<LEN, TYPE>::predict(const VCTR::Math::Matrix<TYPE, LEN, CCOLS>& inputModel, const VCTR::Math::Vector<TYPE, CCOLS>& input) {
            x_ = f_ * x_ + inputModel * input;
            p_ = f_ * p_ * f_.transpose() + q_;
        }

        template<size_t LEN, typename TYPE>
        template<size_t UCOLS> 
        void KalmanFilter<LEN, TYPE>::update(const VCTR::Math::Matrix<TYPE, LEN, UCOLS>& observationModel, const ValueCov<TYPE, UCOLS>& input) {
            auto y = input.val - observationModel * x_;
            auto s = observationModel * p_ * observationModel.transpose() + input.cov;
            auto k = p_ * observationModel.transpose() * s.inverse();

            x_ = x_ + k * y;
            p_ = (VCTR::Math::Matrix<TYPE, LEN, LEN>::eye() - k * observationModel) * p_;

            p_ = (p_ + p_.transpose())*0.5f;
        }

        template<size_t LEN, typename TYPE>
        void KalmanFilter<LEN, TYPE>::setTransitionMatrix(const VCTR::Math::Matrix<TYPE, LEN, LEN>& transition) {
            f_ = transition;
        }

        template<size_t LEN, typename TYPE>
        void KalmanFilter<LEN, TYPE>::setProcessNoise(const VCTR::Math::Matrix<TYPE, LEN, LEN>& noise) {
            q_ = noise;
        }

        template<size_t LEN, typename TYPE>
        const VCTR::Math::Vector<TYPE, LEN>& KalmanFilter<LEN, TYPE>::getState() {
            return x_;
        }

        template<size_t LEN, typename TYPE>
        const VCTR::Math::Matrix<TYPE, LEN, LEN>& KalmanFilter<LEN, TYPE>::getCovariance() {
            return p_;
        }

        template<size_t LEN, typename TYPE>
        void KalmanFilter<LEN, TYPE>::setState(const ValueCov<TYPE, LEN>& state) {
            x_ = state.val;
            p_ = state.cov;
        }

    }
}

#endif