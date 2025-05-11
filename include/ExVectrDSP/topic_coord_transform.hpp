#ifndef EXVECTRDSP_COORDINATETRANSFORM_HPP
#define EXVECTRDSP_COORDINATETRANSFORM_HPP

#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"
#include "ExVectrMath/matrix_base.hpp"

#include "value_covariance.hpp"

namespace VCTR
{

    namespace DSP
    {

        /**
         * @brief Class to transform a 3D vector using a scale and bias matrix. Used primarily for transforming sensor data from sensor frame to body frame and also for applying sensor calibration. 
         * @tparam TYPE Type of the data. (float, double, etc.)
         * @tparam DIM Dimension of the data. (1, 3, 6, etc.)
         */
        template<typename TYPE, size_t DIM = 3>
        class TopicCoordTransform
        {
        private:
            Math::Matrix<TYPE, DIM, DIM> transformMat_;
            Math::Matrix<TYPE, DIM, 1> biasVec_;

            Core::Callback_Subscriber<Core::Timestamped<ValueCov<TYPE, DIM>>, TopicCoordTransform> subr_;
            Core::Topic<Core::Timestamped<ValueCov<TYPE, DIM>>> outputTopic_;

        public:

            /**
             * @brief Constructor that sets the transform matrix and bias vector for the transformation.
             * @note The calibration is done via this formula: val = transformMat * (val - biasVec).
             */
            TopicCoordTransform(const Math::Matrix<TYPE, DIM, DIM> &transformMat = 1, const Math::Matrix<TYPE, DIM, 1> &biasVec = 0);

            /**
             * @brief Constructor that sets the transform matrix and bias vector for the transformation.
             * @note The calibration is done via this formula: val = transformMat * (val - biasVec).
             */
            TopicCoordTransform(const Math::Matrix<TYPE, DIM, DIM> &transformMat, const Math::Matrix<TYPE, DIM, 1> &biasVec, const Core::Topic<Core::Timestamped<ValueCov<TYPE, DIM>>>& inputTopic);

            /**
             * @brief Sets the transform matrix and bias vector for the transformation.
             * @note The calibration is done via this formula: val = transformMat * (val - biasVec).
             */
            void setTransform(const Math::Matrix<TYPE, DIM, DIM> &transformMat, const Math::Matrix<TYPE, DIM, 1> &biasVec);

            void setInputTopic(Core::Topic<Core::Timestamped<ValueCov<TYPE, DIM>>>& topic);

            Core::Topic<Core::Timestamped<ValueCov<TYPE, DIM>>>& getOutputTopic();

        private:

            void receive(const Core::Timestamped<ValueCov<TYPE, DIM>>& item);

        };

        template<typename TYPE, size_t DIM>
        TopicCoordTransform<TYPE, DIM>::TopicCoordTransform(const Math::Matrix<TYPE, DIM, DIM> &transformMat, const Math::Matrix<TYPE, DIM, 1> &biasVec) {
            setTransform(transformMat, biasVec);
            subr_.setCallback(this, &TopicCoordTransform::receive);
        }

        template<typename TYPE, size_t DIM>
        TopicCoordTransform<TYPE, DIM>::TopicCoordTransform(const Math::Matrix<TYPE, DIM, DIM> &transformMat, const Math::Matrix<TYPE, DIM, 1> &biasVec, const Core::Topic<Core::Timestamped<ValueCov<TYPE, DIM>>>& inputTopic)
         : TopicCoordTransform(transformMat, biasVec) {
            transformMat_ = transformMat;
            biasVec_ = biasVec;
            setInputTopic(inputTopic);
        }

        
        template<typename TYPE, size_t DIM>
        void TopicCoordTransform<TYPE, DIM>::setTransform(const Math::Matrix<TYPE, DIM, DIM> &transformMat, const Math::Matrix<TYPE, DIM, 1> &biasVec) {
            transformMat_ = transformMat;
            biasVec_ = biasVec;
        }

        template<typename TYPE, size_t DIM>
        void TopicCoordTransform<TYPE, DIM>::setInputTopic(Core::Topic<Core::Timestamped<ValueCov<TYPE, DIM>>>& topic) {
            //subr_.unsubcribe();
            subr_.subscribe(topic);
        }

        template<typename TYPE, size_t DIM>
        Core::Topic<Core::Timestamped<ValueCov<TYPE, DIM>>>& TopicCoordTransform<TYPE, DIM>::getOutputTopic() {
            return outputTopic_;
        }
        
        template<typename TYPE, size_t DIM>
        void TopicCoordTransform<TYPE, DIM>::receive(const Core::Timestamped<ValueCov<TYPE, DIM>>& item) {
            auto val = item.data;
            val.val = transformMat_ * (val.val - biasVec_);
            outputTopic_.publish(Core::Timestamped<ValueCov<TYPE, DIM>>(val, item.timestamp));
        }

    }
}

#endif