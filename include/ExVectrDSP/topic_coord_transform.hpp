#ifndef ExVectrDSP_COORDINATETRANSFORM_H
#define ExVectrDSP_COORDINATETRANSFORM_H

#include "ExVectrCore/timestamped.hpp"
#include "ExVectrCore/topic.hpp"
#include "ExVectrCore/topic_subscribers.hpp"
#include "ExVectrMath/matrix_base.hpp"

#include "value_covariance.hpp"

namespace VCTR
{

    namespace DSP
    {

        template<typename TYPE>
        class TopicCoordTransform
        {
        private:
            Math::Matrix<TYPE, 3, 3> scaleMat_;
            Math::Matrix<TYPE, 3, 1> biasVec_;

            Core::Callback_Subscriber<Core::Timestamped<ValueCov<TYPE, 3>>, TopicCoordTransform> subr_;
            Core::Topic<Core::Timestamped<ValueCov<TYPE, 3>>> outputTopic_;

        public:
            TopicCoordTransform(const Math::Matrix<TYPE, 3, 3> &scaleMat = 1, const Math::Matrix<TYPE, 3, 1> &biasVec = 0);

            TopicCoordTransform(const Math::Matrix<TYPE, 3, 3> &scaleMat, const Math::Matrix<TYPE, 3, 1> &biasVec, const Core::Topic<Core::Timestamped<ValueCov<TYPE, 3>>>& inputTopic);

            void setTransform(const Math::Matrix<TYPE, 3, 3> &scaleMat, const Math::Matrix<TYPE, 3, 1> &biasVec);

            void setInputTopic(Core::Topic<Core::Timestamped<ValueCov<TYPE, 3>>>& topic);

            Core::Topic<Core::Timestamped<ValueCov<TYPE, 3>>>& getOutputTopic();

        private:

            void receive(const Core::Timestamped<ValueCov<TYPE, 3>>& item);

        };

        template<typename TYPE>
        TopicCoordTransform<TYPE>::TopicCoordTransform(const Math::Matrix<TYPE, 3, 3> &scaleMat, const Math::Matrix<TYPE, 3, 1> &biasVec) {
            setTransform(scaleMat, biasVec);
            subr_.setCallback(this, &TopicCoordTransform::receive);
        }

        template<typename TYPE>
        TopicCoordTransform<TYPE>::TopicCoordTransform(const Math::Matrix<TYPE, 3, 3> &scaleMat, const Math::Matrix<TYPE, 3, 1> &biasVec, const Core::Topic<Core::Timestamped<ValueCov<TYPE, 3>>>& inputTopic)
         : TopicCoordTransform(scaleMat, biasVec) {
            scaleMat_ = scaleMat;
            biasVec_ = biasVec;
            setInputTopic(inputTopic);
        }

        template<typename TYPE>
        void TopicCoordTransform<TYPE>::setTransform(const Math::Matrix<TYPE, 3, 3> &scaleMat, const Math::Matrix<TYPE, 3, 1> &biasVec) {
            scaleMat_ = scaleMat;
            biasVec_ = biasVec;
        }

        template<typename TYPE>
        void TopicCoordTransform<TYPE>::setInputTopic(Core::Topic<Core::Timestamped<ValueCov<TYPE, 3>>>& topic) {
            //subr_.unsubcribe();
            subr_.subscribe(topic);
        }

        template<typename TYPE>
        Core::Topic<Core::Timestamped<ValueCov<TYPE, 3>>>& TopicCoordTransform<TYPE>::getOutputTopic() {
            return outputTopic_;
        }
        
        template<typename TYPE>
        void TopicCoordTransform<TYPE>::receive(const Core::Timestamped<ValueCov<TYPE, 3>>& item) {
            auto val = item.data;
            val.val = scaleMat_ * (val.val - biasVec_);
            outputTopic_.publish(Core::Timestamped<ValueCov<TYPE, 3>>(val, item.timestamp));
        }

    }
}

#endif