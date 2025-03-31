#ifndef ExVectrDSP_HPF_H
#define ExVectrDSP_HPF_H

#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/timestamped.hpp"

namespace VCTR
{

    namespace Data
    {

        /**
         * A first order high-pass filter.
         */
        template <typename TYPE = float>
        class HighPassFilter
        {
        private:
            /// @brief Current filter value.
            TYPE value_;
            /// @brief Current filter value.
            TYPE lastInput_;
            /// @brief Crossover frequency.
            float alpha_;
            /// @brief rc value.
            float rc_;
            /// @brief Last sample time.
            int64_t lastSampleTime_ = 0;

        public:
            /**
             * Constructor.
             * @param cutoffFrequency Frequencies above this will be damped.
             * @param value Initial value to be used.
             */
            HighPassFilter(float cutoffFrequency, const TYPE &value = 0);

            /**
             * Updates the values with new given sample.
             * @param value New sample.
             * @returns Filter output.
             */
            const TYPE &update(const TYPE &value, int64_t sampleTime = VCTR::Core::NOW());

            /**
             * Updates the values with new given sample.
             * @param value New sample with timestamp.
             * @returns Filter output.
             */
            const TYPE &update(const VCTR::Core::Timestamped<TYPE> &value);

            /**
             * Changes the cut off frequency.
             */
            void setCutoff(float cutoffFrequency);

            /**
             * @returns current filter output.
             */
            const TYPE &getOutput();

            /**
             * Sets the filter value to given value.
             * @param value The value to force output to be.
             */
            void setValue(const TYPE &value);
        };

        template <typename TYPE>
        HighPassFilter<TYPE>::HighPassFilter(float cutoffFrequency = 1.0f, const TYPE &value = 0)
        {
            value_ = value;
            rc_ = 1.0 / (cutoffFrequency * 2 * 3.14);
        }

        template <typename TYPE>
        const TYPE &HighPassFilter<TYPE>::update(const VCTR::Core::Timestamped<TYPE> &value)
        {
            return update(value.data, value.timestamp);
        }

        template <typename TYPE>
        const TYPE &HighPassFilter<TYPE>::update(const TYPE &value, int64_t sampleTime = VCTR::Core::NOW())
        {

            float dt = (double)(sampleTime - lastSampleTime_) / VCTR::Core::SECONDS;
            lastSampleTime_ = sampleTime;
            float a = rc_ / (rc_ + dt);

            value_ = (value_ + value - lastInput_) * a;
            lastInput_ = value;

            return value_;
        }

        template <typename TYPE>
        void HighPassFilter<TYPE>::setCutoff(float cutoffFrequency) {
            rc_ = 1.0 / (cutoffFrequency * 2 * 3.14);
        }

        template <typename TYPE>
        const TYPE &HighPassFilter<TYPE>::getOutput()
        {
            return value_;
        }

        template <typename TYPE>
        void HighPassFilter<TYPE>::setValue(const TYPE &value)
        {
            lastInput_ = value_ = value;
        }



    }
}

#endif