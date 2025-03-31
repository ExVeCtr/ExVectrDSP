#ifndef ExVectrDSP_LPF_H
#define ExVectrDSP_LPF_H

#include "ExVectrCore/time_definitions.hpp"
#include "ExVectrCore/timestamped.hpp"

namespace VCTR
{

    namespace DSP
    {

        /**
         * A first order low-pass filter
         */
        template <typename TYPE = float>
        class LowPassFilter
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
            LowPassFilter(float cutoffFrequency, const TYPE &value = 0);

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
        LowPassFilter<TYPE>::LowPassFilter(float cutoffFrequency = 1.0f, const TYPE &value = 0)
        {
            value_ = value;
            rc_ = 1.0 / (cutoffFrequency * 2 * 3.14);
        }

        template <typename TYPE>
        const TYPE &LowPassFilter<TYPE>::update(const VCTR::Core::Timestamped<TYPE> &value)
        {
            return update(value.data, value.timestamp);
        }

        template <typename TYPE>
        const TYPE &LowPassFilter<TYPE>::update(const TYPE &value, int64_t sampleTime = VCTR::Core::NOW())
        {

            float dt = (double)(sampleTime - lastSampleTime_) / VCTR::Core::SECONDS;
            lastSampleTime_ = sampleTime;
            float a = dt / (rc_ + dt);

            value_ = lastInput_ + (value - lastInput_) * a;
            lastInput_ = value;

            return value_;
        }

        template <typename TYPE>
        void LowPassFilter<TYPE>::setCutoff(float cutoffFrequency) {
            rc_ = 1.0 / (cutoffFrequency * 2 * 3.14);
        }

        template <typename TYPE>
        const TYPE &LowPassFilter<TYPE>::getOutput()
        {
            return value_;
        }

        template <typename TYPE>
        void LowPassFilter<TYPE>::setValue(const TYPE &value)
        {
            lastInput_ = value_ = value;
        }


        

    }
}

#endif