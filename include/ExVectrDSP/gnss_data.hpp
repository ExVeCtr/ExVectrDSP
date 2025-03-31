#ifndef ExVectrDSP_GNSSDATA_H
#define ExVectrDSP_GNSSDATA_H

#include "ExVectrMath/matrix_base.hpp"
#include "ExVectrMath/matrix_vector.hpp"

namespace VCTR
{

    namespace Data
    {

        /**
         * Simple container class for containing GNSS data. This should be used as its designed to maximise position accuracy (double) and reduce memory usage.
         * @note The position and velocity are in NED frame. 
         */
        struct GNSSData
        {
        public:

            /// @brief in form [latitude, longitude, altitude] lat, lon in radians altitude in meters above ellipsoid.
            Math::Vector_D position = 0;
            Math::Vector_F velocity = 0;

            Math::Vector_F positionCov = 1;
            Math::Vector_F velocityCov = 1;

            bool positionValid = false;
            bool velocityValid = false;

            uint8_t numSats = 0;

        };

    }
}

#endif