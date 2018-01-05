#ifndef CZI_EXPRESSION_MATRIX2_NORMALIZATION_METHOD_HPP
#define CZI_EXPRESSION_MATRIX2_NORMALIZATION_METHOD_HPP

#include "string.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {



        enum class NormalizationMethod {
            none,       // Using None conflicts with Python use of None as a language keyword
            L1,
            L2,
            Invalid
        };

        // This can be used to loop over the valid values of NormalizationMethod.
        // The unused attribute is necessary to suppress compilation warnings (due to -Wall).
        const auto validNormalizationMethods __attribute__((unused)) =
        {
            NormalizationMethod::none,
            NormalizationMethod::L1,
            NormalizationMethod::L2
        };



    	// Convert a NormalizationMethod to a short string and vice versa.
        inline string normalizationMethodToShortString(NormalizationMethod m)
        {
            switch(m) {
            case NormalizationMethod::none:
                return "none";
            case NormalizationMethod::L1:
                return "L1";
            case NormalizationMethod::L2:
                return "L2";
            default:
                return "Invalid";
            }
        }
        inline NormalizationMethod normalizationMethodFromShortString(const string& s)
        {
            if(s == "none") {
                return NormalizationMethod::none;
            } else if(s == "L1") {
                return NormalizationMethod::L1;
            } else if(s == "L2") {
                return NormalizationMethod::L2;
            } else {
                return NormalizationMethod::Invalid;
            }
        }



        // Convert a NormalizationMethod to a long descriptive string.
        inline string normalizationMethodToLongString(NormalizationMethod m)
        {
            switch(m) {
            case NormalizationMethod::none:
                return "no normalization (raw read count)";
            case NormalizationMethod::L1:
                return "L1 normalization (fractional read count)";
            case NormalizationMethod::L2:
                return "L2 normalization";
            default:
                return "Invalid normalization";
            }
        }


    }
}

#endif
