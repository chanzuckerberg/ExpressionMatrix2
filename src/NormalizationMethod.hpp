#ifndef CZI_EXPRESSION_MATRIX2_NORMALIZATION_METHOD_HPP
#define CZI_EXPRESSION_MATRIX2_NORMALIZATION_METHOD_HPP

#include "string.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {



    	enum class NormalizationMethod {
    		None,
			L1,
			L2,
			Invalid
    	};



    	// Convert a NormalizationMethod to a short string and vice versa.
    	inline string normalizationMethodToShortString(NormalizationMethod m)
    	{
    		switch(m) {
    		case NormalizationMethod::None:
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
    		if(s=="none") {
    			return NormalizationMethod::None;
    		} else if(s=="L1") {
        		return NormalizationMethod::L1;
    		} else if(s=="L2") {
        		return NormalizationMethod::L2;
    		} else {
    			return NormalizationMethod::Invalid;
    		}
    	}



    	// Convert a NormalizationMethod to a long descriptive string.
    	inline string normalizationMethodToLongString(NormalizationMethod m)
    	{
    		switch(m) {
    		case NormalizationMethod::None:
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
