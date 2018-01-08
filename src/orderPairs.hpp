#ifndef CZI_EXPRESSION_MATRIX2_ORDER_PAIRS_HPP
#define CZI_EXPRESSION_MATRIX2_ORDER_PAIRS_HPP



namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        // Function object to order pairs in increasing order
        // of the second item in the pair, breaking ties
        // using the first item in the pair, in increasing order.
        template<class Pair> class OrderPairsBySecondThenByFirst;

        // Function object to order pairs in decreasing order
        // of the second item in the pair, breaking ties
        // using the first item in the pair, in increasing order.
        template<class Pair> class OrderPairsBySecondGreaterThenByFirstLess;

        // Function object to order pairs in decreasing order
        // of the second item in the pair, neglecting the first item in the pair.
        template<class Pair> class OrderPairsBySecondGreater;

        // Function object to order pairs in increasing order
        // of the first item in the pair, neglecting the second item in the pair.
        // (two pairs are considered equal if they differ only in the second item).
        template<class Pair> class OrderPairsByFirstOnly;
    }
}



template<class Pair> class ChanZuckerberg::ExpressionMatrix2::OrderPairsBySecondThenByFirst {
public:
    bool operator()(const Pair& x, const Pair& y) const
    {
        if(x.second < y.second) return true;
        if(y.second < x.second) return false;
        return x.first < y.first;
    }
};



template<class Pair> class ChanZuckerberg::ExpressionMatrix2::OrderPairsBySecondGreaterThenByFirstLess {
public:
    bool operator()(const Pair& x, const Pair& y) const
    {
        if(x.second > y.second) return true;
        if(y.second > x.second) return false;
        return x.first < y.first;
    }
};



template<class Pair> class ChanZuckerberg::ExpressionMatrix2::OrderPairsBySecondGreater {
public:
    bool operator()(const Pair& x, const Pair& y) const
    {
         return x.second > y.second;
    }
};



template<class Pair> class ChanZuckerberg::ExpressionMatrix2::OrderPairsByFirstOnly {
public:
    bool operator()(const Pair& x, const Pair& y) const
    {
         return x.first < y.first;
    }
};

#endif

