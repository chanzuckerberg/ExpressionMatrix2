#ifndef CZI_EXPRESSION_MATRIX2_SHORT_STATIC_STRING_HPP
#define CZI_EXPRESSION_MATRIX2_SHORT_STATIC_STRING_HPP

// Class to store a short string without using any allocated memory,
// so it can be used in shared memory.
// We could use boost::static_vector<char>, but
// this first appeared in Boost 1.54, and we want to keep the ability
// to build on CentOS 7 which uses Boost 1.53.

#include "algorithm.hpp"
#include "array.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"
#include "string.hpp"


namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        template<size_t capacity> class ShortStaticString;
        void testShortStaticString();
        using StaticString255 = ShortStaticString<255>;
    }
}



template<size_t capacity> class ChanZuckerberg::ExpressionMatrix2::ShortStaticString {
public:

    // The number of characters actually stored in the string.
    // This can be at most equal to the specified capacity.
    using CapacityType = uint8_t;
    CapacityType n;

    // Sanity check that we can represent the number of characters.
    static_assert(capacity <= std::numeric_limits<CapacityType>::max(),
        "ShortStaticString capacity exceeded.");

    // The characters of the string.
    // Only the first n are significant.
    array<char, capacity> s;



    // Begin/end iterators.
    using iterator = typename array<char, capacity>::iterator;
    using const_iterator = typename array<char, capacity>::const_iterator;
    iterator begin()
    {
        return s.begin();
    }
    const_iterator begin() const
    {
        return s.begin();
    }
    iterator end()
    {
        return begin() + n;
    }
    const_iterator end() const
    {
        return begin() + n;
    }



    // The default constructor creates an empty string.
    ShortStaticString() : n(0)
    {
    }

    // Constructor from a string.
    ShortStaticString(const string& x)
    {
        setSize(x);
        copy(x.begin(), x.end(), begin());
    }

    // Assignment from a string.
    ShortStaticString<capacity>& operator=(const string& x) {
        setSize(x);
        copy(x.begin(), x.end(), begin());
        return *this;
    }

    // Conversion to a string.
    operator string() const
    {
        string x;
        x.resize(size());
        copy(begin(), end(), x.begin());
        return x;
    }

    size_t size() const
    {
        return n;
    }

    size_t empty() const
    {
        return n==0;
    }

    // Compare with another ShortStaticString.
    template<size_t thatCapacity> bool operator==(const ShortStaticString<thatCapacity>& that) const
    {
        return
            (size() == that.size()) &&
            (std::equal(begin(), end(), that.begin()));
    }
    template<size_t thatCapacity> bool operator!=(const ShortStaticString<thatCapacity>& that) const
    {
        return !(*this == that);
    }

    // Compare with a string.
    bool operator==(const string& that) const
    {
        return
            (size() == that.size()) &&
            (std::equal(begin(), end(), that.begin()));
    }


private:

    void setSize(size_t N)
    {
        if(N > capacity) {
            throw runtime_error("ShortStaticString capacity exceeded.");
        }
        n = CapacityType(N);
    }
    void setSize(const string& x)
    {
        setSize(x.size());
    }


};



// Output operator.
template<size_t capacity> std::ostream& operator<<(
    std::ostream& s,
    const ChanZuckerberg::ExpressionMatrix2::ShortStaticString<capacity>& x)
{
    for(char c: x) {
        s << c;
    }
    return s;
}


inline void ChanZuckerberg::ExpressionMatrix2::testShortStaticString()
{

    const size_t capacity = 100;

    // Default construction.
    ShortStaticString<capacity> x;
    CZI_ASSERT(x.size() == 0);
    CZI_ASSERT(x.empty());

    // Construction from an array of char.
    const ShortStaticString<capacity> y("abc");
    CZI_ASSERT(y.size() == 3);

    // Construction from a std::string.
    const ShortStaticString<capacity> z(string("abcdef"));
    CZI_ASSERT(z.size() == 6);

    // Construction from a std::string that is too long.
    // This must cause an exception.
    bool caughtException = false;
    try {
        const ShortStaticString<capacity> u(string(1000, 'a'));
    } catch(std::runtime_error&) {
        caughtException = true;
    }
    CZI_ASSERT(caughtException);

    // Assignment from an array of characters.
    x = "abcd";
    CZI_ASSERT(x.size() == 4);

    // Assignment from a std::string.
    x = string("abcde");
    CZI_ASSERT(x.size() == 5);

    // Assignment from another ShortStaticString.
    x = y;
    CZI_ASSERT(x.size() == y.size());

    // Comparison with another ShortStaticString.
    CZI_ASSERT(x == y);
    CZI_ASSERT(x != z);

    // Comparison with a std::string.
    CZI_ASSERT(x == string("abc"));

    // Conversion to a std::string.
    const string X = x;
    CZI_ASSERT(X == "abc");

    cout << x << endl;

}

#endif
