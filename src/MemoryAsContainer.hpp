// Simple class to make a range objects in memory appear as a container.

#ifndef CZI_EXPRESSION_MATRIX2_MEMORY_AS_CONTAINER_HPP
#define CZI_EXPRESSION_MATRIX2_MEMORY_AS_CONTAINER_HPP

#include "cstddef.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
    template<class T> class MemoryAsContainer;
    // template<class T> class ConstMemoryAsContainer;
    }
}



template<class T> class ChanZuckerberg::ExpressionMatrix2::MemoryAsContainer {
public:

    MemoryAsContainer(T* begin, T* end) :
        dataBegin(begin),
        dataEnd(end)
    {
    }

    size_t size() const
    {
        return dataEnd - dataBegin;
    }
    T* begin() const
    {
        return dataBegin;
    }
    T* end() const
    {
        return dataEnd;
    }
    T& operator[](size_t i) const
    {
        return dataBegin[i];
    }

private:
    T* dataBegin;
    T* dataEnd;
};


#if 0
template<class T> class ChanZuckerberg::ExpressionMatrix2::ConstMemoryAsContainer {
public:

    ConstMemoryAsContainer(const T* begin, const T* end) :
        dataBegin(begin),
        dataEnd(end)
    {
    }

    size_t size() const
    {
        return dataEnd - dataBegin;
    }
    const T* begin() const
    {
        return dataBegin;
    }
    const T* end() const
    {
        return dataEnd;
    }
    const T& operator[](size_t i) const
    {
        return dataBegin[i];
    }

private:
    const T* dataBegin;
    const T* dataEnd;
};
#endif


#endif
