#ifndef CZI_EXPRESSION_MATRIX2_HEAP_HPP
#define CZI_EXPRESSION_MATRIX2_HEAP_HPP

// Non -recursive implementation of some heap functionality.

// Adapted from code found at the following locations:
// https://www.codeproject.com/Tips/816934/Min-Binary-Heap-Implementation-in-Cplusplus
// https://gist.github.com/truncs/1810804


#include "cstddef.hpp"
// #include "iostream.hpp"
#include <algorithm>
#include <utility>
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        // Replacement for std::pop_heap followed by std::push_heap
        // used in the standard way.
        template<class Iterator, class T, class Comparator> inline void popAndPushHeap(
            Iterator begin,
            Iterator end,
            const T&,
            const Comparator&);


        // Return true if the given range is a heap using the given comparator.
        template<class Iterator, class Comparator> bool isHeap(
            Iterator begin,
            Iterator end,
            const Comparator&);

        // Given a vector, keep only the k best items in no particular order.
        // Best is defined by the given comparator.
        template<class T, class Comparator> void keepBest(vector<T>&, size_t k, const Comparator&);

        // Unit tests.
        void testHeap();
        void testKeepBest();
    }
}







// Return true if the given range is a min-heap.
template<class Iterator, class Comparator> bool ChanZuckerberg::ExpressionMatrix2::isHeap(
    Iterator begin,
    Iterator end,
    const Comparator& comparator)
{
    const size_t n = end - begin;
    for(size_t childIndex=1; childIndex<n; childIndex++) {
        const auto& child = *(begin + childIndex);
        const size_t parentIndex = (childIndex-1) / 2;
        const auto& parent = *(begin + parentIndex);
        if(!comparator(child, parent)) {
            return false;
        }
    }
    return true;
}



// Replacement for std::pop_heap followed by std::push_heap
// used in the standard way.
// Adapted from https://gist.github.com/truncs/1810804
template<class Iterator, class T, class Comparator>
    inline void ChanZuckerberg::ExpressionMatrix2::popAndPushHeap(
    Iterator begin,
    Iterator end,
    const T& t,
    const Comparator& comparator)
{
    // Put the new item in place of the top item.
    *begin = t;

    // Recursive swaps down the tree.
    const size_t n = end - begin;
    size_t parent = 0;

    while (true) {
        const size_t left = 2*parent + 1;
        const size_t right = left + 1;
        size_t best = parent;

        if(left<n && !comparator(*(begin+left), *(begin+best))) {
            best = left;
        }

        if(right<n && !comparator(*(begin+right), *(begin+best))) {
            best = right;
        }

        if(best != parent) {
            // cout << "Swap " << *(begin+parent) << " " << *(begin+best) << endl;
            std::swap(*(begin+parent), *(begin+best));
            parent = best;
        }
        else {
            break;
        }
  }
}



// Given a vector, keep only the k best items in no particular order.
// "Best" is as defined by the given comparator.
template<class T, class Comparator>
    void ChanZuckerberg::ExpressionMatrix2::keepBest(
    vector<T>& v,
    size_t k,
    const Comparator& comparator)
{
    std::nth_element(v.begin(), v.begin()+k, v.end(), comparator);
    v.resize(k);
}


#endif
