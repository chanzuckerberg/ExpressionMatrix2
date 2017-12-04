#ifndef CZI_EXPRESSION_MATRIX2_MULTIPLE_SET_UNION_HPP
#define CZI_EXPRESSION_MATRIX2_MULTIPLE_SET_UNION_HPP

// Efficient computation of the union of multiple sets.
// Each of the input sets is stored in a sorted container,
// and so is the result set.

#include "CZI_ASSERT.hpp"
#include "cstddef.hpp"
#include <queue>
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        template<class Container> inline void multipleSetUnion(
            const vector<const Container*>& inputSets,
            Container& outputSet);

        template<class Container> class MultipleSetUnionHelper {
        public:
            size_t containerId;
            typename Container::const_iterator it;
            typename Container::value_type value;
            MultipleSetUnionHelper(
                size_t containerId,
                typename Container::const_iterator it,
                typename Container::value_type value
                ) :
                containerId(containerId), it(it), value(value)
            {
            }
            bool operator<(const MultipleSetUnionHelper& that) const
            {
                return that.value < value;  // Reversed, because this is the way a priority queue works.
            }

        };

        void multipleSetUnionTest();
    }
}


template<class Container> inline void ChanZuckerberg::ExpressionMatrix2::multipleSetUnion(
    const vector<const Container*>& inputSets,
    Container& outputSet)
{
    CZI_ASSERT(outputSet.empty());

    // We use a priority q in which the top element is the
    // iterator pointing to the next element to be inserted.
    std::priority_queue< MultipleSetUnionHelper<Container> > q;
    for(size_t containerId=0; containerId<inputSets.size(); containerId++) {
        const Container* p = inputSets[containerId];
        CZI_ASSERT(p);
        const Container& container = *p;
        if(!container.empty()) {
            q.push(MultipleSetUnionHelper<Container>(containerId, container.begin(), container.front()));
        }
    }

    while(!q.empty()) {
        MultipleSetUnionHelper<Container> m = q.top();
        q.pop();
        if(outputSet.empty() || m.value != outputSet.back()) {
            outputSet.push_back(m.value);
        }
        const Container& container = *(inputSets[m.containerId]);
        ++m.it;
        if(m.it != container.end()) {
            m.value = *(m.it);
            q.push(m);
        }
    }

}




#endif
