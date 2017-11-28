#include "heap.hpp"
#include "CZI_ASSERT.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

#include "iostream.hpp"
#include "iterator.hpp"
#include "string.hpp"
#include "vector.hpp"

static void dumpMinHeap(const vector<int>& v, const string& s)
{
    cout << s << endl;
    for(size_t i=0; i<v.size(); i++) {
        cout << i << " " << v[i] << endl;
    }
    CZI_ASSERT(isHeap(v.begin(), v.end(), std::greater<int>()));

}

void ChanZuckerberg::ExpressionMatrix2::testHeap()
{
    vector<int> v = {35, 9, 14, 39, 17, 10, 18, 28, 19, 36, 7, 43, 16};
    std::make_heap(v.begin(), v.end(), std::greater<int>());
    dumpMinHeap(v, "Initial:");

    popAndPushHeap(v.begin(), v.end(), 20, std::greater<int>());
    dumpMinHeap(v, "After popAndPush:");
}


void ChanZuckerberg::ExpressionMatrix2::testKeepBest()
{
    vector<int> v = {35, 9, 14, 39, 17, 10, 18, 28, 19, 36, 7, 43, 16};
    keepBest(v, 6, std::greater<int>());
    copy(v.begin(), v.end(), ostream_iterator<int>(cout, " "));
    cout << endl;
}
