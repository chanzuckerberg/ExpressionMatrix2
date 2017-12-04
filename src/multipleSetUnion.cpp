#include "multipleSetUnion.hpp"
#include "iostream.hpp"
#include "iterator.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



void ChanZuckerberg::ExpressionMatrix2::multipleSetUnionTest()
{
    const vector<int> x0 = {3, 7, 10};
    const vector<int> x1 = {2, 7, 25};
    const vector<int> x2 = {7, 10};
    const vector<int> x3 = {3, 8 ,25, 40};

    vector<const vector<int>* > pointers = {&x0, &x1, &x2, &x3};
    vector<int> y;
    multipleSetUnion(pointers, y);

    copy(y.begin(), y.end(), ostream_iterator<int>(cout, " "));
    cout << endl;

}
