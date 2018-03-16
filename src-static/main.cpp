#include "ExpressionMatrix.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

int main()
{
    try {
        ExpressionMatrix e("data");
        e.explore(17100, "");
    } catch(runtime_error& e) {
        cout << e.what();
        return 1;
    }
    
    return 0;
}
