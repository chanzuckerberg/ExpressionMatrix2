#include "MemoryMappedStringTable.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

void ChanZuckerberg::ExpressionMatrix2::testMemoryMappedStringTable()
{

    MemoryMapped::StringTable<uint32_t> stringTable;
    stringTable.createNew("StringTableTest", 6);

    const string s0 = "abcde";
    const uint32_t i0 = stringTable[s0];
    const string s1 = "12345";
    const uint32_t i1 = stringTable[s1];
    cout << i0 << " " << i1 << endl;
    const string s0Check = stringTable[i0];
    const string s1Check = stringTable[i1];
    CZI_ASSERT(s0Check == s0);
    CZI_ASSERT(s1Check == s1);
    CZI_ASSERT(stringTable[s0] == i0);
    CZI_ASSERT(stringTable[s1] == i1);

    cout << stringTable["aaa"] << endl;
    cout << stringTable["bbb"] << endl;
    cout << stringTable["ccc"] << endl;
    cout << stringTable["aaa"] << endl;
    cout << stringTable["bbb"] << endl;
    cout << stringTable("ccc") << endl;

}
