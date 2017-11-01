#ifndef CZI_EXPRESSION_MATRIX2_MEMORY_MAPPED_STRING_TABLE_HPP
#define CZI_EXPRESSION_MATRIX2_MEMORY_MAPPED_STRING_TABLE_HPP



/*******************************************************************************

Class StringTable maintains a bidirectional correspondence between
n distinct strings (a set of strings) and the integers 0 to n-1.

The table is persistent. This is achieved by storing it on memory mapped files.

Strings can only be added to the table, never removed.

Each of the n strings is assigned a string id of type StringId, which
should be an unsigned integer type. The choice of this type
determines the maximum number of strings that can be stored.

The class provides two versions of operator[] that can be used
to find the string corresponding to a given id and vice versa.

*******************************************************************************/



// CZI.
#include "MemoryAsContainer.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "nextPowerOfTwo.hpp"

// MurmurHash.
#include "MurmurHash2.hpp"

// Forward declarations.
namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        namespace MemoryMapped {
            template<class StringId> class StringTable;
        }
        void testMemoryMappedStringTable();
    }
}



template<class StringId> class ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable {
public:

    // Create a new StringTable.
    // The name will be used as a base name for memory mapped files.
    // The specified capacity should be at least double the expected number of
    // strings we are going to put in the table.
    void createNew(const string& name, size_t capacity);

    // Access an existing StringTable.
    void accessExistingReadOnly(const string& name);
    void accessExistingReadWrite(const string& name, bool allowReadOnly);

    // Return the StringId corresponding to a given string,
    // adding it to the table if not already present.
    // (this  requires write access).
    StringId operator[](const string&);

    // Return the StringId corresponding to a given string.
    // Returns invalidStringId if the table does not contain the given string
    static const StringId invalidStringId = std::numeric_limits<StringId>::max();
    StringId operator()(const string&) const;

    // Return the string corresponding to a given StringId.
    string operator[](StringId) const;

    // Return the memory range containing the string with a given StringId.
    MemoryAsContainer<const char> operator()(StringId) const;

    // Return true if the given StringId corresponds to the specified string.
    bool equal(StringId, const string&) const;

    // Return the number of strings currently stored in the table.
    size_t size() const;

    // Return the maximum number of strings that can be stored.
    // For performance, we should stay below half this value.
    size_t capacity() const;

    // The strings are stored using a MemoryMapped::VectorOfVectors.
    // The i-th string is stored in the open range of characters
    // strings.begin(i) through strings.end(i).
    // This is kept public to simplify iterating over the strings in the table.
    VectorOfVectors<char, StringId> strings;

private:

    // To be able to find the string id corresponding to a given string,
    // we use a hash table with open addressing. Each slot stores a string id.
    // An unused slot is filled with InvalidStringId.
    // See the Wikipedia article for background information:
    // https://en.wikipedia.org/wiki/Open_addressing.
    Vector<StringId> hashTable;
    uint64_t mask;

};



template<class StringId> const StringId ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::invalidStringId;



template<class StringId> inline
    void ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::createNew(
        const string& name,
        uint64_t capacity)
{
    // The actual capacity must be a power of two, for speed (this way modulo
    // operations are not needed).
    const uint64_t n = nextPowerOfTwoGreaterThanOrEqual(capacity);
    mask = n - 1ULL;

    // Start with no strings.
    strings.createNew(name + "-strings");

    // Initially fill the hash table with InvalidStringId.
    hashTable.createNew(name + "-hashTable", n);
    fill(hashTable.begin(), hashTable.end(), invalidStringId);
}



template<class StringId> inline
    void ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::accessExistingReadOnly(
        const string& name)
{
    strings.accessExistingReadOnly(name + "-strings");
    hashTable.accessExistingReadOnly(name + "-hashTable");
    mask = hashTable.size() - 1ULL;
}



template<class StringId> inline
    void ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::accessExistingReadWrite(
        const string& name,
        bool allowReadOnly)
{
    strings.accessExistingReadWrite(name + "-strings", allowReadOnly);
    hashTable.accessExistingReadWrite(name + "-hashTable", allowReadOnly);
    mask = hashTable.size() - 1ULL;
}



template<class StringId> inline
    size_t ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::size() const
{
    return strings.size();
}



template<class StringId> inline
    size_t ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::capacity() const
{
    return hashTable.size();
}



// Return the StringId corresponding to a given string, creating it if necessary.
template<class StringId> inline
    StringId ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::operator[](const string& s)
{
    uint64_t bucketIndex = MurmurHash64A(s.data(), int(s.size()), 237) & mask;

    while(true) {
        StringId& stringId = hashTable[bucketIndex];
        if(stringId == invalidStringId) {

            // The bucket is empty. This means that this string is not already in the table. Add it.
            stringId = StringId(strings.size());
            if(stringId > hashTable.size()/2) {
                throw runtime_error(
                    hashTable.fileName +
                    ": string table capacity " +
                    lexical_cast<string>(hashTable.size()) +
                    " exceeded.");
            }
            strings.appendVector(s.begin(), s.end());
            return stringId;

        } else {

            // The bucket is not empty. See if it contains the string we are looking for.

            if((s.size() == strings.size(stringId)) && std::equal(s.begin(), s.end(), strings.begin(stringId))) {
                // The bucket contains this string.
                return stringId;
            } else {
                // The bucket contains another string. Try the next bucket.
                ++bucketIndex;
                bucketIndex &= mask;
            }

        }
    }

    // Note that this will loop forever if the capacity of the hash table is exceeded.

}


// Return the StringId corresponding to a given string.
// Returns invalidStringId if the table does not contain the given string
template<class StringId> inline
    StringId ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::operator()(const string& s) const
{
    uint64_t bucketIndex = MurmurHash64A(s.data(), int(s.size()), 237) & mask;

    while(true) {
        const StringId& stringId = hashTable[bucketIndex];
        if(stringId == invalidStringId) {

            // The bucket is empty. This means that this string is not already in the table.
            return invalidStringId;

        } else {

            // The bucket is not empty. See if it contains the string we are looking for.

            if((s.size() == strings.size(stringId)) && std::equal(s.begin(), s.end(), strings.begin(stringId))) {
                // The bucket contains this string.
                return stringId;
            } else {
                // The bucket contains another string. Try the next bucket.
                ++bucketIndex;
                bucketIndex &= mask;
            }

        }
    }

    // Note that this will loop forever if the hash table is full and does not
    // contain the string we are looking for.

}



// Return the string corresponding to a given StringId.
template<class StringId> inline
    std::string ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::operator[](StringId stringId) const
{
    CZI_ASSERT(stringId < strings.size());
    const auto& s = strings[stringId];
    return string(s.begin(), s.end());
}


// Return the memory range containing the string with a given StringId.
template<class StringId> inline
    ChanZuckerberg::ExpressionMatrix2::MemoryAsContainer<const char>
    ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::operator()(StringId stringId) const
{
    return strings[stringId];
}


// Return true if the given StringId corresponds to the specified string.
template<class StringId> inline
    bool ChanZuckerberg::ExpressionMatrix2::MemoryMapped::StringTable<StringId>::equal(
        StringId stringId,
        const string& s) const
    {
        return
            strings.size(stringId) == s.size() &&
            std::equal(s.begin(), s.end(), strings.begin(stringId));
    }


#endif

