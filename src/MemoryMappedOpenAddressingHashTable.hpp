#ifndef CZI_EXPRESSION_MATRIX2_MEMORY_MAPPED_OPEN_ADDRESSING_HASH_TABLE_HPP
#define CZI_EXPRESSION_MATRIX2_MEMORY_MAPPED_OPEN_ADDRESSING_HASH_TABLE_HPP

// Template class MemoryMapped::OpenAddressingHashTable implements a simple hash table
// with open addressing and linear probing stored in a memory mapped file
// (using class MemoryMapped::Vector).
// See the Wikipedia article for background information:
// https://en.wikipedia.org/wiki/Open_addressing

// For good performance, the load factor should be kept well before 0.5.
// Note that attempting to add an element will result in an endless
// loop if the table is already full.
// For performance, we don't check for that.

// THIS IS UNTESTED CODE. IT IS CURRENTLY UNUSED. *****************************************

#include "nextPowerOfTwo.hpp"

#include "cstddef.hpp"
#include "vector.hpp"



namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        namespace MemoryMapped {
            template<class T, class Hasher> class OpenAddressingHashTable;
        }
    }
}


// The hashes type must have an operator()(const T&) that computes the hash function to be used.
template<class T, class Hasher> class ChanZuckerberg::ExpressionMatrix2::MemoryMapped::OpenAddressingHashTable {
public:

    // Create an empty table with at least n slots.
    // The actual number of slots is the lowest power of 2
    // greater or equal to n.
    void createNew(const string& name, size_t n, const Hasher*);

    // Access a previously created table.
    // The actual number of slots is the lowest power of 2
    // greater or equal to n.
    void accessExistingReadOnly(const string& name, const Hasher*);
    void accessExistingReadWrite(const string& name, const Hasher*);

    // Insert an element in the table.
    void insert(const T&);

    // Insert an element in the table, if it is not already present.
    // Return true if the element was inserted.
    bool insertIfNotPresent(const T&);

    // A slot of the table contains an object of type T, plus a bool
    // flag to indicate is the slot is occupied.
    class Slot {
    public:
        T data;
        bool isOccupied = false;
    };

    // The table where the data are stored.
    // To iterate over elements in the table, iterate over the slots in this vector
    // that have isOccupied set to true.
    Vector<Slot> data;

    // Mark all slots as unoccupied.
    void clear();

    // Set the in use count to zero, but without clearing any in-use flags.
    // The caller is also responsible for clearing all the in use flags.
    void clearInUseCount();

    // Return the number of elements currently stored in the table.
    // This is expensive as it loops over all slots in the table.
    size_t size() const
    {
        size_t inUseCount = 0;
        for(const Slot& slot: data) {
            if(slot.isOccupied) {
                ++inUseCount;
            }
        }
        return inUseCount;
    }

    // Return the maximum number of elements that can be stored in the table.
    // In practice it is good to stay below half this value,
    // because otherwise open addressing becomes notoriously slow.
    size_t capacity() const
    {
        return data.size();
    }

    private:

    // The mask used to convert a hash value to a slot index.
    uint64_t mask;

    // The number of slots currently occupied, which also
    // equals the number of elements currently stored in the table.
    uint64_t inUseCount;

    // The object used to compute hash functions.
    const Hasher* hasher;
};




// Create an empty OpenAddressingHashTable with at least n slots.
// The actual number of slots is the lowest power of 2
// greater or equal to n.
// Note that attempting to add an element will result in an endless
// loop if the tabke is already full.
// For performance, we don't check for that.
template<class T, class Hasher> inline void ChanZuckerberg::ExpressionMatrix2::MemoryMapped::OpenAddressingHashTable<T, Hasher>::createNew(
    const string& name, size_t n, const Hasher* hasherArgument)
{
    hasher = hasherArgument;
    n = nextPowerOfTwoGreaterThanOrEqual(n);
    data.createNew(name, n);
    data.resize(n);
    mask = n - 1ULL;
}



// Insert an element in the table, if it is not already present.
// Return true if the element was inserted.
template<class T, class Hasher> inline
    bool ChanZuckerberg::ExpressionMatrix2::MemoryMapped::OpenAddressingHashTable<T, Hasher>::insertIfNotPresent(const T& t)
{
    const uint64_t hashFunction = hasher(t);
    uint64_t slotIndex = hashFunction & mask;


    while(true) {
        Slot& slot = data[slotIndex];
        if(slot.isOccupied) {

            if(t == slot.data) {
                // This slot is occupied, and contains an element identical to t.
                return false;
            } else {

                // This slot is occupied, and contains an element different from t.
                // Try the next slot.
                ++slotIndex;
                slotIndex &= mask;
            }
        } else {

            // This slot is not occupied. Store this t here.
            slot.data = t;
            slot.isOccupied = true;
            ++inUseCount;
            return true;
        }
    }

}




// Mark all slots as unoccupied.
template<class T, class Hasher> inline
    void ChanZuckerberg::ExpressionMatrix2::MemoryMapped::OpenAddressingHashTable<T, Hasher>::clear()
{
    for(Slot& slot: data) {
        slot.isOccupied = false;
    }
    inUseCount = 0ULL;
}



// Set the in use count to zero.
// This means that the caller is also responsible for cleraring all the in use flags.
template<class T, class Hasher> inline
    void ChanZuckerberg::ExpressionMatrix2::MemoryMapped::OpenAddressingHashTable<T, Hasher>::clearInUseCount()
{
    inUseCount = 0ULL;
}



#endif
