// A vector of lists stored in mapped memory.

#ifndef CZI_EXPRESSION_MATRIX2_MEMORY_MAPPED_VECTOR_OF_LISTS_HPP
#define CZI_EXPRESSION_MATRIX2_MEMORY_MAPPED_VECTOR_OF_LISTS_HPP

// CZI.
#include "MemoryMappedVector.hpp"

// Standard libraries.
#include "algorithm.hpp"
#include "cstdint.hpp"
#include <limits>
#include "string.hpp"
#include "utility.hpp"

// Forward declarations.
namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {
        namespace MemoryMapped {
            template<class T> class VectorOfLists;
        }
        inline void testMemoryMappedVectorOfLists();
    }
}



template<class T> class ChanZuckerberg::ExpressionMatrix2::MemoryMapped::VectorOfLists {
public:

    void createNew(const string& name)
    {
        toc.createNew(name + ".toc");
        data.createNew(name + ".data");
        freeSlots.createNew(name + ".freeSlots");
    }

    void accessExisting(const string& name, bool readWriteAccess)
    {
        toc.accessExisting(name + ".toc", readWriteAccess);
        data.accessExisting(name + ".data", readWriteAccess);
        freeSlots.accessExisting(name + ".freeSlots", readWriteAccess);
    }
    void accessExistingReadOnly(const string& name)
    {
        accessExisting(name, false);
    }

    void accessExistingReadWrite(const string& name)
    {
        accessExisting(name, true);
    }

    void close()
    {
        toc.close();
        data.close();
        freeSlots.close();
    }

    // Return the number of lists.
    size_t size() const
    {
        return toc.size();
    }

    // Return true if there are no lists.
    bool empty() const
    {
        return toc.empty();
    }

    // A list node stores an element of the list and an index offset, in the data vector,
    // to the previous and next node in the same list.
    // For each list we also store a "one past the end" Node.
    // The list is actually stored as a circular list.
    class Node {
    public:
        T t;                // User data (default constructed for the "one past the end" element of each list.
        size_t previous;    // Index in the data vector of the previous element in the list.
        size_t next;        // Index in the data vector of the next element in the list.
    };



    // An iterator stores a pointer to the VectorOfLists and the index of the item pointed to.
    // We don't store pointers so the iterator remains valid if a reallocation occurs.
    // It comes in const and non-const versions.
    class iterator {
    public:
        iterator(
            VectorOfLists<T>* container=0,
            size_t nodeIndex=std::numeric_limits<size_t>::max()) :
            container(container), nodeIndex(nodeIndex) {}
        void operator++()
        {
            nodeIndex = node().next;
        }
        void operator--()
        {
            nodeIndex = node().previous;
        }
        bool operator==(iterator that) const
        {
            return container==that.container && nodeIndex==that.nodeIndex;
        }
        bool operator!=(iterator that) const
        {
            return container!=that.container || nodeIndex!=that.nodeIndex;
        }
        T& operator*() const
        {
            return node().t;
        }
        friend class VectorOfLists<T>;
    private:
        VectorOfLists<T>* container;
        size_t nodeIndex;    // Index of the node in the data vector.
        Node& node() const
        {
            return (*container).data[nodeIndex];
        }
    };
    class const_iterator {
    public:
        const_iterator(
            const VectorOfLists<T>* container=0,
            size_t nodeIndex=std::numeric_limits<size_t>::max()) :
            container(container), nodeIndex(nodeIndex) {}
        void operator++()
        {
            nodeIndex = node().next;
        }
        void operator--()
        {
            nodeIndex = node().previous;
        }
        bool operator==(const_iterator that) const
        {
            return container==that.container && nodeIndex==that.nodeIndex;
        }
        bool operator!=(const_iterator that) const
        {
            return container!=that.container || nodeIndex!=that.nodeIndex;
        }
        const T& operator*() const
        {
            return node().t;
        }
    private:
        const VectorOfLists<T>* container;
        size_t nodeIndex;    // Index of the node in the data vector.
        const Node& node() const
        {
            return (*container).data[nodeIndex];
        }
    };




    // Return begin/end iterators for one of the lists.
    iterator begin(size_t i)
    {
        iterator it = end(i);
        ++it;
        return it;
    }
    iterator end(size_t i)
    {
        return iterator(this, toc[i]);
    }
    const_iterator begin(size_t i) const
    {
        const_iterator it = end(i);
        ++it;
        return it;
    }
    const_iterator end(size_t i) const
    {
        return const_iterator(this, toc[i]);
    }



   // Add an empty list at the end of the vector of lists.
    void push_back()
    {
        // Create the end element for this list.
        const size_t slot = allocateSlot();
        toc.push_back(slot);
        Node& node = data[slot];
        node.previous = slot;
        node.next = slot;
    }

    // Insert a T before a given iterator position.
    iterator insert(iterator it, const T& t)
    {

        // Store this t in an available slot.
        const size_t slot = allocateSlot();
        Node& node = data[slot];
        node.t = t;

        // Splice it in.
        iterator itPrevious = it;
        --itPrevious;
        itPrevious.node().next = slot;
        node.previous = itPrevious.nodeIndex;
        it.node().previous = slot;
        node.next = it.nodeIndex;

        // Return an iterator pointing to the element we just inserted.
        return iterator(this, slot);
    }

    // Insert a T at the beginning of one of the lists.
    iterator push_front(size_t i, const T& t)
    {
        return insert(begin(i), t);
    }

    // Insert a T at the end of one of the lists.
    iterator push_back(size_t i, const T& t)
    {
        return insert(end(i), t);
    }

    // Insert a T at the end of the last of the lists.
    iterator push_back(const T& t)
    {
        return insert(end(size()-1), t);
    }



    // Remove the element pointed to by an iterator.
    void erase(iterator it)
    {
        // Add this slot to the list of free slots.
        it.container->freeSlots.push_back(it.nodeIndex);

        // Make the previous and next list elements point to each other.
        iterator itPrevious = it;
        --itPrevious;
        Node& previousNode = itPrevious.node();
        iterator itNext = it;
        ++itNext;
        Node& nextNode = itNext.node();
        previousNode.next = itNext.nodeIndex;
        nextNode.previous = itPrevious.nodeIndex;
    }


    // Make each of the lists appear as its own container.
    // This only provides iteration over the list.
    // For operations that add or remove elements, the VectorOfList API
    // must be used.
    class List {
    public:
        iterator beginIterator;
        iterator endIterator;
        iterator begin() const
        {
            return beginIterator;
        }
        iterator end() const
        {
            return endIterator;
        }
        List(iterator begin, iterator end) :
            beginIterator(begin), endIterator(end)
        {
        }
    };
    List operator[] (size_t i)
    {
        return List(begin(i), end(i));
    }
    class ConstList {
    public:
        const_iterator beginIterator;
        const_iterator endIterator;
        const_iterator begin() const
        {
            return beginIterator;
        }
        const_iterator end() const
        {
            return endIterator;
        }
        ConstList(const_iterator begin, const_iterator end) :
            beginIterator(begin), endIterator(end)
        {
        }
    };
    ConstList operator[] (size_t i) const
    {
        return ConstList(begin(i), end(i));
    }



    // Dump all internal data structures.
    void dump(ostream& s) const
    {
        s << "toc" << "\n";
        for(size_t i=0; i<toc.size(); i++) {
            s << i << ": " << toc[i] << "\n";
        }
        s << "data" << "\n";
        for(size_t i=0; i<data.size(); i++) {
            s << i << ":";
            s << ", previous " << data[i].previous;
            s << ", next " << data[i].next << "\n";
        }
        s << "freeSlots" << "\n";
        for(size_t i=0; i<freeSlots.size(); i++) {
            s << i << ": " << freeSlots[i] << "\n";
        }
        s << flush;

    }



private:

    // Indexes in the data vector end node of each list.
    Vector<size_t> toc;

    // The list nodes.
    Vector<Node> data;

    // Offsets in the data vector to elements that are currently not used.
    Vector<size_t> freeSlots;



    // Allocate a slot to be added to an existing list.
    // Use a free slot if we have one. Otherwise, increase by one the
    // size of the data vector.
    size_t allocateSlot()
    {
        // If we have any free slots, use one of them to satisfy the request.
        if(!freeSlots.empty()) {
            const size_t slot = freeSlots.back();
            freeSlots.resize(freeSlots.size() - 1);
            return slot;
        }

        // Otherwise, increase by one the size of the data vector.
        size_t slot = data.size();
        data.resize(slot + 1);
        return slot;
    }

};




// Unit test.
inline void ChanZuckerberg::ExpressionMatrix2::testMemoryMappedVectorOfLists()
{
    // Write a message to begin.
    cout << "testVectorOfLists begins." << endl;

    // Create the vector of lists.
    MemoryMapped::VectorOfLists<int> v;
    v.createNew("VectorOfLists");
    CZI_ASSERT(v.empty());

    // Add three empty lists.
    v.push_back();
    v.push_back();
    v.push_back();
    CZI_ASSERT(!v.empty());
    CZI_ASSERT(v.size() == 3);

    // Add some elements to one of the lists.
    v.push_back(1,10);
    const auto it = v.push_back(1,20);
    v.push_back(1,30);
    v.push_back(1,40);
    v.erase(it);
    v.push_back(2, 50);



    // Dump the entire vector of lists.
    cout << "Container contents:" << endl;
    for(size_t i=0; i<v.size(); i++) {
        cout << i;
        for(auto it=v.begin(i); it!=v.end(i); ++it) {
            cout << " " << *it;
        }
        cout << endl;
    }

    cout << "Iteration over v[1]:";
    for(int i: v[1]) {
        cout << " " << i;
    }
    cout << endl;


    // v.dump(cout);

    // Done.
    cout << "testVectorOfLists ends." << endl;
}


#endif


