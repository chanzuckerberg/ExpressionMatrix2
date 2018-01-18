// This file defines OpenCL kernel code
// used by class Lsh (see Lsh.hpp and LshGpu.cpp).

const string code = R"###(



// Compute the number of mismatches between the signature
// of cell cellId0 and the signatures of all other cells.
// This is parallelized over cellId1.
// This is a simple but working kernel.
// It has high overhead and low performance
// because of the small amount of work done by each instance.
void kernel kernel0(
    uint cellCount,
    ulong lshWordCount,
    global const ulong* signatures, 
    uint cellId0, 
    global ushort* mismatchCounts)
{ 
    uint cellId1 = get_global_id(0);
    global const ulong* signature0 = signatures + cellId0 * lshWordCount;
    global const ulong* signature1 = signatures + cellId1 * lshWordCount;
    ulong mismatchCount = 0;
    uint k;
    for(k=0; k<lshWordCount; k++) {
        mismatchCount += popcount((*signature0++) ^ (*signature1++));
    }
    mismatchCounts[cellId1] = mismatchCount;
}



)###";
