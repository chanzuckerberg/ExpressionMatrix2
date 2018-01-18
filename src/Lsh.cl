// This file defines OpenCL kernel code
// used by class Lsh (see Lsh.hpp and LshGpu.cpp).

const string code = R"###(



// Kernel 0: compute the number of mismatches between the signature
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


// Kernel 1: compute the number of mismatches between the signature
// of cells in range [cellId0Begin, cellId0End]
// and the signatures of all other cells.
// This is parallelized over cellId1.
// This has lower overhead than kernel 0
// because each kernel instance does more work
// (processes n0 values of cellId0 instead of just one).
// For efficient memory access, the blockSize mismatch counts 
// for each cellId1 are stored contiguously.
void kernel kernel1(
    uint cellCount,
    ulong lshWordCount,
    global const ulong* signatures, 
    uint cellId0Begin, 
    uint cellId0End,
    uint blockSize,
    global ushort* mismatchCounts)
{ 
    uint cellId1 = get_global_id(0);
    global const ulong* signature0Begin = signatures + cellId0Begin * lshWordCount;
    global const ulong* signature1Begin = signatures + cellId1 * lshWordCount;
    global ushort* mismatchCount1 = mismatchCounts + blockSize*cellId1;
    
    // Loop over the values of cellId0 assigned to this instance of the kernel.
    uint cellId0;
    for(cellId0=cellId0Begin; cellId0!=cellId0End; ++cellId0, signature0Begin+=lshWordCount) {
    
        // Compute the number of mismatches between the signatures
        // of cellId0 and cellId1.
        global const ulong* signature0 = signature0Begin;
        global const ulong* signature1 = signature1Begin;
        ulong mismatchCount = 0;
        uint k;
        for(k=0; k<lshWordCount; k++) {
            mismatchCount += popcount((*signature0++) ^ (*signature1++));
        }
        *mismatchCount1++ = mismatchCount;
    }
}

)###";
