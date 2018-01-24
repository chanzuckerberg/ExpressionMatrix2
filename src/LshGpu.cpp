// This files contains the implementation
// of portions of class Lsh which use OpenCL
// to perform computations on GPUs.

#if CZI_EXPRESSION_MATRIX2_BUILD_FOR_GPU


#include "Lsh.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;



void Lsh::initializeGpu()
{
    gpu.initialize();
}


// Initialize the OpenCL platform, device, context, and queue.
void Lsh::Gpu::initialize()
{
    // If already initialized, do nothing.
    if(isInitialized) {
        return;
    }

    try {
        choosePlatform();
        chooseDevice();
        context = cl::Context(device);
        queue = cl::CommandQueue(context, device);
        buildProgram();

        // Mark it as initialized.
        isInitialized = true;

    } catch(cl::Error e) {
        const string message =
            "OpenCL error " + std::to_string(e.err()) +
            " from " + e.what() + " during GPU initialization.";
        throw runtime_error(message);
    }
}


void Lsh::Gpu::choosePlatform()
{
    // Get the platforms supported on the current system.
     vector<cl::Platform> platforms;
     cl::Platform::get(&platforms);
     if(platforms.empty()){
         throw runtime_error("No OpenCL plaform found.");
     }

     // Use the first platform we found.
     platform = platforms.front();

}

void Lsh::Gpu::chooseDevice()
{
    // Get the devices of this platform.
    vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
    if(devices.empty()){
        throw runtime_error("No OpenCL device found.");
    }

    // Use the first device we found.
    device = devices.front();
}


void Lsh::Gpu::writeDeviceInformation(ostream& s) const
{
    CZI_ASSERT(isInitialized);

    s << "CL_DEVICE_NAME "<< device.getInfo<CL_DEVICE_NAME>() << endl;
    s << "CL_DEVICE_VENDOR "<< device.getInfo<CL_DEVICE_VENDOR>() << endl;
    s << "CL_DRIVER_VERSION "<< device.getInfo<CL_DRIVER_VERSION>() << endl;
    s << "CL_DEVICE_PROFILE "<< device.getInfo<CL_DEVICE_PROFILE>() << endl;
    s << "CL_DEVICE_VERSION "<< device.getInfo<CL_DEVICE_VERSION>() << endl;
    s << "CL_DEVICE_EXTENSIONS "<< device.getInfo<CL_DEVICE_EXTENSIONS>() << endl;
    s << "CL_DEVICE_MAX_COMPUTE_UNITS " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << endl;
    s << "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS " << device.getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>() << endl;
    s << "CL_DEVICE_MAX_WORK_GROUP_SIZE " << device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << endl;
    s << "CL_DEVICE_MAX_CLOCK_FREQUENCY " << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << endl;
    s << "CL_DEVICE_ADDRESS_BITS " << device.getInfo<CL_DEVICE_ADDRESS_BITS>() << endl;
    s << "CL_DEVICE_MAX_MEM_ALLOC_SIZE " << device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << endl;

    const auto maxWorkItemSizes = device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
    s << "CL_DEVICE_MAX_WORK_ITEM_SIZES ";
    for(const auto maxWorkItemSize: maxWorkItemSizes) {
        s << maxWorkItemSize << " ";
    }
    s << endl;
}


string Lsh::getGpuName() const
{
    return gpu.deviceName();
}


string Lsh::Gpu::deviceName() const
{
    CZI_ASSERT(isInitialized);
    return device.getInfo<CL_DEVICE_NAME>();
}



void Lsh::Gpu::buildProgram()
{
    // Define a C++ string named "code" containing the OpenCL code.
    #include "Lsh.cl"

    // Build it.
    cl::Program::Sources sources = {code};
    program = cl::Program(context, sources);
    try {
        program.build({device});
    } catch(cl::BuildError) {
        string buildLog;
        program.getBuildInfo(device, CL_PROGRAM_BUILD_LOG, &buildLog);
        cout << "OpenCL build failed:" << endl;
        cout << buildLog << endl;
        throw;
    }

    // Store our kernels.
    kernel0 = cl::Kernel(program, "kernel0");
    kernel1 = cl::Kernel(program, "kernel1");
    kernel2 = cl::Kernel(program, "kernel2");
    kernel3 = cl::Kernel(program, "kernel3");

}



// Load the LSH signatures to the GPU.
void Lsh::loadSignaturesToGpu()
{
    CZI_ASSERT(gpu.isInitialized);

    const uint64_t signatureBufferSize = cellCount()*wordCount()*sizeof(uint64_t);
    gpu.signatureBuffer =
        cl::Buffer(gpu.context, CL_MEM_READ_WRITE, signatureBufferSize);

    gpu.queue.enqueueWriteBuffer(
        gpu.signatureBuffer, CL_TRUE, 0,
        signatureBufferSize,
        signatures.begin());
    gpu.queue.finish();

}


// Kernel 0: compute the number of mismatches between the signature
// of cell cellId0 and the signatures of all other cells.
// This is a simple but working kernel.
// It has high overhead and low performance
// because of the small amount of work done by each instance.
void Lsh::setupGpuKernel0(vector<uint16_t>& mismatchCounts)
{
    CZI_ASSERT(gpu.isInitialized);

    // Set up the buffer to hold the computed number of mismatches
    // for each cell.
    mismatchCounts.resize(cellCount());
    gpu.mismatchBufferSize = cellCount() * sizeof(uint16_t);
    gpu.mismatchBufferHostPointer = mismatchCounts.data();
    gpu.mismatchBuffer = std::make_shared<cl::Buffer>(
        gpu.context, CL_MEM_READ_WRITE, gpu.mismatchBufferSize);

    // Set the arguments that don't change.
    gpu.kernel0.setArg(0, cellCount());
    gpu.kernel0.setArg(1, wordCount());
    gpu.kernel0.setArg(2, gpu.signatureBuffer);
    gpu.kernel0.setArg(4, *gpu.mismatchBuffer);

}
void Lsh::gpuKernel0(CellId cellId0)
{

    try {

        // Run the kernel.
        gpu.kernel0.setArg(3, cellId0);
        gpu.queue.enqueueNDRangeKernel(gpu.kernel0,
            cl::NullRange, cl::NDRange(cellCount()), cl::NullRange);

        // Get back the results.
        gpu.queue.enqueueReadBuffer(*(gpu.mismatchBuffer), CL_TRUE, 0,
            gpu.mismatchBufferSize, gpu.mismatchBufferHostPointer);
        gpu.queue.finish();

    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " running GPU kernel 0.";
        throw runtime_error(message);
    }

}
void Lsh::cleanupGpuKernel0()
{
    gpu.mismatchBuffer.reset();
}



// Kernel 1: compute the number of mismatches between the signature
// of cells in range [cellId0Begin, cellId0End]
// and the signatures of all other cells.
// This is parallelized over cellId1.
// This has lower overhead than kernel 0
// because each kernel instance does more work
// (processes n0 values of cellId0 instead of just one).
// For efficient memory access, the blockSize mismatch counts
// for each cellId1 are stored contiguously, with stride blockSize
// betwene successive values of cellId1.
void Lsh::setupGpuKernel1(vector<uint16_t>& mismatchCounts, CellId blockSize)
{
    CZI_ASSERT(gpu.isInitialized);

    // Set up the buffer to hold the computed number of mismatches
    // for each cell.
    mismatchCounts.resize(cellCount() * blockSize);
    gpu.mismatchBufferSize = cellCount() * blockSize * sizeof(uint16_t);
    gpu.mismatchBufferHostPointer = mismatchCounts.data();
    gpu.mismatchBuffer = std::make_shared<cl::Buffer>(
        gpu.context, CL_MEM_READ_WRITE, gpu.mismatchBufferSize);

    // Set the arguments that don't change.
    gpu.kernel1.setArg(0, cellCount());
    gpu.kernel1.setArg(1, wordCount());
    gpu.kernel1.setArg(2, gpu.signatureBuffer);
    gpu.kernel1.setArg(5, blockSize);
    gpu.kernel1.setArg(6, *gpu.mismatchBuffer);

}
void Lsh::gpuKernel1(CellId cellId0Begin, CellId cellId0End)
{

    try {

        // Run the kernel.
        gpu.kernel1.setArg(3, cellId0Begin);
        gpu.kernel1.setArg(4, cellId0End);
        gpu.queue.enqueueNDRangeKernel(gpu.kernel1,
            cl::NullRange, cl::NDRange(cellCount()), cl::NullRange);

        // Get back the results.
        gpu.queue.enqueueReadBuffer(*(gpu.mismatchBuffer), CL_TRUE, 0,
            gpu.mismatchBufferSize, gpu.mismatchBufferHostPointer);
        gpu.queue.finish();

    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " running GPU kernel 1.";
        throw runtime_error(message);
    }

}
void Lsh::cleanupGpuKernel1()
{
    gpu.mismatchBuffer.reset();
}



// Set up the neighbors buffer on the host and on the GPU.
void Lsh::Gpu::setupNeighborsBuffer(
    vector<uint32_t>& neighbors,
    size_t blockSize,
    size_t k,
    size_t mismatchCountThreshold)
{
    try {
        neighbors.resize(k * blockSize * (mismatchCountThreshold+1));
        neighborsBufferSize = neighbors.size() * sizeof(uint32_t);
        // cout << "Neighbors buffer size is " << neighborsBufferSize << endl;
        neighborsBufferHostPointer = neighbors.data();
        neighborsBuffer = std::make_shared<cl::Buffer>(
            context, CL_MEM_READ_WRITE, neighborsBufferSize);
    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " setting up neighbors buffer.";
        throw runtime_error(message);
    }

}

// Set up the neighborCounts buffer on the host and on the GPU.
void Lsh::Gpu::setupNeighborCountsBuffer (
    vector<uint32_t>& neighborCounts,
    size_t blockSize,
    size_t mismatchCountThreshold)
{
    try {
        neighborCounts.resize(blockSize * (mismatchCountThreshold+1));
        neighborCountsBufferSize = neighborCounts.size() * sizeof(uint32_t);
        // cout << "Neighbor counts buffer size is " << neighborCountsBufferSize << endl;
        neighborCountsBufferHostPointer = neighborCounts.data();
        neighborCountsBuffer = std::make_shared<cl::Buffer>(
            context, CL_MEM_READ_WRITE, neighborCountsBufferSize);
    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " setting up neighbor counts buffer.";
        throw runtime_error(message);
    }
}



// Kernel2. See Lsh.cl for details.
void Lsh::setupGpuKernel2(
    vector<uint32_t>& neighbors,
    vector<uint32_t>& neighborCounts,
    CellId blockSize,
    size_t k,
    size_t mismatchCountThreshold)
{
    try {
        gpu.setupNeighborsBuffer(neighbors, blockSize, k, mismatchCountThreshold);
        gpu.setupNeighborCountsBuffer(neighborCounts, blockSize, mismatchCountThreshold);

        // Set the arguments that don't change.
        gpu.kernel2.setArg(0, uint64_t(cellCount()));
        gpu.kernel2.setArg(1, wordCount());
        gpu.kernel2.setArg(3, k);
        gpu.kernel2.setArg(4, mismatchCountThreshold);
        gpu.kernel2.setArg(5, gpu.signatureBuffer);
        gpu.kernel2.setArg(6, *gpu.neighborsBuffer);
        gpu.kernel2.setArg(7, *gpu.neighborCountsBuffer);

    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " setting up GPU kernel 2.";
        throw runtime_error(message);
    }
}



void Lsh::gpuKernel2(CellId cellId0Begin, CellId cellId0End)
{
    try {

        // Run the kernel.
        gpu.kernel2.setArg(2, uint64_t(cellId0Begin));
        gpu.queue.enqueueNDRangeKernel(gpu.kernel2,
            cl::NullRange, cl::NDRange(cellId0End-cellId0Begin), cl::NullRange);

        // Get back the results.
        gpu.queue.enqueueReadBuffer(*(gpu.neighborsBuffer), CL_TRUE, 0,
            gpu.neighborsBufferSize, gpu.neighborsBufferHostPointer);
        gpu.queue.enqueueReadBuffer(*(gpu.neighborCountsBuffer), CL_TRUE, 0,
            gpu.neighborCountsBufferSize, gpu.neighborCountsBufferHostPointer);
        gpu.queue.finish();

    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " running GPU kernel 2.";
        throw runtime_error(message);
    }
}



void Lsh::cleanupGpuKernel2()
{
    gpu.neighborsBuffer.reset();
    gpu.neighborCountsBuffer.reset();
}



// Set up the cellIds1 buffer on the host and on the GPU.
void Lsh::Gpu::setupCellId1sBuffer(
    vector<uint32_t>& cellId1s,
    size_t blockSize,
    size_t maxCheck
    )
{
    try {
        cellId1s.resize(blockSize * maxCheck);
        cellId1sBufferSize = cellId1s.size() * sizeof(uint32_t);
        // cout << "cellId1s buffer size is " << cellId1sBufferSize << endl;
        cellId1sBufferHostPointer = cellId1s.data();
        cellId1sBuffer = std::make_shared<cl::Buffer>(
            context, CL_MEM_READ_WRITE, cellId1sBufferSize);
    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " setting up cellIds1 buffer.";
        throw runtime_error(message);
    }

}


// Set up the cellId1Counts buffer on the host and on the GPU.
void Lsh::Gpu::setupCellId1CountsBuffer(
    vector<uint32_t>& cellId1Counts,
    size_t blockSize
    )
{
    try {
        cellId1Counts.resize(blockSize);
        cellId1CountsBufferSize = cellId1Counts.size() * sizeof(uint32_t);
        // cout << "cellId1Counts buffer size is " << cellId1CountsBufferSize << endl;
        cellId1CountsBufferHostPointer = cellId1Counts.data();
        cellId1CountsBuffer = std::make_shared<cl::Buffer>(
            context, CL_MEM_READ_WRITE, cellId1CountsBufferSize);
    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " setting up cellId1Counts buffer.";
        throw runtime_error(message);
    }

}


// Kernel3. See Lsh.cl for details.
void Lsh::setupGpuKernel3(
    size_t k,
    size_t mismatchCountThreshold,
    size_t maxCheck,
    size_t blockSize,
    vector<CellId>& cellIds1,
    vector<CellId>& cellId1Counts,
    vector<CellId>& neighbors,
    vector<CellId>& neighborCounts)
{
    try {
        gpu.setupNeighborsBuffer(neighbors, blockSize, k, mismatchCountThreshold);
        gpu.setupNeighborCountsBuffer(neighborCounts, blockSize, mismatchCountThreshold);
        gpu.setupCellId1sBuffer(cellIds1, blockSize, maxCheck);
        gpu.setupCellId1CountsBuffer(cellId1Counts, blockSize);

        // Set the arguments that don't change.
        gpu.kernel3.setArg(0, wordCount());
        // gpu.kernel3.setArg(1, cellBegin0);   // This is set differently for each block, by gpuGpuKernel3.
        gpu.kernel3.setArg(2, k);
        gpu.kernel3.setArg(3, mismatchCountThreshold);
        gpu.kernel3.setArg(4, maxCheck);
        gpu.kernel3.setArg(5, gpu.signatureBuffer);
        gpu.kernel3.setArg(6, *gpu.cellId1CountsBuffer);
        gpu.kernel3.setArg(7, *gpu.cellId1sBuffer);
        gpu.kernel3.setArg(8, *gpu.neighborCountsBuffer);
        gpu.kernel3.setArg(9, *gpu.neighborsBuffer);

    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " setting up GPU kernel 2.";
        throw runtime_error(message);
    }
}



// Kernel3. See Lsh.cl for details.
void Lsh::gpuKernel3(CellId cellId0Begin, CellId cellId0End)
{
    // cout << "gpuKernel3 begins "<< cellId0Begin << " " << cellId0End << endl;
    try {

        // Copy the cellId1s and cellId1Counts to the GPU.
        gpu.queue.enqueueWriteBuffer(*(gpu.cellId1CountsBuffer), CL_TRUE, 0,
            gpu.cellId1CountsBufferSize, gpu.cellId1CountsBufferHostPointer);
        gpu.queue.enqueueWriteBuffer(*(gpu.cellId1sBuffer), CL_TRUE, 0,
            gpu.cellId1sBufferSize, gpu.cellId1sBufferHostPointer);

        // Run the kernel.
        gpu.kernel3.setArg(1, uint64_t(cellId0Begin));
        gpu.queue.enqueueNDRangeKernel(gpu.kernel3,
            cl::NullRange, cl::NDRange(cellId0End-cellId0Begin), cl::NullRange);

        // Get back the results.
        gpu.queue.enqueueReadBuffer(*(gpu.neighborsBuffer), CL_TRUE, 0,
            gpu.neighborsBufferSize, gpu.neighborsBufferHostPointer);
        gpu.queue.enqueueReadBuffer(*(gpu.neighborCountsBuffer), CL_TRUE, 0,
            gpu.neighborCountsBufferSize, gpu.neighborCountsBufferHostPointer);
        gpu.queue.finish();

    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " running GPU kernel 3.";
        throw runtime_error(message);
    }

}



void Lsh::cleanupGpuKernel3()
{
    gpu.neighborsBuffer.reset();
    gpu.neighborCountsBuffer.reset();
    gpu.cellId1sBuffer.reset();
    gpu.cellId1CountsBuffer.reset();
}

#endif

