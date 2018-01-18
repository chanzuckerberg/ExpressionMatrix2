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
        cout << "OpenCL error " << e.err() <<  " from " << e.what();
        cout << " during GPU initialization." << endl;
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
void Lsh::setupGpuKernel0()
{
    CZI_ASSERT(gpu.isInitialized);
    gpu.kernel0MismatchBuffer = std::make_shared<cl::Buffer>(
        gpu.context, CL_MEM_READ_WRITE, cellCount() * sizeof(uint16_t));

}
void Lsh::gpuKernel0(
    CellId cellId0,
    vector<uint16_t>& mismatchCounts)
{

    try {
        // Resize the output vector.
        mismatchCounts.resize(cellCount());

        // Run the kernel.
        gpu.kernel0.setArg(0, cellCount());
        gpu.kernel0.setArg(1, wordCount());
        gpu.kernel0.setArg(2, gpu.signatureBuffer);
        gpu.kernel0.setArg(3, cellId0);
        gpu.kernel0.setArg(4, *gpu.kernel0MismatchBuffer);
        gpu.queue.enqueueNDRangeKernel(gpu.kernel0,
            cl::NullRange, cl::NDRange(cellCount()), cl::NullRange);

        // Get back the results.
        gpu.queue.enqueueReadBuffer(*(gpu.kernel0MismatchBuffer), CL_TRUE, 0,
            cellCount() * sizeof(uint16_t), mismatchCounts.data());
        gpu.queue.finish();

    } catch(cl::Error e) {
        const string message = "OpenCL error " + std::to_string(e.err()) + " from " +
            e.what() + " running GPU kernel 0.";
        throw runtime_error(message);
    }

}
void Lsh::cleanupGpuKernel0()
{
    gpu.kernel0MismatchBuffer.reset();
}



#endif

