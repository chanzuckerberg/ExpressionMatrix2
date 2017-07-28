################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CellSets.cpp \
../src/CellSimilarityGraph.cpp \
../src/ClusterGraph.cpp \
../src/ExpressionMatrix.cpp \
../src/ExpressionMatrixBioHub.cpp \
../src/ExpressionMatrixFindSimilarPairs.cpp \
../src/ExpressionMatrixHdf5.cpp \
../src/ExpressionMatrixHttpServer.cpp \
../src/ExpressionMatrixLsh.cpp \
../src/GeneSet.cpp \
../src/HttpServer.cpp \
../src/MemoryMappedStringTable.cpp \
../src/MurmurHash2.cpp \
../src/PythonModule.cpp \
../src/SimilarPairs.cpp \
../src/color.cpp \
../src/hdf5.cpp \
../src/regressionCoefficient.cpp \
../src/tokenize.cpp \
../src/touchMemory.cpp 

OBJS += \
./src/CellSets.o \
./src/CellSimilarityGraph.o \
./src/ClusterGraph.o \
./src/ExpressionMatrix.o \
./src/ExpressionMatrixBioHub.o \
./src/ExpressionMatrixFindSimilarPairs.o \
./src/ExpressionMatrixHdf5.o \
./src/ExpressionMatrixHttpServer.o \
./src/ExpressionMatrixLsh.o \
./src/GeneSet.o \
./src/HttpServer.o \
./src/MemoryMappedStringTable.o \
./src/MurmurHash2.o \
./src/PythonModule.o \
./src/SimilarPairs.o \
./src/color.o \
./src/hdf5.o \
./src/regressionCoefficient.o \
./src/tokenize.o \
./src/touchMemory.o 

CPP_DEPS += \
./src/CellSets.d \
./src/CellSimilarityGraph.d \
./src/ClusterGraph.d \
./src/ExpressionMatrix.d \
./src/ExpressionMatrixBioHub.d \
./src/ExpressionMatrixFindSimilarPairs.d \
./src/ExpressionMatrixHdf5.d \
./src/ExpressionMatrixHttpServer.d \
./src/ExpressionMatrixLsh.d \
./src/GeneSet.d \
./src/HttpServer.d \
./src/MemoryMappedStringTable.d \
./src/MurmurHash2.d \
./src/PythonModule.d \
./src/SimilarPairs.d \
./src/color.d \
./src/hdf5.d \
./src/regressionCoefficient.d \
./src/tokenize.d \
./src/touchMemory.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -I/usr/include/python3.5m -I/usr/include/hdf5/serial -O3 -msse4.2 -ggdb3 -Wall -Wconversion -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


