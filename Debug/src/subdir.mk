################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CellSets.cpp \
../src/CellSimilarityGraph.cpp \
../src/ExpressionMatrix.cpp \
../src/ExpressionMatrixFindSimilarPairs.cpp \
../src/ExpressionMatrixHttpServer.cpp \
../src/ExpressionMatrixLsh.cpp \
../src/HttpServer.cpp \
../src/MemoryMappedStringTable.cpp \
../src/MurmurHash2.cpp \
../src/PythonModule.cpp \
../src/SimilarPairs.cpp \
../src/color.cpp \
../src/touchMemory.cpp 

OBJS += \
./src/CellSets.o \
./src/CellSimilarityGraph.o \
./src/ExpressionMatrix.o \
./src/ExpressionMatrixFindSimilarPairs.o \
./src/ExpressionMatrixHttpServer.o \
./src/ExpressionMatrixLsh.o \
./src/HttpServer.o \
./src/MemoryMappedStringTable.o \
./src/MurmurHash2.o \
./src/PythonModule.o \
./src/SimilarPairs.o \
./src/color.o \
./src/touchMemory.o 

CPP_DEPS += \
./src/CellSets.d \
./src/CellSimilarityGraph.d \
./src/ExpressionMatrix.d \
./src/ExpressionMatrixFindSimilarPairs.d \
./src/ExpressionMatrixHttpServer.d \
./src/ExpressionMatrixLsh.d \
./src/HttpServer.d \
./src/MemoryMappedStringTable.d \
./src/MurmurHash2.d \
./src/PythonModule.d \
./src/SimilarPairs.d \
./src/color.d \
./src/touchMemory.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -I/usr/include/python3.5m -O0 -g3 -Wall -Wconversion -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


