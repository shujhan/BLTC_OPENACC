cmake_minimum_required(VERSION 3.9...3.14)
project(FARRSIGHT LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the C++ compiler (e.g., NVIDIA's pgc++)
set(CMAKE_CXX_COMPILER nvc++)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "None")
endif()

set(default_include /usr/local/include)
if(NOT INCLUDE_HEADERS)
    message(STATUS "Setting include directory to '${default_include}' as none was specified.")
    set(INCLUDE_HEADERS ${default_include} CACHE FILEPATH "Choose include directory")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

option(ENABLE_GPU_BUILD "Build OpenACC GPU version" ON)

# Enable OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Interface for OpenACC
add_library(OpenACC_CXX INTERFACE)
# set_property(TARGET OpenACC_CXX PROPERTY INTERFACE_COMPILE_OPTIONS -acc -Minfo=accel)

# Add the source files
set(SOURCES
    toy_code.cpp
    FieldStructure.cpp
)

# Add CPU executable target
add_executable(farrsight_cpu ${SOURCES})
target_link_libraries(farrsight_cpu PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(farrsight_cpu PUBLIC ${INCLUDE_HEADERS} inc)

# Install CPU executable
install(TARGETS farrsight_cpu RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

# Set additional flags for NVIDIA compilers if ENABLE_GPU_BUILD is on
if(ENABLE_GPU_BUILD AND CMAKE_CXX_COMPILER_ID MATCHES NVHPC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -acc -Minfo=accel -gpu=managed") # -gpu=managed -ta=tesla:managed -gpu=nordc
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -acc")
else()
    message(STATUS "GPU build requires PGI compilers! Turning GPU build off.")
    set(ENABLE_GPU_BUILD OFF CACHE BOOL "Build OpenACC GPU version" FORCE)
endif()

if(ENABLE_GPU_BUILD)
    # GPU executable target with OpenACC enabled
    add_executable(farrsight_gpu ${SOURCES})
    
    # # Apply OpenACC compile and link options
    # target_compile_options(farrsight_gpu PUBLIC -fast -acc -Minfo=accel) # -gpu=managed
    # target_link_options(farrsight_gpu PUBLIC -acc)

    # Link with OpenMP and OpenACC libraries
    target_link_libraries(farrsight_gpu PUBLIC OpenMP::OpenMP_CXX OpenACC_CXX)
    target_include_directories(farrsight_gpu PUBLIC ${INCLUDE_HEADERS} inc)

    # Define OPENACC_ENABLED for conditional compilation
    target_compile_definitions(farrsight_gpu PUBLIC OPENACC_ENABLED)

    # Install GPU executable
    install(TARGETS farrsight_gpu RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
endif()
