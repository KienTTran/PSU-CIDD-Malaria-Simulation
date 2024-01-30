cmake_minimum_required(VERSION 3.12)

set(CMAKE_TOOLCHAIN_FILE
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/vcpkg/scripts/buildsystems/vcpkg.cmake CACHE
  STRING "Vcpkg toolchain file")

if(WIN32)
  message(STATUS "Build in Windows")
  set(CMAKE_CUDA_COMPILER "D:/SDK/CUDA/12.2.1/Toolkit/bin/nvcc.exe")
elseif(UNIX)
  if(DOCKER)
    message(STATUS "Build in Docker")
    set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
  else()
    message(STATUS "Build in Linux")
    set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
  endif()
endif()

if (MSVC)
  add_definitions(-Dgtest_force_shared_crt=ON)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    MESSAGE(STATUS "Test Build in Windows Debug")
    add_definitions(-D_ITERATOR_DEBUG_LEVEL=2)
  else()
    MESSAGE(STATUS "Test Build in Windows Release")
    add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)
  endif()
endif ()

project(POMS LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

if (MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNOMINMAX -D__DISABLE_CRIT_ERR=1")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DNOMINMAX -D__DISABLE_CRIT_ERR=1")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-D_ITERATOR_DEBUG_LEVEL=2)
  else()
    add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)
  endif()
  add_definitions(-Dgtest_force_shared_crt=ON)
endif ()
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --generate-line-info")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DELPP_NO_DEFAULT_LOG_FILE")

# the next 2 flags is used for profiler with perf
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

IF (APPLE)
  # Fix linking on 10.14+. See https://stackoverflow.com/questions/54068035
  LINK_DIRECTORIES(/usr/local/lib)
  INCLUDE_DIRECTORIES(/usr/local/include)
ENDIF ()

set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src)

#Setup test project
option(BUILD_TESTS "Build all tests." ON) # Makes boolean 'test' available.
if (BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif ()

#enable_testing()
#add_test(test1 ${CMAKE_BINARY_DIR}/bin/${PROJECT_TEST_NAME} --gtest_color=yes)

