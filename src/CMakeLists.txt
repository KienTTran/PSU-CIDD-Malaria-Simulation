# Enable if building on the Windows Subsystem for Linux (WSL), via the command line use -DBUILD_WSL:BOOL=true
set(BUILD_WSL false CACHE BOOL "Enable for builds on WSL")

include_directories(${PROJECT_SOURCE_DIR}/src ${EXT_PROJECTS_DIR}/easyloggingpp-9.97.0)

file(GLOB SRC_FILES *.cpp
        Core/*.cpp 
        Core/Config/*.cpp 
        Core/Index/*.cpp 
        Events/*.cpp
        Events/Population/*.cpp
        Helpers/*.cpp
        Helpers/*.h
        Helpers/*.hxx
        Malaria/*.cpp
        MDC/*.cpp
        Parasites/*.cpp GeneInfo.h
        Population/*.cpp
        Population/Properties/*.cpp
        Population/Properties/*.hxx
        Reporters/*.cpp
        Spatial/*.cpp
        Spatial/*.hxx
        Environment/*.cpp
        Strategies/*.cpp
        Therapies/*.cpp
        Mosquito/*.cpp
        GIS/*.cpp
        Gpu/*.cu
        Gpu/*.cuh
        Gpu/Core/*.cu
        Gpu/Core/*.cuh
        Gpu/Renderer/*.cu
        Gpu/Renderer/*.cuh
        Gpu/Population/*.cu
        Gpu/Population/*.cuh
        Gpu/Population/Properties/*.cu
        Gpu/Population/Properties/*.cuh
        Gpu/Population/Properties/*.hxx
        Gpu/Parasites/*.cu
        Gpu/Parasites/*.cuh
        Gpu/Mosquito/*.cu
        Gpu/Mosquito/*.cuh
        Gpu/Events/*.cu
        Gpu/Events/*.cuh
        Gpu/MDC/*.cu
        Gpu/MDC/*.cuh
        Gpu/Events/Population/*.cu
        Gpu/Events/Population/*.cuh
        Gpu/Strategies/*.cu
        Gpu/Strategies/*.cuh
        Gpu/Therapies/*.cu
        Gpu/Therapies/*.cuh
        Gpu/Reporters/*.cu
        Gpu/Reporters/*.cuh
        Gpu/Utils/*.cu
        Gpu/Utils/*.cuh
        Gpu/Utils/*.hpp
)

set(SRC_FILES ${SRC_FILES}
        ${EXT_PROJECTS_DIR}/easyloggingpp-9.97.0/easylogging++.h
        ${EXT_PROJECTS_DIR}/easyloggingpp-9.97.0/easylogging++.cc
        Gpu/Renderer/Plot/Plot.cu
        Gpu/Renderer/Plot/Plot.cuh
        Gpu/Renderer/Plot/PlotLines.cu
        Gpu/Renderer/Plot/PlotLines.cuh
)

add_library(MaSimCore STATIC ${SRC_FILES})

set_target_properties(
        MaSimCore
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(MaSimCore PROPERTIES CUDA_ARCHITECTURES "86")


#GLFW3 glfw3
find_package(glfw3 CONFIG REQUIRED)
target_link_libraries(MaSimCore PRIVATE glfw)

#GLEW glew
find_package(GLEW REQUIRED)
target_link_libraries(MaSimCore PRIVATE GLEW::GLEW)

#GLM glm
find_package(glm CONFIG REQUIRED)
target_link_libraries(MaSimCore PRIVATE glm::glm)

#IMGUI imgui[core,glfw-binding,opengl3-binding,docking-experimental]
find_package(imgui CONFIG REQUIRED)
target_link_libraries(MaSimCore PRIVATE imgui::imgui)

#IMPlot implot
find_package(implot CONFIG REQUIRED)
target_link_libraries(MaSimCore PRIVATE implot::implot)

#Date date
find_package(date CONFIG REQUIRED)
target_link_libraries(MaSimCore PRIVATE date::date date::date-tz)

#Yaml yaml-cpp
find_package(yaml-cpp CONFIG REQUIRED)
target_link_libraries(MaSimCore PRIVATE yaml-cpp)

#fmt fmt
find_package(fmt CONFIG REQUIRED)
target_link_libraries(MaSimCore PRIVATE fmt::fmt)

#gsl gsl
find_package(GSL REQUIRED)
target_link_libraries(MaSimCore PRIVATE GSL::gsl GSL::gslcblas)

#args
find_package(args CONFIG REQUIRED)
target_link_libraries(MaSimCore PRIVATE taywee::args)

#CLI11
find_package(CLI11 CONFIG REQUIRED)
target_link_libraries(MaSimCore PRIVATE CLI11::CLI11)

if(WIN32)
    message(STATUS "Build in Windows")
    target_include_directories(MaSimCore PUBLIC "D:/SDK/vcpkg/installed/x64-windows/include")
    target_link_directories(MaSimCore PUBLIC "D:/SDK/vcpkg/installed/x64-windows/lib")
    #CUDA - NVIDIA
    target_include_directories(MaSimCore PUBLIC "D:/SDK/CUDA/12.2.1/Toolkit/include")
    target_link_directories(MaSimCore PUBLIC "D:/SDK/CUDA/12.2.1/Toolkit/lib/x64")
    target_link_libraries(MaSimCore PRIVATE cuda)

    #OpenGL opengl
    find_package(OpenGL REQUIRED)
    target_link_libraries(MaSimCore PRIVATE opengl32)
elseif(UNIX)
    if(DOCKER)
        message(STATUS "Build in Docker")
        target_include_directories(MaSimCore PUBLIC "/root/vcpkg/installed/x64-linux/include")
        target_link_directories(MaSimCore PUBLIC "/root/vcpkg/installed/x64-linux/lib")
    elseif(WSL)
        message(STATUS "Build in WSL")
        target_include_directories(MaSimCore PUBLIC "/mnt/d/SDK/vcpkg/installed/x64-linux/include")
        target_link_directories(MaSimCore PUBLIC "/mnt/d/SDK/vcpkg/installed/x64-linux/lib")
    else()
        message(STATUS "Build in Linux")
        target_include_directories(MaSimCore PUBLIC "/home/ktt/SDK/vcpkg/installed/x64-linux/include")
        target_link_directories(MaSimCore PUBLIC "/home/ktt/SDK/vcpkg/installed/x64-linux/lib")
    endif()
    set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
    include_directories("/usr/local/cuda/include")
    include_directories("/usr/local/include")
    include_directories("/usr/include")
    #OpenGL opengl
    find_package(OpenGL REQUIRED)
    target_link_libraries(MaSimCore PRIVATE OpenGL::GL)
endif()

add_executable(MaSim MaSim/main.cpp)
add_dependencies(MaSim MaSimCore)

target_link_libraries(MaSim PRIVATE MaSimCore)
target_compile_features(MaSim PRIVATE cxx_range_for)
set_target_properties(MaSim PROPERTIES CUDA_ARCHITECTURES "86")

install(TARGETS MaSim DESTINATION ${PROJECT_SOURCE_DIR}/bin)
install(FILES ${PROJECT_SOURCE_DIR}/input/input.yml DESTINATION ${PROJECT_SOURCE_DIR}/bin)

add_subdirectory(DxGGenerator)
add_subdirectory(EfficacyEstimator)

#install(TARGETS DxGGenerator DESTINATION ${PROJECT_SOURCE_DIR}/bin)
#install(FILES ${PROJECT_SOURCE_DIR}/misc/input_DxG.yml DESTINATION ${PROJECT_SOURCE_DIR}/bin)
