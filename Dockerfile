# Use the NVIDIA CUDA base image
FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

# Set the working directory
WORKDIR /root

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies for vcpkg
RUN apt-get update && \
    apt-get install -y \
        git \
        build-essential \
		curl \
		zip \
		gdb \
		cmake \
        vim

RUN apt install -y libxmu-dev libxi-dev libgl-dev
RUN apt install -y libxinerama-dev libxcursor-dev xorg-dev libglu1-mesa-dev pkg-config

# Clone vcpkg repository
RUN git clone https://github.com/microsoft/vcpkg.git

# Bootstrap and integrate vcpkg
RUN cd vcpkg && \
    ./bootstrap-vcpkg.sh && \
	./vcpkg install glfw3 opengl glew glm imgui[core,glfw-binding,opengl3-binding] \
	easyloggingpp date yaml-cpp fmt gsl args CLI11 gtest && \
    ./vcpkg integrate install

# Add vcpkg to the PATH
ENV PATH="/root/vcpkg:${PATH}"

# Create a 'build' directory within MaSimGPU
RUN mkdir /root/MaSimGPU
RUN mkdir /root/MaSimGPU/build
RUN mkdir /root/MaSimGPU/build/output

# Copy the MaSimGPU folder into the Docker image
COPY . /root/MaSimGPU

# Set the working directory to the 'build' directory
WORKDIR /root/MaSimGPU/build

# Run CMake to build the project
RUN cmake \
    -DCMAKE_BUILD_TYPE=Release \
	-DDOCKER=true \
	-DCMAKE_TOOLCHAIN_FILE="/root/vcpkg/scripts/buildsystems/vcpkg.cmake" \
	.. && \
    make MaSim -j8

# Cleanup unnecessary packages
#RUN apt-get autoremove -y && \
#    apt-get clean && \
#    rm -rf /var/lib/apt/lists/*

WORKDIR /root/MaSimGPU/build/bin

# Set the default command to run when the container starts
CMD ["./MaSim", "-i", "/root/MaSimGPU/input/input_pfcrt_avg_crs_raster_gpu_dev.yml", "-r", "ValidationReporter", "-o", "/root/MaSimGPU/build/output/", "--v=1"]

#build cmd
#docker build -t ktt878/masim_gpu .
#docker build -t ktt878/masim_gpu .

#run cmd
#docker pull ktt878/masim_gpu
#docker run --gpus all user/image:tag
#docker run --gpus all user/image:tag bash
#docker exec -it container_id nvidia-smi
#sudo docker run --rm --runtime=nvidia --gpus all -it ktt878/masim_gpu bash
#sudo docker run --rm --runtime=nvidia --gpus all ktt878/masim_gpu