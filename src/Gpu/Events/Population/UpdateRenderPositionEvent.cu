//
// Created by kient on 8/4/2023.
//

#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/remove.h>
#include "easylogging++.h"

#include "Gpu/Core/Scheduler.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/Population/Properties/PersonIndexGPU.cuh"

#include "Gpu/Core/Random.cuh"
#include "Gpu/Renderer/RenderEntity.cuh"
#include "Gpu/Utils/Utils.cuh"
#include "Gpu/Population/Population.cuh"

#include "UpdateRenderPositionEvent.cuh"

GPU::UpdateRenderPositionEvent::UpdateRenderPositionEvent() = default;

GPU::UpdateRenderPositionEvent::~UpdateRenderPositionEvent() = default;

void GPU::UpdateRenderPositionEvent::schedule_event(GPU::Scheduler *scheduler, const int &time) {
    if (scheduler!=nullptr) {
        auto *person_update_event = new UpdateRenderPositionEvent();
        person_update_event->dispatcher = nullptr;
        person_update_event->time = time;

        scheduler->schedule_population_event(person_update_event);
    }
}

std::string GPU::UpdateRenderPositionEvent::name() {
    return "PersonUpdateRenderPositionEvent";
}

__global__ void update_person_position(int work_from, int work_to, int work_batch, float width, float height,
                                       glm::mat4 *buffer_person_models, glm::vec4 *buffer_person_colors,
                                       curandState *state){
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    curandState local_state = state[thread_index];
    for (int index = thread_index; index < work_batch; index += stride) {
//        if(curand_uniform(&local_state) > 0.9f)
        {
            glm::mat4 model = buffer_person_models[index];
            glm::vec4 color = buffer_person_colors[index];
            float velocity = 0.00012;
            float x_n1_1 = (curand_uniform(&local_state) - 0.5f) * 2.0f;
            float y_n1_1 = (curand_uniform(&local_state) - 0.5f) * 2.0f;
            float x = x_n1_1 * width * velocity;
            float y = y_n1_1 * height * velocity;
            float rot = curand_uniform(&local_state) * 360.0f * velocity;
            model = translate(model, glm::vec3(x, y, 0.0f));
            model = translate(model, glm::vec3(0.0f, 0.0f, 1.0f));
            model = rotate(model, rot, glm::vec3(0.0f, 0.0f, 1.0f));
            model = translate(model, glm::vec3(0.0f, 0.0f, -1.0f));
            buffer_person_models[index] = model;
            buffer_person_colors[index] = color;
        }
        __syncthreads();
    }
    state[thread_index] = local_state;
}


__global__ void update_person_position_stream(int offset, int size, float width, float height,
                                       glm::mat4 *buffer_person_models, glm::vec4 *buffer_person_colors,
                                       curandState *state){
    int index = offset + threadIdx.x + blockIdx.x * blockDim.x;
    if (index < offset + size) {
        curandState local_state = state[index];
//        if(curand_uniform(&local_state) > 0.9f)
        {
            glm::mat4 model = buffer_person_models[index];
            glm::vec4 color = buffer_person_colors[index];
            float velocity = 0.00012;
            float x_n1_1 = (curand_uniform(&local_state) - 0.5f) * 2.0f;
            float y_n1_1 = (curand_uniform(&local_state) - 0.5f) * 2.0f;
            float x = x_n1_1 * width * velocity;
            float y = y_n1_1 * height * velocity;
            float rot = curand_uniform(&local_state) * 360.0f * velocity;
            model = translate(model, glm::vec3(x, y, 0.0f));
            model = translate(model, glm::vec3(0.0f, 0.0f, 1.0f));
            model = rotate(model, rot, glm::vec3(0.0f, 0.0f, 1.0f));
            model = translate(model, glm::vec3(0.0f, 0.0f, -1.0f));
            buffer_person_models[index] = model;
            buffer_person_colors[index] = color;
        }
        __syncthreads();
        state[index] = local_state;
    }
}

void GPU::UpdateRenderPositionEvent::execute() {
    auto *pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexGPU>();
    if(pi->h_persons().size() >= pi->h_person_models().size()){
        Model::CONFIG->render_config().display_gui = false;
        this->executable = false;
        return;
    }

    //Update population here
//    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO) << "[Person Update Event] executed at time: " << time;
    if(Model::CONFIG->debug_config().enable_random_position_update){
        glm::mat4 *d_person_models;
        glm::vec4 *d_person_colors;
        auto *d_streams = new cudaStream_t[Model::CONFIG->gpu_config().n_streams];
        for (int i = 0; i < Model::CONFIG->gpu_config().n_streams; ++i) {
            check_cuda_error(cudaStreamCreate(&d_streams[i]));
        }
        auto tp_start = std::chrono::high_resolution_clock::now();

        float width = Model::CONFIG->debug_config().width > 0 ? Model::CONFIG->debug_config().width : Model::CONFIG->render_config().window_width;
        float height = Model::CONFIG->debug_config().height > 0 ? Model::CONFIG->debug_config().height : Model::CONFIG->render_config().window_height;
        auto& location_db = Model::CONFIG->location_db();

        tp_start = std::chrono::high_resolution_clock::now();

        //Update person position
        int batch_size = (pi->h_persons().size() < Model::CONFIG->gpu_config().n_people_1_batch)
                     ? pi->h_persons().size() : Model::CONFIG->gpu_config().n_people_1_batch;
        int batch_offset = 0;
        for (auto [batch_remain, b_index] = std::tuple{pi->h_persons().size(), 0}; batch_remain > 0; batch_remain -= batch_size, b_index++) {
            batch_size = (batch_remain < batch_size) ? batch_remain : batch_size;
            int batch_from = pi->h_persons().size() - batch_remain;
            int batch_to = batch_from + batch_size;
            const int batch_bytes_model = batch_size * sizeof(glm::mat4);
            const int batch_bytes_color = batch_size * sizeof(glm::vec4);
            check_cuda_error(cudaMalloc(&d_person_models, batch_bytes_model));
            check_cuda_error(cudaMalloc(&d_person_colors, batch_bytes_color));
//            pi->buffer_person_models().resize(batch_size);
//            pi->buffer_person_colors().resize(batch_size);
            /*
             * Check if batch_size > Model::CONFIG->gpu_config().n_people_1_batch / 2
             * If true, then use streams
             * */
            if(batch_size > (Model::CONFIG->gpu_config().n_people_1_batch / 2)){
                int stream_size = batch_size / Model::CONFIG->gpu_config().n_streams;
                int stream_offset = 0;
                for (auto [stream_remain, s_index] = std::tuple{batch_size, 0}; stream_remain > 0; stream_remain -= stream_size, s_index++) {
                    if (s_index == Model::CONFIG->gpu_config().n_streams - 1) {
                        stream_size = stream_remain;
                    } else {
                        stream_size = (stream_remain < stream_size) ? stream_remain : stream_size;
                    }
                    const int stream_bytes_model = stream_size * sizeof(glm::mat4);
                    const int stream_bytes_color= stream_size * sizeof(glm::vec4);
                    const int n_threads = (stream_size < Model::CONFIG->gpu_config().n_threads) ? stream_size : Model::CONFIG->gpu_config().n_threads;
                    const int n_blocks = (stream_size + n_threads - 1) / n_threads;
                    check_cuda_error(cudaMemcpyAsync(&d_person_models[stream_offset],
                                                     thrust::raw_pointer_cast(pi->h_person_models().data()) + batch_offset,
                                                     stream_bytes_model,
                                                     cudaMemcpyHostToDevice,
                                                     d_streams[s_index]));
                    check_cuda_error(cudaMemcpyAsync(&d_person_colors[stream_offset],
                                                     thrust::raw_pointer_cast(pi->h_person_colors().data()) + batch_offset,
                                                     stream_bytes_color,
                                                     cudaMemcpyHostToDevice,
                                                     d_streams[s_index]));
                    update_person_position_stream<<<n_blocks, n_threads, 0, d_streams[s_index]>>>(stream_offset,
                                                                                                  stream_size,
                                                                                                  width,height,
                                                                                                  d_person_models,
                                                                                                  d_person_colors,
                                                                                                  Model::GPU_RANDOM->d_states
                    );
                    check_cuda_error(cudaDeviceSynchronize());
                    check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_models().data()) + batch_offset,
                                                     &d_person_models[stream_offset],
                                                     stream_bytes_model,
                                                     cudaMemcpyDeviceToHost,
                                                     d_streams[s_index]));
                    check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_colors().data()) + batch_offset,
                                                     &d_person_colors[stream_offset],
                                                     stream_bytes_color,
                                                     cudaMemcpyDeviceToHost,
                                                     d_streams[s_index]));
                    check_cuda_error(cudaGetLastError());
                    stream_offset += stream_size;
                    batch_offset += stream_size;
                }
                for (int s_index = 0; s_index < Model::CONFIG->gpu_config().n_streams; s_index++) {
                    check_cuda_error(cudaStreamSynchronize(d_streams[s_index]));
                }
            }
            else{
                const int n_threads = (batch_size < Model::CONFIG->gpu_config().n_threads) ? batch_size : Model::CONFIG->gpu_config().n_threads;
                const int n_blocks = (batch_size + n_threads - 1) / n_threads;
                int s_index = 0;
                check_cuda_error(cudaMemcpyAsync(&d_person_models[0],
                                                 thrust::raw_pointer_cast(pi->h_person_models().data()) + batch_offset,
                                                 batch_bytes_model,
                                                 cudaMemcpyHostToDevice,
                                                 d_streams[s_index]));
                check_cuda_error(cudaMemcpyAsync(&d_person_colors[0],
                                                 thrust::raw_pointer_cast(pi->h_person_colors().data()) + batch_offset,
                                                 batch_bytes_color,
                                                 cudaMemcpyHostToDevice,
                                                 d_streams[s_index]));
                update_person_position_stream<<<n_blocks, n_threads, 0, d_streams[s_index]>>>(0,
                                                                                              batch_size,
                                                                                              width,height,
                                                                                              d_person_models,
                                                                                              d_person_colors,
                                                                                              Model::GPU_RANDOM->d_states
                );
                check_cuda_error(cudaDeviceSynchronize());
                check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_models().data()) + batch_offset,
                                                 &d_person_models[0],
                                                 batch_bytes_model,
                                                 cudaMemcpyDeviceToHost,
                                                 d_streams[s_index]));
                check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_colors().data()) + batch_offset,
                                                 &d_person_colors[0],
                                                 batch_bytes_color,
                                                 cudaMemcpyDeviceToHost,
                                                 d_streams[s_index]));
                check_cuda_error(cudaGetLastError());
                batch_offset += batch_size;
            }
            cudaFree(d_person_models);
            cudaFree(d_person_colors);
//            /* H2D */
//            thrust::copy(pi->h_person_models().begin() + batch_from, pi->h_person_models().begin() + batch_to,pi->buffer_person_models().begin());
//            thrust::copy(pi->h_person_colors().begin() + batch_from, pi->h_person_colors().begin() + batch_to,pi->buffer_person_colors().begin());
//            const int n_threads = (batch_size < Model::CONFIG->gpu_config().n_threads) ? batch_size : Model::CONFIG->gpu_config().n_threads;
//            const int n_blocks = (batch_size + n_threads - 1) / n_threads;
//            update_person_position<<<n_blocks, n_threads>>>(batch_from,batch_to,batch_size,width,height,
//                                                                                            thrust::raw_pointer_cast(pi->buffer_person_models().data()),
//                                                                                            thrust::raw_pointer_cast(pi->buffer_person_colors().data()),
//                                                                                            Model::GPU_RANDOM->d_states);
//            check_cuda_error(cudaDeviceSynchronize());
//            check_cuda_error(cudaGetLastError());
//            /* D2H */
//            thrust::copy(pi->buffer_person_models().begin(), pi->buffer_person_models().end(), pi->h_person_models().begin() + batch_from);
//            thrust::copy(pi->buffer_person_colors().begin(), pi->buffer_person_colors().end(), pi->h_person_colors().begin() + batch_from);
//            check_cuda_error(cudaGetLastError());
        }
        cudaFree(d_streams);
//        pi->buffer_person_models().clear();
//        pi->buffer_person_colors().clear();
//        thrust::device_vector<glm::mat4>().swap(pi->buffer_person_models());
//        thrust::device_vector<glm::vec4>().swap(pi->buffer_person_colors());
        auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
        if(Model::CONFIG->debug_config().enable_debug_text){
            LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
            << fmt::format("[GPU Person Update Event] Update population movement time: {} ms",std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count());
        }
    }
}
