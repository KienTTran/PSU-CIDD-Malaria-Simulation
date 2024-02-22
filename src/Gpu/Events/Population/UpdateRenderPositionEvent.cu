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
        auto *person_update_event = new GPU::UpdateRenderPositionEvent();
        person_update_event->dispatcher = nullptr;
        person_update_event->time = time;

        scheduler->schedule_population_event(person_update_event);
    }
}

std::string GPU::UpdateRenderPositionEvent::name() {
    return "PersonUpdateRenderPositionEvent";
}

__global__ void update_person_position_model_stream(int offset, int size, float width, float height,
                                              glm::mat4 *buffer_person_models,
                                              curandState *state){
  int index = offset + threadIdx.x + blockIdx.x * blockDim.x;
  if (index < offset + size) {
    curandState local_state = state[index];
//        if(curand_uniform(&local_state) > 0.9f)
    {
      glm::mat4 model = buffer_person_models[index];
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
    }
    __syncthreads();
    state[index] = local_state;
  }
}


__global__ void update_person_position_color_stream(int offset, int size, float width, float height,
                                              glm::vec4 *buffer_person_colors,
                                              curandState *state){
  int index = offset + threadIdx.x + blockIdx.x * blockDim.x;
  if (index < offset + size) {
    curandState local_state = state[index];
//        if(curand_uniform(&local_state) > 0.9f)
    {
      glm::vec4 color = buffer_person_colors[index];
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

    /*
     * Update population here
     * Stream size is not changed
     * Only number of stream is doubled to handle model and color
     * */
//    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO) << "[Person Update Event] executed at time: " << time;
    if(Model::CONFIG->render_config().display_gui && Model::CONFIG->debug_config().enable_random_position_update){
        auto tp_start = std::chrono::high_resolution_clock::now();
        int n_streams_x2 = Model::CONFIG->gpu_config().n_streams*2;
        auto *d_streams = new cudaStream_t[n_streams_x2];
        for (int i = 0; i < n_streams_x2; ++i) {
            check_cuda_error(cudaStreamCreate(&d_streams[i]));
        }

        float width = Model::CONFIG->debug_config().width > 0 ? Model::CONFIG->debug_config().width : Model::CONFIG->render_config().window_width;
        float height = Model::CONFIG->debug_config().height > 0 ? Model::CONFIG->debug_config().height : Model::CONFIG->render_config().window_height;

        //Update person position
        int batch_size = (pi->h_persons().size() < Model::CONFIG->gpu_config().n_people_1_batch)
                        ? pi->h_persons().size() : Model::CONFIG->gpu_config().n_people_1_batch;
        int batch_offset = 0;
//        Model::GPU_RANDOM->init(pi->h_persons().size());
        for (auto [batch_remain, b_index] = std::tuple{pi->h_persons().size(), 0}; batch_remain > 0; batch_remain -= batch_size, b_index++) {
            batch_size = (batch_remain < batch_size) ? batch_remain : batch_size;
            int batch_from = pi->h_persons().size() - batch_remain;
            int batch_to = batch_from + batch_size;
//            printf("##### BATCH %d size %d remain %d from %d to %d #####\n", b_index, batch_size, batch_remain, batch_from, batch_to);
            glm::mat4 *d_person_models;
            glm::vec4 *d_person_colors;
            const int batch_bytes_model = batch_size * sizeof(glm::mat4);
            const int batch_bytes_color = batch_size * sizeof(glm::vec4);
            check_cuda_error(cudaMalloc(&d_person_models, batch_bytes_model));
            check_cuda_error(cudaMalloc(&d_person_colors, batch_bytes_color));
            if(batch_size > (Model::CONFIG->gpu_config().n_people_1_batch / 2)){
//              printf("  ##### STREAMS #####\n");
              int stream_size = batch_size / (Model::CONFIG->gpu_config().n_streams);
              int stream_offset = 0;
              for (auto [stream_remain, s_index] = std::tuple{batch_size, 0}; stream_remain > 0; stream_remain -= stream_size, s_index+=2) {
                  if (s_index == n_streams_x2 - 2) {
                      stream_size = stream_remain;
                  } else {
                      stream_size = (stream_remain < stream_size) ? stream_remain : stream_size;
                  }
//                  printf("  ##### STREAM %d %d size %d remain %d from %d to %d #####\n", s_index, s_index+1,stream_size, stream_remain, stream_offset, stream_offset + stream_size);
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
                                                   d_streams[s_index+1]));
                  update_person_position_model_stream<<<n_blocks, n_threads, 0, d_streams[s_index]>>>(stream_offset,
                                                                                                stream_size,
                                                                                                width,height,
                                                                                                d_person_models,
                                                                                                Model::GPU_RANDOM->d_states);
                  check_cuda_error(cudaDeviceSynchronize());
                  update_person_position_color_stream<<<n_blocks, n_threads, 0, d_streams[s_index+1]>>>(stream_offset,
                                                                                                stream_size,
                                                                                                width,height,
                                                                                                d_person_colors,
                                                                                                Model::GPU_RANDOM->d_states);
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
                                                   d_streams[s_index+1]));
                  check_cuda_error(cudaGetLastError());
                  stream_offset += stream_size;
                  batch_offset += stream_size;
              }
              for (int s_index = 0; s_index < n_streams_x2; s_index++) {
                  check_cuda_error(cudaStreamSynchronize(d_streams[s_index]));
              }
            }
            else{
//              printf("##### LAST BATCH %d size %d remain %d from %d to %d #####\n", b_index, batch_size, batch_remain, batch_from, batch_to);
              const int n_threads = (batch_size < Model::CONFIG->gpu_config().n_threads) ? batch_size : Model::CONFIG->gpu_config().n_threads;
              const int n_blocks = (batch_size + n_threads - 1) / n_threads;
              check_cuda_error(cudaMemcpyAsync(&d_person_models[0],
                                               thrust::raw_pointer_cast(pi->h_person_models().data()) + batch_offset,
                                               batch_bytes_model,
                                               cudaMemcpyHostToDevice,
                                               d_streams[0]));
              check_cuda_error(cudaMemcpyAsync(&d_person_colors[0],
                                               thrust::raw_pointer_cast(pi->h_person_colors().data()) + batch_offset,
                                               batch_bytes_color,
                                               cudaMemcpyHostToDevice,
                                               d_streams[1]));
              update_person_position_model_stream<<<n_blocks, n_threads, 0, d_streams[0]>>>(0,
                                                                                            batch_size,
                                                                                            width,height,
                                                                                            d_person_models,
                                                                                            Model::GPU_RANDOM->d_states);
              check_cuda_error(cudaDeviceSynchronize());
              update_person_position_color_stream<<<n_blocks, n_threads, 0, d_streams[1]>>>(0,
                                                                                            batch_size,
                                                                                            width,height,
                                                                                            d_person_colors,
                                                                                            Model::GPU_RANDOM->d_states);
              check_cuda_error(cudaDeviceSynchronize());
              check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_models().data()) + batch_offset,
                                               &d_person_models[0],
                                               batch_bytes_model,
                                               cudaMemcpyDeviceToHost,
                                               d_streams[0]));
              check_cuda_error(cudaMemcpyAsync(thrust::raw_pointer_cast(pi->h_person_colors().data()) + batch_offset,
                                               &d_person_colors[0],
                                               batch_bytes_color,
                                               cudaMemcpyDeviceToHost,
                                               d_streams[1]));
              check_cuda_error(cudaGetLastError());
              batch_offset += batch_size;
            }
            cudaFree(d_person_models);
            cudaFree(d_person_colors);

//          LOG(INFO) << fmt::format("[Population update POS] Work batch size {} remain {}, from {} to {} (of {} {})", batch_size, batch_remain, batch_from, batch_to,
//                   pi->h_person_models().size(),pi->h_person_colors().size());
//          pi->buffer_person_models().resize(batch_size);
//          pi->buffer_person_colors().resize(batch_size);
//          /* H2D */
//          thrust::copy(pi->h_person_models().begin() + batch_from, pi->h_person_models().begin() + batch_to,pi->buffer_person_models().begin());
//          thrust::copy(pi->h_person_colors().begin() + batch_from, pi->h_person_colors().begin() + batch_to,pi->buffer_person_colors().begin());
//          const int n_threads = (batch_size < Model::CONFIG->gpu_config().n_threads) ? batch_size : Model::CONFIG->gpu_config().n_threads;
//          const int n_blocks = (batch_size + n_threads - 1) / n_threads;
//          update_person_position<<<n_blocks, n_threads>>>(batch_from,batch_to,batch_size,width,height,
//                                                          thrust::raw_pointer_cast(pi->buffer_person_models().data()),
//                                                          thrust::raw_pointer_cast(pi->buffer_person_colors().data()),
//                                                          Model::GPU_RANDOM->d_states);
//          check_cuda_error(cudaDeviceSynchronize());
//          check_cuda_error(cudaGetLastError());
//          /* D2H */
//          thrust::copy(pi->buffer_person_models().begin(), pi->buffer_person_models().end(), pi->h_person_models().begin() + batch_from);
//          thrust::copy(pi->buffer_person_colors().begin(), pi->buffer_person_colors().end(), pi->h_person_colors().begin() + batch_from);
//          check_cuda_error(cudaGetLastError());
        }
        for (int s_index = 0; s_index < n_streams_x2; s_index++) {
            check_cuda_error(cudaStreamDestroy(d_streams[s_index]));
        }
        auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
        if(Model::CONFIG->debug_config().enable_debug_text){
            LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
            << fmt::format("[GPU Person Update Event] Update population movement time: {} ms",std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count());
        }
    }
}
