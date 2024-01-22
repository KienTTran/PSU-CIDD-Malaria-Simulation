//
// Created by kient on 8/4/2023.
//
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include "easylogging++.h"

#include "Gpu/Core/Scheduler.cuh"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/Population/Properties/PersonIndexGPU.cuh"

#include "Gpu/Core/Random.cuh"
#include "Gpu/Renderer/RenderEntity.cuh"
#include "Gpu/Utils/Utils.cuh"
#include "Gpu/Population/Population.cuh"

#include "UpdateRenderOGLEvent.cuh"

GPU::UpdateRenderOGLEvent::UpdateRenderOGLEvent() = default;

GPU::UpdateRenderOGLEvent::~UpdateRenderOGLEvent() = default;

void GPU::UpdateRenderOGLEvent::schedule_event(GPU::Scheduler *scheduler, const int &time) {
    if (scheduler!=nullptr) {
        auto *person_update_event = new UpdateRenderOGLEvent();
        person_update_event->dispatcher = nullptr;
        person_update_event->time = time;

        scheduler->schedule_population_event(person_update_event);
    }
}

std::string GPU::UpdateRenderOGLEvent::name() {
    return "PersonUpdateRenderOGLEvent";
}

__global__ void update_ogl_buffer(int work_from, int work_to, int work_batch,
                                  glm::mat4 buffer_person_models[],
                                  glm::vec4 buffer_person_colors[],
                                  glm::mat4 ogl_person_models[],
                                  glm::vec4 ogl_person_colors[]){
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int index = thread_index; index < work_batch; index += stride) {
        ogl_person_models[index] = buffer_person_models[index];
        ogl_person_colors[index] = buffer_person_colors[index];
        __syncthreads();
    }
}

void GPU::UpdateRenderOGLEvent::execute() {
    auto *pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexGPU>();
    if(pi->h_persons().size() >= pi->h_person_models().size()){
        Model::CONFIG->render_config().display_gui = false;
        this->executable = false;
        return;
    }

    //Update population here
//    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO) << "[Person Update Render] Event executed at time: " << time;
    auto tp_start = std::chrono::high_resolution_clock::now();

    int n_threads = Model::CONFIG->gpu_config().n_threads;
    int batch_size = (pi->h_persons().size() < Model::CONFIG->gpu_config().people_1_batch)
                     ? pi->h_persons().size() : Model::CONFIG->gpu_config().people_1_batch;
    //This is to make sure threads fit all people in population
    n_threads = (batch_size < n_threads) ? batch_size : n_threads;
    for (int remain = pi->h_persons().size(); remain > 0; remain -= batch_size) {
        batch_size = (remain < batch_size) ? remain : batch_size;
        int batch_from = pi->h_persons().size() - remain;
        int batch_to = batch_from + batch_size;
//        printf("[GPUBUffer update render] Work batch size %d remain %d, from %d to %d\n", batch_size, remain, batch_from, batch_to);
        pi->buffer_person_models().resize(batch_size);
        pi->buffer_person_colors().resize(batch_size);
        /* H2D */
        thrust::copy(pi->h_person_models().begin() + batch_from, pi->h_person_models().begin() + batch_to,pi->buffer_person_models().begin());
        thrust::copy(pi->h_person_colors().begin() + batch_from, pi->h_person_colors().begin() + batch_to,pi->buffer_person_colors().begin());
        //Update OGL buffer, this must be in final step
        update_ogl_buffer<<<((batch_size + n_threads - 1)/n_threads), n_threads>>>(batch_from,batch_to,batch_size,
                                                                                   thrust::raw_pointer_cast(pi->buffer_person_models().data()),
                                                                                   thrust::raw_pointer_cast(pi->buffer_person_colors().data()),
                                                                                   Model::GPU_RENDER_ENTITY->d_ogl_buffer_model_ptr,
                                                                                   Model::GPU_RENDER_ENTITY->d_ogl_buffer_color_ptr);
        cudaDeviceSynchronize();
        check_cuda_error(cudaGetLastError());
    }
    pi->buffer_person_models().clear();
    pi->buffer_person_colors().clear();
    thrust::device_vector<glm::mat4>().swap(pi->buffer_person_models());
    thrust::device_vector<glm::vec4>().swap(pi->buffer_person_colors());

    auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
    if(Model::CONFIG->debug_config().enable_debug_text) {
        LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
        << fmt::format("[GPU Person Update Render] Update population render ({} {}) time: {} ms",
               pi->h_persons().size(),pi->h_person_models().size(),
               std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count());
    }
}