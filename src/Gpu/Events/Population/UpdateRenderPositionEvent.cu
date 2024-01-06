//
// Created by kient on 8/4/2023.
//

#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/remove.h>
#include "easylogging++.h"

#include "Core/Scheduler.h"
#include "Model.h"
#include "Core/Config/Config.h"
#include "Gpu/Population/Properties/PersonIndexGPU.h"

#include "Gpu/Random.cuh"
#include "Gpu/RenderEntity.cuh"
#include "Gpu/Utils.cuh"
#include "Population/Population.h"

#include "UpdateRenderPositionEvent.cuh"

UpdateRenderPositionEvent::UpdateRenderPositionEvent() = default;

UpdateRenderPositionEvent::~UpdateRenderPositionEvent() = default;

void UpdateRenderPositionEvent::schedule_event(Scheduler *scheduler, const int &time) {
    if (scheduler!=nullptr) {
        auto *person_update_event = new UpdateRenderPositionEvent();
        person_update_event->dispatcher = nullptr;
        person_update_event->time = time;

        scheduler->schedule_population_event(person_update_event);
    }
}

std::string UpdateRenderPositionEvent::name() {
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

void UpdateRenderPositionEvent::execute() {
    auto *pi = Model::POPULATION->get_person_index<PersonIndexGPU>();
    if(pi->h_persons().size() >= pi->h_person_models().size()){
        Model::CONFIG->render_config().display_gui = false;
        this->executable = false;
        return;
    }

    //Update population here
//    LOG_IF(Model::SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO) << "[Person Update Event] executed at time: " << time;
    if(Model::CONFIG->debug_config().enable_update){
        auto tp_start = std::chrono::high_resolution_clock::now();

        float width = Model::CONFIG->debug_config().width > 0 ? Model::CONFIG->debug_config().width : Model::CONFIG->render_config().window_width;
        float height = Model::CONFIG->debug_config().height > 0 ? Model::CONFIG->debug_config().height : Model::CONFIG->render_config().window_height;
        int n_threads = Model::CONFIG->gpu_config().n_threads;
        auto& location_db = Model::CONFIG->location_db();
        int batch_size;

        tp_start = std::chrono::high_resolution_clock::now();

        //Update person position
        batch_size = (pi->h_persons().size() < Model::CONFIG->gpu_config().people_1_batch)
                     ? pi->h_persons().size() : Model::CONFIG->gpu_config().people_1_batch;
        //This is to make sure threads fit all people in population
        n_threads = (batch_size < n_threads) ? batch_size : n_threads;
        for (int remain = pi->h_persons().size(); remain > 0; remain -= batch_size) {
            batch_size = (remain < batch_size) ? remain : batch_size;
            int batch_from = pi->h_persons().size() - remain;
            int batch_to = batch_from + batch_size;
//            LOG(INFO) << fmt::format("[Population update] Work batch size %d remain %d, from %d to %d (of %d %d)", batch_size, remain, batch_from, batch_to,
//                   pi->h_person_models().size(),pi->h_person_colors().size());
            pi->buffer_person_models().resize(batch_size);
            pi->buffer_person_colors().resize(batch_size);
            /* H2D */
            thrust::copy(pi->h_person_models().begin() + batch_from, pi->h_person_models().begin() + batch_to,pi->buffer_person_models().begin());
            thrust::copy(pi->h_person_colors().begin() + batch_from, pi->h_person_colors().begin() + batch_to,pi->buffer_person_colors().begin());
            update_person_position<<<((batch_size + n_threads - 1)/n_threads), n_threads>>>(batch_from,batch_to,batch_size,width,height,
                    thrust::raw_pointer_cast(pi->buffer_person_models().data()),
                    thrust::raw_pointer_cast(pi->buffer_person_colors().data()),
                    Model::GPU_RANDOM->d_states);
            check_cuda_error(cudaDeviceSynchronize());
            check_cuda_error(cudaGetLastError());
            /* D2H */
            thrust::copy(pi->buffer_person_models().begin(), pi->buffer_person_models().end(), pi->h_person_models().begin() + batch_from);
            thrust::copy(pi->buffer_person_colors().begin(), pi->buffer_person_colors().end(), pi->h_person_colors().begin() + batch_from);
            check_cuda_error(cudaGetLastError());
        }
        pi->buffer_person_models().clear();
        pi->buffer_person_colors().clear();
        thrust::device_vector<glm::mat4>().swap(pi->buffer_person_models());
        thrust::device_vector<glm::vec4>().swap(pi->buffer_person_colors());

        auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
        if(Model::CONFIG->debug_config().enable_debug_text){
            LOG_IF(Model::SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
            << fmt::format("[GPU Person Update Event] Update population movement time: {} ms",std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count());
        }
    }
}
