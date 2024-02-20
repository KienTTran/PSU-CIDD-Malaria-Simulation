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
  if (scheduler != nullptr) {
    auto *person_update_event = new UpdateRenderOGLEvent();
    person_update_event->dispatcher = nullptr;
    person_update_event->time = time;

    scheduler->schedule_population_event(person_update_event);
  }
}

std::string GPU::UpdateRenderOGLEvent::name() {
  return "PersonUpdateRenderOGLEvent";
}

void GPU::UpdateRenderOGLEvent::execute() {
  auto *pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexGPU>();
  if (pi->h_persons().size() >= pi->h_person_models().size()) {
    Model::CONFIG->render_config().display_gui = false;
    this->executable = false;
    return;
  }

//    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO) << "[Person Update Render] Event executed at time: " << time;
  if (Model::CONFIG->render_config().display_gui){
    auto tp_start = std::chrono::high_resolution_clock::now();

  check_cuda_error(cudaMemcpy(Model::GPU_RENDER_ENTITY->d_ogl_buffer_model_ptr,
                              thrust::raw_pointer_cast(pi->h_person_models().data()), pi->h_person_models().size()*sizeof(glm::mat4), cudaMemcpyHostToDevice));
  check_cuda_error(cudaMemcpy(Model::GPU_RENDER_ENTITY->d_ogl_buffer_color_ptr,
                              thrust::raw_pointer_cast(pi->h_person_colors().data()), pi->h_person_models().size()*sizeof(glm::vec4), cudaMemcpyHostToDevice));
  check_cuda_error(cudaDeviceSynchronize());
  check_cuda_error(cudaGetLastError());

    auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
    if (Model::CONFIG->debug_config().enable_debug_text) {
      LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
        << fmt::format("[GPU Person Update Render] Update population render ({} {}) time: {} ms",
                       pi->h_persons().size(), pi->h_person_models().size(),
                       std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count());
    }
  }
}