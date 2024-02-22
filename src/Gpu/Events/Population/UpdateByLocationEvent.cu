//
// Created by kient on 12/24/2023.
//

#include "UpdateByLocationEvent.cuh"
#include "easylogging++.h"

#include "Gpu/Core/Scheduler.cuh"
#include "Model.h"
#include "Gpu/Utils/Utils.cuh"
#include "Core/Config/Config.h"
#include "Gpu/MDC/ModelDataCollector.cuh"
#include "Gpu/Population/Properties/PersonIndexGPU.cuh"
#include "Gpu/Population/Population.cuh"

GPU::UpdateByLocationEvent::UpdateByLocationEvent() = default;

GPU::UpdateByLocationEvent::~UpdateByLocationEvent() = default;

void GPU::UpdateByLocationEvent::schedule_event(GPU::Scheduler *scheduler, const int &time) {
    if (scheduler!=nullptr) {
        auto *person_update_event = new GPU::UpdateByLocationEvent();
        person_update_event->dispatcher = nullptr;
        person_update_event->time = time;

        scheduler->schedule_population_event(person_update_event);
    }
}

std::string GPU::UpdateByLocationEvent::name() {
    return "PersonUpdateByLocationEvent";
}


void GPU::UpdateByLocationEvent::execute() {
    //Update population here
//    LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO) << "[Person Update Event] executed at time: " << time;
    auto tp_start = std::chrono::high_resolution_clock::now();

    auto *pi = Model::GPU_POPULATION->get_person_index<GPU::PersonIndexGPU>();

    Model::GPU_DATA_COLLECTOR->set_popsize_residence_by_location(Model::GPU_UTILS->count_by_1key<int>(pi->h_person_residence_locations(),
                                                                              Model::CONFIG->location_db().size()));

    if(Model::CONFIG->debug_config().enable_debug_text){
        auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
        LOG_IF(Model::GPU_SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
        << fmt::format("[GPU Person Update Event] Update population by location time: {} ms",std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count());
    }
}