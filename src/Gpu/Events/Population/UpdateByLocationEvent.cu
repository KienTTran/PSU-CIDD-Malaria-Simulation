//
// Created by kient on 12/24/2023.
//

#include "UpdateByLocationEvent.cuh"
#include "easylogging++.h"

#include "Core/Scheduler.h"
#include "Model.h"
#include "Gpu/Utils.cuh"
#include "Core/Config/Config.h"
#include "MDC/ModelDataCollector.h"
#include "Gpu/Population/Properties/PersonIndexGPU.h"
#include "Population/Population.h"

UpdateByLocationEvent::UpdateByLocationEvent() = default;

UpdateByLocationEvent::~UpdateByLocationEvent() = default;

void UpdateByLocationEvent::schedule_event(Scheduler *scheduler, const int &time) {
    if (scheduler!=nullptr) {
        auto *person_update_event = new UpdateByLocationEvent();
        person_update_event->dispatcher = nullptr;
        person_update_event->time = time;

        scheduler->schedule_population_event(person_update_event);
    }
}

std::string UpdateByLocationEvent::name() {
    return "PersonUpdateByLocationEvent";
}


void UpdateByLocationEvent::execute() {
    //Update population here
//    LOG_IF(Model::SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO) << "[Person Update Event] executed at time: " << time;
    auto tp_start = std::chrono::high_resolution_clock::now();

    auto *pi = Model::POPULATION->get_person_index<PersonIndexGPU>();

    Model::DATA_COLLECTOR->set_popsize_residence_by_location_gpu(Model::GPU_UTILS->count_by_1key<int>(pi->h_person_residence_locations(),
                                                                              Model::CONFIG->location_db().size()));

    if(Model::CONFIG->debug_config().enable_debug_text){
        auto lapse = std::chrono::high_resolution_clock::now() - tp_start;
        LOG_IF(Model::SCHEDULER->current_time() % Model::CONFIG->debug_config().log_interval == 0, INFO)
        << fmt::format("[GPU Person Update Event] Update population by location time: {} ms",std::chrono::duration_cast<std::chrono::milliseconds>(lapse).count());
    }
}