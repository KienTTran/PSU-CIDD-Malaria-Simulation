#include "IStrategy.cuh"

std::map<std::string, GPU::IStrategy::StrategyType> GPU::IStrategy::StrategyTypeMap{
    {"SFT",                    SFT},
    {"Cycling",                Cycling},
    {"AdaptiveCycling",        AdaptiveCycling},
    {"MFT",                    MFT},
    {"MFTRebalancing",         MFTRebalancing},
    {"NestedMFT",              NestedMFT},
    {"MFTMultiLocation",       MFTMultiLocation},
    {"NestedMFTMultiLocation", NestedMFTMultiLocation},
    {"NovelDrugIntroduction",     NovelDrugIntroduction},
};
