find_package(CLI11 CONFIG REQUIRED)
file(GLOB DXG_SRC_FILES *.cpp *.h)

add_executable(DxGGenerator_PfCRT ${DXG_SRC_FILES})
add_dependencies(DxGGenerator_PfCRT  MaSimCore)
target_link_libraries(DxGGenerator_PfCRT PRIVATE MaSimCore CLI11::CLI11)
target_compile_features(DxGGenerator_PfCRT PRIVATE cxx_range_for)