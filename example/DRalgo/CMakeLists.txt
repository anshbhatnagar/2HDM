add_executable(run_DRalgo_ah run_DRalgo_ah.cpp)
add_executable(run_DRalgo_xSM run_DRalgo_xSM.cpp)
add_executable(run_DRalgo_2HDM run_DRalgo_2HDM.cpp)
target_link_libraries(run_DRalgo_2HDM PRIVATE nlohmann_json::nlohmann_json)
add_executable(run_2HDM_ST run_2HDM_ST.cpp)
target_link_libraries(run_2HDM_ST PRIVATE nlohmann_json::nlohmann_json)
add_executable(2HDM_runs 2HDM_runs.cpp ../ProgressBar/progress.cpp)
target_link_libraries(2HDM_runs PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(2HDM_runs PRIVATE nlohmann_json::nlohmann_json)
add_executable(2HDM_vt_runs 2HDM_vt_runs.cpp ../ProgressBar/progress.cpp)
target_link_libraries(2HDM_vt_runs PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(2HDM_vt_runs PRIVATE nlohmann_json::nlohmann_json)
