add_executable(mumemto pfp_mum.cpp ref_builder.cpp)
target_link_libraries(mumemto common pfp gsacak64 sdsl)# ${PROFILER_LIB})
target_include_directories(mumemto PUBLIC "../include/")
target_compile_options(mumemto PUBLIC "-std=c++17" "-DM64" "-march=native" "-w")