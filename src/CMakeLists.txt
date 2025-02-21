add_executable(mumemto_exec pfp_mum.cpp ref_builder.cpp)
target_link_libraries(mumemto_exec common pfp gsacak64 sdsl)# ${PROFILER_LIB})
target_include_directories(mumemto_exec PUBLIC "../include/")
target_compile_options(mumemto_exec PUBLIC "-std=c++17" "-DM64" "-march=native" "-w")

# Add new compute_lengths executable
add_executable(compute_lengths compute_lengths.cpp)
target_include_directories(compute_lengths PUBLIC "../include/")
target_compile_options(compute_lengths PUBLIC "-std=c++17" "-DM64" "-march=native" "-w")