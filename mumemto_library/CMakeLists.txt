cmake_minimum_required(VERSION 3.13)

# Small C++ library wrapper around the in-repo mumemto API.

add_library(mumemto_library STATIC
    mumemto_api.cpp
)

target_include_directories(mumemto_library
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
        ${CMAKE_CURRENT_LIST_DIR}/../include
)

target_link_libraries(mumemto_library
    PUBLIC
        common
        pfp
        gsacak64
        sdsl
        z
)

target_compile_options(mumemto_library PUBLIC "-std=c++17" "-DM64" "-w")

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    target_compile_options(mumemto_library PUBLIC "-march=native")
endif()

