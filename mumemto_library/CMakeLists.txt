cmake_minimum_required(VERSION 3.13)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Parent project may not initialize GNUInstallDirs before this subdirectory; an empty
# CMAKE_INSTALL_LIBDIR yields broken absolute paths like /cmake/Mumemto on install.
if(NOT CMAKE_INSTALL_LIBDIR)
    set(CMAKE_INSTALL_LIBDIR "lib")
endif()
if(NOT CMAKE_INSTALL_BINDIR)
    set(CMAKE_INSTALL_BINDIR "bin")
endif()
if(NOT CMAKE_INSTALL_INCLUDEDIR)
    set(CMAKE_INSTALL_INCLUDEDIR "include")
endif()

find_package(ZLIB REQUIRED)

# ------------------------------------------------------------------------------
# Shared library (`mumemto`) is the installable SDK artifact (conda).
# ------------------------------------------------------------------------------

set(MUMEMTO_CXX_STD cxx_std_17)

# Code linked into a shared library must be PIC (sdsl + our object libs).
foreach(_mumemto_pic_target common pfp gsacak64 sdsl)
    if(TARGET ${_mumemto_pic_target})
        set_target_properties(${_mumemto_pic_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif()
endforeach()

add_library(mumemto SHARED
    mumemto_api.cpp
    mumemto_util.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../src/ref_builder.cpp
    $<TARGET_OBJECTS:common>
    $<TARGET_OBJECTS:pfp>
    $<TARGET_OBJECTS:gsacak64>
)

set_target_properties(mumemto PROPERTIES
    OUTPUT_NAME mumemto
    VERSION ${VERSION}
    SOVERSION ${VERSION_MAJOR}
)

if(MSVC)
    set_target_properties(mumemto PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

target_include_directories(mumemto PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_include_directories(mumemto PRIVATE
    $<TARGET_PROPERTY:gsacak64,INTERFACE_INCLUDE_DIRECTORIES>
)

# zlib is shared on most installs; sdsl is linked in from the static libsdsl target.
target_link_libraries(mumemto PRIVATE sdsl ZLIB::ZLIB)

target_compile_features(mumemto PUBLIC ${MUMEMTO_CXX_STD})
target_compile_definitions(mumemto PRIVATE M64)
target_compile_options(mumemto PRIVATE "-std=c++17" "-w")

if(MUMEMTO_SHARED_EXTRA_CXX_FLAGS)
    target_compile_options(mumemto PRIVATE ${MUMEMTO_SHARED_EXTRA_CXX_FLAGS})
endif()

if(UNIX AND NOT APPLE)
    set_target_properties(mumemto PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        BUILD_WITH_INSTALL_RPATH FALSE
    )
elseif(APPLE)
    set_target_properties(mumemto PROPERTIES
        INSTALL_RPATH "@loader_path"
        BUILD_WITH_INSTALL_RPATH FALSE
    )
endif()

# ------------------------------------------------------------------------------
# Install: library, public headers, CMake package (find_package Mumemto).
# ------------------------------------------------------------------------------

install(TARGETS mumemto
    EXPORT MumemtoTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES
    "${CMAKE_CURRENT_LIST_DIR}/mumemto.h"
    "${CMAKE_CURRENT_LIST_DIR}/mumemto_api.hpp"
    "${CMAKE_CURRENT_LIST_DIR}/mumemto.hpp"
    "${CMAKE_CURRENT_LIST_DIR}/../include/mumsio.hpp"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT MumemtoTargets
    FILE MumemtoTargets.cmake
    NAMESPACE Mumemto::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mumemto
)

configure_package_config_file(
    "${CMAKE_CURRENT_LIST_DIR}/MumemtoConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/MumemtoConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mumemto
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/MumemtoConfigVersion.cmake"
    VERSION ${VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/MumemtoConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/MumemtoConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mumemto
)
