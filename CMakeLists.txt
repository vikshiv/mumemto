CMAKE_MINIMUM_REQUIRED(VERSION 3.13)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

# ------------------------------------------------------------------------------
# Set a default build type if none was specified
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# About this project
# ------------------------------------------------------------------------------
project(pfp-mum)
SET(VERSION_MAJOR "1")
SET(VERSION_MINOR "0")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

option(COMPILE_BENCHMARKS "Set ON to compile the benchmarks" OFF)
option(COMPILE_TESTS "Set ON to compile the tests" OFF)

# ------------------------------------------------------------------------------
# Set environment
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install)

find_package(Git)
if(GIT_FOUND)
    message("git found: ${GIT_EXECUTABLE}")
else()
    message(WARNING "git not found. Cloning of submodules will not work.")
endif()

# ------------------------------------------------------------------------------
# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include") # This is an hack because include(GUIInstallDirs) doesn't work
add_subdirectory(thirdparty)

# ------------------------------------------------------------------------------
# Configure the compiler with the appropriate flags
# ------------------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  include(ConfigureCompilerClang)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  include(ConfigureCompilerGcc)
else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message(FATAL_ERROR "Only the compiler gcc and clang are supported")
endif()

add_subdirectory(include)
add_subdirectory(src)

install(TARGETS mumemto DESTINATION ${PROJECT_BINARY_DIR} RENAME mumemto_exec)
#install(TARGETS viz DESTINATION ${PROJECT_BINARY_DIR})
# install(TARGETS newscanNT.x pscan.x newscan.x 
#         DESTINATION ${PROJECT_BINARY_DIR}/bin)
install(TARGETS newscanNT.x DESTINATION ${PROJECT_BINARY_DIR})

# Copy and rename mumemto Python package to build directory
file(COPY ${CMAKE_SOURCE_DIR}/mumemto/ DESTINATION ${PROJECT_BINARY_DIR}/_mumemto)
# Move mumemto executable from package directory to build directory root
file(COPY ${CMAKE_SOURCE_DIR}/mumemto/mumemto DESTINATION ${PROJECT_BINARY_DIR})
# Make mumemto script executable
file(CHMOD ${PROJECT_BINARY_DIR}/mumemto PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
