cmake_minimum_required(VERSION 3.18)
project(mumemto_python LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Python COMPONENTS NumPy REQUIRED)

# We intentionally do NOT modify the parent project. This binding expects the
# Mumemto CMake package to be discoverable (e.g. via CMAKE_PREFIX_PATH).
find_package(Mumemto CONFIG REQUIRED)

pybind11_add_module(_mumemto_core src/mumemto_pybind.cpp)
target_link_libraries(_mumemto_core PRIVATE Mumemto::mumemto)
target_compile_features(_mumemto_core PRIVATE cxx_std_17)

target_link_libraries(_mumemto_core PRIVATE Python::NumPy)

# Ensure the link step actually records DT_NEEDED for libmumemto.
# With --as-needed (common in some toolchains), a shared lib may be dropped when
# producing a Python extension, leaving symbols unresolved until import time.
if(UNIX AND NOT APPLE)
  target_link_options(_mumemto_core PRIVATE "LINKER:--no-as-needed")
endif()

# Name the extension as mumemto._mumemto_core
set_target_properties(_mumemto_core PROPERTIES
  OUTPUT_NAME "_mumemto_core"
)

# For local testing, embed an rpath to the resolved Mumemto shared library dir.
# This avoids requiring users to manually set LD_LIBRARY_PATH.
set(_mumemto_libdir "$<TARGET_FILE_DIR:Mumemto::mumemto>")
if(UNIX AND NOT APPLE)
  set_target_properties(_mumemto_core PROPERTIES
    BUILD_RPATH "${_mumemto_libdir}"
    INSTALL_RPATH "${_mumemto_libdir}"
  )
elseif(APPLE)
  set_target_properties(_mumemto_core PROPERTIES
    BUILD_RPATH "${_mumemto_libdir}"
    INSTALL_RPATH "${_mumemto_libdir}"
  )
endif()

# Ensure the extension is installed inside the Python package directory.
install(TARGETS _mumemto_core
  LIBRARY DESTINATION mumemto
  RUNTIME DESTINATION mumemto
  ARCHIVE DESTINATION mumemto
)
