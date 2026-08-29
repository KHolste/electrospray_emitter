cmake_minimum_required(VERSION 3.20)
project(electrospray LANGUAGES CXX
        DESCRIPTION "Simulation of colloid / electrospray emitters")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

find_package(OpenMP)

# Host-CPU codegen is opt-in and deliberately NOT the default.  Besides
# producing binaries that only run on the machine that built them,
# -march=native breaks this build on MinGW-w64 GCC 16.1 as soon as OpenMP is
# also enabled: -mavx2 together with libgomp segfaults during teardown.
# Neither alone does, and -mstackrealign does not help, so it is a toolchain
# defect rather than a code bug.  The numbers agreed with the scalar build in
# every case checked -- but a toolchain that crashes on exit is not one to
# trust with silent numerics.
option(ES_NATIVE_ARCH "Build for the host CPU (-march=native).  See the comment above." OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
  # -ffast-math is deliberately NOT used: the elliptic-integral AGM and the
  # singularity subtraction rely on well-behaved rounding.
  set(CMAKE_CXX_FLAGS_RELEASE "-O2")
  if(ES_NATIVE_ARCH)
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native")
    if(OpenMP_CXX_FOUND)
      message(WARNING "ES_NATIVE_ARCH with OpenMP is known to crash on MinGW GCC 16.1.")
    endif()
  endif()
endif()

add_library(es_core STATIC
  src/status.cpp
  src/elliptic.cpp
  src/linalg.cpp
  src/fluid.cpp
  src/liquid.cpp
  src/feed.cpp
  src/transport.cpp
  src/geometry.cpp
  src/device_geometry.cpp
  src/boundary_mesh.cpp
  src/materials.cpp
  src/material_data.cpp
  src/material_data_emibf4.cpp
  src/axisym_fem.cpp
  src/volume_mesh.cpp
  src/dielectric_device.cpp
  src/bem.cpp
  src/vacuum_bem.cpp
  src/emission.cpp
  src/meniscus.cpp
  src/capillary.cpp
  src/electrocapillary.cpp
  src/load_projection.cpp
  src/beam.cpp
  src/config.cpp
  src/io.cpp
)
target_include_directories(es_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(OpenMP_CXX_FOUND)
  target_link_libraries(es_core PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(es_core PUBLIC ES_HAVE_OPENMP)
endif()

# ---------------------------------------------------------------- applications
foreach(app es_field es_meniscus es_operate es_beam es_vacuum es_dielectric es_reservoir
             es_capillary es_electrocapillary es_p3b_audit es_feed es_material es_transport)
  add_executable(${app} apps/${app}.cpp)
  target_link_libraries(${app} PRIVATE es_core)
endforeach()

# ----------------------------------------------------------------- diagnostics
# Not part of the physics: regenerates the data behind docs/figures.
add_executable(branch_figure tools/branch_figure.cpp)
target_link_libraries(branch_figure PRIVATE es_core)

add_executable(device_figure tools/device_figure.cpp)
target_link_libraries(device_figure PRIVATE es_core)

add_executable(mesh_figure tools/mesh_figure.cpp)
target_link_libraries(mesh_figure PRIVATE es_core)

# ----------------------------------------------------------------------- tests
enable_testing()
foreach(t test_elliptic test_bem test_meniscus test_emission test_beam test_config
           test_regressions test_device_geometry test_boundary_mesh test_vacuum_bem
           test_axisym_fem test_dielectric_device test_reservoir test_capillary test_electrocapillary test_load_projection test_feed test_material_data test_transport)
  add_executable(${t} tests/${t}.cpp)
  target_link_libraries(${t} PRIVATE es_core)
  add_test(NAME ${t} COMMAND ${t})
endforeach()

# The provenance stamp on every figure is produced by python/provenance.py, so
# its test is a Python one.  It builds throw-away git repositories and checks
# that a dirty tree is detected -- a mocked git would test nothing.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME test_provenance
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_provenance.py)
else()
  message(WARNING "Kein Python3 gefunden: test_provenance wird nicht ausgefuehrt.")
endif()
