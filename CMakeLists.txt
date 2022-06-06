cmake_minimum_required(VERSION 3.14)

project(BatchProcessorCmd LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core)

add_compile_definitions(QT_NO_KEYWORDS)

add_executable(BatchProcessorCmd
  main.cpp
  batchconfiguration.h
  batchconfiguration.cpp
  geometry/ConcurrentPool.h
  geometry/mesh.cpp
  geometry/mesh.h
  geometry/node.cpp
  geometry/node.h
  geometry/element.cpp
  geometry/element.h
  geometry/czinsertiontool.cpp
  geometry/czinsertiontool.h
  geometry/cohesivezone.cpp
  geometry/cohesivezone.h


)
target_link_libraries(BatchProcessorCmd Qt${QT_VERSION_MAJOR}::Core)
target_include_directories(BatchProcessorCmd PRIVATE geometry)


# Gmsh
target_include_directories(BatchProcessorCmd PRIVATE "C:\\Users\\s\\Projects\\libraries\\gmsh-4.9.5-source\\api")
target_link_directories(BatchProcessorCmd PRIVATE "C:\\Users\\s\\Projects\\libraries\\gmsh-build\\Release")
target_link_libraries(BatchProcessorCmd gmsh)


# TBB and MKL
if (WIN32)
    target_include_directories(BatchProcessorCmd PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\tbb\\2021.5.2\\include")
    target_link_directories(BatchProcessorCmd PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\tbb\\2021.5.2\\lib\\intel64\\vc14")
    target_link_libraries(BatchProcessorCmd tbb)

else()
    target_include_directories(BatchProcessorCmd PRIVATE /usr/include/mkl)
    target_link_libraries(BatchProcessorCmd tbb mkl_intel_lp64 mkl_gnu_thread mkl_core gomp m dl)
    target_link_libraries(BatchProcessorCmd ${TBB_IMPORTED_TARGETS})
endif ()

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE HINTS "C:\\Users\\s\\Projects\\libraries\\eigen-build")
target_link_libraries(BatchProcessorCmd Eigen3::Eigen)

# spdlog
if (WIN32)
    target_link_libraries(BatchProcessorCmd spdlog)
    target_include_directories(BatchProcessorCmd PRIVATE "C:\\Users\\s\\Projects\\libraries\\spdlog\\include")
    target_link_directories(BatchProcessorCmd PRIVATE "C:\\Users\\s\\Projects\\libraries\\spdlog-build\\Release")
else()
    target_link_libraries(BatchProcessorCmd spdlog::spdlog)
endif()
