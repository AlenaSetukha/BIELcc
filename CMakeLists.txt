cmake_minimum_required(VERSION 3.20)

project(BIELCC_LIB CXX)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "release")
endif()

#compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -Wall -Wextra -pedantic -fsanitize=leak, -fsanitize=undefined, -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-flto -march=native")



# Найдите LAPACK и BLAS библиотеки
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


set(BIELCCSOURCE
    electrodynamics/currents.cpp
    electrodynamics/material.cpp
    electrodynamics/num_param.cpp
    electrodynamics/radiation_source.cpp
    geometry/element_geom.cpp
    integrators/integral_param.cpp
    integrators/integrals_analytic.cpp
    io/io_aereco.cpp
    io/io_vtk.cpp
    mesh/cell.cpp
    mesh/point.cpp
    mesh/segment.cpp
    mesh/surf_mesh.cpp)
add_library(bielcc STATIC ${BIELCCSOURCE})



# Указываем, где лежат хедеры, чтобы другие проекты их видели
target_include_directories(bielcc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/electrodynamics)
target_include_directories(bielcc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geometry)
target_include_directories(bielcc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/integrators)
target_include_directories(bielcc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/io)
target_include_directories(bielcc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mesh)



target_link_libraries(bielcc PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}) # Линкуем LAPACK и BLAS библиотеки
if(OPENMP_FOUND)
    target_link_libraries(bielcc PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(bielcc PUBLIC MPI::MPI_CXX) # Связывание MPI с вашим исполняемым файлом
