cmake_minimum_required(VERSION 3.20)

project(PCAC_LIB CXX)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "release")
endif()

#compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -Wall -Wextra -pedantic -fsanitize=leak, -fsanitize=undefined, -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-flto -march=native")


### Компиляция с флагами: $ cmake -DCMAKE_BUILD_TYPE=Debug ..
### set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -flto -march=native")



# Найдите LAPACK и BLAS библиотеки
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
include_directories(${LAPACK_INCLUDE_DIRS} ${BLAS_INCLUDE_DIRS})






include_directories(electrodynamics) 
include_directories(geometry)
include_directories(integrators)
include_directories(io)
include_directories(mesh)






set(SOURCE
    electrodynamics/currents.cpp
    electrodynamics/material.cpp
    electrodynamics/num_param.cpp
    electrodynamics/radiation_source.cpp
    geometry/element_geom.cpp
    integrators/integral_param.cpp
    integrators/integrals_analytic.cpp
    io/io_vtk.cpp
    io/io_aereco.cpp
    mesh/cell.cpp
    mesh/point.cpp
    mesh/segment.cpp
    mesh/surf_mesh.cpp)
add_library(pcac_lib STATIC ${SOURCE})





add_executable(main0 ../examples/ED_Operators/main0_R.cpp)


target_link_libraries(main0 PUBLIC pcac_lib)
target_link_libraries(main0 PRIVATE lapacke lapack blas)
target_link_libraries(main0 PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}) # Линкуем LAPACK и BLAS библиотеки
target_link_libraries(main0 PUBLIC MPI::MPI_CXX) # Связывание MPI с вашим исполняемым файлом

