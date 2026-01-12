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
#find_package(LAPACK REQUIRED)
#find_package(BLAS REQUIRED)

# Поиск пакета MPI
find_package(MPI REQUIRED)



# Find OpenMP package
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()



include_directories(electrodynamics)
include_directories(em_project)
include_directories(geometry)
include_directories(integrators)
include_directories(io)
include_directories(mesh)
include_directories(parsers)

include_directories(extern/yaml-cpp/include)




add_subdirectory(extern/yaml-cpp)


set(SOURCE
    electrodynamics/currents.cpp
    electrodynamics/get_j.cpp
    electrodynamics/material.cpp
    electrodynamics/num_param.cpp
    electrodynamics/radiation_source.cpp
    em_project/project.cpp
    geometry/element_geom.cpp
    integrators/integral_param.cpp
    integrators/integrals_analytic.cpp
    io/io_aereco.cpp
    io/io_vtk.cpp
    mesh/cell.cpp
    mesh/point.cpp
    mesh/segment.cpp
    parsers/num_param_parser.cpp
    parsers/parse_refr_index.cpp
    parsers/pp_parser.cpp
    parsers/project_parser.cpp)

add_library(pcac_lib STATIC ${SOURCE})



#add_executable(test_PCAC ../examples/mesh/mesh_vtk_test_main0.cpp)
###add_executable(test_PCAC ../examples/parsing/main_parsing_test.cpp)
###add_executable(test_PCAC ../examples/ED_Operators/main_R.cpp)
add_executable(test_PCAC ../examples/field/main_ideal.cpp)


if(OPENMP_FOUND)
    target_link_libraries(test_PCAC PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(test_PCAC PUBLIC pcac_lib)
target_link_libraries(
    test_PCAC PRIVATE yaml-cpp)

#target_link_libraries(test_PCAC PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}) # Линкуем LAPACK и BLAS библиотеки
target_link_libraries(test_PCAC PUBLIC MPI::MPI_CXX) # Связывание MPI с вашим исполняемым файлом
