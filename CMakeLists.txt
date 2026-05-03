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


include_directories(extensions/epr/)
include_directories(extensions/em_problems/)


### include_directories(extern/yaml-cpp/include)
### add_subdirectory(extern/yaml-cpp)


set(SOURCE
    electrodynamics/currents.cpp
    electrodynamics/material.cpp
    electrodynamics/num_param.cpp
    electrodynamics/radiation_source.cpp
    extensions/em_problems/get_j.cpp
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


add_executable(main0_j ../examples/PEC/main0_PEC_BLAS.cpp)
add_executable(main1_field ../examples/PEC/pp_fieldCalc.cpp)
add_executable(main2_epr ../examples/PEC/pp_epr.cpp)



foreach(TARGET main0_j main1_field main2_epr)
    if(OPENMP_FOUND)
        target_link_libraries(${TARGET} PUBLIC OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(${TARGET} PUBLIC pcac_lib)
    target_link_libraries(${TARGET} PRIVATE lapacke lapack blas)
    target_link_libraries(${TARGET} PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}) # Линкуем LAPACK и BLAS библиотеки
    target_link_libraries(${TARGET} PUBLIC MPI::MPI_CXX) # Связывание MPI с вашим исполняемым файлом
endforeach()


# Custom target для последовательного запуска
add_custom_target(run_all
    COMMAND ${CMAKE_COMMAND} -E echo "Running main0..."
    COMMAND main0_j
    COMMAND ${CMAKE_COMMAND} -E echo "Running main1..."
    COMMAND main1_field
    COMMAND ${CMAKE_COMMAND} -E echo "Running main2..."
    COMMAND main2_epr
    DEPENDS main0_j main1_field main2_epr
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

