cmake_minimum_required(VERSION 3.19)

project(
        CFD
        VERSION 1.0
        LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)

find_package(PythonLibs REQUIRED)
find_package(Python3)
find_package(PkgConfig REQUIRED)

set(Ascent_DIR /home/samrusso/ascent/install/ascent-checkout/lib/cmake/ascent)
set(CONDUIT_DIR /home/samrusso/ascent/install/conduit-v0.9.2/)
set(VTKM_DIR /home/samrusso/ascent/install/vtk-m-v2.1.0/)
set(CAMP_DIR /home/samrusso/ascent/install/camp-v2024.02.1)
set(UMPIRE_DIR /home/samrusso/ascent/install/umpire-2024.02.1)
set(RAJA_DIR /home/samrusso/ascent/install/raja-v2024.02.1)
set(MFEM_DIR /home/samrusso/ascent/install/mfem-4.6)


find_package(Ascent REQUIRED)

include_directories(
        ${PYTHON3_INCLUDE_DIRS}
        "/usr/include/python3.8"
)

find_package(Eigen3 REQUIRED)
#find_package(Boost REQUIRED)
#add_subdirectory(src)

#include_directories( ${Boost_INCLUDE_DIR} )

add_executable(${PROJECT_NAME}_1 assignment_1.cpp)
add_executable(${PROJECT_NAME}_3 assignment_3.cpp)
add_executable(${PROJECT_NAME}_4 assignment_4.cpp)
add_executable(${PROJECT_NAME}_7 assignment_7.cpp)
add_executable(${PROJECT_NAME}_8 assignment_8.cpp include/insitu.hpp)
add_executable(${PROJECT_NAME}_9 assignment_9.cpp)

add_executable(${PROJECT_NAME}_test ascent_test.cpp)


target_link_libraries(
        ${PROJECT_NAME}_1
        python3.8
)
target_link_libraries(
        ${PROJECT_NAME}_3
        python3.8
        Eigen3::Eigen
)

target_link_libraries(
        ${PROJECT_NAME}_4
        python3.8
        Eigen3::Eigen
)

target_link_libraries(
        ${PROJECT_NAME}_7
        python3.8
        Eigen3::Eigen
)

target_link_libraries(
        ${PROJECT_NAME}_8
        python3.8
        Eigen3::Eigen
        ascent::ascent
)

target_link_libraries(
        ${PROJECT_NAME}_test
        PRIVATE
        ascent::ascent
)

target_link_libraries(
        ${PROJECT_NAME}_9
        python3.8
        Eigen3::Eigen
)
#target_link_libraries(main PRIVATE ifd_lib Eigen3::Eigen ${Boost_LIBRARIES})