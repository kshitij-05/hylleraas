cmake_minimum_required(VERSION 3.25)
project(hylleraas)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)
set(Boost_USE_MULTITHREADED ON)

find_package(Eigen3 REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(hylleraas
            hylleraas.cpp)
    target_link_libraries(hylleraas ${Boost_LIBRARIES})
endif()

target_link_libraries(hylleraas Eigen3::Eigen)

