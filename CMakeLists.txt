cmake_minimum_required(VERSION 3.10)
project(MyProject)

set(CMAKE_CXX_STANDARD 17)
include_directories(${CMAKE_SOURCE_DIR}/eigen-3.4.0)

file(GLOB SOURCES "src/*.cpp")


foreach(SOURCE ${SOURCES})
    get_filename_component(EXEC_NAME ${SOURCE} NAME_WE)
    add_executable(${EXEC_NAME} ${SOURCE})
endforeach()
