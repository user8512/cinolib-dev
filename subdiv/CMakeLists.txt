cmake_minimum_required(VERSION 3.2)

set (CMAKE_C_COMPILER "D:/Coding/mingw64/bin/gcc.exe")
set (CMAKE_CXX_COMPILER "D:/Coding/mingw64/bin/g++.exe")

project(CCsubdiv)

add_executable(${PROJECT_NAME} main.cpp)

set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)

set(cinolib_DIR "${PROJECT_SOURCE_DIR}/..")

find_package(cinolib REQUIRED)

target_link_libraries(${PROJECT_NAME} cinolib)