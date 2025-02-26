cmake_minimum_required(VERSION 3.2)

project(trimesh_viewer)

add_executable(${PROJECT_NAME} main.cpp)

set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)

find_package(cinolib REQUIRED)

target_link_libraries(${PROJECT_NAME} cinolib)
