cmake_minimum_required(VERSION 3.2)

project(patch)

add_executable(${PROJECT_NAME} patch.cpp)

set(cinolib_DIR "${PROJECT_SOURCE_DIR}/..")

set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)

find_package(cinolib REQUIRED)

target_link_libraries(${PROJECT_NAME} cinolib)

target_compile_definitions(patch PRIVATE
    DATA_PATH="../../../data"
)