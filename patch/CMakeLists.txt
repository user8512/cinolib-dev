cmake_minimum_required(VERSION 3.18)

project(patch LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

set(CMAKE_CUDA_ARCHITECTURES 86)

add_executable(${PROJECT_NAME}
    main.cpp
    patch.cpp
    kmeans.cpp
    cc_subdiv_gpu.cu   
    #cc_subdiv_gpu_prototype.cu       
 )

# 对多翻译单元/模板友好
set_target_properties(${PROJECT_NAME} PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

set(cinolib_DIR "${PROJECT_SOURCE_DIR}/..")
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
find_package(cinolib REQUIRED)
find_package(CUDAToolkit REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE cinolib CUDA::cudart)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")

target_compile_definitions(${PROJECT_NAME} PRIVATE
    DATA_PATH="../../../data"
)
