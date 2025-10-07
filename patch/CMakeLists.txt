cmake_minimum_required(VERSION 3.18)

project(patch LANGUAGES CXX CUDA)

# 语言标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# 按你的 GPU 修改/追加架构（示例：75/86/89/90 等；可写多个）
set(CMAKE_CUDA_ARCHITECTURES 86)

set(CMAKE_CUDA_FLAGS_INIT "-Xcompiler=/D_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH")

# 你的源码 + 新增的 CUDA 源文件
add_executable(${PROJECT_NAME}
    patch.cpp
    kmeans.cpp
    cc_subdiv.cu        # <--- 在这里加入你的 .cu 文件
)

# 分离编译/设备符号解析对多 TU/模板更稳
set_target_properties(${PROJECT_NAME} PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

# 可选：需要 GPU 端 lambda 时（老些的 CUDA 版本）
target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
)

# 你原有的 cinolib 配置保持不变
set(cinolib_DIR "${PROJECT_SOURCE_DIR}/..")
set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
find_package(cinolib REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE cinolib CUDA::cudart)

target_compile_definitions(${PROJECT_NAME} PRIVATE
    DATA_PATH="../../../data"
)
