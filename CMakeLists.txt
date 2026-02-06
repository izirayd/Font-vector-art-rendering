cmake_minimum_required(VERSION 3.20)
project(MeshShaderFont LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
include(FetchContent)

# Dear ImGui (Win32 + DX12 backends)
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.92.5
)
FetchContent_GetProperties(imgui)
if(NOT imgui_POPULATED)
    FetchContent_Populate(imgui)
endif()

add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_win32.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_dx12.cpp
)

target_include_directories(imgui PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)

target_compile_definitions(imgui PUBLIC
    IMGUI_IMPL_WIN32_DISABLE_GAMEPAD
)

target_link_libraries(imgui PUBLIC
    d3d12
    dxgi
    dxguid
    imm32
)

# Find DXC compiler
find_program(DXC_COMPILER dxc
    HINTS
        "$ENV{WindowsSdkDir}/bin/$ENV{WindowsSDKVersion}/x64"
        "$ENV{VULKAN_SDK}/Bin"
        "C:/Program Files (x86)/Windows Kits/10/bin/10.0.22621.0/x64"
        "C:/Program Files (x86)/Windows Kits/10/bin/10.0.22000.0/x64"
        "C:/Program Files (x86)/Windows Kits/10/bin/10.0.26100.0/x64"
)

if(NOT DXC_COMPILER)
    message(FATAL_ERROR "DXC shader compiler not found. Install Windows SDK or set path manually.")
endif()
message(STATUS "DXC compiler found: ${DXC_COMPILER}")

# Shader compilation
set(SHADER_SOURCE "${CMAKE_SOURCE_DIR}/shaders/font_shader.hlsl")
set(SHADER_OUTPUT_DIR "${CMAKE_BINARY_DIR}/shaders")
file(MAKE_DIRECTORY ${SHADER_OUTPUT_DIR})

set(MS_HEADER "${SHADER_OUTPUT_DIR}/font_ms.h")
set(PS_HEADER "${SHADER_OUTPUT_DIR}/font_ps.h")

add_custom_command(
    OUTPUT ${MS_HEADER}
    COMMAND ${DXC_COMPILER} -T ms_6_5 -E MSMain "${SHADER_SOURCE}" -Fh "${MS_HEADER}" -Vn g_MSMain -O3
    DEPENDS ${SHADER_SOURCE}
    COMMENT "Compiling mesh shader..."
)

add_custom_command(
    OUTPUT ${PS_HEADER}
    COMMAND ${DXC_COMPILER} -T ps_6_5 -E PSMain "${SHADER_SOURCE}" -Fh "${PS_HEADER}" -Vn g_PSMain -O3
    DEPENDS ${SHADER_SOURCE}
    COMMENT "Compiling pixel shader..."
)

add_custom_target(CompileShaders DEPENDS ${MS_HEADER} ${PS_HEADER})

# Main executable
add_executable(MeshShaderFont WIN32
    src/main.cpp
    src/font_processor.cpp
    src/camera2d.cpp
)

add_dependencies(MeshShaderFont CompileShaders)

target_include_directories(MeshShaderFont PRIVATE
    ${CMAKE_SOURCE_DIR}/third_party
    ${CMAKE_SOURCE_DIR}/src
    ${SHADER_OUTPUT_DIR}
)

target_link_libraries(MeshShaderFont PRIVATE
    d3d12
    dxgi
    dxguid
    comdlg32
    imgui
)
