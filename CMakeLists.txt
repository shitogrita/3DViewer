cmake_minimum_required(VERSION 3.16)
project(CpuStaticViewer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 QUIET COMPONENTS Widgets Gui OpenGLWidgets)
if(Qt6_FOUND)
    add_executable(cpu_viewer
            src/main.cpp
            src/core/matrix.cpp
            src/core/affine_transformation.cpp
            src/core/projection.cpp
            src/view/gl_widget.cpp
            src/render/opengl/gl_render.cpp
            src/render/opengl/gl_shader_loader.cpp
            src/view/test_meshes.cpp
            src/view/interface.cpp

    )
    target_include_directories(cpu_viewer PRIVATE src)
    target_link_libraries(cpu_viewer PRIVATE Qt6::Widgets Qt6::Gui Qt6::OpenGLWidgets)
else()
    find_package(Qt5 REQUIRED COMPONENTS Widgets)
    add_executable(cpu_viewer
            src/main.cpp
            src/core/matrix.cpp
            src/core/affine_transformation.cpp
    )
    target_include_directories(cpu_viewer PRIVATE src view)
    target_link_libraries(cpu_viewer PRIVATE Qt5::Widgets)
endif()
