add_subdirectory(core)
add_subdirectory(render)
add_subdirectory(view)
add_subdirectory(resources)

add_executable(app
        main.cpp
)

target_include_directories(app
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(app
        PRIVATE
        core
        Qt6::Widget
        view

)
