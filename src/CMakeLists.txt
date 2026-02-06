# Подмодули
add_subdirectory(core)

# Главный исполняемый файл
add_executable(app
        main.cpp
)

# Общие include-пути (чтобы можно было писать #include "core/matrix.h" и т.п.)
target_include_directories(app
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Линковка с библиотеками модулей
target_link_libraries(app
        PRIVATE
        core
)
