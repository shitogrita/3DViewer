add_library(core STATIC
        matrix.cpp
        affine_transformation.cpp

)

target_include_directories(core
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}         # чтобы было #include "matrix.h"
        ${CMAKE_CURRENT_SOURCE_DIR}/..      # чтобы было #include "core/matrix.h"
)