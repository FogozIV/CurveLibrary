cmake_minimum_required(VERSION 3.20)
project(CurveLibrary)

set(CMAKE_CXX_STANDARD 20)

# --------------------------
# PYBIND11 + CURVE LIBRARY
# --------------------------

add_subdirectory(lib/EurobotCoreMath)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
add_subdirectory(extern/pybind11)  # <- you can mirror PacketHandling structure here

# Gather your source files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS include/*.h)

# Create the static library
add_library(CurveLibrary STATIC ${SOURCES} ${HEADERS})
target_link_libraries(CurveLibrary
        PUBLIC EurobotCoreMath
)
target_include_directories(CurveLibrary
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# --------------------------
# PYTHON BINDINGS MODULE
# --------------------------

# Bindings sources
file(GLOB_RECURSE BINDINGS CONFIGURE_DEPENDS bindings/*.cpp)

pybind11_add_module(curve_library ${BINDINGS})

target_link_libraries(curve_library PRIVATE CurveLibrary)

target_include_directories(curve_library
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Windows .pyd setup
set_target_properties(curve_library PROPERTIES PREFIX "" SUFFIX ".pyd")

# Install compiled module & __init__.py
install(TARGETS curve_library
        LIBRARY DESTINATION curve_library
)
install(FILES curve_library/__init__.py DESTINATION curve_library)

# --------------------------
# STUB GENERATION
# --------------------------

find_program(PYBIND11_STUBGEN_EXECUTABLE pybind11-stubgen)
message(STATUS "Using Python executable: ${Python3_EXECUTABLE}")
if (PYBIND11_STUBGEN_EXECUTABLE)
    add_custom_command(
            TARGET curve_library
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:curve_library>
            ${CMAKE_CURRENT_SOURCE_DIR}/curve_library/curve_library.pyd
            COMMENT "Copying .pyd to curve_library/ for stub generation"
    )
    add_custom_command(
            TARGET curve_library
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}
            ${Python3_EXECUTABLE} -m pybind11_stubgen curve_library
            --output-dir ${CMAKE_CURRENT_BINARY_DIR}/generated_stubs
            COMMENT "Generating .pyi stub file for curve_library"
    )
    #add_custom_command(
    #       TARGET curve_library
    #       POST_BUILD
    #       COMMAND ${CMAKE_COMMAND} -E remove
    #       ${CMAKE_CURRENT_SOURCE_DIR}/curve_library/curve_library.pyd
    #       COMMENT "Cleaning up copied curve_library.pyd"
    #)

    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/generated_stubs/curve_library/curve_library.pyi
            DESTINATION curve_library
    )
else()
    message(WARNING "pybind11-stubgen not found; .pyi stub will not be generated")
endif()
