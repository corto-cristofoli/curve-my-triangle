CMAKE_MINIMUM_REQUIRED(VERSION 3.14)

project(curve-my-triangle)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # add by me

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include_directories(${PROJECT_SOURCE_DIR}/deps)
message(STATUS "Deps: ${PROJECT_SOURCE_DIR}/deps")

## MODULES #################################################
include(geometry-central)
include(polyscope)

file(GLOB CCP_SOURCE
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/*/*.cpp"
)

# add_executable(main CCP_SOURCE)
# target_link_libraries(main polyscope geometry-central)
