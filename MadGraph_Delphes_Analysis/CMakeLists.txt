cmake_minimum_required(VERSION 3.16)

project(MadGraph_Delphes_Analysis LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

## Application is always built from App/main.cpp

# Find ROOT with required components
find_package(ROOT REQUIRED COMPONENTS
	Core RIO Hist Tree ROOTDataFrame ROOTVecOps Physics MathCore)

message(STATUS "Found ROOT: ${ROOT_DIR}")
message(STATUS "ROOT include dirs: ${ROOT_INCLUDE_DIRS}")
# message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")

include(${ROOT_USE_FILE})

# Sources and headers
set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SRC_DIR ${PROJECT_ROOT}/Analysis)
set(INC_DIR ${PROJECT_ROOT}/Include)

file(GLOB ANALYSIS_SOURCES
	${SRC_DIR}/*.cpp
)

# Create a library from analysis code
add_library(analysis_core ${ANALYSIS_SOURCES})
target_include_directories(analysis_core
	PUBLIC
		${INC_DIR}
		${ROOT_INCLUDE_DIRS}
)
target_link_libraries(analysis_core
	PUBLIC
		${ROOT_LIBRARIES}
)

# If there are ROOT dictionaries needed later, they can be added here using root_generate_dictionary

# Add the app target (always expected to exist)
set(APP_MAIN ${PROJECT_ROOT}/App/main.cpp)
add_executable(analysis_app ${APP_MAIN})
target_link_libraries(analysis_app PRIVATE analysis_core)
target_include_directories(analysis_app PRIVATE ${INC_DIR})

# Installation rules
# Install the library and executable
install(TARGETS analysis_core analysis_app
    EXPORT MyProjectTargets
    RUNTIME DESTINATION bin           # for executables
    LIBRARY DESTINATION lib           # for shared libraries
    ARCHIVE DESTINATION lib           # for static libraries
    INCLUDES DESTINATION include
)

# Install headers
install(DIRECTORY Include/ DESTINATION include)