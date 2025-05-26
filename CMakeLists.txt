# CMakeLists.txt

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

#----------------------------------------------------------------------------  
# Setup the project
project(VictreAndDetector)

#----------------------------------------------------------------------------  
# Find Geant4 package, activating all available UI and Vis drivers by default
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------  
# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------  
# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------  
# Add the executable and link it to the Geant4 libraries
add_executable(VictreAndDetector main.cc ${sources} ${headers})
target_link_libraries(VictreAndDetector ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------  
# Copy all .txt and .mac files to the build directory
file(GLOB TXT_FILES "${PROJECT_SOURCE_DIR}/*.txt")
file(GLOB MAC_FILES "${PROJECT_SOURCE_DIR}/*.mac")
set(SCRIPT_FILES ${TXT_FILES} ${MAC_FILES})

if(NOT SCRIPT_FILES)
  message(WARNING "No .txt or .mac files found in ${PROJECT_SOURCE_DIR}.")
else()
  foreach(_script ${SCRIPT_FILES})
    get_filename_component(filename ${_script} NAME)
    configure_file(${_script} ${PROJECT_BINARY_DIR}/${filename} COPYONLY)
    message(STATUS "Copied script: ${filename}")
  endforeach()
endif()

#----------------------------------------------------------------------------  
# Copy ExampleVictrePhantom folder recursively to the build directory
file(COPY ${PROJECT_SOURCE_DIR}/ExampleVictrePhantom DESTINATION ${PROJECT_BINARY_DIR})

#----------------------------------------------------------------------------  
# Add custom build target
add_custom_target(VictreAndDetectorTarget DEPENDS VictreAndDetector)

#----------------------------------------------------------------------------  
# Install the executable
install(TARGETS VictreAndDetector DESTINATION bin)

